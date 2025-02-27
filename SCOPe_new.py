import os
import requests
import tarfile
import shutil
from io import BytesIO, StringIO
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from MMseqs import MMseqs2

PDBSTYLE_URL_TEMPLATE = "https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-2.08-{}.tgz"
FASTA_URL = "https://scop.berkeley.edu/downloads/scopeseq-2.08/astral-scopedom-seqres-gd-all-2.08-stable.fa"

def create_directories(base_dir: Path) -> dict:
    """
    Create and return a dictionary of required directories.
    """
    dirs = {
        "main": base_dir,
        "input": base_dir / "input",
        "pdbs": base_dir / "pdbs",
        "temp": base_dir / "temp",
        "tarballs": base_dir / "tarballs",
        "data": base_dir / "data",
        "final_output": base_dir / "final_output",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs

def download_and_extract_tarball(url: str, tarballs_dir: Path) -> dict:
    """
    Download a tarball from a URL (caching if already present) and extract text files.
    Returns a dict mapping member filename to file content.
    """
    tarball_filename = url.split("/")[-1]
    tarball_path = tarballs_dir / tarball_filename

    if tarball_path.exists():
        print(f"Using cached tarball: {tarball_filename}")
        content = tarball_path.read_bytes()
    else:
        print(f"Downloading {tarball_filename}")
        try:
            response = requests.get(url, verify=False)
            response.raise_for_status()
        except Exception as e:
            print(f"Error downloading {url}: {e}")
            return {}
        content = response.content
        tarball_path.write_bytes(content)
        print(f"Saved tarball to {tarball_path}")

    extracted_files = {}
    with tarfile.open(fileobj=BytesIO(content), mode="r:gz") as tar:
        for member in tar.getmembers():
            if member.isfile():
                file_obj = tar.extractfile(member)
                if file_obj:
                    try:
                        file_content = file_obj.read().decode('utf-8')
                        extracted_files[member.name] = file_content
                    except Exception as e:
                        print(f"Error decoding file {member.name}: {e}")
    return extracted_files

def write_pdbstyle_data_to_disk(pdbstyle_data: dict, output_dir: Path) -> None:
    """
    Write extracted pdbstyle files to the specified output directory.
    """
    for name, content in pdbstyle_data.items():
        file_path = output_dir / name
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_path.write_text(content)
    print(f"All files have been written to {output_dir}")

def download_pdbstyle(dirs: dict, write: bool = True) -> list:
    """
    Download multiple pdbstyle tarballs, extract .ent files,
    optionally write them to disk, and return a list of pdbstyle file paths.
    """
    pdbstyle_tree = []
    data_dir = dirs["data"]
    for i in tqdm(range(1, 9), desc="Downloading pdbstyle tarballs"):
        print(f"Processing pdbstyle-2.08-{i}")
        url = PDBSTYLE_URL_TEMPLATE.format(i)
        extracted_files = download_and_extract_tarball(url, dirs["tarballs"])

        ent_files = [name for name in extracted_files if name.endswith('.ent')]
        pdbstyle_tree.extend(ent_files)
        
        pdb_output_dir = data_dir / "pdbs" / "SCOPe" / f"pdbstyle-2.08-{i}"
        pdb_output_dir.mkdir(parents=True, exist_ok=True)
        if write:
            write_pdbstyle_data_to_disk(extracted_files, output_dir=pdb_output_dir)

    tree_file = dirs["main"] / "pdbstyle_tree.txt"
    if write:
        with tree_file.open("w") as f:
            for line in pdbstyle_tree:
                f.write(f"{line}\n")
        print(f"Pdbstyle tree saved to {tree_file}")
    return pdbstyle_tree

def download_and_extract_class_sequences(url: str, min_length: int, max_length: int,
                                           dirs: dict, save_to_disk: bool = False) -> dict:
    """
    Download a FASTA file from the given URL, filter sequences by standard amino acids
    and length, and optionally save the FASTA locally.
    Returns a dict mapping class letter to a list of SeqRecord objects.
    """
    print("Downloading and extracting sequences ...")
    response = requests.get(url, verify=False)
    response.raise_for_status()
    fasta_content = response.text

    fasta_file_path = dirs["input"] / "astral-scopedom-seqres-gd-all-2.08-stable.fa"
    fasta_file_path.parent.mkdir(parents=True, exist_ok=True)
    if save_to_disk:
        fasta_file_path.write_text(fasta_content)
        print(f"FASTA file saved to {fasta_file_path}")

    fasta_io = StringIO(fasta_content)
    classes = {}
    valid_AAs = set("ACDEFGHIKLMNPQRSTVWY")
    expected_classes = "abcdeg"

    for record in SeqIO.parse(fasta_io, "fasta"):
        try:
            class_letter = record.description.split()[1][0].lower()
        except IndexError:
            print(f"Skipping record {record.id} due to unexpected description format.")
            continue

        seq = str(record.seq).upper()
        if set(seq).issubset(valid_AAs) and min_length <= len(seq) <= max_length and class_letter in expected_classes:
            if class_letter not in classes:
                classes[class_letter] = []
                print(f"Adding class {class_letter}")
            classes[class_letter].append(SeqRecord(Seq(seq), id=record.id, description=record.description))
    return classes

def correlate_structures_and_sequences(pdbstyle_tree: list, sequences: dict) -> dict:
    """
    Correlate structure file names (from pdbstyle_tree) with sequence records by matching identifiers.
    Returns a dict mapping structure id to a dict with "structure" (the file name) and "sequence" (SeqRecord).
    """

    id_to_path = {fname.split('/')[-1].split('.')[0]: Path(fname) for fname in pdbstyle_tree}
    correlated_data = {}

    for seq_id, seq_record in sequences.items():
        if seq_id in id_to_path:
            correlated_data[seq_id] = {
                "structure": id_to_path[seq_id],
                "sequence": seq_record
            }
        else:
            print(f"No matching structure for sequence ID {seq_id}")
    return correlated_data

def cleanup_intermediate_dirs(dirs: dict) -> None:
    """
    Remove all intermediate directories and files, keeping only final outputs.
    """
    targets = [
        dirs["input"],
        dirs["temp"],
        dirs["data"],
        dirs["main"] / "pdbstyle_tree.txt"
    ]
    for path in targets:
        if path.exists():
            try:
                if path.is_file():
                    path.unlink()
                    print(f"Deleted file: {path}")
                else:
                    shutil.rmtree(path)
                    print(f"Deleted directory: {path}")
            except Exception as e:
                print(f"Error deleting {path}: {e}")

def run_pipeline(min_length: int = 20, max_length: int = 100,
                 base_dir: Path = Path("SCOPe_database"),
                 cleanup: bool = False) -> int:
    """
    Run the full pipeline:
      - Set up directories.
      - Download and extract pdbstyle tarballs.
      - Download the FASTA file and filter sequences by min/max length.
      - Separate sequences into groups (abcde and g).
      - Write intermediate FASTA files.
      - Use MMseqs2 to get representative sequences.
      - Correlate sequences with structure files.
      - Write final output files.
      - Optionally clean up intermediate files.
    
    Parameters:
      min_length: Minimum sequence length to keep.
      max_length: Maximum sequence length to keep.
      base_dir: Base directory for the pipeline files.
      cleanup: If True, remove intermediate directories after processing.
    
    Returns:
      0 upon completion.
    """
    dirs = create_directories(base_dir)

    pdbstyle_tree = download_pdbstyle(dirs, write=True)
    tree_file = dirs["main"] / "pdbstyle_tree.txt"
    if tree_file.exists():
        with tree_file.open("r") as f:
            pdbstyle_tree = [line.strip() for line in f if line.strip()]

    classes = download_and_extract_class_sequences(FASTA_URL, min_length, max_length, dirs, save_to_disk=True)

    abcde_sequences = []
    g_sequences = []
    for class_letter, seq_list in classes.items():
        if class_letter in "abcde":
            abcde_sequences.extend(seq_list)
        elif class_letter == "g":
            g_sequences.extend(seq_list)
        else:
            print(f"Skipping unexpected class: {class_letter}")

    abcde_temp_fasta = dirs["temp"] / "abcde_sequences.fasta"
    abcde_temp_fasta.parent.mkdir(parents=True, exist_ok=True)
    with abcde_temp_fasta.open("w") as f:
        SeqIO.write(abcde_sequences, f, "fasta")
    print(f"abcde sequences written to {abcde_temp_fasta}")

    g_temp_fasta = dirs["temp"] / "g_sequences.fasta"
    g_temp_fasta.parent.mkdir(parents=True, exist_ok=True)
    with g_temp_fasta.open("w") as f:
        SeqIO.write(g_sequences, f, "fasta")
    print(f"g sequences written to {g_temp_fasta}")

    mmseqs2_api = MMseqs2(threads=os.cpu_count() - 2, cleanup=True)
    representatives_abcde = mmseqs2_api.fasta2representativeseq(
        str(abcde_temp_fasta), cov=0.7, iden=0.3, cov_mode=0
    )
    print(f"Found {len(representatives_abcde)} representative sequences for abcde.")

    representatives_g = mmseqs2_api.fasta2representativeseq(
        str(g_temp_fasta), cov=0.7, iden=0.3, cov_mode=0
    )
    print(f"Found {len(representatives_g)} representative sequences for g.")

    correlated_data_abcde = correlate_structures_and_sequences(pdbstyle_tree, representatives_abcde)
    print(f"Correlated {len(correlated_data_abcde)} abcde sequences with structure files.")
    correlated_data_g = correlate_structures_and_sequences(pdbstyle_tree, representatives_g)
    print(f"Correlated {len(correlated_data_g)} g sequences with structure files.")

    output_base = dirs["final_output"]
    output_base.mkdir(parents=True, exist_ok=True)

    pdb_output_dir_abcde = output_base / "pdb_files" / "abcde"
    for struct_id, data in correlated_data_abcde.items():
        out_dir = pdb_output_dir_abcde / struct_id
        out_dir.mkdir(parents=True, exist_ok=True)
        structure_file = out_dir / f"{struct_id}.ent"
        structure_file.write_text(f"Structure file path: {data['structure']}")
    
    pdb_output_dir_g = output_base / "pdb_files" / "g"
    for struct_id, data in correlated_data_g.items():
        out_dir = pdb_output_dir_g / struct_id
        out_dir.mkdir(parents=True, exist_ok=True)
        structure_file = out_dir / f"{struct_id}.ent"
        structure_file.write_text(f"Structure file path: {data['structure']}")

    abcde_final_file = output_base / "abcde_representatives.fasta"
    g_final_file = output_base / "g_representatives.fasta"
    abcde_final_seqs = [data["sequence"] for data in correlated_data_abcde.values()]
    g_final_seqs = [data["sequence"] for data in correlated_data_g.values()]

    with abcde_final_file.open("w") as f:
        SeqIO.write(abcde_final_seqs, f, "fasta")
    with g_final_file.open("w") as f:
        SeqIO.write(g_final_seqs, f, "fasta")

    print("Multifasta files for abcde and g classes have been written.")
    print("Pipeline complete.")

    if cleanup:
        cleanup_intermediate_dirs(dirs)
    return 0

if __name__ == "__main__":
    run_pipeline(min_length=20, max_length=100)
