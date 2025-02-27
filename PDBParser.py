from pathlib import Path
from collections import defaultdict, OrderedDict
from typing import Optional, Any, Dict
import numpy as np
import logging
import functools
import inspect

import re
import os 

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import argparse
from tqdm import tqdm
import concurrent.futures

import polars as pl 
from glob import glob 



def log_on_exception(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            # Try to extract a local variable called 'log_messages' from one of the frames in the traceback.
            captured_logs = []
            tb = e.__traceback__
            for frame_info in inspect.getinnerframes(tb):
                local_vars = frame_info.frame.f_locals
                if 'log_messages' in local_vars:
                    captured_logs = local_vars['log_messages']
                    break
            # Build an error message including the captured logs.
            log_text = "\n".join(str(msg) for msg in captured_logs) if captured_logs else "No log messages captured."
            new_msg = f"Exception in {func.__name__}:\n{log_text}\nOriginal error: {e}"
            # Option 1: Modify the exception message (if the exception accepts string args).
            e.args = (new_msg,) + e.args[1:]
            # Option 2: You could also print it or send it elsewhere here.
            raise e
    return wrapper

    
def strip_dict(d):

    keys = list(d.keys())
    for key in keys:
        if d[key] == "X":
            del d[key]
        else:
            break
    for key in reversed(keys):
        if key in d and d[key] == "X":
            del d[key]
        else:
            break


class SlicedOrderedDict(OrderedDict):
    """
    Just a custom OD that supports slicing 
    """
    def __getitem__(self, key):
        if isinstance(key, slice):
            keys = list(self.keys())[key]
            return SlicedOrderedDict((k, self[k]) for k in keys)
        else:
            return OrderedDict.__getitem__(self, key)

class OPMPDBStruct:

    def __init__(self, 
                 logger,
                 min_length:int = 20,
                 max_length:int = 100,
                 margin:float = 5.0, 
                 inner_margin:float = 0.0,
                 min_segment_length:int = 15, 
                 gaps:int = 1
                ):

        self.protein_name: str = ""
        self.full: Dict[str, Any] = {}
        self.CA: Dict[str, Any] = {}
        self.sequences: Dict[str, str] = {}  
        self.membrane_coord: Optional[np.ndarray] = None
        self.protein_length: Dict[str, int] = {}
        self.in_membrane: Dict[str, str] = {}  
        self.in_margin: Dict[str, str] = {}
        self.tm_segments: Dict[str, list] = {} 

        self.logger = logger
        self.min_length = min_length
        self.max_length = max_length 
        self.margin = margin 
        self.inner_margin = inner_margin 
        self.min_segment_length = min_segment_length 
        self.gaps = gaps

    def __repr__(self):
        return (f"PDBStruct(protein_name={self.protein_name}, sequences={self.sequences}, "
                f"protein_length={self.protein_length})")
    
    def analyze(self):

        self._binarize_transmembrane()
        self._define_tm_segments()

        return self._extract_elongated_segments()
    
    @log_on_exception
    def _binarize_transmembrane(self):
        """
        Create binary strings for each chain indicating if a residue is within the membrane region.
        
        A residue is flagged as '1' (in the membrane) if its x and y coordinates are within the 
        membrane boundaries and its z coordinate is between min_z+inner_margin and max_z-inner_margin.
        An expanded condition with an additional margin flag is also computed.
        
        Parameters:
            margin (float): Additional margin in angstroms for the z-coordinate.
            inner_margin (float): Reduction in z-boundaries for the stricter condition,
                                better left untouched

        Returns:
            tuple: Two dictionaries mapping chain IDs to their binary strings for:
                   - strict membrane inclusion (in_membrane_binaries)
                   - inclusion within the expanded margin (in_margin_binaries)
        """
        if self.membrane_coord is None or self.membrane_coord.size == 0:
            raise ValueError("No membrane coordinates available to binarize.")
        
        log_messages = []
        
        min_z_membrane = np.min(self.membrane_coord[:, 2])
        max_z_membrane = np.max(self.membrane_coord[:, 2])
        min_x_membrane = np.min(self.membrane_coord[:, 0])
        max_x_membrane = np.max(self.membrane_coord[:, 0])
        min_y_membrane = np.min(self.membrane_coord[:, 1])
        max_y_membrane = np.max(self.membrane_coord[:, 1])
        
        in_membrane_binaries = {chain_id: "" for chain_id in self.CA.keys()}
        in_margin_binaries = {chain_id: "" for chain_id in self.CA.keys()}
        
        for chain_id, residues in self.CA.items():
            log_messages.append(f"Binarizing chain {chain_id} ...")
            sorted_resnums = sorted(residues.keys())
            if not sorted_resnums:
                continue
            last_res_number = sorted_resnums[0]
            log_messages.append(f"First residue number: {last_res_number}")
            # Subtract one to handle potential gaps 
            last_res_number = last_res_number - 1
            
            for current_res_number in sorted_resnums:
                data = residues[current_res_number]
                x, y, z = data["coord"]
                
                # Check for missing residues and fill with zeros
                if int(current_res_number) != last_res_number + 1:
                    missing_count = int(current_res_number) - last_res_number - 1
                    in_membrane_binaries[chain_id] += "0" * missing_count
                    in_margin_binaries[chain_id] += "0" * missing_count
                    log_messages.append(f"Found missing amino acids from {last_res_number} to {current_res_number}")
                
                # Determine if the coordinates fall within the membrane boundaries
                # x and y just allow to detect pathological pdbs where AA are out 
                # of the membrane 
                is_x = (min_x_membrane <= x <= max_x_membrane)
                is_y = (min_y_membrane <= y <= max_y_membrane)
                is_z = (min_z_membrane + self.inner_margin <= z <= max_z_membrane - self.inner_margin)
                is_z_plus_margin = (min_z_membrane - self.margin <= z <= max_z_membrane + self.margin)
                
                in_membrane = "1" if is_x and is_y and is_z else "0"
                in_margin = "1" if is_x and is_y and is_z_plus_margin else "0"
                
                # Update the residue data with binary flags
                data.update({"in_membrane": in_membrane, "in_margin": in_margin})
                in_membrane_binaries[chain_id] += in_membrane
                in_margin_binaries[chain_id] += in_margin
                
                last_res_number = int(current_res_number)
            
            log_messages.append(f"Binary membrane for chain {chain_id}: {in_membrane_binaries[chain_id]}")
            log_messages.append(f"Binary margin for chain {chain_id}: {in_margin_binaries[chain_id]}")
        
        if self.logger:
            self.logger.info("\n".join(log_messages))
            
        self.in_membrane = in_membrane_binaries
        self.in_margin = in_margin_binaries

    @log_on_exception
    def _define_tm_segments(self) -> Dict[str, list]:
        """
        Identify transmembrane segments from the self.in_membrane binary strings.

        This method assumes:
        - self.in_membrane is a dict mapping chain IDs to binary strings (each character is either "1" or "0")
            where "1" indicates a residue is predicted to be in the membrane.
        - self.CA is a dict mapping chain IDs to another dict whose keys are the actual (biological)
            residue numbers (typically 1-based) of CA atoms.
        
        The binary string is iterated using 1-based indices (matching the residue order in the binary string).
        For each chain, a segment of contiguous "1"s is considered a transmembrane region if its length is at least
        min_segment_length. The segment's start and end indices (which are computed relative to the binary string)
        are then adjusted using the first residue number from the CA data so that the reported indices correspond
        to the PDB residue numbering.

        Parameters:
        min_segment_length (int): The minimum number of consecutive '1's to qualify as a TM segment.

        Returns:
        A dict mapping each chain to a list of tuples (start_res, end_res, length), where start_res and end_res are
        the adjusted (biological) residue numbers for the transmembrane segments.
        """
        log_messages = []
        tm_segments = {chain_id: [] for chain_id in self.in_membrane.keys()}
        log_messages.append("Searching for transmembrane segments ...")
        log_messages.append(f"Chains in input: {list(self.in_membrane.keys())}")

        for chain_id, binary_str in self.in_membrane.items():
            if chain_id not in self.CA or not self.CA[chain_id]:
                log_messages.append(f"No CA data for chain {chain_id}, skipping.")
                continue

            # Residues do not always start at 1
            first_residue = sorted(self.CA[chain_id].keys())[0]
            start_index = None  

            for pos, bit in enumerate(binary_str, start=1):
                if bit == "1": 
                    if start_index is None:
                        start_index = pos 
                else: 
                    if start_index is not None:
                        segment_length = pos - start_index
                        if segment_length >= self.min_segment_length:
                            tm_segments[chain_id].append((start_index, pos - 1, segment_length))
                        start_index = None  # Reset for the next segment

            if start_index is not None:
                segment_length = len(binary_str) - start_index + 1
                if segment_length >= self.min_segment_length:
                    tm_segments[chain_id].append((start_index, len(binary_str), segment_length))

 
            adjusted_segments = []
            for seg_start, seg_end, seg_length in tm_segments[chain_id]:
                adjusted_start = first_residue + seg_start - 1
                adjusted_end = first_residue + seg_end - 1
                adjusted_segments.append((adjusted_start, adjusted_end, seg_length))
            tm_segments[chain_id] = adjusted_segments

            log_messages.append(f"TM segments for chain {chain_id}: {tm_segments[chain_id]}")

        if self.logger:
            self.logger.info("\n".join(log_messages))

        self.tm_segments = tm_segments

    @log_on_exception
    def _extract_elongated_segments(self) -> dict:

        log_messages = []
        records = []
        records_shorts = []
        structures = {chain_id: {} for chain_id in self.tm_segments}
        structures_shorts = {chain_id: {} for chain_id in self.tm_segments}
        protein_name = self.protein_name

        log_messages.append(f"Extracting elongated sequences with random desired lengths...\nChains: {list(self.tm_segments.keys())}")

        for chain_id in self.tm_segments:

            if not self.tm_segments[chain_id]:
                log_messages.append(f"Chain {chain_id} has no TM segments; skipping.")
                continue

            sorted_res = sorted(self.CA[chain_id].keys())
            margin_mask = "".join('1' if self.CA[chain_id][res]["in_margin"] == "1" else '0' for res in sorted_res)
            margin_mask = self._smooth_margin(margin_mask)
            chain_seq = "".join(self.CA[chain_id][res]["res_name"] for res in sorted_res)
            index_to_res = {i+1: res for i, res in enumerate(sorted_res)}
            res_to_index = {res: i+1 for i, res in enumerate(sorted_res)}
            
            log_messages.append(f"Chain {chain_id}: margin_mask = {margin_mask}")
            
            num_segments = len(self.tm_segments[chain_id])
            for i, segment in enumerate(self.tm_segments[chain_id]):
                # Convert core TM boundaries (tm_start, tm_end) from residue numbers to sequence positions.

                log_messages.append(f"Segments : {segment}")

                tm_start, tm_end, _ = segment
                try:
                    pos_tm_start = res_to_index[tm_start]
                    pos_tm_end = res_to_index[tm_end]
                except KeyError:
                    log_messages.append(f"Chain {chain_id}: TM segment {tm_start}-{tm_end} not found in CA data.")
                    continue
                
                current_length = pos_tm_end - pos_tm_start + 1
                
                desired_length = np.random.randint(self.min_length, self.max_length)
                log_messages.append(f"Chain {chain_id}, segment {i}: core from {tm_start} to {tm_end} (pos {pos_tm_start}-{pos_tm_end}, length {current_length}); desired length: {desired_length}")
                
                if desired_length <= current_length:
                    pos_elong_start = pos_tm_start
                    pos_elong_end = pos_tm_end
                else:
                    extra_needed = desired_length - current_length
                    upstream_extra = np.random.randint(0, extra_needed)
                    downstream_extra = extra_needed - upstream_extra

                    # --- BACKWARD EXTENSION ---
                    # Which residues are in the margin before the segment ?
                    left_mask = margin_mask[:pos_tm_start] 
                    left_block = re.search(r"(1+)$", left_mask)
                    if left_block:
                        available_left = pos_tm_start - left_block.start()  # how many ones are available
                        actual_upstream = min(available_left, upstream_extra)
                    else:
                        actual_upstream = 0
                    pos_elong_start = pos_tm_start - actual_upstream
                    
                    # --- FORWARD EXTENSION ---
                    # Whihc ones are after ?
                    right_mask = margin_mask[pos_tm_end:]  
                    right_block = re.match(r"(1+)", right_mask)
                    if right_block:
                        available_right = right_block.end()  
                        actual_downstream = min(available_right, downstream_extra)
                    else:
                        actual_downstream = 0
                    pos_elong_end = pos_tm_end + actual_downstream

                # Map back from position in sequence index to absolute value of res
                new_start = index_to_res.get(pos_elong_start, tm_start)
                new_end = index_to_res.get(pos_elong_end, tm_end)
                
                # Extract sequences.
                full_seq = chain_seq[pos_elong_start - 1: pos_elong_end]  
                core_seq = chain_seq[pos_tm_start - 1: pos_tm_end]        
                
                fasta_key = f"{protein_name}_{chain_id}_{i+1}_elong"
                fasta_key_short = f"{protein_name}_{chain_id}_{i+1}_short"
                
                # Extract pdb lines 
                full_struct = {res: self.full[chain_id][res] for res in sorted_res if pos_elong_start <= res_to_index[res] <= pos_elong_end}
                full_struct_short = {res: self.full[chain_id][res] for res in sorted_res if pos_tm_start <= res_to_index[res] <= pos_tm_end}

                # Trim Xss
                full_seq = full_seq.strip("X")
                core_seq = core_seq.strip("X")
                
                strip_dict(full_struct)
                strip_dict(full_struct_short)
                
                if ("X" * self.gaps not in full_seq) and full_seq:
                    record = SeqRecord(Seq(full_seq), id=fasta_key, description="")
                    records.append(record)
                    if chain_id not in structures:
                        structures[chain_id] = {}
                    structures[chain_id][fasta_key] = full_struct
                
                if ("X" * self.gaps not in core_seq) and core_seq:
                    record_short = SeqRecord(Seq(core_seq), id=fasta_key_short, description="")
                    records_shorts.append(record_short)
                    if chain_id not in structures_shorts:
                        structures_shorts[chain_id] = {}
                    structures_shorts[chain_id][fasta_key_short] = full_struct_short
                
                log_messages.append(
                    f"Chain {chain_id}, segment {i}: pos_tm_start={pos_tm_start}, pos_tm_end={pos_tm_end}, "
                    f"pos_elong_start={pos_elong_start}, pos_elong_end={pos_elong_end}, new_start={new_start}, new_end={new_end}, "
                    f"full_seq={full_seq}, core_seq={core_seq}"
                )
        
        if self.logger:
            self.logger.info("\n".join(log_messages))
        
        return {
            "sequences": records,
            "sequences_shorts": records_shorts,
            "structures": structures,
            "structures_shorts": structures_shorts
        }

    def _smooth_margin(self, binary_str) -> str:

        # 0s flanked by two 1s
        # Allow stretches of 0s up to tolerance parameter
        pattern = r'(?<=1)(0{1,' + str(self.gaps) + r'})(?=1)'
        # Replace the matched block of zeros by 1s
        smoothed = re.sub(pattern, lambda m: '1' * len(m.group(0)), binary_str)
        return smoothed

class OPMPDBParser:
    aa_dict = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    def __init__(self, 
                 file_path: Path, 
                 secondary_structure_path: Optional[Path] = None, 
                 logger: Optional[Any] = None, 
                 verbose: bool = False,
                 min_length:int = 20,
                 max_length:int = 100,
                 margin:float = 5.0, 
                 inner_margin:float = 0.0,
                 min_segment_length:int = 15, 
                 gaps:int = 1
                 ):
        
        if not isinstance(file_path, Path):
            raise TypeError("file_path must be a pathlib.Path instance")
        
        self.file_path = file_path
        self.secondary_structure_path = secondary_structure_path
        self.logger = logger
        self.verbose = verbose
        self.log_messages = []
        self.pdb_struct = OPMPDBStruct(
            logger=logger,
            min_length = min_length,
            max_length = max_length,
            margin = margin, 
            inner_margin = inner_margin,
            min_segment_length = min_segment_length, 
            gaps = gaps
        )
        self.mb_coords = []


    def parse(self) -> OPMPDBStruct:
        """
        Parse the PDB file (and the secondary structure file, if provided) 
        and return a PDBStruct object with the parsed data.
        """
        self._initialize_structure()
        self._parse_pdb_file()
        self._compute_protein_length()
        if self.secondary_structure_path:
            self._parse_secondary_structure()
        self._finalize()
        if self.verbose and self.logger:
            self.logger.info("\n".join(self.log_messages))
        return self.pdb_struct

    def _initialize_structure(self):
        self.pdb_struct.protein_name = self.file_path.stem
        self.pdb_struct.full = defaultdict(lambda: SlicedOrderedDict())
        self.pdb_struct.CA = defaultdict(lambda: SlicedOrderedDict())
        self.pdb_struct.sequences = defaultdict(str)
        self.log_messages.append(f"Initialized parsing for {self.file_path.stem}")
        
    def _parse_pdb_file(self):
        with open(self.file_path, "r") as f:
            for line in f:
                suspected_membrane_line = False
                if line.startswith("ATOM"):
                    try:
                        atom_number = int(line[6:11].strip())
                        atom_name = line[12:16].strip()
                        res_number = int(line[22:26].strip())
                        res_name = line[17:20].strip()
                        chain_id = line[21].strip() or "A"
                        x, y, z = map(float, (
                            line[30:38].strip(),
                            line[38:46].strip(),
                            line[46:54].strip()
                        ))
                    except Exception as e:
                        self.log_messages.append(f"Error parsing line: {line.strip()} ({e})")
                        continue
                    
                    if res_name == "DUM":
                        suspected_membrane_line = True
                        self.log_messages.append(f"Membrane line found: {line.strip()}")
                    
                    if not suspected_membrane_line:
                        if chain_id not in self.pdb_struct.full:
                            self.pdb_struct.full[chain_id] = SlicedOrderedDict()
                            self.pdb_struct.CA[chain_id] = SlicedOrderedDict()
                            self.log_messages.append(f"Chain {chain_id} found")
                        
                        if res_number not in self.pdb_struct.full[chain_id]:
                            self.pdb_struct.full[chain_id][res_number] = SlicedOrderedDict()
                        self.pdb_struct.full[chain_id][res_number][atom_number] = line
                        
                        if atom_name == "CA":
                            self.pdb_struct.CA[chain_id][res_number] = {
                                "coord": [x, y, z],
                                "res_name": self.aa_dict.get(res_name, "X"),
                                "res_number": res_number,
                            }
                            self.pdb_struct.sequences[chain_id] += self.aa_dict.get(res_name, "X")
                elif line.startswith("HETATM") and "DUM" in line:
                    try:
                        x, y, z = map(float, (
                            line[30:38].strip(),
                            line[38:46].strip(),
                            line[46:54].strip()
                        ))
                        self.mb_coords.append([x, y, z])
                    except Exception as e:
                        self.log_messages.append(f"Error parsing HETATM line: {line.strip()} ({e})")
                        continue
        
        self.pdb_struct.membrane_coord = np.array(self.mb_coords)
        self.log_messages.append("Finished parsing PDB atoms.")

    def _compute_protein_length(self):
        self.pdb_struct.protein_length = {}
        for chain_id, chain in self.pdb_struct.CA.items():
            self.pdb_struct.protein_length[chain_id] = len(chain)
        self.log_messages.append("Computed protein lengths for all chains.")

    def _parse_secondary_structure(self):
        ss_dict = defaultdict(dict)
        self.log_messages.append("Secondary structure found")
        # Format of the file should be parseable with .split() and
        # have chain_id / res_num / secondary_structure
        with open(self.secondary_structure_path, "r") as f:
            for line in f:
                parts = line.split()
                if len(parts) < 3:
                    continue
                chain_id = parts[0]
                try:
                    res_number = int(parts[1])
                except ValueError:
                    self.log_messages.append(f"Invalid residue number in line: {line.strip()}")
                    continue
                secondary_structure = parts[2]
                self.log_messages.append(f"Ss: {secondary_structure}, Chain: {chain_id}, Res: {res_number}")
                ss_dict[chain_id][res_number] = secondary_structure

        # Update CA records with secondary structure informations 
        for chain_id, residues in self.pdb_struct.CA.items():
            for res_number, data in residues.items():
                sec_struct = ss_dict.get(chain_id, {}).get(res_number, "O")
                data.update({"secondary_structure": sec_struct})
                data.update({"folded": sec_struct in ("H", "S")})
                if sec_struct in ("H", "S"):
                    self.log_messages.append(f"Res {res_number} in chain {chain_id} is in a folded region")
        self.log_messages.append("Secondary structure parsing complete.")

    def _finalize(self):
        self.log_messages.append(f"{self.pdb_struct.protein_name} processed")
        self.log_messages.append(f"Chains found: {list(self.pdb_struct.CA.keys())}")
        self.log_messages.append(f"Pdb structure: {self.pdb_struct}")

def write_fasta_records(seq_records, filename):
    """
    Write a list of SeqRecord objects to a FASTA file.
    """
    with open(filename, "w") as f:
        SeqIO.write(seq_records, f, "fasta")


def write_structures_to_pdb(structures, output_dir):
    """
    Write structure dictionaries to individual PDB files in the given directory.
    
    Parameters:
      structures (dict): Dictionary mapping chain IDs to dictionaries mapping a unique key (e.g. a fasta key)
                         to a residue-level structure dictionary.
                         The residue-level structure dictionary should have residue numbers as keys,
                         and each value should be a dictionary mapping atom numbers to the atom record (a string).
      output_dir (str): Directory to write the PDB files into.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # Iterate over each chain.
    for chain_id, chain_struct in structures.items():
        for key, residue_dict in chain_struct.items():
            pdb_lines = []
            # Sort residues by number
            for res in sorted(residue_dict.keys()):
                value = residue_dict[res]
                # If the value is a dictionary, then we assume it maps atom numbers to atom lines.
                if isinstance(value, dict):
                    # Sort by atom number so that the atom lines are in order.
                    for atom_number in sorted(value.keys()):
                        pdb_lines.append(value[atom_number])
                else:
                    # If it's not a dictionary (e.g. a placeholder like "X"), we skip it.
                    continue
            # Write the concatenated atom lines to a PDB file named after the key.
            output_file = os.path.join(output_dir, f"{key}.pdb")
            with open(output_file, "w") as f:
                f.write("".join(pdb_lines))


def write_output_files(output_dict, fasta_output_dir, pdb_output_dir):
    """
    Write the output dictionary to FASTA and PDB files.
    
    Parameters:
      output_dict (dict): Dictionary with keys:
                          - "sequences": list of full elongated SeqRecord objects,
                          - "sequences_shorts": list of core (short) SeqRecord objects,
                          - "structures": dict mapping chain IDs to structure dicts for full sequences,
                          - "structures_shorts": dict mapping chain IDs to structure dicts for core sequences.
      fasta_output_dir (str): Directory where FASTA files will be written.
      pdb_output_dir (str): Directory where PDB files will be written.
    """
    # Create output directories if they don't exist.
    if not os.path.exists(fasta_output_dir):
        os.makedirs(fasta_output_dir)
    if not os.path.exists(pdb_output_dir):
        os.makedirs(pdb_output_dir)
    
    # Write FASTA files.
    full_fasta_file = os.path.join(fasta_output_dir, "full_sequences.fasta")
    short_fasta_file = os.path.join(fasta_output_dir, "short_sequences.fasta")
    write_fasta_records(output_dict["sequences"], full_fasta_file)
    write_fasta_records(output_dict["sequences_shorts"], short_fasta_file)
    
    # Write PDB files into separate subdirectories for full and short structures.
    full_pdb_dir = os.path.join(pdb_output_dir, "full_structures")
    short_pdb_dir = os.path.join(pdb_output_dir, "short_structures")
    write_structures_to_pdb(output_dict["structures"], full_pdb_dir)
    write_structures_to_pdb(output_dict["structures_shorts"], short_pdb_dir)


def process_pdb_file(pdb_path):
    """
    Process a single PDB file and return the output dictionary.
    """
    try:
        parser_obj = OPMPDBParser(pdb_path, min_length=60, max_length=100, margin=10)
        OPMStruct = parser_obj.parse()
        data = OPMStruct.analyze()
        return data
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return None

def merge_output(master_output, data):
    """
    Merge a single data dictionary into the master_output.
    """
    if data is None:
        return
    master_output["sequences"].extend(data["sequences"])
    master_output["sequences_shorts"].extend(data["sequences_shorts"])

    for chain_id, struct_dict in data["structures"].items():
        if chain_id not in master_output["structures"]:
            master_output["structures"][chain_id] = {}
        master_output["structures"][chain_id].update(struct_dict)

    for chain_id, struct_dict in data["structures_shorts"].items():
        if chain_id not in master_output["structures_shorts"]:
            master_output["structures_shorts"][chain_id] = {}
        master_output["structures_shorts"][chain_id].update(struct_dict)



def search_main_directory():

    """
    Iteratively searches for the peptide directory and returns its absolute path.
    """

    global main_directory
    main_directory = None
    for i in range(3): 
        backward = "../" * i
        main_directory = Path(f"{backward}Peptides").resolve()
        if main_directory.exists():
            break

    if main_directory is None:
        raise FileNotFoundError("Peptide directory not found")
    
    print(f"Working on the main directory : {main_directory}")

def get_paths():

    """
    ## type_id correspondance : 

    #   1 : Multitopic transmembrane
    #   2 : Monotopic / Peripherals
    #   3 : Peptides 

    ## classtype_id correspondance :

    ###### Multitopic transmembrane ######
    #  1  : Alpha-helical polytopic 
    #  11 : Bitopic proteins 
    #  2  : Beta barrel 

    ###### Monotopic ######
    #  4  : All alpha 
    #  3  : All beta
    #  5  : Alpha / beta
    #  6  : Alpha + beta

    ###### Peptides #######
    #  7  : Alpha-helical peptides
    #  8  : Beta-helical peptides
    #  9  : Beta hairpins 
    #  10 : Non-regular peptides



    Returns paths for each category of peptides/proteins that are present in the OPM and Membranome databases.

    """

    metadata = (

        pl.read_csv("/Users/simonherman/Documents/I2BC/Peptides/input/proteins-2024-05-07.csv", separator = ",", infer_schema_length = 20000)
            .with_columns(
                pl.concat_str([
                        pl.lit("/Users/simonherman/Documents/I2BC/Peptides/input/OPM"),
                        pl.concat_str([pl.col("pdbid"), pl.lit(".pdb")], separator = "")           
                    ], separator = "/",
                ).alias("pdb_path")
            )
    )

    ## FILTERS ## 

    # Transmembranes

    # Annotated bitopic proteins, put threshold of 20 to avoid some false positives
    bitopic_proteins = ((pl.col("classtype_id") == 11) & (pl.col("thickness") >= 20)) 
    # All peptides that are crossing the membrane ( > 20 ), regardless of their folding type
    bitopic_peptides = ((pl.col("type_id") == 3) & (pl.col("thickness") >= 20)) 
    # Multitopic proteins to be cut 
    polytopic_proteins = (pl.col("classtype_id") == 1) & (pl.col("thickness") >= 20)

    ## CREATING PATHS ##

    membranome = [ main_directory / path for path in list(glob("/Users/simonherman/Documents/I2BC/Peptides/input/Membranome/*.pdb")) ]
    bitopic_proteins = [ main_directory / path for path in metadata.filter(bitopic_proteins)["pdb_path"].to_list() ]
    bitopic_peptides = [ main_directory / path for path in metadata.filter(bitopic_peptides)["pdb_path"].to_list() ]
    polytopics = [ main_directory / path for path in metadata.filter(polytopic_proteins)["pdb_path"].to_list() ]

    tm_paths = {
        
        "bitopic" : membranome + bitopic_proteins + bitopic_peptides,
        "polytopic" : polytopics

    }

    return tm_paths

def main():

    master_output = {
        "sequences": [],
        "sequences_shorts": [],
        "structures": {},
        "structures_shorts": {}
    }

    search_main_directory()

    paths = get_paths()
    bitopics = paths["bitopic"]
    polytopics = paths["polytopic"]

    pdbs_list = bitopics + polytopics
    total_files = len(pdbs_list)


    with concurrent.futures.ProcessPoolExecutor(max_workers=10, max_tasks_per_child=1) as executor:

        results = list(tqdm(executor.map(process_pdb_file, pdbs_list), total=total_files, desc="Processing PDBs"))

    for data in results:
        merge_output(master_output, data)

    fasta_output_dir = "all_fasta_outputs"
    pdb_output_dir = "all_pdb_outputs"
    write_output_files(master_output, fasta_output_dir, pdb_output_dir)

if __name__ == "__main__":
    main()


    