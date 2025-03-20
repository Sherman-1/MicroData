import requests
import polars as pl 
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from io import StringIO

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from MMseqs import MMseqs2

import random

def print_banner(msg):

    if msg == "disprot":

        logo = r"""
            .----------------.  .----------------.  .----------------.  .----------------.  .----------------.  .----------------.  .----------------. 
            | .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |
            | |  ________    | || |     _____    | || |    _______   | || |   ______     | || |  _______     | || |     ____     | || |  _________   | |
            | | |_   ___ `.  | || |    |_   _|   | || |   /  ___  |  | || |  |_   __ \   | || | |_   __ \    | || |   .'    `.   | || | |  _   _  |  | |
            | |   | |   `. \ | || |      | |     | || |  |  (__ \_|  | || |    | |__) |  | || |   | |__) |   | || |  /  .--.  \  | || | |_/ | | \_|  | |
            | |   | |    | | | || |      | |     | || |   '.___`-.   | || |    |  ___/   | || |   |  __ /    | || |  | |    | |  | || |     | |      | |
            | |  _| |___.' / | || |     _| |_    | || |  |`\____) |  | || |   _| |_      | || |  _| |  \ \_  | || |  \  `--'  /  | || |    _| |_     | |
            | | |________.'  | || |    |_____|   | || |  |_______.'  | || |  |_____|     | || | |____| |___| | || |   `.____.'   | || |   |_____|    | |
            | |              | || |              | || |              | || |              | || |              | || |              | || |              | |
            | '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |
            '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------' 
        """

    elif msg == "mobid":

        logo = r"""
            .----------------.  .----------------.  .----------------.  .----------------.  .----------------.  .----------------. 
            | .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |
            | | ____    ____ | || |     ____     | || |   ______     | || |     _____    | || |  ________    | || |   ______     | |
            | ||_   \  /   _|| || |   .'    `.   | || |  |_   _ \    | || |    |_   _|   | || | |_   ___ `.  | || |  |_   _ \    | |
            | |  |   \/   |  | || |  /  .--.  \  | || |    | |_) |   | || |      | |     | || |   | |   `. \ | || |    | |_) |   | |
            | |  | |\  /| |  | || |  | |    | |  | || |    |  __'.   | || |      | |     | || |   | |    | | | || |    |  __'.   | |
            | | _| |_\/_| |_ | || |  \  `--'  /  | || |   _| |__) |  | || |     _| |_    | || |  _| |___.' / | || |   _| |__) |  | |
            | ||_____||_____|| || |   `.____.'   | || |  |_______/   | || |    |_____|   | || | |________.'  | || |  |_______/   | |
            | |              | || |              | || |              | || |              | || |              | || |              | |
            | '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |
            '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------' 
        """

    print(logo)

def merge_intervals_with_type(intervals):
    """
    Merge a list of intervals (each as a tuple (start, end, prediction_type)) that overlap.
    If overlapping intervals have different prediction types, keep the one of highest fidelity.
    Returns a list of merged intervals with their associated prediction type.
    """
    if not intervals:
        return []
    
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    
    priority = {"curated": 3, "derived": 2, "homology": 1, "predicted" : 0} # Curated is best 
    
    merged = [sorted_intervals[0]]
    
    for current in sorted_intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            new_start = last[0]
            new_end = max(last[1], current[1])
            best_type = last[2] if priority[last[2]] >= priority[current[2]] else current[2]
            merged[-1] = (new_start, new_end, best_type)
        else:
            merged.append(current)
    
    return merged

def merge_intervals(intervals):
    """
    Merge a list of intervals (each as a tuple (start, end)) that overlap.
    Returns a list of merged intervals.
    """
    if not intervals:
        return []
    
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_intervals[0]]
    
    for current in sorted_intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    
    return merged

def random_subregion_slices(start, stop, min_length, max_length):
    """
    Takes a region [start, stop] (1-based indices) and returns a list
    of random sub-regions, each of length in [min_length, max_length].
    """
    slices = []
    region_len = stop - start + 1
    
    current_start = start
    
    while True:
        leftover = stop - current_start + 1
        
        if leftover < min_length:

            if slices and ( (slices[-1][1] + leftover) - slices[-1][0] + 1 ) <= max_length:
 
                old_start, old_stop = slices[-1]
                new_stop = old_stop + leftover
                slices[-1] = (old_start, new_stop)

            break
        

        if leftover <= max_length:
            slices.append((current_start, stop))
            break
        
        max_chunk_possible = leftover - min_length

        chunk_size = random.randint(min_length, min(max_chunk_possible, max_length))
        
        chunk_end = current_start + chunk_size - 1
        slices.append((current_start, chunk_end))
        
        current_start = chunk_end + 1

    return slices


def add_subregion(
    results, acc, pred_type, full_sequence,
    sub_region, merged_lip, merged_extended, merged_compact,
    merged_dto, merged_dtd, merged_context, merged_phase_separation
):
    """
    Helper to handle overlap calculations and append data to 'results'.
    """
    start, stop = sub_region
    region_length = stop - start + 1
    sequence = full_sequence[start - 1:stop]  # 1-based indexing
    header = f"{acc}|{pred_type}|{start}-{stop}"
    
    results["header"].append(header)
    results["sequence"].append(sequence)
    results["length"].append(region_length)
    results["prediction_type"].append(pred_type)

    lip_overlap = sum(overlap((start, stop), (lip_start, lip_stop)) 
                      for lip_start, lip_stop in merged_lip)
    results["lip_overlap"].append(float(lip_overlap))

    extended_overlap = sum(overlap((start, stop), (ext_start, ext_stop)) 
                           for ext_start, ext_stop in merged_extended)
    results["extended_overlap"].append(float(extended_overlap))

    compact_overlap = sum(overlap((start, stop), (comp_start, comp_stop)) 
                          for comp_start, comp_stop in merged_compact)
    results["compact_overlap"].append(float(compact_overlap))

    dto_overlap = sum(overlap((start, stop), (dto_start, dto_stop)) 
                      for dto_start, dto_stop in merged_dto)
    results["dto_overlap"].append(float(dto_overlap))

    dtd_overlap = sum(overlap((start, stop), (dtd_start, dtd_stop)) 
                      for dtd_start, dtd_stop in merged_dtd)
    results["dtd_overlap"].append(float(dtd_overlap))

    context_overlap = sum(overlap((start, stop), (ctx_start, ctx_stop)) 
                          for ctx_start, ctx_stop in merged_context)
    results["context_overlap"].append(float(context_overlap))

    phase_overlap = sum(overlap((start, stop), (ps_start, ps_stop))
                        for ps_start, ps_stop in merged_phase_separation)
    results["phase_separation_overlap"].append(float(phase_overlap))

    # neutral overlap is simply what's left from the region 
    # after subtracting extended and compact
    results["neutral_overlap"].append(1.0 - extended_overlap - compact_overlap)


def overlap(interval_a, interval_b):
    """
    Example overlap function, returns fraction of interval_a
    that overlaps with interval_b.
    """
    start_a, end_a = interval_a
    start_b, end_b = interval_b
    overlap_start = max(start_a, start_b)
    overlap_end = min(end_a, end_b)
    
    if overlap_start > overlap_end:
        return 0.0
    
    overlap_len = (overlap_end - overlap_start + 1)
    length_a = (end_a - start_a + 1)
    return overlap_len / length_a

def extract_disordered_regions(entry, min_length=20, max_length=100):
    """
    Extracts non-overlapping disordered regions (after merging overlaps)
    whose length is between min_length and max_length (inclusive) from the protein sequence.
    
    For each region, the header includes the prediction type (curated, derived, homology)
    with a fidelity ranking: curated > derived > homology.
    
    Returns:
        A tuple (results, flag) where results is a dictionary containing region data,
        and flag is True if the protein length is less than 100, False otherwise.
    """
    
    try:
        acc = entry["acc"]
        full_sequence = entry["sequence"]
        protein_length = entry["length"]
    except KeyError:
        print("Entry is missing required fields.")
        return None, None
    
    disorder_intervals = []  #  (start, end, prediction_type)
    lip_intervals = []
    extended_intervals = []
    compact_intervals = []
    disorder_to_order_intervals = []
    disorder_to_disorder_intervals = []
    context_dependent_intervals = []
    phase_separation_intervals = []
    
    for key, value in entry.items():
        if not isinstance(value, dict):
            continue
        
        regions = value.get("regions", [])
        
        disord_match = re.match(r"^(curated|homology|derive|predicted)-disorder-(merge|mobidb_lite)", key)
        if disord_match:
            prediction_type = disord_match.group(1)
            for region in regions:
                disorder_intervals.append((region[0], region[1], prediction_type))
        elif re.match(r"^(?:curated|homology|derived)-lip-merge$", key):
            lip_intervals.extend(regions)
        elif re.match(r"^(?:curated|homology|derived)-phase_separation-merge$", key):
            phase_separation_intervals.extend(regions)
        elif key == "prediction-extended-mobidb_lite_sub":
            extended_intervals.extend(regions)
        elif key == "prediction-compact-mobidb_lite_sub":
            compact_intervals.extend(regions)
        elif re.match(r"^(?:curated|homology|derived)-binding_mode_disorder_to_order-merge$", key):
            disorder_to_order_intervals.extend(regions)
        elif re.match(r"^(?:curated|homology|derived)-binding_mode_disorder_to_disorder-merge$", key):
            disorder_to_disorder_intervals.extend(regions)
        elif re.match(r'^(?:curated|homology|derived)-binding_mode_context_dependent-.*$', key):
            context_dependent_intervals.extend(regions)
    
    merged_disorder = merge_intervals_with_type(disorder_intervals) if disorder_intervals else []
    merged_lip = merge_intervals(lip_intervals) if lip_intervals else []
    merged_extended = merge_intervals(extended_intervals) if extended_intervals else []
    merged_compact = merge_intervals(compact_intervals) if compact_intervals else []
    merged_dto = merge_intervals(disorder_to_order_intervals) if disorder_to_order_intervals else []
    merged_dtd = merge_intervals(disorder_to_disorder_intervals) if disorder_to_disorder_intervals else []
    merged_context = merge_intervals(context_dependent_intervals) if context_dependent_intervals else []
    merged_phase_separation = merge_intervals(phase_separation_intervals) if phase_separation_intervals else []
    
    results = {
        "header": [],
        "sequence": [],
        "length" : [],
        "prediction_type": [],
        "lip_overlap": [],
        "extended_overlap": [],
        "compact_overlap": [],
        "neutral_overlap": [],
        "dto_overlap": [],
        "dtd_overlap": [],
        "context_overlap": [],
        "phase_separation_overlap": []  
    }
    
    for idr_start, idr_stop, pred_type in merged_disorder:
        region_length = idr_stop - idr_start + 1
        
        if region_length < min_length:
            continue
        
        if region_length <= max_length:
            add_subregion(
                results, acc, pred_type, full_sequence,
                (idr_start, idr_stop),
                merged_lip, merged_extended, merged_compact,
                merged_dto, merged_dtd, merged_context, merged_phase_separation
            )
        else:

            slices = random_subregion_slices(idr_start, idr_stop, min_length, max_length)
            for (sub_start, sub_stop) in slices:
                add_subregion(
                    results, acc, pred_type, full_sequence,
                    (sub_start, sub_stop),
                    merged_lip, merged_extended, merged_compact,
                    merged_dto, merged_dtd, merged_context, merged_phase_separation
                )

    # Check that extended + compact + neutral = 1.0
    for i in range(len(results["header"])):
        total = results["extended_overlap"][i] + results["compact_overlap"][i] + results["neutral_overlap"][i]
        if abs(total - 1.0) > 1e-6:
            print(f"Warning: extended + compact + neutral != 1.0 for acc {acc} in region {results['header'][i]}")

    flag = True if protein_length < 100 else False
    return results, flag


def write_fasta_append(region_dict, output_file):
    """
    Appends disordered regions in FASTA format to an existing file.
    Each key in `region_dict` is used as the FASTA header.
    Each value is the protein subsequence.
    """
    with open(output_file, "a") as f:
        for header, sequence in region_dict.items():
            f.write(f">{header}\n{sequence}\n")

def fetch_page(last_id=None):
    """
    Fetches a single page from the MobiDB API using pagination.
    If last_id is provided, it will be used in the '_id_gt' parameter.
    """
    base_url = "https://mobidb.org/api/download_page"

    params = {
        "sort_asc": "_id",
        "limit": 1000,
        "format": "json",
        "level": 2
    }

    if last_id:
        params["_id_gt"] = last_id

    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Request failed with status code {response.status_code}")
        print(response.text)
        return None

def update_dict_of_lists(target : dict, source : dict):
    """
    Updates the target dict of lists by appending the values from the source dict.
    If a key does not exist in target, it is created.
    """

    try:
        for key, values in source.items():
            if key not in target:
                target[key] = list(values)  
            else:
                target[key].extend(values)

    except Exception as e:
        print("Error updating dict of lists:", e)

         
def download_disprot_fasta(link):

    """
    Download the DisProt FASTA file and save it to the output_file.
    """
    
    response = requests.get(link)

    if response.status_code == 200:
        fasta_content = response.text
        fasta_io = StringIO(fasta_content)
        fasta_iterable = SeqIO.parse(fasta_io, "fasta")
    else:
        print("Failed to retrieve the file.")

    return fasta_iterable


def cut_fasta_random_chunks(fasta_iterable, min_length=20, max_length=100):
    """
    Split sequences into random-sized chunks (each between min_length and max_length).
    
    Rules:
      1. If length < min_length, skip.
      2. If min_length <= length <= max_length, keep whole.
      3. If length > max_length, cut into random-length chunks in [min_length, max_length].
         - Each chunk is at least min_length in length.
         - If the leftover is shorter than min_length, discard it.
    """
    
    results = []
    for record in fasta_iterable:
        full_len = len(record.seq)
        
        if full_len < min_length:
            continue
        
        if full_len <= max_length:
            results.append(
                SeqRecord(
                    record.seq,
                    id=f"{str(record.id).split('|')[1]}|curated|1-{full_len}",
                    description=""
                )
            )
            continue
        
        chunk_start = 0
        while chunk_start < full_len:
            leftover = full_len - chunk_start
            if leftover < min_length:
                break
            
            chunk_size = random.randint(min_length, min(max_length, leftover))
            chunk_end = chunk_start + chunk_size
            
            sub_seq = record.seq[chunk_start:chunk_end]
            new_record = SeqRecord(
                Seq(str(sub_seq)),
                id=f"{str(record.id).split('|')[1]}|curated|{chunk_start+1}-{chunk_end}",
                description=""
            )
            results.append(new_record)
            
            chunk_start = chunk_end
            
    return results

def main():

    proteins_output = "mobidb_disordered_regions.csv"
    microproteins_output = "mobidb_disordered_microproteins.csv"
    mmseqs_path = "/home/simon.herman/.local/bin/mmseqs"
    disprot_url = "https://disprot.org/api/search?release=2024_06&show_ambiguous=false&show_obsolete=false&format=fasta&namespace=all&get_consensus=false"

    mmseqs2_api = MMseqs2(threads=6, cleanup=True, mmseqs2_path=mmseqs_path)
    
    for output in [proteins_output, microproteins_output]:
        if os.path.exists(output):
            os.remove(output)

    ###### MOBIDB ######
    print_banner("mobid")
    
    proteins = {}
    microproteins = {}
    last_id = None 
    page_count = 0

    while True:

        data_json = fetch_page(last_id=last_id)
        print("Page fetched!")
        
        if not data_json:
            print("No data present in the fetched page.")
            print("Metadata:", data_json.get("metadata"))
            break
        
        entries = data_json.get("data", [])
        if not entries:
            print("No entries found.")
            print("Metadata:", data_json.get("metadata"))
            break
        
        for entry in entries:
            disordered, micro_flag = extract_disordered_regions(entry)
            if micro_flag:
                update_dict_of_lists(microproteins, disordered)
            else:
                update_dict_of_lists(proteins, disordered)
    
        metadata = data_json.get("metadata", {})
        last_id = metadata.get("last_id")
        if not last_id:
            print("No more pages to fetch.")
            break

        page_count += 1

        print(f"Found {len(proteins.get('header', []))} proteins and {len(microproteins.get('header', []))} microproteins so far.")

    print(f"Fetched {page_count} pages.")
    print(f"Found {len(proteins.get('header', []))} proteins and {len(microproteins.get('header', []))} microproteins.")

    # Turn dicts into fasta writeable format
    proteins_sequences = [ SeqRecord(Seq(sequence), id=header, description="", name="") for header, sequence in zip(proteins["header"], proteins["sequence"]) ]
    microproteins_sequences = [ SeqRecord(Seq(sequence), id=header, description="", name="") for header, sequence in zip(microproteins["header"], microproteins["sequence"]) ]


    ###### DISPROT ######
    print_banner("disprot")

    disprot_iterable = download_disprot_fasta(disprot_url)
    records = cut_fasta_random_chunks(disprot_iterable)


    full_sequences = proteins_sequences + microproteins_sequences + records

    SeqIO.write(full_sequences, "disordered_sequences.fasta", "fasta")

    representatives = mmseqs2_api.fasta2representative_mobidb("disordered_sequences.fasta", cov = 0.7, iden = 0.3)

    SeqIO.write(representatives.values(), "representative_disordered_sequences.fasta", "fasta")

if __name__ == "__main__":

    main()