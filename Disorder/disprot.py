#!/usr/bin/env python3

from pathlib import Path 
import shutil
import requests
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
import random
from MMseqs import MMseqs2

MIN_LENGTH = 20
MAX_LENGTH = 100

            
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


def cut_fasta_mobidb_style(fasta_iterable, min_length=20, max_length=100):
    """
    Replicates the 'MobiDB-style' cutting approach on an iterable of SeqRecords.

    For each sequence:
      1. If length < min_length, skip it.
      2. If min_length <= length <= max_length, keep it whole.
      3. If length > max_length, split into consecutive chunks of size up to max_length.
         - Each chunk must be at least min_length in length, or it is discarded.
    """

    results = []
    for record in fasta_iterable:
        full_len = len(record.seq)
        
        if full_len < min_length:
            continue
        
        if full_len <= max_length:
            results.append(record)
            continue
        
        chunk_start = 0
        while chunk_start < full_len:
            chunk_end = chunk_start + max_length
            if chunk_end > full_len:
                chunk_end = full_len  
            chunk_len = chunk_end - chunk_start
            
            if chunk_len >= min_length:
                sub_seq = record.seq[chunk_start:chunk_end]
                new_record = SeqRecord(
                    Seq(str(sub_seq)),
                    id=f"{record.id}_{chunk_start+1}-{chunk_end}",
                    description=""
                )
                results.append(new_record)
            
            chunk_start = chunk_end
        
    return results

def main():

    disprot_link = "https://disprot.org/api/search?release=2024_06&show_ambiguous=false&show_obsolete=false&format=fasta&namespace=all&get_consensus=false"

    disprot_iterable = download_disprot_fasta(disprot_link)

    records = cut_fasta_mobidb_style(disprot_iterable, MIN_LENGTH, MAX_LENGTH)

    print(f"Extracted {len(records)} sequences.")

    SeqIO.write(records, "disprot.fasta", "fasta")


if __name__ == "__main__":

    main()