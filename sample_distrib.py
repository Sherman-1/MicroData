#!/usr/bin/env python3

import sys
import random
from collections import defaultdict
from Bio import SeqIO

import matplotlib.pyplot as plt
import seaborn as sns

import polars as pl 

###############################################################################
DESIRED_DISTRIBUTION = {
    (0, 24): 0.25,
    (25, 49): 0.25,
    (50, 74): 0.25,
    (75, 100): 0.25  
}

desired_total = 10000
random_seed = 66


def length_to_bin(seq_len, bins):
    """
    Given a sequence length and a dictionary of bins (range tuples) -> probability,
    find which bin the length belongs to. If no bin is found, return None.
    """
    for (low, high) in bins:
        if low <= seq_len <= high:
            return (low, high)
    return None

def check_distribution_probs(dist_dict):
    """
    Check that probabilities sum to something near 1.0 (or exactly).
    Raise a warning if not.
    """
    s = sum(dist_dict.values())
    if abs(s - 1.0) > 1e-6:
        print(f"WARNING: distribution sums to {s:.4f}, not 1.0!", file=sys.stderr)

def sample_fasta_custom_dist(input_fasta, output_fasta):
    
    random.seed(random_seed)
    check_distribution_probs(DESIRED_DISTRIBUTION)
    
    bin_dict = defaultdict(list)
    
    with open(input_fasta, "r") as f_in:
        for record in SeqIO.parse(f_in, "fasta"):
            seq_len = len(record.seq)
            b = length_to_bin(seq_len, DESIRED_DISTRIBUTION)
            if b is not None:
                bin_dict[b].append(record)

    print(bin_dict.keys())
    
    sampled_records = []
    
    for b, probability in DESIRED_DISTRIBUTION.items():

        target_count = int(round(probability * desired_total))
        
        bin_seqs = bin_dict.get(b, [])
        print(f"Bin {b}: {len(bin_seqs)} sequences, target count: {target_count}, target probability: {probability}")
        
        # If more or same number of sequences in the bin as the target count, sample without replacement
        if len(bin_seqs) <= target_count:
            sampled_records.extend(bin_seqs)
        else:
            sampled_records.extend(random.sample(bin_seqs, target_count))
    
    with open(output_fasta, "w") as f_out:
        SeqIO.write(sampled_records, f_out, "fasta")

    return sampled_records

def plot_output_fasta_distribution(input_records, sampled_records) -> None:

    resampled_lengths = [len(record.seq) for record in sampled_records]
    input_lengths = [len(record.seq) for record in input_records]

    data = pl.DataFrame(
        {
            "length": input_lengths + resampled_lengths,
            "dataset": ["input"] * len(input_lengths) + ["resampled"] * len(resampled_lengths)
        }
    )

    plt.figure(figsize=(10, 6))
    
    sns.violinplot(
        data = data,
        y = "length",
        split = True,
        hue = "dataset",
        density_norm = "count",
        inner="quart", 
        fill=False,
        palette={"input": "g", "resampled": ".35"}
    )

    plt.title("Sequence length distribution")
    plt.savefig("resampled_vs_original_length_distribution.png")
    plt.close()

    return 0

if __name__ == "__main__":
    """
    Usage:
        python sample_fasta_custom_dist.py input.fasta output.fasta
    """

    import argparse
    
    parser = argparse.ArgumentParser(description="Sample sequences from a FASTA file according to a custom length distribution.")
    parser.add_argument("--input_fasta", help="Input FASTA file")
    parser.add_argument("--output_fasta", help="Output FASTA file")

    args = parser.parse_args()
    
    input_fasta = args.input_fasta
    output_fasta = args.output_fasta
    resample_records = sample_fasta_custom_dist(input_fasta, output_fasta)
    plot_output_fasta_distribution(list(SeqIO.parse(input_fasta, "fasta")), resample_records)
