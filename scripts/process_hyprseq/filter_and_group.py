"""
Run the script with the following command:
python /broad/thechenlab/stag_seq/hypr_pipeline/filter_and_group.py \
    -p /broad/thechenlab/stag_seq/matrix_screen_v2/hypr_samples/donor2/Probe_sets.txt \
    -f /broad/thechenlab/stag_seq/test/Donor1_Ctrl_HyPR_S1_R1_001_work_dir/0_output.txt \
    -b /broad/thechenlab/stag_seq/hypr_pipeline/barcode_lookup.pkl

Or use extracted.txt as output. 
"""

import numpy as np
import pandas as pd
import pickle
import sys
import argparse


# create the mismatch table from the probe set
def one_mismatch_sequences(sequence):
    """Generate all sequences that are one mismatch away from the input sequence."""
    nucleotides = ['A', 'C', 'T', 'G']
    mismatch_seqs = []
    for i in range(len(sequence)):
        for nuc in nucleotides:
            if sequence[i] != nuc:
                mismatch_seq = sequence[:i] + nuc + sequence[i+1:]
                mismatch_seqs.append(mismatch_seq)
                
    return mismatch_seqs

def parse_args():
    parser = argparse.ArgumentParser(description="Filter and group reads")
    parser.add_argument(
        "-p", "--probe_file", type=str, required=True, help="Probe set file",
    )
    parser.add_argument(
        "-f", "--in_file", type=str, required=True, help="Concatenated output file",
    )
    parser.add_argument(
        "-b", "--barcode_hash", type=str, required=True, help="Barcode hash table",
    )
    parser.add_argument(
        "-s", "--save_name_index", type=str, required=False, help="Save name index", default="test",
    )
    args = parser.parse_args()
    return args

def create_probe_lookup(probe_file):
    """Create a lookup table for probes and their mismatches."""
    # Read the file of correct sequences into a dictionary
    correct_sequences = {}
    with open(probe_file, 'r') as f:
        for line in f:
            try:
                seq_id, seq = line.strip().split('\t')
                correct_sequences[seq] = seq_id
            except:
                print(line)
    print(str(len(correct_sequences.keys())) + " probe sequences identified")

    # Create probe lookup table
    lookup_table = {}
    for seq, seq_id in correct_sequences.items():
        lookup_table[seq] = seq_id
        for mismatch_seq in one_mismatch_sequences(seq):
            lookup_table[mismatch_seq] = seq_id
    return lookup_table

if __name__ == "__main__":
    # Parse arguments.
    args = parse_args()
    probe_file = args.probe_file
    in_file = args.in_file
    barcode_hash = args.barcode_hash
    save_name_index = args.save_name_index
    
    print("probe_file: " + probe_file)
    print("in_file: " + in_file)
    print("barcode_hash: " + barcode_hash)
    

    # load in the data
    sys.stdout.write("reading concatenated output into dataframe\n")
    df = pd.read_csv(in_file, sep=' ', names=['barcode', 'umi', 'probe'])

    # Probe lookup table
    sys.stdout.write("creating probe lookup table\n")
    lookup_table = create_probe_lookup(probe_file)

    # Use the probe lookup to replace mismatches and discard anything with more than 1 mismatch
    sys.stdout.write("replacing probe mismatches\n")
    sys.stdout.write(str(len(df)) + " reads before filtering\n")
    df['probe'] = df['probe'].map(lookup_table)
    df.dropna(subset=['probe'], inplace=True)
    sys.stdout.write(str(len(df)) + " reads have valid probes\n")

    # Drop duplicate from dataframe
    # Print number of rows before and after
    sys.stdout.write("dropping duplicates\n")
    sys.stdout.write(str(len(df)) + " reads before dropping duplicates\n")
    df = df.drop_duplicates(subset=['probe', 'umi', 'barcode'])
    sys.stdout.write(str(len(df)) + " reads after dropping duplicates\n")
    # grouped_barcode = df.groupby(['barcode']).size().reset_index(name='count')
    
    # GROUPBY
    grouped_df = df.groupby(['probe', 'barcode']).size().reset_index(name='count') 
    filtered_df = grouped_df[grouped_df['count'] > 0]
    # filtered_df.to_csv(f'filtered_df_{save_name_index}.txt', sep='\t', index=True, mode='w', chunksize=1000)
    
    # to matrix df
    # Every column is a probe, every row is a barcode
    # Each row (cell) must have at least 100 unique reads
    matrix_df = filtered_df.pivot(index='barcode', columns='probe', values='count').fillna(0).astype(int)
    matrix_df = matrix_df[matrix_df.sum(axis=1) >= 100]
    # write to .csv
    # matrix_df.to_csv(f'count_matrix_f100_{save_name_index}.txt', sep='\t', index=True, mode='w', chunksize=1000)
    matrix_df.to_csv('count_matrix_f100.txt', sep='\t', index=True, mode='w', chunksize=1000)