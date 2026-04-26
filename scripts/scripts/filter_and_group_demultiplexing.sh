import numpy as np
import pandas as pd
import sys
import argparse
import anndata as ad
import scipy.sparse as sp

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
    parser.add_argument("-p", "--probe_file", type=str, required=True, help="Probe set file")
    parser.add_argument("-f", "--in_file", type=str, required=True, help="Concatenated output file")
    parser.add_argument("-b", "--barcode_hash", type=str, required=True, help="Barcode hash table")
    parser.add_argument("-s", "--save_name_index", type=str, required=False, help="Save name index", default="test")
    parser.add_argument("-m", "--mapping_file", type=str, required=False, help="Barcode to batch mapping CSV")
    
    # NEW: Argument for the multiplex definition table
    parser.add_argument("-x", "--multiplex_file", type=str, required=True, help="Multiplex lookup TSV (Seq in last column)")
    
    args = parser.parse_args()
    return args

def create_probe_lookup(probe_file):
    """Create a lookup table for probes and their mismatches."""
    correct_sequences = {}
    with open(probe_file, 'r') as f:
        for line in f:
            try:
                seq_id, seq = line.strip().split('\t')
                correct_sequences[seq] = seq_id
            except:
                print(line)
    print(str(len(correct_sequences.keys())) + " probe sequences identified")

    lookup_table = {}
    for seq, seq_id in correct_sequences.items():
        lookup_table[seq] = seq_id
        for mismatch_seq in one_mismatch_sequences(seq):
            lookup_table[mismatch_seq] = seq_id
    return lookup_table

def generate_multiplex_lookup(multiplex_file):
    """Dynamic lookup for multiplex indices from a TSV file."""
    raw_indices = {}
    
    # Read the file dynamically: assumes the actual DNA sequence is the LAST column, 
    # and joins the previous columns to make a descriptive name.
    with open(multiplex_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines or header lines if they contain 'Need to be'
            if not line or "Need to be" in line:
                continue
                
            parts = line.split('\t')
            if len(parts) >= 2:
                raw_seq = parts[-1].strip()
                # Join the remaining columns with underscores to make a clean condition name
                condition_name = "_".join([p.strip() for p in parts[:-1]]).replace(" ", "_")
                raw_indices[raw_seq] = condition_name

    # Helper to generate reverse complement
    def rc(seq):
        comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        return "".join(comp.get(b, b) for b in reversed(seq))
        
    mux_lookup = {}
    for raw_seq, condition_name in raw_indices.items():
        # Get the actual sequence we expect to see in the read (Reverse Complement)
        true_seq = rc(raw_seq)
        mux_lookup[true_seq] = condition_name
        
        # Add all 1-bp mismatch permutations
        for mismatch in one_mismatch_sequences(true_seq):
            if mismatch not in mux_lookup:
                mux_lookup[mismatch] = condition_name
                
    print(f"{len(raw_indices.keys())} multiplex conditions loaded and reverse-complemented.")
    return mux_lookup

if __name__ == "__main__":
    # Parse arguments
    args = parse_args()
    probe_file = args.probe_file
    in_file = args.in_file
    barcode_hash = args.barcode_hash
    save_name_index = args.save_name_index
    mapping_file = args.mapping_file
    multiplex_file = args.multiplex_file # NEW
    
    print("probe_file: " + probe_file)
    print("in_file: " + in_file)
    print("barcode_hash: " + barcode_hash)
    if mapping_file:
        print("mapping_file: " + mapping_file)
    if multiplex_file:
        print("multiplex_file: " + multiplex_file)
    
    # Load in the main data
    sys.stdout.write("reading concatenated output into dataframe\n")
    df = pd.read_csv(in_file, sep=' ', names=['barcode', 'umi', 'probe'])

    # Probe lookup table
    sys.stdout.write("creating probe lookup table\n")
    lookup_table = create_probe_lookup(probe_file)

    # Use the probe lookup to replace mismatches
    sys.stdout.write("replacing probe mismatches\n")
    sys.stdout.write(str(len(df)) + " reads before filtering\n")
    df['probe'] = df['probe'].map(lookup_table)
    df.dropna(subset=['probe'], inplace=True)
    sys.stdout.write(str(len(df)) + " reads have valid probes\n")

    # Drop duplicate from dataframe
    sys.stdout.write("dropping duplicates\n")
    sys.stdout.write(str(len(df)) + " reads before dropping duplicates\n")
    df = df.drop_duplicates(subset=['probe', 'umi', 'barcode'])
    sys.stdout.write(str(len(df)) + " reads after dropping duplicates\n")
    
    # GROUPBY
    grouped_df = df.groupby(['probe', 'barcode']).size().reset_index(name='count') 
    filtered_df = grouped_df[grouped_df['count'] > 0]
    
    # To matrix df
    matrix_df = filtered_df.pivot(index='barcode', columns='probe', values='count').fillna(0).astype(int)
    sparse_matrix = sp.csr_matrix(matrix_df.values, dtype='int32')
    adata = ad.AnnData(X=sparse_matrix, obs=pd.DataFrame(index=matrix_df.index), var=pd.DataFrame(index=matrix_df.columns))
    
    # ==========================================
    # Multiplexing Side-Channel Error Correction
    # ==========================================
    if mapping_file and multiplex_file:
        sys.stdout.write("processing multiplexing batch IDs with 1bp error correction\n")
        
        batch_df = pd.read_csv(mapping_file, sep=',', names=['barcode', 'raw_idx'])
        mux_lookup = generate_multiplex_lookup(multiplex_file)
        
        batch_df['condition'] = batch_df['raw_idx'].map(mux_lookup)
        batch_df.dropna(subset=['condition'], inplace=True)
        
        majority_condition = batch_df.groupby('barcode')['condition'].agg(
            lambda x: x.value_counts().index[0]
        )
        
        adata.obs['experiment_condition'] = adata.obs_names.map(majority_condition)
        adata.obs['experiment_condition'] = adata.obs['experiment_condition'].fillna('Unknown')
        
        sys.stdout.write("successfully attached experiment_condition to adata.obs\n")
    elif mapping_file and not multiplex_file:
        sys.stdout.write("WARNING: mapping_file provided but no multiplex_file (-x) provided. Skipping multiplex annotation.\n")
    
    # Save the final file
    adata.write("count_matrix_f100.h5ad")