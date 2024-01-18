import numpy as np
import pandas as pd
import pickle
from matplotlib import pyplot as plt
import sys

probe_file = sys.argv[1]
in_file =  sys.argv[2]
barcode_hash = sys.argv[3]

print("probe file:")
print(probe_file)
print("in file")
print(in_file)
print("hash table:")
print(barcode_hash)


# load in the data

sys.stdout.write("reading concatenated output into dataframe\n")
df = pd.read_csv(in_file, sep=' ', names=['barcode', 'umi', 'probe'])

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

# Read the file of correct sequences into a dictionary
sys.stdout.write("create hash table of probe mismatches \n")

correct_sequences = {}
with open(probe_file, 'r') as f:
    for line in f:
        try:
            seq_id, seq = line.strip().split('\t')
            correct_sequences[seq] = seq_id
        except:
            print(line)
sys.stdout.write(str(len(correct_sequences.keys())) + " probe sequences identified \n")


# Create a lookup table
lookup_table = {}
for seq, seq_id in correct_sequences.items():
    lookup_table[seq] = seq_id
    #for mismatch_seq in one_mismatch_sequences(seq):
        #lookup_table[mismatch_seq] = seq_id

# Use the probe lookup to replace mismatches and discard anything with more than 1 mismatch

sys.stdout.write("replacing probe mismatches\n")

sys.stdout.write(str(len(df)) + " reads before filtering\n")

df['probe'] = df['probe'].map(lookup_table)
df.dropna(subset=['probe'], inplace=True)

sys.stdout.write(str(len(df)) + " reads have valid probes\n")




sys.stdout.write("loading the whitelist hash table\n")
# load the barcode lookup table
with open(barcode_hash, 'rb') as f:
    barcode_table = pickle.load(f)

sys.stdout.write("replacing barcode mismatches\n")
# Use the barcode lookup to replace mismatches and discard anything with more than 1 mismatch
df['barcode'] = df['barcode'].map(barcode_table)
df.dropna(subset=['barcode'], inplace=True)





# sys.stdout.write("triplets before dropping duplicates " + str(len(df)) + "\n")

# sys.stdout.write("drop duplicates\n")
# Drop duplicate from dataframe
df = df.drop_duplicates(subset=['probe', 'umi', 'barcode'])

sys.stdout.write("triplets after dropping duplicates " + str(len(df)) + "\n")

sys.stdout.write("identifying valid barcodes\n")
# check the histogram for local minima
grouped_barcode = df.groupby(['barcode']).size().reset_index(name='count')
counts, bin_edges = np.histogram(grouped_barcode['count'].values, bins=75)

print("sum of counts before bc filtering: " + str(np.sum(grouped_barcode['count'].values)))

local_minima_indices = []
for i in range(2, len(counts) - 2):
    if counts[i] < counts[i - 1] and counts[i] < counts[i + 1]:
        if counts[i] < counts[i - 2] and counts[i] < counts[i + 2]:
            local_minima_indices.append(i)
local_minima_bin_edges = [bin_edges[i] for i in local_minima_indices]

# could plot the histogram here if you want
#plt.hist(grouped_barcode['count'].values, bins=75, alpha=0.5, edgecolor='black', log=True)
#plt.axvline(x=local_minima_bin_edges[0], color='r', linestyle='--')
#plt.show()
sys.stdout.write(str(len(grouped_barcode)) + " barcodes identified\n")


# # filter out barcodes that have a number of UMIs below the threshold
# #valid_barcodes = grouped_barcode[grouped_barcode['count'] >= local_minima_bin_edges[0]]['barcode']
# valid_barcodes = grouped_barcode[grouped_barcode['count'] >= 2000]['barcode']
# sys.stdout.write(str(len(valid_barcodes)) + " barcodes passed your threshold\n")
# # filtered_df = df[df['barcode'].isin(valid_barcodes)] #qiyu modify
# sys.stdout.write("count UMIs per barcode per probe\n")


# GROUPBY
grouped_df = df.groupby(['probe', 'barcode']).size().reset_index(name='count') #qiyu modify

# to matrix df
matrix_df = grouped_df.pivot(index='barcode', columns='probe', values='count').fillna(0).astype(int)

sys.stdout.write("write the matrix out\n")
# write to .csv
matrix_df.to_csv('count_matrix_all.txt', sep='\t', index=True, mode='w', chunksize=1000)


