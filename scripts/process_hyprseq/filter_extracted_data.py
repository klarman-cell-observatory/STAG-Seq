import pandas as pd
import os
import sys
import random


def filter_extracted_data(extracted_file, barcodes_file, output_file):
    """
    Filters the extracted data based on the barcodes and writes the result to a new file,
    allowing for duplicate barcodes.
    """
    if os.path.exists(output_file):
        print(f"{output_file} already exists. Skipping file generation")
        return
    
    barcodes = read_barcodes(barcodes_file)

    with open(extracted_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.split()
            if parts[0] in barcodes:
                outfile.write(line)


def read_barcodes(file_path):
    """
    Reads the barcode file and returns a set of barcodes.
    """
    with open(file_path, 'r') as file:
        barcodes = set(line.strip() for line in file)
    return barcodes


def sample_lines(input_file, sample_size, output, seed=1234):
    """
    Samples a given number of lines from the input file and writes them to a new file.
    Sets a random seed for reproducibility.
    """
    random.seed(seed) 
    
    if sample_size <= 0:
        return

    with open(input_file, 'r') as file:
        lines = file.readlines()

    sampled_lines = random.sample(lines, min(sample_size, len(lines)))

    with open(output, 'w') as file:
        file.writelines(sampled_lines)
        

def get_info_from_command_line():
    index = int(sys.argv[1])
    paths = [
        'MDM_MixingP500_HyPR_S1_R1_001_work_dir/', 
        'MDM_P1000_HyPR_S2_R1_001_work_dir/', 
        'MDM_P2000_HyPR_S3_R1_001_work_dir/'
    ]
    # sample_reads = [0, 56299520, 112599040]
    # sample_reads = [28149760, 28655483, 104566108]
    sample_reads = [28149760, 14283107, 52094261]
    
    command_line_1 = [
        'python /broad/thechenlab/stag_seq/hypr_pipeline/filter_and_group.py -p /broad/thechenlab/stag_seq/figure_1/sample_fastqs/p500_MDM/p500_probes.txt -f /broad/thechenlab/stag_seq/hypr_pipeline/MDM_MixingP500_HyPR_S1_R1_001_work_dir/extracted_intersect_sample.txt -b /broad/thechenlab/stag_seq/hypr_pipeline/barcode_lookup.pkl',
        'python /broad/thechenlab/stag_seq/hypr_pipeline/filter_and_group.py -p /broad/thechenlab/stag_seq/figure_1/sample_fastqs/p1000_MDM/p1000_probes.txt -f /broad/thechenlab/stag_seq/hypr_pipeline/MDM_P1000_HyPR_S2_R1_001_work_dir/extracted_intersect_sample.txt -b /broad/thechenlab/stag_seq/hypr_pipeline/barcode_lookup.pkl',
        'python /broad/thechenlab/stag_seq/hypr_pipeline/filter_and_group.py -p /broad/thechenlab/stag_seq/figure_1/sample_fastqs/p2000_MDM/p2000_probes.txt -f /broad/thechenlab/stag_seq/hypr_pipeline/MDM_P2000_HyPR_S3_R1_001_work_dir/extracted_intersect_sample.txt -b /broad/thechenlab/stag_seq/hypr_pipeline/barcode_lookup.pkl'
    ]
    command_line_2 = [
        'mv count_matrix_f100.txt MDM_MixingP500_HyPR_S1_R1_001_work_dir/count_matrix_sample_p500_new.txt',
        'mv count_matrix_f100.txt MDM_P1000_HyPR_S2_R1_001_work_dir/count_matrix_sample_p1000_new.txt',
        'mv count_matrix_f100.txt MDM_P2000_HyPR_S3_R1_001_work_dir/count_matrix_sample_p2000_new.txt'
    ]
    return paths[index], sample_reads[index], command_line_1[index], command_line_2[index]


os.chdir("/broad/thechenlab/stag_seq/hypr_pipeline/")
path, sample_size,  command_line_1, command_line_2 = get_info_from_command_line()
print(path)
print('sample_size:',sample_size)
filter_extracted_data(path+'extracted.txt', path+'barcode_to_sample.tsv', path+'extracted_intersect.txt')
sample_lines(path+'extracted_intersect.txt', sample_size, path+'extracted_intersect_sample.txt')
os.system(command_line_1)
os.system(command_line_2)