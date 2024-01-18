#!/bin/bash
#$ -N chunk_job
#$ -pe smp 1
#$ -l h_rt=08:00:00  # adjust walltime as needed

FASTQ_R1_PREFIX=$FASTQ_R1_PREFIX
FASTQ_R2_PREFIX=$FASTQ_R2_PREFIX

work_dir=$work_dir
cd $work_dir

ADJUSTED_TASK_ID=$((SGE_TASK_ID - 1))

INPUT_FILE_R1="${FASTQ_R1_PREFIX}${ADJUSTED_TASK_ID}.fastq"
INPUT_FILE_R2="${FASTQ_R2_PREFIX}${ADJUSTED_TASK_ID}.fastq"

echo "$SGE_TASK_ID"

# activate python environment

. /home/unix/arasmuss/anaconda3/etc/profile.d/conda.sh
conda activate /home/unix/arasmuss/anaconda3/envs/HyPR_env > activate_output.txt

# first we use umi_tools extract to pull the cell barcode and place it in the header of an output .fastq the cell barcode is 18bp (bc1 + bc2). bc1 is located at the first 9bp of the read. following this there is a constant region AGTAAGTACGAGTC  with a variable offset between 0 and 3 bps. Follwing this region is bc2. 
# we allow up to two mismatchs in the constant region 

echo "extracting barcodes..."

umi_tools extract \
    --stdin "$INPUT_FILE_R1" \
    --read2-in "$INPUT_FILE_R2" \
    --whitelist /broad/thechenlab/Dawn/HyPR_Tap/mix_HyPR_fastqs/whitelist_missionbio.txt \
    --extract-method "regex" \
    --bc-pattern="(?P<cell_1>.{9})(?P<discard_1>.{0,3}AGTAAGTACGAGTC){s<=2}(?P<cell_2>.{9}).*" \
    --bc-pattern2="[ACGT]{3}(?P<umi_1>.{35})" \
    --stdout "bc_extracted_${ADJUSTED_TASK_ID}.fastq" \
    -L "bc_extract_${ADJUSTED_TASK_ID}.log" \
    --subset-reads=300000000 || exit 1

awk 'NR%4==1 {split($0,a,"_"); print a[length(a)-1], substr(a[length(a)], 26, 10), substr(a[length(a)], 1, 25)}' bc_extracted_${ADJUSTED_TASK_ID}.fastq > "${ADJUSTED_TASK_ID}_output.txt"

#awk 'NR%4==1 {barcode=substr($0, 2, 9); umi=substr($0, 12, 35); last_10=substr(umi, 26, 10); first_25=substr(umi, 1, 25); print barcode "\t" last_10 "\t" first_25}' bc_extracted_${ADJUSTED_TASK_ID}.fastq > "${ADJUSTED_TASK_ID}_output.txt"
