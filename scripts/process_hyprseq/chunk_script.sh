#!/bin/bash
#$ -pe smp 1
#$ -l h_rt=08:00:00  # adjust walltime as needed

FASTQ_PREFIX=$FASTQ_PREFIX
work_dir=$work_dir
cd $work_dir

ADJUSTED_TASK_ID=$((SGE_TASK_ID - 1))

INPUT_FILE="${FASTQ_PREFIX}${ADJUSTED_TASK_ID}.fastq"

if [[ ! -f "$INPUT_FILE" ]]; then
	echo "Input file $INPUT_FILE does not exist. Exiting."
	exit 1
fi

echo "$SGE_TASK_ID"

# activate python environment

. /home/unix/gongqiyu/miniconda3/etc/profile.d/conda.sh
conda activate /home/unix/gongqiyu/miniconda3 >activate_output.txt

# first we use umi_tools extract to pull the cell barcode and place it in the header of an output .fastq the cell barcode is 18bp (bc1 + bc2). bc1 is located at the first 9bp of the read. following this there is a constant region AGTAAGTACGAGTC  with a variable offset between 0 and 3 bps. Follwing this region is bc2.
# we allow up to two mismatchs in the constant region

echo "extracting barcodes..."

umi_tools extract --stdin "$INPUT_FILE" \
	--stdout "bc_extracted_${ADJUSTED_TASK_ID}.fastq" \
    --whitelist /broad/thechenlab/stag_seq/hypr_pipeline/whitelist_missionbio.txt \
	--extract-method "regex" \
	--bc-pattern="(?P<cell_1>.{9})(?P<discard_1>.{0,3}AGTAAGTACGAGTC){s<=2}(?P<cell_2>.{9})(?P<umi_dummy>.{2}).*" \
	--filtered-out "no_bc_match_${ADJUSTED_TASK_ID}.fastq" \
	-L "bc_extract_${ADJUSTED_TASK_ID}.log" \
	--subset-reads=300000000 || exit 1

# we now need to extract the probe and UMI. These are in a 35bp region that follows (towards 3' end) the sequence GTACTCGCAGTAGTCTCG. We extract by taking the 35bps after that constant sequence (I should change this to allow mismatches).

echo "extracting UMI + PROBE..."

awk 'NR % 4 == 2 {
        start = index($0, "CGCAGTAGTCTCG")
        if (start) {
            seq = substr($0, start+13, 35)
            getline plus
            getline qscore
            print header"\n"seq"\n"plus"\n"substr(qscore, start+13, 35)
        }
    }
    NR % 4 == 1 {
        header = $0
    }' "bc_extracted_${ADJUSTED_TASK_ID}.fastq" >"bc_umi_probe_extracted_${ADJUSTED_TASK_ID}.fastq"

# since we are reading R1 we need to reverse complement the 35bp sequence to match it to our known probes. A fast way to do this is with seqtk

echo "taking reverse complement..."

/broad/thechenlab/Dawn/software/seqtk/seqtk seq -r "bc_umi_probe_extracted_${ADJUSTED_TASK_ID}.fastq" >"rev_bc_umi_probe_extracted_${ADJUSTED_TASK_ID}.fastq"

# this output file is barcode, umi, probe

echo "creating output file..."

awk -F'_' 'NR % 4 == 1 {
        barcode = $2
        getline seq
        print barcode, substr(seq, length(seq)-9), substr(seq, 1, 25)
    }' "rev_bc_umi_probe_extracted_${ADJUSTED_TASK_ID}.fastq" >"${ADJUSTED_TASK_ID}_output.txt"
