#!/bin/bash
INPUT_FILE=$1
RAW_FASTQ_NAME=$2
BASE_NAME=$3

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Input file $INPUT_FILE does not exist. Exiting."
    exit 1
fi

echo "Extracting barcodes for ${BASE_NAME}..."

zcat "${INPUT_FILE}" | umi_tools extract \
    --stdout="/mnt/disks/cromwell_root/bc_extracted_${BASE_NAME}.fastq" \
    --whitelist=/ref/whitelist_missionbio.txt \
    --extract-method=regex \
    --bc-pattern="(?P<cell_1>.{9})(?P<discard_1>.{0,3}AGTAAGTACGAGTC){s<=2}(?P<cell_2>.{9})(?P<umi_dummy>.{2}).*" \
    --filtered-out="/mnt/disks/cromwell_root/no_bc_match_${BASE_NAME}.fastq" \
    -L "bc_extract_${BASE_NAME}.log" \
    -v 2 \
    --subset-reads=300000000

# we now need to extract the probe and UMI. These are in a 35bp region that follows (towards 3' end) the sequence GTACTCGCAGTAGTCTCG. We extract by taking the 35bps after that constant sequence (I should change this to allow mismatches).

if [[ ! -f "/mnt/disks/cromwell_root/bc_extracted_${BASE_NAME}.fastq" ]]; then
    echo "Extraction step did not complete successfully, barcode extracted file was not found."
    exit 1
fi

echo "Extracting UMI + PROBE..."

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
    }' "/mnt/disks/cromwell_root/bc_extracted_${BASE_NAME}.fastq" >"/mnt/disks/cromwell_root/bc_umi_probe_extracted_${BASE_NAME}.fastq"

# since we are reading R1 we need to reverse complement the 35bp sequence to match it to our known probes. A fast way to do this is with seqtk

echo "Taking reverse complement..."

seqtk seq -r "/mnt/disks/cromwell_root/bc_umi_probe_extracted_${BASE_NAME}.fastq" >"/mnt/disks/cromwell_root/rev_bc_umi_probe_extracted_${BASE_NAME}.fastq"

# this output file is barcode, umi, probe

echo "Creating output file..."

awk -F'_' 'NR % 4 == 1 {
        barcode = $2
        getline seq
        print barcode, substr(seq, length(seq)-9), substr(seq, 1, 25)
    }' "/mnt/disks/cromwell_root/rev_bc_umi_probe_extracted_${BASE_NAME}.fastq" >"/mnt/disks/cromwell_root/${BASE_NAME}_output.txt"
