#!/bin/bash
#$ -N chunk_job
#$ -pe smp 1
#$ -l h_rt=08:00:00

probe_file=$probe_file
work_dir=$work_dir
cd $work_dir

. /home/unix/gongqiyu/miniconda3/etc/profile.d/conda.sh
conda activate /home/unix/gongqiyu/miniconda3 >activate_output.txt

python /broad/thechenlab/stag_seq/hypr_pipeline/filter_and_group.py -p $probe_file \
    -f extracted.txt \
    -b /broad/thechenlab/stag_seq/hypr_pipeline/barcode_lookup.pkl
