#!/bin/bash

help_me_obi1() {
  echo "Usage: $0 [R1.fastq] [probeset.txt] [optional arguments]"
  echo
  echo "  R1.fastq          path to R1 reads from sample"
  echo "  probeset.txt      path to the probeset file"
  echo
  echo "Optional Arguments:"
  echo "  -w  working directory"
  echo "  -n  number of jobs to launch"
  echo "  -c  (not implemented yet) cluster / computer type. {local, slurm, sge, auto}. automatically detects by default."
  echo "  -p  amount of memory to request per extraction job (in GB)"
  echo "  -m  amount of memory to request for mismatch and groupby script (in GB)"
  echo "  -u  (not implemented yet) filter out barcodes with fewer than N UMIs. Default is 'auto' which selects N. Set to 0 to do no filtering. "
  echo "  -r  R2.fastq if R1 does not have 150 cycles"
  echo "  -o  directory where stdout and stderr are streamed to"
  echo "  -s  {1,2,3,4}. 1 runs everything. 2 runs extract onwards. 3 runs cat onwards. 4 runs Python launcher only."
  echo "  -h  help"
  echo
}

if [ "$#" -lt 2 ]; then
  help_me_obi1
  exit 1
fi

# CURRENTLY IVE ONLY WRITTEN THIS FOR UGER / SUN GRID ENGINE

FASTQ_FILE="$1"
PROBE_FILE="$2"
shift 2

fastq_name=$(basename $FASTQ_FILE | rev | cut -d '.' -f2- | rev)

# set the working dir to be fastq_name by default
s=1
work_dir="${PWD}/${fastq_name}_work_dir"
n_jobs=5
mem_extract=10
mem_final=100
R2=""
output_dir="${PWD}/${fastq_name}_output_dir"

# Function to unzip and get the new filename
gunzip_file() {
  INPUT_FILE=$1
  # Check if the file is gzipped
  if [[ "$INPUT_FILE" == *.gz ]]; then
    # Unzip the file
    gunzip "$INPUT_FILE"
    # Remove the .gz from the name to get the new name
    INPUT_FILE=${INPUT_FILE%.gz}
  fi
  # Return the (possibly new) filename
  echo "$INPUT_FILE"
}

# parse the optional arguments
while getopts "w:n:s:o:m:r:h" opt; do
  case $opt in

  w)
    work_dir="$OPTARG"
    ;;
  n)
    n_jobs="$OPTARG"
    ;;
  s)
    s="$OPTARG"
    ;;
  p)
    mem_extract="$OPTARG"
    ;;
  m)
    mem_final="$OPTARG"
    ;;
  r)
    R2="$OPTARG"
    ;;
  o)
    output_dir="$OPTARG"
    ;;
  h)
    display_help
    exit 0
    ;;
  \?)
    echo "Invalid option: -$OPTARG" >&2
    exit 1
    ;;
  :)
    echo "Option -$OPTARG requires an argument." >&2
    exit 1
    ;;
  esac
done

echo "parsed arguments"

# make sure the .fastq and probe.txt were provided
if [ -z "$FASTQ_FILE" ] || [ -z "$PROBE_FILE" ]; then
  echo "Error: Both R1.fastq and probe .txt file must be provided."
  exit 1
fi

echo "extracting .fastqs"
# Unzip files and get new names
FASTQ_FILE=$(gunzip_file "$FASTQ_FILE")
# Only unzip R2 if it is provided
if [ -n "$R2" ]; then
  R2=$(gunzip_file "$R2")
fi

current_dir=$(pwd)
source /broad/software/scripts/useuse
use UGER

if [[ "$FASTQ_FILE" == *.gz ]]; then
  echo "unzipping the .fastq file..."
  gunzip "$FASTQ_FILE"
  echo "unzipped the file"
fi

echo "loaded UGER successfully"

mkdir "$work_dir"

# we will do the umi, probe, barcode extraction differently if we have R2 or just R1:

if [[ "$s" -le 1 ]]; then
  echo "split the .fastq into ${n_jobs} chunks"

  if [ -z "$R2" ]; then
    echo "No read 2 .fastq was provided. Extracting everything from R1."

    # submit the job that splits the file
    split_job_name="split_${fastq_name}"
    qsub -N "$split_job_name" -l h_vmem=${mem_extract}G -o ${output_dir}/split.out -e ${output_dir}/split.err -v N="$n_jobs",FASTQ_FILE="$FASTQ_FILE",work_dir="$work_dir" split_script.sh
  else
    echo "Read 2 .fastq was provided. Extracting umi + probe from R2."

    # submit the job that splits the file
    split_job_name="split_${fastq_name}"

    qsub -N "$split_job_name" -l h_vmem=${mem_extract}G -o ${output_dir}/split.out -e ${output_dir}/split.err -v N="$n_jobs",FASTQ_FILE_R1="$FASTQ_FILE",FASTQ_FILE_R2="$R2",work_dir="$work_dir" split_script_two_reads.sh
  fi
fi

if [[ "$s" -le 2 ]]; then
  if [ -z "$R2" ]; then
    # submit the jobs that do the extraction
    JOB_RANGE="1-$n_jobs"
    ARRAY_JOB_ID=$(qsub -N "${fastq_name}_extract" -t $JOB_RANGE -hold_jid "$split_job_name" -o ${output_dir}/extract.out -e ${output_dir}/extract.err -v FASTQ_PREFIX="fastq_chunk_",work_dir="$work_dir" chunk_script.sh | awk '{print $3}')
  else
    # submit the jobs that do the extraction
    JOB_RANGE="1-$n_jobs"
    ARRAY_JOB_ID=$(qsub -N "${fastq_name}_extract" -t $JOB_RANGE -hold_jid "$split_job_name" -o ${output_dir}/extract.out -e ${output_dir}/extract.err -v FASTQ_R1_PREFIX="fastq_1_chunk_",FASTQ_R2_PREFIX="fastq_2_chunk_",work_dir="$work_dir" chunk_script_two_reads.sh | awk '{print $3}')
  fi
fi

if [[ "$s" -le 3 ]]; then
  # submit the job that concatenates the outputs from the above jobs
  qsub -N "${fastq_name}_cat" -hold_jid "${fastq_name}_extract" -o ${output_dir}/cat.out -e ${output_dir}/cat.err -v work_dir="$work_dir" cat_job.sh
fi

if [[ "$s" -le 4 ]]; then
  # submit the job that does the filtering and groupby in Python
  qsub -N "${fastq_name}_python" -l h_vmem=${mem_final}G -hold_jid "${fastq_name}_cat" -o ${output_dir}/python_run.out -e ${output_dir}/python_run.err -v work_dir="$work_dir",probe_file="$PROBE_FILE" python_launcher.sh
fi
