#!/bin/bash

#$ -l h_vmem=5G

# Cores
#$ -pe smp 1
#$ -binding linear:1

# I like single output files
#$ -j y

# Runtime request.  Usually 30 minutes is plenty for me and helps me get backfilled into reserved slots.
#$ -l h_rt=48:00:00

# I don't like the top level of my homedir filling up.
#$ -o $HOME/outputs/

# Job name.
#$ -N z_test

WORKING_DIR=$work_dir
FASTQ_FILE_R1=$FASTQ_FILE_R1
FASTQ_FILE_R2=$FASTQ_FILE_R2

N=$N

cd "$WORKING_DIR"

echo "Determining lines per chunk..."

head -n 1000 "$FASTQ_FILE_R1" > sample.fastq
sample_lines=$(wc -l < sample.fastq)
sample_size=$(wc -c < sample.fastq)
total_size=$(wc -c < "$FASTQ_FILE_R1")
estimated_total_lines=$(($total_size * $sample_lines / $sample_size))
total_lines=$(($estimated_total_lines / 4 * 4))  # Adjust to multiple of 4

# Calculate lines per chunk
lines_per_chunk=$(($total_lines / ${N}))  # Replace N with the number of chunks you want
lines_per_chunk=$(($lines_per_chunk / 4 * 4))  # Make sure it's a multiple of 4

echo "total lines: ${total_lines}"
echo "lines per chunk: ${lines_per_chunk}"
echo "Splitting into ${N} CHUNKS"

split -l $lines_per_chunk $FASTQ_FILE_R1 fastq_1_chunk_
split -l $lines_per_chunk $FASTQ_FILE_R2 fastq_2_chunk_

counter=0
for f in fastq_1_chunk_*; do
    mv "$f" "fastq_1_chunk_${counter}.fastq"
    counter=$((counter + 1))
done

counter=0
for f in fastq_2_chunk_*; do
    mv "$f" "fastq_2_chunk_${counter}.fastq"
    counter=$((counter + 1))
done

