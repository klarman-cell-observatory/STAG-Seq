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
FASTQ_FILE=$FASTQ_FILE
N=$N

cd "$WORKING_DIR"

#############################
/broad/thechenlab/Dawn/software/seqkit split2 $FASTQ_FILE -p $N -o fastq_chunk_ -f -O ./
# Rename the files so that they are of the form: fastq_chunk_0.fastq, fastq_chunk_1.fastq, etc.
# The original files are named fastq_chunk_.part_001.fastq etc.
for file in fastq_chunk_.part_*.fastq; do
    # Extract the part number from the file name and decrement by 1
    part=$(echo $file | sed 's/fastq_chunk_.part_0*\([1-9][0-9]*\).fastq/\1/')
    part=$((part - 1)) # decrement part by 1

    # Construct the new file name using the decremented part number
    newname="fastq_chunk_$part.fastq"
    
    # Rename the file
    mv "$file" "$newname"
done
#############################

# echo "Determining lines per chunk..."

# head -n 1000 "$FASTQ_FILE" > sample.fastq
# sample_lines=$(wc -l < sample.fastq)
# sample_size=$(wc -c < sample.fastq)
# total_size=$(wc -c < "$FASTQ_FILE")
# estimated_total_lines=$(($total_size * $sample_lines / $sample_size))
# total_lines=$(($estimated_total_lines / 4 * 4))  # Adjust to multiple of 4

# # Calculate lines per chunk
# lines_per_chunk=$(($total_lines / ${N}))  # Replace N with the number of chunks you want
# lines_per_chunk=$(($lines_per_chunk / 4 * 4))  # Make sure it's a multiple of 4

# echo "total lines: ${total_lines}"
# echo "lines per chunk: ${lines_per_chunk}"
# echo "Splitting into ${N} CHUNKS"

# split -l $lines_per_chunk $FASTQ_FILE fastq_chunk_

# counter=0
# for f in fastq_chunk_*; do
#     mv "$f" "fastq_chunk_${counter}.fastq"
#     counter=$((counter + 1))
# done

# awk -v total_lines="$total_lines" -v N="${N}" 'BEGIN { lines_per_chunk = int(total_lines / N / 4) * 4; print "Lines per chunk: " lines_per_chunk; } NR % 4 == 1 { x="fastq_chunk_" int((NR/4)/lines_per_chunk) ".fastq"; print "Writing to " x; } { print > x; }' "$FASTQ_FILE"

