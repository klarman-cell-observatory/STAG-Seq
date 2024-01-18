#!/bin/bash
#$ -N cat_job
#$ -pe smp 1
#$ -l h_rt=08:00:00  # adjust walltime as needed

work_dir=$work_dir
cd $work_dir
cat *output.txt > extracted.txt
