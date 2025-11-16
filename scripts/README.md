## 1. How to process HyPR-seq data
This code is adapted to processed on Terra, the original code was a product of the STAG-Seq repo:

Dockerfile is included in this repo, a file that is not included for external users is the `barcode_lookup.pkl` file, which is included for internal Broad users within the Docker image.

Other required bioinformatics tools: seqtk, seqkit, umi_tools

`environment.yml` has been provided for HPC or local usage.

## To run the pipeline:

Pipeline monitoring and runtime will be managed within a terra.bio (Cromwell) workspace.

Import the stag-seq-rna pipeline from the Terra WDL repository into a workspace of choice:

https://app.terra.bio/#workflows/STAG_Seq/STAG_Seq_RNA/10 (Not currently publically readable due to final polishing)

Fill in the provided `STAG_Seq_TerraInput_Example.csv` with each row in the csv corresponding to the paths of a sample fastq and its corresponding probe file. Upload the csv onto the cloud as well, that will serve as the `input_fastq_probe_table`.

Within the Terra workflow interface fill out the required fields (`copy_intermediates`, `docker_registry`, `input_fastq_probe_table`, `output_directory`, `run_filter_group_only` and scale the computational specifications  (`disk_space`, `memory`, `num_cpu`) to the size of your data.

<img width="1133" height="380" alt="image" src="https://github.com/user-attachments/assets/8ab074b3-13c2-46a1-bf42-eed3ab35b896" />

`FASTQ_PATH` could be either in .fastq or .fastq.gz

`FASTQ_PATH`: `/path/in/cloud/Th17_ROR_Acti_Probe4_S1_R1_001.fastq.gz` 

`PROBE_PATH`: `/path/in/cloud/5-9-2025_probe_RORC.txt` (THIS WILL BE EXPERIMENT SPECIFIC)

Example csv provided here:
<img width="973" height="87" alt="image" src="https://github.com/user-attachments/assets/958b4477-b502-43ca-8ccf-287ec60b6ac0" />

The `copy_intermediates` flag is helpful if you want to export the intermediary files generated within the scripts (barcode/umi/probe extracted fastqs), exported to the output_directory.

The `run_filter_group_only` flag is useful for skipping the lengthy extraction step, granted you have copied the intermediate outputs, it will look for files in the export location specified above in the variable.

This will generate a sample labeled `count_matrix_f100.h5ad` file that is the count of all probes in HyPR-seq dataset per cell barcode. 
