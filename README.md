STAG-Seq
Advancements in genomics have revolutionized human genetics by uncovering the genetic architecture of disease. However, identifying causal variants and their mechanisms of action remains a challenge for therapeutic translation. Here, we integrate variant editing with simultaneous measurements of genomic DNA and RNA to phenotype variants associated with inborn errors of immunity and autoimmune disease. Sensitive Transcriptomics And Genotyping by sequencing (STAG-seq) enabled high-content variant functionalization in primary human immune cells while controlling for allele-dose effects and editing precision.

Accordingly, we systematically defined functions for rare coding variants associated with monogenic diseases of immune dysregulation that impact type II interferon signaling in macrophages, thus defining mechanisms and explaining modes of inheritance. Additionally, we elucidated the function of an autoimmunity-linked pleiotropic noncoding variant by defining its cis-regulatory role in TNRC18 expression in T cells and ensuing trans effects on endogenous retroviral element derepression. Finally, we functionally dissected a credible set of noncoding variants in linkage disequilibrium associated with autoimmunity, implicating a causal variant that increases TCF7 expression in T cells and modulates cell fate decisions skewing memory versus effector status.

Collectively, we provide a framework to bridge the gap between genetic associations and functional mechanisms to advance our understanding of disease biology.

This repository contains the scripts and functions necessary to process, merge, and analyze this data, allowing for comprehensive insights into cellular heterogeneity.


🧬 Data Processing
HyPR-seq Data
To process your HyPR-seq data, you will need to be in a UGER environment, such as the Broad cluster.

Activate your Conda Environment: Edit the scripts/process_hyprseq/python_launcher.sh file to point to your Anaconda environment.

.~/miniconda3/etc/profile.d/conda.sh
conda activate ~/miniconda3 > activate_output.txt

Run the Processing Script: Execute the following command to process your FASTQ files.

sh process_hyprseq/launch_script_v2.sh [FASTQ_PATH] [PROBE_SET_FILE]

This script will generate a count_matrix_f100.txt file, which contains the probe counts for each cell barcode. This matrix can be used to construct a Scanpy or Seurat object for downstream analysis.

Tapestri Data and Merging
Once you have the .loom file from the Tapestri pipeline, you can call variants and merge the data with your HyPR-seq results.

Use the following command to run the merging script:

python merge/merge_hyprseq_and_tapestri_output.py \
    -r example/HyPR_matrix.txt \
    -d example/Tapestri_cells.loom \
    -o example/outs \
    -g germline_mutations.txt

For more information on the script's arguments, you can use the -h flag.

📊 Analysis Notebooks for Paper Figures
The notebooks directory contains the Jupyter and R notebooks used to generate the figures for our paper. Each subdirectory corresponds to a specific figure.

Figure 1: Sensitivity
fig1_e-g.ipynb

fig1_h-k.ipynb

STAGseq_10X_Fig1.Rmd

THP1_Jurkat_Fig1.Rmd

Figure 2: Matrix Screen
preprocess_matrix_screen.ipynb

matrix_screen_DE.ipynb

make_quadrant_plot.ipynb

Figure 3: DSIFNG
anders_version_preprocess_IFNG.ipynb

figure3b_umap.ipynb

figure3c_boxplot.ipynb

figure3d-i_ifng_score.ipynb

figure3k.ipynb

figure3k.R

Figure 4: TNRC18
preprocess_TNRC18.ipynb

make_figures.ipynb

figure4d.ipynb

figure4fgh.ipynb

Figure 5: TCF7
fig5_c_efg.ipynb