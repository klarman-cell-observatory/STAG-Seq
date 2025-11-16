# STAG-seq

Here, we integrate variant editing with simultaneous measurements of genomic DNA and RNA to phenotype variants associated with inborn errors of immunity and autoimmune disease.

Sensitive Transcriptomics And Genotyping by sequencing (**STAG-seq**) enabled high-content variant functionalization in primary human immune cells while controlling for allele-dose effects and editing precision.

Collectively, we provide a framework to bridge the gap between genetic associations and functional mechanisms to advance our understanding of disease biology.

This repository contains the scripts and functions necessary to **process, merge, and analyze** this data, allowing for comprehensive insights into cellular heterogeneity.

---

## 🧬 Workflow Overview

<p align="center">
  <img src="https://github.com/klarman-cell-observatory/STAG-Seq/blob/main/python_functions/Github.jpg" alt="STAG-seq workflow" width="600"/>
  <br/>
  <b>Fig 1.</b> A brief illustration of the STAG-seq workflow.
</p>
---

## 🧬 Data Processing

Raw STAG-seq RNA modality FASTQ data (`.fastq.gz`) was processed into cell-by-probe UMI count matrices (`count_matrix_f100.h5ad`) using a reproducible, containerized WDL pipeline executed on the Terra.bio platform.

The complete pipeline, including the WDL workflow, Docker environment, and detailed, step-by-step instructions for execution, is available in the `scripts` directory of this repository. Please refer to the `README.md` file within that folder for the full data processing guide.

---

## 📊 Analysis Notebooks for Paper Figures

The `notebooks/` directory contains the **Jupyter** and **R** notebooks used to generate the figures for our paper. Each subdirectory corresponds to a specific figure.

### Figure 1: Sensitivity

* `fig1_e-g.ipynb`
* `fig1_h-k.ipynb`
* `STAGseq_10X_Fig1.Rmd`
* `THP1_Jurkat_Fig1.Rmd`
* `count_tapestri_S7.csv`
* `count_tapestri_S8.csv`


### Figure 2: Matrix Screen

* `preprocess_matrix_screen.ipynb`
* `matrix_screen_DE.ipynb`
* `make_quadrant_plot.ipynb`

### Figure 3: DS IFNγ

* `figure3_preprocess.ipynb`
* `figure3b_umap.ipynb`
* `figure3c_boxplot.ipynb`
* `figure3d-i_ifng_score.ipynb`
* `figure3k.ipynb`
* `figure3k.R`
* `variants_config.ini`

### Figure 4: TNRC18

* `preprocess_TNRC18.ipynb`
* `make_figures.ipynb`
* `figure4d.ipynb`
* `figure4fgh.ipynb`
* `figure4_preprocess.ipynb`
* `TNRC18.ini`

### Figure 5: TCF7

* `fig5_c_efg.ipynb`
* `fig5_preprocess.ipynb`
* `TCF7.ini`

