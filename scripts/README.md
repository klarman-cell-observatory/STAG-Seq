# Scripts used to process HyPR-seq and Tapestri data

## 1. How to process HyPR-seq data
This code is currently written to be processed on the UGER in Broad cluster. To run: 

Edit to the `scripts/process_hyprseq/python_launcher.sh` to activate your own Anaconda environment:

```
.~/miniconda3/etc/profile.d/conda.sh
conda activate ~/miniconda3 > activate_output.txt
```

Run by using:
```
sh process_hyprseq/launch_script_v2.sh  FASTQ_PATH  PROBE_SET_FILE
```

This will generate a `count_matrix_f100.txt` file that is the count of all probes in HyPR-seq dataset per cell barcode. 

## 2. How to call variants from Tapestri data and merge with HyPR-seq data
Have your `.loom` file from the Tapestri pipeline available. 

Run the variant calling script using the following command. This will also merge the HyPR-seq and Tapestri outputs:

```
Usage: 
	python merge/merge_hyprseq_and_tapestri.py -h
Run: 
	python merge/merge_hyprseq_and_tapestri_output.py \
		-r example/HyPR_matrix.txt \  
  		-d example/Tapestri_cells.loom \
  		-o example/outs \
  		-g germline_mutations.txt
```
