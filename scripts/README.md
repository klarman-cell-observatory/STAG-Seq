# Scripts used to process HyPR-seq and Tapestri data

## How to process HyPR-seq data
This code is currently written to be processed on the UGER in Broad cluster. To run: 

Edit to the `scripts/process_hyprseq/python_launcher.sh` to activate your own Anaconda environment:

```
. /home/unix/gongqiyu/miniconda3/etc/profile.d/conda.sh
conda activate /home/unix/gongqiyu/miniconda3 >activate_output.txt
```

Run by using:
```
sh launch_script_v2.sh  FASTQ_PATH  PROBE_SET_FILE
```

This will generate a `count_matrix_f100.txt` file that is the count of all probes in HyPR-seq dataset per cell barcode. 

## How to call variants from Tapestri data and merge with HyPR-seq data
Have your `.loom` file from the Tapestri pipeline available. 

Run the variant calling script using the following command. This will also merge the HyPR-seq and Tapestri outputs:

```
python merge_hyprseq_and_tapestri_output.py \
  -r /broad/thechenlab/stag_seq/test/MDM_PolyIC2_HyPR_S4_R1_001_work_dir/count_matrix_f100.txt \  
  -d /broad/thechenlab/stag_seq/test/matrix_v1_loom/MDM_PolyIC2.cells.loom \
  -g /broad/thechenlab/stag_seq/test/germline_mutations.txt \
  -o /broad/thechenlab/stag_seq/test/MDM_PolyIC2_HyPR_S4_R1_001_work_dir
```
