version 1.0

workflow stag_seq_rna {
    input {
        File input_fastq_probe_table
        String output_directory
        Boolean save_barcodes_table = true
        Boolean skip_barcode_extraction = false
        Boolean demultiplex_mode = false
        Int num_cpu = 8
        String memory = "16G"
        Int disk_space = 1000
        String docker_registry
    }

    parameter_meta {
        skip_barcode_extraction: "false = full run (FASTQ -> barcode extraction -> count matrix). true = reuse saved SAMPLE_output.txt (and SAMPLE_barcode_batch_mapping.csv when demultiplex_mode=true) from output_directory and run only the filter/group step. Requires an earlier run with save_barcodes_table=true."
        demultiplex_mode: "true = demultiplexing pipeline (column 3 of the input table must be filled). false = standard pipeline (column 3 must be empty)."
        save_barcodes_table: "Upload barcode tables to output_directory after extraction. Needed for later skip_barcode_extraction runs."
    }
    # Read all lines
    Array[String] all_rows = read_lines(input_fastq_probe_table)

    # Remove the header manually by scattering with an index check
    scatter (i in range(length(all_rows))) {
        String line = all_rows[i]

        if (i > 0) {
            call process_sample {
                input:
                    line = line,
                    save_barcodes_table = save_barcodes_table,
                    skip_barcode_extraction = skip_barcode_extraction,
                    demultiplex_mode = demultiplex_mode,
                    output_directory = output_directory,
                    docker_registry = docker_registry,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space = disk_space
            }
        }
    }

    output {
        # Array[File?] filtered_outputs = process_sample.result_h5
    }
}

task process_sample {
    input {
        String line
        Boolean save_barcodes_table
        Boolean skip_barcode_extraction
        Boolean demultiplex_mode
        String output_directory
        String docker_registry
        Int num_cpu
        String memory
        Int disk_space
    }

    command <<<
        set -euo pipefail
        chmod -R 777 /mnt/disks/cromwell_root/

        # Split CSV line manually using awk (since WDL 1.0 lacks split())
        FASTQ_PATH=$(echo "~{line}" | awk -F, '{print $1}')
        PROBE_PATH=$(echo "~{line}" | awk -F, '{print $2}')
        MULTIPLEX_PATH=$(echo "~{line}" | awk -F, '{print $3}')
        RAW_NAME=$(basename "$FASTQ_PATH")

        MNT_PATH="/mnt/disks/cromwell_root/"

        # Strip extensions and Illumina _S# suffix
        ROOT_NAME=$(echo "$RAW_NAME" | sed -E 's/(\.fastq(\.gz)?|\.fq(\.gz)?)?(_S[0-9]+.*)?$//')

        echo "Processing sample: $ROOT_NAME"
        mkdir -p "$MNT_PATH"/raw
        mkdir -p "$MNT_PATH"/out

        OUT_DIR="~{output_directory}"
        OUT_DIR="${OUT_DIR%/}/"

        echo "Sample:               $ROOT_NAME"
        echo "FASTQ:                $FASTQ_PATH"
        echo "Probe file:           $PROBE_PATH"
        echo "Multiplex file:       ${MULTIPLEX_PATH:-<none>}"
        echo "Output directory:     $OUT_DIR"
        echo "skip_barcode_extraction=~{skip_barcode_extraction}  demultiplex_mode=~{demultiplex_mode}  save_barcodes_table=~{save_barcodes_table}"

        check() {
            if ! gcloud storage ls "$1" > /dev/null 2>&1; then
                echo "ERROR: $2 not found or not accessible: $1" >&2
                exit 1
            fi
        }

        # Table vs. demultiplex flag
        if [ "~{demultiplex_mode}" == true ] && [ -z "$MULTIPLEX_PATH" ]; then
            echo "ERROR: demultiplex_mode=true but column 3 (multiplex file) is empty" >&2; exit 1
        fi
        if [ "~{demultiplex_mode}" == false ] && [ -n "$MULTIPLEX_PATH" ]; then
            echo "ERROR: demultiplex_mode=false but column 3 is filled ($MULTIPLEX_PATH)" >&2; exit 1
        fi

        # Files every run localizes
        check "$PROBE_PATH" "Probe file"
        if [ "~{demultiplex_mode}" == true ]; then
            check "$MULTIPLEX_PATH" "Multiplex file"
        fi

        if [ "~{skip_barcode_extraction}" == true ]; then
            # Partial run: saved tables must already exist
            check "${OUT_DIR}${ROOT_NAME}_output.txt" "Saved barcode table (run once with skip_barcode_extraction=false and save_barcodes_tables=true)"
            if [ "~{demultiplex_mode}" == true ]; then
                check "${OUT_DIR}${ROOT_NAME}_barcode_batch_mapping.csv" "Saved barcode batch mapping"
            fi
        else
            # Full run: the FASTQ is localized
            check "$FASTQ_PATH" "FASTQ"
        fi
        echo "Pre-flight checks passed"
 
        #All runs (partial or not) will pull the probe file
        echo "Localizing probe file"
        gcloud storage cp "$PROBE_PATH" "$MNT_PATH"/raw/

        #[PARTIAL]
        if [ "~{skip_barcode_extraction}" == true ]; then
            echo "[PARTIAL] Skipping barcode extraction"
            echo "[PARTIAL] Localizing existing output file from GCS"
            gcloud storage cp "${OUT_DIR}${ROOT_NAME}_output.txt" "$MNT_PATH"/out/

            if [ "~{demultiplex_mode}" == true ]; then
                echo "[PARTIAL] Localizing multiplex file"
                gcloud storage cp "$MULTIPLEX_PATH" "$MNT_PATH"/raw/

                echo "[PARTIAL] Localizing barcode batch mapping"
                   gcloud storage cp "${OUT_DIR}${ROOT_NAME}_barcode_batch_mapping.csv" "$MNT_PATH"/barcode_batch_mapping.csv

                echo "[PARTIAL] Running filter_and_group_demultiplexing.py"
                python /STAG_Seq_RNA/filter_and_group_demultiplexing.py \
                    -p "${MNT_PATH}/raw/$(basename $PROBE_PATH)" \
                    -f "${MNT_PATH}/out/${ROOT_NAME}_output.txt" \
                    -b /ref/barcode_lookup.pkl \
                    -m "${MNT_PATH}/barcode_batch_mapping.csv" \
                    -x "${MNT_PATH}/raw/$(basename $MULTIPLEX_PATH)"

            else
                echo "[PARTIAL] Running filter_and_group.py"   
                python /STAG_Seq_RNA/filter_and_group.py \
                    -p "${MNT_PATH}/raw/$(basename $PROBE_PATH)" \
                    -f "${MNT_PATH}/out/${ROOT_NAME}_output.txt" \
                    -b /ref/barcode_lookup.pkl
            

            fi
            echo "[PARTIAL] Uploading results to GCS"
            gcloud storage cp count_matrix_f100.h5ad "${OUT_DIR}${ROOT_NAME}_count_matrix.h5ad"
        #[FULL]
        else
            echo "[FULL] Localizing FASTQ from $FASTQ_PATH"
            gcloud storage cp "$FASTQ_PATH" "$MNT_PATH"/raw/

            echo "[FULL] Running chunk_script.sh"
            if [ "~{demultiplex_mode}" == true ]; then
                bash /STAG_Seq_RNA/chunk_script_demultiplexing.sh "$MNT_PATH/raw/$RAW_NAME" "$RAW_NAME" "$ROOT_NAME"
            else
                bash /STAG_Seq_RNA/chunk_script.sh "$MNT_PATH/raw/$RAW_NAME" "$RAW_NAME" "$ROOT_NAME"
            fi
            
            cat bc_extract_${ROOT_NAME}.log
            ls /mnt/disks/cromwell_root/

            echo "[FULL] Uploading intermediate output files (for future partial runs and troubleshooting)"
            if [ "~{save_barcodes_table}" == true ]; then
                gcloud storage cp "${MNT_PATH}${ROOT_NAME}_output.txt" "$OUT_DIR"
                if [ "~{demultiplex_mode}" == true ]; then
                    gcloud storage cp "${MNT_PATH}barcode_batch_mapping.csv" "${OUT_DIR}${ROOT_NAME}_barcode_batch_mapping.csv"
                fi
            fi
            
            mv ${MNT_PATH}${ROOT_NAME}_output.txt ${MNT_PATH}/out/

            if [ "~{demultiplex_mode}" == true ]; then
                echo "[FULL] Localizing multiplex file"
                gcloud storage cp "$MULTIPLEX_PATH" "$MNT_PATH"/raw/

                echo "[FULL] Running filter_and_group_demultiplexing.py"
                python /STAG_Seq_RNA/filter_and_group_demultiplexing.py \
                    -p "${MNT_PATH}/raw/$(basename $PROBE_PATH)" \
                    -f "${MNT_PATH}/out/${ROOT_NAME}_output.txt" \
                    -b /ref/barcode_lookup.pkl \
                    -m "${MNT_PATH}/barcode_batch_mapping.csv" \
                    -x "${MNT_PATH}/raw/$(basename $MULTIPLEX_PATH)"
            else
                echo "[FULL] Running filter_and_group.py"
                python /STAG_Seq_RNA/filter_and_group.py \
                    -p "${MNT_PATH}/raw/$(basename $PROBE_PATH)" \
                    -f "${MNT_PATH}/out/${ROOT_NAME}_output.txt" \
                    -b /ref/barcode_lookup.pkl
            
            fi
        echo "[FULL] Uploading results to GCS"
        gcloud storage cp count_matrix_f100.h5ad "${OUT_DIR}${ROOT_NAME}_count_matrix.h5ad"
        fi
    >>>

    output {
        # File result_h5 = "count_matrix_f100.h5ad"
    }

    runtime {
        docker: docker_registry
        cpu: num_cpu
        memory: memory
        bootDiskSizeGB: 25
        disks: "local-disk ${disk_space} HDD"
        }
    }
