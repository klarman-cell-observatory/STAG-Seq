version 1.0

workflow stag_seq_rna {
    input {
        File input_fastq_probe_table
        String output_directory
        Boolean copy_intermediates
        Boolean run_filter_group_only
        Int num_cpu = 8
        String memory = "16G"
        Int disk_space = 1000
        String docker_registry
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
                    copy_intermediates = copy_intermediates,
                    run_filter_group_only = run_filter_group_only,
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
        Boolean copy_intermediates
        Boolean run_filter_group_only
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
        RAW_NAME=$(basename "$FASTQ_PATH")

        MNT_PATH="/mnt/disks/cromwell_root/"

        # Strip extensions and Illumina _S# suffix
        ROOT_NAME=$(echo "$RAW_NAME" | sed -E 's/(\.fastq(\.gz)?|\.fq(\.gz)?)?(_S[0-9]+.*)?$//')

        echo "Processing sample: $ROOT_NAME"
        mkdir -p "$MNT_PATH"/raw
        mkdir -p "$MNT_PATH"/out

        if [ "~{run_filter_group_only}" == "true" ]; then
            echo "Copying probe file"
            gcloud storage cp "$PROBE_PATH" "$MNT_PATH"/raw/
            
            echo "Retrieving existing bc file"
            gcloud storage cp ~{output_directory}"$ROOT_NAME"_output.txt "$MNT_PATH"/out/

            echo "Running filter_and_group.py"
            python /STAG_Seq_RNA/filter_and_group.py \
                -p "${MNT_PATH}/raw/$(basename $PROBE_PATH)" \
                -f "${MNT_PATH}/out/${ROOT_NAME}_output.txt" \
                -b /ref/barcode_lookup.pkl

            echo "Uploading results to GCS"
            gcloud storage cp count_matrix_f100.h5ad "~{output_directory}${ROOT_NAME}_count_matrix.h5ad"
        else
            echo "Copying FASTQ from $FASTQ_PATH"
            gcloud storage cp "$FASTQ_PATH" "$MNT_PATH"/raw/

            echo "Running chunk_script.sh"
            bash /STAG_Seq_RNA/chunk_script.sh "$MNT_PATH/raw/$RAW_NAME" "$RAW_NAME" "$ROOT_NAME"
            
            cat bc_extract_${ROOM_NAME}.log
            ls /mnt/disks/cromwell_root/

            if [ "~{copy_intermediates}" == "true" ]; then
                echo "Copying intermediate file output from chunking script"
                gcloud storage cp "$MNT_PATH$ROOT_NAME_output.txt" ~{output_directory}
            fi
            
            mv ${MNT_PATH}${ROOT_NAME}_output.txt ${MNT_PATH}/out/

            echo "Copying probe file"
            gcloud storage cp "$PROBE_PATH" "$MNT_PATH"/raw/

            echo "Running filter_and_group.py"
            python /STAG_Seq_RNA/filter_and_group.py \
                -p "${MNT_PATH}/raw/$(basename $PROBE_PATH)" \
                -f "${MNT_PATH}/out/${ROOT_NAME}_output.txt" \
                -b /ref/barcode_lookup.pkl

            echo "Uploading results to GCS"
            gcloud storage cp count_matrix_f100.h5ad "~{output_directory}${ROOT_NAME}_count_matrix.h5ad"
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
