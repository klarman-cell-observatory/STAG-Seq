FROM continuumio/miniconda3:latest

WORKDIR /STAG_Seq_RNA

RUN mkdir -p /ref && chmod 777 /ref
RUN mkdir -p /raw && chmod 777 /raw

# Pull the repo
RUN apk add --no-cache git
RUN git clone https://github.com/klarman-cell-observatory/STAG-Seq /STAG-Seq-main

# Copy required files
#COPY STAG-Seq-main/barcode_lookup.pkl /ref/
COPY STAG-Seq-main/scripts/environment.yml .
COPY STAG-Seq-main/scripts/scripts/chunk_script.sh .
COPY STAG-Seq-main/scripts/scripts/chunk_script_demultiplexing.sh .
COPY STAG-Seq-main/scripts/scripts/filter_and_group.py .
COPY STAG-Seq-main/scripts/scripts/filter_and_group_demultiplexing.py .
COPY STAG-Seq-main/scripts/ref/whitelist_missionbio.txt /ref/

RUN chmod +x chunk_script.sh
RUN chmod +x chunk_script_demultiplexing.sh
RUN chmod +x filter_and_group.py
RUN chmod +x filter_and_group_demultiplexing.py

RUN rm -rf STAG-Seq-main

# Install system dependencies
RUN apt-get update && apt-get install -y \
        curl \
        gnupg \
        build-essential \
        python3-dev && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
        | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update && apt-get install -y google-cloud-sdk && \
    rm -rf /var/lib/apt/lists/*

# Create conda environment from environment.yml
RUN conda env create -f environment.yml -n stag

# Install seqtk into the same environment via bioconda
RUN conda install -n stag -c bioconda seqtk -y

# Ensure conda environment is active by default
ENV CONDA_DEFAULT_ENV=stag
ENV PATH=/opt/conda/envs/stag/bin:$PATH

# Set working directory
WORKDIR /STAG_Seq_RNA

# Optional: verify installations
RUN /bin/bash -c "source activate stag && python --version && seqtk 2>&1 | head -n 5"

# Default command
CMD ["/bin/bash"]
