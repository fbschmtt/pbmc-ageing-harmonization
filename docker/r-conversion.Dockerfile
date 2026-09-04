FROM satijalab/seurat@sha256:048d965487afc2961c450879692281a73e24f6de47360ed2f11b31112ebd844e

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libhdf5-dev \
        libssl-dev \
        libxml2-dev \
        python3 \
        python3-dev \
        python3-pip \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

ENV RETICULATE_PYTHON=/opt/venv/bin/python
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir anndata==0.12.2 h5py==3.14.0

WORKDIR /opt/pbmc-pipeline
COPY scripts/convert_rds.R /opt/pbmc-pipeline/scripts/convert_rds.R

CMD ["Rscript", "scripts/convert_rds.R", "--help"]
