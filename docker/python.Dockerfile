FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLCONFIGDIR=/tmp/matplotlib \
    CELLTYPIST_FOLDER=/tmp/celltypist

WORKDIR /opt/pbmc-pipeline

COPY pyproject.toml README.md ./
COPY src ./src
RUN python -m pip install --no-cache-dir '.[qc]'

COPY config ./config
COPY scripts ./scripts
COPY aifi_models ./aifi_models
COPY reports ./reports

CMD ["pbmc-harmonize", "--help"]
