FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLCONFIGDIR=/tmp/matplotlib \
    CELLTYPIST_FOLDER=/tmp/celltypist

WORKDIR /opt/pbmc-pipeline

COPY pyproject.toml README.md ./
# Install dependencies before application source so source-only changes reuse this layer.
# A minimal placeholder package lets pip resolve the project dependencies without copying src.
RUN mkdir -p src/pbmc_pipeline \
    && touch src/pbmc_pipeline/__init__.py \
    && python -m pip install --no-cache-dir '.[qc]'

COPY src ./src
# Remove the placeholder build output and metadata, then install the real
# package without re-resolving its already installed dependencies.
RUN rm -rf build src/*.egg-info \
    && python -m pip install --no-cache-dir --no-deps '.[qc]'

COPY config ./config
COPY scripts ./scripts
COPY aifi_models ./aifi_models
COPY reports ./reports

CMD ["pbmc-harmonize", "--help"]
