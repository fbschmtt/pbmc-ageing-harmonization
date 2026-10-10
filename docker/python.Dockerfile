FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLCONFIGDIR=/tmp/matplotlib \
    CELLTYPIST_FOLDER=/tmp/celltypist

WORKDIR /opt/pbmc-pipeline

COPY pyproject.toml ./
COPY requirements.lock ./
# Resolve the application environment once from the reviewed lock before copying
# source. Source-only changes then reuse this layer without selecting new wheels.
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && python -m pip install --no-cache-dir -r requirements.lock \
    && apt-get purge -y --auto-remove build-essential \
    && rm -rf /var/lib/apt/lists/*

COPY src ./src
# Install the project itself without re-resolving locked dependencies.
RUN python -m pip install --no-cache-dir --no-deps --no-build-isolation '.[qc]'

COPY config ./config

CMD ["pbmc-harmonize", "--help"]
