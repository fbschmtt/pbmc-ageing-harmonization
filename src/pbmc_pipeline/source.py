"""Read configured expression inputs, including 10x matrices held in tar archives."""
from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy.io import mmread


def read_study_input(root: Path, path: Path, study: dict) -> ad.AnnData:
    """Return an AnnData object without extracting an archive to disk."""
    if path.suffix != ".tar":
        return ad.read_h5ad(path)
    return _read_10x_tar(root, path, study)


def _read_10x_tar(root: Path, path: Path, study: dict) -> ad.AnnData:
    samples = study.get("raw_samples")
    if not isinstance(samples, list) or not samples:
        raise ValueError("10x tar input requires a non-empty raw_samples configuration")
    parts = []
    reference_var = None
    with tarfile.open(path) as archive:
        for sample in samples:
            prefix = sample["archive_prefix"]
            barcodes = _read_tsv_gz(archive, f"{prefix}_barcodes.tsv.gz").iloc[:, 0].astype(str)
            features = _read_tsv_gz(archive, f"{prefix}_features.tsv.gz")
            matrix_member = archive.extractfile(f"{prefix}_matrix.mtx.gz")
            if matrix_member is None:
                raise ValueError(f"Archive member is not a file: {prefix}_matrix.mtx.gz")
            with matrix_member, gzip.GzipFile(fileobj=matrix_member) as compressed:
                matrix = mmread(compressed).tocsr().transpose().tocsr()
            var = pd.DataFrame(
                {"gene_id": features.iloc[:, 0].astype(str).to_numpy(), "gene_symbol": features.iloc[:, 1].astype(str).to_numpy(), "feature_type": features.iloc[:, 2].astype(str).to_numpy()},
                index=pd.Index(features.iloc[:, 0].astype(str), name="gene_id"),
            )
            if reference_var is None:
                reference_var = var
            elif not reference_var.index.equals(var.index):
                raise ValueError(f"{prefix}: feature order differs from the other archive members")
            cell_ids = pd.Index([f"{sample['subject_id']}__{barcode}" for barcode in barcodes], name="cell_id")
            obs = pd.DataFrame({"sample_id": sample["sample_id"], "subject_id": sample["subject_id"]}, index=cell_ids)
            parts.append(ad.AnnData(X=matrix, obs=obs, var=var.copy()))
    result = ad.concat(parts, axis=0, join="inner", merge="same", index_unique=None)
    _attach_matching_h5ad_labels(root, result, study)
    return result


def _read_tsv_gz(archive: tarfile.TarFile, member: str) -> pd.DataFrame:
    source = archive.extractfile(member)
    if source is None:
        raise ValueError(f"Archive member is not a file: {member}")
    with source, gzip.GzipFile(fileobj=source) as compressed:
        return pd.read_csv(compressed, sep="\t", header=None)


def _attach_matching_h5ad_labels(root: Path, adata: ad.AnnData, study: dict) -> None:
    """Attach published labels only where barcode identity can be established."""
    label_path = study.get("label_h5ad")
    label_columns = study.get("label_columns", [])
    if not label_path or not label_columns:
        return
    source_path = root / label_path
    if not source_path.exists():
        raise FileNotFoundError(f"Nehar-Belaid label H5AD does not exist: {source_path}")
    labels = ad.read_h5ad(source_path, backed="r")
    try:
        table = labels.obs.loc[:, [column for column in label_columns if column in labels.obs]].copy()
    finally:
        labels.file.close()
    table["_barcode"] = table.index.astype(str).str.extract(r"([ACGT]+-\d+)", expand=False)
    if "sample_id" not in labels.obs:
        return
    table["_key"] = labels.obs["sample_id"].astype(str) + "\0" + table["_barcode"]
    table = table.dropna(subset=["_barcode"]).drop_duplicates("_key", keep=False).set_index("_key")
    raw_barcodes = adata.obs_names.astype(str).str.rsplit("__", n=1).str[-1]
    raw_keys = adata.obs["sample_id"].astype(str) + "\0" + raw_barcodes
    for column in label_columns:
        values = table[column] if column in table else pd.Series(dtype="string")
        transferred = raw_keys.map(values).astype("string").fillna("not_provided")
        adata.obs[column] = transferred.to_numpy()
