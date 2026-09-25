"""Read configured expression inputs, including 10x matrices held in tar archives."""
from __future__ import annotations

import gzip
import re
import tarfile
import warnings
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy.io import mmread


def read_study_input(root: Path, path: Path, study: dict) -> ad.AnnData:
    """Return an AnnData object without extracting an archive to disk."""
    if path.suffix != ".tar":
        return ad.read_h5ad(path)
    return _read_10x_tar(path)


def _read_10x_tar(path: Path) -> ad.AnnData:
    """Read every complete 10x triplet in a GEO tar archive in memory."""
    parts = []
    with tarfile.open(path) as archive:
        prefixes = _archive_prefixes(archive)
        for prefix, feature_member in prefixes:
            subject_id = _subject_id(prefix)
            barcodes = _read_tsv_gz(archive, f"{prefix}_barcodes.tsv.gz").iloc[:, 0].astype(str)
            features = _read_tsv_gz(archive, feature_member)
            matrix_member = archive.extractfile(f"{prefix}_matrix.mtx.gz")
            if matrix_member is None:
                raise ValueError(f"Archive member is not a file: {prefix}_matrix.mtx.gz")
            with matrix_member, gzip.GzipFile(fileobj=matrix_member) as compressed:
                matrix = mmread(compressed).tocsr().transpose().tocsr()
            var = pd.DataFrame(
                {
                    "gene_id": features.iloc[:, 0].astype(str).to_numpy(),
                    "gene_symbol": features.iloc[:, 1].astype(str).to_numpy(),
                    "feature_type": (
                        features.iloc[:, 2].astype(str).to_numpy()
                        if features.shape[1] > 2
                        else ["Gene Expression"] * len(features)
                    ),
                },
                index=pd.Index(features.iloc[:, 0].astype(str), name="gene_id"),
            )
            cell_ids = pd.Index([f"{subject_id}__{barcode}" for barcode in barcodes], name="cell_id")
            obs = pd.DataFrame({"subject_id": subject_id, "barcode": barcodes.to_numpy()}, index=cell_ids)
            parts.append(ad.AnnData(X=matrix, obs=obs, var=var.copy()))
    if not parts:
        raise ValueError(f"No complete 10x Matrix Market triplets found in {path}")
    result = ad.concat(parts, axis=0, join="outer", merge="first", index_unique=None)
    # Preserve original IDs in ``var["gene_id"]`` but avoid naming the feature
    # index ``gene_id``: outer concatenation can make its values diverge from
    # that column before the generated test fixture is written.
    result.var_names.name = None
    return result


def _archive_prefixes(archive: tarfile.TarFile) -> list[tuple[str, str]]:
    """Find complete GEX triplets and warn about non-GEX/incomplete candidates."""
    barcode_suffix = "_barcodes.tsv.gz"
    matrix_suffix = "_matrix.mtx.gz"
    feature_suffixes = ("_features.tsv.gz", "_genes.tsv.gz")
    names = {member.name for member in archive.getmembers() if member.isfile()}
    prefixes = sorted(name.removesuffix(barcode_suffix) for name in names if name.endswith(barcode_suffix))
    complete = []
    skipped = []
    for prefix in prefixes:
        feature_member = next((f"{prefix}{suffix}" for suffix in feature_suffixes if f"{prefix}{suffix}" in names), None)
        if f"{prefix}{matrix_suffix}" not in names or feature_member is None:
            skipped.append(prefix)
        else:
            complete.append((prefix, feature_member))
    if skipped:
        warnings.warn(
            f"Skipping {len(skipped)} incomplete or non-GEX 10x candidates: {', '.join(skipped)}",
            UserWarning,
            stacklevel=2,
        )
    return complete


def _subject_id(prefix: str) -> str:
    match = re.search(r"_(JB\d+)$", Path(prefix).name)
    if match is None:
        raise ValueError(f"Cannot derive Nehar-Belaid subject ID from archive member prefix: {prefix}")
    return match.group(1)


def _read_tsv_gz(archive: tarfile.TarFile, member: str) -> pd.DataFrame:
    source = archive.extractfile(member)
    if source is None:
        raise ValueError(f"Archive member is not a file: {member}")
    with source, gzip.GzipFile(fileobj=source) as compressed:
        return pd.read_csv(compressed, sep="\t", header=None)
