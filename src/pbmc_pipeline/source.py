"""Read configured expression inputs, including 10x matrices held in tar archives."""
from __future__ import annotations

import gzip
import re
import tarfile
import tempfile
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


def read_study_obs(root: Path, path: Path, study: dict) -> pd.DataFrame:
    """Read cell metadata without materializing an H5AD expression matrix."""
    del root, study
    if path.suffix == ".tar":
        return _read_10x_tar_obs(path)
    import h5py

    with h5py.File(path, "r") as source:
        return ad.io.read_elem(source["obs"])


def read_study_counts_as_x(root: Path, path: Path, study: dict) -> ad.AnnData:
    """Read an H5AD with only its configured count representation as X.

    A temporary H5AD uses HDF5 external links for the original metadata and
    selected count matrix. AnnData then loads one expression matrix instead of
    loading X, raw, and every layer before selecting counts.
    """
    counts_source = study["counts_source"]
    if path.suffix == ".tar":
        return read_study_input(root, path, study)

    import h5py

    source_path = path.resolve()
    with h5py.File(source_path, "r") as source:
        if counts_source == "X":
            matrix_path, var_path, varm_path = "/X", "/var", "/varm"
            keep_source_annotations = True
        elif counts_source == "raw":
            if "raw" not in source or "X" not in source["raw"]:
                raise ValueError("Configured raw counts are absent")
            matrix_path, var_path = "/raw/X", "/raw/var"
            varm_path = "/raw/varm" if "varm" in source["raw"] else None
            keep_source_annotations = False
        elif counts_source.startswith("layer:"):
            layer = counts_source.removeprefix("layer:")
            if not layer or "layers" not in source or layer not in source["layers"]:
                raise ValueError(f"Configured count layer is absent: {layer!r}")
            matrix_path, var_path, varm_path = f"/layers/{layer}", "/var", "/varm"
            keep_source_annotations = True
        else:
            raise ValueError(f"Unsupported counts_source: {counts_source}")

    with tempfile.TemporaryDirectory(prefix="pbmc-count-view-") as temporary_dir:
        view_path = Path(temporary_dir) / "counts-as-x.h5ad"
        with h5py.File(source_path, "r") as source, h5py.File(view_path, "w") as view:
            for key, value in source.attrs.items():
                view.attrs[key] = value
            view["X"] = h5py.ExternalLink(str(source_path), matrix_path)
            view["var"] = h5py.ExternalLink(str(source_path), var_path)
            if varm_path is not None and varm_path.lstrip("/") in source:
                view["varm"] = h5py.ExternalLink(str(source_path), varm_path)
            for key in source:
                if key in {"X", "layers", "raw", "var", "varm"}:
                    continue
                if not keep_source_annotations and key != "obs":
                    continue
                view[key] = h5py.ExternalLink(str(source_path), f"/{key}")
        return ad.read_h5ad(view_path)


def _read_10x_tar_obs(path: Path) -> pd.DataFrame:
    """Read barcodes and source metadata from a 10x archive without its matrix."""
    frames = []
    with tarfile.open(path) as archive:
        for prefix, _feature_member in _archive_prefixes(archive):
            subject_id = _subject_id(prefix)
            barcodes = _read_tsv_gz(archive, f"{prefix}_barcodes.tsv.gz").iloc[:, 0].astype(str)
            cell_ids = pd.Index([f"{subject_id}__{barcode}" for barcode in barcodes], name="cell_id")
            frames.append(pd.DataFrame(
                {"subject_id": subject_id, "barcode": barcodes.to_numpy()}, index=cell_ids
            ))
    if not frames:
        raise ValueError(f"No complete 10x Matrix Market triplets found in {path}")
    return pd.concat(frames, axis=0)


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
