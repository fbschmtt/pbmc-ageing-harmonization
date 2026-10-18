"""Shared feature selection and sparse PCA preparation for embeddings."""
from __future__ import annotations

import logging
import time

import numpy as np

LOGGER = logging.getLogger(__name__)


def prepare_sparse_pca(adata, pipeline: dict, study_id: str):
    """Select non-V(D)J HVGs, scale sparsely, and compute configured PCs.

    Harmonization and per-cell-type analysis use the same mean-dispersion HVG
    selection, V(D)J exclusion, sparse scaling, and implicitly centered ARPACK
    PCA. The returned boolean masks are aligned to ``adata.var_names``.
    """
    import scanpy as sc

    started = time.perf_counter()
    hvg_settings = dict(pipeline["processing"]["hvg"])
    # Scanpy defaults max_mean to 3.0; pass infinity to disable that upper
    # expression cutoff while retaining the configured lower mean/dispersion
    # thresholds.
    hvg_settings["max_mean"] = np.inf
    sc.pp.highly_variable_genes(adata, **hvg_settings)
    excluded = adata.var_names.to_series().str.contains(
        pipeline["processing"]["exclude_vdj_regex"], case=False, regex=True
    ).to_numpy()
    selected = adata.var["highly_variable"].to_numpy(dtype=bool) & ~excluded
    if selected.sum() < 2:
        # Keep the embedding defined for studies with fewer than two passing
        # mean-dispersion HVGs while applying the same fallback in each stage.
        selected = ~excluded
    if selected.sum() < 2:
        raise ValueError(f"{study_id}: fewer than two non-V(D)J genes are available for PCA")
    adata.var["highly_variable"] = selected

    working = adata[:, selected].copy()
    working.X = working.X.astype(np.float32, copy=False)
    LOGGER.info(
        "study=%s step=select_embedding_features duration_seconds=%.2f input_features=%d selected_features=%d",
        study_id,
        time.perf_counter() - started,
        adata.n_vars,
        int(selected.sum()),
    )

    started = time.perf_counter()
    sc.pp.scale(working, zero_center=False)
    LOGGER.info(
        "study=%s step=scale_embedding duration_seconds=%.2f",
        study_id,
        time.perf_counter() - started,
    )

    neighbors = pipeline["processing"]["neighbors"]
    n_comps = min(int(neighbors["n_pcs"]), working.n_obs - 1, working.n_vars - 1)
    if n_comps < 2:
        raise ValueError(f"{study_id}: too few cells or genes for PCA")
    started = time.perf_counter()
    sc.tl.pca(
        working,
        n_comps=n_comps,
        zero_center=True,
        svd_solver="arpack",
        random_state=pipeline["processing"]["random_seed"],
    )
    LOGGER.info(
        "study=%s step=pca duration_seconds=%.2f components=%d",
        study_id,
        time.perf_counter() - started,
        int(working.obsm["X_pca"].shape[1]),
    )
    return working, selected, excluded, n_comps
