"""Importable feature transformers persisted inside age-clock model bundles."""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.stats import rankdata
from sklearn.base import BaseEstimator, TransformerMixin


class WithinSampleRankTransformer(TransformerMixin, BaseEstimator):
    """Replace expression magnitudes with within-sample gene percentile ranks."""

    def fit(self, x, y=None):
        shape = x.shape if hasattr(x, "shape") else np.asarray(x).shape
        if len(shape) != 2:
            raise ValueError("Rank input must be a sample-by-gene matrix")
        self.n_features_in_ = shape[1]
        return self

    def transform(self, x):
        values = x.toarray() if sparse.issparse(x) else np.asarray(x)
        if values.ndim != 2 or values.shape[1] != self.n_features_in_:
            raise ValueError("Rank input does not match the fitted gene universe")
        if values.shape[1] < 2:
            return np.zeros(values.shape, dtype=float)
        ranks = rankdata(values, axis=1, method="average")
        return (ranks - 1.0) / (values.shape[1] - 1.0)


def within_study_gene_zscore(x, studies) -> tuple[np.ndarray, list[dict]]:
    """Z-score each gene across samples separately within each study.

    Study-specific means and standard deviations are computed from the samples
    being transformed. This is intended for cohort-level held-out study
    prediction, where the target study contributes an unlabeled expression cohort.
    """
    values = x.toarray() if sparse.issparse(x) else np.asarray(x, dtype=float)
    labels = np.asarray(studies, dtype=object).astype(str)
    if values.ndim != 2 or len(labels) != values.shape[0]:
        raise ValueError("Study z-score input must align samples with study labels")
    if not np.isfinite(values).all():
        raise ValueError("Study z-score input must contain finite values")
    transformed = np.empty_like(values, dtype=float)
    diagnostics = []
    for study in sorted(np.unique(labels)):
        mask = labels == study
        study_values = values[mask]
        means = study_values.mean(axis=0)
        scales = study_values.std(axis=0, ddof=0)
        constant = scales <= np.finfo(float).eps
        scales[constant] = 1.0
        transformed[mask] = (study_values - means) / scales
        diagnostics.append({
            "study": study,
            "n_samples": int(mask.sum()),
            "n_constant_genes": int(constant.sum()),
        })
    return transformed, diagnostics
