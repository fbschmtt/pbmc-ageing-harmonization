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
