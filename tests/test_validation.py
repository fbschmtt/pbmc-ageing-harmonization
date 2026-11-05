import anndata as ad
import numpy as np
import pytest

from pbmc_pipeline.validation import ValidationError, validate_counts


def test_validate_counts_accepts_integer_valued_float_counts():
    validate_counts(ad.AnnData(X=np.array([[1.0, 2.0]])))


def test_validate_counts_rejects_any_non_integer_value():
    with pytest.raises(ValidationError, match="maximum deviation"):
        validate_counts(ad.AnnData(X=np.array([[1.009, 2.0]])))
