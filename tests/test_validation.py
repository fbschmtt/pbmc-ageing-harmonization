import anndata as ad
import numpy as np
import pytest

from pbmc_pipeline.validation import ValidationError, validate_counts


def test_validate_counts_accepts_small_floating_point_rounding_error():
    validate_counts(ad.AnnData(X=np.array([[1.009, 2.0]])))


def test_validate_counts_rejects_material_non_integer_values():
    with pytest.raises(ValidationError, match="maximum deviation"):
        validate_counts(ad.AnnData(X=np.array([[1.02, 2.0]])))
