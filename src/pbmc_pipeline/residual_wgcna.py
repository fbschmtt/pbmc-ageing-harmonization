"""WGCNA utilities for per-cell-type pseudobulk Pearson residuals."""
from __future__ import annotations

import contextlib
import io
from dataclasses import dataclass
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, pointbiserialr
from statsmodels.stats.multitest import multipletests

from .config import ResidualWGCNASettings
from .covariates import normalize_cmv_status, normalize_smoking_status


@dataclass(frozen=True)
class ResidualWGCNAResult:
    """Published tables and status for one residual-network analysis."""

    metadata: dict[str, Any]
    soft_threshold: pd.DataFrame
    module_assignments: pd.DataFrame
    eigengenes: pd.DataFrame
    module_summary: pd.DataFrame


def _empty_result(*, settings: ResidualWGCNASettings, reason: str, n_samples: int, n_genes: int):
    return ResidualWGCNAResult(
        metadata={
            "status": "skipped",
            "reason": reason,
            "settings": settings.to_mapping(),
            "n_samples": n_samples,
            "n_genes_input": n_genes,
            "n_genes_selected": 0,
        },
        soft_threshold=pd.DataFrame(),
        module_assignments=pd.DataFrame(),
        eigengenes=pd.DataFrame(),
        module_summary=pd.DataFrame(),
    )


def _associated_genes(results: pd.DataFrame | None, *, fdr_threshold: float) -> set[str]:
    if results is None or not {"gene", "padj"}.issubset(results):
        return set()
    padj = pd.to_numeric(results["padj"], errors="coerce")
    return set(results.loc[padj.le(fdr_threshold), "gene"].astype(str))


def _module_trait_associations(
    eigengenes: pd.DataFrame, sample_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Relate eigengenes to observed yes/no CMV and smoking statuses."""
    associations = pd.DataFrame({
        "module": [column.removeprefix("ME") for column in eigengenes.columns]
    })
    for covariate, normalizer in (
        ("cmv", normalize_cmv_status),
        ("smoking_status", normalize_smoking_status),
    ):
        suffix = "smoking" if covariate == "smoking_status" else "cmv"
        columns = [
            f"{suffix}_eigengene_correlation", f"{suffix}_eigengene_pvalue",
            f"{suffix}_eigengene_fdr", f"n_{suffix}_samples",
            f"n_{suffix}_yes", f"n_{suffix}_no",
        ]
        if covariate not in sample_metadata or eigengenes.empty:
            for column in columns:
                associations[column] = np.nan
            continue
        trait = normalizer(sample_metadata[covariate]).map({"no": 0, "yes": 1})
        trait = pd.to_numeric(trait, errors="coerce")
        rows = []
        for column in eigengenes.columns:
            values = pd.to_numeric(eigengenes[column], errors="coerce")
            valid = values.notna() & trait.notna()
            observed_trait = trait.loc[valid]
            if len(observed_trait) < 3 or observed_trait.nunique() < 2:
                correlation, pvalue = np.nan, np.nan
            else:
                correlation, pvalue = pointbiserialr(observed_trait, values.loc[valid])
            rows.append({
                "module": column.removeprefix("ME"),
                f"{suffix}_eigengene_correlation": correlation,
                f"{suffix}_eigengene_pvalue": pvalue,
                f"n_{suffix}_samples": int(valid.sum()),
                f"n_{suffix}_yes": int(observed_trait.eq(1).sum()),
                f"n_{suffix}_no": int(observed_trait.eq(0).sum()),
            })
        trait_result = pd.DataFrame(rows)
        fdr_column = f"{suffix}_eigengene_fdr"
        trait_result[fdr_column] = np.nan
        valid_pvalues = trait_result[f"{suffix}_eigengene_pvalue"].notna()
        if valid_pvalues.any():
            trait_result.loc[valid_pvalues, fdr_column] = multipletests(
                trait_result.loc[valid_pvalues, f"{suffix}_eigengene_pvalue"], method="fdr_bh"
            )[1]
        associations = associations.merge(
            trait_result[["module", *columns]], on="module", how="left", validate="one_to_one"
        )
    return associations


def run_residual_wgcna(
    residuals: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    *,
    settings: ResidualWGCNASettings,
    cmv_results: pd.DataFrame | None = None,
    smoking_results: pd.DataFrame | None = None,
) -> ResidualWGCNAResult:
    """Run PyWGCNA on one type's age-bin-model Pearson residuals.

    The residuals already remove the fitted age-bin, sex, study-site, and depth
    effects. WGCNA is descriptive here: optional CMV and smoking comparisons
    are downstream checks against traits omitted from that residual model.
    """
    numeric = residuals.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    numeric = numeric.dropna(axis=0, how="any")
    numeric = numeric.loc[:, numeric.var(axis=0).gt(0)]
    if len(numeric) < settings.minimum_samples:
        return _empty_result(
            settings=settings,
            reason=(f"{len(numeric)} residual samples; at least {settings.minimum_samples} are required"),
            n_samples=len(numeric), n_genes=residuals.shape[1],
        )
    if numeric.shape[1] < settings.minimum_module_size:
        return _empty_result(
            settings=settings,
            reason=(
                f"{numeric.shape[1]} variable residual genes; at least "
                f"{settings.minimum_module_size} are required"
            ),
            n_samples=len(numeric), n_genes=residuals.shape[1],
        )

    variance = numeric.var(axis=0).sort_values(ascending=False)
    selected_genes = variance.head(settings.maximum_genes).index
    selected = numeric.loc[:, selected_genes]
    metadata = sample_metadata.copy()
    if "sample_key" in metadata:
        metadata = metadata.drop_duplicates("sample_key").set_index("sample_key", drop=False)
    metadata = metadata.reindex(selected.index).copy()
    cmv_genes = _associated_genes(cmv_results, fdr_threshold=settings.cmv_fdr_threshold)
    smoking_genes = _associated_genes(
        smoking_results, fdr_threshold=settings.cmv_fdr_threshold,
    )

    try:
        from PyWGCNA import WGCNA

        expression = ad.AnnData(
            X=selected.to_numpy(dtype=float),
            obs=metadata,
            var=pd.DataFrame(index=selected.columns.astype(str)),
        )
        network = WGCNA(
            name="residual_wgcna",
            anndata=expression,
            powers=settings.soft_threshold_powers,
            RsquaredCut=settings.scale_free_fit_threshold,
            MeanCut=settings.mean_connectivity_threshold,
            networkType="signed",
            TOMType="signed",
            minModuleSize=settings.minimum_module_size,
            MEDissThres=settings.module_eigengene_dissimilarity_threshold,
            save=False,
        )
        # PyWGCNA's preprocessing assumes non-negative TPM values. Pearson
        # residuals are centered signed values, so inputs have already been
        # cleaned above and we run its network/module step directly.
        existing_figures = set(plt.get_fignums())
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                network.findModules({
                    "cutreeHybrid": {
                        "deepSplit": settings.dynamic_tree_deep_split,
                        "pamRespectsDendro": False,
                    }
                })
        finally:
            for figure_number in set(plt.get_fignums()) - existing_figures:
                plt.close(figure_number)
    except (Exception, SystemExit) as exc:  # noqa: BLE001 - preserve reports on WGCNA failure
        return ResidualWGCNAResult(
            metadata={
                "status": "failed",
                "reason": f"{type(exc).__name__}: {exc}",
                "settings": settings.to_mapping(),
                "n_samples": len(selected),
                "n_genes_input": residuals.shape[1],
                "n_genes_selected": selected.shape[1],
            },
            soft_threshold=pd.DataFrame(),
            module_assignments=pd.DataFrame(),
            eigengenes=pd.DataFrame(),
            module_summary=pd.DataFrame(),
        )

    soft_threshold = network.sft.copy()
    module_assignments = network.datExpr.var[["moduleColors", "moduleLabels"]].copy()
    module_assignments.index.name = "gene"
    module_assignments = module_assignments.reset_index().rename(columns={
        "moduleColors": "module", "moduleLabels": "module_label",
    })
    module_assignments["residual_variance"] = module_assignments["gene"].map(variance)
    module_assignments["is_cmv_associated"] = module_assignments["gene"].isin(cmv_genes)
    module_assignments["is_smoking_associated"] = module_assignments["gene"].isin(smoking_genes)
    eigengenes = network.MEs.copy()
    eigengenes.index = selected.index.astype(str)
    eigengenes.index.name = "sample_key"

    trait_associations = _module_trait_associations(eigengenes, metadata)
    network_genes = set(module_assignments["gene"])
    module_rows = []
    for module, group in module_assignments.groupby("module", sort=True, observed=True):
        overlap = set(group["gene"]) & cmv_genes
        smoking_overlap = set(group["gene"]) & smoking_genes
        n_background_cmv = len(network_genes & cmv_genes)
        n_background_smoking = len(network_genes & smoking_genes)
        if n_background_cmv:
            enrichment_pvalue = hypergeom.sf(
                len(overlap) - 1, len(network_genes), n_background_cmv, len(group)
            )
        else:
            enrichment_pvalue = np.nan
        if n_background_smoking:
            smoking_enrichment_pvalue = hypergeom.sf(
                len(smoking_overlap) - 1, len(network_genes), n_background_smoking, len(group)
            )
        else:
            smoking_enrichment_pvalue = np.nan
        module_rows.append({
            "module": module,
            "n_genes": len(group),
            "n_cmv_associated_genes": len(overlap),
            "cmv_gene_overlap_pvalue": enrichment_pvalue,
            "n_smoking_associated_genes": len(smoking_overlap),
            "smoking_gene_overlap_pvalue": smoking_enrichment_pvalue,
        })
    module_summary = pd.DataFrame(module_rows)
    module_summary["cmv_gene_overlap_fdr"] = np.nan
    valid_enrichment = module_summary["cmv_gene_overlap_pvalue"].notna()
    if valid_enrichment.any():
        module_summary.loc[valid_enrichment, "cmv_gene_overlap_fdr"] = multipletests(
            module_summary.loc[valid_enrichment, "cmv_gene_overlap_pvalue"], method="fdr_bh"
        )[1]
    module_summary["smoking_gene_overlap_fdr"] = np.nan
    valid_smoking_enrichment = module_summary["smoking_gene_overlap_pvalue"].notna()
    if valid_smoking_enrichment.any():
        module_summary.loc[valid_smoking_enrichment, "smoking_gene_overlap_fdr"] = multipletests(
            module_summary.loc[valid_smoking_enrichment, "smoking_gene_overlap_pvalue"],
            method="fdr_bh",
        )[1]
    module_summary = module_summary.merge(trait_associations, on="module", how="left")
    module_summary = module_summary.sort_values(
        ["n_cmv_associated_genes", "n_genes", "module"], ascending=[False, False, True]
    ).reset_index(drop=True)
    return ResidualWGCNAResult(
        metadata={
            "status": "complete",
            "reason": "",
            "settings": settings.to_mapping(),
            "n_samples": len(selected),
            "n_genes_input": residuals.shape[1],
            "n_genes_variable": numeric.shape[1],
            "n_genes_selected": selected.shape[1],
            "n_modules": int(module_assignments["module"].nunique()),
            "selected_soft_threshold_power": int(network.power),
            "n_cmv_associated_genes_in_network": len(network_genes & cmv_genes),
            "n_smoking_associated_genes_in_network": len(network_genes & smoking_genes),
            "cmv_de_fdr_threshold": settings.cmv_fdr_threshold,
            "network_type": "signed",
            "tom_type": "signed",
            "implementation": "PyWGCNA",
        },
        soft_threshold=soft_threshold,
        module_assignments=module_assignments,
        eigengenes=eigengenes.reset_index(),
        module_summary=module_summary,
    )


def trait_associated_module_scores(
    result: ResidualWGCNAResult, sample_metadata: pd.DataFrame, *, covariate: str,
) -> pd.DataFrame:
    """Score the module most associated with CMV or smoking status."""
    if covariate not in {"cmv", "smoking_status"}:
        raise ValueError(f"Unsupported module-score covariate: {covariate}")
    suffix = "smoking" if covariate == "smoking_status" else "cmv"
    normalizer = normalize_smoking_status if suffix == "smoking" else normalize_cmv_status
    status_column = "trait_status"
    correlation_column = f"{suffix}_eigengene_correlation"
    columns = ["sample_key", "module", "module_score", status_column, "trait_correlation"]
    summary = result.module_summary.dropna(subset=[correlation_column])
    summary = summary.loc[summary["module"].astype(str).str.lower().ne("grey")]
    if summary.empty or result.eigengenes.empty or covariate not in sample_metadata:
        return pd.DataFrame(columns=columns)
    selected = summary.loc[summary[correlation_column].abs().idxmax()]
    module = str(selected["module"])
    eigengene_column = f"ME{module}"
    if eigengene_column not in result.eigengenes:
        return pd.DataFrame(columns=columns)
    correlation = float(selected[correlation_column])
    metadata = sample_metadata.drop_duplicates("sample_key").copy()
    metadata[status_column] = normalizer(metadata[covariate])
    metadata[status_column] = metadata[status_column].where(
        metadata[status_column].isin({"yes", "no"}), "not provided"
    )
    scores = result.eigengenes[["sample_key", eigengene_column]].rename(
        columns={eigengene_column: "module_score"}
    )
    scores["module_score"] = pd.to_numeric(scores["module_score"], errors="coerce")
    if correlation < 0:
        scores["module_score"] = -scores["module_score"]
    scores["module"] = module
    scores["trait_correlation"] = correlation
    return scores.merge(
        metadata[["sample_key", status_column]],
        on="sample_key", how="left", validate="one_to_one",
    )[columns]
