# %%
import base64
import json
import os
import re
import warnings
from html import escape
from io import BytesIO
from pathlib import Path
from types import SimpleNamespace

# The report runs in a notebook kernel where ipywidgets is intentionally not a
# runtime dependency; tqdm's optional-progress warning is not analytically
# relevant.
warnings.filterwarnings("ignore", message="IProgress not found.*")

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from IPython.display import HTML, Markdown, display
from matplotlib.colors import Normalize
from scipy.stats import f_oneway, pointbiserialr, spearmanr
from statsmodels.stats.multitest import multipletests

from pbmc_pipeline.report_theme import render_report_header

sns.set_theme(style="whitegrid")
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 110})
from pbmc_pipeline.config import AgeTrajectorySettings, ResidualWGCNASettings, read_json
from pbmc_pipeline.covariates import (
    categorical_contrast_label,
    categorical_levels,
    normalize_cmv_status,
    normalize_smoking_status,
)
from pbmc_pipeline.go_enrichment import (
    GeneSetEnrichmentRun,
    leading_edge_gene_preview,
    run_preranked_gsea,
    save_go_specificity_outputs,
    signed_de_gene_ranking,
)
from pbmc_pipeline.reporting import (
    plot_trajectory_merge_diagnostic,
    sample_cell_type_fractions,
    signed_value_at_largest_absolute_magnitude,
)
from pbmc_pipeline.residual_wgcna import run_residual_wgcna, trait_associated_module_scores
from pbmc_pipeline.trajectory_analysis import (
    cluster_sample_residuals,
    cluster_trajectory_profiles,
    summarize_residual_cluster_markers,
    summarize_residual_clustering,
    trajectory_merge_diagnostics,
)

# %%
de_dir = Path(os.environ["DE_REPORT_DIR"])
output_dir = Path(os.environ["DE_REPORT_OUTPUT_DIR"])
output_dir.mkdir(parents=True, exist_ok=True)
pipeline = read_json(Path(os.environ["DE_REPORT_CONFIG"]))
de_metadata_path = de_dir / "run_metadata.json"
de_run_metadata = read_json(de_metadata_path) if de_metadata_path.is_file() else {}
result_record_path = de_dir / "cell_type_result.json"
result_record = read_json(result_record_path) if result_record_path.is_file() else {}
cell_type_name = str(
    de_run_metadata.get("cell_type") or result_record.get("cell_type")
    or de_dir.parent.name.replace("-", " ").title()
)
abundance_path = de_dir / "abundance_input.tsv"
if not abundance_path.is_file():
    raise FileNotFoundError(f"{abundance_path}: rerun the per-type DE task to create abundance inputs")
adata = SimpleNamespace(obs=pd.read_csv(abundance_path, sep="	"))
parent_values = adata.obs["aifi_l1_parent_for_l2"].unique()
if len(parent_values) != 1:
    raise ValueError(f"{cell_type_name!r} maps to multiple configured L1 parents")
parent = str(parent_values[0])
trajectory_settings = AgeTrajectorySettings.from_mapping(
    pipeline["differential_expression"]["age_trajectory"]
)
completed_models = [
    model for model in de_run_metadata.get("models", [])
    if model.get("status") == "complete"
]
summary_metrics = [
    ("Pseudobulk samples", str(result_record.get("n_pseudobulks", "—"))),
    ("Completed models", str(len(completed_models))),
    ("Analysis mode", str(de_run_metadata.get("analysis_mode", "unknown")).replace("_", " ").title()),
]
display(HTML(render_report_header(
    title=cell_type_name,
    eyebrow="PBMC ageing · pseudobulk differential expression",
    subtitle=(
        "Sample-level pseudobulk differential-expression models, age-bin trajectories, "
        "and residual structure for this cell type."
    ),
    metrics=summary_metrics,
)))
if de_run_metadata.get("analysis_mode") == "test_only":
    display(Markdown(
        "> **TEST OUTPUT — NOT FOR BIOLOGICAL INTERPRETATION.** "
        + de_run_metadata["interpretation_warning"]
    ))

def display_collapsible_table(data, summary, *, index=False, max_height=None):
    """Render a table inside a closed-by-default details panel."""
    table_html = data.to_html(index=index, escape=True, border=0)
    table_html = (
        '<div class="report-table-scroll" role="region" '
        f'aria-label="{escape(summary, quote=True)}" tabindex="0">'
        + table_html + "</div>"
    )
    if max_height is not None:
        table_html = (
            f'<div style="max-height: {int(max_height)}px; overflow: auto;" '
            'tabindex="0">' + table_html + "</div>"
        )
    display(HTML(
        f'<details class="report-details"><summary>{escape(summary)}</summary>'
        f"{table_html}</details>"
    ))


def display_collapsible_figure(figure, summary, *, open_by_default=False):
    """Render a figure in a details panel, optionally open on load."""
    image_buffer = BytesIO()
    figure.savefig(image_buffer, format="png", bbox_inches="tight")
    image_data = base64.b64encode(image_buffer.getvalue()).decode("ascii")
    plt.close(figure)
    escaped_summary = escape(summary, quote=True)
    display(HTML(
        '<details class="report-details report-details-figure"'
        + (" open" if open_by_default else "")
        + ">"
        f"<summary>{escape(summary)}</summary>"
        f'<img alt="{escaped_summary}" src="data:image/png;base64,{image_data}">'
        "</details>"
    ))


WGCNA_MAX_DISPLAY_MODULES = 25


def wgcna_modules_for_display(summary, *, limit=WGCNA_MAX_DISPLAY_MODULES):
    """Return the largest modules retained in compact report figures."""
    return summary.sort_values(
        ["n_genes", "module"], ascending=[False, True], kind="stable"
    )["module"].head(limit).tolist()


def wgcna_display_caption(n_shown, n_total):
    """Describe when a WGCNA figure omits smaller modules for legibility."""
    if n_shown == n_total:
        return ""
    return (
        f" Figures show the {n_shown} largest modules of {n_total}; "
        "complete module-level results remain available in the downloadable tables."
    )


def load_residual_sample_clusters(de_dir, fit, settings):
    """Cluster one cell type's fitted age-bin residuals for its report."""
    if fit.get("status") != "complete":
        return None, None, 0, None
    residual_path = Path(de_dir) / fit.get("pearson_residuals", "")
    metadata_path = Path(de_dir) / fit.get("pearson_residual_sample_metadata", "")
    if not residual_path.is_file() or not metadata_path.is_file():
        return None, None, 0, None
    with np.load(residual_path, allow_pickle=False) as stored:
        residuals = pd.DataFrame(
            stored["residuals"],
            index=stored["pseudobulk_id"].astype(str),
            columns=stored["genes"].astype(str),
        )
    samples = pd.read_csv(metadata_path, dtype={"pseudobulk_id": str})
    samples["sample_key"] = samples["study"].astype(str) + "|" + samples["sample"].astype(str)
    sample_by_id = samples.drop_duplicates("pseudobulk_id").set_index("pseudobulk_id")
    shared_ids = residuals.index.intersection(sample_by_id.index, sort=False)
    residuals = residuals.loc[shared_ids]
    sample_by_id = sample_by_id.loc[shared_ids]
    residuals.index = sample_by_id["sample_key"].to_numpy()
    keep = ~residuals.index.duplicated(keep="first")
    residuals = residuals.loc[keep]
    sample_by_key = sample_by_id.loc[~sample_by_id["sample_key"].duplicated()].set_index("sample_key")
    annotations, pc_scores, n_components = cluster_sample_residuals(
        residuals, sample_by_key, settings
    )
    residuals = residuals.loc[annotations["sample_key"].astype(str)].copy()
    residuals.index = annotations["sample_key"].astype(str)
    return annotations, pc_scores, n_components, residuals


def load_merged_cmv_results(de_dir):
    """Load the per-cell-type merged CMV model when it completed."""
    if not de_dir:
        return None
    path = Path(de_dir) / "combined" / "cmv.csv"
    if not path.is_file():
        return None
    return pd.read_csv(path)


def load_merged_smoking_results(de_dir):
    """Load the per-cell-type merged smoking-status model when it completed."""
    if not de_dir:
        return None
    path = Path(de_dir) / "combined" / "smoking_status.csv"
    if not path.is_file():
        return None
    return pd.read_csv(path)


def publish_residual_wgcna(result):
    """Publish residual-network artifacts beside the per-type report."""
    wgcna_dir = output_dir / "residual_wgcna"
    wgcna_dir.mkdir(parents=True, exist_ok=True)
    (wgcna_dir / "metadata.json").write_text(json.dumps(result.metadata, indent=2) + "\n")
    result.soft_threshold.to_csv(wgcna_dir / "soft_threshold.csv", index=False)
    result.module_assignments.to_csv(wgcna_dir / "module_assignments.csv", index=False)
    result.eigengenes.to_csv(wgcna_dir / "eigengenes.csv", index=False)
    result.module_summary.to_csv(wgcna_dir / "module_summary.csv", index=False)
    return wgcna_dir


def wgcna_sample_connectivity(residuals, annotations):
    """Summarize each sample's signed-residual-network connectivity."""
    numeric = residuals.apply(pd.to_numeric, errors="coerce").dropna(axis=0, how="any")
    if len(numeric) < 3:
        return pd.DataFrame()
    correlations = numeric.T.corr().clip(-1, 1)
    adjacency = (correlations + 1) / 2
    np.fill_diagonal(adjacency.values, 0)
    connectivity = adjacency.sum(axis=1)
    standard_deviation = connectivity.std(ddof=0)
    standardized = (
        (connectivity - connectivity.mean()) / standard_deviation
        if standard_deviation > 0 else pd.Series(0.0, index=connectivity.index)
    )
    result = pd.DataFrame({
        "sample_key": connectivity.index.astype(str),
        "sample_connectivity": connectivity.to_numpy(),
        "standardized_connectivity": standardized.to_numpy(),
    })
    metadata = annotations.drop_duplicates("sample_key")
    columns = [column for column in ("sample_key", "study", "age", "sex") if column in metadata]
    return result.merge(metadata[columns], on="sample_key", how="left", validate="one_to_one")


def wgcna_module_membership(result, residuals):
    """Calculate correlations between each gene and its assigned eigengene."""
    assignments = result.module_assignments.copy()
    eigengenes = result.eigengenes.set_index("sample_key")
    shared_samples = residuals.index.intersection(eigengenes.index, sort=False)
    assignments["module_membership"] = np.nan
    for module, group in assignments.groupby("module", observed=True):
        eigengene_column = f"ME{module}"
        genes = [gene for gene in group["gene"] if gene in residuals.columns]
        if eigengene_column not in eigengenes or not genes:
            continue
        eigengene = pd.to_numeric(eigengenes.loc[shared_samples, eigengene_column], errors="coerce")
        values = residuals.loc[shared_samples, genes].apply(pd.to_numeric, errors="coerce")
        assignments.loc[group.index, "module_membership"] = values.corrwith(eigengene).to_numpy()
    return assignments


def wgcna_module_rankings(result, residuals):
    """Rank all network genes by correlation with each non-grey eigengene."""
    eigengenes = result.eigengenes.set_index("sample_key")
    shared_samples = residuals.index.intersection(eigengenes.index, sort=False)
    rankings = {}
    for module in result.module_summary["module"].astype(str):
        if module.lower() == "grey":
            continue
        eigengene_column = f"ME{module}"
        if eigengene_column not in eigengenes:
            continue
        eigengene = pd.to_numeric(eigengenes.loc[shared_samples, eigengene_column], errors="coerce")
        values = residuals.loc[shared_samples].apply(pd.to_numeric, errors="coerce")
        membership = values.corrwith(eigengene).replace([np.inf, -np.inf], np.nan).dropna()
        membership = membership.sort_values(ascending=False, kind="stable")
        if not membership.empty:
            rankings[module] = membership
    return rankings


def run_wgcna_module_gsea(result, residuals):
    """Run preranked GO BP GSEA for every retained WGCNA module."""
    run = run_preranked_gsea(wgcna_module_rankings(result, residuals))
    return GeneSetEnrichmentRun(
        results=run.results.rename(columns={
            "analysis_set": "module", "analysis_set_fdr": "module_fdr",
        }),
        gene_sets=run.gene_sets,
        provenance=run.provenance,
        failures=run.failures,
    )


def display_merged_age_go_enrichment(age_results, contrast):
    """Run and display signed GO BP GSEA for one merged age-DE contrast."""
    ranking = signed_de_gene_ranking(age_results)
    ranking_label = (
        "signed Wald statistic" if ranking.name == "wald_statistic"
        else "signed −log10(p-value)"
    )
    display(Markdown("##### GO Biological Process enrichment of merged age effects"))
    if ranking.dropna().shape[0] < 10:
        display(Markdown("_Too few finite age statistics for preranked enrichment._"))
        return

    gsea_run = run_preranked_gsea({
        f"Merged age ({contrast})": ranking,
    })
    if gsea_run.results.empty and gsea_run.failures:
        display(Markdown(
            "_GO enrichment could not run: "
            f"{escape('; '.join(gsea_run.failures.values()))}._"
        ))
        return
    if gsea_run.results.empty:
        display(Markdown("_No GO Biological Process terms met the GSEA size limits._"))
        return

    enrichment_dir = output_dir / "go_enrichment"
    enrichment_dir.mkdir(parents=True, exist_ok=True)
    gsea_run.results.to_csv(enrichment_dir / "merged_age_gsea.csv", index=False)
    if gsea_run.failures:
        display(Markdown(
            "_Some enrichment rankings were skipped: "
            f"{escape('; '.join(gsea_run.failures.values()))}. Details are recorded in provenance._"
        ))
    try:
        summary, _ = save_go_specificity_outputs(
            gsea_run, output_dir=enrichment_dir, prefix="merged_age_gsea",
        )
    except Exception as exc:  # noqa: BLE001 - keep DE report usable if ontology is unavailable
        display(Markdown(
            "_The complete GO results were saved, but the hierarchy-aware summary could not be "
            f"built: {escape(f'{type(exc).__name__}: {exc}')}. "
            "See [full GSEA results](go_enrichment/merged_age_gsea.csv)._"
        ))
        return

    if summary.empty:
        display(Markdown(
            "_No GO terms passed the native GSEA FDR threshold of 0.05. See [full ranked results]"
            "(go_enrichment/merged_age_gsea.csv)._"
        ))
        return
    display(Markdown(
        f"All genes tested in the merged model are ranked by the {ranking_label}; positive NES "
        "indicates enrichment toward genes increasing with age and negative NES toward genes "
        "decreasing with age. The table includes GSEA FDR-significant terms after collapsing "
        "redundant parent terms when a same-direction descendant has at least 0.5 leading-edge "
        "Jaccard overlap. Full results and the term audit are available as [full GSEA results]"
        "(go_enrichment/merged_age_gsea.csv), [specific summary]"
        "(go_enrichment/merged_age_gsea_specific_summary.csv), and [term audit]"
        "(go_enrichment/merged_age_gsea_specificity_audit.csv). The exact library, ontology, "
        "parameters, and any failed rankings are recorded in [provenance]"
        "(go_enrichment/merged_age_gsea_provenance.json)."
    ))
    age_summary = summary.copy()
    leading_edge = age_summary["leading_edge_genes"].map(leading_edge_gene_preview)
    age_summary["leading_edge_gene_count"] = leading_edge.map(lambda value: value[0])
    age_summary["leading_edge_gene_preview"] = leading_edge.map(lambda value: value[1])
    age_summary = age_summary.loc[:, [
        "gene_set", "normalized_enrichment_score", "analysis_set_fdr",
        "leading_edge_gene_count", "leading_edge_gene_preview",
    ]].rename(columns={
        "gene_set": "GO Biological Process",
        "normalized_enrichment_score": "NES",
        "analysis_set_fdr": "GSEA FDR q-value",
        "leading_edge_gene_count": "Leading-edge genes (n)",
        "leading_edge_gene_preview": "Leading-edge genes",
    })
    display_collapsible_table(
        age_summary, "Show specific GO terms and their leading-edge genes", max_height=560,
    )


def wgcna_covariate_associations(result, annotations):
    """Associate module eigengenes with recorded continuous and categorical covariates."""
    eigengenes = result.eigengenes.copy()
    metadata = annotations.drop_duplicates("sample_key").copy()
    data = eigengenes.merge(metadata, on="sample_key", how="inner", validate="one_to_one")
    continuous = [
        column for column in ("age", "log10_total_counts", "bmi") if column in data
    ]
    categorical = [
        column for column in (
            "sex", "cmv", "smoking_status", "study", "study_site", "technology",
            "batch_single_cell",
            "include_intronic",
        ) if column in data
    ]
    rows = []
    eigengene_columns = [column for column in eigengenes if column.startswith("ME")]
    for eigengene_column in eigengene_columns:
        module = eigengene_column.removeprefix("ME")
        for covariate in continuous:
            values = data[[eigengene_column, covariate]].apply(pd.to_numeric, errors="coerce").dropna()
            if len(values) < 3 or values[covariate].nunique() < 2:
                continue
            effect_size, pvalue = spearmanr(values[eigengene_column], values[covariate])
            rows.append({
                "module": module, "covariate": covariate, "association": "Spearman correlation",
                "effect_size": effect_size, "pvalue": pvalue, "n_samples": len(values),
                "n_levels": np.nan,
            })
        for covariate in categorical:
            values = data[[eigengene_column, covariate]].dropna().copy()
            counts = values[covariate].astype(str).value_counts()
            counts = counts.loc[counts >= 2]
            values = values.loc[values[covariate].astype(str).isin(counts.index)]
            if len(values) < 3 or len(counts) < 2:
                continue
            codes, _ = pd.factorize(values[covariate].astype(str), sort=True)
            if len(counts) == 2:
                effect_size, pvalue = pointbiserialr(codes, values[eigengene_column])
                association = "Point-biserial correlation"
            else:
                groups = [
                    values.loc[codes == code, eigengene_column].to_numpy()
                    for code in range(len(counts))
                ]
                _, pvalue = f_oneway(*groups)
                overall_mean = values[eigengene_column].mean()
                total_sum_squares = ((values[eigengene_column] - overall_mean) ** 2).sum()
                between_sum_squares = sum(
                    len(group) * (np.mean(group) - overall_mean) ** 2 for group in groups
                )
                effect_size = between_sum_squares / total_sum_squares if total_sum_squares else np.nan
                association = "One-way ANOVA η²"
            rows.append({
                "module": module, "covariate": covariate, "association": association,
                "effect_size": effect_size, "pvalue": pvalue, "n_samples": len(values),
                "n_levels": len(counts),
            })
    associations = pd.DataFrame(rows)
    if associations.empty:
        return associations
    associations["fdr"] = np.nan
    valid = associations["pvalue"].notna()
    associations.loc[valid, "fdr"] = multipletests(
        associations.loc[valid, "pvalue"], method="fdr_bh",
    )[1]
    return associations.sort_values(["fdr", "module", "covariate"], na_position="last")


def wgcna_leave_one_study_out_stability(result, residuals, annotations, settings):
    """Compare full-data modules with modules refitted after omitting each study."""
    if "study" not in annotations:
        return pd.DataFrame()
    studies = sorted(annotations["study"].dropna().astype(str).unique())
    if len(studies) < 2:
        return pd.DataFrame()
    full_modules = {
        str(module): set(group["gene"].astype(str))
        for module, group in result.module_assignments.groupby("module", observed=True)
    }
    rows = []
    for held_out_study in studies:
        retained = annotations.loc[annotations["study"].astype(str) != held_out_study].copy()
        retained_residuals = residuals.loc[residuals.index.intersection(
            retained["sample_key"].astype(str), sort=False
        )]
        refit = run_residual_wgcna(retained_residuals, retained, settings=settings)
        if refit.metadata["status"] != "complete":
            rows.append({
                "held_out_study": held_out_study,
                "module": "—",
                "best_refit_module": "—",
                "jaccard_overlap": np.nan,
                "n_retained_samples": len(retained_residuals),
                "status": refit.metadata.get("reason", "refit failed"),
            })
            continue
        refit_modules = [
            (str(module), set(group["gene"].astype(str)))
            for module, group in refit.module_assignments.groupby("module", observed=True)
        ]
        for module, genes in full_modules.items():
            best_module, best_overlap = max(
                (
                    (refit_module, len(genes & refit_genes) / len(genes | refit_genes))
                    for refit_module, refit_genes in refit_modules
                    if genes | refit_genes
                ),
                key=lambda item: item[1], default=("—", np.nan),
            )
            rows.append({
                "held_out_study": held_out_study,
                "module": module,
                "best_refit_module": best_module,
                "jaccard_overlap": best_overlap,
                "n_retained_samples": len(retained_residuals),
                "status": "complete",
            })
    return pd.DataFrame(rows)


def plot_residual_wgcna(result, residuals, annotations, settings, cmv_results, smoking_results):
    """Show WGCNA diagnostics, covariates, hubs, and cross-study stability."""
    if result.metadata["status"] != "complete":
        status = result.metadata["status"]
        message = result.metadata.get("reason", "unknown reason")
        if status == "failed":
            display(Markdown(f"_Residual WGCNA could not be completed: {message}._"))
        else:
            display(Markdown(f"_Residual WGCNA was not run: {message}._"))
        return
    wgcna_dir = output_dir / "residual_wgcna"

    display(Markdown("#### Network diagnostics"))
    soft_threshold = result.soft_threshold
    if not soft_threshold.empty and {"Power", "SFT.R.sq", "mean(k)"}.issubset(soft_threshold):
        figure, axes = plt.subplots(1, 2, figsize=(12, 4.3))
        signed_fit = -np.sign(soft_threshold["slope"]) * soft_threshold["SFT.R.sq"]
        axes[0].plot(soft_threshold["Power"], signed_fit, marker="o", color="#4c72b0")
        axes[0].axhline(
            result.metadata["settings"]["scale_free_fit_threshold"], color="#c44e52",
            linestyle="--", linewidth=1,
        )
        axes[0].axvline(result.metadata["selected_soft_threshold_power"], color="#555555", linewidth=1)
        axes[0].set(xlabel="Soft-threshold power", ylabel="Signed scale-free topology fit", title="Soft-threshold diagnostic")
        axes[1].plot(soft_threshold["Power"], soft_threshold["mean(k)"], marker="o", color="#55a868")
        axes[1].axvline(result.metadata["selected_soft_threshold_power"], color="#555555", linewidth=1)
        axes[1].set(xlabel="Soft-threshold power", ylabel="Mean connectivity", title="Network connectivity")
        figure.tight_layout()
        display_collapsible_figure(figure, "Residual WGCNA soft-threshold diagnostics")

    display(Markdown("#### Sample connectivity"))
    connectivity = wgcna_sample_connectivity(residuals, annotations)
    connectivity.to_csv(wgcna_dir / "sample_connectivity.csv", index=False)
    if connectivity.empty:
        display(Markdown("_Sample connectivity is unavailable for fewer than three residual samples._"))
    else:
        figure, axis = plt.subplots(figsize=(11, 4.4))
        plot_data = connectivity.sort_values("standardized_connectivity").reset_index(drop=True)
        sns.scatterplot(
            data=plot_data, x=plot_data.index, y="standardized_connectivity",
            hue="study" if "study" in plot_data else None, palette="tab20", s=46, ax=axis,
        )
        axis.axhline(-2, color="#c44e52", linestyle="--", linewidth=1, label="−2 SD")
        axis.axhline(0, color="#666666", linewidth=0.8)
        axis.set(xlabel="Samples ranked by connectivity", ylabel="Standardized sample connectivity")
        axis.legend(title="Study", bbox_to_anchor=(1.02, 1), loc="upper left")
        figure.subplots_adjust(right=0.82)
        display_collapsible_figure(figure, "Show sample-connectivity diagnostic")
        display_collapsible_table(
            connectivity.sort_values("standardized_connectivity").head(15),
            "Show lowest-connectivity residual samples",
        )

    summary = result.module_summary.copy()
    if summary.empty:
        display(Markdown("_PyWGCNA did not retain any non-grey modules._"))
        return
    display(Markdown("#### Module structure"))
    display(Markdown(
        "Module color names are WGCNA identifiers, not biological labels, and can change when "
        "the network is refit; the stability check matches modules by shared genes. Grey marks "
        "genes that were not assigned to a module. The `n_cmv_yes/no` and `n_smoking_yes/no` fields "
        "count samples with observed yes/no status used in each eigengene association. They "
        "repeat across modules because eigengenes use the same samples."
    ))
    summary_columns = [
        "module", "n_genes", "n_cmv_associated_genes", "cmv_eigengene_correlation",
        "cmv_eigengene_fdr", "n_smoking_associated_genes",
        "smoking_eigengene_correlation", "smoking_eigengene_fdr",
    ]
    summary_display = summary[[column for column in summary_columns if column in summary]].rename(
        columns={
            "module": "Module", "n_genes": "Genes",
            "n_cmv_associated_genes": "CMV-DE genes",
            "cmv_eigengene_correlation": "CMV–ME r",
            "cmv_eigengene_fdr": "CMV–ME FDR",
            "n_smoking_associated_genes": "Smoking-DE genes",
            "smoking_eigengene_correlation": "Smoking–ME r",
            "smoking_eigengene_fdr": "Smoking–ME FDR",
        }
    )
    display_collapsible_table(
        summary_display, "Show residual WGCNA module summary", max_height=480,
    )
    display_modules = wgcna_modules_for_display(summary)
    shown = summary.loc[summary["module"].isin(display_modules)].sort_values(
        "n_genes", ascending=True
    )
    display_caption = wgcna_display_caption(len(shown), len(summary))
    if display_caption:
        display(Markdown(f"_{display_caption.strip()}_"))
    figure_label_suffix = (
        f" ({len(shown)} largest of {len(summary)})" if display_caption else ""
    )
    figure, axes = plt.subplots(1, 3, figsize=(20, max(4.5, 0.36 * len(shown) + 2)))
    colors = [module if module in mcolors.get_named_colors_mapping() else "#4c72b0" for module in shown["module"]]
    axes[0].barh(
        shown["module"], shown["n_genes"], color=colors,
        edgecolor="#38434a", linewidth=0.55,
    )
    axes[0].set(xlabel="Genes", title="Residual WGCNA module sizes" + figure_label_suffix)
    for axis, covariate, label, de_results in (
        (axes[1], "cmv", "CMV", cmv_results),
        (axes[2], "smoking", "Smoking", smoking_results),
    ):
        correlation_column = f"{covariate}_eigengene_correlation"
        overlap_column = f"{covariate}_gene_overlap_fdr"
        fdr_column = f"{covariate}_eigengene_fdr"
        trait_data = shown.dropna(subset=[correlation_column]).copy()
        if trait_data.empty:
            axis.text(0.5, 0.5, f"{label} status unavailable\nor no within-group variation",
                      ha="center", va="center")
            axis.set_axis_off()
            continue
        overlap_fdr = pd.to_numeric(trait_data[overlap_column], errors="coerce")
        if de_results is None or not overlap_fdr.notna().any():
            axis.text(
                0.5, 0.06, f"{label}-DE gene-overlap evidence unavailable",
                ha="center", va="bottom", transform=axis.transAxes, fontsize=8,
            )
            trait_data["overlap_score"] = 0.0
        else:
            trait_data["overlap_score"] = -np.log10(overlap_fdr.clip(lower=1e-300))
        points = axis.scatter(
            trait_data[correlation_column], trait_data["overlap_score"],
            s=np.maximum(30, trait_data["n_genes"] * 3), c=trait_data[fdr_column],
            cmap="viridis_r", vmin=0, vmax=1, alpha=0.85, linewidths=0,
        )
        for row in trait_data.itertuples(index=False):
            axis.annotate(row.module, (getattr(row, correlation_column), row.overlap_score),
                          xytext=(3, 3), textcoords="offset points", fontsize=8)
        axis.axvline(0, color="#666666", linewidth=0.8)
        axis.set(
            xlabel=f"{label}-status vs eigengene correlation",
            ylabel=f"−log10({label}-DE overlap FDR)",
            title=f"{label} checks of residual modules" + figure_label_suffix,
        )
        figure.colorbar(points, ax=axis, label=f"Eigengene–{label} FDR")
    figure.tight_layout()
    display_collapsible_figure(figure, "Residual WGCNA modules and CMV diagnostics")

    for covariate, label, status_field in (
        ("cmv", "CMV", "cmv_status"),
        ("smoking_status", "smoking", "smoking_status"),
    ):
        score_data = trait_associated_module_scores(result, annotations, covariate=covariate)
        if score_data.empty or not {"umap_1", "umap_2"}.issubset(annotations):
            continue
        score_data = score_data.rename(columns={"trait_status": status_field})
        score_data = score_data.merge(
            annotations[["sample_key", "umap_1", "umap_2"]].drop_duplicates("sample_key"),
            on="sample_key", how="left", validate="one_to_one",
        ).dropna(subset=["module_score", "umap_1", "umap_2"])
        score_data.to_csv(wgcna_dir / f"{label}_associated_module_scores.csv", index=False)
        if not score_data.empty:
            module = score_data["module"].iloc[0]
            figure, axes = plt.subplots(1, 2, figsize=(13, 5.2))
            points = axes[0].scatter(
                score_data["umap_1"], score_data["umap_2"],
                c=score_data["module_score"], cmap="viridis", s=28, alpha=0.9,
                linewidths=0,
            )
            axes[0].set(
                xlabel="UMAP 1", ylabel="UMAP 2",
                title=f"Residual UMAP: {label}-associated module {module}",
            )
            figure.colorbar(points, ax=axes[0], label="Oriented module eigengene score")
            order = [status for status in ("no", "yes", "not provided")
                     if status in score_data[status_field].unique()]
            group_scores = score_data.groupby(status_field, observed=True)["module_score"]
            means = group_scores.mean().reindex(order)
            confidence = (1.96 * group_scores.sem()).reindex(order).fillna(0)
            axes[1].bar(
                order, means, yerr=confidence, color="#9ecae1", edgecolor="#38434a",
                linewidth=0.8, capsize=4,
            )
            sns.stripplot(
                data=score_data, x=status_field, y="module_score", order=order,
                color="#315a73", alpha=0.45, size=3, jitter=0.18, ax=axes[1],
            )
            axes[1].set(
                xlabel=f"{label.title()} status", ylabel="Oriented module eigengene score",
                title=f"Module {module} mean score by {label} status\nBars: mean ± 1.96 SE; points: samples",
            )
            figure.suptitle(
                f"Highest absolute eigengene–{label} correlation; score sign oriented so higher "
                f"values correspond to {label} yes", y=1.02, fontsize=10,
            )
            figure.tight_layout()
            display_collapsible_figure(
                figure, f"Show the strongest {label}-associated module on residual UMAP",
                open_by_default=True,
            )

    eigengenes = result.eigengenes.set_index("sample_key")
    module_columns = [
        f"ME{module}" for module in display_modules if f"ME{module}" in eigengenes
    ]
    if len(module_columns) > 1:
        figure, axis = plt.subplots(figsize=(
            max(5, 0.5 * len(module_columns) + 3),
            max(4, 0.5 * len(module_columns) + 2),
        ))
        sns.heatmap(eigengenes[module_columns].corr(), cmap="vlag", center=0, vmin=-1, vmax=1,
                    annot=True, fmt=".2f", square=True, ax=axis)
        axis.set_title("Module eigengene correlation" + figure_label_suffix)
        figure.tight_layout()
        display_collapsible_figure(figure, "Show module-eigengene correlation heatmap")

    display(Markdown("#### Covariate associations"))
    associations = wgcna_covariate_associations(result, annotations)
    associations.to_csv(wgcna_dir / "module_covariate_associations.csv", index=False)
    if associations.empty:
        display(Markdown("_No recorded covariates had sufficient variation for eigengene association tests._"))
    else:
        display(Markdown(
            "Continuous covariates use Spearman correlation. Binary covariates use point-biserial "
            "correlation (signed r; for CMV and smoking, positive means higher eigengene scores "
            "in yes than no. Other binary covariates use sorted category coding). Multi-level covariates "
            "use one-way ANOVA η², the fraction of eigengene variance explained by category "
            "membership (0 to 1, unsigned). These are unadjusted, descriptive associations. FDR "
            "is controlled across all module–covariate tests."
        ))
        plotted_associations = associations.loc[associations["module"].isin(display_modules)]
        figure, axes = plt.subplots(
            1, 2, figsize=(15, max(4, 0.42 * plotted_associations["module"].nunique() + 2)),
        )
        continuous = plotted_associations.loc[
            plotted_associations["association"] == "Spearman correlation"
        ]
        categorical = plotted_associations.loc[
            plotted_associations["association"] != "Spearman correlation"
        ]
        for axis, subset, title, cmap, center in (
            (axes[0], continuous, "Continuous covariates" + figure_label_suffix, "vlag", 0),
            (axes[1], categorical, "Categorical covariates" + figure_label_suffix, "YlGnBu", None),
        ):
            if subset.empty:
                axis.text(0.5, 0.5, "No eligible covariates", ha="center", va="center")
                axis.set_axis_off()
                continue
            matrix = subset.pivot(index="module", columns="covariate", values="effect_size")
            sns.heatmap(matrix, cmap=cmap, center=center, annot=True, fmt=".2f", linewidths=0.5,
                        cbar_kws={"label": "Effect size"}, ax=axis)
            axis.set_title(title)
            axis.set_xlabel("")
            axis.set_ylabel("Module")
        figure.tight_layout()
        display_collapsible_figure(
            figure, "Show module-eigengene covariate association heatmaps", open_by_default=True,
        )
        display_collapsible_table(
            associations, "Show module-eigengene covariate association statistics", max_height=480,
        )

    display(Markdown("#### Hub genes and CMV evidence"))
    membership = wgcna_module_membership(result, residuals)
    hub_genes = membership.assign(
        absolute_module_membership=membership["module_membership"].abs()
    ).sort_values(["module", "absolute_module_membership"], ascending=[True, False])
    hub_genes["module_rank"] = hub_genes.groupby("module", observed=True).cumcount() + 1
    hub_genes.to_csv(wgcna_dir / "hub_genes.csv", index=False)
    display_collapsible_table(
        hub_genes.loc[hub_genes["module_rank"] <= 20], "Show top 20 hub genes per module", max_height=560,
    )
    display(Markdown("#### Gene-set enrichment (preranked GSEA)"))
    gsea_run = run_wgcna_module_gsea(result, residuals)
    gsea_results = gsea_run.results
    if gsea_results.empty and gsea_run.failures:
        display(Markdown(
            "_GO Biological Process GSEA could not run: "
            f"{escape('; '.join(gsea_run.failures.values()))}._"
        ))
    elif gsea_results.empty:
        display(Markdown("_No GO Biological Process terms met the configured GSEA size limits._"))
    else:
        gsea_results.to_csv(wgcna_dir / "module_gsea.csv", index=False)
        if gsea_run.failures:
            display(Markdown(
                "_Some modules could not be enriched: "
                f"{escape('; '.join(gsea_run.failures.values()))}. Details are recorded in provenance._"
            ))
        display(Markdown(
            "Genes are ranked by their Pearson correlation with each module eigengene, using all "
            "selected network genes. Grey/unassigned genes are excluded as modules. Results use "
            "GO Biological Process 2023; the full table includes the native GSEA FDR q-value "
            "within each module."
        ))
        summary_error = None
        try:
            specific_gsea, _ = save_go_specificity_outputs(
                gsea_run, output_dir=wgcna_dir, prefix="module_gsea", group_column="module",
                fdr_column="module_fdr",
            )
        except Exception as exc:  # noqa: BLE001 - preserve the full GSEA output if GO is unavailable
            summary_error = exc
            specific_gsea = pd.DataFrame()
            display(Markdown(
                "_The GO-specific summary could not be built from the bundled ontology: "
                f"{escape(f'{type(exc).__name__}: {exc}')}. "
                "The complete preranked results remain available in "
                "[module_gsea.csv](residual_wgcna/module_gsea.csv)._"
            ))
        if specific_gsea.empty:
            if summary_error is None:
                display(Markdown(
                    "_No GO terms passed the native GSEA FDR threshold of 0.05. The full ranked "
                    "results remain available in [module_gsea.csv]"
                    "(residual_wgcna/module_gsea.csv)._"
                ))
        else:
            display(Markdown(
                "This summary shows terms with native GSEA FDR ≤ 0.05. A broader GO parent is "
                "collapsed only when a same-direction descendant has at least 0.5 leading-edge "
                "Jaccard overlap; terms without such a supported child remain. The complete "
                "results and summary audit are available as [full GSEA results]"
                "(residual_wgcna/module_gsea.csv), [specific summary]"
                "(residual_wgcna/module_gsea_specific_summary.csv), and [term audit]"
                "(residual_wgcna/module_gsea_specificity_audit.csv). The exact library, ontology, "
                "parameters, and any failed modules are recorded in [provenance]"
                "(residual_wgcna/module_gsea_provenance.json)."
            ))
            display_collapsible_table(
                specific_gsea, "Show specific, nonredundant GO terms per module", max_height=560,
            )
        plotted_gsea = specific_gsea.head(15).copy()
        if not plotted_gsea.empty:
            plotted_gsea["minus_log10_module_fdr"] = -np.log10(
                plotted_gsea["module_fdr"].clip(lower=1e-300)
            )
            plotted_gsea["label"] = plotted_gsea["module"] + ": " + plotted_gsea["gene_set"]
            figure, axis = plt.subplots(figsize=(10, max(4.5, 0.34 * len(plotted_gsea) + 1.5)))
            axis.barh(
                plotted_gsea["label"].iloc[::-1],
                plotted_gsea["minus_log10_module_fdr"].iloc[::-1],
                color="#7aa6c2", edgecolor="#38434a", linewidth=0.55,
            )
            axis.set(xlabel="−log10(module GSEA FDR)", title="Specific module GO enrichments")
            figure.tight_layout()
            display_collapsible_figure(figure, "Show top module GO enrichments")
    if cmv_results is not None and {"gene", "log2FoldChange", "padj"}.issubset(cmv_results):
        cmv = cmv_results[["gene", "log2FoldChange", "padj"]].copy()
        cmv["log2FoldChange"] = pd.to_numeric(cmv["log2FoldChange"], errors="coerce")
        cmv["padj"] = pd.to_numeric(cmv["padj"], errors="coerce")
        evidence = membership.merge(cmv, on="gene", how="left", validate="one_to_one").dropna(
            subset=["module_membership", "log2FoldChange", "padj"]
        )
        if not evidence.empty:
            figure, axis = plt.subplots(figsize=(8.5, 5.5))
            points = axis.scatter(
                evidence["module_membership"], evidence["log2FoldChange"],
                c=-np.log10(evidence["padj"].clip(lower=np.finfo(float).tiny)), cmap="viridis",
                s=28, alpha=0.8, linewidths=0,
            )
            axis.axhline(0, color="#666666", linewidth=0.8)
            axis.axvline(0, color="#666666", linewidth=0.8)
            axis.set(xlabel="Gene membership in assigned module", ylabel="CMV log2 fold change",
                     title="Module membership and CMV-DE effect")
            figure.colorbar(points, ax=axis, label="−log10(CMV adjusted p-value)")
            figure.tight_layout()
            display_collapsible_figure(figure, "Show module-membership versus CMV-DE evidence")

    display(Markdown("#### Leave-one-study-out stability"))
    stability = wgcna_leave_one_study_out_stability(result, residuals, annotations, settings)
    stability.to_csv(wgcna_dir / "leave_one_study_out_stability.csv", index=False)
    completed = stability.loc[stability.get("status", pd.Series(dtype=str)).eq("complete")].copy()
    if completed.empty:
        display(Markdown("_Leave-one-study-out module refits require at least two studies with sufficient retained samples._"))
    else:
        plotted_stability = completed.loc[completed["module"].isin(display_modules)]
        figure, axis = plt.subplots(
            figsize=(
                max(6, 0.9 * plotted_stability["held_out_study"].nunique() + 4),
                max(4, 0.36 * plotted_stability["module"].nunique() + 2),
            )
        )
        matrix = plotted_stability.pivot(
            index="module", columns="held_out_study", values="jaccard_overlap"
        )
        best_matches = plotted_stability.pivot(
            index="module", columns="held_out_study", values="best_refit_module"
        ).reindex(index=matrix.index, columns=matrix.columns)
        labels = matrix.copy().astype(object)
        for module_name in matrix.index:
            for study_name in matrix.columns:
                overlap = matrix.loc[module_name, study_name]
                match = best_matches.loc[module_name, study_name]
                labels.loc[module_name, study_name] = (
                    f"{match}\n{overlap:.2f}" if pd.notna(overlap) and pd.notna(match) else ""
                )
        sns.heatmap(matrix, vmin=0, vmax=1, cmap="mako", annot=labels, fmt="",
                    annot_kws={"fontsize": 7}, ax=axis,
                    cbar_kws={"label": "Best refit-module Jaccard overlap"})
        axis.set(
            xlabel="Held-out study", ylabel="Full-data module",
            title="Leave-one-study-out module stability" + figure_label_suffix,
        )
        axis.tick_params(axis="x", labelrotation=35)
        figure.tight_layout()
        display_collapsible_figure(figure, "Show leave-one-study-out module stability", open_by_default=True)
        display_collapsible_table(
            stability, "Show leave-one-study-out module matches", max_height=480,
        )


def annotate_cluster_centroids(axis, table, label_column):
    labels = pd.to_numeric(table[label_column], errors="coerce")
    points = table.assign(_cluster=labels).dropna(
        subset=["_cluster", "umap_1", "umap_2"]
    )
    for cluster, group in points.loc[points["_cluster"] >= 0].groupby("_cluster"):
        axis.annotate(
            str(int(cluster)),
            (group["umap_1"].mean(), group["umap_2"].mean()),
            ha="center", va="center", fontsize=9, fontweight="bold",
            color="#111111",
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.8,
                  "edgecolor": "none"},
            zorder=6,
        )


def plot_residual_cluster_markers(residuals, annotations, settings):
    """Show the residual UMAP clusters and their largest mean-residual contrasts."""
    display(Markdown("### Residual Cluster Markers"))
    coordinates = annotations.dropna(subset=["umap_1", "umap_2"]).copy()
    if coordinates.empty:
        display(Markdown("_Residual UMAP coordinates are unavailable for this cell type._"))
    else:
        cluster_counts = annotations["residual_cluster"].value_counts()
        coordinates["residual_cluster_label"] = coordinates["residual_cluster"].map(
            lambda cluster: (
                f"Noise (n={cluster_counts[cluster]})" if int(cluster) == -1
                else f"Cluster {int(cluster)} (n={cluster_counts[cluster]})"
            )
        )
        figure, axis = plt.subplots(figsize=(6.6, 4.8))
        sns.scatterplot(
            data=coordinates, x="umap_1", y="umap_2", hue="residual_cluster_label",
            palette="tab20", s=40, alpha=0.84, linewidth=0, ax=axis,
        )
        annotate_cluster_centroids(axis, coordinates, "residual_cluster")
        axis.set_title("Residual UMAP: sample clusters")
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        if axis.legend_ is not None:
            handles, labels = axis.get_legend_handles_labels()
            axis.legend(
                handles, labels, title="Residual cluster", loc="best",
                ncol=2 if len(labels) > 6 else 1,
                fontsize=7, title_fontsize=8, framealpha=0.9,
                borderpad=0.35, labelspacing=0.25, handletextpad=0.35,
            )
        figure.tight_layout()
        plt.show()

    markers, _, reference_cluster = summarize_residual_cluster_markers(
        residuals, annotations, settings.residual_marker_max_clusters
    )
    if markers.empty:
        display(Markdown(
            "_No between-cluster marker contrasts are available; at least two residual "
            "clusters are needed._"
        ))
        return

    display(Markdown(
        f"For up to the {settings.residual_marker_max_clusters} largest clusters, these are "
        f"the genes with the largest lower and higher mean "
        f"Pearson residuals versus the largest cluster (cluster {reference_cluster}). "
        "The top 10 in each direction are shown. These are descriptive residual contrasts, "
        "not differential-expression significance tests."
    ))
    marker_clusters = sorted(markers["residual_cluster"].unique())
    figure, axes = plt.subplots(
        len(marker_clusters), 2, figsize=(13, 3.4 * len(marker_clusters)),
        squeeze=False,
    )
    for row_index, cluster in enumerate(marker_clusters):
        for column_index, direction in enumerate(("lower", "higher")):
            axis = axes[row_index, column_index]
            selected = markers.loc[
                (markers["residual_cluster"] == cluster)
                & (markers["direction"] == direction)
            ].sort_values("rank").head(10)
            axis.barh(
                selected["feature"].iloc[::-1],
                selected["mean_residual_difference_vs_reference"].iloc[::-1],
                color="#4c72b0", edgecolor="#38434a", linewidth=0.55,
            )
            axis.axvline(0, color="#555555", linewidth=0.8)
            axis.set_title(f"Cluster {cluster}: {direction} vs. largest cluster")
            axis.set_xlabel("Mean Pearson-residual difference")
    figure.tight_layout()
    display_collapsible_figure(
        figure, "Show residual-cluster marker contrasts"
    )


def plot_residual_sample_umaps(table, *, minimum_umap_samples):
    """Show residual clusters and sample metadata over the residual UMAP."""
    coordinates = table.dropna(subset=["umap_1", "umap_2"])
    if coordinates.empty:
        display(Markdown(
            f"_At least {minimum_umap_samples} samples are required for residual UMAPs; "
            "the cluster assignments are still available below._"
        ))
        return
    def plot_covariate_grid(
        covariates, title, *, categorical_covariates, figsize, open_by_default=False
    ):
        n_columns = 2
        n_rows = max(1, int(np.ceil(len(covariates) / n_columns)))
        figure_size = (figsize[0], max(figsize[1], figsize[1] * n_rows / 2))
        figure, axes = plt.subplots(
            n_rows, n_columns, figsize=figure_size, squeeze=False,
            layout="constrained",
        )
        for axis, covariate in zip(axes.flat, covariates):
            if covariate not in coordinates:
                axis.text(0.5, 0.5, "Not recorded", ha="center", va="center")
                axis.set_axis_off()
                continue
            values = coordinates.dropna(subset=[covariate]).copy()
            present = len(values)
            if values.empty:
                axis.text(0.5, 0.5, "No recorded values", ha="center", va="center")
                axis.set_axis_off()
                continue
            categorical = covariate in categorical_covariates
            if covariate == "batch_single_cell":
                plot_values = values.assign(_batch_label=values[covariate].astype(str))
                batch_order = sorted(plot_values["_batch_label"].unique())
                batch_palette = dict(zip(
                    batch_order,
                    sns.color_palette("husl", n_colors=len(batch_order)),
                ))
                sns.scatterplot(
                    data=plot_values, x="umap_1", y="umap_2", hue="_batch_label",
                    hue_order=batch_order, palette=batch_palette, s=50, alpha=0.85,
                    linewidth=0, ax=axis, legend=False,
                )
            elif categorical:
                plot_values = values.assign(_category_label=values[covariate].astype(str))
                category_order = sorted(plot_values["_category_label"].unique())
                category_palette = sns.color_palette(
                    "husl" if len(category_order) > 20 else "tab20",
                    n_colors=len(category_order),
                )
                sns.scatterplot(
                    data=plot_values, x="umap_1", y="umap_2", hue="_category_label",
                    hue_order=category_order, palette=dict(zip(category_order, category_palette)),
                    s=50, alpha=0.85, linewidth=0, ax=axis, legend="brief",
                )
                if axis.legend_ is not None:
                    handles, labels = axis.get_legend_handles_labels()
                    axis.legend(
                        handles, labels, loc="best",
                        ncol=2 if len(labels) > 6 else 1,
                        fontsize=6, title_fontsize=7, framealpha=0.9,
                        borderpad=0.35, labelspacing=0.25, handletextpad=0.35,
                    )
            else:
                sns.scatterplot(
                    data=values, x="umap_1", y="umap_2", hue=covariate,
                    palette="viridis", s=50, alpha=0.85, linewidth=0,
                    ax=axis, legend=False,
                )
                if axis.collections:
                    colorbar = figure.colorbar(
                        axis.collections[-1], ax=axis, orientation="horizontal",
                        pad=0.08, fraction=0.05, aspect=32,
                    )
                    colorbar.set_label(covariate.replace("_", " "), fontsize=7)
                    colorbar.ax.tick_params(labelsize=6)
            axis.set_title(f"{covariate.replace('_', ' ')} · {present}/{len(coordinates)} present")
            axis.set_xlabel("UMAP 1")
            axis.set_ylabel("UMAP 2")
        for axis in axes.flat[len(covariates):]:
            axis.set_axis_off()
        figure.suptitle(title)
        display_collapsible_figure(
            figure, title, open_by_default=open_by_default,
        )

    optional_covariates = [
        covariate for covariate in (
            "bmi", "cmv", "smoking_status", "ethnicity", "technology",
            "batch_single_cell"
        ) if covariate in coordinates
    ]
    coordinates["residual_cluster"] = coordinates["residual_cluster"].map(
        lambda cluster: "Noise" if pd.notna(cluster) and int(cluster) == -1
        else f"Cluster {int(cluster)}" if pd.notna(cluster) else np.nan
    )
    if optional_covariates or "residual_cluster" in coordinates:
        if "batch_single_cell" in optional_covariates:
            display(Markdown(
                "_Technical-batch colors are assigned per batch; its legend is omitted "
                "because there are many batch levels._"
            ))
        plot_covariate_grid(
            ["residual_cluster", *optional_covariates],
            "Residual clusters and non-model covariates over residual UMAP",
            categorical_covariates={
                "residual_cluster", "cmv", "smoking_status", "ethnicity", "technology",
                "batch_single_cell"
            },
            figsize=(18, 10.2),
            open_by_default=True,
        )
    plot_covariate_grid(
        ["sex", "age", "log10_total_counts", "study"],
        "Model covariates over residual UMAP",
        categorical_covariates={"sex", "study"},
        figsize=(18, 10.2),
    )


def plot_residual_pca_and_pc_umaps(pc_scores, coordinates):
    """Show residual PCA axes and their scores over the residual UMAP."""
    if pc_scores is None or not {"PC1", "PC2", "residual_cluster"}.issubset(pc_scores):
        display(Markdown("_Residual PCA scores are unavailable for these samples._"))
        return
    plot_data = coordinates.merge(
        pc_scores[["sample_key", "residual_cluster", "PC1", "PC2"]],
        on=["sample_key", "residual_cluster"],
        how="inner",
        validate="one_to_one",
        suffixes=("", "_pc"),
    ).dropna(subset=["PC1", "PC2"])
    if plot_data.empty:
        display(Markdown("_At least two residual PCs are needed for PCA plots._"))
        return

    figure, axes = plt.subplots(
        1, 3, figsize=(18, 5.2), layout="constrained",
    )
    sns.scatterplot(
        data=plot_data, x="PC1", y="PC2", hue="residual_cluster",
        palette="tab20", s=42, alpha=0.85, linewidth=0, ax=axes[0],
    )
    axes[0].set_title("Residual PCA: PC1 vs. PC2")
    if axes[0].legend_ is not None:
        handles, labels = axes[0].get_legend_handles_labels()
        axes[0].legend(
            handles, labels, title="Residual cluster", loc="best",
            ncol=2 if len(labels) > 6 else 1,
            fontsize=7, title_fontsize=8, framealpha=0.9,
            borderpad=0.35, labelspacing=0.25, handletextpad=0.35,
        )
    for axis, component in zip(axes[1:], ("PC1", "PC2")):
        values = plot_data.dropna(subset=["umap_1", "umap_2", component])
        points = axis.scatter(
            values["umap_1"], values["umap_2"], c=values[component],
            cmap="viridis", s=42, alpha=0.85, linewidths=0,
        )
        axis.set_title(f"Residual UMAP colored by {component}")
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        figure.colorbar(points, ax=axis, label=f"{component} score")
    display_collapsible_figure(figure, "Residual PCA and PC scores over UMAP")


def plot_residual_alternative_clusterings(annotations, settings):
    """Compare UMAP and PC-based residual cluster assignments."""
    summary = summarize_residual_clustering(annotations, settings)
    shown = summary[[
        "method", "n_clusters", "n_assigned", "n_noise", "noise_fraction", "cluster_sizes",
    ]].copy()
    shown["noise_fraction"] = shown["noise_fraction"].map(
        lambda value: f"{value:.1%}" if pd.notna(value) else "—"
    )
    display_collapsible_table(shown, "Residual clustering summary")

    plot_summary = summary.loc[
        summary["label_column"].isin([
            "residual_cluster_kmeans",
            "residual_cluster_hdbscan_umap",
            "residual_cluster_hdbscan_pcs",
        ])
    ]
    coordinates = annotations.dropna(subset=["umap_1", "umap_2"])
    n_rows = max(1, (len(plot_summary) + 1) // 2)
    figure, axes = plt.subplots(n_rows, 2, figsize=(14, 5.5 * n_rows), squeeze=False)
    for axis, row in zip(axes.flat, plot_summary.to_dict("records")):
        column = row["label_column"]
        if coordinates.empty or column not in coordinates:
            axis.text(0.5, 0.5, "UMAP or cluster labels unavailable", ha="center", va="center")
            axis.set_axis_off()
            continue
        values = coordinates.dropna(subset=[column]).copy()
        if values.empty:
            axis.text(0.5, 0.5, "Cluster labels unavailable", ha="center", va="center")
            axis.set_axis_off()
            continue
        labels = pd.to_numeric(values[column], errors="coerce")
        noise = labels.eq(-1)
        if noise.any():
            axis.scatter(
                values.loc[noise, "umap_1"], values.loc[noise, "umap_2"],
                color="#8a8a8a", s=35, alpha=0.75, linewidths=0,
            )
        clusters = sorted(labels.loc[~noise].dropna().astype(int).unique())
        for index, cluster in enumerate(clusters):
            selected = labels.eq(cluster)
            axis.scatter(
                values.loc[selected, "umap_1"], values.loc[selected, "umap_2"],
                color=plt.get_cmap("tab20")(index % 20), s=35, alpha=0.82, linewidths=0,
            )
        annotate_cluster_centroids(axis, values, column)
        axis.set_title(
            f"{row['method']}\n{row['n_clusters']} clusters; {row['n_noise']} noise"
        )
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
    for axis in axes.flat[len(plot_summary):]:
        axis.set_axis_off()
    figure.tight_layout()
    display_collapsible_figure(figure, "Residual cluster method comparison")



"""Existing plotting and model helpers for per-sample abundance reports."""


def has_age_span(data, *, min_samples=2):
    return len(data) >= min_samples and data["age"].nunique() >= 2


def facet_grid_shape(n_panels, *, max_columns=4):
    """Keep study panels compact while growing to any number of studies."""
    n_columns = min(max_columns, n_panels)
    return (n_panels + n_columns - 1) // n_columns, n_columns


def plot_study_linear_trends(
    data, *, fraction_column, title, ylabel, palette, show_fit_legend=True,
    open_by_default=True,
):
    """Show within-study trends with and without samples younger than 20."""
    if data.empty:
        display(Markdown(f"_No samples meet the denominator cutoff for {title}._"))
        return
    studies = sorted(data["study"].unique())
    n_rows, n_columns = facet_grid_shape(len(studies))
    figure, axes = plt.subplots(
        n_rows, n_columns, figsize=(4.8 * n_columns, 3.9 * n_rows),
        squeeze=False, sharex=True, sharey=True,
    )
    for axis, study in zip(axes.flat, studies):
        subset = data.loc[data["study"] == study]
        sns.scatterplot(
            data=subset, x="age", y=fraction_column, color=palette[study],
            s=24, alpha=0.7, edgecolor="none", ax=axis,
        )
        for include_all_ages, color, linestyle in [
            (False, "#111111", "-"),
            (True, "#777777", "--"),
        ]:
            fit_data = subset if include_all_ages else subset.loc[subset["age"] >= 20]
            if has_age_span(fit_data):
                sns.regplot(
                    data=fit_data, x="age", y=fraction_column, scatter=False, ci=None,
                    color=color,
                    line_kws={"linewidth": 2.4, "alpha": 0.95, "linestyle": linestyle},
                    ax=axis,
                )
        if not has_age_span(subset):
            axis.text(0.5, 0.08, "No age span for trend", transform=axis.transAxes,
                      ha="center", va="bottom", fontsize=9)
        axis.set_title(f"{study} (n={len(subset)})")
        axis.set_ylabel(ylabel)
    for axis in axes.flat[len(studies):]:
        axis.remove()
    figure.suptitle(title, y=0.99)
    if show_fit_legend:
        fit_handles = [
            plt.Line2D([], [], color="#111111", linewidth=2.4,
                       label="Linear fit: age ≥20"),
            plt.Line2D([], [], color="#777777", linewidth=2.4, linestyle="--",
                       label="Linear fit: all samples"),
        ]
        figure.legend(handles=fit_handles, loc="upper center", ncol=2,
                      bbox_to_anchor=(0.5, 0.94))
        figure.tight_layout(rect=(0, 0, 1, 0.87))
    else:
        figure.tight_layout(rect=(0, 0, 1, 0.94))
    display_collapsible_figure(
        figure, f"Within-study trends: {title}", open_by_default=open_by_default,
    )


def plot_study_fraction_curves(fraction_specs, *, palette):
    """Plot each study's descriptive lowess curve without a pooled fit."""
    figure, axes = plt.subplots(
        1, len(fraction_specs), figsize=(7.8 * len(fraction_specs), 5.5), squeeze=False,
    )
    for axis, (data, fraction_column, label) in zip(axes.flat, fraction_specs):
        plot_data = data[["study", "age", fraction_column]].copy()
        plot_data["study"] = plot_data["study"].astype(str)
        plot_data["age"] = pd.to_numeric(plot_data["age"], errors="coerce")
        plot_data[fraction_column] = pd.to_numeric(plot_data[fraction_column], errors="coerce")
        plot_data = plot_data.dropna().reset_index(drop=True)
        if plot_data.empty:
            axis.text(
                0.5, 0.5, "No samples available for descriptive curves",
                transform=axis.transAxes, ha="center", va="center",
            )
        else:
            studies = sorted(plot_data["study"].unique())
            sns.scatterplot(
                data=plot_data, x="age", y=fraction_column, hue="study", palette=palette,
                s=26, alpha=0.2, edgecolor="none", ax=axis,
            )
            for study in studies:
                study_data = plot_data.loc[plot_data["study"] == study]
                if has_age_span(study_data):
                    sns.regplot(
                        data=study_data, x="age", y=fraction_column, scatter=False,
                        lowess=True, ci=None, color=palette[study],
                        line_kws={"linewidth": 2.4, "alpha": 0.95}, ax=axis,
                    )
            axis.legend(title="Study", loc="best", fontsize=8)
        axis.set_title(label)
        axis.set_xlabel("Age (years)")
        axis.set_ylabel(label)
    figure.suptitle("Study-specific descriptive fraction curves", y=0.98)
    figure.text(
        0.5, 0.015,
        "Points are study × sample (alpha = 0.2). Each coloured curve is a separate study lowess fit; no combined model is fitted.",
        ha="center", fontsize=9,
    )
    figure.tight_layout(rect=(0, 0.06, 1, 0.92))
    display_collapsible_figure(figure, "Show study-specific fraction curves")


def _bootstrap_binomial_residual_sd(design, denominators, fitted_probabilities, *, rng, n_bootstrap):
    """Estimate post-model residual SD expected from binomial cell sampling alone.

    The fitted probabilities and design are fixed inputs.  Thus this is a
    conditional measurement-noise reference: only finite cell sampling varies
    between replicates, while the observed covariate structure is retained.
    """
    n_samples, n_terms = design.shape
    if n_samples <= n_terms:
        return np.nan
    denominators = np.asarray(denominators, dtype=np.int64)
    fitted_probabilities = np.asarray(fitted_probabilities, dtype=float)
    simulated = rng.binomial(
        denominators, fitted_probabilities, size=(n_bootstrap, n_samples),
    ) / denominators
    orthonormal_design, _ = np.linalg.qr(design, mode="reduced")
    # Refit the same fixed design because the observed quantity is also a
    # post-covariate-model residual SD.  This removes chance alignment between
    # simulated binomial noise and the covariates, just as in the observed fit.
    residuals = simulated - (simulated @ orthonormal_design) @ orthonormal_design.T
    return float(np.sqrt(np.mean(np.sum(residuals**2, axis=1) / (n_samples - n_terms))))


def fit_parent_fraction_models(
    fractions, adata, *, fraction_column, bootstrap_replicates, bootstrap_seed,
):
    """Fit adult-only, study-specific OLS models and residual diagnostics."""
    candidate_covariates = ["age", "sex", "bmi", "cmv", "smoking_status"]
    present = [column for column in candidate_covariates if column in adata.obs]
    metadata_columns = [column for column in present if column != "age"]
    metadata = adata.obs.groupby(["study", "sample"], observed=True)[metadata_columns].first().reset_index()
    # Repeated cell metadata must describe one sample consistently.
    for column in present:
        n_unique = adata.obs.groupby(["study", "sample"], observed=True)[column].nunique(dropna=False)
        if (n_unique > 1).any():
            raise ValueError(f"{column} must be constant within each study × sample")
    data = fractions.merge(metadata, on=["study", "sample"], how="left", validate="one_to_one")
    data["age"] = pd.to_numeric(data["age"], errors="coerce")
    if "cmv" in present:
        data["cmv"] = normalize_cmv_status(data["cmv"])
    if "smoking_status" in present:
        data["smoking_status"] = normalize_smoking_status(data["smoking_status"])
    data = data.loc[data["age"] >= 20].copy()
    results = []
    residual_diagnostics = []
    unavailable = {
        "", "not_provided", "heterogeneous during bulk", "unknown", "nan", "none",
    }
    bootstrap_rng = np.random.default_rng(bootstrap_seed)

    def valid_values(values, column):
        if column == "smoking_status":
            return normalize_smoking_status(values).notna()
        if column in {"age", "bmi"}:
            numeric = pd.to_numeric(values, errors="coerce")
            return numeric.notna() & np.isfinite(numeric)
        normalized = values.astype("string").str.strip().str.lower()
        return normalized.notna() & ~normalized.isin(unavailable)

    for study, subset in data.groupby("study", sort=True, observed=True):
        eligible = subset.copy()
        selected = []
        for column in candidate_covariates:
            if column not in present:
                continue
            valid = valid_values(eligible[column], column)
            coverage = float(valid.mean()) if len(eligible) else np.nan
            if coverage >= 0.95:
                if column in {"sex", "cmv", "smoking_status"}:
                    levels = categorical_levels(column, eligible.loc[valid, column])
                    variable = len(levels) >= 2
                else:
                    variable = pd.to_numeric(eligible.loc[valid, column], errors="coerce").nunique() >= 2
                if variable:
                    selected.append(column)

        if not selected:
            continue
        valid_rows = pd.Series(True, index=eligible.index)
        for column in selected:
            valid_rows &= valid_values(eligible[column], column)
        model_data = eligible.loc[valid_rows].copy()
        if model_data.empty:
            continue
        pieces = [pd.Series(1.0, index=model_data.index, name="intercept")]
        term_covariates = {"intercept": "intercept"}
        for column in selected:
            if column in {"age", "bmi"}:
                pieces.append(pd.to_numeric(model_data[column], errors="coerce").rename(column))
                term_covariates[column] = column
            else:
                values = model_data[column].astype("string").str.strip().str.lower()
                levels = categorical_levels(column, values)
                values = pd.Series(pd.Categorical(values, categories=levels), index=model_data.index)
                dummies = pd.get_dummies(values, prefix=column, drop_first=True, dtype=float)
                for term in dummies:
                    pieces.append(dummies[term].rename(term))
                    term_covariates[term] = column
        design = pd.concat(pieces, axis=1).astype(float)
        if len(model_data) <= design.shape[1] or np.linalg.matrix_rank(design.to_numpy()) < design.shape[1]:
            continue
        design_array = design.to_numpy()
        fit = sm.OLS(model_data[fraction_column].to_numpy(dtype=float), design_array).fit()
        confidence = fit.conf_int()
        mean_fraction = float(model_data[fraction_column].mean())
        unbounded_fitted_probabilities = np.asarray(fit.predict(design_array), dtype=float)
        n_clipped_probabilities = int(
            ((unbounded_fitted_probabilities < 0) | (unbounded_fitted_probabilities > 1)).sum()
        )
        # OLS is retained as the reported fraction model.  Its fitted values
        # provide donor-specific conditional means for the binomial null; clip
        # only because binomial probabilities must lie in [0, 1].  The count is
        # published so any material boundary correction is visible.
        fitted_probabilities = np.clip(unbounded_fitted_probabilities, 0, 1)
        residual_diagnostics.append({
            "study": str(study),
            "residual_sd": float(np.sqrt(fit.mse_resid)),
            "binomial_sampling_sd": _bootstrap_binomial_residual_sd(
                design_array,
                model_data["denominator_cells"].to_numpy(),
                fitted_probabilities,
                rng=bootstrap_rng,
                n_bootstrap=bootstrap_replicates,
            ),
            "n_samples": len(model_data),
            "mean_fraction": mean_fraction,
            "n_clipped_binomial_probabilities": n_clipped_probabilities,
            "fraction_clipped_binomial_probabilities": n_clipped_probabilities / len(model_data),
            "covariates": ";".join(selected),
            "model_terms": ";".join(design.columns.drop("intercept")),
        })
        for index, term in enumerate(design.columns):
            if term == "intercept":
                continue
            covariate = term_covariates[term]
            if covariate in {"sex", "cmv", "smoking_status"}:
                levels = categorical_levels(
                    covariate,
                    model_data[covariate].astype("string").str.strip().str.lower(),
                )
                contrast = f"{term.removeprefix(f'{covariate}_')} vs {levels[0]}"
            else:
                contrast = "per decade" if covariate == "age" else "per BMI unit"
            effect_scale = 10 if covariate == "age" else 1
            results.append({
                "study": str(study), "covariate": covariate, "term": term,
                "contrast": contrast, "estimate": fit.params[index] * effect_scale,
                "ci_low": confidence[index, 0] * effect_scale,
                "ci_high": confidence[index, 1] * effect_scale,
                "n_samples": len(model_data),
            })
    coefficient_columns = [
        "study", "covariate", "term", "contrast", "estimate", "ci_low", "ci_high",
        "n_samples",
    ]
    diagnostic_columns = [
        "study", "residual_sd", "binomial_sampling_sd", "n_samples", "mean_fraction",
        "n_clipped_binomial_probabilities", "fraction_clipped_binomial_probabilities",
        "covariates", "model_terms",
    ]
    # Keep the report-facing schema stable even when no study has enough adult
    # samples or covariate variation for an estimable model.
    return (
        pd.DataFrame(results, columns=coefficient_columns),
        pd.DataFrame(residual_diagnostics, columns=diagnostic_columns),
    )


def plot_parent_fraction_coefficients(
    coefficients, residual_diagnostics, *, palette, cell_type, parent,
):
    """Plot adjusted coefficients plus observed and binomial residual variation."""
    if coefficients.empty:
        display(Markdown("_No estimable adult-only study models for this fraction._"))
        return
    covariate_order = [name for name in ["age", "sex", "bmi", "cmv", "smoking_status"]
                       if name in coefficients["covariate"].unique()]
    studies = sorted(coefficients["study"].unique())
    residual_label = "residual_sd"
    fraction_to_percentage_points = 100.0
    coefficient_rows = []
    for covariate in covariate_order:
        subset = coefficients.loc[coefficients["covariate"] == covariate]
        for contrast in sorted(subset["contrast"].astype(str).unique()):
            if covariate in {"sex", "cmv", "smoking_status"}:
                label = categorical_contrast_label(
                    covariate, contrast.split(" vs ", maxsplit=1)[0],
                    contrast.split(" vs ", maxsplit=1)[1],
                )
            elif covariate == "age":
                label = "AGE (per decade)"
            else:
                label = "BMI (per BMI unit)"
            coefficient_rows.append((covariate, contrast, label))
    y_positions = {
        (covariate, contrast): index
        for index, (covariate, contrast, _) in enumerate(coefficient_rows)
    }
    residual_y = len(coefficient_rows)
    y_positions[residual_label] = residual_y

    # Keep broad confidence intervals from setting the scale. Estimates and
    # residual variation define the visible range; truncated CI ends are marked
    # explicitly at the panel boundary below.
    scale_values = [0.0]
    scale_values.extend(
        pd.to_numeric(coefficients["estimate"], errors="coerce")
        * fraction_to_percentage_points
    )
    for column in ["residual_sd", "binomial_sampling_sd"]:
        if column in residual_diagnostics:
            scale_values.extend(
                pd.to_numeric(residual_diagnostics[column], errors="coerce")
                * fraction_to_percentage_points
            )
    scale_values = np.asarray(scale_values, dtype=float)
    scale_values = scale_values[np.isfinite(scale_values)]
    data_min, data_max = float(scale_values.min()), float(scale_values.max())
    padding = max((data_max - data_min) * 0.08, max(abs(data_min), abs(data_max), 1.0) * 0.02)
    x_limits = (data_min - padding, data_max + padding)

    figure, axis = plt.subplots(figsize=(10, max(4.8, 1.15 * (len(coefficient_rows) + 1))))
    clipped_ci_ends = []
    for study_index, study in enumerate(studies):
        subset = coefficients.loc[coefficients["study"] == study]
        for _, row in subset.iterrows():
            offset = (study_index - (len(studies) - 1) / 2) * 0.055
            y = y_positions[(row["covariate"], row["contrast"])] + offset
            estimate = row["estimate"] * fraction_to_percentage_points
            ci_low = row["ci_low"] * fraction_to_percentage_points
            ci_high = row["ci_high"] * fraction_to_percentage_points
            axis.errorbar(
                estimate, y, xerr=[[estimate - ci_low], [ci_high - estimate]],
                fmt="o", color=palette[study], capsize=3, markersize=5,
            )
            if ci_low < x_limits[0]:
                clipped_ci_ends.append((x_limits[0], y, "<", palette[study]))
            if ci_high > x_limits[1]:
                clipped_ci_ends.append((x_limits[1], y, ">", palette[study]))
    for study_index, study in enumerate(studies):
        subset = residual_diagnostics.loc[residual_diagnostics["study"] == study]
        if subset.empty:
            continue
        row = subset.iloc[0]
        # Interpretation check: assess whether residual variation is too large
        # for the coefficient estimates above to support useful interpretation.
        offset = (study_index - (len(studies) - 1) / 2) * 0.055
        y = residual_y + offset
        axis.plot(
            row["residual_sd"] * fraction_to_percentage_points,
            y, "o", color=palette[study], markersize=6,
        )
        axis.plot(
            row["binomial_sampling_sd"] * fraction_to_percentage_points,
            y, marker="D", linestyle="", markersize=5,
            markerfacecolor="none", markeredgecolor=palette[study], markeredgewidth=1.25,
        )
    axis.axvline(0, color="#333333", linewidth=1, linestyle="--")
    axis.set_xlim(*x_limits)
    for boundary, y, direction, color in clipped_ci_ends:
        axis.plot(
            boundary, y, marker=direction, linestyle="", color=color,
            markersize=5, clip_on=False, zorder=5,
        )
    axis.set_yticks(
        list(range(residual_y + 1)),
        [*(label for _, _, label in coefficient_rows), "RESIDUAL SD"],
    )
    axis.set_xlabel("Fraction change (percentage points; effects with 95% CI)")
    axis.set_ylabel("Covariate / diagnostic")
    axis.invert_yaxis()
    study_legend = axis.legend(
        handles=[plt.Line2D([], [], marker="o", linestyle="", color=palette[study], label=study)
                 for study in studies],
        title="Study", bbox_to_anchor=(1.02, 1), loc="upper left",
    )
    axis.add_artist(study_legend)
    axis.legend(
        handles=[
            plt.Line2D([], [], marker="o", linestyle="", color="#333333", label="Observed residual SD"),
            plt.Line2D([], [], marker="D", linestyle="", color="#333333", markerfacecolor="none",
                       label="Simulated sampling SD"),
        ],
        title="Residual diagnostic", bbox_to_anchor=(1.02, 0.53), loc="upper left",
    )
    axis.set_title(
        f"Adult-only covariate effects on fraction of {cell_type}\nwithin {parent}"
    )
    figure.tight_layout()
    plt.show()



# %% [markdown]
# <a id="sample-fractions"></a>
# ## Sample-Level Fraction Across Age
#
# Each point is one study × sample. The left panel uses all retained PBMCs as
# denominator; the right panel uses the manually configured AIFI-L1 parent
# derived from every cell's merged L2 call, never from a separate L1 classifier.
# Samples with fewer than the configured denominator-cell cutoff are excluded.
# Facets show descriptive within-study linear fits for all samples and for age
# ≥20. The all-study view places all studies in one panel per denominator, with
# semi-transparent points and a separate lowess curve per study. It is purely
# descriptive: no combined model or pooled age effect is fitted.

# %%
min_fraction_denominator = pipeline["differential_expression"]["abundance"]["min_fraction_denominator_cells"]
fraction_specs = []
for denominator, fraction_column, title, ylabel in [
    ("n_cells_in_sample", "fraction_of_retained_pbmc", "Fraction of all retained PBMCs", "Fraction of retained PBMCs"),
    ("n_cells_in_sample_l1_parent", "fraction_within_aifi_l1_parent", None, None),
]:
    fraction_specs.append((
        sample_cell_type_fractions(
            adata, denominator=denominator, fraction_name=fraction_column,
            min_denominator_cells=min_fraction_denominator,
            count_column="n_cells",
        ),
        fraction_column,
        title,
        ylabel,
    ))
parent = adata.obs["aifi_l1_parent_for_l2"].iloc[0]
fraction_specs[1] = (
    fraction_specs[1][0], fraction_specs[1][1], f"Fraction within {parent}", f"Fraction of {parent}"
)
all_sample_count = adata.obs[["study", "sample"]].drop_duplicates().shape[0]
all_fraction_samples = [all_sample_count] * len(fraction_specs)
display_collapsible_table(pd.DataFrame({
    "fraction": [spec[2] for spec in fraction_specs],
    "eligible_samples": [len(spec[0]) for spec in fraction_specs],
    "excluded_below_denominator_cutoff": [
        total - len(spec[0]) for total, spec in zip(all_fraction_samples, fraction_specs)
    ],
    "minimum_denominator_cells": min_fraction_denominator,
}), "Show sample eligibility by fraction denominator")

studies = sorted(pd.concat([spec[0] for spec in fraction_specs])["study"].unique())
palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))

# %%
display(Markdown(
    f"### Adjusted Fraction Effects: {cell_type_name} Within {parent}"
))

# %% [markdown]
# Each study has a separate ordinary least-squares model with one row per
# eligible sample aged **20 years or older**. The response is this cell type's
# fraction of its derived AIFI-L1 parent compartment. Age and BMI are numeric;
# sex, CMV, and smoking status are categorical treatment-coded terms. A
# covariate is included
# when it is present and observed in at least 95% of adult samples in that
# study. Models use complete cases for the included covariates, and categorical
# contrasts use female as the sex reference and the no CMV and no smoking-status
# categories as references (negative CMV labels are normalized to no). Age coefficients and
# intervals are
# scaled to a 10-year change. The bottom row compares each model's
# residual SD with a conditional parametric bootstrap of binomial cell-sampling
# noise. For donor *i*, the reported OLS model supplies its fitted fraction
# \(p_i\), clipped to [0, 1], and its observed parent-cell count is \(n_i\).
# Each replicate draws \(K_i \sim Binomial(n_i, p_i)\), analyzes
# \(K_i / n_i\) with the same covariate design, and records the residual SD.
# This holds the fitted covariate structure fixed and varies only finite
# cell-count sampling; it is therefore directly comparable with the observed
# post-model residual SD. The number and fraction of clipped OLS predictions
# are published in the diagnostics TSV. More bootstrap replicates improve Monte-Carlo
# precision, but do not change the target quantity. The residuals need review
# before interpreting the other terms: if they are too large relative to their
# scale, those effects may not support a reasonable interpretation. In the forest
# plot, “Observed residual SD” is the fitted OLS residual SD; “Simulated sampling
# SD” is the mean residual SD across the conditional binomial bootstrap replicates.
# The response is a fraction on the 0–1 scale. The plot multiplies effects, their
# 95% confidence limits, and residual SDs by 100 for display in percentage points
# (a 0.01 fraction change is 1 percentage point); these are absolute percentage-
# point changes, not relative percent changes.

# %%
display(Markdown("**Model sample restriction:** adults aged ≥20 years only."))
parent_fraction_data = fraction_specs[1][0]
parent_fraction_column = fraction_specs[1][1]
fraction_model_settings = pipeline["differential_expression"]["abundance"]["fraction_model"]
fraction_coefficients, fraction_residual_diagnostics = fit_parent_fraction_models(
    parent_fraction_data, adata, fraction_column=parent_fraction_column,
    bootstrap_replicates=fraction_model_settings["binomial_bootstrap_replicates"],
    bootstrap_seed=fraction_model_settings["binomial_bootstrap_seed"],
)
fraction_residual_diagnostics.insert(0, "cell_type", cell_type_name)
fraction_residual_diagnostics.insert(1, "aifi_l1_parent", parent)
fraction_residual_diagnostics.insert(2, "response", parent_fraction_column)
fraction_residual_diagnostics["binomial_bootstrap_replicates"] = (
    fraction_model_settings["binomial_bootstrap_replicates"]
)
fraction_residual_diagnostics["binomial_bootstrap_seed"] = fraction_model_settings["binomial_bootstrap_seed"]
fraction_diagnostics_path = output_dir / "fraction_model_diagnostics.tsv"
fraction_residual_diagnostics.to_csv(fraction_diagnostics_path, sep="\t", index=False)
display(Markdown("#### Fraction-Model Residual and Sampling Diagnostics"))
fraction_diagnostics_display = fraction_residual_diagnostics.rename(columns={
    "cell_type": "Cell type",
    "aifi_l1_parent": "L1 parent",
    "response": "Response",
    "study": "Study",
    "residual_sd": "Residual SD",
    "binomial_sampling_sd": "Binomial sampling SD",
    "n_samples": "Adult samples",
    "mean_fraction": "Mean fraction",
    "n_clipped_binomial_probabilities": "Clipped predictions (n)",
    "fraction_clipped_binomial_probabilities": "Clipped prediction fraction",
    "covariates": "Covariates",
    "model_terms": "Model terms",
    "binomial_bootstrap_replicates": "Bootstrap replicates",
    "binomial_bootstrap_seed": "Bootstrap seed",
})
display_collapsible_table(
    fraction_diagnostics_display, "Show full per-study residual diagnostics",
)

age_effects = fraction_coefficients.loc[fraction_coefficients["covariate"] == "age"]
n_positive_age_effects = int((age_effects["estimate"] > 0).sum())
n_negative_age_effects = int((age_effects["estimate"] < 0).sum())
if age_effects.empty:
    age_direction = "No estimable age coefficient"
elif len(age_effects) == 1:
    age_direction = "1 positive estimate" if n_positive_age_effects else "1 negative estimate"
elif n_positive_age_effects and n_negative_age_effects:
    age_direction = f"Mixed: {n_positive_age_effects} positive, {n_negative_age_effects} negative"
else:
    age_direction = (
        f"{n_positive_age_effects} positive estimates"
        if n_positive_age_effects else f"{n_negative_age_effects} negative estimates"
    )
sampling_ratio = (
    fraction_residual_diagnostics["residual_sd"]
    / fraction_residual_diagnostics["binomial_sampling_sd"].replace(0, np.nan)
)
sampling_ratio = sampling_ratio.replace([np.inf, -np.inf], np.nan).dropna()
residual_summary = (
    f"Median {sampling_ratio.median():.1f}× the conditional sampling reference"
    if not sampling_ratio.empty else "Sampling comparison unavailable"
)
model_covariates = sorted({
    covariate
    for terms in fraction_residual_diagnostics["covariates"].dropna()
    for covariate in str(terms).split(";") if covariate
})
display(HTML(
    "<section id=\"fraction-model-evidence\" class=\"report-evidence\">"
    "<h3>Fraction-Model Evidence at a Glance</h3>"
    "<p>These are study-specific adult-sample OLS fits, not a pooled effect estimate. "
    "Use them to assess direction, precision, and the residual diagnostic before interpreting individual coefficients.</p>"
    "<div class=\"report-evidence-grid\">"
    f"<div class=\"report-evidence-item\"><p class=\"report-evidence-label\">Estimable study models</p><p class=\"report-evidence-value\">{len(fraction_residual_diagnostics)}</p></div>"
    f"<div class=\"report-evidence-item\"><p class=\"report-evidence-label\">Age-effect direction</p><p class=\"report-evidence-value\">{escape(age_direction)}</p></div>"
    f"<div class=\"report-evidence-item\"><p class=\"report-evidence-label\">Residual check</p><p class=\"report-evidence-value\">{escape(residual_summary)}</p></div>"
    f"<div class=\"report-evidence-item\"><p class=\"report-evidence-label\">Covariates represented</p><p class=\"report-evidence-value\">{escape(', '.join(model_covariates) or 'None')}</p></div>"
    "</div></section>"
))
plot_parent_fraction_coefficients(
    fraction_coefficients, fraction_residual_diagnostics, palette=palette,
    cell_type=cell_type_name, parent=parent,
)
display(Markdown(
    "_How to read this plot: each colour is a study and intervals are its 95% OLS confidence "
    "intervals. Fraction effects and both residual SD markers are shown in percentage points. "
    "The filled dot on the residual row is the observed residual SD; the open diamond is the "
    "simulated sampling SD described above. Compare studies for consistency and precision, "
    "rather than treating their collection as a pooled estimate. Arrowheads mark 95% CIs "
    "extending beyond the displayed range._"
))

display(Markdown(
    "<a id=\"additional-fraction-trends\"></a>\n### Additional Descriptive Fraction Trends"
))
display(Markdown(
    "These age-trend views are secondary to the adult-only adjusted estimates above. "
    "The within-L1 within-study panel opens on load; the all-retained-PBMC panel starts "
    "closed. Both show all-age and adult-only linear fits. "
    "The all-study panels show one descriptive lowess curve per study."
))
for index, (fractions, fraction_column, title, ylabel) in enumerate(fraction_specs):
    plot_study_linear_trends(
        fractions, fraction_column=fraction_column, title=title, ylabel=ylabel, palette=palette,
        show_fit_legend=index == 0,
        open_by_default=index != 0,
    )
plot_study_fraction_curves(
    [(fractions, fraction_column, ylabel) for fractions, fraction_column, _, ylabel in fraction_specs],
    palette=palette,
)


# %%
if de_dir:
    de_dir = Path(de_dir)
    indexed_result_paths = {
        (de_dir / relative_path).resolve()
        for model in de_run_metadata.get("models", [])
        if model.get("status") == "complete"
        for relative_path in model.get("results_by_covariate", {}).values()
    }
    if "models" in de_run_metadata:
        per_study_files = sorted(
            path for path in (de_dir / "per_study").rglob("*.csv")
            if path.resolve() in indexed_result_paths
        )
        merged_file = de_dir / "merged.csv"
        if merged_file.resolve() not in indexed_result_paths:
            merged_file = None
        combined_files = sorted(
            path for path in (de_dir / "combined").glob("*.csv")
            if path.resolve() in indexed_result_paths
        )
    else:
        per_study_files = sorted((de_dir / "per_study").rglob("*.csv"))
        merged_file = de_dir / "merged.csv"
        combined_files = sorted((de_dir / "combined").glob("*.csv"))
    per_study_results = []
    for result_path in per_study_files:
        result = pd.read_csv(result_path)
        if "covariate" not in result:
            result["covariate"] = "age"
        if "contrast" not in result:
            result["contrast"] = "per year"
        if "study_sample_counts" not in result:
            result["study_sample_counts"] = "{}"
        per_study_results.append(result)
    has_per_study = bool(per_study_results)
    combined_results_by_covariate = {}
    age_trajectory_results = None
    if merged_file is not None and merged_file.is_file():
        result = pd.read_csv(merged_file)
        if "covariate" not in result:
            result["covariate"] = "age"
        if "contrast" not in result:
            result["contrast"] = "per year"
        combined_results_by_covariate["age"] = result
    for result_path in combined_files:
        result = pd.read_csv(result_path)
        if result_path.stem == "age_bin":
            age_trajectory_results = result
            continue
        if "covariate" not in result:
            result["covariate"] = result_path.stem
        if "contrast" not in result:
            result["contrast"] = "per BMI unit" if result_path.stem == "bmi" else result_path.stem
        combined_results_by_covariate[result_path.stem] = result
    has_combined = bool(combined_results_by_covariate)
    if has_per_study or has_combined or de_run_metadata:
        display(Markdown("<a id=\"differential-expression\"></a>\n## Differential Gene Expression"))
        display(Markdown(
            "Pseudobulk models aggregate counts at the study × sample × cell-type level. "
            "They provide the sample-level complement to the descriptive single-cell views above."
        ))
        display(Markdown(
            "**Gene universes:** per-study fits use genes observed in that study; each combined "
            "fit uses the intersection across its included studies, excluding outer-join zeros for "
            "study-absent genes. The cross-cell-type expression atlas uses the all-study intersection."
        ))
        if de_run_metadata.get("analysis_mode") == "test_only":
            display(Markdown(
                "> **TEST OUTPUT — NOT FOR BIOLOGICAL INTERPRETATION.** "
                + de_run_metadata["interpretation_warning"]
            ))
        if not has_per_study and not has_combined:
            display(Markdown(
                "_No estimable PyDESeq2 result files were produced for this cell type._"
            ))
    age_trajectory_fit = next(
        (model for model in de_run_metadata.get("models", [])
         if model.get("purpose") == "age_trajectory"),
        {},
    )
    age_trajectory_info = de_run_metadata.get("age_trajectory", {})
    if has_per_study or has_combined or age_trajectory_fit or age_trajectory_results is not None:
        de_alpha = pipeline["differential_expression"]["alpha"]
        per_study_by_covariate = {}
        for result in per_study_results:
            for covariate, subset in result.groupby("covariate", observed=True):
                per_study_by_covariate.setdefault(str(covariate), []).append(subset)
        age_per_study = (
            pd.concat(per_study_by_covariate.get("age", []), ignore_index=True)
            if per_study_by_covariate.get("age") else pd.DataFrame()
        )
        age_study_results = (
            [group for _, group in age_per_study.groupby("study", observed=True)]
            if not age_per_study.empty else []
        )
        has_age_per_study = bool(age_study_results)
        if has_age_per_study:
            gene_universes = [set(result["gene"].astype(str)) for result in age_study_results]
            common_per_study_genes = set.intersection(*gene_universes)
            per_study_significant = age_per_study.loc[
                pd.to_numeric(age_per_study["padj"], errors="coerce") < de_alpha
            ].copy()
            per_study_hits = per_study_significant.loc[
                per_study_significant["gene"].astype(str).isin(common_per_study_genes)
            ].copy()
            recurrence = (
                per_study_hits.groupby("gene", as_index=False)
                .agg(
                    studies_associated=("study", "nunique"),
                    studies=("study", lambda values: ", ".join(sorted(set(values)))),
                    maximum_absolute_log2_fold_change=(
                        "log2FoldChange", signed_value_at_largest_absolute_magnitude
                    ),
                )
                .assign(_absolute_effect=lambda frame: frame[
                    "maximum_absolute_log2_fold_change"
                ].abs())
                .sort_values(
                    ["studies_associated", "_absolute_effect", "gene"],
                    ascending=[False, False, True],
                )
                .drop(columns="_absolute_effect")
            )
        else:
            recurrence = pd.DataFrame()

        diagnostic_rows = []
        for covariate, model_results in per_study_by_covariate.items():
            covariate_results = pd.concat(model_results, ignore_index=True)
            study_gene_sets = [
                set(group["gene"].astype(str))
                for _, group in covariate_results.groupby("study", observed=True)
            ]
            if not study_gene_sets:
                continue
            common_genes = set.intersection(*study_gene_sets)
            for (study, contrast), study_results in covariate_results.groupby(
                ["study", "contrast"], sort=True, observed=True
            ):
                significant = study_results.loc[
                    pd.to_numeric(study_results["padj"], errors="coerce") < de_alpha
                ]
                in_intersection = significant.loc[
                    significant["gene"].astype(str).isin(common_genes), "gene"
                ].nunique()
                total_significant = int(significant["gene"].nunique())
                outside_intersection = total_significant - int(in_intersection)
                diagnostic_rows.append({
                    "covariate": covariate,
                    "contrast": contrast,
                    "study": study,
                    "design": (
                        study_results["design"].iloc[0]
                        if "design" in study_results else "—"
                    ),
                    "complete-case samples": (
                        int(study_results["n_samples_model"].iloc[0])
                        if "n_samples_model" in study_results else "—"
                    ),
                    "excluded for missing design values": (
                        int(study_results["n_samples_excluded_missing_design_covariates"].iloc[0])
                        if "n_samples_excluded_missing_design_covariates" in study_results else "—"
                    ),
                    "significant genes": total_significant,
                    "significant genes in study intersection": int(in_intersection),
                    "significant genes outside intersection": outside_intersection,
                    "outside intersection (%)": (
                        100 * outside_intersection / total_significant if total_significant else 0.0
                    ),
                })
        intersection_diagnostics = pd.DataFrame(diagnostic_rows)

        summary_rows = []
        if has_age_per_study:
            summary_rows.append(
                {
                    "model": "Per-study age significant-gene union",
                    "FDR-significant age-associated genes": int(
                        per_study_significant["gene"].nunique()
                    ),
                }
            )
            summary_rows.append(
                {
                    "model": "Per-study age significant genes in all-study intersection",
                    "FDR-significant age-associated genes": int(per_study_hits["gene"].nunique()),
                }
            )
        combined_fit_records = [
            model for model in de_run_metadata.get("models", [])
            if model.get("model") in {"merged", "combined"}
        ]
        age_diagnostics_path = de_dir / "age_model_diagnostics.csv"
        age_combined_results = combined_results_by_covariate.get("age")
        if age_diagnostics_path.is_file() and age_combined_results is not None:
            age_diagnostics = pd.read_csv(age_diagnostics_path)
            expected_columns = {
                "study", "age", "sex", "counts_in_age_gene_intersection",
                "genes_in_age_intersection",
            }
            if expected_columns.issubset(age_diagnostics):
                display(Markdown(
                    "<a id=\"shared-age-model-diagnostics\"></a>\n### Shared Age-Model Diagnostics"
                ))
                display(Markdown(
                    "This single figure checks study-site/age support, the relation between age and "
                    "the raw pseudobulk counts used in the shared age gene intersection, and "
                    "the unadjusted age-model p-value calibration. Q-Q points are the raw "
                    "p-values (`pvalue`, before multiple-testing correction) from the combined "
                    "age model fitted across all eligible studies for this cell type, adjusted "
                    "for study site, sex, age, and log10(total_counts), and using the model's shared "
                    "gene intersection. "
                    "The dashed y=x line is the null reference; "
                    "the axes retain independent scales so departures remain readable. This is "
                    "descriptive: it does not diagnose a specific gene or replace model checks."
                ))
                figure, axes = plt.subplots(1, 3, figsize=(18, 4.8))
                support_column = "study_site" if "study_site" in age_diagnostics else "study"
                sns.stripplot(
                    data=age_diagnostics, x=support_column, y="age", hue="sex", dodge=True,
                    jitter=0.18, alpha=0.8, ax=axes[0],
                )
                axes[0].set_xlabel("Study site" if support_column == "study_site" else "Study")
                axes[0].set_ylabel("Age (years)")
                axes[0].set_title("Age support by study site")
                axes[0].tick_params(axis="x", rotation=75, labelsize=7)
                axes[0].legend(title="Sex", fontsize=8, title_fontsize=8)

                plot_samples = age_diagnostics.loc[
                    age_diagnostics["counts_in_age_gene_intersection"] > 0
                ].copy()
                plot_samples["log10_counts"] = np.log10(
                    plot_samples["counts_in_age_gene_intersection"]
                )
                sns.scatterplot(
                    data=plot_samples, x="age", y="log10_counts", hue=support_column, style="sex",
                    s=42, alpha=0.85, ax=axes[1],
                )
                axes[1].set_xlabel("Age (years)")
                axes[1].set_ylabel("log10 raw counts in age gene intersection")
                axes[1].set_title("Age and pseudobulk depth")
                axes[1].legend(fontsize=7, title_fontsize=8)

                pvalue_source = (
                    age_combined_results["pvalue"]
                    if "pvalue" in age_combined_results else pd.Series(dtype=float)
                )
                pvalues = pd.to_numeric(pvalue_source, errors="coerce")
                pvalues = pvalues.loc[(pvalues > 0) & (pvalues <= 1)].sort_values().to_numpy()
                if len(pvalues):
                    expected = (np.arange(1, len(pvalues) + 1) - 0.5) / len(pvalues)
                    axes[2].scatter(
                        -np.log10(expected), -np.log10(pvalues), s=9, alpha=0.65,
                        color="#4c72b0", linewidths=0,
                    )
                    axes[2].axline(
                        (0, 0), slope=1, color="black", linestyle="--",
                        linewidth=1.5, zorder=4, label="Null (y = x)",
                    )
                    axes[2].set_xlabel("Expected −log10(p)")
                    axes[2].set_ylabel("Observed −log10(p)")
                    axes[2].set_title("Age-model p-value QQ plot")
                    axes[2].legend(fontsize=8)
                else:
                    axes[2].text(
                        0.5, 0.5, "No finite unadjusted p-values", ha="center", va="center",
                        transform=axes[2].transAxes,
                    )
                    axes[2].set_axis_off()
                figure.tight_layout()
                plt.show()

        display(Markdown("### Covariate Effects"))
        display(Markdown(
            f"Each combined fit estimates one covariate effect from eligible studies using "
            f"the complete cases in its displayed design. **FDR threshold:** {de_alpha:g}. "
            "Model details and sample counts are listed below each volcano plot."
        ))
        if not combined_results_by_covariate:
            display(Markdown("_No combined covariate results are available for volcano plots._"))
        covariate_order = {
            "age": 0,
            "sex": 1,
            "bmi": 2,
            "cmv": 3,
            "smoking_status": 4,
            "log10_total_counts": 5,
        }
        ordered_combined_covariates = sorted(
            combined_results_by_covariate.items(),
            key=lambda item: (covariate_order.get(item[0], 99), item[0]),
        )
        for covariate, combined_results in ordered_combined_covariates:
            covariate_label = {
                "age": "Age",
                "sex": "Sex",
                "bmi": "BMI",
                "cmv": "CMV",
                "smoking_status": "Smoking status",
                "log10_total_counts": "log10(total_counts)",
            }.get(covariate, covariate.title())
            covariate_anchor_suffix = re.sub(
                r"[^a-z0-9]+", "-", covariate_label.lower()
            ).strip("-")
            display(Markdown(
                f'<a id="combined-{covariate_anchor_suffix}"></a>\n#### {covariate_label}'
            ))
            fit_record = next(
                (model for model in combined_fit_records
                 if covariate in model.get("results_by_covariate", {})),
                {},
            )
            for contrast, contrast_results in combined_results.groupby("contrast", sort=True):
                volcano = contrast_results.copy()
                merged_gene_universe = set(volcano["gene"].astype(str))
                volcano["padj"] = pd.to_numeric(volcano["padj"], errors="coerce")
                volcano["log2FoldChange"] = pd.to_numeric(
                    volcano["log2FoldChange"], errors="coerce"
                )
                volcano = volcano.dropna(subset=["padj", "log2FoldChange"]).copy()
                if volcano.empty:
                    display(Markdown(
                        f"_{covariate} ({contrast}) has no finite adjusted p-values to plot._"
                    ))
                    continue
                volcano["minus_log10_padj"] = -np.log10(
                    volcano["padj"].clip(lower=np.finfo(float).tiny)
                )
                study_results = per_study_by_covariate.get(covariate, [])
                matching_study_results = pd.DataFrame()
                if study_results:
                    matching_study_results = pd.concat(study_results, ignore_index=True)
                    matching_study_results = matching_study_results.loc[
                        matching_study_results["contrast"].astype(str) == str(contrast)
                    ].copy()
                    matching_study_results["gene"] = matching_study_results[
                        "gene"
                    ].astype(str)
                    matching_study_results["padj"] = pd.to_numeric(
                        matching_study_results["padj"], errors="coerce"
                    )
                    study_hits = (
                        matching_study_results.loc[
                            matching_study_results["padj"] < de_alpha
                        ]
                        .groupby("gene")["study"]
                        .nunique()
                    )
                    volcano["studies_associated"] = (
                        volcano["gene"].astype(str).map(study_hits).fillna(0).astype(int)
                    )
                else:
                    volcano["studies_associated"] = 0
                studies_used = (
                    str(volcano["studies_included"].iloc[0]).split(";")
                    if "studies_included" in volcano else fit_record.get("studies", [])
                )
                sample_counts = (
                    json.loads(str(volcano["study_sample_counts"].iloc[0]))
                    if "study_sample_counts" in volcano
                    else fit_record.get("study_sample_counts", {})
                )
                study_summary = ", ".join(
                    f"{study} (n={sample_counts.get(study, 0)})" for study in studies_used if study
                )
                if not studies_used:
                    studies_used = sorted(sample_counts)
                if not study_summary and sample_counts:
                    study_summary = ", ".join(
                        f"{study} (n={count})" for study, count in sample_counts.items()
                    )
                nonintersection_age_markers = pd.DataFrame()
                if covariate == "age" and not matching_study_results.empty:
                    marker_results = matching_study_results.copy()
                    if studies_used:
                        marker_results = marker_results.loc[
                            marker_results["study"].astype(str).isin(studies_used)
                        ].copy()
                    marker_results["log2FoldChange"] = pd.to_numeric(
                        marker_results["log2FoldChange"], errors="coerce"
                    )
                    marker_results = marker_results.dropna(
                        subset=["gene", "padj", "log2FoldChange"]
                    )
                    # Compare study-only hits to the actual combined-model
                    # result universe. The intersection of per-study outputs
                    # can differ from the combined fit's availability mask.
                    nonintersection_age_markers = marker_results.loc[
                        (marker_results["padj"] < de_alpha)
                        & ~marker_results["gene"].astype(str).isin(merged_gene_universe)
                    ].copy()
                    if not nonintersection_age_markers.empty:
                        nonintersection_age_markers["minus_log10_padj"] = -np.log10(
                            nonintersection_age_markers["padj"].clip(
                                lower=np.finfo(float).tiny
                            )
                        )
                        nonintersection_age_markers = (
                            nonintersection_age_markers.sort_values(
                                ["study", "padj", "gene"], kind="stable"
                            )
                            .groupby("study", as_index=False, sort=True)
                            .head(1)
                        )
                max_associated_studies = max(
                    int(volcano["studies_associated"].max()), 1
                )
                # Leave room for the study-count legend while keeping these panels
                # comparable in plot-area width to the baseMean-colored panels.
                figure, axis = plt.subplots(figsize=(12.6, 6.5))
                sns.scatterplot(
                    data=volcano, x="log2FoldChange", y="minus_log10_padj",
                    hue="studies_associated", palette="viridis",
                    hue_norm=(0, max_associated_studies), legend="brief",
                    s=16, alpha=0.7, linewidth=0, ax=axis,
                )
                # Fix the plotting window from the merged results before adding
                # per-study markers, so study-specific extremes cannot compress
                # the merged volcano cloud.
                merged_xlim = axis.get_xlim()
                merged_ylim = axis.get_ylim()
                axis.set_xlim(merged_xlim)
                axis.set_ylim(merged_ylim)
                study_count_legend = axis.legend(
                    title="Per-study FDR-significant\nassociations", fontsize=8,
                    bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0,
                )
                if not nonintersection_age_markers.empty:
                    study_colors = dict(zip(
                        sorted(nonintersection_age_markers["study"].astype(str).unique()),
                        sns.color_palette("tab10", n_colors=nonintersection_age_markers["study"].nunique()),
                    ))
                    for _, marker in nonintersection_age_markers.iterrows():
                        study = str(marker["study"])
                        effect = float(marker["log2FoldChange"])
                        significance = float(marker["minus_log10_padj"])
                        in_x = merged_xlim[0] <= effect <= merged_xlim[1]
                        in_y = merged_ylim[0] <= significance <= merged_ylim[1]
                        if in_x and in_y:
                            axis.scatter(
                                effect, significance, marker="D", s=72,
                                facecolors="none", edgecolors=study_colors[study],
                                linewidths=1.8, label=study, zorder=4,
                            )
                            axis.annotate(
                                f"{marker['gene']} ({study})",
                                (effect, significance), xytext=(4, -10),
                                textcoords="offset points", fontsize=7,
                                color=study_colors[study],
                            )
                        else:
                            x_fraction = np.clip(
                                (effect - merged_xlim[0]) / (merged_xlim[1] - merged_xlim[0]),
                                0, 1,
                            )
                            y_fraction = np.clip(
                                (significance - merged_ylim[0]) / (merged_ylim[1] - merged_ylim[0]),
                                0, 1,
                            )
                            dx = np.sign(effect - np.clip(effect, *merged_xlim))
                            dy = np.sign(significance - np.clip(significance, *merged_ylim))
                            length = np.hypot(dx, dy) or 1
                            inset = 0.035
                            axis.annotate(
                                "", xy=(x_fraction, y_fraction),
                                xytext=(
                                    np.clip(x_fraction - inset * dx / length, 0.01, 0.99),
                                    np.clip(y_fraction - inset * dy / length, 0.01, 0.99),
                                ), xycoords="axes fraction",
                                arrowprops={
                                    "arrowstyle": "-|>", "color": study_colors[study],
                                    "lw": 1.5, "mutation_scale": 9,
                                }, zorder=5,
                            )
                    marker_studies = sorted(
                        nonintersection_age_markers["study"].astype(str).unique()
                    )
                    marker_handles = [
                        plt.Line2D(
                            [], [], marker="D", linestyle="", markersize=7,
                            markerfacecolor="none", markeredgecolor=study_colors[study],
                            label=study,
                        )
                        for study in marker_studies
                    ]
                    axis.legend(
                        handles=marker_handles,
                        title="Top FDR-significant\nstudy-only age gene",
                        fontsize=8, bbox_to_anchor=(1.02, 0.46),
                        loc="upper left", borderaxespad=0,
                    )
                    axis.add_artist(study_count_legend)
                axis.axhline(-np.log10(de_alpha), color="#555555", linestyle="--", linewidth=1)
                axis.axvline(0, color="#555555", linewidth=0.8)
                if covariate == "age":
                    effect_label = "log2 fold change per year"
                elif covariate == "bmi":
                    effect_label = "log2 fold change per BMI unit"
                elif covariate == "log10_total_counts":
                    effect_label = "log2 fold change per log10(total_counts) unit"
                else:
                    effect_label = f"log2 fold change ({contrast})"
                axis.set_xlabel(f"{covariate_label} effect ({effect_label})")
                axis.set_ylabel("−log10(adjusted p-value)")
                if len(studies_used) > 1:
                    title = f"Combined site-adjusted {covariate_label} association: {contrast}"
                else:
                    study_name = studies_used[0] if studies_used else "single study"
                    title = f"{covariate_label} association in {study_name}: {contrast}"
                axis.set_title(title)
                finite_label_rows = volcano.loc[
                    np.isfinite(volcano["padj"])
                    & np.isfinite(volcano["log2FoldChange"])
                ].copy()
                top_significance = finite_label_rows.sort_values(
                    ["padj", "gene"], kind="stable"
                ).head(20)
                top_effect = finite_label_rows.assign(
                    _absolute_log2_fold_change=(
                        finite_label_rows["log2FoldChange"].abs()
                    )
                ).sort_values(
                    ["_absolute_log2_fold_change", "padj", "gene"],
                    ascending=[False, True, True],
                    kind="stable",
                ).head(20)
                labels = pd.concat([top_significance, top_effect]).drop_duplicates(
                    "gene", keep="first"
                )
                for _, row in labels.iterrows():
                    axis.annotate(
                        str(row["gene"]),
                        (row["log2FoldChange"], row["minus_log10_padj"]),
                        xytext=(3, 3), textcoords="offset points", fontsize=7,
                    )
                figure.tight_layout(rect=(0, 0, 0.68 if not nonintersection_age_markers.empty else 0.78, 1))
                plt.show()
                if not nonintersection_age_markers.empty:
                    display(Markdown(
                        "Outlined diamonds use each study's own age-model estimate and adjusted "
                        "p-value. They select the smallest adjusted p-value among FDR-significant "
                        "genes outside the combined model's gene universe."
                    ))
                    display(nonintersection_age_markers[[
                        "study", "gene", "log2FoldChange", "padj",
                    ]].rename(columns={
                        "log2FoldChange": "per-study log2 fold change per year",
                        "padj": "per-study adjusted p-value",
                    }).sort_values("study", kind="stable"))
                if covariate == "log10_total_counts" and "baseMean" in volcano.columns:
                    base_mean = pd.to_numeric(volcano["baseMean"], errors="coerce")
                    finite = np.isfinite(base_mean) & (base_mean > 0)
                    if finite.any():
                        log_base_mean = np.log10(base_mean.loc[finite])
                        base_mean_volcano = volcano.loc[finite]
                        low = float(log_base_mean.min())
                        high = float(log_base_mean.max())
                        norm = Normalize(vmin=low, vmax=high if high > low else low + 0.01)
                        figure, axis = plt.subplots(figsize=(10.5, 6.5))
                        points = axis.scatter(
                            base_mean_volcano["log2FoldChange"],
                            base_mean_volcano["minus_log10_padj"],
                            c=log_base_mean,
                            cmap="viridis",
                            norm=norm,
                            s=16,
                            alpha=0.7,
                            linewidths=0,
                        )
                        colorbar = figure.colorbar(points, ax=axis, pad=0.02)
                        colorbar.set_label("log10(combined-fit baseMean)")
                        axis.axhline(
                            -np.log10(de_alpha), color="#555555", linestyle="--", linewidth=1
                        )
                        axis.axvline(0, color="#555555", linewidth=0.8)
                        axis.set_xlabel(f"{covariate_label} effect ({effect_label})")
                        axis.set_ylabel("−log10(adjusted p-value)")
                        axis.set_title(f"{title} (colored by log10(baseMean))")
                        for _, row in labels.iterrows():
                            axis.annotate(
                                str(row["gene"]),
                                (row["log2FoldChange"], row["minus_log10_padj"]),
                                xytext=(3, 3), textcoords="offset points", fontsize=7,
                            )
                        figure.tight_layout()
                        plt.show()
                fit_design = fit_record.get(
                    "design",
                    str(volcano["design"].iloc[0]) if "design" in volcano else "—",
                )
                study_sites = (
                    str(volcano["study_sites_included"].iloc[0]).split(";")
                    if "study_sites_included" in volcano
                    else fit_record.get("study_sites", [])
                )
                details = pd.DataFrame(
                    {
                        "Value": [
                            contrast,
                            fit_design,
                            f"{de_alpha:g}",
                            int((volcano["padj"] < de_alpha).sum()),
                            study_summary or "unavailable in this legacy result set",
                            ", ".join(site for site in study_sites if site) or "—",
                            sum(sample_counts.values()) if sample_counts else "unavailable",
                            fit_record.get(
                                "n_samples_excluded_missing_design_covariates", "—"
                            ),
                        ]
                    },
                    index=[
                        "Contrast",
                        "Model design",
                        "FDR threshold",
                        "FDR-significant genes",
                        "Studies included (samples per study)",
                        "Study sites included",
                        "Total samples used",
                        "Samples excluded for missing design values",
                    ],
                )
                display(details)
                if covariate == "age":
                    display_merged_age_go_enrichment(contrast_results, contrast)
        display(Markdown(
            "_DE interpretation: the recurrence panels emphasize agreement among available per-study "
            "models; each combined volcano summarizes eligible samples with study-site, age, sex, "
            "and log-total-count adjustment, plus the target optional covariate when applicable. "
            "In the primary volcano, point color gives the number of available per-study fits "
            "for the same covariate and contrast with FDR-significant association for that gene "
            "(using the displayed FDR threshold); the plotted effect and adjusted p-value come "
            "from the combined fit. The log10(total_counts) volcano also has a companion colored "
            "by combined-fit log10(baseMean). "
            "For age, outlined diamonds mark the smallest per-study adjusted p-value among "
            "FDR-significant genes outside the combined model's gene universe; their position instead uses that study's "
            "own effect and adjusted p-value, so they are a visibility check rather than combined-fit points. "
            "The site term is omitted when only one site remains estimable. Combined fits do not "
            "by themselves establish replication across studies._"
        ))
        per_study_covariate_summary = []
        for result in per_study_results:
            for (covariate, contrast), subset in result.groupby(
                ["covariate", "contrast"], observed=True
            ):
                sample_counts = json.loads(str(subset["study_sample_counts"].iloc[0]))
                study = str(subset["study"].iloc[0])
                fit_record = next(
                    (model for model in de_run_metadata.get("models", [])
                     if model.get("model") == "per_study"
                     and model.get("study") == study
                     and covariate in model.get("results_by_covariate", {})),
                    {},
                )
                per_study_covariate_summary.append({
                    "study": study,
                    "covariate": covariate,
                    "contrast": contrast,
                    "design": fit_record.get("design", "—"),
                    "complete-case samples used": (
                        sum(sample_counts.values()) if sample_counts else "—"
                    ),
                    "excluded for missing design values": fit_record.get(
                        "n_samples_excluded_missing_design_covariates", "—"
                    ),
                    "genes tested": subset["gene"].nunique(),
                    "FDR-significant genes": int(
                        (pd.to_numeric(subset["padj"], errors="coerce") < de_alpha).sum()
                    ),
                })
        if per_study_covariate_summary:
            display(Markdown(
                "<a id=\"per-study-covariate-models\"></a>\n### Per-Study Covariate Models"
            ))
            display(Markdown(
                "Each study is fitted once with its maximal available design. Every coefficient for "
                "that study comes from the same complete-case sample set and formula."
            ))
            per_study_summary = pd.DataFrame(per_study_covariate_summary)
            per_study_summary["_covariate_order"] = per_study_summary["covariate"].map({
                "age": 0, "sex": 1, "bmi": 2, "cmv": 3,
                "smoking_status": 4, "log10_total_counts": 5,
            }).fillna(99)
            per_study_summary = (
                per_study_summary
                .sort_values(
                    ["study", "_covariate_order", "contrast"], kind="stable",
                )
                .drop(columns="_covariate_order")
            )
            display_collapsible_table(
                per_study_summary,
                "Show per-study fits, sorted by study then coefficient",
            )

        display(Markdown(
            "Per-study fits use one maximal available design and also adjust for site when "
            "multiple sites contribute. Missing values exclude a sample only from models whose "
            "design includes that variable. Categorical effects use female as the sex reference, "
            "no CMV (negative where that is the source label), and no smoking status as references."
        ))
        if summary_rows:
            display_collapsible_table(
                pd.DataFrame(summary_rows), "Show per-study age gene-count summary",
            )

        if has_age_per_study:
            display(Markdown(
                "### Per-Study Age-Association Recurrence in the Common Tested-Gene Set"
            ))
            display(Markdown(
                "These plots and the list below include genes tested by every available "
                "per-study age fit and FDR-significant in at least one study. They summarize "
                "overlap among per-study fits; they do not show the combined-model volcano "
                "gene universe."
            ))
            figure, axes = plt.subplots(1, 2, figsize=(14.5, 8.7))
            if recurrence.empty:
                empty_message = (
                    "No genes were tested in every per-study model"
                    if not common_per_study_genes
                    else "No common-universe genes pass the FDR threshold"
                )
                axes[0].text(0.5, 0.5, empty_message,
                             transform=axes[0].transAxes, ha="center", va="center")
                axes[1].text(0.5, 0.5, empty_message,
                             transform=axes[1].transAxes, ha="center", va="center")
            else:
                top_genes = recurrence.head(50).sort_values("studies_associated")
                sns.barplot(
                    data=top_genes, x="studies_associated", y="gene", color="#4c72b0",
                    ax=axes[0],
                )
                axes[0].set_xlabel("Studies with FDR-significant age association")
                axes[0].set_ylabel("")
                axes[0].set_title("Genes recurring across studies (top 50)")
                recurrence_counts = recurrence["studies_associated"].value_counts().sort_index()
                sns.barplot(
                    x=recurrence_counts.index.astype(str), y=recurrence_counts.values,
                    color="#55a868", ax=axes[1],
                )
                axes[1].set_xlabel("Number of studies")
                axes[1].set_ylabel("Unique associated genes")
                axes[1].set_title(
                    f"Common-universe union: {recurrence['gene'].nunique():,} genes"
                )
            figure.tight_layout()
            plt.show()
            if not recurrence.empty:
                display(Markdown(
                    "The effect shown for each gene is the **signed** age log2 fold change "
                    "with the largest absolute magnitude among its per-study estimates."
                ))
                display_collapsible_table(
                    recurrence.head(50),
                    "Show top 50 recurrent per-study age genes in the common tested-gene set",
                    max_height=480,
                )

        if not intersection_diagnostics.empty:
            display(Markdown(
                "<a id=\"per-study-gene-intersection\"></a>\n"
                "### Per-Study Gene-Intersection Diagnostic"
            ))
            display(Markdown(
                "This table shows how much restricting a cross-study comparison to genes tested "
                "by every available study would remove from each study's significant-gene list. "
                "It is a gene-universe sensitivity check, not a replication result. Each row uses "
                "the displayed single-study design and complete-case sample count; age, sex, BMI, "
                "CMV, and smoking-status coefficients from the same study therefore come from that one maximal model."
            ))
            display_collapsible_table(intersection_diagnostics.assign(
                **{"outside intersection (%)": intersection_diagnostics[
                    "outside intersection (%)"
                ].map(lambda value: f"{value:.1f}%")}
            ), "Show per-study gene-intersection sensitivity details")

        if de_run_metadata:
            adult_age_fit = next(
                (model for model in de_run_metadata.get("models", [])
                 if model.get("model") == "merged"),
                {},
            )
            developmental_age_fit = next(
                (model for model in de_run_metadata.get("models", [])
                 if model.get("purpose") == "developmental_adult_effect_comparison"),
                {},
            )
            if (
                adult_age_fit.get("status") == "complete"
                and developmental_age_fit.get("status") == "complete"
            ):
                adult_age_path = Path(de_dir) / adult_age_fit.get(
                    "results_by_covariate", {}
                ).get("age", "")
                developmental_age_path = Path(de_dir) / developmental_age_fit.get(
                    "results_by_covariate", {}
                ).get("age", "")
                required_age_columns = {"gene", "log2FoldChange", "padj"}
                if adult_age_path.is_file() and developmental_age_path.is_file():
                    adult_age = pd.read_csv(adult_age_path)
                    developmental_age = pd.read_csv(developmental_age_path)
                    if (
                        required_age_columns.issubset(adult_age.columns)
                        and required_age_columns.issubset(developmental_age.columns)
                    ):
                        developmental_alpha = float(
                            pipeline["differential_expression"]["alpha"]
                        )
                        developmental_volcano = developmental_age.copy()
                        developmental_volcano["log2FoldChange"] = pd.to_numeric(
                            developmental_volcano["log2FoldChange"], errors="coerce"
                        )
                        developmental_volcano["padj"] = pd.to_numeric(
                            developmental_volcano["padj"], errors="coerce"
                        )
                        developmental_volcano = developmental_volcano.replace(
                            [np.inf, -np.inf], np.nan
                        ).dropna(subset=["log2FoldChange", "padj"])
                        developmental_volcano = developmental_volcano.loc[
                            developmental_volcano["padj"].between(0, 1)
                        ].copy()
                        if not developmental_volcano.empty:
                            developmental_volcano["minus_log10_padj"] = -np.log10(
                                developmental_volcano["padj"].clip(
                                    lower=np.finfo(float).tiny, upper=1.0
                                )
                            )
                            developmental_significant = (
                                developmental_volcano["padj"] < developmental_alpha
                            )
                            figure, axis = plt.subplots(figsize=(9.5, 5.8))
                            axis.scatter(
                                developmental_volcano.loc[~developmental_significant, "log2FoldChange"],
                                developmental_volcano.loc[~developmental_significant, "minus_log10_padj"],
                                color="#aab4bf", s=15, alpha=0.55, linewidths=0,
                                label="Not FDR significant",
                            )
                            axis.scatter(
                                developmental_volcano.loc[developmental_significant, "log2FoldChange"],
                                developmental_volcano.loc[developmental_significant, "minus_log10_padj"],
                                color="#c44e52", s=20, alpha=0.78, linewidths=0,
                                label=f"FDR < {developmental_alpha:g}",
                            )
                            axis.axvline(0, color="#666666", linewidth=0.8)
                            axis.axhline(
                                -np.log10(developmental_alpha), color="#666666",
                                linewidth=0.9, linestyle="--",
                            )
                            axis.set(
                                xlabel="Juvenile age effect ([0,20) years)\n"
                                "raw log2 fold change per year",
                                ylabel="−log10(adjusted p-value)",
                                title=(
                                    f"Juvenile age-associated genes "
                                    f"({int(developmental_significant.sum()):,} FDR significant)"
                                ),
                            )
                            axis.legend(frameon=True, fontsize=8)
                            figure.tight_layout()
                            plt.show()

                            juvenile_gene_rows = (
                                developmental_volcano.loc[developmental_significant]
                                if developmental_significant.any()
                                else developmental_volcano
                            )
                            top_juvenile_genes = juvenile_gene_rows.sort_values(
                                ["padj", "gene"], kind="stable"
                            ).head(25)
                            display_columns = [
                                column for column in (
                                    "gene", "log2FoldChange", "padj", "baseMean"
                                ) if column in top_juvenile_genes
                            ]
                            top_juvenile_genes = top_juvenile_genes[display_columns].rename(
                                columns={
                                    "gene": "Gene",
                                    "log2FoldChange": "log2 fold change per year",
                                    "padj": "Adjusted p-value",
                                    "baseMean": "Mean expression",
                                }
                            )
                            n_developmental_significant = int(developmental_significant.sum())
                            summary = (
                                f"Top juvenile age-DE genes (FDR < {developmental_alpha:g})"
                                if n_developmental_significant
                                else "No juvenile genes pass the FDR threshold; show top 25 tested genes"
                            )
                            display_collapsible_table(
                                top_juvenile_genes,
                                summary,
                            )
                        else:
                            display(Markdown(
                                "_No juvenile genes have finite age effects and adjusted "
                                "p-values for the volcano plot._"
                            ))
                        paired_age_effects = developmental_age[
                            ["gene", "log2FoldChange", "padj"]
                        ].rename(columns={
                            "log2FoldChange": "developmental_effect",
                            "padj": "developmental_padj",
                        }).merge(
                            adult_age[["gene", "log2FoldChange", "padj"]].rename(
                                columns={
                                    "log2FoldChange": "adult_effect",
                                    "padj": "adult_padj",
                                }
                            ),
                            on="gene", how="inner", validate="one_to_one",
                        )
                        for column in [
                            "developmental_effect", "adult_effect",
                            "developmental_padj", "adult_padj",
                        ]:
                            paired_age_effects[column] = pd.to_numeric(
                                paired_age_effects[column], errors="coerce"
                            )
                        plot_age_effects = paired_age_effects.replace(
                            [np.inf, -np.inf], np.nan
                        ).dropna(subset=[
                            "developmental_effect", "adult_effect",
                            "developmental_padj", "adult_padj",
                        ])
                        display(Markdown("### Developmental versus adult age effects"))
                        support_path = de_dir / developmental_age_fit.get(
                            "age_support_samples", "developmental/age_support_samples.csv",
                        )
                        if support_path.is_file():
                            developmental_support_samples = pd.read_csv(support_path)
                            if not developmental_support_samples.empty:
                                support_column = (
                                    "study_site"
                                    if "study_site" in developmental_support_samples
                                    else "study"
                                )
                                support_samples = developmental_support_samples.dropna(
                                    subset=[support_column, "age"],
                                ).copy()
                                support_order = sorted(
                                    support_samples[support_column].astype(str).unique()
                                )
                                figure, axis = plt.subplots(
                                    figsize=(max(10, 0.42 * len(support_order) + 5), 4.8)
                                )
                                sns.stripplot(
                                    data=support_samples,
                                    x=support_column,
                                    y="age",
                                    hue="sex" if "sex" in support_samples else None,
                                    order=support_order,
                                    hue_order=(
                                        sorted(support_samples["sex"].dropna().astype(str).unique())
                                        if "sex" in support_samples else None
                                    ),
                                    dodge="sex" in support_samples,
                                    jitter=0.18,
                                    alpha=0.8,
                                    ax=axis,
                                )
                                axis.set_xlabel(
                                    "Study site" if support_column == "study_site" else "Study"
                                )
                                axis.set_ylabel("Age (years)")
                                axis.set_title("Developmental age support by study site ([0,20) years)")
                                axis.tick_params(axis="x", rotation=75, labelsize=7)
                                if axis.legend_ is not None:
                                    axis.legend(title="Sex", fontsize=8, title_fontsize=8)
                                figure.tight_layout()
                                display_collapsible_figure(
                                    figure,
                                    "Developmental age support by study site",
                                    open_by_default=True,
                                )
                            else:
                                display(Markdown("_No samples were retained for the [0,20) age-support plot._"))
                        else:
                            display(Markdown(
                                "_Developmental sample-level support is unavailable; rerun the DE fit to publish it._"
                            ))
                        display(Markdown(
                            "Each point in the scatterplots is a gene, comparing its raw age "
                            "coefficient from the [0,20)-year fit with its coefficient from the "
                            "adult fit (ages ≥20). Points require a non-missing adjusted p-value "
                            "in both fits, which is PyDESeq2’s independent expression filter; no "
                            "significance cutoff is used. The developmental range is [0,20) and "
                            "the adult range starts at 20, so the age ranges do not overlap. "
                            "Both panels use the same color scale."
                        ))
                        if plot_age_effects.empty:
                            display(Markdown(
                                "_No genes pass the expression filter in both fits with finite "
                                "age coefficients._"
                            ))
                        else:
                            developmental_color = -np.log10(
                                plot_age_effects["developmental_padj"].clip(
                                    lower=np.finfo(float).tiny, upper=1.0
                                )
                            )
                            adult_color = -np.log10(
                                plot_age_effects["adult_padj"].clip(
                                    lower=np.finfo(float).tiny, upper=1.0
                                )
                            )
                            color_max = max(
                                1.0,
                                float(pd.concat(
                                    [developmental_color, adult_color], ignore_index=True
                                ).quantile(0.99)),
                            )
                            effect_limit = max(
                                float(np.nanmax(np.abs(plot_age_effects[[
                                    "developmental_effect", "adult_effect",
                                ]].to_numpy(dtype=float)))) * 1.05,
                                0.01,
                            )
                            figure = plt.figure(figsize=(12.2, 5.8), layout="constrained")
                            grid = figure.add_gridspec(
                                2, 2, height_ratios=[1, 0.075],
                            )
                            axes = [figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 1])]
                            colorbar_axis = figure.add_subplot(grid[1, :])
                            scatter = None
                            for axis, color_values, title in [
                                (axes[0], developmental_color, "Colored by developmental −log10(padj)"),
                                (axes[1], adult_color, "Colored by adult −log10(padj)"),
                            ]:
                                scatter = axis.scatter(
                                    plot_age_effects["developmental_effect"],
                                    plot_age_effects["adult_effect"],
                                    c=color_values.clip(upper=color_max),
                                    s=14, alpha=0.65, linewidths=0,
                                    cmap="viridis", vmin=0, vmax=color_max,
                                )
                                axis.axhline(0, color="#666666", linewidth=0.8)
                                axis.axvline(0, color="#666666", linewidth=0.8)
                                axis.plot(
                                    [-effect_limit, effect_limit],
                                    [-effect_limit, effect_limit],
                                    color="#c44e52", linewidth=0.9, linestyle="--",
                                )
                                axis.set(
                                    xlim=(-effect_limit, effect_limit),
                                    ylim=(-effect_limit, effect_limit),
                                    xlabel="Developmental age effect ([0,20) years)\n"
                                    "raw log2 fold change per year",
                                    ylabel="Adult age effect (ages ≥20)\n"
                                    "raw log2 fold change per year",
                                    title=title,
                                )
                                axis.set_aspect("equal", adjustable="box")
                            colorbar = figure.colorbar(
                                scatter, cax=colorbar_axis, orientation="horizontal",
                            )
                            colorbar.set_label(
                                f"−log10(adjusted p-value; values above {color_max:.1f} capped)"
                            )
                            figure.suptitle(
                                f"{cell_type_name}: developmental versus adult age effects "
                                f"(n={len(plot_age_effects):,} genes)"
                            )
                            plt.show()
        if age_trajectory_fit or age_trajectory_results is not None:
            display(Markdown(
                "<a id=\"age-bin-trajectories\"></a>\n## Age-Bin Trajectories"
            ))
            if de_dir:
                display(Markdown(
                    "See the [cross-cell-type trajectory report](../../../combined/trajectory/"
                    "report.html) for "
                    "shared-pattern clusters and gene recurrence across cell types."
                ))
            width = trajectory_settings.bin_width_years
            minimum_bin_samples = trajectory_settings.minimum_samples_per_bin
            minimum_bins = trajectory_settings.minimum_bins
            strict_age_cutoff = trajectory_settings.strict_age_cutoff_exclusive
            age_cutoff_excluded = age_trajectory_info.get(
                "samples_excluded_by_age_cutoff", 0
            )
            manually_excluded_bins = age_trajectory_info.get(
                "manually_excluded_age_bins", []
            )
            raw_bin_counts = age_trajectory_info.get(
                "age_bin_counts_before_filtering", {}
            )
            manual_bin_counts_text = ", ".join(
                f"{age_bin}: {raw_bin_counts.get(age_bin, 0)} eligible pseudobulks"
                for age_bin in manually_excluded_bins
            ) or "none"
            reference_bin = age_trajectory_info.get(
                "reference_bin",
                f"{trajectory_settings.reference_bin_start_age}-"
                f"{trajectory_settings.reference_bin_start_age + width}",
            )
            retained_bins = age_trajectory_info.get("retained_bins", [])
            display(Markdown(
                f"Trajectory sample filtering: manually excluded age bin(s) "
                f"{manual_bin_counts_text}. The strict age cutoff is < {strict_age_cutoff} years; "
                f"{age_cutoff_excluded} eligible pseudobulks aged at or above it were excluded "
                "(this total includes the manually excluded 90–100 bin). Each retained "
                f"{width}-year bin must contain at least {minimum_bin_samples} eligible pseudobulks."
            ))
            if age_trajectory_fit.get("status") != "complete" or age_trajectory_results is None:
                display(Markdown(
                    f"_Age-bin trajectory fit skipped: "
                    f"{age_trajectory_fit.get('reason', 'no trajectory result was produced')}._"
                ))
            elif len(retained_bins) < trajectory_settings.minimum_bins:
                display(Markdown(
                    f"_Age-bin trajectory results omitted: only {len(retained_bins)} bins "
                    f"meet the sample-count threshold; at least "
                    f"{trajectory_settings.minimum_bins} "
                    "are required._"
                ))
            else:
                dropped_bins = age_trajectory_info.get("dropped_bins", [])
                top_gene_count = trajectory_settings.report_top_n_genes
                trajectory_de_fdr = trajectory_settings.de_fdr_threshold
                display(Markdown(
                    f"The combined model uses {width}-year age bins, adjusts for study and sex, "
                    f"and uses **{reference_bin}** as the zero point. Bins with fewer than "
                    f"{minimum_bin_samples} eligible pseudobulks are excluded; at least "
                    f"{minimum_bins} supported bins are required. Retained bins: "
                    f"{', '.join(retained_bins)}. "
                    f"Dropped bins: {', '.join(dropped_bins) if dropped_bins else 'none'}. "
                    "The omnibus joint Wald test asks whether any retained age bin differs from "
                    "the reference, with FDR controlled across genes. Each trajectory is scaled "
                    "by its population SD across all retained bins, including the zero reference; "
                    "a trajectory with incomplete estimates or zero SD has no standardized shape."
                ))
                omnibus = age_trajectory_results.drop_duplicates("gene").copy()
                omnibus["omnibus_padj"] = pd.to_numeric(
                    omnibus["omnibus_padj"], errors="coerce"
                )
                selected_genes = omnibus.loc[
                    omnibus["omnibus_padj"] < trajectory_de_fdr
                ].sort_values("omnibus_padj").head(top_gene_count)
                if selected_genes.empty:
                    display(Markdown(
                        f"_No genes pass the omnibus age-bin FDR threshold "
                        f"(< {trajectory_de_fdr:g}) for this cell type._"
                    ))
                else:
                    selected_names = selected_genes["gene"].astype(str).tolist()
                    selected_shape = age_trajectory_results.loc[
                        age_trajectory_results["gene"].astype(str).isin(selected_names)
                    ].copy()
                    shape_matrix = selected_shape.pivot(
                        index="gene", columns="age_bin", values="trajectory_scaled"
                    ).reindex(index=selected_names, columns=retained_bins)
                    figure, axis = plt.subplots(
                        figsize=(max(7.0, 1.1 * len(retained_bins)),
                                 max(3.5, 0.28 * len(selected_names) + 1.5))
                    )
                    sns.heatmap(
                        shape_matrix, cmap="vlag", center=0, ax=axis,
                        cbar_kws={"label": "Trajectory shape (SD units)"},
                    )
                    axis.set_xlabel("Age bin (years)")
                    axis.set_ylabel("Gene")
                    axis.set_title("Top omnibus-significant standardized trajectories")
                    figure.tight_layout()
                    display_collapsible_figure(figure, "Show omnibus-significant age trajectories")

                    selected_table = selected_shape.pivot(
                        index="gene", columns="age_bin", values="trajectory_scaled"
                    ).reindex(index=selected_names, columns=retained_bins)
                    selected_table.insert(
                        0, "Omnibus FDR", selected_genes.set_index("gene").loc[
                            selected_table.index, "omnibus_padj"
                        ].to_numpy(),
                    )
                    selected_table.insert(
                        1, "Trajectory SD (log2 fold change)", selected_genes.set_index(
                            "gene"
                        ).loc[selected_table.index, "trajectory_sd"].to_numpy(),
                    )
                    selected_table.columns = [
                        str(column) if column in {"Omnibus FDR", "Trajectory SD (log2 fold change)"}
                        else f"Shape {column}"
                        for column in selected_table.columns
                    ]
                    display_collapsible_table(
                        selected_table.reset_index(),
                        "Show top omnibus-significant age-bin trajectories",
                    )

                cluster_fdr = trajectory_settings.cluster_fdr_threshold
                cluster_significant = omnibus.loc[
                    omnibus["omnibus_padj"] <= cluster_fdr
                ]
                if cluster_significant.empty:
                    display(Markdown(
                        f"_No genes pass the trajectory clustering FDR threshold "
                        f"({cluster_fdr:g})._"
                    ))
                else:
                    profile_data = age_trajectory_results.loc[
                        age_trajectory_results["gene"].astype(str).isin(
                            cluster_significant["gene"].astype(str)
                        )
                    ].copy()
                    profiles = profile_data.pivot(
                        index="gene", columns="age_bin", values="trajectory_scaled"
                    ).reindex(columns=retained_bins).dropna()
                    cluster_labels, cluster_means = cluster_trajectory_profiles(
                        profiles,
                        trajectory_settings,
                    )
                    merge_diagnostics = trajectory_merge_diagnostics(
                        profiles, trajectory_settings, maximum_merges=20
                    )
                    cluster_rows = profiles.join(cluster_labels).reset_index()
                    cluster_rows.insert(1, "cell_type", cell_type_name)
                    cluster_rows["omnibus_padj"] = cluster_rows["gene"].map(
                        omnibus.set_index("gene")["omnibus_padj"]
                    )
                    cluster_rows.to_csv(output_dir / "age_trajectory_clusters.csv", index=False)
                    cluster_means = cluster_means.reset_index()
                    cluster_means.insert(0, "cell_type", cell_type_name)
                    cluster_means.to_csv(
                        output_dir / "age_trajectory_cluster_means.csv", index=False
                    )
                    figure, axes = plt.subplots(1, 3, figsize=(18, 4.8))
                    palette = sns.color_palette(
                        "husl", n_colors=max(len(cluster_means), 1)
                    )
                    cluster_colors = {
                        str(int(row["trajectory_cluster"])): palette[index]
                        for index, row in cluster_means.iterrows()
                    }
                    for index, row in cluster_means.iterrows():
                        cluster = int(row["trajectory_cluster"])
                        label = f"Cluster {cluster} (n={int(row['n_trajectories'])})"
                        axes[0].plot(
                            retained_bins,
                            [row[age_bin] for age_bin in retained_bins],
                            marker="o", linewidth=2, color=cluster_colors[str(cluster)],
                            label=label,
                        )
                    axes[0].axhline(0, color="#555555", linestyle="--", linewidth=0.8)
                    axes[0].set_xlabel("Age bin")
                    axes[0].set_ylabel("Mean standardized trajectory (SD units)")
                    axes[0].set_title("Mean trajectory by hierarchical cluster")
                    axes[0].tick_params(axis="x", rotation=35)
                    axes[0].legend(fontsize=8)
                    embedded = cluster_labels.dropna(subset=["umap_1", "umap_2"])
                    if embedded.empty:
                        axes[1].text(
                            0.5, 0.5,
                            f"At least {trajectory_settings.minimum_umap_trajectories} "
                            "trajectories are needed for UMAP",
                            ha="center", va="center", transform=axes[1].transAxes,
                        )
                    else:
                        plot_points = embedded.reset_index()
                        plot_points["trajectory_cluster"] = plot_points[
                            "trajectory_cluster"
                        ].astype(str)
                        sns.scatterplot(
                            data=plot_points, x="umap_1", y="umap_2",
                            hue="trajectory_cluster", palette=cluster_colors, s=24,
                            alpha=0.8, linewidth=0, ax=axes[1], legend="brief",
                        )
                        axes[1].set_title(
                            f"Trajectory UMAP ({trajectory_settings.distance_metric} "
                            "bin distance)"
                        )
                        axes[1].set_xlabel("UMAP 1")
                        axes[1].set_ylabel("UMAP 2")
                    plot_trajectory_merge_diagnostic(
                        axes[2], merge_diagnostics,
                        max_clusters=trajectory_settings.max_clusters,
                    )
                    figure.tight_layout()
                    display_collapsible_figure(
                        figure,
                        f"Show trajectory clusters, UMAP, and merge diagnostic ({len(profiles):,} genes, "
                        f"FDR ≤ {cluster_fdr:g})",
                        open_by_default=True,
                    )
                    display(Markdown(
                        "The merge-cost panel shows up to the final 20 hierarchical merges. "
                        "A sharp increase at **K → K−1** supports retaining K trajectory groups; "
                        "the red marker shows the configured cluster cap when that transition is available."
                    ))
                    display(Markdown(
                        "Download the [gene cluster assignments and standardized profiles]"
                        "(age_trajectory_clusters.csv) or [cluster mean trajectories]"
                        "(age_trajectory_cluster_means.csv)."
                    ))

        display(Markdown(
            "<a id=\"sample-residual-clustering\"></a>\n"
                "## Per-Sample Residual Clustering"
        ))
        display(Markdown(
            "These sample profiles come from the merged age-bin model. It replaces continuous "
            "age with configured decade bins and includes sex, plus study-site and "
            "log10(total counts) terms when estimable. For each sample and gene, the saved "
            "Pearson residual is the observed count minus the fitted mean, divided by the "
            "fitted standard deviation \\(\\sqrt{\\mu + \\alpha\\mu^2}\\); it is a scaled "
            "deviation from the model fit. PCA centers each gene's residuals across samples, "
            "then retains up to the first 100 principal components, capped by the available "
            "samples and genes. Those PC scores are the input to the two-dimensional UMAP "
            "shown below. Colors overlay sample covariates and residual-cluster labels on "
            "the same UMAP coordinates."
        ))
        display(Markdown("### Residual UMAP and PCA plots"))
        if de_dir and age_trajectory_fit:
            (
                residual_clusters, residual_pc_scores, residual_pcs, residual_matrix,
            ) = load_residual_sample_clusters(de_dir, age_trajectory_fit, trajectory_settings)
        else:
            residual_clusters, residual_pc_scores, residual_pcs, residual_matrix = (
                None, None, 0, None
            )
        if residual_clusters is None:
            residual_reason = age_trajectory_fit.get(
                "reason", "the age-bin model did not save Pearson residuals"
            )
            display(Markdown(f"_Residual clustering unavailable: {residual_reason}._"))
        else:
            residual_clusters.to_csv(output_dir / "residual_sample_clusters.csv", index=False)
            residual_pc_scores.to_csv(
                output_dir / "residual_sample_pc_scores.csv", index=False
            )
            display(Markdown(
                f"The embedding includes **{len(residual_clusters)} samples** and "
                f"**{residual_matrix.shape[1]:,} genes**, with {residual_pcs} "
                "principal components retained. The primary clusters use Leiden on a "
                f"{trajectory_settings.residual_leiden_neighbors}-neighbor graph built from "
                f"up to {trajectory_settings.residual_pc_hdbscan_components} PCs. The "
                f"selected resolution targets at least "
                f"{trajectory_settings.residual_leiden_min_clusters} clusters."
            ))
            plot_residual_sample_umaps(
                residual_clusters,
                minimum_umap_samples=trajectory_settings.minimum_umap_trajectories,
            )
            plot_residual_pca_and_pc_umaps(residual_pc_scores, residual_clusters)
            display(Markdown("### Alternative clusterings"))
            display(Markdown(
                "The primary Leiden resolution is selected from the configured candidates "
                f"{trajectory_settings.residual_leiden_resolutions}: choose the lowest that "
                f"reaches {trajectory_settings.residual_leiden_min_clusters} clusters, or "
                "the candidate with the most clusters if none reach the target. The table "
                "shows candidate cluster counts; plots compare K-means and one configured "
                "HDBSCAN result on UMAP and PCs. Marker contrasts cover at most the "
                f"{trajectory_settings.residual_marker_max_clusters} largest clusters."
            ))
            plot_residual_alternative_clusterings(residual_clusters, trajectory_settings)
            plot_residual_cluster_markers(
                residual_matrix, residual_clusters, trajectory_settings
            )
            display(Markdown(
                "Download the [residual cluster assignments](residual_sample_clusters.csv) "
                "or [sample PCA scores](residual_sample_pc_scores.csv)."
            ))
            display(Markdown("### Residual WGCNA modules"))
            wgcna_settings = ResidualWGCNASettings.from_mapping(
                pipeline["cell_type_analysis"]["residual_wgcna"]
            )
            residual_wgcna = run_residual_wgcna(
                residual_matrix,
                residual_clusters,
                settings=wgcna_settings,
                cmv_results=load_merged_cmv_results(de_dir),
                smoking_results=load_merged_smoking_results(de_dir),
            )
            publish_residual_wgcna(residual_wgcna)
            display(Markdown(
                "PyWGCNA fits a signed correlation network to the age-bin-model Pearson "
                "residuals, after retaining the configured most variable genes. The residual "
                "model already accounts for age bins, sex, study site, and pseudobulk depth. "
                "CMV and smoking are omitted from that model, so these checks focus on two "
                "plausible sources of remaining biological structure in its residual modules. "
                "When its merged model is available, the module gene-overlap check uses "
                "FDR-significant genes from that model and the selected network genes as background; the "
                "eigengene check is a point-biserial association with observed yes/no status. "
                "The selected exposures are examined because they were not adjusted out of the "
                "residual model, not because residual effects are assumed to be present. These "
                "are exploratory, unadjusted checks."
            ))
            display(Markdown(
                f"The network used **{residual_wgcna.metadata['n_samples']:,} samples** and "
                f"**{residual_wgcna.metadata['n_genes_selected']:,} genes**. "
                "Download its [module assignments](residual_wgcna/module_assignments.csv), "
                "[eigengenes](residual_wgcna/eigengenes.csv), [module summary]"
                "(residual_wgcna/module_summary.csv), and [soft-threshold diagnostics]"
                "(residual_wgcna/soft_threshold.csv), [sample connectivity]"
                "(residual_wgcna/sample_connectivity.csv), [covariate associations]"
                "(residual_wgcna/module_covariate_associations.csv), [hub genes]"
                "(residual_wgcna/hub_genes.csv), and [leave-one-study-out stability]"
                "(residual_wgcna/leave_one_study_out_stability.csv)."
            ))
            plot_residual_wgcna(
                residual_wgcna,
                residual_matrix,
                residual_clusters,
                wgcna_settings,
                load_merged_cmv_results(de_dir),
                load_merged_smoking_results(de_dir),
            )
