"""Prepare per-cell-type DE diagnostics independently of notebook rendering."""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import shutil
import tempfile
from importlib.metadata import version
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, pointbiserialr, spearmanr
from statsmodels.stats.multitest import multipletests

from .abundance import fit_parent_fraction_models
from .config import (
    AgeTrajectorySettings,
    ResidualWGCNASettings,
    config_digest,
    load_analysis_configuration,
    read_json,
)
from .fraction_age import fit_fraction_age_bins
from .go_enrichment import (
    DEFAULT_GO_BP_CONFIG,
    GO_ENRICHMENT_DATA_DIR_ENV,
    GeneSetEnrichmentRun,
    run_preranked_gsea,
    save_go_specificity_outputs,
    signed_de_gene_ranking,
)
from .logging_utils import configure_logging
from .reporting import sample_cell_type_fractions
from .residual_wgcna import ResidualWGCNAResult, run_residual_wgcna, trait_associated_module_scores
from .trajectory_analysis import (
    cluster_sample_residuals,
    cluster_trajectory_profiles,
    summarize_residual_cluster_markers,
    trajectory_merge_diagnostics,
)

LOGGER = logging.getLogger(__name__)


def read_table(path: Path, *, optional: bool = False, **kwargs):
    """Read published tables, including valid empty diagnostic results."""
    if optional and not path.is_file():
        return None
    try:
        return pd.read_csv(path, **kwargs)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _sha256(path: Path) -> str:
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def read_diagnostics(directory: Path, pipeline: dict) -> dict:
    """Require diagnostics prepared with the report's configuration."""
    path = directory / "diagnostics.json"
    if not path.is_file():
        raise FileNotFoundError(f"{path}: run the DE diagnostics stage before rendering")
    metadata = read_json(path)
    if metadata.get("schema_version") != 1:
        raise ValueError(f"{path}: unsupported DE diagnostic schema")
    if metadata.get("configuration_sha256") != config_digest(pipeline):
        raise ValueError(f"{path}: configuration changed; regenerate DE diagnostics")
    return metadata


def _analysis_identity() -> dict:
    """Invalidate prepared calculations when their implementation or environment changes."""
    modules = (
        "abundance", "de_diagnostics", "fraction_age", "residual_wgcna",
        "trajectory_analysis", "go_enrichment", "covariates", "age_trajectory",
    )
    return {
        "code": {name: _sha256(Path(__file__).with_name(f"{name}.py")) for name in modules},
        "software": {name: version(name) for name in (
            "numpy", "pandas", "scipy", "scikit-learn", "scanpy", "PyWGCNA", "gseapy",
            "statsmodels", "igraph", "leidenalg", "umap-learn", "numba",
        )},
    }


def _enrichment_resources() -> dict:
    directory = Path(os.environ.get(GO_ENRICHMENT_DATA_DIR_ENV, "input_data/go"))
    return {name: _sha256(directory / name) if (directory / name).is_file() else None
            for name in (DEFAULT_GO_BP_CONFIG.gene_set_resource, DEFAULT_GO_BP_CONFIG.ontology_resource)
            if name is not None}


def _source_hashes(result_dir: Path, metadata: dict) -> dict[str, str]:
    names = {"run_metadata.json", "cell_type_result.json", "abundance_input.tsv"}
    for model in metadata.get("models", []):
        if model.get("status") == "complete":
            names.update(model.get("results_by_covariate", {}).values())
            for key in ("pearson_residuals", "pearson_residual_sample_metadata"):
                if model.get(key):
                    names.add(model[key])
    for name in ("combined/cmv.csv", "combined/smoking_status.csv", "combined/age_bin.csv"):
        if (result_dir / name).is_file():
            names.add(name)
    for name in names:
        if Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError(f"Unsafe DE source path: {name}")
    return {name: _sha256(result_dir / name) for name in sorted(names)}


def validate_diagnostic_sources(result_dir: Path, analysis_dir: Path) -> None:
    metadata = read_json(result_dir / "run_metadata.json")
    prepared = read_json(analysis_dir / "diagnostics.json")
    if prepared.get("sources") != _source_hashes(result_dir, metadata):
        raise ValueError("DE diagnostic inputs differ from the model results; regenerate diagnostics")


def read_wgcna(directory: Path) -> ResidualWGCNAResult:
    return ResidualWGCNAResult(
        metadata=read_json(directory / "metadata.json"),
        **{name: read_table(directory / f"{name}.csv") for name in (
            "soft_threshold", "module_assignments", "eigengenes", "module_summary",
        )},
    )


def copy_diagnostics(source: Path, destination: Path, pipeline: dict) -> None:
    """Validate and copy exactly the analysis artifacts used by a report."""
    metadata = read_diagnostics(source, pipeline)
    for name, expected_digest in metadata["artifacts"].items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe DE diagnostic artifact: {name}")
        path = source / relative
        if _sha256(path) != expected_digest:
            raise ValueError(f"{path}: DE diagnostic checksum mismatch")
    if source.resolve() == destination.resolve():
        return
    destination.mkdir(parents=True, exist_ok=True)
    previous_path = destination / "diagnostics.json"
    if previous_path.is_file():
        previous = read_json(previous_path)
        for name in previous.get("artifacts", {}).keys() - metadata["artifacts"].keys():
            relative = Path(name)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError(f"Unsafe previous DE diagnostic artifact: {name}")
            (destination / relative).unlink(missing_ok=True)
    for name in [*metadata["artifacts"], "diagnostics.json"]:
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / name, target)


def _save_table(table: pd.DataFrame, directory: Path, name: str, **kwargs) -> None:
    path = directory / name
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False, **kwargs)


def _save_enrichment(run, directory: Path, prefix: str, **kwargs) -> dict:
    """Keep full results and diagnostic status even if ontology summarization fails."""
    directory.mkdir(parents=True, exist_ok=True)
    _save_table(run.results, directory, f"{prefix}.csv")
    status = {"failures": run.failures, "summary_error": None}
    if not run.results.empty:
        try:
            save_go_specificity_outputs(run, output_dir=directory, prefix=prefix, **kwargs)
        except Exception as exc:  # noqa: BLE001 - retain the established optional-summary policy
            status["summary_error"] = f"{type(exc).__name__}: {exc}"
            LOGGER.warning("GO summary unavailable: %s", status["summary_error"])
    summary = directory / f"{prefix}_specific_summary.csv"
    if not summary.is_file():
        _save_table(pd.DataFrame(), directory, summary.name)
    (directory / f"{prefix}_status.json").write_text(json.dumps(status, indent=2) + "\n")
    return status


def _prepare_fraction_diagnostics(obs, directory, pipeline, cell_type, parent):
    settings = AgeTrajectorySettings.from_mapping(
        pipeline["differential_expression"]["age_trajectory"]
    )
    abundance = pipeline["differential_expression"]["abundance"]
    minimum = abundance["min_fraction_denominator_cells"]
    estimates, status = fit_fraction_age_bins(
        obs, settings=settings, min_denominator_cells=minimum,
        minimum_samples_per_bin=int(
            pipeline["differential_expression"]["study_comparison_minimum_samples_per_bin"]
        ),
    )
    estimates.insert(0, "cell_type", cell_type)
    estimates.insert(1, "aifi_l1_parent", parent)
    _save_table(estimates, directory, "fraction_age_bins.tsv", sep="\t")
    (directory / "fraction_age_bins_metadata.json").write_text(json.dumps(status, indent=2) + "\n")
    response = "fraction_within_aifi_l1_parent"
    adata = SimpleNamespace(obs=obs)
    fractions = sample_cell_type_fractions(
        adata, denominator="n_cells_in_sample_l1_parent", fraction_name=response,
        min_denominator_cells=minimum, count_column="n_cells",
    )
    model_settings = abundance["fraction_model"]
    coefficients, diagnostics = fit_parent_fraction_models(
        fractions, adata, fraction_column=response,
        bootstrap_replicates=model_settings["binomial_bootstrap_replicates"],
        bootstrap_seed=model_settings["binomial_bootstrap_seed"],
    )
    for table in (coefficients, diagnostics):
        table.insert(0, "cell_type", cell_type)
        table.insert(1, "aifi_l1_parent", parent)
        table.insert(2, "response", response)
    for key in ("binomial_bootstrap_replicates", "binomial_bootstrap_seed"):
        diagnostics[key] = model_settings[key]
    _save_table(coefficients, directory, "fraction_model_coefficients.tsv", sep="\t")
    _save_table(diagnostics, directory, "fraction_model_diagnostics.tsv", sep="\t")


def _prepare_trajectory_diagnostics(result_dir, directory, metadata, settings, cell_type):
    fit = next((model for model in metadata.get("models", [])
                if model.get("purpose") == "age_trajectory"), {})
    bins = metadata.get("age_trajectory", {}).get("retained_bins", [])
    path = result_dir / "combined/age_bin.csv"
    if fit.get("status") != "complete" or not path.is_file() or len(bins) < settings.minimum_bins:
        return fit
    results = pd.read_csv(path)
    omnibus = results.drop_duplicates("gene").copy()
    omnibus["omnibus_padj"] = pd.to_numeric(omnibus["omnibus_padj"], errors="coerce")
    significant = omnibus.loc[omnibus["omnibus_padj"] <= settings.cluster_fdr_threshold]
    if significant.empty:
        return fit
    profiles = results.loc[results["gene"].astype(str).isin(significant["gene"].astype(str))].pivot(
        index="gene", columns="age_bin", values="trajectory_scaled",
    ).reindex(columns=bins).dropna()
    labels, means = cluster_trajectory_profiles(profiles, settings)
    merges = trajectory_merge_diagnostics(profiles, settings, maximum_merges=20)
    rows = profiles.join(labels).reset_index()
    rows.insert(1, "cell_type", cell_type)
    rows["omnibus_padj"] = rows["gene"].map(omnibus.set_index("gene")["omnibus_padj"])
    means = means.reset_index()
    means.insert(0, "cell_type", cell_type)
    _save_table(rows, directory, "age_trajectory_clusters.csv")
    _save_table(means, directory, "age_trajectory_cluster_means.csv")
    _save_table(merges, directory, "age_trajectory_merge_diagnostics.csv")
    return fit


def _prepare_wgcna_diagnostics(result, residuals, annotations, directory):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "metadata.json").write_text(json.dumps(result.metadata, indent=2) + "\n")
    for name in ("soft_threshold", "module_assignments", "eigengenes", "module_summary"):
        _save_table(getattr(result, name), directory, f"{name}.csv")
    if result.metadata["status"] != "complete":
        return
    _save_table(wgcna_sample_connectivity(residuals, annotations), directory, "sample_connectivity.csv")
    if result.module_summary.empty:
        return
    associations = wgcna_covariate_associations(result, annotations)
    _save_table(associations, directory, "module_covariate_associations.csv")
    membership = wgcna_module_membership(result, residuals)
    hubs = membership.assign(
        absolute_module_membership=membership["module_membership"].abs(),
    ).sort_values(["module", "absolute_module_membership"], ascending=[True, False])
    hubs["module_rank"] = hubs.groupby("module", observed=True).cumcount() + 1
    _save_table(hubs, directory, "hub_genes.csv")
    for covariate, label, status_field in (
        ("cmv", "CMV", "cmv_status"), ("smoking_status", "smoking", "smoking_status"),
    ):
        scores = trait_associated_module_scores(result, annotations, covariate=covariate)
        if not scores.empty and {"umap_1", "umap_2"}.issubset(annotations):
            scores = scores.rename(columns={"trait_status": status_field}).merge(
                annotations[["sample_key", "umap_1", "umap_2"]].drop_duplicates("sample_key"),
                on="sample_key", how="left", validate="one_to_one",
            ).dropna(subset=["module_score", "umap_1", "umap_2"])
        else:
            scores = pd.DataFrame()
        _save_table(scores, directory, f"{label}_associated_module_scores.csv")
    _save_enrichment(
        run_wgcna_module_gsea(result, residuals), directory, "module_gsea",
        group_column="module", fdr_column="module_fdr",
    )


def _prepare_merged_age_enrichment(result_dir, stage, metadata, manifest):
    merged = next((model for model in metadata.get("models", [])
                   if model.get("model") == "merged" and model.get("status") == "complete"), None)
    if merged and "age" in merged.get("results_by_covariate", {}):
        results = pd.read_csv(result_dir / merged["results_by_covariate"]["age"])
        if "contrast" not in results:
            results["contrast"] = "per year"
        for index, (contrast, rows) in enumerate(results.groupby("contrast", observed=True)):
            plotted = rows[["padj", "log2FoldChange"]].apply(pd.to_numeric, errors="coerce").dropna()
            if plotted.empty:
                continue
            ranking = signed_de_gene_ranking(rows)
            if ranking.dropna().shape[0] < 10:
                continue
            directory = stage / "go_enrichment"
            prefix = "merged_age_gsea" if index == 0 else f"merged_age_gsea_{index}"
            status = _save_enrichment(
                run_preranked_gsea({f"Merged age ({contrast})": ranking}), directory, prefix,
            )
            status["summary"] = f"go_enrichment/{prefix}_specific_summary.csv"
            manifest["merged_age_enrichment"][str(contrast)] = status


def prepare_de_diagnostics(result_dir: Path, output_dir: Path, config_path: Path) -> dict:
    """Publish the same DE diagnostic calculations formerly executed by the notebook."""
    pipeline = load_analysis_configuration(config_path)
    metadata = read_json(result_dir / "run_metadata.json")
    sources = _source_hashes(result_dir, metadata)
    identity = _analysis_identity()
    resources = _enrichment_resources()
    existing = result_dir / "diagnostics.json"
    if existing.is_file():
        prepared = read_json(existing)
        if (prepared.get("configuration_sha256") == config_digest(pipeline)
                and prepared.get("analysis_identity") == identity
                and prepared.get("enrichment_resources") == resources
                and prepared.get("sources") == sources):
            try:
                copy_diagnostics(result_dir, output_dir, pipeline)
            except (FileNotFoundError, ValueError) as exc:
                LOGGER.warning("Regenerating invalid prepared diagnostics: %s", exc)
            else:
                LOGGER.info("Reusing prepared DE diagnostics for %s", metadata.get("cell_type"))
                return prepared
    record = read_json(result_dir / "cell_type_result.json")
    cell_type = str(metadata.get("cell_type") or record.get("cell_type") or result_dir.parent.name)
    obs = pd.read_csv(result_dir / "abundance_input.tsv", sep="\t")
    parents = obs["aifi_l1_parent_for_l2"].unique()
    if len(parents) != 1:
        raise ValueError(f"{cell_type!r} maps to multiple configured L1 parents")
    settings = AgeTrajectorySettings.from_mapping(pipeline["differential_expression"]["age_trajectory"])
    manifest = {"schema_version": 1, "cell_type": cell_type,
                "configuration_sha256": config_digest(pipeline), "residuals": None,
                "merged_age_enrichment": {}, "sources": sources, "analysis_identity": identity,
                "enrichment_resources": resources}
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".de-diagnostics-", dir=output_dir.parent) as temporary:
        stage = Path(temporary)
        LOGGER.info("Preparing fraction diagnostics for %s", cell_type)
        _prepare_fraction_diagnostics(obs, stage, pipeline, cell_type, str(parents[0]))
        _prepare_merged_age_enrichment(result_dir, stage, metadata, manifest)
        LOGGER.info("Preparing trajectory diagnostics for %s", cell_type)
        fit = _prepare_trajectory_diagnostics(result_dir, stage, metadata, settings, cell_type)
        annotations, pcs, n_components, residuals = load_residual_sample_clusters(result_dir, fit, settings)
        if annotations is not None:
            _save_table(annotations, stage, "residual_sample_clusters.csv")
            _save_table(pcs, stage, "residual_sample_pc_scores.csv")
            markers, _, reference = summarize_residual_cluster_markers(
                residuals, annotations, settings.residual_marker_max_clusters,
                progress=LOGGER.info,
            )
            _save_table(markers, stage, "residual_cluster_markers.csv")
            manifest["residuals"] = {"n_components": n_components, "shape": list(residuals.shape),
                                     "reference_cluster": int(reference) if reference is not None else None}
            LOGGER.info("Preparing residual WGCNA for %s", cell_type)
            network = run_residual_wgcna(
                residuals, annotations,
                settings=ResidualWGCNASettings.from_mapping(pipeline["cell_type_analysis"]["residual_wgcna"]),
                cmv_results=load_merged_cmv_results(result_dir),
                smoking_results=load_merged_smoking_results(result_dir),
            )
            _prepare_wgcna_diagnostics(network, residuals, annotations, stage / "residual_wgcna")
        manifest["artifacts"] = {
            str(path.relative_to(stage)): _sha256(path)
            for path in sorted(stage.rglob("*")) if path.is_file()
        }
        if sources != _source_hashes(result_dir, metadata) or resources != _enrichment_resources():
            raise ValueError("DE inputs changed while preparing diagnostics; retry from stable inputs")
        (stage / "diagnostics.json").write_text(json.dumps(manifest, indent=2) + "\n")
        copy_diagnostics(stage, output_dir, pipeline)
    return manifest


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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    configure_logging("INFO")
    prepare_de_diagnostics(args.result_dir, args.output_dir, args.config)


if __name__ == "__main__":
    main()
