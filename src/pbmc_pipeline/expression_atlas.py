"""Build descriptive, study-balanced expression summaries across cell types."""
from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

from .config import ExpressionAtlasSettings


def _technology_family(value: object) -> str | None:
    normalized = str(value).upper().replace(" ", "")
    if normalized.startswith("10X3'"):
        return "3_prime"
    if normalized.startswith("10X5'"):
        return "5_prime"
    return None


def _write_table(path: Path, table: pd.DataFrame, *, index: bool = False) -> None:
    table.to_csv(path, index=index)


def _empty_atlas_tables(output_dir: Path, expression_columns: list[str] | None = None) -> None:
    """Write header-only atlas artifacts when no group passes the depth gate."""
    columns = expression_columns or []
    _write_table(
        output_dir / "expression_atlas_gene_clusters.csv",
        pd.DataFrame(
            columns=[
                "gene",
                "expression_cluster",
                "profile_sd",
                *[f"z_{column}" for column in columns],
                "n_cell_types_at_minimum_cpm",
            ]
        ),
    )
    _write_table(
        output_dir / "expression_atlas_cluster_means.csv",
        pd.DataFrame(columns=["expression_cluster", *columns, "n_genes"]),
    )
    _write_table(
        output_dir / "expression_atlas_gene_umap.csv",
        pd.DataFrame(columns=["gene", "umap_1", "umap_2", "scaling_sd_log2_cpm"]),
    )
    _write_table(
        output_dir / "expression_atlas_technology_contrast.csv",
        pd.DataFrame(
            columns=[
                "gene",
                "cell_type",
                "log2_cpm_3_prime",
                "log2_cpm_5_prime",
                "log2_cpm_difference_5_prime_minus_3_prime",
                "n_studies_3_prime",
                "n_studies_5_prime",
                "meets_configured_study_replication",
            ]
        ),
    )
    _write_table(
        output_dir / "expression_atlas_technology_summary.csv",
        pd.DataFrame(
            columns=[
                "gene",
                "median_log2_cpm_difference_5_prime_minus_3_prime",
                "mean_log2_cpm_across_technologies",
                "n_cell_types_with_technology_contrast",
                "n_cell_types_meeting_configured_study_replication",
                "spearman_rho_cell_type_profile_3_prime_vs_5_prime",
            ]
        ),
    )
    _write_table(
        output_dir / "expression_atlas_intronic_contrast.csv",
        pd.DataFrame(
            columns=[
                "gene",
                "cell_type",
                "log2_cpm_intronic",
                "log2_cpm_non_intronic",
                "log2_cpm_difference_intronic_minus_non_intronic",
                "n_studies_intronic",
                "n_studies_non_intronic",
                "meets_configured_study_replication",
            ]
        ),
    )
    _write_table(
        output_dir / "expression_atlas_intronic_summary.csv",
        pd.DataFrame(
            columns=[
                "gene",
                "median_log2_cpm_difference_intronic_minus_non_intronic",
                "mean_log2_cpm_across_intronic_status",
                "n_cell_types_with_intronic_contrast",
                "n_cell_types_meeting_configured_study_replication",
            ]
        ),
    )


def _cluster_expression_profiles(
    expression: pd.DataFrame, settings: ExpressionAtlasSettings
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Cluster row-standardized cell-type expression profiles."""
    values = expression.to_numpy(dtype=float)
    centered = values - values.mean(axis=1, keepdims=True)
    standard_deviation = centered.std(axis=1, keepdims=True)
    valid = np.isfinite(values).all(axis=1) & (standard_deviation[:, 0] > 0)
    z_scores = np.full_like(values, np.nan, dtype=float)
    z_scores[valid] = centered[valid] / standard_deviation[valid]
    clusterable = pd.DataFrame(z_scores[valid], index=expression.index[valid], columns=expression.columns)
    if clusterable.empty:
        assignments = pd.DataFrame(columns=["gene", "expression_cluster", "profile_sd"])
        means = pd.DataFrame(columns=["expression_cluster", *expression.columns, "n_genes"])
        return assignments, means
    if len(clusterable) == 1:
        labels = np.ones(1, dtype=int)
    else:
        tree = linkage(clusterable.to_numpy(), method="average", metric="correlation")
        labels = fcluster(tree, t=min(settings.max_clusters, len(clusterable)), criterion="maxclust")
    assignments = pd.DataFrame({
        "gene": clusterable.index.astype(str),
        "expression_cluster": labels,
        "profile_sd": standard_deviation[valid, 0],
        **{f"z_{column}": clusterable[column].to_numpy() for column in clusterable.columns},
    })
    means = (
        assignments.groupby("expression_cluster", sort=True)[[f"z_{column}" for column in expression.columns]]
        .mean()
        .rename(columns={f"z_{column}": column for column in expression.columns})
    )
    means["n_genes"] = assignments.groupby("expression_cluster", sort=True).size()
    return assignments, means.reset_index()


def _expression_umap(
    expression: pd.DataFrame, settings: ExpressionAtlasSettings
) -> tuple[pd.DataFrame, int]:
    """Embed genes after per-gene centering/scaling solely for PCA and UMAP."""
    values = expression.to_numpy(dtype=float)
    centered = values - values.mean(axis=1, keepdims=True)
    standard_deviation = centered.std(axis=1, keepdims=True)
    valid = np.isfinite(values).all(axis=1) & (standard_deviation[:, 0] > 0)
    genes = expression.index[valid].astype(str)
    if len(genes) < 4:
        return pd.DataFrame(
            columns=["gene", "umap_1", "umap_2", "scaling_sd_log2_cpm"]
        ), 0
    normalized = centered[valid] / standard_deviation[valid]
    n_components = min(
        settings.umap_pca_components, normalized.shape[0] - 1, normalized.shape[1]
    )
    scores = PCA(n_components=n_components, random_state=settings.random_state).fit_transform(normalized)
    import umap

    embedding = umap.UMAP(
        n_components=2,
        n_neighbors=min(settings.umap_neighbors, len(genes) - 1),
        min_dist=settings.umap_min_dist,
        metric="euclidean",
        random_state=settings.random_state,
        transform_seed=settings.random_state,
        n_jobs=1,
    ).fit_transform(scores)
    return pd.DataFrame({
        "gene": genes,
        "umap_1": embedding[:, 0],
        "umap_2": embedding[:, 1],
        "scaling_sd_log2_cpm": standard_deviation[valid, 0],
    }), n_components


def analyze_expression_atlas(
    pseudobulk_path: Path,
    output_dir: Path,
    settings: ExpressionAtlasSettings,
    *,
    split_by: str,
    l2_parent_l1: Mapping[str, str] | None = None,
) -> dict:
    """Write a common-gene, depth-gated expression atlas and technology sensitivity tables."""
    import anndata as ad

    adata = ad.read_h5ad(pseudobulk_path)
    required = {"study", split_by}
    missing = sorted(required - set(adata.obs.columns))
    if missing:
        raise ValueError(f"{pseudobulk_path}: expression atlas metadata columns are missing: {missing}")
    studies = sorted(adata.obs["study"].astype(str).unique())
    technology_available = "technology" in adata.obs.columns
    intronic_available = "include_intronic" in adata.obs.columns
    availability_columns = [f"available_in_{study}" for study in studies]
    missing_availability = sorted(set(availability_columns) - set(adata.var.columns))
    if missing_availability:
        raise ValueError(
            f"{pseudobulk_path}: expression atlas gene-availability columns are missing: "
            f"{missing_availability}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_availability = adata.var[availability_columns].fillna(False).astype(bool).copy()
    gene_availability.columns = studies
    gene_availability.index = pd.Index(adata.var_names.astype(str), name="gene")
    _write_table(output_dir / "expression_atlas_gene_availability.csv", gene_availability, index=True)
    shared_mask = gene_availability.all(axis=1).to_numpy()
    genes = adata.var_names[shared_mask].astype(str)
    empty_matrix = pd.DataFrame(index=pd.Index([], name="gene"))
    intersection_rows = []
    total_input_counts = 0
    total_retained_counts = 0
    for study in studies:
        selected = np.flatnonzero(adata.obs["study"].astype(str).to_numpy() == study)
        gene_counts = np.asarray(adata.X[selected].sum(axis=0)).ravel().astype(np.int64)
        available = gene_availability[study].to_numpy(dtype=bool)
        discarded = available & ~shared_mask
        input_counts = int(gene_counts.sum(dtype=np.int64))
        retained_counts = int(gene_counts[shared_mask].sum(dtype=np.int64))
        discarded_counts = input_counts - retained_counts
        top_discarded_gene = None
        top_discarded_gene_counts = 0
        if discarded_counts:
            discarded_positions = np.flatnonzero(discarded)
            top_position = discarded_positions[np.argmax(gene_counts[discarded])]
            top_discarded_gene = str(adata.var_names[top_position])
            top_discarded_gene_counts = int(gene_counts[top_position])
        intersection_rows.append({
            "study": study,
            "input_gene_symbols": int(available.sum()),
            "retained_gene_symbols": int(shared_mask.sum()),
            "discarded_gene_symbols": int(discarded.sum()),
            "input_counts": input_counts,
            "retained_counts": retained_counts,
            "discarded_counts": discarded_counts,
            "discarded_fraction": (
                discarded_counts / input_counts if input_counts else 0.0
            ),
            "top_discarded_gene": top_discarded_gene,
            "top_discarded_gene_counts": top_discarded_gene_counts,
        })
        total_input_counts += input_counts
        total_retained_counts += retained_counts
    intersection_accounting = pd.DataFrame(intersection_rows)
    _write_table(
        output_dir / "expression_atlas_gene_intersection_accounting.csv",
        intersection_accounting,
    )
    total_discarded_counts = total_input_counts - total_retained_counts
    availability_metadata = {
        "n_genes_in_outer_union": int(adata.n_vars),
        "n_shared_genes": len(genes),
        "gene_availability_by_study": {
            study: int(gene_availability[study].sum()) for study in studies
        },
        "gene_intersection_accounting": {
            "input_counts": total_input_counts,
            "retained_counts": total_retained_counts,
            "discarded_counts": total_discarded_counts,
            "discarded_fraction": (
                total_discarded_counts / total_input_counts if total_input_counts else 0.0
            ),
        },
    }
    if not len(genes):
        _write_table(output_dir / "expression_atlas_matrix.csv", empty_matrix, index=True)
        _empty_atlas_tables(output_dir)
        return {
            "status": "skipped",
            "reason": "no genes are shared by every study",
            "technology_available": technology_available,
            "intronic_available": intronic_available,
            **availability_metadata,
            **settings.to_mapping(),
        }

    counts = adata.X[:, shared_mask]
    obs = adata.obs.reset_index(drop=False).copy()
    grouped_records = []
    group_profiles: dict[tuple[str, str], np.ndarray] = {}
    for (study, cell_type), positions in obs.groupby(["study", split_by], observed=True).indices.items():
        selected = np.asarray(positions, dtype=int)
        profile = np.asarray(counts[selected].sum(axis=0)).ravel().astype(float)
        total_counts = float(profile.sum())
        if technology_available:
            technologies = sorted(set(obs.iloc[selected]["technology"].dropna().astype(str)))
            technology = technologies[0] if len(technologies) == 1 else "heterogeneous_during_bulk"
            family = _technology_family(technology)
        else:
            technology = "unknown"
            family = None
        if intronic_available:
            intronic_values = sorted(
                set(obs.iloc[selected]["include_intronic"].dropna().astype(str).str.lower())
            )
            intronic_status = (
                intronic_values[0] if len(intronic_values) == 1 else "heterogeneous_during_bulk"
            )
        else:
            intronic_status = "unknown"
        retained = total_counts >= settings.minimum_study_cell_type_total_counts
        key = (str(study), str(cell_type))
        group_profiles[key] = profile
        grouped_records.append({
            "study": str(study),
            "cell_type": str(cell_type),
            "technology": technology,
            "technology_family": family or "other_or_mixed",
            "include_intronic": intronic_status,
            "total_counts_in_gene_intersection": int(total_counts),
            "retained_by_depth_gate": retained,
        })
    support = pd.DataFrame(grouped_records).sort_values(["cell_type", "study"], kind="stable")
    _write_table(output_dir / "expression_atlas_study_support.csv", support)

    def study_balanced_matrix(records: pd.DataFrame, minimum_studies: int) -> tuple[pd.DataFrame, pd.Series]:
        vectors: dict[str, list[np.ndarray]] = {}
        study_counts: dict[str, int] = {}
        for cell_type, rows in records.groupby("cell_type", observed=True):
            unique_studies = rows["study"].nunique()
            if unique_studies < minimum_studies:
                continue
            vectors[str(cell_type)] = [
                np.log2(group_profiles[(str(row.study), str(row.cell_type))] / row.total_counts_in_gene_intersection * 1_000_000 + settings.cpm_pseudocount)
                for row in rows.itertuples(index=False)
            ]
            study_counts[str(cell_type)] = int(unique_studies)
        if not vectors:
            return empty_matrix.copy(), pd.Series(dtype="int64", name="n_studies")
        matrix = pd.DataFrame(
            {cell_type: np.mean(values, axis=0) for cell_type, values in sorted(vectors.items())},
            index=pd.Index(genes, name="gene"),
        )
        return matrix, pd.Series(study_counts, name="n_studies").reindex(matrix.columns)

    retained = support.loc[support["retained_by_depth_gate"]].copy()
    expression, _study_counts = study_balanced_matrix(
        retained, settings.minimum_studies_per_cell_type
    )
    _write_table(output_dir / "expression_atlas_matrix.csv", expression, index=True)

    if expression.empty:
        _empty_atlas_tables(output_dir)
        return {
            "status": "skipped",
            "reason": "no cell type has the configured number of depth-qualified studies",
            "n_shared_genes": len(genes),
            "n_study_cell_type_groups": len(support),
            "n_depth_qualified_groups": len(retained),
            "technology_available": technology_available,
            "technology_contrast_reason": (
                None if technology_available else "technology metadata was not present in the pseudobulk input"
            ),
            "intronic_contrast_reason": (
                None if intronic_available else "intronic-read metadata was not present in the pseudobulk input"
            ),
            **availability_metadata,
            **settings.to_mapping(),
        }

    cpm = np.exp2(expression) - settings.cpm_pseudocount
    signal = (cpm >= settings.minimum_cpm_for_clustering).sum(axis=1)
    clusters, cluster_means = _cluster_expression_profiles(
        expression.loc[signal > 0], settings
    )
    if not clusters.empty:
        clusters = clusters.merge(signal.rename("n_cell_types_at_minimum_cpm"), left_on="gene", right_index=True)
    _write_table(output_dir / "expression_atlas_gene_clusters.csv", clusters)
    _write_table(output_dir / "expression_atlas_cluster_means.csv", cluster_means)
    gene_umap, n_expression_umap_pcs = _expression_umap(
        expression.loc[signal > 0], settings
    )
    _write_table(output_dir / "expression_atlas_gene_umap.csv", gene_umap)
    parent_map = l2_parent_l1 or {}
    l1_parent_by_l2 = {
        cell_type: str(parent_map[cell_type])
        for cell_type in expression.columns
        if cell_type in parent_map
    }
    unmapped_l2_cell_types = sorted(set(expression.columns) - set(l1_parent_by_l2))

    by_technology = {}
    technology_study_counts = {}
    for family in ("3_prime", "5_prime"):
        by_technology[family], technology_study_counts[family] = study_balanced_matrix(
            retained.loc[retained["technology_family"] == family],
            1,
        )
    shared_types = by_technology["3_prime"].columns.intersection(by_technology["5_prime"].columns)
    contrast_rows = []
    if len(shared_types):
        three_prime = by_technology["3_prime"].loc[:, shared_types]
        five_prime = by_technology["5_prime"].loc[:, shared_types]
        difference = five_prime - three_prime
        for cell_type in shared_types:
            n_three = int(technology_study_counts["3_prime"].loc[cell_type])
            n_five = int(technology_study_counts["5_prime"].loc[cell_type])
            contrast_rows.append(pd.DataFrame({
                "gene": genes,
                "cell_type": cell_type,
                "log2_cpm_3_prime": three_prime[cell_type].to_numpy(),
                "log2_cpm_5_prime": five_prime[cell_type].to_numpy(),
                "log2_cpm_difference_5_prime_minus_3_prime": difference[cell_type].to_numpy(),
                "n_studies_3_prime": n_three,
                "n_studies_5_prime": n_five,
                "meets_configured_study_replication": (
                    n_three >= settings.minimum_studies_per_technology
                    and n_five >= settings.minimum_studies_per_technology
                ),
            }))
        contrast = pd.concat(contrast_rows, ignore_index=True)
        summaries = []
        for gene, values in difference.iterrows():
            paired = pd.DataFrame({"three": three_prime.loc[gene], "five": five_prime.loc[gene]}).dropna()
            correlation = spearmanr(paired["three"], paired["five"]).statistic if len(paired) >= 3 else np.nan
            summaries.append({
                "gene": gene,
                "median_log2_cpm_difference_5_prime_minus_3_prime": float(values.median()),
                "mean_log2_cpm_across_technologies": float(pd.concat([three_prime.loc[gene], five_prime.loc[gene]]).mean()),
                "n_cell_types_with_technology_contrast": len(paired),
                "n_cell_types_meeting_configured_study_replication": int(sum(
                    technology_study_counts["3_prime"].loc[cell_type]
                    >= settings.minimum_studies_per_technology
                    and technology_study_counts["5_prime"].loc[cell_type]
                    >= settings.minimum_studies_per_technology
                    for cell_type in shared_types
                )),
                "spearman_rho_cell_type_profile_3_prime_vs_5_prime": correlation,
            })
        technology_summary = pd.DataFrame(summaries)
    else:
        contrast = pd.DataFrame(columns=[
            "gene", "cell_type", "log2_cpm_3_prime", "log2_cpm_5_prime",
            "log2_cpm_difference_5_prime_minus_3_prime",
            "n_studies_3_prime", "n_studies_5_prime",
            "meets_configured_study_replication",
        ])
        technology_summary = pd.DataFrame(columns=[
            "gene", "median_log2_cpm_difference_5_prime_minus_3_prime",
            "mean_log2_cpm_across_technologies", "n_cell_types_with_technology_contrast",
            "n_cell_types_meeting_configured_study_replication",
            "spearman_rho_cell_type_profile_3_prime_vs_5_prime",
        ])
    _write_table(output_dir / "expression_atlas_technology_contrast.csv", contrast)
    _write_table(output_dir / "expression_atlas_technology_summary.csv", technology_summary)

    by_intronic_status = {}
    intronic_study_counts = {}
    for status in ("yes", "no"):
        by_intronic_status[status], intronic_study_counts[status] = study_balanced_matrix(
            retained.loc[retained["include_intronic"] == status], 1
        )
    shared_intronic_types = (
        by_intronic_status["yes"].columns.intersection(by_intronic_status["no"].columns)
    )
    intronic_rows = []
    if len(shared_intronic_types):
        intronic = by_intronic_status["yes"].loc[:, shared_intronic_types]
        non_intronic = by_intronic_status["no"].loc[:, shared_intronic_types]
        intronic_difference = intronic - non_intronic
        for cell_type in shared_intronic_types:
            n_intronic = int(intronic_study_counts["yes"].loc[cell_type])
            n_non_intronic = int(intronic_study_counts["no"].loc[cell_type])
            intronic_rows.append(pd.DataFrame({
                "gene": genes,
                "cell_type": cell_type,
                "log2_cpm_intronic": intronic[cell_type].to_numpy(),
                "log2_cpm_non_intronic": non_intronic[cell_type].to_numpy(),
                "log2_cpm_difference_intronic_minus_non_intronic": (
                    intronic_difference[cell_type].to_numpy()
                ),
                "n_studies_intronic": n_intronic,
                "n_studies_non_intronic": n_non_intronic,
                "meets_configured_study_replication": (
                    n_intronic >= settings.minimum_studies_per_technology
                    and n_non_intronic >= settings.minimum_studies_per_technology
                ),
            }))
        intronic_contrast = pd.concat(intronic_rows, ignore_index=True)
        intronic_summary = pd.DataFrame([
            {
                "gene": gene,
                "median_log2_cpm_difference_intronic_minus_non_intronic": float(values.median()),
                "mean_log2_cpm_across_intronic_status": float(
                    pd.concat([intronic.loc[gene], non_intronic.loc[gene]]).mean()
                ),
                "n_cell_types_with_intronic_contrast": len(shared_intronic_types),
                "n_cell_types_meeting_configured_study_replication": int(sum(
                    intronic_study_counts["yes"].loc[cell_type]
                    >= settings.minimum_studies_per_technology
                    and intronic_study_counts["no"].loc[cell_type]
                    >= settings.minimum_studies_per_technology
                    for cell_type in shared_intronic_types
                )),
            }
            for gene, values in intronic_difference.iterrows()
        ])
    else:
        intronic_contrast = pd.DataFrame(columns=[
            "gene", "cell_type", "log2_cpm_intronic", "log2_cpm_non_intronic",
            "log2_cpm_difference_intronic_minus_non_intronic", "n_studies_intronic",
            "n_studies_non_intronic", "meets_configured_study_replication",
        ])
        intronic_summary = pd.DataFrame(columns=[
            "gene", "median_log2_cpm_difference_intronic_minus_non_intronic",
            "mean_log2_cpm_across_intronic_status", "n_cell_types_with_intronic_contrast",
            "n_cell_types_meeting_configured_study_replication",
        ])
    _write_table(output_dir / "expression_atlas_intronic_contrast.csv", intronic_contrast)
    _write_table(output_dir / "expression_atlas_intronic_summary.csv", intronic_summary)
    return {
        "status": "complete",
        **availability_metadata,
        "n_study_cell_type_groups": len(support),
        "n_depth_qualified_groups": len(retained),
        "cell_types": expression.columns.tolist(),
        "n_cell_types": int(expression.shape[1]),
        "n_clustered_genes": len(clusters),
        "n_gene_expression_umap_genes": len(gene_umap),
        "n_gene_expression_umap_pcs": n_expression_umap_pcs,
        "l1_parent_by_l2": l1_parent_by_l2,
        "unmapped_l2_cell_types": unmapped_l2_cell_types,
        "technology_contrast_cell_types": shared_types.tolist(),
        "technology_available": technology_available,
        "technology_contrast_reason": (
            None if technology_available else "technology metadata was not present in the pseudobulk input"
        ),
        "intronic_contrast_cell_types": shared_intronic_types.tolist(),
        "intronic_available": intronic_available,
        "intronic_contrast_reason": (
            None if intronic_available else "intronic-read metadata was not present in the pseudobulk input"
        ),
        **settings.to_mapping(),
    }
