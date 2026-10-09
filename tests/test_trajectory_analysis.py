from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from pbmc_pipeline.config import AgeTrajectorySettings, ExpressionAtlasSettings, read_json
from pbmc_pipeline.expression_atlas import analyze_expression_atlas
from pbmc_pipeline.trajectory_analysis import (
    analyze_age_trajectories,
    build_cross_cell_type_trajectories,
    cluster_trajectory_profiles,
)

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = AgeTrajectorySettings.from_mapping(
    read_json(ROOT / "config/pipeline.json")["differential_expression"]["age_trajectory"]
)
ATLAS_SETTINGS = ExpressionAtlasSettings.from_mapping(
    read_json(ROOT / "config/pipeline.json")["differential_expression"]["expression_atlas"]
)


def test_expression_atlas_uses_common_genes_depth_gate_and_technology_contrast(tmp_path) -> None:
    rows = []
    counts = []
    for study, technology, cell_type, values in [
        ("three_a", "10X3'v2", "Type A", [80, 20, 0]),
        ("three_b", "10X3'v2", "Type A", [60, 40, 0]),
        ("five_a", "10X5'v2", "Type A", [20, 80, 0]),
        ("five_b", "10X5'v2", "Type A", [30, 70, 0]),
        ("three_a", "10X3'v2", "Type B", [20, 80, 0]),
        ("three_b", "10X3'v2", "Type B", [30, 70, 0]),
        ("five_a", "10X5'v2", "Type B", [80, 20, 0]),
        ("five_b", "10X5'v2", "Type B", [60, 40, 0]),
    ]:
        rows.append({"study": study, "technology": technology, "aifi_l2_majority": cell_type})
        counts.append(values)
    pseudobulk = ad.AnnData(
        X=sparse.csr_matrix(np.asarray(counts)),
        obs=pd.DataFrame(rows),
        var=pd.DataFrame(
            {f"available_in_{study}": [True, True, True] for study in sorted({row["study"] for row in rows})},
            index=["G1", "G2", "G3"],
        ),
    )
    pseudobulk_path = tmp_path / "pseudobulk.h5ad"
    pseudobulk.write_h5ad(pseudobulk_path)
    output_dir = tmp_path / "atlas"

    metadata = analyze_expression_atlas(
        pseudobulk_path,
        output_dir,
        replace(
            ATLAS_SETTINGS,
            minimum_study_cell_type_total_counts=100,
            report_top_n_genes=10,
        ),
        split_by="aifi_l2_majority",
        l2_parent_l1={"Type A": "L1 A", "Type B": "L1 B"},
    )

    assert metadata["status"] == "complete"
    assert metadata["n_depth_qualified_groups"] == 8
    assert metadata["l1_parent_by_l2"] == {"Type A": "L1 A", "Type B": "L1 B"}
    matrix = pd.read_csv(output_dir / "expression_atlas_matrix.csv", index_col="gene")
    expected = np.mean(np.log2(np.array([800_000, 600_000, 200_000, 300_000]) + 1))
    assert np.isclose(matrix.loc["G1", "Type A"], expected)
    support = pd.read_csv(output_dir / "expression_atlas_study_support.csv")
    assert support["retained_by_depth_gate"].all()
    availability = pd.read_csv(output_dir / "expression_atlas_gene_availability.csv", index_col="gene")
    assert availability.columns.tolist() == ["five_a", "five_b", "three_a", "three_b"]
    assert availability.astype(bool).all(axis=None)
    intersection_accounting = pd.read_csv(
        output_dir / "expression_atlas_gene_intersection_accounting.csv"
    )
    assert intersection_accounting["discarded_counts"].eq(0).all()
    assert intersection_accounting["top_discarded_gene"].isna().all()
    assert metadata["n_genes_in_outer_union"] == 3
    assert metadata["gene_availability_by_study"]["three_a"] == 3
    assert metadata["gene_intersection_accounting"]["discarded_counts"] == 0
    contrast = pd.read_csv(output_dir / "expression_atlas_technology_contrast.csv")
    assert (contrast.loc[(contrast["gene"] == "G1") & (contrast["cell_type"] == "Type A"), "log2_cpm_difference_5_prime_minus_3_prime"] < 0).all()
    assert (output_dir / "expression_atlas_gene_clusters.csv").is_file()

    skipped_output_dir = tmp_path / "skipped_atlas"
    skipped = analyze_expression_atlas(
        pseudobulk_path,
        skipped_output_dir,
        replace(ATLAS_SETTINGS, minimum_study_cell_type_total_counts=10_000),
        split_by="aifi_l2_majority",
    )
    assert skipped["status"] == "skipped"
    assert pd.read_csv(skipped_output_dir / "expression_atlas_gene_clusters.csv").empty
    assert pd.read_csv(skipped_output_dir / "expression_atlas_technology_contrast.csv").empty

    missing_technology = pseudobulk.copy()
    missing_technology.obs = missing_technology.obs.drop(columns=["technology"])
    missing_technology_path = tmp_path / "pseudobulk_without_technology.h5ad"
    missing_technology.write_h5ad(missing_technology_path)
    unknown_technology_output_dir = tmp_path / "unknown_technology_atlas"
    unknown_technology = analyze_expression_atlas(
        missing_technology_path,
        unknown_technology_output_dir,
        replace(ATLAS_SETTINGS, minimum_study_cell_type_total_counts=100),
        split_by="aifi_l2_majority",
    )
    assert unknown_technology["status"] == "complete"
    assert unknown_technology["technology_available"] is False
    assert pd.read_csv(
        unknown_technology_output_dir / "expression_atlas_technology_contrast.csv"
    ).empty


def test_expression_atlas_records_counts_lost_by_shared_gene_intersection(tmp_path) -> None:
    pseudobulk = ad.AnnData(
        X=sparse.csr_matrix(np.array([[10, 4, 0], [20, 0, 8]])),
        obs=pd.DataFrame({
            "study": ["study_a", "study_b"],
            "aifi_l2_majority": ["Type A", "Type A"],
        }),
        var=pd.DataFrame({
            "available_in_study_a": [True, True, False],
            "available_in_study_b": [True, False, True],
        }, index=["SHARED", "ONLY_A", "ONLY_B"]),
    )
    pseudobulk_path = tmp_path / "pseudobulk.h5ad"
    pseudobulk.write_h5ad(pseudobulk_path)
    output_dir = tmp_path / "atlas"

    metadata = analyze_expression_atlas(
        pseudobulk_path,
        output_dir,
        replace(ATLAS_SETTINGS, minimum_study_cell_type_total_counts=1),
        split_by="aifi_l2_majority",
    )

    accounting = pd.read_csv(
        output_dir / "expression_atlas_gene_intersection_accounting.csv"
    ).set_index("study")
    assert metadata["n_shared_genes"] == 1
    assert accounting.loc["study_a", "discarded_counts"] == 4
    assert accounting.loc["study_a", "top_discarded_gene"] == "ONLY_A"
    assert accounting.loc["study_b", "discarded_counts"] == 8
    assert accounting.loc["study_b", "top_discarded_gene"] == "ONLY_B"
    assert metadata["gene_intersection_accounting"]["discarded_counts"] == 12


def test_hierarchical_trajectory_clusters_and_umap_use_complete_profiles() -> None:
    profiles = pd.DataFrame(
        [
            [0.0, 1.0, -1.0, 1.0, -1.0],
            [0.0, 1.1, -1.1, 1.0, -1.0],
            [0.0, -1.0, 1.0, -1.0, 1.0],
            [0.0, -1.1, 1.1, -1.0, 1.0],
        ],
        index=pd.Index(["up_a", "up_b", "down_a", "down_b"], name="gene"),
        columns=["20-30", "30-40", "40-50", "50-60", "60-70"],
    )

    annotations, means = cluster_trajectory_profiles(
        profiles, replace(SETTINGS, max_clusters=2)
    )

    assert annotations["trajectory_cluster"].nunique() == 2
    assert annotations.loc["up_a", "trajectory_cluster"] == annotations.loc[
        "up_b", "trajectory_cluster"
    ]
    assert annotations.loc["down_a", "trajectory_cluster"] == annotations.loc[
        "down_b", "trajectory_cluster"
    ]
    assert np.isfinite(annotations[["umap_1", "umap_2"]].to_numpy()).all()
    assert means["n_trajectories"].sum() == len(profiles)


def test_cross_cell_type_analysis_uses_shared_bins_and_counts_de_recurrence() -> None:
    bins = [
        "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90",
    ]
    results = {}
    metadata = {}
    for cell_type, retained in {
        "Type A": bins,
        "Type B": [*bins, "90-100"],
    }.items():
        metadata[cell_type] = {
            "cell_type": cell_type,
            "age_trajectory": {
                "status": "complete",
                "retained_bins": retained,
            },
        }
        rows = []
        for gene, fdr, amplitude in [
            ("G_CLUSTER", 0.0005, 1.0),
            ("G_RECURRENCE", 0.01, -0.8),
            ("G_NOT_DE", 0.2, 0.5),
        ]:
            for index, age_bin in enumerate(retained):
                rows.append({
                    "gene": gene,
                    "age_bin": age_bin,
                    "log2FoldChange": 0.0 if index == 0 else amplitude * (-1) ** index,
                    "omnibus_padj": fdr,
                })
        results[cell_type] = pd.DataFrame(rows)

    metadata["Low-support type"] = {
        "cell_type": "Low-support type",
        "age_trajectory": {"status": "complete", "retained_bins": bins[:6]},
    }
    results["Low-support type"] = results["Type A"].loc[
        results["Type A"]["age_bin"].isin(bins[:6])
    ].copy()

    analysis = build_cross_cell_type_trajectories(results, metadata, settings=SETTINGS)

    assert analysis["common_bins"] == bins
    assert analysis["supported_cell_types"] == ["Type A", "Type B"]
    assert analysis["eligible_cell_types"] == ["Type A", "Type B"]
    assert "Low-support type" in analysis["skipped_cell_types"]
    assert "Low-support type" not in set(analysis["significance"]["cell_type"])
    assert set(analysis["trajectories"]["gene"]) == {"G_CLUSTER"}
    recurrence = analysis["gene_recurrence"].set_index("gene")
    assert recurrence.loc["G_CLUSTER", "n_cell_types_DE_significant"] == 2
    assert recurrence.loc["G_RECURRENCE", "n_cell_types_DE_significant"] == 2
    assert recurrence.loc["G_NOT_DE", "n_cell_types_DE_significant"] == 0


def test_completed_fit_with_fewer_than_seven_bins_is_omitted_from_report_data(tmp_path) -> None:
    de_dir = tmp_path / "differential_expression"
    cell_type_dir = de_dir / "rare-type"
    combined_dir = cell_type_dir / "combined"
    combined_dir.mkdir(parents=True)
    retained_bins = [
        "20-30", "30-40", "40-50", "50-60", "60-70", "70-80",
    ]
    (cell_type_dir / "run_metadata.json").write_text(json.dumps({
        "cell_type": "Rare type",
        "models": [{"purpose": "age_trajectory", "status": "complete"}],
        "age_trajectory": {
            "status": "complete",
            "retained_bins": retained_bins,
        },
    }))
    pd.DataFrame({
        "gene": ["GENE1"] * len(retained_bins),
        "age_bin": retained_bins,
        "log2FoldChange": [0.0, 0.3, 0.9, 0.4, -0.3, -0.7],
        "trajectory_scaled": [0.0, 0.5, 1.2, 0.4, -0.5, -1.0],
        "omnibus_padj": [0.0005] * len(retained_bins),
    }).to_csv(combined_dir / "age_bin.csv", index=False)

    metadata = analyze_age_trajectories(
        de_dir, tmp_path / "trajectory_analysis", SETTINGS
    )

    assert metadata["cell_types_eligible_by_bin_support"] == []
    assert metadata["cell_types_in_cross_cell_type_analysis"] == []
    assert "Rare type" in metadata["cell_types_excluded_from_cross_cell_type_analysis"]
    assert pd.read_csv(
        tmp_path / "trajectory_analysis" / "gene_recurrence.csv"
    ).empty
    assert pd.read_csv(
        tmp_path / "trajectory_analysis" / "cross_cell_type_trajectory_clusters.csv"
    ).empty
    assert not (
        tmp_path / "trajectory_analysis" / "cell_type_trajectory_clusters.csv"
    ).exists()
