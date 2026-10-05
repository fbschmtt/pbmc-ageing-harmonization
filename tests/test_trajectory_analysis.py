from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from pbmc_pipeline.config import AgeTrajectorySettings, read_json
from pbmc_pipeline.trajectory_analysis import (
    analyze_age_trajectories,
    build_cross_cell_type_trajectories,
    cluster_trajectory_profiles,
)

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = AgeTrajectorySettings.from_mapping(
    read_json(ROOT / "config/pipeline.json")["differential_expression"]["age_trajectory"]
)


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
