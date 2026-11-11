import subprocess

from scripts import image_context_hash as context


def test_context_hash_changes_for_copied_sources_and_ignored_optional_color_table(tmp_path, monkeypatch):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    monkeypatch.setattr(context, "ROOT", tmp_path)
    (tmp_path / ".gitignore").write_text("aifi_models/*.csv\noutput/\n")
    (tmp_path / "src").mkdir()
    source = tmp_path / "src" / "analysis.py"
    source.write_text("first")
    first = context.context_hash("python")
    (tmp_path / "output").mkdir()
    (tmp_path / "output" / "generated.txt").write_text("not copied into the image")
    assert context.context_hash("python") == first
    source.write_text("second")
    second = context.context_hash("python")
    assert second != first
    (tmp_path / "aifi_models").mkdir()
    (tmp_path / "aifi_models" / "AIFI_L2_imm_health_atlas_type_order_colors.csv").write_text("label,color\n")
    assert context.context_hash("python") != second


def test_site_context_excludes_generated_files_and_tracks_hierarchy(tmp_path, monkeypatch):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    monkeypatch.setattr(context, "ROOT", tmp_path)
    (tmp_path / ".gitignore").write_text("site/node_modules/\nsite/data-inputs/\n")
    for directory in ("site/node_modules", "site/data-inputs", "site/src"):
        (tmp_path / directory).mkdir(parents=True)
    (tmp_path / "site/src/App.tsx").write_text("app")
    dependencies = tmp_path / "site/node_modules/package.json"
    dependencies.write_text("first")
    hierarchy = tmp_path / "site/data-inputs/cell-type-hierarchy.json"
    hierarchy.write_text("first")
    first = context.context_hash("site")
    dependencies.write_text("second")
    assert context.context_hash("site") == first
    hierarchy.write_text("second")
    assert context.context_hash("site") != first
