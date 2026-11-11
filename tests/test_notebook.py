import json
import os
import subprocess
from pathlib import Path

import pytest

from pbmc_pipeline.notebook import render_notebook


def test_notebook_lifecycle_isolates_all_tool_phases_without_changing_parent_environment(tmp_path, monkeypatch):
    template = tmp_path / "report.py"
    template.write_text("# %%\nprint('report')\n")
    output = tmp_path / "output"
    environment_before = dict(os.environ)
    calls = []

    def run(command, **kwargs):
        env = kwargs["env"]
        for key in ("MPLCONFIGDIR", "NUMBA_CACHE_DIR", "JUPYTER_DATA_DIR", "IPYTHONDIR",
                    "JUPYTER_CONFIG_DIR", "JUPYTER_RUNTIME_DIR"):
            assert Path(env[key]).is_dir()
            assert env[key].startswith("/tmp/pbmc-notebook-")
        assert env["REPORT_INPUT"] == "test-data"
        calls.append(command)
        if "--execute" in command:
            (output / "executed.ipynb").write_text(json.dumps({"cells": [], "metadata": {},
                                                            "nbformat": 4, "nbformat_minor": 5}))
        elif "html" in command:
            (output / "report.html").write_text('<nav><!-- REPORT_TOC_PLACEHOLDER --></nav><h2 id="results">Results</h2>')
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("pbmc_pipeline.notebook.subprocess.run", run)
    result = render_notebook(template_path=template, output_dir=output, project_root=tmp_path,
                             environment={"REPORT_INPUT": "test-data"}, title="Test report")
    assert len(calls) == 3
    assert "--ExecutePreprocessor.timeout=-1" in calls[1]
    assert '<a href="#results">Results</a>' in result.read_text()
    assert json.loads((output / "executed.ipynb").read_text())["metadata"]["title"] == "Test report"
    assert dict(os.environ) == environment_before


def test_execution_failure_preserves_the_previous_html_and_cleans_isolated_caches(tmp_path, monkeypatch):
    template = tmp_path / "report.py"
    template.write_text("# %%\nraise ValueError('bad report input')\n")
    output = tmp_path / "output"
    output.mkdir()
    (output / "report.html").write_text("previous successful report")
    calls = []
    caches = []

    def run(command, **kwargs):
        calls.append(command)
        caches.append(Path(kwargs["env"]["JUPYTER_RUNTIME_DIR"]))
        assert caches[-1].is_dir()
        failed = "--execute" in command
        return subprocess.CompletedProcess(command, 1 if failed else 0,
                                           stdout="kernel output\n", stderr="bad report input\n")

    monkeypatch.setattr("pbmc_pipeline.notebook.subprocess.run", run)
    with pytest.raises(subprocess.CalledProcessError) as error:
        render_notebook(template_path=template, output_dir=output, project_root=tmp_path,
                        environment={}, title="Report", log_path=output / "render.log")
    assert error.value.stderr == "bad report input\n"
    assert len(calls) == 2  # An execution failure must stop before HTML export.
    assert (output / "report.html").read_text() == "previous successful report"
    assert "bad report input" in (output / "render.log").read_text()
    assert all(not path.exists() for path in caches)
