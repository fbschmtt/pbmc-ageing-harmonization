"""Execute and export the per-cell-type pseudobulk DE report."""
from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

import psutil

from .report_toc import populate_html_toc, set_notebook_title


def _run_report_command(command, *, cwd: Path, env: dict[str, str], log_path: Path) -> None:
    """Stream and timestamp report commands so long renders remain observable."""
    started_at = time.monotonic()
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with log_path.open("a", buffering=1) as log:
        log.write(f"[{timestamp}] START $ {shlex.join(command)}\n")
        process = subprocess.Popen(
            command, cwd=cwd, env=env, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, bufsize=1,
        )
        assert process.stdout is not None
        stop_monitor = threading.Event()
        log_lock = threading.Lock()

        def log_process_memory() -> None:
            root = psutil.Process(process.pid)
            while not stop_monitor.is_set():
                process_rows = []
                try:
                    processes = [root, *root.children(recursive=True)]
                except (psutil.Error, OSError):
                    processes = [root]
                for proc in processes:
                    try:
                        memory = proc.memory_info()
                        process_rows.append((memory.rss, proc.pid, proc.name()))
                    except (psutil.Error, OSError):
                        # Kernel children can exit between the process-tree and RSS reads.
                        continue
                rss_bytes = sum(row[0] for row in process_rows)
                largest = ",".join(
                    f"{name}:{pid}:{rss / 1024**2:.0f}MiB"
                    for rss, pid, name in sorted(process_rows, reverse=True)[:4]
                ) or "unavailable"
                timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
                with log_lock:
                    log.write(
                        f"[{timestamp}] RESOURCE process_tree_rss_mib={rss_bytes / 1024**2:.1f} "
                        f"process_count={len(process_rows)} largest_processes={largest}\n"
                    )
                if stop_monitor.wait(60):
                    return

        monitor = threading.Thread(target=log_process_memory, daemon=True)
        monitor.start()
        tail = ""
        for line in process.stdout:
            with log_lock:
                log.write(line)
            sys.stdout.write(line)
            sys.stdout.flush()
            tail = (tail + line)[-8192:]
        returncode = process.wait()
        stop_monitor.set()
        monitor.join()
        elapsed = time.monotonic() - started_at
        timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
        log.write(
            f"[{timestamp}] END exit_code={returncode} elapsed_seconds={elapsed:.1f}\n"
        )
    if returncode:
        raise subprocess.CalledProcessError(
            returncode, command, output=tail,
        )


def generate_differential_expression_report(
    *, result_dir: Path, output_dir: Path, template_path: Path, config_path: Path,
    project_root: Path,
) -> Path:
    """Execute a standalone Jupytext template and export its HTML report."""
    for path in (result_dir, template_path, config_path, project_root):
        if not path.exists():
            raise FileNotFoundError(path)
    result_dir = result_dir.resolve()
    output_dir = output_dir.resolve()
    template_path = template_path.resolve()
    config_path = config_path.resolve()
    project_root = project_root.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "age_trajectory_clusters.csv",
        "age_trajectory_cluster_means.csv",
        "residual_sample_clusters.csv",
        "residual_sample_pc_scores.csv",
        "residual_wgcna/leave_one_study_out_stability.csv",
    ):
        (output_dir / name).unlink(missing_ok=True)
    shutil.rmtree(output_dir / "go_enrichment", ignore_errors=True)
    render_log = output_dir / "report_render.log"
    render_log.unlink(missing_ok=True)
    progress_log = output_dir / "report_progress.log"
    progress_log.unlink(missing_ok=True)
    output_notebook = output_dir / "executed.ipynb"
    output_html = output_dir / "report.html"
    environment = os.environ.copy()
    environment.update({
        "DE_REPORT_DIR": str(result_dir),
        "DE_REPORT_OUTPUT_DIR": str(output_dir),
        "DE_REPORT_CONFIG": str(config_path),
        "DE_REPORT_PROGRESS_LOG": str(progress_log),
    })
    source_path = str(project_root / "src")
    environment["PYTHONPATH"] = os.pathsep.join(
        value for value in (source_path, environment.get("PYTHONPATH", "")) if value
    )

    with tempfile.TemporaryDirectory(prefix="pbmc-de-notebook-") as temporary_dir:
        temporary_path = Path(temporary_dir)
        matplotlib_cache = temporary_path / "matplotlib"
        jupyter_data = temporary_path / "jupyter-data"
        numba_cache = temporary_path / "numba"
        matplotlib_cache.mkdir()
        jupyter_data.mkdir()
        numba_cache.mkdir()
        environment["MPLCONFIGDIR"] = str(matplotlib_cache)
        environment["JUPYTER_DATA_DIR"] = str(jupyter_data)
        environment["NUMBA_CACHE_DIR"] = str(numba_cache)
        materialized_notebook = temporary_path / "differential_expression_report.ipynb"
        _run_report_command([
            sys.executable, "-m", "jupytext", "--to", "notebook",
            "--output", str(materialized_notebook), str(template_path),
        ], cwd=project_root, env=environment, log_path=render_log)
        _run_report_command([
            sys.executable, "-m", "jupyter", "nbconvert", "--execute", "--to", "notebook",
            "--ExecutePreprocessor.timeout=-1", f"--output={output_notebook.name}",
            f"--output-dir={output_dir}", str(materialized_notebook),
        ], cwd=project_root, env=environment, log_path=render_log)

    metadata_path = result_dir / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    cell_type = str(metadata.get("cell_type") or result_dir.parent.name)
    set_notebook_title(output_notebook, f"{cell_type} pseudobulk differential expression | PBMC ageing")
    _run_report_command([
        sys.executable, "-m", "jupyter", "nbconvert", "--to", "html",
        "--HTMLExporter.exclude_input=True",
        "--HTMLExporter.exclude_input_prompt=True",
        "--HTMLExporter.exclude_output_prompt=True",
        f"--output={output_html.name}", f"--output-dir={output_dir}",
        str(output_notebook),
    ], cwd=project_root, env=environment, log_path=render_log)
    populate_html_toc(output_html)
    return output_html


def render_per_type_reports(
    input_dir: Path, output_root: Path, *, template_path: Path, config_path: Path,
    project_root: Path,
) -> list[Path]:
    """Render a standalone DE report for every completed cell-type fit."""
    rendered = []
    for result_dir in sorted(input_dir.glob("*/differential_expression")):
        if not (result_dir / "run_metadata.json").is_file():
            continue
        output_dir = output_root / result_dir.parent.name / "differential_expression"
        rendered.append(generate_differential_expression_report(
            result_dir=result_dir,
            output_dir=output_dir,
            template_path=template_path,
            config_path=config_path,
            project_root=project_root,
        ))
    return rendered


def main() -> None:
    parser = argparse.ArgumentParser(description="Render per-cell-type pseudobulk DE reports")
    input_mode = parser.add_mutually_exclusive_group(required=True)
    input_mode.add_argument(
        "--input-dir", type=Path,
        help="render all completed per-type fits found under this directory",
    )
    input_mode.add_argument(
        "--result-dir", type=Path,
        help="render one completed per-type fit directory",
    )
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()
    if args.result_dir is not None:
        if args.output_dir is None:
            parser.error("--result-dir requires --output-dir")
        report = generate_differential_expression_report(
            result_dir=args.result_dir,
            output_dir=args.output_dir,
            template_path=args.template,
            config_path=args.config,
            project_root=args.project_root,
        )
        print(report)
        return
    if args.output_root is None:
        parser.error("--input-dir requires --output-root")
    reports = render_per_type_reports(
        args.input_dir,
        args.output_root,
        template_path=args.template,
        config_path=args.config,
        project_root=args.project_root,
    )
    print(f"Rendered {len(reports)} per-type DE reports")


if __name__ == "__main__":
    main()
