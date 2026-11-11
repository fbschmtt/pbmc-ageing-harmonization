"""Shared execution lifecycle for Jupytext reports and their HTML audit artifacts."""
from __future__ import annotations

import os
import shlex
import subprocess
import sys
import tempfile
import threading
import time
from collections.abc import Callable, Mapping
from datetime import datetime, timezone
from pathlib import Path

import psutil

from .report_toc import populate_html_toc, set_notebook_title


def run_command(command, *, cwd: Path, env: dict[str, str], log_path: Path | None,
                progress_path: Path | None = None) -> None:
    """Capture routine output in a log and surface it if execution fails."""
    result = subprocess.run(
        command, cwd=cwd, env=env, capture_output=log_path is not None,
        text=True, check=False,
    )
    if log_path is not None:
        with log_path.open("a") as log:
            log.write(f"$ {shlex.join(command)}\n")
            log.write(result.stdout or "")
            log.write(result.stderr or "")
    if result.returncode:
        if log_path is not None:
            sys.stdout.write(result.stdout or "")
            sys.stderr.write(result.stderr or "")
        raise subprocess.CalledProcessError(
            result.returncode, command, output=result.stdout, stderr=result.stderr,
        )


def render_notebook(
    *, template_path: Path, output_dir: Path, project_root: Path,
    environment: Mapping[str, str], title: str | Callable[[], str],
    log_path: Path | None = None, stream: bool = False,
    progress_variable: str | None = None, minimum_heading_level: int = 2,
) -> Path:
    """Materialize, execute, title and export a report in an isolated environment.

    Report-specific wrappers own their input validation and domain preparation.
    All notebook tools share one temporary cache/configuration directory, kept
    alive until HTML export finishes. The caller's environment is never changed.
    """
    template_path, output_dir, project_root = (
        path.resolve() for path in (template_path, output_dir, project_root)
    )
    if not template_path.is_file():
        raise FileNotFoundError(template_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    if log_path is not None:
        log_path = log_path.resolve()
        log_path.unlink(missing_ok=True)
    executed, html = output_dir / "executed.ipynb", output_dir / "report.html"
    command_runner = stream_command if stream else run_command
    with tempfile.TemporaryDirectory(prefix="pbmc-notebook-") as temporary_dir:
        temporary = Path(temporary_dir)
        env = {**os.environ, **environment}
        source_path = str(project_root / "src")
        env["PYTHONPATH"] = os.pathsep.join(
            value for value in (source_path, env.get("PYTHONPATH", "")) if value
        )
        for key, directory in (
            ("MPLCONFIGDIR", "matplotlib"), ("NUMBA_CACHE_DIR", "numba"),
            ("JUPYTER_DATA_DIR", "jupyter-data"), ("IPYTHONDIR", "ipython"),
            ("JUPYTER_CONFIG_DIR", "jupyter-config"),
            ("JUPYTER_RUNTIME_DIR", "jupyter-runtime"),
        ):
            path = temporary / directory
            path.mkdir()
            env[key] = str(path)
        progress_path = temporary / "progress.log" if progress_variable else None
        if progress_variable:
            env[progress_variable] = str(progress_path)
        notebook = temporary / f"{template_path.stem}.ipynb"
        command_runner([
            sys.executable, "-m", "jupytext", "--to", "notebook",
            "--output", str(notebook), str(template_path),
        ], cwd=project_root, env=env, log_path=log_path)
        command_runner([
            sys.executable, "-m", "jupyter", "nbconvert", "--execute", "--to", "notebook",
            "--ExecutePreprocessor.timeout=-1", "--KernelManager.transport_encryption=required",
            f"--output={executed.name}", f"--output-dir={output_dir}", str(notebook),
        ], cwd=project_root, env=env, log_path=log_path, progress_path=progress_path)
        set_notebook_title(executed, title() if callable(title) else title)
        command_runner([
            sys.executable, "-m", "jupyter", "nbconvert", "--to", "html",
            "--HTMLExporter.exclude_input=True", "--HTMLExporter.exclude_input_prompt=True",
            "--HTMLExporter.exclude_output_prompt=True", f"--output={html.name}",
            f"--output-dir={output_dir}", str(executed),
        ], cwd=project_root, env=env, log_path=log_path)
    populate_html_toc(html, minimum_heading_level=minimum_heading_level)
    return html


def stream_command(
    command, *, cwd: Path, env: dict[str, str], log_path: Path,
    progress_path: Path | None = None,
) -> None:
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

        def tail_progress() -> None:
            if progress_path is None:
                return
            offset = 0
            while True:
                try:
                    with progress_path.open(encoding="utf-8") as progress:
                        progress.seek(offset)
                        lines = progress.readlines()
                        offset = progress.tell()
                except FileNotFoundError:
                    lines = []
                if lines:
                    text = "".join(lines)
                    with log_lock:
                        log.write(text)
                    sys.stdout.write(text)
                    sys.stdout.flush()
                if stop_monitor.is_set():
                    return
                stop_monitor.wait(0.5)

        monitor = threading.Thread(target=log_process_memory, daemon=True)
        monitor.start()
        progress_monitor = threading.Thread(target=tail_progress, daemon=True)
        progress_monitor.start()
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
        progress_monitor.join()
        elapsed = time.monotonic() - started_at
        timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
        log.write(
            f"[{timestamp}] END exit_code={returncode} elapsed_seconds={elapsed:.1f}\n"
        )
    if returncode:
        raise subprocess.CalledProcessError(
            returncode, command, output=tail,
        )
