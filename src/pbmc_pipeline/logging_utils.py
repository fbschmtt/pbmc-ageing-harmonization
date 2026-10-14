"""Logging setup shared by the pipeline command-line entry points."""

from __future__ import annotations

import logging
import sys
import time


def configure_logging(level: str = "INFO", *, stream=None) -> None:
    """Write timestamped pipeline progress messages to the selected stream."""
    formatter = logging.Formatter(
        "%(asctime)sZ %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    formatter.converter = time.gmtime
    handler = logging.StreamHandler(sys.stderr if stream is None else stream)
    handler.setFormatter(formatter)

    logger = logging.getLogger("pbmc_pipeline")
    logger.handlers.clear()
    logger.addHandler(handler)
    logger.setLevel(level.upper())
    logger.propagate = False
