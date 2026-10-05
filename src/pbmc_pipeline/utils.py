"""Small shared helpers used across pipeline stages."""
from __future__ import annotations

import re


def slugify(value: str) -> str:
    """Convert a label to a stable lowercase, filename-friendly slug."""
    slug = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return slug or "unnamed"
