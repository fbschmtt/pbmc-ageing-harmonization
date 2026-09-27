"""Check local Markdown links and documented Make targets without extra tooling."""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCUMENTS = (
    ROOT / "AGENTS.md",
    ROOT / "README.md",
    ROOT / "IMPLEMENTATION.md",
    ROOT / "INPUT_FILES.md",
    ROOT / "PLAN.md",
    ROOT / "input_data/README.md",
    ROOT / "aifi_models/README",
)
MARKDOWN_LINK = re.compile(r"(?<!!)\[[^\]]*\]\(([^)]+)\)")
MAKE_COMMAND = re.compile(r"^\s*make\s+([A-Za-z0-9_-]+)\b", re.MULTILINE)
MAKE_TARGET = re.compile(r"^([A-Za-z0-9_-]+):", re.MULTILINE)


def local_link_errors(document: Path) -> list[str]:
    errors: list[str] = []
    for destination in MARKDOWN_LINK.findall(document.read_text(encoding="utf-8")):
        destination = destination.strip().strip("<>")
        path, _, _fragment = destination.partition("#")
        if not path or "://" in path or path.startswith(("mailto:", "/")):
            continue
        if not (document.parent / path).exists():
            errors.append(f"{document.relative_to(ROOT)}: missing link target {path}")
    return errors


def documented_make_target_errors() -> list[str]:
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    targets = set(MAKE_TARGET.findall(makefile))
    errors: list[str] = []
    for document in DOCUMENTS:
        for target in MAKE_COMMAND.findall(document.read_text(encoding="utf-8")):
            if target not in targets:
                errors.append(
                    f"{document.relative_to(ROOT)}: documented Make target does not exist: {target}"
                )
    return errors


def main() -> None:
    errors = [error for document in DOCUMENTS for error in local_link_errors(document)]
    errors.extend(documented_make_target_errors())
    if errors:
        print("Documentation checks failed:", *errors, sep="\n", file=sys.stderr)
        raise SystemExit(1)
    print(f"Documentation checks passed for {len(DOCUMENTS)} files.")


if __name__ == "__main__":
    main()
