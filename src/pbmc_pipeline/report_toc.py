"""Build a report navigation list from headings in rendered HTML."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from html import escape
from html.parser import HTMLParser
from pathlib import Path


class _HeadingParser(HTMLParser):
    """Collect rendered H2-H4 headings and their anchor IDs."""

    def __init__(self, *, minimum_level: int = 2) -> None:
        super().__init__(convert_charrefs=True)
        self.minimum_level = minimum_level
        self.headings: list[tuple[int, str, str]] = []
        self._current: tuple[int, str, list[str]] | None = None
        self._skip_anchor = False
        self._stack: list[tuple[str, str | None]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        if tag in {f"h{level}" for level in range(self.minimum_level, 5)}:
            inherited_id = next(
                (anchor for _, anchor in reversed(self._stack) if anchor), ""
            )
            self._current = (int(tag[1]), attributes.get("id") or inherited_id, [])
        elif self._current is not None and tag == "a":
            self._skip_anchor = "anchor-link" in (attributes.get("class") or "").split()
        if tag not in {
            "area", "base", "br", "col", "embed", "hr", "img", "input", "link",
            "meta", "param", "source", "track", "wbr",
        }:
            self._stack.append((tag, attributes.get("id")))

    def handle_data(self, data: str) -> None:
        if self._current is not None and not self._skip_anchor:
            self._current[2].append(data)

    def handle_endtag(self, tag: str) -> None:
        if self._current is not None and tag == f"h{self._current[0]}":
            level, anchor, pieces = self._current
            title = " ".join("".join(pieces).split()).replace("¶", "").strip()
            if title:
                if not anchor:
                    anchor = re.sub(r"[^a-z0-9]+", "-", title.casefold()).strip("-")
                self.headings.append((level, anchor, title))
            self._current = None
        elif tag == "a":
            self._skip_anchor = False
        for index in range(len(self._stack) - 1, -1, -1):
            if self._stack[index][0] == tag:
                del self._stack[index:]
                break


def set_notebook_title(notebook_path: Path, title: str) -> None:
    """Set the document title used by nbconvert for the browser tab."""
    import nbformat

    notebook = nbformat.read(notebook_path, as_version=4)
    notebook.metadata["title"] = title
    nbformat.write(notebook, notebook_path)


@dataclass
class _TocNode:
    level: int
    anchor: str
    title: str
    children: list[_TocNode] = field(default_factory=list)


def _render_toc_items(headings: list[tuple[int, str, str]]) -> str:
    sections: list[_TocNode] = []
    stack: list[_TocNode] = []
    for level, anchor, title in headings:
        node = _TocNode(level, anchor, title)
        while stack and stack[-1].level >= level:
            stack.pop()
        if stack:
            stack[-1].children.append(node)
        else:
            sections.append(node)
        stack.append(node)

    def render(node: _TocNode) -> str:
        link = f'<a href="#{escape(node.anchor, quote=True)}">{escape(node.title)}</a>'
        nested = ""
        if node.children:
            nested = "<ul>" + "".join(render(child) for child in node.children) + "</ul>"
        return f"<li>{link}{nested}</li>"

    return "".join(render(section) for section in sections)


def populate_html_toc(
    html_path: Path,
    *,
    placeholder: str = "<!-- REPORT_TOC_PLACEHOLDER -->",
    minimum_heading_level: int = 2,
) -> None:
    """Replace a report TOC placeholder with links derived from HTML headings.

    Place ``placeholder`` inside the report's TOC list. The final rendered HTML
    is parsed so conditional and programmatically emitted notebook headings are
    included without maintaining a second list of section names. By default,
    include H2-H4; reports with semantic H1 model headings can start at H1.
    """
    if minimum_heading_level not in {1, 2, 3, 4}:
        raise ValueError("minimum_heading_level must be between 1 and 4")
    document = html_path.read_text()
    if placeholder not in document:
        raise ValueError(f"Report TOC placeholder is missing from {html_path}")
    parser = _HeadingParser(minimum_level=minimum_heading_level)
    parser.feed(document)
    if not parser.headings:
        raise ValueError(f"No report headings were found in {html_path}")
    html_path.write_text(document.replace(placeholder, _render_toc_items(parser.headings), 1))
