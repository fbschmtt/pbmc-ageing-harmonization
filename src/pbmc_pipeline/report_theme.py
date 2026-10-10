"""Shared visual shell for generated analysis reports."""
from __future__ import annotations

from html import escape

REPORT_CSS = r"""
:root {
  --report-ink: #182230;
  --report-muted: #5f6b7a;
  --report-line: #d8e0e8;
  --report-paper: #ffffff;
  --report-page: #f5f7fa;
  --report-accent: #176b87;
  --report-accent-soft: #e6f3f7;
}
body.jp-Notebook {
  background: var(--report-page);
  color: var(--report-ink);
  font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  font-size: 17px;
}
main { max-width: 1600px; margin: 0 auto; padding: 2rem clamp(1.25rem, 3vw, 2.5rem) 4rem; }
.jp-Cell { margin: 1.4rem 0; }
.jp-Cell-outputWrapper, .jp-OutputArea { width: 100%; }
.jp-RenderedImage img { display: block; max-width: 100%; height: auto; margin: 0 auto; }
.jp-RenderedMarkdown h2 {
  margin-top: 3rem; padding-bottom: .55rem; border-bottom: 2px solid var(--report-line);
  color: var(--report-ink); font-size: clamp(1.45rem, 2.7vw, 2rem);
}
.jp-RenderedMarkdown h3 { color: var(--report-accent); margin-top: 2.25rem; }
.jp-RenderedMarkdown h4 { color: var(--report-accent); margin-top: 1.65rem; }
.jp-RenderedMarkdown p, .jp-RenderedMarkdown li { color: #354152; font-size: 1.04rem; line-height: 1.65; }
.jp-RenderedMarkdown blockquote {
  margin: 1rem 0; padding: .75rem 1rem; border-left: 4px solid #bd6a29;
  background: #fff6eb; color: #6e3e18;
}
.report-hero {
  margin: 0 0 2.5rem; padding: clamp(1.4rem, 4vw, 2.6rem); border-radius: 18px;
  background: linear-gradient(135deg, #315c6d, #527990); color: #fff;
  box-shadow: 0 12px 30px rgba(49, 92, 109, .16);
}
.report-eyebrow { margin: 0 0 .45rem; color: #bde7f0; font-size: .76rem; font-weight: 700; letter-spacing: .11em; text-transform: uppercase; }
.report-hero-title { margin: 0; color: #fff; font-size: clamp(2rem, 5vw, 3.4rem); line-height: 1.08; }
.report-subtitle { max-width: 48rem; margin: .85rem 0 1.4rem; color: #e2f3f7; line-height: 1.55; }
.report-metrics { display: flex; flex-wrap: wrap; align-items: flex-start; gap: .55rem; margin: 0; }
.report-metric { display: flex; align-items: baseline; gap: .55rem; padding: .62rem .78rem; border: 1px solid rgba(255,255,255,.22); border-radius: 10px; background: rgba(255,255,255,.1); }
.report-metric dt { color: #d6edf2; font-size: .9rem; font-weight: 600; }
.report-metric dd { margin: 0; color: #fff; font-size: .9rem; font-weight: 650; text-align: right; }
.report-toc-shell { position: fixed; z-index: 1000; top: 1rem; bottom: 1rem; left: 1rem; width: 260px; }
.report-toc-toggle { position: absolute; width: 1px; height: 1px; margin: -1px; overflow: hidden; clip: rect(0, 0, 0, 0); white-space: nowrap; clip-path: inset(50%); }
.report-toc-toggle-label { display: none; }
.report-toc { height: 100%; overflow-y: auto; padding: .9rem .8rem 1.1rem; border: 1px solid var(--report-line); border-radius: 14px; background: var(--report-paper); box-shadow: 0 10px 28px rgba(24, 34, 48, .12); }
.report-toc-heading { margin: 0 0 .65rem; padding: .25rem .45rem .65rem; border-bottom: 1px solid var(--report-line); color: var(--report-muted); font-size: .72rem; font-weight: 750; letter-spacing: .08em; text-transform: uppercase; }
.report-toc-list, .report-toc-list ul { margin: 0; padding: 0; list-style: none; }
.report-toc-list > li { margin: .12rem 0 .45rem; }
.report-toc-list ul { margin: .2rem 0 .35rem .65rem; padding-left: .65rem; border-left: 1px solid var(--report-line); }
.report-toc-list ul li { margin: .12rem 0; }
.report-toc a, .report-toc a:visited { display: block; padding: .25rem .45rem; border-radius: 6px; color: #344454; font-size: .88rem; font-weight: 600; line-height: 1.35; text-decoration: none; }
.report-toc-list ul a { color: var(--report-muted); font-size: .82rem; font-weight: 500; }
.report-toc a:hover, .report-toc a:focus-visible { background: var(--report-accent-soft); color: var(--report-accent); outline: none; }
.report-toc a:focus-visible { box-shadow: 0 0 0 2px var(--report-accent); }
html { scroll-behavior: smooth; scroll-padding-top: 1.25rem; }
.report-provenance { padding: .75rem 1rem; border: 1px solid var(--report-line); border-radius: 10px; background: var(--report-paper); }
.report-provenance summary { cursor: pointer; color: var(--report-accent); font-weight: 700; }
.report-provenance table { margin-top: .8rem; }
.report-details { margin: .75rem 0; padding: .65rem .85rem; border: 1px solid var(--report-line); border-radius: 9px; background: var(--report-paper); }
.report-details summary { cursor: pointer; color: var(--report-accent); font-weight: 650; }
.report-details table { width: 100%; margin-top: .7rem; border-collapse: collapse; }
.report-details th { background: #edf3f6; color: #263544; font-weight: 700; }
.report-details th, .report-details td { padding: .4rem .6rem; border: 1px solid var(--report-line); vertical-align: top; }
.report-details-figure img { display: block; max-width: 100%; height: auto; margin: .75rem auto .1rem; }
.report-evidence { margin: 2.5rem 0 1.25rem; padding: 1.25rem 1.4rem; border: 1px solid #c8dde5; border-radius: 14px; background: #f1f8fa; }
.report-evidence h3 { margin: 0 0 .45rem; color: var(--report-ink); font-size: 1.45rem; }
.report-evidence p { margin: 0 0 .85rem; color: #354152; line-height: 1.55; }
.report-evidence-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: .65rem; }
.report-evidence-item { padding: .7rem .8rem; border-radius: 9px; background: #fff; border: 1px solid #d7e7ec; }
.report-evidence-label { margin: 0; color: var(--report-muted); font-size: .77rem; font-weight: 700; text-transform: uppercase; letter-spacing: .04em; }
.report-evidence-value { margin: .22rem 0 0; color: var(--report-ink); font-weight: 650; line-height: 1.35; }
.jp-RenderedHTMLCommon { max-width: 100%; overflow-x: auto; }
.jp-RenderedHTMLCommon table { border-collapse: collapse; background: var(--report-paper); }
.jp-RenderedHTMLCommon th { background: #edf3f6; color: #263544; font-weight: 700; }
.jp-RenderedHTMLCommon th, .jp-RenderedHTMLCommon td { padding: .45rem .65rem; border: 1px solid var(--report-line); vertical-align: top; }
@media print {
  body.jp-Notebook { background: #fff; } main { max-width: none; padding: 0; }
  .report-hero { box-shadow: none; }
  .report-toc-shell { display: none; }
}
@media (min-width: 1200px) {
  main { width: calc(100% - 320px); max-width: 1500px; margin: 0 1.5rem 0 295px; }
}
@media (max-width: 1199px) {
  main { padding-top: 5rem; }
  .report-toc-shell { top: .5rem; right: .5rem; bottom: auto; left: .5rem; width: auto; }
  .report-toc-toggle-label { display: flex; align-items: center; justify-content: space-between; min-height: 2.8rem; padding: .6rem .9rem; border: 1px solid var(--report-line); border-radius: 10px; background: var(--report-paper); box-shadow: 0 6px 18px rgba(24, 34, 48, .12); color: var(--report-accent); cursor: pointer; font-weight: 700; }
  .report-toc-toggle-label::after { content: "＋"; font-size: 1.2rem; }
  .report-toc-toggle:checked + .report-toc-toggle-label::after { content: "−"; }
  .report-toc { display: none; height: auto; max-height: min(72vh, 680px); margin-top: .35rem; border-radius: 10px; }
  .report-toc-toggle:checked ~ .report-toc { display: block; }
  .report-toc-list { columns: 2; column-gap: 1rem; }
  .report-toc-list > li { break-inside: avoid; }
}
@media (max-width: 560px) {
  .report-toc-list { columns: 1; }
}
"""


def render_report_header(
    *, title: str, eyebrow: str, subtitle: str, metrics: list[tuple[str, str]],
    title_tag: str = "h1",
) -> str:
    """Render the shared report hero and responsive heading-based TOC shell."""
    if title_tag not in {"div", "h1"}:
        raise ValueError("title_tag must be 'div' or 'h1'")
    metrics_html = "".join(
        f'<div class="report-metric"><dt>{escape(label)}</dt>'
        f'<dd>{escape(value)}</dd></div>'
        for label, value in metrics
    )
    return f'''
<style>\n{REPORT_CSS}\n</style>
<header class="report-hero">
  <p class="report-eyebrow">{escape(eyebrow)}</p>
  <{title_tag} class="report-hero-title">{escape(title)}</{title_tag}>
  <p class="report-subtitle">{escape(subtitle)}</p>
  <dl class="report-metrics">{metrics_html}</dl>
  <div class="report-toc-shell">
    <input class="report-toc-toggle" type="checkbox" id="report-toc-toggle">
    <label class="report-toc-toggle-label" for="report-toc-toggle">Contents</label>
    <nav class="report-toc" aria-label="Report table of contents">
      <p class="report-toc-heading">On this page</p>
      <ul class="report-toc-list">
        <!-- REPORT_TOC_PLACEHOLDER -->
      </ul>
    </nav>
  </div>
</header>
'''
