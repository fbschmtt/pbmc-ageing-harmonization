from datetime import datetime
from html.parser import HTMLParser

from pbmc_pipeline.report_theme import render_report_header


class _GeneratedTimeParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.values = []

    def handle_starttag(self, tag, attrs):
        if tag == "time":
            self.values.append(dict(attrs).get("datetime"))


def test_shared_report_header_includes_utc_generation_timestamp():
    header = render_report_header(
        title="Example report",
        eyebrow="Testing",
        subtitle="Report subtitle",
        metrics=[],
    )
    parser = _GeneratedTimeParser()
    parser.feed(header)

    assert len(parser.values) == 1
    generated_at = datetime.fromisoformat(parser.values[0])
    assert generated_at.tzinfo is not None
    assert generated_at.utcoffset().total_seconds() == 0
    assert " UTC</time>" in header
    assert "Generated " in header
    assert 'class="report-generated-at"' in header
