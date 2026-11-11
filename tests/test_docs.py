from scripts.check_docs import (
    documentation_files,
    documented_make_target_errors,
    heading_anchors,
    local_link_errors,
)


def test_nested_guides_are_checked_and_links_resolve_from_their_new_directory(tmp_path):
    guide = tmp_path / "docs" / "analysis" / "models.md"
    guide.parent.mkdir(parents=True)
    readme = tmp_path / "README.md"
    readme.write_text("# Start here\n\n[Models](docs/analysis/models.md#age-models)\n")
    guide.write_text("# Age models\n\n[Start](../../README.md#start-here)\n")
    assert guide in documentation_files(tmp_path)
    assert not local_link_errors(readme, tmp_path)
    assert not local_link_errors(guide, tmp_path)
    guide.write_text("# Age models\n\n[Old guide](../../README.md#removed-section)\n")
    assert "missing heading" in local_link_errors(guide, tmp_path)[0]
    guide.write_text("# Age models\n\n[Missing](missing.md)\n")
    assert "missing link target" in local_link_errors(guide, tmp_path)[0]


def test_anchors_include_duplicate_headings_and_explicit_ids_but_ignore_code():
    contents = '# Sample `age`\n# Sample `age`\n<a id="stable-link"></a>\n```python\n# Example only\n```\n'
    assert heading_anchors(contents) == {"sample-age", "sample-age-1", "stable-link"}


def test_commands_check_every_target_and_discover_included_makefiles(tmp_path):
    (tmp_path / "make").mkdir()
    (tmp_path / "Makefile").write_text("include make/checks.mk\nhelp:\n")
    (tmp_path / "make/checks.mk").write_text("lint:\ntest-unit:\n")
    (tmp_path / "docs").mkdir()
    guide = tmp_path / "docs" / "testing.md"
    guide.write_text("```bash\nmake lint test-unit STUDIES=cima\n```\nUse `make help`.\n")
    assert not documented_make_target_errors(tmp_path)
    guide.write_text("```bash\nmake lint deleted-target\n```\n")
    assert documented_make_target_errors(tmp_path) == [
        "docs/testing.md: documented Make target does not exist: deleted-target"
    ]


def test_documented_make_targets_can_share_a_recipe(tmp_path):
    from scripts.check_docs import documented_make_target_errors

    (tmp_path / 'README.md').write_text('Run `make first-test second-test`.\n')
    (tmp_path / 'Makefile').write_text('first-test second-test:\n\t@true\n')
    assert documented_make_target_errors(tmp_path) == []
