import json
from pathlib import Path

from pbmc_pipeline.config import load_configuration

ROOT = Path(__file__).resolve().parents[1]


def test_configuration_is_complete():
    pipeline, studies, schema = load_configuration(ROOT, Path("config/pipeline.json"))
    assert set(studies["studies"]) == {"aida25", "aifi", "onek1k", "terekhova23", "wang25"}
    assert pipeline["schema_version"] == schema["schema_version"] == 1


def test_all_json_files_are_valid():
    for path in (ROOT / "config").glob("*.json"):
        with path.open(encoding="utf-8") as handle:
            assert isinstance(json.load(handle), dict)


def test_test_inputs_are_unique():
    _, document, _ = load_configuration(ROOT, Path("config/pipeline.json"))
    names = [study["test_input"] for study in document["studies"].values()]
    assert len(names) == len(set(names))
