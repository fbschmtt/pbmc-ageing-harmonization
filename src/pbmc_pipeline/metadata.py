"""Small shared helpers for adapters that map source metadata to the output contract."""
from __future__ import annotations

from typing import Any


class MetadataError(ValueError):
    """Raised when metadata cannot be transformed without ambiguity."""


def build_homogeneous_obs(source, study_id: str, rules: dict[str, dict]):
    """Apply an adapter-owned mapping from source columns to canonical columns."""
    import pandas as pd

    output = pd.DataFrame(index=source.index.copy())
    output["study"] = study_id
    for target, rule in rules.items():
        value = _evaluate_rule(source, rule)
        if rule.get("lower"):
            value = value.astype("string").str.lower()
        if "replace" in rule:
            for old, new in rule["replace"].items():
                value = value.astype("string").str.replace(old, new, regex=False)
        if "split" in rule:
            split = rule["split"]
            value = value.astype("string").str.split(split["separator"]).str[split["index"]]
        if rule.get("prefix"):
            value = rule["prefix"] + value.astype("string")
        if rule.get("dtype") == "float":
            value = pd.to_numeric(value, errors="raise").astype(float)
        output[target] = value
    return output


def _evaluate_rule(source, rule: dict[str, Any]):
    import pandas as pd

    if "constant" in rule:
        return pd.Series(rule["constant"], index=source.index)
    if rule.get("source_index"):
        return source.index.to_series(index=source.index)
    if "concat" in rule:
        parts = []
        for part in rule["concat"]:
            if isinstance(part, dict) and "constant" in part:
                parts.append(pd.Series(part["constant"], index=source.index, dtype="string"))
            else:
                _require_column(source, part)
                parts.append(source[part].astype("string"))
        result = parts[0]
        for part in parts[1:]:
            result = result + rule.get("separator", "") + part
        return result
    column = rule["source"]
    if column not in source and "fallback_constant" in rule:
        return pd.Series(rule["fallback_constant"], index=source.index)
    _require_column(source, column)
    result = source[column].copy()
    if "map" in rule:
        mapped = result.astype("string").map(rule["map"])
        result = mapped.fillna(rule.get("default")) if "default" in rule else mapped
    return result


def _require_column(frame, column: str) -> None:
    if column not in frame:
        raise MetadataError(f"Required source metadata column is absent: {column}")
