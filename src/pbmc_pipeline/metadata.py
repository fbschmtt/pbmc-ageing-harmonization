from __future__ import annotations

from pathlib import Path
from typing import Any


class MetadataError(ValueError):
    """Raised when metadata cannot be transformed without ambiguity."""


def _read_table(
    path: Path,
    file_format: str,
    columns: list[str] | None = None,
    read_options: dict[str, Any] | None = None,
):
    import pandas as pd

    options = read_options or {}
    if file_format == "csv":
        return pd.read_csv(path, usecols=columns, **options)
    if file_format == "tsv":
        return pd.read_csv(path, sep="\t", usecols=columns, **options)
    if file_format == "excel":
        return pd.read_excel(path, usecols=columns, **options)
    raise MetadataError(f"Unsupported metadata format: {file_format}")


def apply_joins(obs, joins: list[dict[str, Any]], root: Path, warnings: list[str]):
    """Attach external metadata while preserving exactly the original cell index."""
    result = obs.copy()
    original_index = result.index.copy()
    for spec in joins:
        path = root / spec["path"]
        if not path.exists():
            message = f"Auxiliary metadata missing: {path}"
            if spec.get("required", True):
                raise FileNotFoundError(message)
            warnings.append(message)
            continue
        table = _read_table(
            path, spec["format"], spec.get("columns"), spec.get("read_options")
        )
        table = table.rename(columns=spec.get("rename", {}))
        for column, rule in spec.get("derive_left", {}).items():
            value = _evaluate_rule(result, rule)
            for old, new in rule.get("replace", {}).items():
                value = value.astype("string").str.replace(old, new, regex=False)
            result[column] = value
        for column, mapping in spec.get("value_maps", {}).items():
            table[column] = table[column].replace(mapping)
        if duplicate_spec := spec.get("resolve_duplicates"):
            keys = duplicate_spec["keys"]
            duplicated = table.duplicated(keys, keep=False)
            for column, value in duplicate_spec.get("set", {}).items():
                table.loc[duplicated, column] = value
            table = table.drop_duplicates(keys, keep=duplicate_spec.get("keep", "first"))

        mode = spec.get("mode", "merge")
        if mode == "cell_lookup":
            index_column = spec["index_column"]
            if table[index_column].duplicated().any():
                raise MetadataError(f"Cell lookup {path} contains duplicate cell IDs")
            table = table.set_index(index_column)
            missing = original_index.difference(table.index)
            if len(missing):
                raise MetadataError(f"Cell lookup {path} is missing {len(missing)} cells")
            for column in table.columns:
                result[column] = table.loc[original_index, column].to_numpy()
        elif mode == "replace_obs_by_index":
            table = table.rename(columns={table.columns[0]: "cell_id"}).set_index("cell_id")
            result = table.loc[original_index].copy()
        elif mode == "terekhova_visits":
            table = _reshape_terekhova_visits(table)
            result = _safe_merge(
                result,
                table[["Tube_id", "Donor_id", "Gender", "Age", "BMI", "Ethnicity"]],
                ["Tube_id"],
            )
        elif mode == "merge":
            result = _safe_merge(result, table, spec["on"])
        else:
            raise MetadataError(f"Unsupported join mode: {mode}")

        if not result.index.equals(original_index):
            raise MetadataError(f"Join with {path} changed cell order or identity")
    return result


def _safe_merge(obs, table, on: list[str]):
    if table.duplicated(on).any():
        duplicated = int(table.duplicated(on, keep=False).sum())
        raise MetadataError(f"Auxiliary table has {duplicated} rows with duplicate keys {on}")
    index_name = obs.index.name or "cell_id"
    left = obs.reset_index(names=index_name)
    merged = left.merge(table, on=on, how="left", validate="many_to_one", suffixes=("", "_joined"))
    return merged.set_index(index_name)


def _reshape_terekhova_visits(table):
    import pandas as pd

    rename = {}
    for visit in (1, 2, 3):
        rename.update({
            f"Visit_{visit}_Tube_id": f"Tube_id_{visit}",
            f"Visit_{visit}_Age": f"Age_{visit}",
            f"Visit_{visit}_BMI": f"BMI_{visit}",
        })
    table = table.rename(columns=rename)
    return pd.wide_to_long(
        table,
        stubnames=["Tube_id", "Age", "BMI"],
        i=["Donor_id", "Gender", "Ethnicity", "Age_group"],
        j="Visit",
        sep="_",
        suffix=r"\d+",
    ).reset_index().dropna(subset=["Tube_id"])


def build_homogeneous_obs(source, study_id: str, rules: dict[str, dict]):
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


def _evaluate_rule(source, rule: dict):
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
