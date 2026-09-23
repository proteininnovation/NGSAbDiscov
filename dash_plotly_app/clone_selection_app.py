"""Dash app for interactive clone selection from NGSAbDiscov results."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import threading
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_RESULTS_ROOT = "/alphafold/combio/NGS"
APP_CACHE_VERSION = "20260729_table_v1"
APP_CACHE_DIR = Path(os.environ.get("NGSABDISCOV_APP_CACHE", "/tmp/ngsabdiscov_clone_app_cache"))
_CACHE_WARM_LOCK = threading.Lock()
_CACHE_WARMING_FOLDERS: set[str] = set()


DEFAULT_DISPLAY_COLUMNS = [
    "CDR1",
    "CDR2",
    "CDR3",
    "cdr3_aa",
    "vh_scaffold",
    "vl_scaffold",
    "rank",
    "max_freq",
    "log2_fold_change",
    "condition_enrichment_score",
    "condition_enrichment_pass",
    "max_negative_control_freq",
    "negative_control_condition",
    "selected_over_F_N",
    "F_N_over_selected",
    "F_N_ratio_pass",
    "F_N_PSR_like",
    "log2_4nM_vs_20nM",
    "log2_4nM_vs_100nM",
    "log2_20nM_vs_100nM",
    "cdr3_cluster",
    "cdr3_cluster_size",
    "diversity_representative",
    "mean_delphi_score",
    "mean_psr_score",
    "mean_sec_score",
    "ML_SEQUENCE_OK",
    "HSEQ_len",
    "LSEQ_len",
]

AFFINITY_CONDITION_PAIRS = [
    ("4nM", "20nM"),
    ("4nM", "100nM"),
    ("20nM", "100nM"),
]

CONDITION_RULE_OPTIONS = [
    {"label": "Off", "value": "off"},
    {"label": "Affinity ladder: all available", "value": "affinity_all"},
    {"label": "Affinity ladder: any available", "value": "affinity_any"},
    {"label": "4nM > 20nM", "value": "4nM_vs_20nM"},
    {"label": "4nM > 100nM", "value": "4nM_vs_100nM"},
    {"label": "20nM > 100nM", "value": "20nM_vs_100nM"},
]

FOLD_CHANGE_COMPARISON_OPTIONS = [
    {"label": "4nM / 20nM", "value": "4nM_vs_20nM"},
    {"label": "4nM / 100nM", "value": "4nM_vs_100nM"},
    {"label": "20nM / 100nM", "value": "20nM_vs_100nM"},
    {"label": "Manual sample columns", "value": "manual"},
]

NEGATIVE_CONTROL_OPTIONS = [
    {"label": "Off", "value": "off"},
    {"label": "Auto selected / F_N", "value": "auto"},
    {"label": "4nM / F_N", "value": "4nM"},
    {"label": "20nM / F_N", "value": "20nM"},
    {"label": "100nM / F_N", "value": "100nM"},
]


APP_CSS = """
html, body { height: 100%; }
body { margin: 0; background: #f1f4f8; color: #18212f; font-family: Arial, Helvetica, sans-serif; overflow: auto; }
.header { padding: 10px 16px; background: #ffffff; border-bottom: 1px solid #d8dee8; }
h1 { margin: 0 0 4px 0; font-size: 22px; }
.path { color: #667085; font-size: 12px; }
.cards { display: flex; gap: 8px; flex-wrap: wrap; margin-top: 8px; }
.card { border: 1px solid #d8dee8; background: #f8fafc; border-radius: 6px; padding: 6px 9px; min-width: 112px; }
.card strong { display: block; font-size: 16px; color: #176d7a; }
.card span { font-size: 12px; }
.workspace { display: grid; grid-template-columns: 300px minmax(0, 1fr); gap: 10px; padding: 10px; min-height: calc(100vh - 108px); box-sizing: border-box; }
.sidebar { background: #ffffff; border: 1px solid #d8dee8; border-radius: 6px; padding: 10px; align-self: start; max-height: 100%; overflow-y: auto; box-sizing: border-box; }
.sidebar label { display: block; margin: 7px 0 3px; font-weight: 700; font-size: 11px; color: #344054; }
.sidebar input { width: 100%; box-sizing: border-box; padding: 6px; border: 1px solid #ccd3dd; border-radius: 4px; }
.sidebar .action-button { width: 100%; margin-top: 9px; padding: 8px; border: 1px solid #176d7a; border-radius: 4px; background: #176d7a; color: white; font-weight: 700; cursor: pointer; }
.sidebar .secondary-button { background: #ffffff; color: #176d7a; }
.hint { margin-top: 4px; color: #667085; font-size: 11px; line-height: 1.35; }
.main { display: grid; grid-template-rows: auto auto minmax(250px, 1fr); gap: 10px; min-width: 0; min-height: 0; }
.dash-graph, .table-panel, .plot-panel { background: #ffffff; border: 1px solid #d8dee8; border-radius: 6px; padding: 6px; min-height: 0; box-sizing: border-box; }
.selection-toolbar { display: flex; justify-content: space-between; align-items: center; min-height: 34px; gap: 10px; }
.selection-status { color: #667085; font-size: 12px; }
.selection-button { padding: 7px 10px; border: 1px solid #176d7a; border-radius: 4px; background: #ffffff; color: #176d7a; font-weight: 700; cursor: pointer; }
.table-panel { max-height: 340px; overflow: hidden; }
.plot-panel { overflow: hidden; }
.plot-tabs { height: 100%; display: flex; flex-direction: column; }
.plot-tabs .tab-content { flex: 1 1 auto; min-height: 0; overflow: hidden; }
.plot-graph { height: calc(100vh - 475px); min-height: 245px; max-height: 340px; }
@media (max-height: 780px) {
  .workspace { min-height: calc(100vh - 96px); }
  .main { grid-template-rows: auto auto minmax(225px, 1fr); }
  .table-panel { max-height: 260px; }
  .plot-graph { height: calc(100vh - 430px); min-height: 220px; max-height: 300px; }
}
@media (max-width: 1000px) {
  .workspace { grid-template-columns: 1fr; height: auto; overflow: visible; }
  .sidebar { max-height: none; }
  .main { display: grid; grid-template-rows: auto auto auto; overflow: visible; }
  .table-panel { max-height: 320px; }
  .plot-graph { height: 320px; max-height: none; }
}
"""


def _clean_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none"} else text


def _score_columns(df: pd.DataFrame, label: str | None = None) -> list[str]:
    cols = []
    for col in df.columns:
        name = str(col)
        lowered = name.lower()
        if not lowered.endswith("_score"):
            continue
        if lowered in {"mean_delphi_score", "mean_psr_score", "mean_sec_score"}:
            continue
        if label and not lowered.startswith(f"{label.lower()}_"):
            continue
        if any(token in lowered for token in ["xgboost", "transformer", "rf", "delph"]):
            cols.append(name)
    return cols


def _freq_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in df.columns if str(col).startswith("freq ")]


def _count_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in df.columns if str(col).startswith("count ")]


def _ensure_numeric(df: pd.DataFrame, columns: list[str]) -> None:
    for col in columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")


def _to_optional_float(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        value = value.strip()
        if not value:
            return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _to_float(value: Any, default: float = 0.0) -> float:
    parsed = _to_optional_float(value)
    return default if parsed is None else parsed


def _condition_slug(condition: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", condition).strip("_")


def _condition_from_freq_col(freq_col: str) -> str | None:
    label = str(freq_col).removeprefix("freq ").strip()
    lowered = label.lower()
    concentration_match = re.search(r"(?<![A-Za-z0-9])(\d+(?:\.\d+)?)(?:\s*)(nm|um|µm|μm)(?![A-Za-z0-9])", label, re.IGNORECASE)
    is_negative = "__f_n__" in lowered or "negative" in lowered
    if concentration_match:
        amount = concentration_match.group(1)
        if "." in amount:
            amount = amount.rstrip("0").rstrip(".")
        unit = concentration_match.group(2).replace("µ", "u").replace("μ", "u")
        unit = unit[0].lower() + unit[1].upper()
        condition = f"{amount}{unit}"
        if is_negative or condition.lower() == "2um":
            return "negative_control"
        return condition
    if is_negative:
        return "negative_control"
    return None


def _condition_frequency_map(df: pd.DataFrame) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    for col in _freq_columns(df):
        condition = _condition_from_freq_col(str(col))
        if condition:
            mapping.setdefault(condition, []).append(col)
    return mapping


def _negative_control_freq_columns(df: pd.DataFrame) -> list[str]:
    cols = []
    for col in _freq_columns(df):
        text = str(col).lower()
        condition = _condition_from_freq_col(str(col))
        if condition == "negative_control" or "__f_n__" in text or "negative" in text or "2um" in text or "2um" in text.replace("µ", "u").replace("μ", "u"):
            cols.append(col)
    return cols


def _with_negative_control_summary(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    cols = _negative_control_freq_columns(out)
    if cols:
        out["max_negative_control_freq"] = out[cols].apply(pd.to_numeric, errors="coerce").fillna(0).max(axis=1)
    else:
        out["max_negative_control_freq"] = 0.0
    return out


def _condition_freq_col(condition: str) -> str:
    return f"freq_condition_{_condition_slug(condition)}"


def _condition_fc_col(numerator_condition: str, denominator_condition: str) -> str:
    return f"log2_{_condition_slug(numerator_condition)}_vs_{_condition_slug(denominator_condition)}"


def _condition_pair_from_value(value: str | None) -> tuple[str, str] | None:
    if not value or value == "manual" or "_vs_" not in str(value):
        return None
    numerator, denominator = str(value).split("_vs_", 1)
    return numerator, denominator


def _condition_rule_pairs(rule: str) -> list[tuple[str, str]]:
    if rule in {"affinity_all", "affinity_any"}:
        return AFFINITY_CONDITION_PAIRS
    if "_vs_" in str(rule):
        numerator, denominator = str(rule).split("_vs_", 1)
        return [(numerator, denominator)]
    return []


def _available_fold_change_options(df: pd.DataFrame) -> list[dict[str, str]]:
    mapping = _condition_frequency_map(df)
    options = []
    for option in FOLD_CHANGE_COMPARISON_OPTIONS:
        pair = _condition_pair_from_value(option["value"])
        if pair is None or (pair[0] in mapping and pair[1] in mapping):
            options.append(option)
    return options or [{"label": "Manual sample columns", "value": "manual"}]


def _best_freq_column_for_condition(df: pd.DataFrame, condition: str) -> str | None:
    cols = _condition_frequency_map(df).get(condition, [])
    if not cols:
        return None
    return sorted(cols, key=_round_sort_key)[-1]


def _default_fold_change_columns(df: pd.DataFrame) -> tuple[str | None, str | None]:
    for numerator_condition, denominator_condition in AFFINITY_CONDITION_PAIRS:
        numerator = _best_freq_column_for_condition(df, numerator_condition)
        denominator = _best_freq_column_for_condition(df, denominator_condition)
        if numerator and denominator and numerator != denominator:
            return numerator, denominator

    freq_cols = sorted(_freq_columns(df), key=_round_sort_key)
    if len(freq_cols) >= 2:
        return freq_cols[-1], freq_cols[-2]
    if len(freq_cols) == 1:
        return freq_cols[0], None
    return None, None


def _condition_available(df: pd.DataFrame, condition: str) -> bool:
    return _condition_freq_col(condition) in df.columns


def _selected_negative_control_condition(
    df: pd.DataFrame,
    negative_control_mode: str | None,
    fold_change_comparison: str | None,
    numerator: str | None,
) -> str | None:
    mode = negative_control_mode or "off"
    if mode == "off":
        return None
    if mode != "auto":
        return mode if _condition_available(df, mode) else None

    pair = _condition_pair_from_value(fold_change_comparison)
    if pair and _condition_available(df, pair[0]):
        return pair[0]

    if fold_change_comparison == "manual" and numerator:
        condition = _condition_from_freq_col(numerator)
        if condition and condition != "negative_control" and _condition_available(df, condition):
            return condition

    for condition in ["4nM", "20nM", "100nM"]:
        if _condition_available(df, condition):
            return condition
    return None


def _ratio_with_zero_handling(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    numerator = pd.to_numeric(numerator, errors="coerce").fillna(0.0)
    denominator = pd.to_numeric(denominator, errors="coerce").fillna(0.0)
    ratio = pd.Series(np.nan, index=numerator.index, dtype=float)
    ratio.loc[(numerator > 0) & (denominator == 0)] = np.inf
    ratio.loc[(numerator == 0) & (denominator == 0)] = np.nan
    mask = denominator > 0
    ratio.loc[mask] = numerator.loc[mask] / denominator.loc[mask]
    return ratio


def _with_negative_control_ratio(
    df: pd.DataFrame,
    negative_control_mode: str | None,
    ratio_cutoff: float | None,
    psr_cutoff: float | None,
    fold_change_comparison: str | None,
    numerator: str | None,
) -> pd.DataFrame:
    out = df.copy()
    out["negative_control_condition"] = ""
    out["selected_over_F_N"] = np.nan
    out["F_N_over_selected"] = np.nan
    out["F_N_ratio_pass"] = True
    out["F_N_PSR_like"] = False

    if (negative_control_mode or "off") == "off":
        return out

    negative_col = _condition_freq_col("negative_control")
    condition = _selected_negative_control_condition(out, negative_control_mode, fold_change_comparison, numerator)
    if not condition or negative_col not in out.columns:
        return out

    positive_col = _condition_freq_col(condition)
    if positive_col not in out.columns:
        return out

    positive_freq = pd.to_numeric(out[positive_col], errors="coerce").fillna(0.0)
    negative_freq = pd.to_numeric(out[negative_col], errors="coerce").fillna(0.0)
    selected_over_negative = _ratio_with_zero_handling(positive_freq, negative_freq)
    negative_over_selected = _ratio_with_zero_handling(negative_freq, positive_freq)

    min_ratio = _to_float(ratio_cutoff, 0.3)
    psr_ratio = _to_float(psr_cutoff, 3.0)
    out["negative_control_condition"] = condition
    out["selected_over_F_N"] = selected_over_negative
    out["F_N_over_selected"] = negative_over_selected
    out["F_N_PSR_like"] = negative_over_selected > psr_ratio
    out["F_N_ratio_pass"] = (selected_over_negative >= min_ratio) & ~out["F_N_PSR_like"]
    return out


def _with_condition_enrichment(
    df: pd.DataFrame,
    condition_rule: str = "off",
    condition_cutoff: float | None = 0.0,
) -> pd.DataFrame:
    out = df.copy()
    mapping = _condition_frequency_map(out)
    for condition, cols in mapping.items():
        numeric = out[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        out[_condition_freq_col(condition)] = numeric.max(axis=1)

    pseudo = 1e-9
    pair_cols: list[str] = []
    for numerator_condition, denominator_condition in AFFINITY_CONDITION_PAIRS:
        numerator_col = _condition_freq_col(numerator_condition)
        denominator_col = _condition_freq_col(denominator_condition)
        fc_col = _condition_fc_col(numerator_condition, denominator_condition)
        if numerator_col in out.columns and denominator_col in out.columns:
            numerator_values = pd.to_numeric(out[numerator_col], errors="coerce").fillna(0)
            denominator_values = pd.to_numeric(out[denominator_col], errors="coerce").fillna(0)
            out[fc_col] = np.log2((numerator_values + pseudo) / (denominator_values + pseudo))
            pair_cols.append(fc_col)
        elif fc_col not in out.columns:
            out[fc_col] = np.nan

    out["condition_enrichment_score"] = np.nan
    out["condition_enrichment_pass"] = False
    rule = condition_rule or "off"
    cutoff = _to_float(condition_cutoff, 0.0)
    selected_pair_cols = [
        _condition_fc_col(numerator, denominator)
        for numerator, denominator in _condition_rule_pairs(rule)
        if _condition_fc_col(numerator, denominator) in out.columns
        and pd.to_numeric(out[_condition_fc_col(numerator, denominator)], errors="coerce").notna().any()
    ]

    if selected_pair_cols:
        scores = out[selected_pair_cols].apply(pd.to_numeric, errors="coerce")
        if rule == "affinity_any":
            out["condition_enrichment_score"] = scores.max(axis=1)
            out["condition_enrichment_pass"] = out["condition_enrichment_score"] > cutoff
        else:
            out["condition_enrichment_score"] = scores.min(axis=1)
            out["condition_enrichment_pass"] = scores.notna().all(axis=1) & (out["condition_enrichment_score"] > cutoff)
    return out


def _with_selected_log2_fold_change(
    df: pd.DataFrame,
    fold_change_comparison: str | None,
    numerator: str | None,
    denominator: str | None,
) -> pd.DataFrame:
    if fold_change_comparison == "manual":
        return _with_log2_fold_change(df, numerator, denominator)

    pair = _condition_pair_from_value(fold_change_comparison)
    if pair:
        fc_col = _condition_fc_col(pair[0], pair[1])
        if fc_col in df.columns and pd.to_numeric(df[fc_col], errors="coerce").notna().any():
            out = df.copy()
            out["log2_fold_change"] = pd.to_numeric(out[fc_col], errors="coerce")
            return out
    return _with_log2_fold_change(df, numerator, denominator)


def _fold_change_axis_label(
    df: pd.DataFrame,
    fold_change_comparison: str | None,
    numerator: str | None,
    denominator: str | None,
) -> str | None:
    if fold_change_comparison == "manual" and numerator in df.columns and denominator in df.columns and numerator != denominator:
        return f"log2({_short_sample_label(numerator)} / {_short_sample_label(denominator)})"

    pair = _condition_pair_from_value(fold_change_comparison)
    if pair:
        fc_col = _condition_fc_col(pair[0], pair[1])
        if fc_col in df.columns and pd.to_numeric(df[fc_col], errors="coerce").notna().any():
            return f"log2({pair[0]} / {pair[1]})"
    if numerator in df.columns and denominator in df.columns and numerator != denominator:
        return f"log2({_short_sample_label(numerator)} / {_short_sample_label(denominator)})"
    return None


def _mean_score(df: pd.DataFrame, label: str) -> pd.Series:
    mean_col = f"mean_{label}_score"
    if mean_col in df.columns:
        return pd.to_numeric(df[mean_col], errors="coerce")
    cols = _score_columns(df, label)
    if cols:
        return df[cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    return pd.Series(np.nan, index=df.index)


def _cdr3_value(row: pd.Series) -> str:
    for col in ["CDR3", "cdr3_aa"]:
        if col in row.index:
            value = _clean_text(row[col]).upper()
            if value:
                return value
    return ""


def _hamming_distance(left: str, right: str) -> int:
    if len(left) != len(right):
        return max(len(left), len(right))
    return sum(aa != bb for aa, bb in zip(left, right))


def _with_cdr3_clusters(df: pd.DataFrame, max_mismatches: int = 2, max_rows: int = 3000) -> pd.DataFrame:
    out = df.copy()
    out["cdr3_cluster"] = ""
    out["cdr3_cluster_size"] = 0
    out["diversity_representative"] = False
    if out.empty:
        return out

    cluster_centroids: dict[int, list[tuple[str, str]]] = {}
    row_to_cluster: dict[Any, str] = {}
    cluster_members: dict[str, list[Any]] = {}
    cluster_count = 0

    for idx, row in out.head(max_rows).iterrows():
        seq = _cdr3_value(row)
        if not seq:
            continue
        length_bucket = len(seq)
        assigned_cluster = None
        for cluster_id, centroid in cluster_centroids.get(length_bucket, []):
            if _hamming_distance(seq, centroid) <= max_mismatches:
                assigned_cluster = cluster_id
                break
        if assigned_cluster is None:
            cluster_count += 1
            assigned_cluster = f"C{cluster_count:04d}"
            cluster_centroids.setdefault(length_bucket, []).append((assigned_cluster, seq))
            out.at[idx, "diversity_representative"] = True
        row_to_cluster[idx] = assigned_cluster
        cluster_members.setdefault(assigned_cluster, []).append(idx)

    for idx, cluster_id in row_to_cluster.items():
        out.at[idx, "cdr3_cluster"] = cluster_id
        out.at[idx, "cdr3_cluster_size"] = len(cluster_members[cluster_id])
    return out


def _sort_leads(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    ranked = df.copy()
    ranked["_sort_freq"] = pd.to_numeric(ranked.get("max_freq", pd.Series(0, index=ranked.index)), errors="coerce").fillna(0)
    ranked["_sort_rank"] = pd.to_numeric(ranked.get("rank", pd.Series(np.inf, index=ranked.index)), errors="coerce").fillna(np.inf)
    ranked["_sort_score"] = pd.to_numeric(ranked.get("mean_delphi_score", pd.Series(-1, index=ranked.index)), errors="coerce").fillna(-1)
    ranked = ranked.sort_values(["_sort_freq", "_sort_rank", "_sort_score"], ascending=[False, True, False])
    return ranked.drop(columns=["_sort_freq", "_sort_rank", "_sort_score"], errors="ignore")


def _normalize_final_table(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path)
    target = path.stem.replace("_final_leads", "")
    if "target" in df.columns:
        df["target"] = df["target"].apply(_clean_text)
        df.loc[df["target"] == "", "target"] = target
    else:
        df.insert(0, "target", target)
    df["_source_file"] = path.name

    if "CDR3" not in df.columns and "cdr3_aa" in df.columns:
        df["CDR3"] = df["cdr3_aa"]
    if "cdr3_aa" not in df.columns and "CDR3" in df.columns:
        df["cdr3_aa"] = df["CDR3"]

    freq_cols = _freq_columns(df)
    count_cols = _count_columns(df)
    _ensure_numeric(df, [*freq_cols, *count_cols, "rank", "max_freq"])

    if "max_freq" not in df.columns and freq_cols:
        df["max_freq"] = df[freq_cols].max(axis=1)
    elif "max_freq" not in df.columns:
        df["max_freq"] = np.nan

    df["mean_psr_score"] = _mean_score(df, "psr")
    df["mean_sec_score"] = _mean_score(df, "sec")
    if "mean_delphi_score" not in df.columns:
        df["mean_delphi_score"] = df[["mean_psr_score", "mean_sec_score"]].mean(axis=1)
    else:
        df["mean_delphi_score"] = pd.to_numeric(df["mean_delphi_score"], errors="coerce")

    for col in ["HSEQ", "LSEQ", "CDR1", "CDR2", "CDR3", "cdr3_aa", "vh_scaffold", "vl_scaffold", "library"]:
        if col in df.columns:
            df[col] = df[col].apply(_clean_text)
    if "HSEQ" in df.columns:
        df["HSEQ_len"] = df["HSEQ"].astype(str).str.len()
    if "LSEQ" in df.columns:
        df["LSEQ_len"] = df["LSEQ"].astype(str).str.len()

    if "ML_SEQUENCE_OK" not in df.columns:
        has_hseq = df["HSEQ"].astype(str).str.len() > 0 if "HSEQ" in df.columns else False
        if "LSEQ" in df.columns:
            if "library" in df.columns:
                library_text = df["library"].astype(str).str.lower()
            else:
                library_text = pd.Series("", index=df.index)
            has_lseq_or_vhh = (df["LSEQ"].astype(str).str.len() > 0) | library_text.str.contains("vhh")
        else:
            has_lseq_or_vhh = True
        df["ML_SEQUENCE_OK"] = has_hseq & has_lseq_or_vhh

    return df


def scan_results(results_folder: str | Path) -> pd.DataFrame:
    folder = Path(results_folder).expanduser().resolve()
    lead_dir = folder / "by_protein"
    files = sorted(lead_dir.glob("*_final_leads.xlsx"))
    if not files:
        raise FileNotFoundError(f"No by_protein/*_final_leads.xlsx files found in {folder}")
    return pd.DataFrame(
        {
            "target": [path.stem.replace("_final_leads", "") for path in files],
            "path": [str(path) for path in files],
        }
    )


def _result_label(results_folder: str | Path) -> str:
    folder = Path(results_folder)
    run = folder.parent.name
    return f"{run} / {folder.name}"


def scan_project_results(results_root: str | Path = DEFAULT_RESULTS_ROOT) -> pd.DataFrame:
    root = Path(results_root).expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(f"Results root does not exist: {root}")

    records = []
    seen = set()
    for by_protein_dir in root.glob("*/*/by_protein"):
        if not by_protein_dir.is_dir() or by_protein_dir in seen:
            continue
        if not any(by_protein_dir.glob("*_final_leads.xlsx")):
            continue
        seen.add(by_protein_dir)
        results_folder = by_protein_dir.parent
        records.append(
            {
                "label": _result_label(results_folder),
                "results_folder": str(results_folder),
                "run_folder": str(results_folder.parent),
                "results_name": results_folder.name,
            }
        )
    if not records:
        raise FileNotFoundError(f"No NGSAbDiscov result folders found under {root}")
    return pd.DataFrame(records).sort_values(["label"]).reset_index(drop=True)


def _table_cache_path(path: Path) -> Path | None:
    try:
        stat = path.stat()
        cache_text = f"{APP_CACHE_VERSION}|{path.resolve()}|{stat.st_size}|{stat.st_mtime_ns}"
    except OSError:
        return None
    digest = hashlib.sha1(cache_text.encode("utf-8")).hexdigest()
    return APP_CACHE_DIR / f"{digest}.pkl"


def _load_final_table_with_disk_cache(path: Path) -> pd.DataFrame:
    cache_path = _table_cache_path(path)
    if cache_path and cache_path.exists():
        try:
            return pd.read_pickle(cache_path)
        except Exception:
            pass

    df = _normalize_final_table(path)
    if cache_path:
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            tmp_path = cache_path.with_name(f"{cache_path.name}.{os.getpid()}.tmp")
            df.to_pickle(tmp_path)
            os.replace(tmp_path, cache_path)
        except Exception:
            pass
    return df


@lru_cache(maxsize=64)
def _load_final_table_cached(path: str) -> pd.DataFrame:
    return _load_final_table_with_disk_cache(Path(path)).copy()


@lru_cache(maxsize=32)
def _scan_results_cached(results_folder: str) -> pd.DataFrame:
    return scan_results(results_folder).copy()


def _warm_result_folder_cache(results_folder: str | Path) -> None:
    folder = str(Path(results_folder).expanduser().resolve())
    with _CACHE_WARM_LOCK:
        if folder in _CACHE_WARMING_FOLDERS:
            return
        _CACHE_WARMING_FOLDERS.add(folder)

    def worker() -> None:
        try:
            index = _scan_results_cached(folder)
            for path in index["path"].tolist():
                _load_final_table_cached(str(path))
        finally:
            with _CACHE_WARM_LOCK:
                _CACHE_WARMING_FOLDERS.discard(folder)

    threading.Thread(target=worker, name="ngsabdiscov-cache-warm", daemon=True).start()


def _finalize_loaded_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["_row_id"] = np.arange(len(df))
    missing = {
        col: ""
        for col in ["target", "library", "CDR3", "cdr3_aa", "vh_scaffold", "vl_scaffold"]
        if col not in df.columns
    }
    if missing:
        df = df.assign(**missing)
    return df.copy()


def load_result_tables(index: pd.DataFrame, target: str = "__ALL__") -> pd.DataFrame:
    if target and target != "__ALL__":
        paths = index.loc[index["target"] == target, "path"].tolist()
    else:
        paths = index["path"].tolist()
    if not paths:
        return _finalize_loaded_data(pd.DataFrame())

    frames = []
    errors = []
    for path in paths:
        try:
            frames.append(_load_final_table_cached(str(path)).copy())
        except Exception as exc:
            errors.append(f"{Path(path).name}: {exc}")
    if not frames:
        raise RuntimeError("Could not load any final lead workbook. " + " | ".join(errors))

    df = pd.concat(frames, ignore_index=True, sort=False)
    return _finalize_loaded_data(df)


def load_results(results_folder: str | Path) -> pd.DataFrame:
    return load_result_tables(scan_results(results_folder), "__ALL__")


@lru_cache(maxsize=128)
def _load_target_data_cached(results_folder: str, target: str) -> pd.DataFrame:
    current_index = _scan_results_cached(results_folder)
    return load_result_tables(current_index, target).copy()


def _round_sort_key(value: str) -> tuple[int, int, str]:
    text = str(value)
    match = re.search(r"(?:round|rnd|r)[_\s-]*(\d+)", text, flags=re.IGNORECASE)
    if not match:
        match = re.search(r"\b(\d+)\b", text)
    return (0, int(match.group(1)), text) if match else (1, 0, text)


def _short_sample_label(freq_col: str) -> str:
    text = str(freq_col).replace("freq ", "")
    parts = text.split("__")
    if len(parts) >= 5:
        return f"{parts[2]} {parts[4]}"
    return text


def _has_fold_change_pair(df: pd.DataFrame, numerator: str | None, denominator: str | None) -> bool:
    return bool(numerator and denominator and numerator != denominator and numerator in df.columns and denominator in df.columns)


def _with_log2_fold_change(df: pd.DataFrame, numerator: str | None, denominator: str | None) -> pd.DataFrame:
    out = df.copy()
    if _has_fold_change_pair(out, numerator, denominator):
        pseudo = 1e-9
        numerator_values = pd.to_numeric(out[numerator], errors="coerce").fillna(0)
        denominator_values = pd.to_numeric(out[denominator], errors="coerce").fillna(0)
        out["log2_fold_change"] = np.log2((numerator_values + pseudo) / (denominator_values + pseudo))
    elif "log2_fold_change" not in out.columns:
        out["log2_fold_change"] = np.nan
    return out


def _filter_data(
    df: pd.DataFrame,
    *,
    target: str,
    library: str,
    min_freq: float,
    min_delphi: float | None,
    min_psr: float | None,
    min_sec: float | None,
    fold_change_mode: str,
    fold_change_cutoff: float | None,
    condition_rule: str,
    negative_control_mode: str,
    max_rank: float | None,
    require_sequence: bool,
    cdr3_search: str,
) -> pd.DataFrame:
    filtered = df.copy()
    if target and target != "__ALL__":
        filtered = filtered[filtered["target"] == target]
    if library and library != "__ALL__" and "library" in filtered.columns:
        filtered = filtered[filtered["library"].astype(str) == library]

    filtered = filtered[pd.to_numeric(filtered["max_freq"], errors="coerce").fillna(0) >= float(min_freq or 0)]

    if min_delphi is not None:
        score = pd.to_numeric(filtered["mean_delphi_score"], errors="coerce")
        if score.notna().any():
            filtered = filtered[score >= float(min_delphi)]
    if min_psr is not None:
        score = pd.to_numeric(filtered["mean_psr_score"], errors="coerce")
        if score.notna().any():
            filtered = filtered[score >= float(min_psr)]
    if min_sec is not None:
        score = pd.to_numeric(filtered["mean_sec_score"], errors="coerce")
        if score.notna().any():
            filtered = filtered[score >= float(min_sec)]
    if fold_change_mode and fold_change_mode != "all":
        fc = pd.to_numeric(filtered["log2_fold_change"], errors="coerce")
        if fc.notna().any():
            cutoff = abs(float(fold_change_cutoff or 0))
            if fold_change_mode == "enriched":
                filtered = filtered[fc >= cutoff]
            elif fold_change_mode == "depleted":
                filtered = filtered[fc <= -cutoff]
            elif fold_change_mode == "either":
                filtered = filtered[fc.abs() >= cutoff]
    if condition_rule and condition_rule != "off":
        condition_scores = pd.to_numeric(filtered["condition_enrichment_score"], errors="coerce")
        if condition_scores.notna().any():
            filtered = filtered[filtered["condition_enrichment_pass"].astype(bool)]
    if negative_control_mode and negative_control_mode != "off" and "selected_over_F_N" in filtered.columns:
        ratios = pd.to_numeric(filtered["selected_over_F_N"], errors="coerce")
        if ratios.notna().any():
            filtered = filtered[filtered["F_N_ratio_pass"].astype(bool)]
    if max_rank is not None and "rank" in filtered.columns:
        filtered = filtered[pd.to_numeric(filtered["rank"], errors="coerce") <= float(max_rank)]
    if require_sequence and "ML_SEQUENCE_OK" in filtered.columns:
        text = filtered["ML_SEQUENCE_OK"].astype(str).str.lower()
        filtered = filtered[text.isin({"true", "1", "yes", "pass", "passed"})]
    if cdr3_search:
        query = cdr3_search.strip().upper()
        cdr3 = filtered["CDR3"].fillna("").astype(str).str.upper()
        cdr3_alt = filtered["cdr3_aa"].fillna("").astype(str).str.upper()
        filtered = filtered[cdr3.str.contains(query, regex=False) | cdr3_alt.str.contains(query, regex=False)]

    return filtered


def _display_columns(df: pd.DataFrame) -> list[str]:
    cols = [col for col in DEFAULT_DISPLAY_COLUMNS if col in df.columns]
    freq_cols = _freq_columns(df)
    useful_freq = sorted(freq_cols, key=_round_sort_key)
    cols.extend(useful_freq[:12])
    return list(dict.fromkeys(cols))


def _column_label(col: str) -> str:
    return _short_sample_label(col) if str(col).startswith("freq ") else col


def _format_for_table(df: pd.DataFrame, columns: list[str], max_rows: int) -> list[dict[str, Any]]:
    payload_columns = list(columns)
    if "_row_id" in df.columns and "_row_id" not in payload_columns:
        payload_columns.append("_row_id")
    table = df[payload_columns].head(max_rows).copy()
    for col in table.columns:
        if pd.api.types.is_numeric_dtype(table[col]):
            table[col] = table[col].round(6)
    return table.replace({np.nan: "", np.inf: "inf", -np.inf: "-inf"}).to_dict("records")


def _selection_context_key(*values: Any) -> str:
    parts = []
    for value in values:
        if isinstance(value, (list, tuple, set)):
            parts.append(",".join(sorted(map(str, value))))
        else:
            parts.append(_clean_text(value))
    return hashlib.sha1("|".join(parts).encode("utf-8")).hexdigest()


def _row_ids_from_plot_selection(selected_data: dict[str, Any] | None, expected_context: str | None = None) -> list[int]:
    if not selected_data:
        return []
    ids: list[int] = []
    for point in selected_data.get("points", []):
        customdata = point.get("customdata")
        if isinstance(customdata, (list, tuple)) and customdata:
            value = customdata[0]
            point_context = customdata[1] if len(customdata) > 1 else None
        else:
            value = customdata
            point_context = None
        if expected_context and point_context != expected_context:
            continue
        try:
            ids.append(int(float(value)))
        except (TypeError, ValueError):
            continue
    return sorted(set(ids))


def _selection_store(context_key: str, selected_row_ids: list[int] | None = None) -> dict[str, Any]:
    return {"context": context_key, "row_ids": selected_row_ids or []}


def _row_ids_from_selection_store(selection_data: Any, context_key: str) -> list[int]:
    if not isinstance(selection_data, dict) or selection_data.get("context") != context_key:
        return []
    ids = []
    for value in selection_data.get("row_ids", []):
        try:
            ids.append(int(float(value)))
        except (TypeError, ValueError):
            continue
    return sorted(set(ids))


def _apply_plot_selection(df: pd.DataFrame, selection_data: Any, context_key: str) -> pd.DataFrame:
    selected_row_ids = _row_ids_from_selection_store(selection_data, context_key)
    if not selected_row_ids or "_row_id" not in df.columns:
        return df
    return df[df["_row_id"].isin(set(selected_row_ids))]


def _edit_distance(left: str, right: str) -> int:
    if left == right:
        return 0
    if not left:
        return len(right)
    if not right:
        return len(left)
    previous = list(range(len(right) + 1))
    for i, aa in enumerate(left, start=1):
        current = [i]
        for j, bb in enumerate(right, start=1):
            current.append(min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (aa != bb)))
        previous = current
    return previous[-1]


def _tree_label(row: pd.Series) -> str:
    cdr3 = _clean_text(row.get("CDR3", row.get("cdr3_aa", "")))
    rank = _clean_text(row.get("rank", ""))
    freq = _clean_text(row.get("max_freq", ""))
    if freq:
        try:
            freq = f"{float(freq):.3g}"
        except ValueError:
            pass
    return "|".join(part for part in [cdr3[:16], f"r{rank}" if rank else "", f"f{freq}" if freq else ""] if part)


def _build_cdr3_tree_figure(rows: list[dict[str, Any]], max_clones: int, go):
    if not rows:
        fig = go.Figure()
        fig.update_layout(title="No clones available for CDR3 tree", height=300)
        return fig

    tree_rows = []
    seen = set()
    for row in rows:
        cdr3 = _clean_text(row.get("CDR3", row.get("cdr3_aa", ""))).upper()
        if not cdr3 or cdr3 in seen:
            continue
        seen.add(cdr3)
        tree_rows.append(row)
        if len(tree_rows) >= max_clones:
            break

    if len(tree_rows) < 2:
        fig = go.Figure()
        fig.update_layout(title="Need at least two unique CDR3s for a tree", height=300)
        return fig

    try:
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import squareform
    except ImportError:
        fig = go.Figure()
        fig.update_layout(title="CDR3 tree requires scipy", height=300)
        return fig

    sequences = [_clean_text(row.get("CDR3", row.get("cdr3_aa", ""))).upper() for row in tree_rows]
    labels = [_tree_label(pd.Series(row)) for row in tree_rows]
    n = len(sequences)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            denom = max(len(sequences[i]), len(sequences[j]), 1)
            distance = _edit_distance(sequences[i], sequences[j]) / denom
            matrix[i, j] = matrix[j, i] = distance
    condensed = squareform(matrix)
    linkage_matrix = linkage(condensed, method="average")
    dendro = dendrogram(linkage_matrix, labels=labels, no_plot=True)

    fig = go.Figure()
    for xs, ys in zip(dendro["icoord"], dendro["dcoord"]):
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line={"color": "#176d7a", "width": 1.4}, hoverinfo="skip"))
    tickvals = [5 + 10 * i for i in range(len(dendro["ivl"]))]
    fig.update_layout(
        title=f"CDR3 Distance Tree ({len(tree_rows)} unique CDR3s)",
        height=300,
        showlegend=False,
        margin={"l": 50, "r": 15, "t": 42, "b": 78},
        xaxis={
            "tickmode": "array",
            "tickvals": tickvals,
            "ticktext": dendro["ivl"],
            "tickangle": 55,
            "tickfont": {"size": 8},
        },
        yaxis_title="Normalized CDR3 edit distance",
    )
    return fig


def create_app(results_folder: str | Path | None = None, *, results_root: str | Path | None = None, max_table_rows: int = 500):
    try:
        from dash import Dash, Input, Output, State, callback_context, dash_table, dcc, html, no_update
        import plotly.express as px
        import plotly.graph_objects as go
    except ImportError as exc:
        raise SystemExit(
            "Dash clone-selection app requires dash and plotly. Install with: pip install dash plotly"
        ) from exc

    if results_folder:
        resolved_folder = Path(results_folder).expanduser().resolve()
        project_index = pd.DataFrame(
            [
                {
                    "label": _result_label(resolved_folder),
                    "results_folder": str(resolved_folder),
                    "run_folder": str(resolved_folder.parent),
                    "results_name": resolved_folder.name,
                }
            ]
        )
    else:
        project_index = scan_project_results(results_root or DEFAULT_RESULTS_ROOT)
    initial_results_folder = str(project_index.iloc[0]["results_folder"])
    result_index = _scan_results_cached(initial_results_folder)
    target_values = sorted(result_index["target"].dropna().astype(str).unique())
    initial_target = target_values[0] if target_values else "__ALL__"
    initial_data = load_result_tables(result_index, initial_target)

    project_options = [
        {"label": row["label"], "value": row["results_folder"]}
        for _, row in project_index.iterrows()
    ]
    target_options = [{"label": target, "value": target} for target in target_values]
    freq_options = [{"label": _short_sample_label(col), "value": col} for col in sorted(_freq_columns(initial_data), key=_round_sort_key)]
    default_numerator, default_denominator = _default_fold_change_columns(initial_data)

    app = Dash(__name__, title="NGSAbDiscov Clone Selection")
    app.index_string = f"""<!DOCTYPE html>
<html>
    <head>
        {{%metas%}}
        <title>{{%title%}}</title>
        {{%favicon%}}
        {{%css%}}
        <style>{APP_CSS}</style>
    </head>
    <body>
        {{%app_entry%}}
        <footer>
            {{%config%}}
            {{%scripts%}}
            {{%renderer%}}
        </footer>
    </body>
</html>"""
    app.layout = html.Div(
        [
            dcc.Store(id="fold-change-selected-row-ids"),
            dcc.Download(id="download-selected"),
            dcc.Download(id="download-filtered"),
            html.Div(
                [
                    html.H1("NGSAbDiscov Clone Selection"),
                    html.Div(id="results-path", className="path"),
                    html.Div(id="summary-cards", className="cards"),
                ],
                className="header",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("MiSeq / results"),
                            dcc.Dropdown(project_options, initial_results_folder, id="results-folder"),
                            html.Label("Target"),
                            dcc.Dropdown(target_options, initial_target, id="target"),
                            html.Label("Min max_freq"),
                            dcc.Input(id="min-freq", type="text", value="0.00005", debounce=True),
                            html.Label("Min Delphi mean score"),
                            dcc.Input(id="min-delphi", type="text", value="0.3", debounce=True),
                            html.Label("Min PSR mean score"),
                            dcc.Input(id="min-psr", type="text", value="", debounce=True),
                            html.Label("Min SEC mean score"),
                            dcc.Input(id="min-sec", type="text", value="", debounce=True),
                            html.Label("Max rank"),
                            dcc.Input(id="max-rank", type="text", value="", debounce=True),
                            html.Label("F_N / 2uM filter"),
                            dcc.Dropdown(NEGATIVE_CONTROL_OPTIONS, "auto", id="negative-control-mode"),
                            html.Div(
                                "Keeps clones where selected condition frequency divided by F_N frequency passes the cutoff. If no F_N/2uM sample exists, this filter is skipped.",
                                className="hint",
                            ),
                            html.Label("Min selected/F_N ratio"),
                            dcc.Input(id="negative-ratio-cutoff", type="text", value="0.3", debounce=True),
                            html.Label("PSR-like F_N/selected cutoff"),
                            dcc.Input(id="negative-psr-cutoff", type="text", value="3", debounce=True),
                            dcc.Checklist(
                                id="diversity-mode",
                                options=[{"label": "One representative per CDR3 cluster", "value": "representative"}],
                                value=["representative"],
                            ),
                            html.Label("CDR3 cluster max mismatches"),
                            dcc.Input(id="cdr3-cluster-distance", type="text", value="2", debounce=True),
                            html.Label("CDR3 contains"),
                            dcc.Input(id="cdr3-search", type="text", value="", debounce=True),
                            dcc.Checklist(
                                id="require-sequence",
                                options=[{"label": "Require valid ML sequence", "value": "yes"}],
                                value=[],
                            ),
                            html.Label("Fold-change comparison"),
                            dcc.Dropdown(
                                FOLD_CHANGE_COMPARISON_OPTIONS,
                                "4nM_vs_20nM",
                                id="fold-change-comparison",
                            ),
                            html.Div(
                                [
                                    html.Label("Manual numerator"),
                                    dcc.Dropdown(freq_options, default_numerator, id="numerator"),
                                    html.Label("Manual denominator"),
                                    dcc.Dropdown(freq_options, default_denominator, id="denominator"),
                                    html.Div(
                                        "Only used when Fold-change comparison is Manual sample columns.",
                                        className="hint",
                                    ),
                                ],
                                id="manual-fold-change-controls",
                                style={"display": "none"},
                            ),
                            html.Label("Fold-change filter"),
                            dcc.Dropdown(
                                [
                                    {"label": "All", "value": "all"},
                                    {"label": "Enriched", "value": "enriched"},
                                    {"label": "Depleted", "value": "depleted"},
                                    {"label": "Either direction", "value": "either"},
                                ],
                                "all",
                                id="fold-change-mode",
                            ),
                            html.Label("Log2 fold-change cutoff"),
                            dcc.Input(id="fold-change-cutoff", type="text", value="1", debounce=True),
                            html.Label("Condition enrichment"),
                            dcc.Dropdown(CONDITION_RULE_OPTIONS, "off", id="condition-rule"),
                            html.Div(
                                "Uses max frequency across columns detected as the same condition, then tests 4nM/20nM, 4nM/100nM, and 20nM/100nM when available.",
                                className="hint",
                            ),
                            html.Label("Condition log2 cutoff"),
                            dcc.Input(id="condition-cutoff", type="text", value="0", debounce=True),
                            html.Label("Max table rows"),
                            dcc.Input(id="max-rows", type="text", value=str(max_table_rows), debounce=True),
                            html.Label("Tree source"),
                            dcc.Dropdown(
                                [
                                    {"label": "Top representative clusters", "value": "representatives"},
                                    {"label": "All clones currently in table", "value": "displayed"},
                                    {"label": "Selected clones", "value": "selected"},
                                ],
                                "representatives",
                                id="tree-source",
                            ),
                            html.Label("Tree clone cap"),
                            dcc.Input(id="tree-max-clones", type="text", value="100", debounce=True),
                            html.Div(
                                "Default tree uses representative CDR3 clusters from the current table.",
                                className="hint",
                            ),
                            html.Button("Download selected CSV", id="download-selected-button", className="action-button secondary-button"),
                            html.Button("Download filtered CSV", id="download-filtered-button", className="action-button secondary-button"),
                        ],
                        className="sidebar",
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Div("Showing all filtered clones", id="plot-selection-status", className="selection-status"),
                                    html.Button(
                                        "Show all filtered clones",
                                        id="clear-plot-selection-button",
                                        className="selection-button",
                                        title="Clear the lasso/box selection and restore the clone table",
                                    ),
                                ],
                                className="selection-toolbar",
                            ),
                            html.Div(
                                dash_table.DataTable(
                                    id="clone-table",
                                    row_selectable="multi",
                                    selected_rows=[],
                                    hidden_columns=["_row_id"],
                                    page_size=12,
                                    sort_action="native",
                                    style_table={
                                        "overflowX": "auto",
                                        "maxHeight": "320px",
                                        "overflowY": "auto",
                                    },
                                    style_cell={
                                        "fontFamily": "Arial",
                                        "fontSize": 11,
                                        "padding": "4px",
                                        "maxWidth": "220px",
                                        "overflow": "hidden",
                                        "textOverflow": "ellipsis",
                                    },
                                    style_header={"fontWeight": "bold", "backgroundColor": "#e8edf4"},
                                    style_data_conditional=[
                                        {
                                            "if": {"filter_query": "{mean_delphi_score} >= 0.5"},
                                            "backgroundColor": "#edf8ef",
                                        }
                                    ],
                                ),
                                className="table-panel",
                            ),
                            html.Div(
                                [
                                    dcc.Tabs(
                                        id="plot-tabs",
                                        value="fold-change",
                                        className="plot-tabs",
                                        children=[
                                            dcc.Tab(
                                                label="Fold Change",
                                                value="fold-change",
                                                children=[
                                                    dcc.Graph(
                                                        id="fold-change-plot",
                                                        className="plot-graph",
                                                        config={"responsive": True},
                                                    )
                                                ],
                                            ),
                                            dcc.Tab(
                                                label="Frequency Profile",
                                                value="frequency-profile",
                                                children=[
                                                    dcc.Graph(
                                                        id="selected-frequency-plot",
                                                        className="plot-graph",
                                                        config={"responsive": True},
                                                    )
                                                ],
                                            ),
                                            dcc.Tab(
                                                label="CDR3 Tree",
                                                value="cdr3-tree",
                                                children=[
                                                    dcc.Graph(
                                                        id="phylogenetic-tree",
                                                        className="plot-graph",
                                                        config={"responsive": True},
                                                    )
                                                ],
                                            ),
                                        ],
                                    ),
                                ],
                                className="plot-panel",
                            ),
                        ],
                        className="main",
                    ),
                ],
                className="workspace",
            ),
        ]
    )

    def data_for_target(selected_results_folder: str, target: str) -> pd.DataFrame:
        selected_results_folder = str(Path(selected_results_folder or initial_results_folder).expanduser().resolve())
        current_index = _scan_results_cached(selected_results_folder)
        target_values_for_folder = set(current_index["target"].dropna().astype(str))
        if not target or target not in target_values_for_folder:
            target = str(current_index.iloc[0]["target"]) if len(current_index) else initial_target
        return _load_target_data_cached(selected_results_folder, target)

    def filtered_from_inputs(
        selected_results_folder,
        target,
        min_freq,
        min_delphi,
        min_psr,
        min_sec,
        fold_change_mode,
        fold_change_cutoff,
        condition_rule,
        condition_cutoff,
        fold_change_comparison,
        numerator,
        denominator,
        negative_control_mode,
        negative_ratio_cutoff,
        negative_psr_cutoff,
        diversity_mode,
        cdr3_cluster_distance,
        max_rank,
        require_sequence,
        cdr3_search,
    ):
        target_data = data_for_target(selected_results_folder, target)
        source_data = _with_condition_enrichment(
            _with_negative_control_summary(
                target_data
            ),
            condition_rule=condition_rule or "off",
            condition_cutoff=_to_optional_float(condition_cutoff),
        )
        source_data = _with_selected_log2_fold_change(
            source_data,
            fold_change_comparison,
            numerator,
            denominator,
        )
        source_data = _with_negative_control_ratio(
            source_data,
            negative_control_mode,
            _to_optional_float(negative_ratio_cutoff),
            _to_optional_float(negative_psr_cutoff),
            fold_change_comparison,
            numerator,
        )
        filtered = _filter_data(
            source_data,
            target="__ALL__",
            library="__ALL__",
            min_freq=_to_float(min_freq, 0.0),
            min_delphi=_to_optional_float(min_delphi),
            min_psr=_to_optional_float(min_psr),
            min_sec=_to_optional_float(min_sec),
            fold_change_mode=fold_change_mode or "all",
            fold_change_cutoff=_to_optional_float(fold_change_cutoff),
            condition_rule=condition_rule or "off",
            negative_control_mode=negative_control_mode or "off",
            max_rank=_to_optional_float(max_rank),
            require_sequence="yes" in (require_sequence or []),
            cdr3_search=cdr3_search or "",
        )
        filtered = _sort_leads(filtered)
        if "representative" in (diversity_mode or []):
            filtered = _with_cdr3_clusters(filtered, max_mismatches=int(_to_float(cdr3_cluster_distance, 2)), max_rows=1000)
            filtered = filtered[filtered["diversity_representative"].astype(bool)]
        return _sort_leads(filtered)

    def plot_selection_context(
        selected_results_folder,
        target,
        min_freq,
        min_delphi,
        min_psr,
        min_sec,
        max_rank,
        require_sequence,
        cdr3_search,
        numerator,
        denominator,
        fold_change_comparison,
        fold_change_mode,
        fold_change_cutoff,
        condition_rule,
        condition_cutoff,
        negative_control_mode,
        negative_ratio_cutoff,
        negative_psr_cutoff,
        diversity_mode,
        cdr3_cluster_distance,
    ):
        return _selection_context_key(
            selected_results_folder,
            target,
            min_freq,
            min_delphi,
            min_psr,
            min_sec,
            max_rank,
            require_sequence,
            cdr3_search,
            numerator,
            denominator,
            fold_change_comparison,
            fold_change_mode,
            fold_change_cutoff,
            condition_rule,
            condition_cutoff,
            negative_control_mode,
            negative_ratio_cutoff,
            negative_psr_cutoff,
            diversity_mode,
            cdr3_cluster_distance,
        )

    def build_fold_change_figure(filtered, fold_change_comparison, numerator, denominator, plot_context_key):
        plot_df = filtered.copy()
        if len(plot_df) > 3000:
            plot_df = plot_df.sample(3000, random_state=11)
        hover_cols = [
            col
            for col in [
                "target",
                "CDR3",
                "vh_scaffold",
                "vl_scaffold",
                "max_freq",
                "log2_fold_change",
                "negative_control_condition",
                "selected_over_F_N",
                "F_N_over_selected",
                "F_N_PSR_like",
                "mean_delphi_score",
                "mean_psr_score",
                "mean_sec_score",
            ]
            if col in plot_df.columns
        ]

        x_label = _fold_change_axis_label(plot_df, fold_change_comparison, numerator, denominator)
        if x_label:
            plot_df["_x"] = pd.to_numeric(plot_df["log2_fold_change"], errors="coerce")
        else:
            plot_df["_x"] = pd.to_numeric(plot_df["max_freq"], errors="coerce")
            x_label = "max_freq"
        plot_df["_rank_log"] = np.log10(pd.to_numeric(plot_df.get("rank", pd.Series(1, index=plot_df.index)), errors="coerce").fillna(1).clip(lower=1))
        plot_df["_score"] = pd.to_numeric(plot_df["mean_delphi_score"], errors="coerce")

        if plot_df.empty:
            fig = go.Figure()
            fig.update_layout(title="No clones match the current filters")
        else:
            scatter_kwargs = {
                "data_frame": plot_df,
                "x": "_x",
                "y": "_rank_log",
                "color": "_score",
                "color_continuous_scale": "Viridis",
                "range_color": (0, 1),
                "hover_data": hover_cols,
                "labels": {"_x": x_label, "_rank_log": "log10(rank)", "_score": "Delphi mean"},
                "title": "Clone Fold Change / Frequency View",
            }
            if "_row_id" in plot_df.columns:
                plot_df["_selection_context"] = plot_context_key
                scatter_kwargs["custom_data"] = ["_row_id", "_selection_context"]
            fig = px.scatter(**scatter_kwargs)
            fig.update_traces(marker={"size": 7, "opacity": 0.78})
            if x_label != "max_freq":
                fig.add_vline(x=0, line_dash="dash", line_color="#888888")
        fig.update_layout(
            height=320,
            margin={"l": 55, "r": 15, "t": 45, "b": 45},
            dragmode="select",
            clickmode="event+select",
            uirevision=plot_context_key,
            selectionrevision=plot_context_key,
        )
        return fig

    @app.callback(
        Output("target", "options"),
        Output("target", "value"),
        Output("results-path", "children"),
        Input("results-folder", "value"),
    )
    def update_targets_for_results(selected_results_folder):
        selected_results_folder = selected_results_folder or initial_results_folder
        _warm_result_folder_cache(selected_results_folder)
        current_index = _scan_results_cached(selected_results_folder)
        targets = sorted(current_index["target"].dropna().astype(str).unique())
        options = [{"label": target, "value": target} for target in targets]
        value = targets[0] if targets else None
        return options, value, selected_results_folder

    @app.callback(
        Output("fold-change-comparison", "options"),
        Output("fold-change-comparison", "value"),
        Output("numerator", "options"),
        Output("numerator", "value"),
        Output("denominator", "options"),
        Output("denominator", "value"),
        Input("results-folder", "value"),
        Input("target", "value"),
    )
    def update_target_options(selected_results_folder, target):
        source_data = data_for_target(selected_results_folder, target)
        comparison_options_for_target = _available_fold_change_options(source_data)
        comparison_value = comparison_options_for_target[0]["value"] if comparison_options_for_target else "manual"
        freq_options_for_target = [
            {"label": _short_sample_label(col), "value": col}
            for col in sorted(_freq_columns(source_data), key=_round_sort_key)
        ]
        numerator_value, denominator_value = _default_fold_change_columns(source_data)

        return (
            comparison_options_for_target,
            comparison_value,
            freq_options_for_target,
            numerator_value,
            freq_options_for_target,
            denominator_value,
        )

    @app.callback(
        Output("manual-fold-change-controls", "style"),
        Input("fold-change-comparison", "value"),
    )
    def toggle_manual_fold_change_controls(fold_change_comparison):
        if fold_change_comparison == "manual":
            return {"display": "block"}
        return {"display": "none"}

    @app.callback(
        Output("fold-change-selected-row-ids", "data"),
        Input("fold-change-plot", "selectedData"),
        Input("clear-plot-selection-button", "n_clicks"),
        Input("results-folder", "value"),
        Input("target", "value"),
        Input("min-freq", "value"),
        Input("min-delphi", "value"),
        Input("min-psr", "value"),
        Input("min-sec", "value"),
        Input("max-rank", "value"),
        Input("require-sequence", "value"),
        Input("cdr3-search", "value"),
        Input("numerator", "value"),
        Input("denominator", "value"),
        Input("fold-change-comparison", "value"),
        Input("fold-change-mode", "value"),
        Input("fold-change-cutoff", "value"),
        Input("condition-rule", "value"),
        Input("condition-cutoff", "value"),
        Input("negative-control-mode", "value"),
        Input("negative-ratio-cutoff", "value"),
        Input("negative-psr-cutoff", "value"),
        Input("diversity-mode", "value"),
        Input("cdr3-cluster-distance", "value"),
    )
    def update_fold_change_selection(selected_data, _clear_clicks, *_filter_values):
        context_key = plot_selection_context(*_filter_values)
        triggered = callback_context.triggered[0]["prop_id"].split(".")[0] if callback_context.triggered else ""
        if triggered == "fold-change-plot":
            return _selection_store(context_key, _row_ids_from_plot_selection(selected_data, context_key))
        return _selection_store(context_key)

    @app.callback(
        Output("fold-change-plot", "selectedData"),
        Input("clear-plot-selection-button", "n_clicks"),
        Input("results-folder", "value"),
        Input("target", "value"),
        Input("min-freq", "value"),
        Input("min-delphi", "value"),
        Input("min-psr", "value"),
        Input("min-sec", "value"),
        Input("max-rank", "value"),
        Input("require-sequence", "value"),
        Input("cdr3-search", "value"),
        Input("numerator", "value"),
        Input("denominator", "value"),
        Input("fold-change-comparison", "value"),
        Input("fold-change-mode", "value"),
        Input("fold-change-cutoff", "value"),
        Input("condition-rule", "value"),
        Input("condition-cutoff", "value"),
        Input("negative-control-mode", "value"),
        Input("negative-ratio-cutoff", "value"),
        Input("negative-psr-cutoff", "value"),
        Input("diversity-mode", "value"),
        Input("cdr3-cluster-distance", "value"),
        prevent_initial_call=True,
    )
    def clear_fold_change_plot_selection(*_inputs):
        return None

    @app.callback(
        Output("summary-cards", "children"),
        Output("plot-selection-status", "children"),
        Output("clone-table", "data"),
        Output("clone-table", "columns"),
        Output("clone-table", "selected_rows"),
        Input("results-folder", "value"),
        Input("target", "value"),
        Input("min-freq", "value"),
        Input("min-delphi", "value"),
        Input("min-psr", "value"),
        Input("min-sec", "value"),
        Input("max-rank", "value"),
        Input("require-sequence", "value"),
        Input("cdr3-search", "value"),
        Input("numerator", "value"),
        Input("denominator", "value"),
        Input("fold-change-comparison", "value"),
        Input("fold-change-mode", "value"),
        Input("fold-change-cutoff", "value"),
        Input("condition-rule", "value"),
        Input("condition-cutoff", "value"),
        Input("negative-control-mode", "value"),
        Input("negative-ratio-cutoff", "value"),
        Input("negative-psr-cutoff", "value"),
        Input("diversity-mode", "value"),
        Input("cdr3-cluster-distance", "value"),
        Input("max-rows", "value"),
        Input("fold-change-selected-row-ids", "data"),
    )
    def update_table_and_plot(
        selected_results_folder,
        target,
        min_freq,
        min_delphi,
        min_psr,
        min_sec,
        max_rank,
        require_sequence,
        cdr3_search,
        numerator,
        denominator,
        fold_change_comparison,
        fold_change_mode,
        fold_change_cutoff,
        condition_rule,
        condition_cutoff,
        negative_control_mode,
        negative_ratio_cutoff,
        negative_psr_cutoff,
        diversity_mode,
        cdr3_cluster_distance,
        max_rows,
        plot_selected_row_ids,
    ):
        filtered = filtered_from_inputs(
            selected_results_folder,
            target,
            min_freq,
            min_delphi,
            min_psr,
            min_sec,
            fold_change_mode,
            fold_change_cutoff,
            condition_rule,
            condition_cutoff,
            fold_change_comparison,
            numerator,
            denominator,
            negative_control_mode,
            negative_ratio_cutoff,
            negative_psr_cutoff,
            diversity_mode,
            cdr3_cluster_distance,
            max_rank,
            require_sequence,
            cdr3_search,
        )
        context_key = plot_selection_context(
            selected_results_folder,
            target,
            min_freq,
            min_delphi,
            min_psr,
            min_sec,
            max_rank,
            require_sequence,
            cdr3_search,
            numerator,
            denominator,
            fold_change_comparison,
            fold_change_mode,
            fold_change_cutoff,
            condition_rule,
            condition_cutoff,
            negative_control_mode,
            negative_ratio_cutoff,
            negative_psr_cutoff,
            diversity_mode,
            cdr3_cluster_distance,
        )
        selected_row_ids = _row_ids_from_selection_store(plot_selected_row_ids, context_key)
        table_filtered = _apply_plot_selection(filtered, plot_selected_row_ids, context_key)
        max_rows = min(int(_to_float(max_rows, float(max_table_rows))), 500)
        column_source = table_filtered if not table_filtered.empty else filtered
        columns = _display_columns(column_source)
        table_data = _format_for_table(table_filtered, columns, max_rows)
        dash_column_ids = [*columns, "_row_id"] if "_row_id" in column_source.columns else columns
        dash_columns = [{"name": _column_label(col), "id": col} for col in dash_column_ids]

        scored = pd.to_numeric(table_filtered["mean_delphi_score"], errors="coerce").notna().sum()
        candidate_label = "Plot-selected candidates" if selected_row_ids else "Lead candidates"
        cards = [
            html.Div([html.Strong(f"{len(table_filtered):,}"), html.Span(candidate_label)], className="card"),
            html.Div([html.Strong(f"{table_filtered['cdr3_cluster'].replace('', np.nan).nunique():,}" if "cdr3_cluster" in table_filtered.columns else "0"), html.Span("CDR3 clusters")], className="card"),
            html.Div([html.Strong(f"{pd.to_numeric(table_filtered['max_freq'], errors='coerce').max(skipna=True):.4g}" if len(table_filtered) else "0"), html.Span("Top max_freq")], className="card"),
            html.Div([html.Strong(f"{scored:,}"), html.Span("ML-scored rows")], className="card"),
        ]

        if selected_row_ids:
            selection_status = f"Lasso selection active: {len(table_filtered):,} clones. Click Show all filtered clones to undo."
        else:
            selection_status = "Showing all filtered clones"

        return cards, selection_status, table_data, dash_columns, []

    @app.callback(
        Output("fold-change-plot", "figure"),
        Input("results-folder", "value"),
        Input("target", "value"),
        Input("min-freq", "value"),
        Input("min-delphi", "value"),
        Input("min-psr", "value"),
        Input("min-sec", "value"),
        Input("max-rank", "value"),
        Input("require-sequence", "value"),
        Input("cdr3-search", "value"),
        Input("numerator", "value"),
        Input("denominator", "value"),
        Input("fold-change-comparison", "value"),
        Input("fold-change-mode", "value"),
        Input("fold-change-cutoff", "value"),
        Input("condition-rule", "value"),
        Input("condition-cutoff", "value"),
        Input("negative-control-mode", "value"),
        Input("negative-ratio-cutoff", "value"),
        Input("negative-psr-cutoff", "value"),
        Input("diversity-mode", "value"),
        Input("cdr3-cluster-distance", "value"),
    )
    def update_fold_change_plot(
        selected_results_folder,
        target,
        min_freq,
        min_delphi,
        min_psr,
        min_sec,
        max_rank,
        require_sequence,
        cdr3_search,
        numerator,
        denominator,
        fold_change_comparison,
        fold_change_mode,
        fold_change_cutoff,
        condition_rule,
        condition_cutoff,
        negative_control_mode,
        negative_ratio_cutoff,
        negative_psr_cutoff,
        diversity_mode,
        cdr3_cluster_distance,
    ):
        filtered = filtered_from_inputs(
            selected_results_folder,
            target,
            min_freq,
            min_delphi,
            min_psr,
            min_sec,
            fold_change_mode,
            fold_change_cutoff,
            condition_rule,
            condition_cutoff,
            fold_change_comparison,
            numerator,
            denominator,
            negative_control_mode,
            negative_ratio_cutoff,
            negative_psr_cutoff,
            diversity_mode,
            cdr3_cluster_distance,
            max_rank,
            require_sequence,
            cdr3_search,
        )
        context_key = plot_selection_context(
            selected_results_folder,
            target,
            min_freq,
            min_delphi,
            min_psr,
            min_sec,
            max_rank,
            require_sequence,
            cdr3_search,
            numerator,
            denominator,
            fold_change_comparison,
            fold_change_mode,
            fold_change_cutoff,
            condition_rule,
            condition_cutoff,
            negative_control_mode,
            negative_ratio_cutoff,
            negative_psr_cutoff,
            diversity_mode,
            cdr3_cluster_distance,
        )
        return build_fold_change_figure(filtered, fold_change_comparison, numerator, denominator, context_key)

    @app.callback(
        Output("selected-frequency-plot", "figure"),
        Input("clone-table", "data"),
        Input("clone-table", "selected_rows"),
    )
    def update_selected_profile(table_rows, selected_rows):
        if not table_rows or not selected_rows:
            fig = go.Figure()
            fig.update_layout(title="Select one or more table rows to view frequency profiles", height=300)
            return fig

        selected = pd.DataFrame([table_rows[idx] for idx in selected_rows if idx < len(table_rows)])
        freq_cols = sorted(_freq_columns(selected), key=_round_sort_key)
        if selected.empty or not freq_cols:
            fig = go.Figure()
            fig.update_layout(title="No frequency columns available for selected clones", height=300)
            return fig

        rows = []
        for row_idx, row in selected.head(10).iterrows():
            label_parts = [row.get("vh_scaffold", ""), row.get("vl_scaffold", ""), row.get("CDR3", "")]
            label = "|".join(part for part in map(str, label_parts) if part) or f"row {row_idx}"
            for col in freq_cols:
                rows.append({"clone": label[:100], "sample": _short_sample_label(col), "freq": row.get(col, 0), "order": _round_sort_key(col)})
        profile = pd.DataFrame(rows)
        profile["freq"] = pd.to_numeric(profile["freq"], errors="coerce").fillna(0)
        profile = profile.sort_values("order")
        fig = px.line(profile, x="sample", y="freq", color="clone", markers=True, title="Selected Clone Frequency Profiles")
        fig.update_layout(height=300, xaxis_tickangle=45, margin={"l": 55, "r": 15, "t": 42, "b": 72})
        return fig

    @app.callback(
        Output("phylogenetic-tree", "figure"),
        Input("plot-tabs", "value"),
        Input("clone-table", "data"),
        Input("clone-table", "selected_rows"),
        Input("tree-source", "value"),
        Input("tree-max-clones", "value"),
    )
    def update_phylogenetic_tree(active_tab, table_rows, selected_rows, tree_source, tree_max_clones):
        if active_tab != "cdr3-tree":
            fig = go.Figure()
            fig.update_layout(title="Open the CDR3 Tree tab to build the tree", height=300)
            return fig
        if not table_rows:
            fig = go.Figure()
            fig.update_layout(title="No clones available to build a CDR3 tree", height=300)
            return fig
        max_clones = int(_to_float(tree_max_clones, 100))
        if tree_source == "selected":
            rows = [table_rows[idx] for idx in (selected_rows or []) if idx < len(table_rows)]
            if not rows:
                fig = go.Figure()
                fig.update_layout(title="Select clones to build a CDR3 tree", height=300)
                return fig
        elif tree_source == "representatives":
            rows = [row for row in table_rows if row.get("diversity_representative") in {True, "True", "true", "1", 1}]
            if not rows:
                rows = table_rows
        else:
            rows = table_rows
        return _build_cdr3_tree_figure(rows, max_clones=max_clones, go=go)

    @app.callback(
        Output("download-selected", "data"),
        Input("download-selected-button", "n_clicks"),
        State("clone-table", "data"),
        State("clone-table", "selected_rows"),
        State("results-folder", "value"),
        State("target", "value"),
        State("min-freq", "value"),
        State("min-delphi", "value"),
        State("min-psr", "value"),
        State("min-sec", "value"),
        State("max-rank", "value"),
        State("require-sequence", "value"),
        State("cdr3-search", "value"),
        State("numerator", "value"),
        State("denominator", "value"),
        State("fold-change-comparison", "value"),
        State("fold-change-mode", "value"),
        State("fold-change-cutoff", "value"),
        State("condition-rule", "value"),
        State("condition-cutoff", "value"),
        State("negative-control-mode", "value"),
        State("negative-ratio-cutoff", "value"),
        State("negative-psr-cutoff", "value"),
        State("diversity-mode", "value"),
        State("cdr3-cluster-distance", "value"),
        prevent_initial_call=True,
    )
    def download_selected(
        _clicks,
        table_rows,
        selected_rows,
        selected_results_folder,
        target,
        min_freq,
        min_delphi,
        min_psr,
        min_sec,
        max_rank,
        require_sequence,
        cdr3_search,
        numerator,
        denominator,
        fold_change_comparison,
        fold_change_mode,
        fold_change_cutoff,
        condition_rule,
        condition_cutoff,
        negative_control_mode,
        negative_ratio_cutoff,
        negative_psr_cutoff,
        diversity_mode,
        cdr3_cluster_distance,
    ):
        if not table_rows or not selected_rows:
            return no_update
        selected_ids = {
            int(table_rows[idx]["_row_id"])
            for idx in selected_rows
            if idx < len(table_rows) and "_row_id" in table_rows[idx]
        }
        if not selected_ids:
            return no_update
        filtered = filtered_from_inputs(
            selected_results_folder,
            target,
            min_freq,
            min_delphi,
            min_psr,
            min_sec,
            fold_change_mode,
            fold_change_cutoff,
            condition_rule,
            condition_cutoff,
            fold_change_comparison,
            numerator,
            denominator,
            negative_control_mode,
            negative_ratio_cutoff,
            negative_psr_cutoff,
            diversity_mode,
            cdr3_cluster_distance,
            max_rank,
            require_sequence,
            cdr3_search,
        )
        selected = filtered[filtered["_row_id"].isin(selected_ids)]
        return dcc.send_data_frame(selected.to_csv, "selected_clones.csv", index=False)

    @app.callback(
        Output("download-filtered", "data"),
        Input("download-filtered-button", "n_clicks"),
        State("results-folder", "value"),
        State("target", "value"),
        State("min-freq", "value"),
        State("min-delphi", "value"),
        State("min-psr", "value"),
        State("min-sec", "value"),
        State("max-rank", "value"),
        State("require-sequence", "value"),
        State("cdr3-search", "value"),
        State("numerator", "value"),
        State("denominator", "value"),
        State("fold-change-comparison", "value"),
        State("fold-change-mode", "value"),
        State("fold-change-cutoff", "value"),
        State("condition-rule", "value"),
        State("condition-cutoff", "value"),
        State("negative-control-mode", "value"),
        State("negative-ratio-cutoff", "value"),
        State("negative-psr-cutoff", "value"),
        State("diversity-mode", "value"),
        State("cdr3-cluster-distance", "value"),
        State("fold-change-selected-row-ids", "data"),
        prevent_initial_call=True,
    )
    def download_filtered(
        _clicks,
        selected_results_folder,
        target,
        min_freq,
        min_delphi,
        min_psr,
        min_sec,
        max_rank,
        require_sequence,
        cdr3_search,
        numerator,
        denominator,
        fold_change_comparison,
        fold_change_mode,
        fold_change_cutoff,
        condition_rule,
        condition_cutoff,
        negative_control_mode,
        negative_ratio_cutoff,
        negative_psr_cutoff,
        diversity_mode,
        cdr3_cluster_distance,
        plot_selected_row_ids,
    ):
        filtered = filtered_from_inputs(
            selected_results_folder,
            target,
            min_freq,
            min_delphi,
            min_psr,
            min_sec,
            fold_change_mode,
            fold_change_cutoff,
            condition_rule,
            condition_cutoff,
            fold_change_comparison,
            numerator,
            denominator,
            negative_control_mode,
            negative_ratio_cutoff,
            negative_psr_cutoff,
            diversity_mode,
            cdr3_cluster_distance,
            max_rank,
            require_sequence,
            cdr3_search,
        )
        context_key = plot_selection_context(
            selected_results_folder,
            target,
            min_freq,
            min_delphi,
            min_psr,
            min_sec,
            max_rank,
            require_sequence,
            cdr3_search,
            numerator,
            denominator,
            fold_change_comparison,
            fold_change_mode,
            fold_change_cutoff,
            condition_rule,
            condition_cutoff,
            negative_control_mode,
            negative_ratio_cutoff,
            negative_psr_cutoff,
            diversity_mode,
            cdr3_cluster_distance,
        )
        filtered = _apply_plot_selection(filtered, plot_selected_row_ids, context_key)
        if filtered.empty:
            return no_update
        return dcc.send_data_frame(filtered.to_csv, "filtered_clones.csv", index=False)

    return app


def run_app(
    results_folder: str | Path | None = None,
    results_root: str | Path | None = None,
    host: str = "127.0.0.1",
    port: int = 8050,
    debug: bool = False,
    max_table_rows: int = 500,
) -> None:
    app = create_app(results_folder, results_root=results_root, max_table_rows=max_table_rows)
    try:
        app.run(host=host, port=port, debug=debug)
    except AttributeError:
        app.run_server(host=host, port=port, debug=debug)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the NGSAbDiscov clone-selection Dash app.")
    parser.add_argument("--results-folder", default=None, help="Single NGSAbDiscov results folder")
    parser.add_argument("--results-root", default=None, help="Root folder containing multiple MiSeq NGSAbDiscov results")
    parser.add_argument("--host", default="127.0.0.1", help="Bind host")
    parser.add_argument("--port", type=int, default=8050, help="Bind port")
    parser.add_argument("--debug", action="store_true", help="Run Dash in debug mode")
    parser.add_argument("--max-table-rows", type=int, default=500, help="Maximum rows sent to the browser table")
    args = parser.parse_args()
    run_app(
        args.results_folder,
        results_root=args.results_root,
        host=args.host,
        port=args.port,
        debug=args.debug,
        max_table_rows=args.max_table_rows,
    )


if __name__ == "__main__":
    main()
