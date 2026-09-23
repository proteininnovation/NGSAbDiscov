"""Stage B bridge between NGSAbDiscov lead tables and AbForge/OpenDDE.

Stage A owns NGS processing, enrichment, PSR/SEC, and first-pass lead picking.
This module prepares those selected clones for structural scoring, then merges
AbForge/OpenDDE scores back into selection-ready tables.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


FINAL_LEAD_GLOB = "by_protein/*_final_leads.xlsx"
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYX")

HSEQ_COLUMNS = [
    "HSEQ",
    "hseq",
    "VH",
    "VH_sequence",
    "vh_sequence",
    "vhh_sequence",
    "heavy_sequence",
    "heavy",
]
LSEQ_COLUMNS = ["LSEQ", "lseq", "VL", "VL_sequence", "vl_sequence", "light_sequence", "light"]
ANTIGEN_SEQUENCE_COLUMNS = [
    "antigen_aa",
    "antigen_sequence",
    "ag_sequence",
    "target_sequence",
    "Target Sequence",
    "protein_sequence",
    "sequence_antigen",
]
TARGET_COLUMNS = ["target", "Target", "antigen", "Antigen", "protein", "Protein"]
CLONE_ID_COLUMNS = ["clone_id", "Clone", "clone", "Antibody", "antibody", "BARCODE", "TAB_ID", "id", "ID"]
CDR3_COLUMNS = ["CDR3", "cdr3_aa", "HCDR3", "hcdr3"]

CORE_EXPORT_COLUMNS = [
    "stage_b_example_id",
    "stage_b_ready",
    "stage_b_issues",
    "target",
    "clone_id",
    "antibody_format",
    "hseq",
    "lseq",
    "antigen_sequence",
    "cdr3",
    "stage_a_rank",
    "max_freq",
    "mean_delphi_score",
    "mean_psr_score",
    "mean_sec_score",
    "source_file",
    "source_row",
    "source_sheet",
]

SCORE_ID_COLUMNS = ["stage_b_example_id", "example_id", "input_id", "job_name", "design_id", "record_id"]
SCORE_COLUMN_ALIASES = {
    "interface_viability_label_pred": "abforge_interface_viability",
    "interface_viability_pred": "abforge_interface_viability",
    "interface_viability_probability": "abforge_interface_viability",
    "p_interface_viability": "abforge_interface_viability",
    "abforge_interface_viability": "abforge_interface_viability",
    "binder_label_pred": "abforge_binder_probability",
    "binder_pred": "abforge_binder_probability",
    "binder_probability": "abforge_binder_probability",
    "p_binder": "abforge_binder_probability",
    "abforge_binder_probability": "abforge_binder_probability",
    "pose_success_label_pred": "abforge_pose_success",
    "pose_success_pred": "abforge_pose_success",
    "pose_success_probability": "abforge_pose_success",
    "p_pose_success": "abforge_pose_success",
    "abforge_pose_success": "abforge_pose_success",
    "dockq_pred": "abforge_dockq_pred",
    "pkd_pred": "abforge_pkd_pred",
    "iptm": "abforge_iptm",
    "best_iptm": "abforge_iptm",
    "ptm": "abforge_ptm",
    "ranking_score": "abforge_ranking_score",
    "min_ipsae_20": "abforge_min_ipsae_20",
    "cif_path": "abforge_model_path",
    "model_path": "abforge_model_path",
    "prediction_cif": "abforge_model_path",
    "confidence_json": "abforge_confidence_json",
}

SELECTION_WEIGHTS = {
    "abforge_interface_viability": 0.45,
    "abforge_pose_success": 0.25,
    "abforge_binder_probability": 0.20,
    "abforge_iptm": 0.10,
}


@dataclass(frozen=True)
class StageBExportResult:
    output_csv: Path
    ready_csv: Path
    summary_json: Path
    rows: int
    ready_rows: int
    missing_hseq_rows: int
    missing_antigen_rows: int
    targets: int


@dataclass(frozen=True)
class StageBImportResult:
    output_dir: Path
    combined_csv: Path
    combined_xlsx: Path
    summary_json: Path
    rows: int
    scored_rows: int
    files_written: int


def export_stage_b_candidates(
    results_folder: str | Path,
    *,
    output_csv: str | Path | None = None,
    target_sequence_table: str | Path | None = None,
    top_per_target: int | None = None,
    max_total: int | None = None,
) -> StageBExportResult:
    """Export Stage A final leads into an AbForge-ready candidate CSV."""

    folder = Path(results_folder)
    output_csv = Path(output_csv) if output_csv else folder / "stage_b" / "abforge_stage_b_candidates.csv"
    ready_csv = output_csv.with_name(output_csv.stem.replace("_candidates", "_ready") + output_csv.suffix)
    summary_json = output_csv.with_suffix(".summary.json")

    target_sequences = load_target_sequences(target_sequence_table) if target_sequence_table else {}
    rows = build_stage_b_candidate_rows(folder, target_sequences=target_sequences, top_per_target=top_per_target)
    if max_total is not None and max_total > 0:
        rows = _sort_candidate_rows(rows).head(max_total).to_dict("records")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    candidate_df = pd.DataFrame(rows)
    candidate_df = _ensure_columns(candidate_df, CORE_EXPORT_COLUMNS)
    candidate_df.to_csv(output_csv, index=False)
    ready_df = candidate_df[candidate_df["stage_b_ready"].astype(bool)].copy()
    ready_df.to_csv(ready_csv, index=False)

    summary = summarize_export(candidate_df, folder, output_csv, ready_csv, target_sequence_table)
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return StageBExportResult(
        output_csv=output_csv,
        ready_csv=ready_csv,
        summary_json=summary_json,
        rows=int(summary["rows"]),
        ready_rows=int(summary["ready_rows"]),
        missing_hseq_rows=int(summary["missing_hseq_rows"]),
        missing_antigen_rows=int(summary["missing_antigen_rows"]),
        targets=int(summary["targets"]),
    )


def import_stage_b_scores(
    results_folder: str | Path,
    scores_csv: str | Path,
    *,
    output_dir: str | Path | None = None,
    target_sequence_table: str | Path | None = None,
    write_back: bool = False,
) -> StageBImportResult:
    """Merge AbForge/OpenDDE Stage B scores back into final lead tables."""

    folder = Path(results_folder)
    scores_path = Path(scores_csv)
    output_dir = Path(output_dir) if output_dir else folder / "stage_b"
    enriched_dir = output_dir / "enriched_by_protein"
    enriched_dir.mkdir(parents=True, exist_ok=True)

    target_sequences = load_target_sequences(target_sequence_table) if target_sequence_table else {}
    scores = normalize_scores(pd.read_csv(scores_path))
    score_columns = [col for col in scores.columns if col != "stage_b_example_id"]

    all_enriched: list[pd.DataFrame] = []
    files_written = 0
    for lead_path in final_lead_files(folder):
        lead_df = pd.read_excel(lead_path)
        target = target_from_lead_path(lead_path)
        annotated = annotate_existing_lead_table(lead_df, lead_path.name, target, target_sequences)
        enriched = annotated.merge(scores, on="stage_b_example_id", how="left")
        for col in score_columns:
            if col not in enriched.columns:
                enriched[col] = pd.NA
        enriched = add_structural_selection_score(enriched)
        enriched["stage_b_score_source"] = scores_path.name

        out_path = lead_path if write_back else enriched_dir / f"{lead_path.stem}_stage_b.xlsx"
        enriched.to_excel(out_path, index=False)
        files_written += 1
        all_enriched.append(enriched)

    combined = pd.concat(all_enriched, ignore_index=True) if all_enriched else pd.DataFrame()
    combined = _sort_enriched_rows(combined)
    combined_csv = output_dir / "abforge_stage_b_enriched_leads.csv"
    combined_xlsx = output_dir / "abforge_stage_b_enriched_leads.xlsx"
    combined.to_csv(combined_csv, index=False)
    combined.to_excel(combined_xlsx, index=False)

    scored_rows = int(combined["abforge_structural_selection_score"].notna().sum()) if "abforge_structural_selection_score" in combined else 0
    summary = {
        "results_folder": str(folder),
        "scores_csv": str(scores_path),
        "output_dir": str(output_dir),
        "write_back": bool(write_back),
        "rows": int(len(combined)),
        "scored_rows": scored_rows,
        "unscored_rows": int(len(combined) - scored_rows),
        "files_written": files_written,
        "combined_csv": str(combined_csv),
        "combined_xlsx": str(combined_xlsx),
    }
    summary_json = output_dir / "abforge_stage_b_import_summary.json"
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return StageBImportResult(
        output_dir=output_dir,
        combined_csv=combined_csv,
        combined_xlsx=combined_xlsx,
        summary_json=summary_json,
        rows=int(len(combined)),
        scored_rows=scored_rows,
        files_written=files_written,
    )


def build_stage_b_candidate_rows(
    folder: Path,
    *,
    target_sequences: dict[str, str] | None = None,
    top_per_target: int | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    target_sequences = target_sequences or {}
    for path in final_lead_files(folder):
        df = pd.read_excel(path)
        if top_per_target is not None and top_per_target > 0:
            df = _sort_lead_rows(df).head(top_per_target).copy()
        target = target_from_lead_path(path)
        annotated = annotate_existing_lead_table(df, path.name, target, target_sequences)
        rows.extend(annotated.to_dict("records"))
    return rows


def annotate_existing_lead_table(
    df: pd.DataFrame,
    source_file: str,
    target: str,
    target_sequences: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Add deterministic Stage B fields to an existing final-leads table."""

    target_sequences = target_sequences or {}
    records = []
    for source_row, (_, row) in enumerate(df.iterrows(), start=1):
        record = dict(row)
        target_value = first_text(row, TARGET_COLUMNS) or target
        hseq = clean_sequence(first_text(row, HSEQ_COLUMNS))
        lseq = clean_sequence(first_text(row, LSEQ_COLUMNS))
        antigen_sequence = clean_sequence(first_text(row, ANTIGEN_SEQUENCE_COLUMNS))
        if not antigen_sequence:
            antigen_sequence = target_sequences.get(normalize_key(target_value), "")
        cdr3 = clean_text(first_text(row, CDR3_COLUMNS))
        clone_id = clean_text(first_text(row, CLONE_ID_COLUMNS)) or f"{safe_token(target_value)}_{source_row:05d}"
        antibody_format = infer_antibody_format(row, lseq=lseq)
        example_id = stage_b_example_id(target_value, clone_id, hseq, lseq, antigen_sequence, source_row)
        issues = stage_b_issues(hseq, antigen_sequence)

        record.update(
            {
                "stage_b_example_id": example_id,
                "stage_b_ready": len(issues) == 0,
                "stage_b_issues": ";".join(issues),
                "target": target_value,
                "clone_id": clone_id,
                "antibody_format": antibody_format,
                "hseq": hseq,
                "lseq": lseq,
                "antigen_sequence": antigen_sequence,
                "cdr3": cdr3,
                "stage_a_rank": numeric_or_blank(first_text(row, ["rank", "Rank", "stage_a_rank"])),
                "source_file": source_file,
                "source_row": source_row,
                "source_sheet": "Sheet1",
            }
        )
        records.append(record)
    out = pd.DataFrame(records)
    return _ensure_columns(out, CORE_EXPORT_COLUMNS)


def final_lead_files(folder: Path) -> list[Path]:
    files = sorted(folder.glob(FINAL_LEAD_GLOB))
    if not files:
        raise FileNotFoundError(f"No {FINAL_LEAD_GLOB} files found in {folder}")
    return files


def load_target_sequences(path: str | Path | None) -> dict[str, str]:
    if not path:
        return {}
    table_path = Path(path)
    if not table_path.exists():
        raise FileNotFoundError(f"Stage B target sequence table not found: {table_path}")
    df = read_table(table_path)
    target_col = first_existing_column(df, TARGET_COLUMNS)
    sequence_col = first_existing_column(df, ANTIGEN_SEQUENCE_COLUMNS + ["sequence", "aa"])
    if target_col is None or sequence_col is None:
        raise ValueError(
            f"Target sequence table must contain target and sequence columns. Found columns: {', '.join(map(str, df.columns))}"
        )
    mapping: dict[str, str] = {}
    for _, row in df.iterrows():
        key = normalize_key(row.get(target_col))
        sequence = clean_sequence(row.get(sequence_col))
        if key and sequence:
            mapping[key] = sequence
    return mapping


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        return pd.read_excel(path)
    if suffix in {".tsv", ".tab"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def normalize_scores(scores: pd.DataFrame) -> pd.DataFrame:
    id_col = first_existing_column(scores, SCORE_ID_COLUMNS)
    if id_col is None:
        raise ValueError(f"Score table needs one ID column from: {', '.join(SCORE_ID_COLUMNS)}")

    normalized = pd.DataFrame({"stage_b_example_id": scores[id_col].astype(str).map(clean_text)})
    for col in scores.columns:
        if col == id_col:
            continue
        mapped = SCORE_COLUMN_ALIASES.get(str(col), f"abforge_{safe_column_name(str(col))}")
        if mapped not in normalized.columns:
            normalized[mapped] = scores[col]
    normalized = normalized.drop_duplicates("stage_b_example_id", keep="first")
    return normalized


def add_structural_selection_score(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    scores: list[float | None] = []
    for _, row in out.iterrows():
        weighted = 0.0
        weight_sum = 0.0
        for col, weight in SELECTION_WEIGHTS.items():
            value = probability_value(row.get(col))
            if value is None:
                continue
            weighted += weight * value
            weight_sum += weight
        scores.append(weighted / weight_sum if weight_sum else None)
    out["abforge_structural_selection_score"] = scores
    out["abforge_selection_band"] = out["abforge_structural_selection_score"].apply(selection_band)
    return out


def selection_band(value: Any) -> str:
    number = float_or_none(value)
    if number is None:
        return "missing"
    if number >= 0.75:
        return "strong"
    if number >= 0.60:
        return "support"
    if number >= 0.45:
        return "review"
    return "low"


def summarize_export(
    df: pd.DataFrame,
    folder: Path,
    output_csv: Path,
    ready_csv: Path,
    target_sequence_table: str | Path | None,
) -> dict[str, Any]:
    ready = df["stage_b_ready"].astype(bool) if "stage_b_ready" in df else pd.Series([], dtype=bool)
    missing_hseq = df["hseq"].fillna("").astype(str).eq("") if "hseq" in df else pd.Series([], dtype=bool)
    missing_antigen = (
        df["antigen_sequence"].fillna("").astype(str).eq("") if "antigen_sequence" in df else pd.Series([], dtype=bool)
    )
    target_counts = df["target"].fillna("").astype(str).value_counts().to_dict() if "target" in df else {}
    return {
        "results_folder": str(folder),
        "output_csv": str(output_csv),
        "ready_csv": str(ready_csv),
        "target_sequence_table": str(target_sequence_table) if target_sequence_table else "",
        "rows": int(len(df)),
        "ready_rows": int(ready.sum()) if len(df) else 0,
        "not_ready_rows": int(len(df) - ready.sum()) if len(df) else 0,
        "missing_hseq_rows": int(missing_hseq.sum()) if len(df) else 0,
        "missing_antigen_rows": int(missing_antigen.sum()) if len(df) else 0,
        "targets": int(len(target_counts)),
        "rows_by_target": target_counts,
    }


def _sort_lead_rows(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    ranked = df.copy()
    ranked["_stage_b_sort_freq"] = pd.to_numeric(ranked.get("max_freq", pd.Series(0, index=ranked.index)), errors="coerce").fillna(0)
    ranked["_stage_b_sort_rank"] = pd.to_numeric(ranked.get("rank", pd.Series(math.inf, index=ranked.index)), errors="coerce").fillna(math.inf)
    ranked["_stage_b_sort_delphi"] = pd.to_numeric(ranked.get("mean_delphi_score", pd.Series(-1, index=ranked.index)), errors="coerce").fillna(-1)
    ranked = ranked.sort_values(["_stage_b_sort_freq", "_stage_b_sort_rank", "_stage_b_sort_delphi"], ascending=[False, True, False])
    return ranked.drop(columns=["_stage_b_sort_freq", "_stage_b_sort_rank", "_stage_b_sort_delphi"], errors="ignore")


def _sort_candidate_rows(rows: Iterable[dict[str, Any]]) -> pd.DataFrame:
    return _sort_lead_rows(pd.DataFrame(rows))


def _sort_enriched_rows(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out["_stage_b_score_sort"] = pd.to_numeric(out.get("abforge_structural_selection_score", pd.Series(-1, index=out.index)), errors="coerce").fillna(-1)
    out["_stage_b_freq_sort"] = pd.to_numeric(out.get("max_freq", pd.Series(0, index=out.index)), errors="coerce").fillna(0)
    out = out.sort_values(["_stage_b_score_sort", "_stage_b_freq_sort"], ascending=[False, False])
    return out.drop(columns=["_stage_b_score_sort", "_stage_b_freq_sort"], errors="ignore")


def _ensure_columns(df: pd.DataFrame, first_columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in first_columns:
        if col not in out.columns:
            out[col] = ""
    ordered = [col for col in first_columns if col in out.columns]
    ordered.extend([col for col in out.columns if col not in ordered])
    return out[ordered]


def first_existing_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower_to_original = {str(col).lower(): str(col) for col in df.columns}
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
        lowered = candidate.lower()
        if lowered in lower_to_original:
            return lower_to_original[lowered]
    return None


def first_text(row: Any, columns: list[str]) -> str:
    for col in columns:
        if col in row:
            value = clean_text(row.get(col))
            if value:
                return value
    lower_to_original = {str(key).lower(): key for key in getattr(row, "index", [])}
    for col in columns:
        original = lower_to_original.get(col.lower())
        if original is None:
            continue
        value = clean_text(row.get(original))
        if value:
            return value
    return ""


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    if text.lower() in {"nan", "none", "null", "na", "n/a"}:
        return ""
    return text


def clean_sequence(value: Any) -> str:
    text = clean_text(value).upper()
    text = re.sub(r"\s+", "", text)
    text = text.replace("-", "")
    if not text:
        return ""
    return "".join(ch for ch in text if ch.isalpha() or ch == "*").rstrip("*")


def invalid_sequence_letters(sequence: str) -> set[str]:
    return {letter for letter in sequence if letter not in AA_ALPHABET}


def stage_b_issues(hseq: str, antigen_sequence: str) -> list[str]:
    issues: list[str] = []
    if not hseq:
        issues.append("missing_hseq")
    elif len(hseq) < 80:
        issues.append("hseq_too_short")
    invalid_h = invalid_sequence_letters(hseq)
    if invalid_h:
        issues.append(f"invalid_hseq_letters:{''.join(sorted(invalid_h))}")

    if not antigen_sequence:
        issues.append("missing_antigen_sequence")
    elif len(antigen_sequence) < 20:
        issues.append("antigen_too_short")
    invalid_ag = invalid_sequence_letters(antigen_sequence)
    if invalid_ag:
        issues.append(f"invalid_antigen_letters:{''.join(sorted(invalid_ag))}")
    return issues


def infer_antibody_format(row: Any, *, lseq: str) -> str:
    text = " ".join(clean_text(row.get(col)) for col in ["library", "library_type", "antibody_format"] if col in row).lower()
    if "scfv" in text:
        return "scfv"
    if "fab" in text:
        return "fab"
    if "vhh" in text or "nanobody" in text:
        return "vhh"
    return "fab" if lseq else "vhh"


def stage_b_example_id(target: str, clone_id: str, hseq: str, lseq: str, antigen_sequence: str, source_row: int) -> str:
    digest = hashlib.sha1(f"{target}|{clone_id}|{hseq}|{lseq}|{antigen_sequence}|{source_row}".encode("utf-8")).hexdigest()[:12]
    return f"ngs_{safe_token(target)}_{source_row:05d}_{digest}"


def target_from_lead_path(path: Path) -> str:
    return path.stem.replace("_final_leads", "")


def safe_token(value: Any) -> str:
    text = clean_text(value).lower()
    text = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    return text or "unknown"


def safe_column_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_").lower() or "column"


def normalize_key(value: Any) -> str:
    return safe_token(value)


def numeric_or_blank(value: Any) -> float | str:
    number = float_or_none(value)
    return "" if number is None else number


def float_or_none(value: Any) -> float | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number


def probability_value(value: Any) -> float | None:
    number = float_or_none(value)
    if number is None:
        return None
    if number > 1.0 and number <= 100.0:
        number = number / 100.0
    return min(max(number, 0.0), 1.0)
