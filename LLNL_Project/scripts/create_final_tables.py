#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pipeline.predict_ml import Generate_FullVHVL


DEFAULT_SELECTED_REPORT = Path("LLNL_Project/results/selected_sample_reports/matched_remote_csvs.csv")
DEFAULT_SOURCE_DIR = Path("LLNL_Project/results/extracted_selected_samples")
DEFAULT_OUTPUT_DIR = Path("LLNL_Project/results/final_tables")
DEFAULT_REPORT_DIR = Path("LLNL_Project/results/final_table_reports")
DEFAULT_LIBRARY_CSV = Path("data/IPI_VLVH_LIB_ALL.csv")
DEFAULT_TAB_ID_DB = Path("data/All_mAb_20260626_FACS_BLI.xlsx")
DEFAULT_MIN_CDR3_LEN = 5
FINAL_COLUMNS = ["cdr3_aa", "vh_scaffold", "vl_scaffold", "HSEQ", "LSEQ", "count", "freq", "TAB_ID"]
K4_1_LSEQ_OVERRIDE = (
    "DIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNNKNYLAWYQQKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQYYSTPLTFGQGTKVEIK"
)


def is_true_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.fillna("").astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def source_path(remote_path: str, source_dir: Path) -> Path:
    return source_dir / remote_path.lstrip("/")


def final_filename(filename: str) -> str:
    if filename.endswith(".csv.gz"):
        return filename[:-3]
    return f"{Path(filename).stem}.csv"


def is_k4_1_scaffold(value: object) -> bool:
    if pd.isna(value):
        return False
    normalized = str(value).strip().upper().replace("_", "-")
    if normalized.startswith("V"):
        normalized = normalized[1:]
    return normalized == "K4-1"


def normalize_cdr3_for_tab_id(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip().upper().replace("_", "")
    return text[1:] if text.startswith("C") else text


def normalize_heavy_for_tab_id(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip().upper()
    return text[1:] if text.startswith("V") else text


def normalize_light_for_tab_id(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip().upper()
    if text.startswith("V"):
        text = text[1:]
    text = text.replace("K4_1_C", "K4-1")
    text = text.replace("K4-1_C", "K4-1")
    text = text.replace("K4_1", "K4-1")
    return text


def tab_id_key(cdr3: object, vh_scaffold: object, vl_scaffold: object) -> str:
    return "|".join(
        [
            normalize_cdr3_for_tab_id(cdr3),
            normalize_heavy_for_tab_id(vh_scaffold),
            normalize_light_for_tab_id(vl_scaffold),
        ]
    )


def normalize_antigen_for_tab_id(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"[^A-Z0-9]+", "", str(value).strip().upper())


def tab_id_antigen_key(antigen: object, block: object) -> str:
    antigen_text = "" if pd.isna(antigen) else str(antigen).strip()
    block_text = "" if pd.isna(block) else str(block).strip()
    if block_text and normalize_antigen_for_tab_id(block_text) not in normalize_antigen_for_tab_id(antigen_text):
        antigen_text = f"{antigen_text}_{block_text}"
    return normalize_antigen_for_tab_id(antigen_text)


def load_tab_id_lookup(db_path: Path) -> dict[str, str]:
    if not db_path.exists():
        raise FileNotFoundError(f"TAB_ID database not found: {db_path}")

    db = pd.read_excel(db_path, usecols=["BARCODE", "CDR3", "heavy", "light", "antigen"])
    db = db.dropna(subset=["BARCODE", "CDR3", "heavy", "light", "antigen"]).copy()
    db["sequence_key"] = [
        tab_id_key(cdr3, heavy, light)
        for cdr3, heavy, light in zip(db["CDR3"], db["heavy"], db["light"], strict=False)
    ]
    db["antigen_key"] = db["antigen"].apply(normalize_antigen_for_tab_id)
    db["BARCODE"] = db["BARCODE"].fillna("").astype(str).str.strip()
    db = db[(db["sequence_key"].str.len() > 2) & db["antigen_key"].ne("") & db["BARCODE"].ne("")]
    db["key"] = db["antigen_key"] + "|" + db["sequence_key"]
    return db.groupby("key")["BARCODE"].apply(lambda values: ";".join(dict.fromkeys(values))).to_dict()


def add_tab_ids(
    df: pd.DataFrame,
    tab_id_lookup: dict[str, str],
    antigen: object,
    block: object,
) -> tuple[pd.DataFrame, str]:
    antigen_key = tab_id_antigen_key(antigen, block)
    keys = [
        f"{antigen_key}|{tab_id_key(cdr3, vh_scaffold, vl_scaffold)}"
        for cdr3, vh_scaffold, vl_scaffold in zip(
            df["cdr3_aa"], df["vh_scaffold"], df["vl_scaffold"], strict=False
        )
    ]
    df = df.copy()
    df["TAB_ID"] = [tab_id_lookup.get(key, "") for key in keys]
    return df, antigen_key


def build_final_table(
    input_path: Path,
    library_csv: Path,
    min_cdr3_len: int,
    tab_id_lookup: dict[str, str],
    antigen: object,
    block: object,
) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(input_path, compression="gzip")
    input_rows = len(df)
    input_count = int(pd.to_numeric(df.get("count", pd.Series(dtype=float)), errors="coerce").fillna(0).sum())
    antigen_key = tab_id_antigen_key(antigen, block)

    df = Generate_FullVHVL(df, library_csv=library_csv)

    if "cdr3_functional" in df.columns:
        df = df[is_true_series(df["cdr3_functional"])].copy()

    for column in ["cdr3_aa", "vh_scaffold", "vl_scaffold", "HSEQ", "LSEQ", "count"]:
        if column not in df.columns:
            df[column] = ""

    text_cols = ["cdr3_aa", "vh_scaffold", "vl_scaffold", "HSEQ", "LSEQ"]
    for column in text_cols:
        df[column] = df[column].fillna("").astype(str).str.strip()

    k4_1_mask = df["vl_scaffold"].apply(is_k4_1_scaffold)
    df.loc[k4_1_mask, "LSEQ"] = K4_1_LSEQ_OVERRIDE
    k4_1_override_rows = int(k4_1_mask.sum())

    cdr3_len_mask = df["cdr3_aa"].str.len() >= min_cdr3_len
    cdr3_len_filtered_rows = int((~cdr3_len_mask).sum())

    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(int)
    df = df[
        (df["count"] > 0)
        & cdr3_len_mask
        & df["cdr3_aa"].ne("")
        & df["vh_scaffold"].ne("")
        & df["vl_scaffold"].ne("")
        & df["vh_scaffold"].ne("UNK")
        & df["vl_scaffold"].ne("UNK")
        & df["HSEQ"].ne("")
        & df["LSEQ"].ne("")
    ].copy()

    total_count = int(df["count"].sum())
    if total_count > 0:
        df = (
            df.groupby(["cdr3_aa", "vh_scaffold", "vl_scaffold", "HSEQ", "LSEQ"], as_index=False)
            .agg({"count": "sum"})
            .sort_values("count", ascending=False)
            .reset_index(drop=True)
        )
        df["freq"] = df["count"] / total_count
        df, antigen_key = add_tab_ids(df, tab_id_lookup, antigen, block)
    else:
        df = pd.DataFrame(columns=FINAL_COLUMNS)

    df = df[FINAL_COLUMNS]
    summary = {
        "input_rows": input_rows,
        "input_count": input_count,
        "final_rows": len(df),
        "final_count": int(df["count"].sum()) if not df.empty else 0,
        "min_cdr3_len": min_cdr3_len,
        "cdr3_len_filtered_rows": cdr3_len_filtered_rows,
        "tab_id_antigen_key": antigen_key,
        "tab_id_matched_rows": int(df["TAB_ID"].fillna("").astype(str).str.len().gt(0).sum()) if not df.empty else 0,
        "tab_id_unmatched_rows": int(df["TAB_ID"].fillna("").astype(str).str.len().eq(0).sum()) if not df.empty else 0,
        "k4_1_override_rows": k4_1_override_rows,
        "missing_hseq_rows": 0,
        "missing_lseq_rows": 0,
    }
    return df, summary


def write_summary(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "antigen",
        "miseq",
        "block",
        "matched_sample_column",
        "matched_sample",
        "source_file",
        "output_file",
        "input_rows",
        "input_count",
        "min_cdr3_len",
        "cdr3_len_filtered_rows",
        "tab_id_antigen_key",
        "final_rows",
        "final_count",
        "tab_id_matched_rows",
        "tab_id_unmatched_rows",
        "k4_1_override_rows",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Create final uncompressed LLNL CSV tables.")
    parser.add_argument("--selected-report", type=Path, default=DEFAULT_SELECTED_REPORT)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--library-csv", type=Path, default=DEFAULT_LIBRARY_CSV)
    parser.add_argument("--tab-id-db", type=Path, default=DEFAULT_TAB_ID_DB)
    parser.add_argument("--min-cdr3-len", type=int, default=DEFAULT_MIN_CDR3_LEN)
    args = parser.parse_args()

    rows = list(csv.DictReader(args.selected_report.open()))
    tab_id_lookup = load_tab_id_lookup(args.tab_id_db)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict] = []
    missing_sources: list[str] = []
    for row in rows:
        source = source_path(row["remote_path"], args.source_dir)
        if not source.exists():
            missing_sources.append(str(source))
            continue

        out_path = args.output_dir / final_filename(row["filename"])
        final_df, summary = build_final_table(
            source,
            args.library_csv,
            args.min_cdr3_len,
            tab_id_lookup,
            row["antigen"],
            row["block"],
        )
        final_df.to_csv(out_path, index=False)

        summary_rows.append(
            {
                "antigen": row["antigen"],
                "miseq": row["miseq"],
                "block": row["block"],
                "matched_sample_column": row["matched_sample_column"],
                "matched_sample": row["matched_sample"],
                "source_file": str(source),
                "output_file": str(out_path),
                "input_rows": summary["input_rows"],
                "input_count": summary["input_count"],
                "min_cdr3_len": summary["min_cdr3_len"],
                "cdr3_len_filtered_rows": summary["cdr3_len_filtered_rows"],
                "tab_id_antigen_key": summary["tab_id_antigen_key"],
                "final_rows": summary["final_rows"],
                "final_count": summary["final_count"],
                "tab_id_matched_rows": summary["tab_id_matched_rows"],
                "tab_id_unmatched_rows": summary["tab_id_unmatched_rows"],
                "k4_1_override_rows": summary["k4_1_override_rows"],
            }
        )

    summary_path = args.report_dir / "final_table_summary.csv"
    write_summary(summary_path, summary_rows)
    if missing_sources:
        missing_path = args.report_dir / "missing_sources.txt"
        missing_path.write_text("\n".join(missing_sources) + "\n")
        raise FileNotFoundError(f"{len(missing_sources)} selected source file(s) were missing")

    print(f"Final CSV tables written: {len(summary_rows)}")
    print(f"Output folder: {args.output_dir}")
    print(f"Summary report: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
