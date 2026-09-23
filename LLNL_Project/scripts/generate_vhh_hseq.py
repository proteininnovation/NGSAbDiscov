#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_INPUT = Path("LLNL_Project/sample_sheets/ipi_vhh_selected.xlsx")
DEFAULT_LIBRARY = Path("data/IPI_VHH_LIB_SEQ.csv")
DEFAULT_OUTPUT_DIR = Path("LLNL_Project/results/vhh_hseq")
DEFAULT_SHEET = "Sheet1"

DERIVED_COLUMNS = [
    "vh_scaffold_normalized",
    "CDR1",
    "CDR2",
    "CDR3",
    "CDR3_for_HSEQ",
    "HSEQ",
    "HSEQ_STATUS",
    "HSEQ_SOURCE",
]


def clean_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "null", "n/a", "na"} else text


def clean_aa(value: Any) -> str:
    text = clean_text(value)
    if not text:
        return ""
    return re.sub(r"[^A-Za-z*]", "", text).upper()


def scaffold_candidates(value: Any) -> list[str]:
    raw = clean_text(value)
    if not raw:
        return []
    compact = re.sub(r"[^A-Za-z0-9]+", "", raw).upper()
    candidates = [raw, raw.upper(), compact]

    match = re.fullmatch(r"(?:IPI)?VHH(?:IPI)?(\d+)", compact)
    if match:
        number = match.group(1)
        candidates.extend([f"VHH-{number}", f"IPI-VHH-{number}", f"VHH{number}"])

    if compact.startswith("VHH") and compact[3:].isdigit():
        number = compact[3:]
        candidates.extend([f"VHH-{number}", f"IPI-VHH-{number}"])

    if compact.startswith("IPIVHH") and compact[6:].isdigit():
        number = compact[6:]
        candidates.extend([f"IPI-VHH-{number}", f"VHH-{number}"])

    return list(dict.fromkeys(candidates))


def load_vhh_library(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"VHH scaffold library not found: {path}")

    library = pd.read_csv(path, encoding="utf-8-sig").fillna("")
    required = {"ID", "FR1", "CDR1", "FR2", "CDR2", "FR3", "FR4"}
    missing = required.difference(library.columns)
    if missing:
        raise ValueError(f"VHH scaffold library is missing columns: {', '.join(sorted(missing))}")

    lookup: dict[str, dict[str, str]] = {}
    for _, row in library.iterrows():
        record = {column: clean_aa(row[column]) for column in ["FR1", "CDR1", "FR2", "CDR2", "FR3", "FR4"]}
        record["ID"] = clean_text(row["ID"])
        for candidate in scaffold_candidates(row["ID"]):
            lookup[candidate.upper()] = record
            lookup[re.sub(r"[^A-Za-z0-9]+", "", candidate).upper()] = record
    return lookup


def lookup_scaffold(library: dict[str, dict[str, str]], value: Any) -> dict[str, str] | None:
    for candidate in scaffold_candidates(value):
        for key in [candidate.upper(), re.sub(r"[^A-Za-z0-9]+", "", candidate).upper()]:
            if key in library:
                return library[key]
    return None


def parse_cdrs(row: pd.Series, scaffold: dict[str, str] | None) -> tuple[str, str, str, str, str]:
    cdr_h3 = clean_text(row.get("CDR H3", ""))
    cdr_h1 = clean_aa(row.get("CDR H1", ""))
    cdr_h2 = clean_aa(row.get("CDR H2", ""))

    if not cdr_h3:
        return "", "", "", "", "missing_cdr_h3"

    if ":" in cdr_h3:
        parts = [clean_aa(part) for part in cdr_h3.split(":")]
        if len(parts) != 3 or not all(parts):
            return "", "", "", "", "invalid_cdr_h3_parts"
        cdr1, cdr2, cdr3 = parts
        return cdr1, cdr2, cdr3, "cdr_h3_triplet", "ok"

    cdr3 = clean_aa(cdr_h3)
    if not cdr3:
        return "", "", "", "", "missing_cdr3"
    if not cdr_h1 and scaffold:
        cdr_h1 = clean_aa(scaffold.get("CDR1", ""))
    if not cdr_h2 and scaffold:
        cdr_h2 = clean_aa(scaffold.get("CDR2", ""))
    source = "cdr_h3_only_with_source_cdrh1_cdrh2" if cdr_h1 and cdr_h2 else "cdr_h3_only"
    return cdr_h1, cdr_h2, cdr3, source, "ok"


def cdr3_for_hseq(fr3: str, cdr3: str) -> str:
    if fr3.endswith("C") and cdr3.startswith("C"):
        return cdr3[1:]
    return cdr3


def build_hseq(scaffold: dict[str, str], cdr1: str, cdr2: str, cdr3: str) -> tuple[str, str]:
    pieces = {
        "FR1": clean_aa(scaffold.get("FR1", "")),
        "FR2": clean_aa(scaffold.get("FR2", "")),
        "FR3": clean_aa(scaffold.get("FR3", "")),
        "FR4": clean_aa(scaffold.get("FR4", "")),
        "CDR1": cdr1,
        "CDR2": cdr2,
        "CDR3": cdr3,
    }
    if any("*" in value for value in pieces.values()):
        return "", "stop_codon"
    missing = [name for name, value in pieces.items() if not value]
    if missing:
        return "", f"missing_{'_'.join(missing).lower()}"

    inserted_cdr3 = cdr3_for_hseq(pieces["FR3"], cdr3)
    hseq = "".join([pieces["FR1"], cdr1, pieces["FR2"], cdr2, pieces["FR3"], inserted_cdr3, pieces["FR4"]])
    if len(hseq) < 105 or "WGQGTLVTVSS" not in hseq:
        return "", "incomplete_hseq"
    return hseq, "ok"


def generate(input_path: Path, library_path: Path, sheet_name: str) -> tuple[pd.DataFrame, dict]:
    df = pd.read_excel(input_path, sheet_name=sheet_name)
    library = load_vhh_library(library_path)

    derived_rows: list[dict[str, str]] = []
    for _, row in df.iterrows():
        scaffold = lookup_scaffold(library, row.get("VH Scaffold", ""))
        cdr1, cdr2, cdr3, cdr_source, parse_status = parse_cdrs(row, scaffold)

        hseq = ""
        hseq_status = parse_status
        cdr3_insert = ""
        normalized_scaffold = scaffold["ID"] if scaffold else ""
        if not scaffold:
            hseq_status = "unknown_scaffold"
        elif parse_status == "ok":
            cdr3_insert = cdr3_for_hseq(clean_aa(scaffold.get("FR3", "")), cdr3)
            hseq, hseq_status = build_hseq(scaffold, cdr1, cdr2, cdr3)

        derived_rows.append(
            {
                "vh_scaffold_normalized": normalized_scaffold,
                "CDR1": cdr1,
                "CDR2": cdr2,
                "CDR3": cdr3,
                "CDR3_for_HSEQ": cdr3_insert,
                "HSEQ": hseq,
                "HSEQ_STATUS": hseq_status,
                "HSEQ_SOURCE": cdr_source if hseq_status == "ok" else "",
            }
        )

    derived = pd.DataFrame(derived_rows)
    output = pd.concat([df, derived], axis=1)

    cdr_h3_text = df.get("CDR H3", pd.Series(dtype=str)).fillna("").astype(str)
    summary = {
        "input_file": str(input_path),
        "sheet": sheet_name,
        "library_file": str(library_path),
        "rows": int(len(output)),
        "hseq_ok_rows": int(output["HSEQ_STATUS"].eq("ok").sum()),
        "hseq_missing_rows": int(output["HSEQ_STATUS"].ne("ok").sum()),
        "cdr_h3_triplet_rows": int(cdr_h3_text.str.contains(":", regex=False).sum()),
        "cdr_h3_only_rows": int((~cdr_h3_text.str.contains(":", regex=False)).sum()),
        "status_counts": output["HSEQ_STATUS"].value_counts(dropna=False).to_dict(),
        "scaffold_counts": output["vh_scaffold_normalized"].value_counts(dropna=False).to_dict(),
    }
    return output, summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate full VHH HSEQ values from selected IPI VHH spreadsheet rows.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--library", type=Path, default=DEFAULT_LIBRARY)
    parser.add_argument("--sheet", default=DEFAULT_SHEET)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    output, summary = generate(args.input, args.library, args.sheet)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = args.output_dir / f"{args.input.stem}_hseq.csv"
    summary_path = args.output_dir / f"{args.input.stem}_hseq_summary.json"
    output.to_csv(csv_path, index=False)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    print(f"Rows processed: {summary['rows']}")
    print(f"HSEQ ok rows: {summary['hseq_ok_rows']}")
    print(f"HSEQ missing rows: {summary['hseq_missing_rows']}")
    print(f"CDR H3 triplet rows: {summary['cdr_h3_triplet_rows']}")
    print(f"CDR H3-only rows: {summary['cdr_h3_only_rows']}")
    print(f"Output CSV: {csv_path}")
    print(f"Summary JSON: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
