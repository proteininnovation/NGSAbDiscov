#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_INPUT = Path("LLNL_Project/results/vhh_hseq/ipi_vhh_selected_hseq.xlsx")
DEFAULT_OUTPUT_DIR = Path("LLNL_Project/results/vhh_hseq")
DEFAULT_TARGET_TABLE = Path("/Users/Hoan.Nguyen/Downloads/2026-08-13 04_43_32 PM test - Block 1.csv")
DEFAULT_TARGET_AA_META_TABLE = Path("/Users/Hoan.Nguyen/Downloads/2026-08-13 04_42_41 PM test - Block 1.csv")
DEFAULT_AA_TABLE = Path("/Users/Hoan.Nguyen/Downloads/2026-08-13 04_41_35 PM test - Block 1.csv")

ANTIGEN_COLUMNS = [
    "antigen_aa_raw",
    "antigen_aa",
    "antigen_aa_removed_tags",
    "antigen_aa_removed_length",
    "antigen_aa_clean_status",
    "antigen_aa_source_name",
    "antigen_aa_source_id",
    "antigen_aa_target",
    "antigen_aa_match_status",
    "antigen_aa_candidate_count",
    "antigen_aa_raw_length",
    "antigen_aa_length",
]

BAD_CONSTRUCT_PATTERNS = [
    "C_V5TM",
    "V5TM",
    "EGFP",
    "_FL",
    "-FL",
    "FULL",
    "ANTITNP",
    "TAG_H",
    "TRANSLATION",
]
TAG_PATTERNS = [
    ("StrepTag", re.compile(r"WSHPQFEK")),
    ("AviTag", re.compile(r"GLNDIFEAQKIEWHE")),
    ("FLAG", re.compile(r"DYKDDDDK")),
    ("HA", re.compile(r"YPYDVPDYA")),
    ("HisTag", re.compile(r"H{6,}")),
]
TERMINAL_TAG_WINDOW = 180


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


def terminal_tag_trim(sequence: Any) -> tuple[str, str, str]:
    raw = clean_aa(sequence)
    if not raw:
        return "", "", "missing_raw_sequence"

    search_start = max(0, len(raw) - TERMINAL_TAG_WINDOW)
    matches: list[tuple[int, str]] = []
    for tag_name, pattern in TAG_PATTERNS:
        for match in pattern.finditer(raw):
            if match.start() >= search_start:
                matches.append((match.start(), tag_name))
    if not matches:
        return raw, "", "no_terminal_tag_detected"

    first_tag_pos = min(position for position, _tag_name in matches)
    cut_pos = first_tag_pos
    while cut_pos > 0 and raw[cut_pos - 1] in {"G", "S"}:
        cut_pos -= 1

    cleaned = raw[:cut_pos]
    removed = raw[cut_pos:]
    removed_tags = [tag_name for position, tag_name in sorted(matches) if position >= cut_pos]
    status = "removed_terminal_tags:" + ";".join(dict.fromkeys(removed_tags))
    return cleaned, removed, status


def normalize_key(value: Any) -> str:
    text = clean_text(value).upper()
    text = text.replace("Α", "ALPHA").replace("Β", "BETA")
    text = text.replace("α", "ALPHA").replace("β", "BETA")
    return re.sub(r"[^A-Z0-9]+", "", text)


def row_species(value: Any) -> str:
    text = clean_text(value)
    return text if text else ""


def selected_lookup_keys(row: pd.Series) -> list[tuple[str, str]]:
    target = clean_text(row.get("Target", ""))
    gene = clean_text(row.get("Gene", ""))
    species = row_species(row.get("Species", ""))
    keys: list[tuple[str, str]] = []
    for label, value in [
        ("target", target),
        ("gene_species", f"{gene}_{species}" if gene and species else ""),
        ("gene", gene),
    ]:
        normalized = normalize_key(value)
        if normalized:
            keys.append((label, normalized))
    return list(dict.fromkeys(keys))


def load_candidates(target_aa_meta_path: Path, aa_path: Path) -> pd.DataFrame:
    meta = pd.read_csv(target_aa_meta_path, low_memory=False).fillna("")
    aa = pd.read_csv(aa_path, low_memory=False).fillna("")
    joined = meta.merge(
        aa[["_sync_key", "id", "name", "amino_acids"]],
        on="_sync_key",
        how="left",
        suffixes=("", "_aa"),
    )
    joined["amino_acids"] = joined["amino_acids"].apply(clean_aa)
    joined["aa_len"] = joined["amino_acids"].str.len()
    joined = joined[joined["aa_len"] > 0].copy()
    for column in ["target", "name$", "id", "name", "ipi_code"]:
        if column not in joined.columns:
            joined[column] = ""
        joined[f"norm_{column}"] = joined[column].apply(normalize_key)
    joined["has_antigenic_region"] = joined.get("antigenic_region_residues", "").apply(clean_text).ne("")
    joined["construct_text"] = (
        joined["id"].astype(str) + " " + joined["name$"].astype(str) + " " + joined["name"].astype(str)
    ).str.upper()
    joined["bad_construct"] = joined["construct_text"].apply(
        lambda text: any(pattern in text for pattern in BAD_CONSTRUCT_PATTERNS)
    )
    return joined


def candidate_subset(candidates: pd.DataFrame, row: pd.Series) -> tuple[pd.DataFrame, str]:
    keys = selected_lookup_keys(row)
    for label, key in keys:
        if label == "target":
            mask = candidates["norm_target"].eq(key)
        elif label == "gene_species":
            mask = candidates["norm_target"].eq(key)
        else:
            mask = candidates["norm_target"].str.contains(key, regex=False)
        subset = candidates[mask].copy()
        if not subset.empty:
            return subset, label

    target_key = keys[0][1] if keys else ""
    if target_key:
        name_cols = ["norm_name$", "norm_id", "norm_name"]
        mask = False
        for column in name_cols:
            mask = mask | candidates[column].eq(target_key)
        subset = candidates[mask].copy()
        if not subset.empty:
            return subset, "construct_name"
    return candidates.iloc[0:0].copy(), "unmatched"


def rank_candidates(subset: pd.DataFrame) -> pd.DataFrame:
    ranked = subset.copy()
    ranked["rank_score"] = 0
    ranked.loc[ranked["has_antigenic_region"], "rank_score"] += 500
    ranked.loc[~ranked["bad_construct"], "rank_score"] += 250
    ranked.loc[ranked["bad_construct"], "rank_score"] -= 250
    ranked["rank_length"] = pd.to_numeric(ranked["aa_len"], errors="coerce").fillna(10**9)
    return ranked.sort_values(
        ["rank_score", "rank_length", "name$"],
        ascending=[False, True, True],
        kind="mergesort",
    )


def choose_candidate(candidates: pd.DataFrame, row: pd.Series) -> dict[str, Any]:
    subset, match_by = candidate_subset(candidates, row)
    if subset.empty:
        return {
            "antigen_aa_raw": "",
            "antigen_aa": "",
            "antigen_aa_removed_tags": "",
            "antigen_aa_removed_length": 0,
            "antigen_aa_clean_status": "unmatched",
            "antigen_aa_source_name": "",
            "antigen_aa_source_id": "",
            "antigen_aa_target": clean_text(row.get("Target", "")),
            "antigen_aa_match_status": "unmatched",
            "antigen_aa_candidate_count": 0,
            "antigen_aa_raw_length": 0,
            "antigen_aa_length": 0,
            "_match_by": match_by,
            "_candidate_summary": [],
        }

    ranked = rank_candidates(subset)
    best = ranked.iloc[0]
    status = "unique_match" if len(ranked) == 1 else "multiple_candidates_ranked"
    raw_sequence = best["amino_acids"]
    cleaned_sequence, removed_tags, clean_status = terminal_tag_trim(raw_sequence)
    return {
        "antigen_aa_raw": raw_sequence,
        "antigen_aa": cleaned_sequence,
        "antigen_aa_removed_tags": removed_tags,
        "antigen_aa_removed_length": len(removed_tags),
        "antigen_aa_clean_status": clean_status,
        "antigen_aa_source_name": clean_text(best.get("name$", "")) or clean_text(best.get("name", "")),
        "antigen_aa_source_id": clean_text(best.get("id", "")),
        "antigen_aa_target": clean_text(best.get("target", "")),
        "antigen_aa_match_status": f"{status}:{match_by}",
        "antigen_aa_candidate_count": int(len(ranked)),
        "antigen_aa_raw_length": int(best["aa_len"]),
        "antigen_aa_length": len(cleaned_sequence),
        "_match_by": match_by,
        "_candidate_summary": [
            {
                "name": clean_text(candidate.get("name$", "")) or clean_text(candidate.get("name", "")),
                "id": clean_text(candidate.get("id", "")),
                "target": clean_text(candidate.get("target", "")),
                "aa_len": int(candidate["aa_len"]),
                "has_antigenic_region": bool(candidate["has_antigenic_region"]),
                "bad_construct": bool(candidate["bad_construct"]),
            }
            for _, candidate in ranked.head(10).iterrows()
        ],
    }


def generate(input_path: Path, target_table_path: Path, target_aa_meta_path: Path, aa_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    workbook = pd.read_excel(input_path, sheet_name="VHH_HSEQ")
    candidates = load_candidates(target_aa_meta_path, aa_path)
    _targets = pd.read_csv(target_table_path, low_memory=False)

    added = [choose_candidate(candidates, row) for _, row in workbook.iterrows()]
    added_df = pd.DataFrame([{column: row[column] for column in ANTIGEN_COLUMNS} for row in added])
    output = pd.concat([workbook.drop(columns=[col for col in ANTIGEN_COLUMNS if col in workbook.columns]), added_df], axis=1)

    target_summary_rows = []
    seen_targets = set()
    for source_row, added_row in zip(workbook.to_dict("records"), added, strict=False):
        target = clean_text(source_row.get("Target", ""))
        if target in seen_targets:
            continue
        seen_targets.add(target)
        target_summary_rows.append(
            {
                "Target": target,
                "Gene": clean_text(source_row.get("Gene", "")),
                "Species": clean_text(source_row.get("Species", "")),
                "antigen_aa_source_name": added_row["antigen_aa_source_name"],
                "antigen_aa_source_id": added_row["antigen_aa_source_id"],
                "antigen_aa_target": added_row["antigen_aa_target"],
                "antigen_aa_match_status": added_row["antigen_aa_match_status"],
                "antigen_aa_candidate_count": added_row["antigen_aa_candidate_count"],
                "antigen_aa_raw_length": added_row["antigen_aa_raw_length"],
                "antigen_aa_length": added_row["antigen_aa_length"],
                "antigen_aa_clean_status": added_row["antigen_aa_clean_status"],
                "antigen_aa_removed_length": added_row["antigen_aa_removed_length"],
                "top_candidates": json.dumps(added_row["_candidate_summary"], ensure_ascii=False),
            }
        )
    target_summary = pd.DataFrame(target_summary_rows).sort_values("Target")

    summary = {
        "input_workbook": str(input_path),
        "target_table": str(target_table_path),
        "target_aa_meta_table": str(target_aa_meta_path),
        "aa_table": str(aa_path),
        "rows": int(len(output)),
        "unique_targets": int(output["Target"].nunique()),
        "rows_with_antigen_aa": int(output["antigen_aa"].fillna("").astype(str).str.len().gt(0).sum()),
        "rows_without_antigen_aa": int(output["antigen_aa"].fillna("").astype(str).str.len().eq(0).sum()),
        "rows_with_tags_removed": int(output["antigen_aa_removed_tags"].fillna("").astype(str).str.len().gt(0).sum()),
        "rows_without_tags_removed": int(output["antigen_aa_removed_tags"].fillna("").astype(str).str.len().eq(0).sum()),
        "targets_with_antigen_aa": int(target_summary["antigen_aa_length"].gt(0).sum()),
        "targets_without_antigen_aa": int(target_summary["antigen_aa_length"].eq(0).sum()),
        "targets_with_tags_removed": int(target_summary["antigen_aa_removed_length"].gt(0).sum()),
        "targets_without_tags_removed": int(target_summary["antigen_aa_removed_length"].eq(0).sum()),
        "match_status_counts": output["antigen_aa_match_status"].value_counts(dropna=False).to_dict(),
        "clean_status_counts": output["antigen_aa_clean_status"].value_counts(dropna=False).to_dict(),
        "source_unique_sequences": int(candidates["amino_acids"].nunique()),
    }
    return output, target_summary, summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Add target antigen amino-acid sequences to the VHH selected HSEQ table.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--target-table", type=Path, default=DEFAULT_TARGET_TABLE)
    parser.add_argument("--target-aa-meta-table", type=Path, default=DEFAULT_TARGET_AA_META_TABLE)
    parser.add_argument("--aa-table", type=Path, default=DEFAULT_AA_TABLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--suffix", default="antigen_aa")
    args = parser.parse_args()

    output, target_summary, summary = generate(args.input, args.target_table, args.target_aa_meta_table, args.aa_table)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    base = args.input.stem.replace("_antigen_aa", "")
    output_csv = args.output_dir / f"{base}_{args.suffix}.csv"
    lookup_csv = args.output_dir / f"{base}_{args.suffix}_lookup.csv"
    payload_json = args.output_dir / f"{base}_{args.suffix}_payload.json"
    summary_json = args.output_dir / f"{base}_{args.suffix}_summary.json"

    output.to_csv(output_csv, index=False)
    target_summary.to_csv(lookup_csv, index=False)
    payload = {
        "summary": summary,
        "columns": ANTIGEN_COLUMNS,
        "values": output[ANTIGEN_COLUMNS].fillna("").values.tolist(),
        "lookup_columns": list(target_summary.columns),
        "lookup_values": target_summary.fillna("").values.tolist(),
    }
    payload_json.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    print(f"Rows: {summary['rows']}")
    print(f"Rows with antigen_aa: {summary['rows_with_antigen_aa']}")
    print(f"Rows without antigen_aa: {summary['rows_without_antigen_aa']}")
    print(f"Rows with terminal tags removed: {summary['rows_with_tags_removed']}")
    print(f"Targets with antigen_aa: {summary['targets_with_antigen_aa']}")
    print(f"Targets without antigen_aa: {summary['targets_without_antigen_aa']}")
    print(f"Payload JSON: {payload_json}")
    print(f"Lookup CSV: {lookup_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
