# pipeline/step4_repeat_check.py
"""
Step 4: Check for repeats against previously ordered antibodies
Adds repeat annotations + 'is_repeat' column to final leads in by_protein/
Updated per your request:
- In previous DB: only remove first 'V' from heavy/light if it starts with 'V'
- CDR3: only remove first 'C' if it starts with 'C'
- Cluster sequences remain plain (cdr3_aa, vh_scaffold, vl_scaffold)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio.Align import substitution_matrices
import shutil
import warnings
warnings.filterwarnings("ignore")

BLOSUM62 = substitution_matrices.load("BLOSUM62")
BLOSUM_ALPHABET = set(BLOSUM62.alphabet)
BLOSUM62_THRESHOLD = 0.8
MAX_CDR3_LEN_FOR_BLOSUM = 30
REPEAT_RESULT_COLUMNS = [
    "Cluster_match", "Cluster_names", "CDR3_match", "CDR3_names",
    "HC_match", "HC_names", "Ab_match", "Ab_names", "Ab_FACS", "Ab_BLI",
]


def _build_blosum_lookup() -> np.ndarray:
    lookup = np.zeros((256, 256), dtype=np.int16)
    for aa in BLOSUM_ALPHABET:
        for bb in BLOSUM_ALPHABET:
            lookup[ord(aa), ord(bb)] = int(BLOSUM62[aa, bb])
    return lookup


BLOSUM_LOOKUP = _build_blosum_lookup()


def _as_text(value) -> str:
    if pd.isna(value):
        return ""
    return str(value)


def _join_values(series: pd.Series) -> str:
    return ";".join(series.fillna("").astype(str))


def _sequence_codes(seq: str) -> np.ndarray:
    return np.fromiter(
        (ord(char) if ord(char) < 256 else 0 for char in seq),
        dtype=np.uint8,
        count=len(seq),
    )


def _build_repeat_index(ordered: pd.DataFrame) -> dict:
    indexed = ordered.copy()
    indexed["CDR3_clean"] = indexed["CDR3_clean"].fillna("").astype(str)
    indexed["cdr3_len"] = indexed["CDR3_clean"].str.len()
    by_len = {}
    for length, group in indexed.groupby("cdr3_len", dropna=False):
        group = group.reset_index(drop=True)
        if int(length) > 0:
            seq_codes = np.vstack([_sequence_codes(seq) for seq in group["CDR3_clean"]])
        else:
            seq_codes = np.empty((len(group), 0), dtype=np.uint8)
        by_len[int(length)] = {
            "df": group,
            "seq_codes": seq_codes,
            "self_scores": group["Self_BLOSUM"].fillna(0).to_numpy(dtype=float),
        }
    return {
        "by_len": by_len,
        "by_cdr3": {
            str(cdr3): group.reset_index(drop=True)
            for cdr3, group in indexed.groupby("CDR3_clean", dropna=False)
        },
    }


def load_previous_antibodies(antibody_list_path: Path) -> pd.DataFrame:
    print("Loading previously ordered antibodies...")
    ordered = pd.read_excel(antibody_list_path)
    required = ["BARCODE", "CDR3", "heavy"]
    missing = [col for col in required if col not in ordered.columns]
    if missing:
        raise ValueError(f"Previous antibodies DB missing required columns: {missing}")

    optional = [col for col in ["light", "name", "FACS", "BLI"] if col in ordered.columns]
    ordered = ordered[required + optional].dropna(subset=["CDR3", "heavy"]).copy()
    for col in ["light", "name", "FACS", "BLI"]:
        if col not in ordered.columns:
            ordered[col] = ""
    ordered = ordered[~ordered["CDR3"].str.contains(":", na=False)]

    # CDR3: remove first 'C' only if starts with 'C'
    ordered["CDR3_clean"] = ordered["CDR3"].apply(lambda x: x[1:] if str(x).startswith('C') else str(x))
    ordered["CDR3_clean"] = ordered["CDR3_clean"].str.replace("_", "")

    # Heavy/light: remove first 'V' only if starts with 'V'
    ordered["heavy_clean"] = ordered["heavy"].apply(lambda x: x[1:] if str(x).startswith('V') else str(x))
    ordered["light_clean"] = ordered["light"].apply(lambda x: x[1:] if str(x).startswith('V') else str(x))
    ordered["light_clean"] = ordered["light_clean"].str.replace("4-1_C", "4-1")

    ordered["Self_BLOSUM"] = [
        sum(BLOSUM62[a, a] for a in seq if a in BLOSUM_ALPHABET)
        for seq in ordered["CDR3_clean"]
    ]

    print(f"Loaded and normalized {len(ordered)} previous antibodies.")
    return ordered


def generate_cluster_files(main_dir: Path, targets: list[str]):
    cluster_dir = main_dir / "cluster"
    cluster_dir.mkdir(exist_ok=True)

    for target in targets:
        full_file = main_dir / "by_protein" / f"{target}_final_leads.xlsx"
        if not full_file.exists():
            continue

        df = pd.read_excel(full_file)

        has_light_chain = "vl_scaffold" in df.columns

        cluster_data = {
            "CDR3": df["cdr3_aa"].fillna("").astype(str),
            "heavy": df["vh_scaffold"].fillna("").astype(str),
            "light": df["vl_scaffold"].fillna("").astype(str) if has_light_chain else "",
            "Aff3_Combined": 1,
            "CDR3_ClustNum": range(1, len(df) + 1)
        }
        for col in ["CDR1", "CDR2"]:
            if col in df.columns:
                cluster_data[col] = df[col].fillna("").astype(str)

        cluster_df = pd.DataFrame(cluster_data)
        cluster_cols = [col for col in ["CDR1", "CDR2", "CDR3", "heavy", "light", "Aff3_Combined", "CDR3_ClustNum"] if col in cluster_df.columns]

        out = cluster_dir / f"{target}.xlsx"
        cluster_df[cluster_cols].to_excel(
            out, sheet_name="cluster", index=False
        )


def check_repeats_for_target(target: str, main_dir: Path, repeat_index: dict, repeat_dir: Path):
    cluster_file = main_dir / "cluster" / f"{target}.xlsx"
    if not cluster_file.exists():
        print(f"  → No cluster file for {target} — creating empty repeat CSV")
        empty_df = pd.DataFrame(columns=["CDR1", "CDR2", "CDR3", "heavy", "light", "Aff3_Combined", "CDR3_ClustNum", *REPEAT_RESULT_COLUMNS])
        empty_df.to_csv(repeat_dir / f"{target}.csv", index=False)
        return

    df = pd.read_excel(cluster_file, sheet_name="cluster")
    for col in ["CDR1", "CDR2", "CDR3", "heavy", "light"]:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].fillna("").astype(str)

    print(f"  → Processing {target}: {len(df)} clones for repeats")

    results = []
    for row in df.itertuples(index=False):
        cdr3_new = _as_text(getattr(row, "CDR3", ""))
        heavy_new = _as_text(getattr(row, "heavy", ""))
        light_new = _as_text(getattr(row, "light", ""))

        candidates = repeat_index["by_len"].get(len(cdr3_new))

        row_result = {
            "Cluster_match": "", "Cluster_names": "",
            "CDR3_match": "", "CDR3_names": "",
            "HC_match": "", "HC_names": "",
            "Ab_match": "", "Ab_names": "", "Ab_FACS": "", "Ab_BLI": ""
        }

        if candidates is not None and len(cdr3_new) <= MAX_CDR3_LEN_FOR_BLOSUM:
            self_scores = candidates["self_scores"]
            valid_self = self_scores > 0
            if valid_self.any():
                query_codes = _sequence_codes(cdr3_new)
                scores = BLOSUM_LOOKUP[query_codes[None, :], candidates["seq_codes"]].sum(axis=1)
                matched = valid_self & ((scores / self_scores) > BLOSUM62_THRESHOLD)
                if matched.any():
                    clustered = candidates["df"].iloc[np.flatnonzero(matched)]
                    row_result["Cluster_match"] = _join_values(clustered["BARCODE"])
                    row_result["Cluster_names"] = _join_values(clustered["name"])

        exact = repeat_index["by_cdr3"].get(cdr3_new)
        if exact is not None and not exact.empty:
            row_result["CDR3_match"] = _join_values(exact["BARCODE"])
            row_result["CDR3_names"] = _join_values(exact["name"])

            hc = exact[exact["heavy_clean"] == heavy_new]
            if len(hc) > 0:
                row_result["HC_match"] = _join_values(hc["BARCODE"])
                row_result["HC_names"] = _join_values(hc["name"])

                if light_new.strip():
                    ab = hc[hc["light_clean"] == light_new]
                    if len(ab) > 0:
                        row_result["Ab_match"] = _join_values(ab["BARCODE"])
                        row_result["Ab_names"] = _join_values(ab["name"])
                        row_result["Ab_FACS"] = _join_values(ab["FACS"])
                        row_result["Ab_BLI"] = _join_values(ab["BLI"])

        results.append(row_result)

    # Preserve the repeat-result schema when there are no rows. Without explicit
    # columns, DataFrame([]) has no HC_match/Ab_match columns and later
    # enrichment cannot distinguish an empty result from a malformed one.
    repeat_results = pd.DataFrame(results, columns=REPEAT_RESULT_COLUMNS)
    result_df = pd.concat([df.reset_index(drop=True), repeat_results], axis=1)
    for col in ["CDR1", "CDR2", "CDR3", "heavy", "light"]:
        result_df[col] = result_df[col].fillna("").astype(str)
    result_df["CDR_key"] = (
        result_df["CDR1"] + ";" + result_df["CDR2"] + ";" + result_df["CDR3"]
        + ";" + result_df["heavy"] + ";" + result_df["light"]
    )
    result_df = result_df.drop_duplicates("CDR_key")
    result_df.drop(columns=["CDR_key"], errors="ignore", inplace=True)

    csv_path = repeat_dir / f"{target}.csv"
    result_df.to_csv(csv_path, index=False)
    print(f"  → Wrote cluster_repeat CSV: {csv_path.name} ({len(result_df)} rows)")


def enrich_with_repeats(main_dir: Path, targets: list[str]):
    print("\nEnriching leads with repeat info + adding 'is_repeat' column...")
    repeat_dir = main_dir / "cluster_repeat"
    dedup_dir = main_dir / "by_protein" / "final_leads_dedup_bytopfreq"
    dedup_dir.mkdir(parents=True, exist_ok=True)

    for target in targets:
        full_file = main_dir / "by_protein" / f"{target}_final_leads.xlsx"
        if not full_file.exists():
            continue

        leads_df = pd.read_excel(full_file)
        if "CDR3_x" in leads_df.columns and "CDR3" not in leads_df.columns:
            leads_df.rename(columns={"CDR3_x": "CDR3"}, inplace=True)

        stale_cols = [
            "is_repeat", "CDR1_y", "CDR2_y", "CDR3_y",
            "repeat_CDR1", "repeat_CDR2", "repeat_CDR3", "repeat_heavy", "repeat_light",
            "heavy", "light", "Aff3_Combined", "CDR3_ClustNum", *REPEAT_RESULT_COLUMNS,
        ]
        leads_df.drop(columns=[col for col in stale_cols if col in leads_df.columns], inplace=True)
        has_light_chain = "vl_scaffold" in leads_df.columns

        repeat_file = repeat_dir / f"{target}.csv"
        if repeat_file.exists():
            repeat_df = pd.read_csv(repeat_file).fillna("")
            repeat_df.rename(
                columns={
                    "CDR1": "repeat_CDR1",
                    "CDR2": "repeat_CDR2",
                    "CDR3": "repeat_CDR3",
                    "heavy": "repeat_heavy",
                    "light": "repeat_light",
                },
                inplace=True,
            )

            left_keys = ["cdr3_aa", "vh_scaffold"]
            right_keys = ["repeat_CDR3", "repeat_heavy"]
            for left_col, right_col in [("CDR1", "repeat_CDR1"), ("CDR2", "repeat_CDR2")]:
                if left_col in leads_df.columns and right_col in repeat_df.columns:
                    left_keys.append(left_col)
                    right_keys.append(right_col)
            if has_light_chain and "repeat_light" in repeat_df.columns:
                left_keys.append("vl_scaffold")
                right_keys.append("repeat_light")

            enriched = leads_df.merge(
                repeat_df,
                left_on=left_keys,
                right_on=right_keys,
                how="left"
            ).fillna("")

            enriched.drop(
                columns=[
                    c for c in [
                        "repeat_CDR1", "repeat_CDR2", "repeat_CDR3",
                        "repeat_heavy", "repeat_light", "Aff3_Combined", "CDR3_ClustNum"
                    ]
                    if c in enriched.columns
                ],
                inplace=True,
            )
        else:
            enriched = leads_df.copy()
            for col in REPEAT_RESULT_COLUMNS:
                enriched[col] = ""

        repeat_source = "Ab_match" if has_light_chain and "Ab_match" in enriched.columns else "HC_match"
        if repeat_source in enriched.columns:
            repeat_values = enriched[repeat_source]
        else:
            # Empty or legacy repeat files may not contain annotation columns.
            # Use an index-aligned Series so the empty-table case remains valid.
            repeat_values = pd.Series("", index=enriched.index, dtype="object")
        enriched["is_repeat"] = repeat_values.fillna("").astype(str).str.len() > 0

        cols = enriched.columns.tolist()
        if "is_repeat" in cols:
            cols.remove("is_repeat")
            insert_pos = cols.index("max_freq") + 1 if "max_freq" in cols else 5
            cols.insert(insert_pos, "is_repeat")
        enriched = enriched[cols]

        enriched.to_excel(full_file, index=False, engine="openpyxl")
        print(f"  → Updated and saved: {full_file.name} (is_repeat=True for {enriched['is_repeat'].sum()} clones)")

        if "max_freq" in enriched.columns and not enriched.empty:
            dedup = enriched.loc[enriched.groupby("cdr3_aa")["max_freq"].idxmax()]
        else:
            dedup = enriched.drop_duplicates("cdr3_aa")
        dedup_file = dedup_dir / f"{target}_final_leads_dedup_bytopfreq.xlsx"
        dedup.to_excel(dedup_file, index=False, engine="openpyxl")

    print("\nAll files updated with 'is_repeat' column.")


def run_repeat_check(cfg, folder: Path, antibody_list_path: Path | None = None):
    folder = Path(folder)

    if antibody_list_path is None:
        antibody_list_path = Path(cfg["general"]["previous_antibodies_db"])

    ordered = load_previous_antibodies(antibody_list_path)
    repeat_index = _build_repeat_index(ordered)

    files = list((folder / "by_protein").glob("*_final_leads.xlsx"))
    targets = [f.stem.replace("_final_leads", "") for f in files]

    generate_cluster_files(folder, targets)

    repeat_dir = folder / "cluster_repeat"
    repeat_dir.mkdir(exist_ok=True)
    print(f"Created cluster_repeat folder at {repeat_dir}")

    print("Running repeat check sequentially for debugging...")
    for target in targets:
        check_repeats_for_target(target, folder, repeat_index, repeat_dir)

    enrich_with_repeats(folder, targets)

    print("\nRepeat check complete! All files updated with 'is_repeat' column.")
