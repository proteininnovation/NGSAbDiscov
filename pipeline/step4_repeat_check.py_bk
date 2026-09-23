# pipeline/step4_repeat_check.py
"""
Step 4: Check for repeats against previously ordered antibodies
Adds repeat annotations + 'is_repeat' column to final leads in by_protein/
"""

import pandas as pd
from pathlib import Path
import re
from Bio.Align import substitution_matrices
from concurrent.futures import ProcessPoolExecutor
import shutil
import warnings
warnings.filterwarnings("ignore")

BLOSUM62 = substitution_matrices.load("BLOSUM62")
BLOSUM62_THRESHOLD = 0.8
MAX_CDR3_LEN_FOR_BLOSUM = 25
N_WORKERS = 8


def load_previous_antibodies(antibody_path: Path) -> pd.DataFrame:
    print("Loading previously ordered antibodies...")
    ordered = pd.read_excel(antibody_path)
    required = ["BARCODE", "CDR3", "heavy", "light", "name", "FACS", "BLI"]
    ordered = ordered[required].dropna(subset=["CDR3", "heavy", "light"]).copy()
    ordered = ordered[~ordered["CDR3"].str.contains(":", na=False)]

    ordered["CDR3_clean"] = ordered["CDR3"].str.replace("_", "")
    ordered["Self_BLOSUM"] = [
        sum(BLOSUM62[a, a] for a in seq if a in BLOSUM62.alphabet)
        for seq in ordered["CDR3_clean"]
    ]
    print(f"Loaded {len(ordered)} previous antibodies.")
    return ordered


def generate_cluster_files(main_dir: Path, targets: list[str]):
    cluster_dir = main_dir / "cluster"
    cluster_dir.mkdir(exist_ok=True)

    for target in targets:
        full_file = main_dir / "by_protein" / f"{target}_final_leads.xlsx"
        if not full_file.exists():
            continue
        print(full_file)
        df = pd.read_excel(full_file)

        cluster_df = pd.DataFrame({
            "CDR3": "C" + df["cdr3_aa"],
            "heavy": "V" + df["vh_scaffold"],
            "light": "V" + df["vl_scaffold"].replace("K4-1_C", "K4-1"),
            "Aff3_Combined": 1,
            "CDR3_ClustNum": range(1, len(df) + 1)
        })

        out = cluster_dir / f"{target}.xlsx"
        cluster_df[["CDR3", "heavy", "light", "key" not in cluster_df.columns and "Aff3_Combined", "CDR3_ClustNum"]].to_excel(
            out, sheet_name="cluster", index=False
        )


def check_repeats_for_target(target: str, main_dir: Path, ordered: pd.DataFrame):
    cluster_file = main_dir / "cluster" / f"{target}.xlsx"
    if not cluster_file.exists():
        return

    df = pd.read_excel(cluster_file, sheet_name="cluster")
    df = df[df["CDR3"].str.len() < MAX_CDR3_LEN_FOR_BLOSUM]

    results = []
    for _, row in df.iterrows():
        cdr3_new = row["CDR3"][1:]
        heavy_new = row["heavy"][1:]
        light_new = row["light"][1:]

        candidates = ordered[ordered["CDR3_clean"].str.len() == len(cdr3_new)]

        row_result = {
            "Cluster_match": "", "Cluster_names": "",
            "CDR3_match": "", "CDR3_names": "",
            "HC_match": "", "HC_names": "",
            "Ab_match": "", "Ab_names": "", "Ab_FACS": "", "Ab_BLI": ""
        }

        if len(candidates) > 0:
            scores = [sum(BLOSUM62[a, b] for a, b in zip(cdr3_new, seq) if a in BLOSUM62.alphabet and b in BLOSUM62.alphabet)
                      for seq in candidates["CDR3_clean"]]
            candidates = candidates.copy()
            candidates["norm_score"] = scores / candidates["Self_BLOSUM"]

            clustered = candidates[candidates["norm_score"] > BLOSUM62_THRESHOLD]
            row_result["Cluster_match"] = ";".join(clustered["BARCODE"]) if len(clustered) > 0 else ""
            row_result["Cluster_names"] = ";".join(clustered["name"]) if len(clustered) > 0 else ""

            exact = candidates[candidates["CDR3_clean"] == cdr3_new]
            if len(exact) > 0:
                row_result["CDR3_match"] = ";".join(exact["BARCODE"])
                row_result["CDR3_names"] = ";".join(exact["name"])

                hc = exact[exact["Heavy"] == heavy_new]
                if len(hc) > 0:
                    row_result["HC_match"] = ";".join(hc["BARCODE"])
                    row_result["HC_names"] = ";".join(hc["name"])

                    ab = hc[hc["Light"] == light_new]
                    if len(ab) > 0:
                        row_result["Ab_match"] = ";".join(ab["BARCODE"])
                        row_result["Ab_names"] = ";".join(ab["name"])
                        row_result["Ab_FACS"] = ";".join(ab["FACS"].astype(str))
                        row_result["Ab_BLI"] = ";".join(ab["BLI"].astype(str))

        results.append(row_result)

    result_df = pd.concat([df.reset_index(drop=True), pd.DataFrame(results)], axis=1)
    result_df["CDR_key"] = result_df["CDR3"] + ";" + result_df["heavy"] + ";" + result_df["light"]
    result_df = result_df.drop_duplicates("CDR_key")
    result_df.drop(columns=["key", "CDR_key"], errors="ignore", inplace=True)

    (main_dir / "cluster_repeat").mkdir(exist_ok=True)
    result_df.to_csv(main_dir / "cluster_repeat" / f"{target}.csv", index=False)


def enrich_with_repeats(main_dir: Path, targets: list[str]):
    print("\nEnriching leads with repeat info + adding 'is_repeat' column...")
    for target in targets:
        repeat_file = main_dir / "cluster_repeat" / f"{target}.csv"
        if not repeat_file.exists():
            continue

        repeat_df = pd.read_csv(repeat_file).fillna("")
        repeat_df["CDR3"] = repeat_df["CDR3"].astype(str).str[1:]
        repeat_df["heavy"] = repeat_df["heavy"].astype(str).str[1:]
        repeat_df["light"] = repeat_df["light"].astype(str).str[1:]

        full_file = main_dir / "by_protein" / f"{target}_final_leads.xlsx"
        if not full_file.exists():
            continue

        leads_df = pd.read_excel(full_file)

        enriched = leads_df.merge(
            repeat_df,
            left_on=["cdr3_aa", "vh_scaffold", "vl_scaffold"],
            right_on=["CDR3", "heavy", "light"],
            how="left"
        ).fillna("")

        enriched.drop(columns=[c for c in ["CDR3", "heavy", "light", "Aff3_Combined", "CDR3_ClustNum"] if c in enriched.columns],
                      inplace=True)

        # Add is_repeat: True if full antibody match
        enriched["is_repeat"] = enriched["Ab_match"].str.len() > 0

        # Reorder: put is_repeat early
        cols = enriched.columns.tolist()
        if "is_repeat" in cols:
            cols.remove("is_repeat")
            insert_pos = cols.index("max_freq") if "max_freq" in cols else 5
            cols.insert(insert_pos, "is_repeat")
        enriched = enriched[cols]

        # Overwrite
        enriched.to_excel(full_file, index=False)

        # Update dedup
        dedup = enriched.loc[enriched.groupby("cdr3_aa")["max_freq"].idxmax()]
        dedup_file = main_dir / "by_protein/final_leads_dedup" / f"{target}_final_leads_dedup_bytopfreq.xlsx"
        dedup.to_excel(dedup_file, index=False)

        print(f"  → Updated: {target}")


def run_repeat_check(cfg, folder: Path, antibody_list_path: Path | None = None):
    folder = Path(folder)

    if antibody_list_path is None:
        antibody_list_path = Path(cfg["general"]["previous_antibodies_db"])

    ordered = load_previous_antibodies(antibody_list_path)

    # Find targets from by_protein files
    files = list((folder / "by_protein").glob("*_final_leads.xlsx"))
    targets = [f.stem.replace("_final_leads", "") for f in files]

    generate_cluster_files(folder, targets)

    print("Running repeat check in parallel...")
    with ProcessPoolExecutor(max_workers=N_WORKERS) as exec:
        exec.map(lambda t: check_repeats_for_target(t, folder, ordered), targets)
    enrich_with_repeats(folder, targets)

    # Cleanup
    shutil.rmtree(folder / "cluster", ignore_errors=True)
    shutil.rmtree(folder / "cluster_repeat", ignore_errors=True)

    print("\nRepeat check complete! All files updated with 'is_repeat' column.")
