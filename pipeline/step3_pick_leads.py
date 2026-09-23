# pipeline/step3_pick_leads.py
"""
Step 3: Global lead selection across all targets
Includes:
- low-frequency filtering (MIN_FREQ_1 on max_freq)
- negative control analysis (Biotin, hIgG1_Fc, PSR, Streptavidin)
- cross-contamination dedup (VH+VL+CDR3, then CDR3)
"""

import pandas as pd
from pathlib import Path
import re


def _base_dir(cfg) -> Path:
    return Path(cfg["general"].get("base_dir", "."))


def _clean_cdr3(value) -> str:
    value = str(value)
    return value[1:] if value.startswith("C") else value


DEFAULT_FINAL_EXCLUDE_COLUMNS = [
    "prevalent_targets_vh_cdr3",
    "contaminant_vh_cdr3",
    "cross_target_reads",
    "cross_target_samples",
    "library",
    "library_type",
    "CHAIN",
    "aa",
    "FR1",
    "FR2",
    "FR3",
    "FR4",
    "ML_SEQUENCE_OK",
    "HSEQ_SOURCE",
    "LSEQ_SOURCE",
]


def _final_exclude_columns(cfg: dict) -> list[str]:
    configured = cfg.get("pick_leads", {}).get("final_exclude_columns")
    if configured is None:
        return DEFAULT_FINAL_EXCLUDE_COLUMNS
    return [str(col) for col in configured]


def _clean_final_output_columns(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    drop_cols = [col for col in _final_exclude_columns(cfg) if col in df.columns]
    return df.drop(columns=drop_cols, errors="ignore")


def _cleanup_intermediate_lead_files(folder: Path) -> None:
    deleted = 0
    for name in ["leads.xlsx", "all_ranked_leads.xlsx"]:
        path = folder / name
        if not path.exists():
            continue
        try:
            path.unlink()
            deleted += 1
        except Exception as exc:
            print(f"Warning: could not remove {name}: {exc}")
    if deleted:
        print(f"Removed {deleted} intermediate lead workbook(s)")


def _available_key_columns(df: pd.DataFrame, candidates: list[str]) -> bool:
    return all(col in df.columns for col in candidates)


def _lead_key_columns(clones: pd.DataFrame, leads: pd.DataFrame) -> list[str]:
    vhh_key = ["CDR1", "CDR2", "CDR3", "vh_scaffold"]
    fab_key = ["cdr3_aa", "vh_scaffold", "vl_scaffold"]
    heavy_key = ["cdr3_aa", "vh_scaffold"]

    if _available_key_columns(clones, fab_key) and _available_key_columns(leads, fab_key):
        return fab_key
    if _available_key_columns(clones, vhh_key) and _available_key_columns(leads, vhh_key):
        return vhh_key
    return heavy_key


def _key_set(df: pd.DataFrame, columns: list[str]) -> set[tuple]:
    keyed = df[columns].fillna("").astype(str)
    return set(map(tuple, keyed.itertuples(index=False, name=None)))


def _row_keys(df: pd.DataFrame, columns: list[str]) -> list[tuple]:
    keyed = df[columns].fillna("").astype(str)
    return list(map(tuple, keyed.itertuples(index=False, name=None)))


def run_pick_leads(cfg, folder: Path):
    folder = Path(folder)

    c = cfg["pick_leads"]

    n_leads = c["n_leads"]
    min_freq_first = c["min_freq_first"]
    min_freq_sum = c["min_freq_sum"]
    min_freq_table = c["min_freq_table"]  # your MIN_FREQ_1
    critical_filtering = c["critical_filtering"]
    priority_concs = c["priority_concentrations"]
    dont_order = set(c["dont_order_antigens"])
    write_intermediate_lead_files = bool(c.get("write_intermediate_lead_files", False))
    cleanup_intermediate_lead_files = bool(c.get("cleanup_intermediate_lead_files", True))

    # Negative control files 
    base_dir = _base_dir(cfg)
    negative_tables = {
        "Biotin-90N_1": base_dir / "data/Biotin-90N_Strep_Biotin_HCDR3.txt",
        "Biotin-90N_2": base_dir / "data/Biotin-90N_Strep_HCDR3.txt",
        "hIgG1_Fc": base_dir / "data/HIS-AVI-hIgG1_Fc_pk1_HCDR3.txt",
        "PSR_reagent": base_dir / "data/PSR_reagent_HCDR3.txt",
        "Streptavidin": base_dir / "data/Streptavidin_beads_HCDR3.txt",
    }

    # Load negative controls
    negative_sets = {}
    for name, path in negative_tables.items():
        try:
            with open(path) as f:
                seqs = {line.strip() for line in f if line.strip()}
            negative_sets[name] = seqs
            print(f"Loaded {name}: {len(seqs)} negative CDR3s")
        except Exception as e:
            print(f"Warning: Could not load {name}: {e}")

    # Load previous antibodies for exact repeat removal
    prev_db_path = cfg["general"]["previous_antibodies_db"]
    old_cdr3 = set()
    if Path(prev_db_path).exists():
        old_ab = pd.read_excel(prev_db_path)
        if "CDR3" in old_ab.columns:
            old_cdr3 = {_clean_cdr3(value) for value in old_ab["CDR3"].dropna()}
        print(f"Loaded {len(old_cdr3)} previous CDR3s for repeat removal")

    # Load all clones files
    clone_files = list(folder.glob("*_clones.csv"))

    if not clone_files:
        print("No _clones.csv files found — Step 3 skipped.")
        return

    all_clones = []
    for f in clone_files:
        target = f.stem.replace("_clones", "")
        if target in dont_order:
            print(f"Skipping dont_order target: {target}")
            continue

        try:
            df = pd.read_csv(f)
            df["target"] = target
            all_clones.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}")

    if not all_clones:
        print("No data loaded — Step 3 skipped.")
        return

    leads = pd.concat(all_clones, ignore_index=True)

    # Exact CDR3 repeat removal
    before = len(leads)
    leads = leads[~leads["cdr3_aa"].astype(str).map(_clean_cdr3).isin(old_cdr3)]
    print(f"Removed {before - len(leads)} exact CDR3 repeats from previous antibodies")

    # Critical filtering
    if critical_filtering:
        critical_cols = [col for col in leads.columns if col.startswith("l_")]
        if critical_cols:
            before = len(leads)
            leads["critical"] = leads[critical_cols].any(axis=1)
            leads = leads[~leads["critical"]]
            print(f"Removed {before - len(leads)} clones with critical liabilities")

    # Frequency columns
    freq_cols = [col for col in leads.columns if col.startswith("freq ")]
    if freq_cols:
        leads["max_freq"] = leads[freq_cols].max(axis=1)
        leads["sum_freq"] = leads[freq_cols].sum(axis=1)

        leads = leads[
            (leads["max_freq"] >= min_freq_first) |
            (leads["sum_freq"] >= min_freq_sum)
        ]

    # === ORIGINAL LOW-FREQUENCY FILTERING ===
    print(f"Before low-freq filter: {len(leads)} clones")
    leads = leads[leads["max_freq"] >= min_freq_table]
    print(f"After low-freq filter (max_freq >= {min_freq_table}): {len(leads)} clones")

    # === ORIGINAL CROSS-CONTAMINATION DEDUP ===
    print(f"Before cross-contamination removal: {len(leads)} clones")

    leads.sort_values("max_freq", ascending=False, inplace=True)

    before = len(leads)
    temp=["vh_scaffold", "cdr3_aa"]

    if 'vl_scaffold' in leads.columns:
        temp=["vh_scaffold","vl_scaffold", "cdr3_aa"]
        
    leads = leads.drop_duplicates(subset=temp, keep="first")
    print(f"After VH+VL+CDR3 dedup: {len(leads)} clones (removed {before - len(leads)})")

    before = len(leads)
    leads = leads.drop_duplicates(subset=["cdr3_aa"], keep="first")
    print(f"After CDR3 dedup: {len(leads)} clones (removed {before - len(leads)})")

    # === ORIGINAL NEGATIVE CONTROL ANALYSIS ===
    leads["negative"] = ""
    for name, seqs in negative_sets.items():
        leads["negative"] = leads[["cdr3_aa", "negative"]].apply(
            lambda r: r["negative"] + name + " " if r["cdr3_aa"] in seqs else r["negative"],
            axis=1
        )

    negative_hits = (leads["negative"] != "").sum()
    print(f"Negative control hits: {negative_hits} clones")

    # Priority concentrations
    if priority_concs:
        priority_samples = [col for col in leads.columns if any(conc in col for conc in priority_concs)]
        if priority_samples:
            leads["priority_freq"] = leads[priority_samples].max(axis=1)
            leads = leads.sort_values("priority_freq", ascending=False)

    # Final selection
    leads = leads.sort_values("max_freq", ascending=False)
    final_leads = leads.head(n_leads)

    if write_intermediate_lead_files:
        out = folder / "leads.xlsx"
        final_leads.to_excel(out, index=False)
        print(f"\nSaved intermediate global leads: {out.name} ({len(final_leads)} clones)")

        ranked_out = folder / "all_ranked_leads.xlsx"
        leads.to_excel(ranked_out, index=False)
        print(f"Saved intermediate ranked table: {ranked_out.name} ({len(leads)} clones)")
    elif cleanup_intermediate_lead_files:
        _cleanup_intermediate_lead_files(folder)


    ## Per-target final leads ##
    print("\nGenerating per-target final leads...")
    leads_for_targets = final_leads.copy()
    print(f"Using {len(leads_for_targets)} selected global leads in memory")

    # Create by_protein directory
    by_protein = folder / "by_protein"
    by_protein.mkdir(exist_ok=True)

    # Group leads by target for fast lookup
    leads_by_target = leads_for_targets.groupby("target")["cdr3_aa"].apply(set).to_dict()

    # Process each target
    for target, valid_cdr3_set in leads_by_target.items():
        print(f"\nProcessing target: {target} ({len(valid_cdr3_set)} leads)")
        target= target.replace(".csv", "")
        clones_file = folder / f"{target}_clones.csv"
        if not clones_file.exists():
            print(f"  Warning: {clones_file.name} not found — skipping")
            continue

        df = pd.read_csv(clones_file)

        # Keep only CDR3s in leads (your contamination/low-freq removal)
        before = len(df)
        df = df[df["cdr3_aa"].isin(valid_cdr3_set)]
        print(f"  Kept {len(df)} clones present in leads (removed {before - len(df)})")

        # Remove UNK scaffolds
        before = len(df)
        df = df[df["vh_scaffold"] != "UNK"]
        if 'vl_scaffold' in df.columns:
            df = df[df["vl_scaffold"] != "UNK"]
        print(f"  Removed UNK scaffolds: {len(df)} remaining (removed {before - len(df)})")

        # === Charge calculation in CDR3 ===
        df["pos_charge"] = df["cdr3_aa"].apply(lambda x: len(re.findall(r"[KRH]", x)))
        df["neg_charge"] = df["cdr3_aa"].apply(lambda x: len(re.findall(r"[ED]", x)))
        df["neg_pos"] = df["neg_charge"] - df["pos_charge"]

        # === Concentration ratios ===
        freq_cols = [col for col in df.columns if col.startswith("freq ")]
        freq_2 = [col for col in freq_cols if "2uM" in col]
        freq_4 = [col for col in freq_cols if "4nM" in col]
        freq_20 = [col for col in freq_cols if "20nM" in col]
        freq_100 = [col for col in freq_cols if "100nM" in col]

        if freq_2 and freq_20:
            df["ratio_2uM_20nM"] = df[freq_2[0]] / (df[freq_20[0]] + 1e-8)  # avoid div by zero
        if freq_2 and freq_4:
            df["ratio_2uM_4nM"] = df[freq_2[0]] / (df[freq_4[0]] + 1e-8)  # avoid div by zero
        if freq_2 and freq_100:
            df["ratio_2uM_100nM"] = df[freq_2[0]] / (df[freq_100[0]] + 1e-8)  # avoid div by zero

        if freq_4 and freq_20:
            df["ratio_4nM_20nM"] = df[freq_4[0]] / (df[freq_20[0]] + 1e-8)  # avoid div by zero
        if freq_20 and freq_100:
            df["ratio_20nM_100nM"] = df[freq_20[0]] / (df[freq_100[0]] + 1e-8)
        if freq_4 and freq_100:
            df["ratio_4nM_100nM"] = df[freq_4[0]] / (df[freq_100[0]] + 1e-8)

        # Save final leads for this target
        out_file = by_protein / f"{target}_final_leads.xlsx"
        _clean_final_output_columns(df, cfg).to_excel(out_file, index=False)
        print(f"  Saved: {out_file.name} ({len(df)} clones)")

    print("\nStep 3 complete! Per-target final leads saved in by_protein/")


    ##### update *clones.csv.gz files to only have final leads ##

    leads = leads_for_targets
    print(f"Using {len(leads)} best leads for clone LEAD flags")

    # Process each clones file
    clone_files = list(folder.glob("*_clones.csv"))

    for f in clone_files:
        target = f.stem.replace("_clones", "")
        print(f"\nProcessing target: {target}")

        df = pd.read_csv(f)

        if len(df) == 0:
            print("  No clones — skipping")
            continue

        key_columns = _lead_key_columns(df, leads)
        lead_keys = _key_set(leads, key_columns)
        df_keys = _row_keys(df, key_columns)
        print(f"  Matching leads by: {', '.join(key_columns)}")
        
        # Add LEAD flag
        df["LEAD"] = [key in lead_keys for key in df_keys]

        leads_count = df["LEAD"].sum()
        print(f"  Flagged {leads_count} leads out of {len(df)} clones")


        # === Concentration ratios ===
        freq_cols = [col for col in df.columns if col.startswith("freq ")]
        freq_2 = [col for col in freq_cols if "2uM" in col]
        freq_4 = [col for col in freq_cols if "4nM" in col]
        freq_20 = [col for col in freq_cols if "20nM" in col]
        freq_100 = [col for col in freq_cols if "100nM" in col]

        if freq_4 and freq_20:
            df["ratio_4nM_20nM"] = df[freq_4[0]] / (df[freq_20[0]] + 1e-8)  # avoid div by zero
        if freq_2 and freq_20:
            df["ratio_2uM_20nM"] = df[freq_2[0]] / (df[freq_20[0]] + 1e-8)  # avoid div by zero
        if freq_2 and freq_4:
            df["ratio_2uM_4nM"] = df[freq_2[0]] / (df[freq_4[0]] + 1e-8)  # avoid div by zero
        if freq_2 and freq_100:
            df["ratio_2uM_100nM"] = df[freq_2[0]] / (df[freq_100[0]] + 1e-8)  # avoid div by zero
        if freq_20 and freq_100:
            df["ratio_20nM_100nM"] = df[freq_20[0]] / (df[freq_100[0]] + 1e-8)
        if freq_4 and freq_100:
            df["ratio_4nM_100nM"] = df[freq_4[0]] / (df[freq_100[0]] + 1e-8)


        # Save back (overwrite with LEAD column)
        df.to_csv(f, index=False)

    if cleanup_intermediate_lead_files:
        _cleanup_intermediate_lead_files(folder)

    print("\nStep 3 complete! LEAD = True added to all _clones.csv files")
