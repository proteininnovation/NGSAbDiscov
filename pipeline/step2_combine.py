# pipeline/step2_combine.py
"""
Step 2: Combine per-sample tables per target
Now includes exact repeat checking from previous antibodies DB
Updated:
- Filenames use .csv (no gzip compression)
- Block number automatically detected:
  - If already in the common prefix (target), preserved
  - Else, if all samples share the same BlockXXX in their sample part (after __), add it as _BlockXXX
- Dynamic greedy cluster column and filename
- Greedy representatives saved
"""

import pandas as pd
from pathlib import Path
from collections import Counter, defaultdict
import re  # Added for block detection
from utilities.clustering import greedy_clustering_by_levenshtein
from utilities.liabilities import annotate_liabilities


def _sample_stem(path: Path) -> str:
    return path.name.removesuffix(".csv.gz")


def _safe_name_part(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip()).strip("_")


def _target_block_from_sample_file(path: Path) -> tuple[str, str]:
    parts = _sample_stem(path).split("__", 4)
    target = parts[0].strip()
    block = parts[1].strip() if len(parts) > 1 else ""
    if not block:
        match = re.search(r"Block\d+(?:_[A-Za-z0-9.-]+)?", target, re.IGNORECASE)
        block = match.group(0) if match else ""
    return target, block


def _target_name_for_group(target: str, block: str) -> str:
    block = _safe_name_part(block)
    if block and block.lower() not in target.lower():
        return f"{target}_{block}"
    return target


def _library_from_sample_file(path: Path, default_library: str) -> str:
    try:
        first_row = pd.read_csv(path, usecols=lambda col: col == "library", nrows=1)
    except Exception:
        return default_library
    if "library" not in first_row.columns or first_row.empty:
        return default_library
    value = str(first_row["library"].iloc[0]).strip()
    return value or default_library


def _combine_setting_for_library(cfg: dict, key: str, library: str, default):
    combine_cfg = cfg.get("combine", {})
    library_cfg = cfg.get("libraries", {}).get(library, {})
    library_type = library_cfg.get("library_type", "")
    by_library = combine_cfg.get(f"{key}_by_library", {})
    by_library_type = combine_cfg.get(f"{key}_by_library_type", {})
    return by_library.get(library) or by_library_type.get(library_type) or combine_cfg.get(key, default)


def _normalize_grouping_columns(pivot_cols: list[str], library_type: str, cdr3_col: str) -> list[str]:
    requested = [str(col).strip() for col in pivot_cols if str(col).strip()]
    library_type = str(library_type or "").lower()
    if library_type == "fab":
        requested = ["CDR3", "vh_scaffold", "vl_scaffold", "anarci_anno"]
    elif library_type == "vhh":
        requested = ["CDR1", "CDR2", "CDR3", "vh_scaffold"]

    if cdr3_col == "CDR3" and "CDR3" in requested:
        requested = [col for col in requested if col != "cdr3_aa"]
    return list(dict.fromkeys(requested))


def _sync_cdr3_aliases(df: pd.DataFrame, preferred_col: str = "CDR3") -> tuple[pd.DataFrame, str]:
    if preferred_col not in df.columns and preferred_col == "CDR3" and "cdr3_aa" in df.columns:
        df["CDR3"] = df["cdr3_aa"]
    if preferred_col not in df.columns and "CDR3" in df.columns:
        preferred_col = "CDR3"
    if preferred_col not in df.columns and "cdr3_aa" in df.columns:
        preferred_col = "cdr3_aa"
    if preferred_col in df.columns:
        df[preferred_col] = df[preferred_col].fillna("").astype(str)
        if preferred_col == "CDR3":
            df["cdr3_aa"] = df["CDR3"]
        elif "CDR3" not in df.columns:
            df["CDR3"] = df[preferred_col]
    return df, preferred_col


def _ensure_clone_cdr3_alias(p: pd.DataFrame, cdr3_col: str) -> pd.DataFrame:
    if "CDR3" not in p.columns and cdr3_col in p.columns:
        insert_at = list(p.columns).index(cdr3_col) + 1
        p.insert(insert_at, "CDR3", p[cdr3_col])
    if "cdr3_aa" not in p.columns and "CDR3" in p.columns:
        insert_at = list(p.columns).index("CDR3") + 1
        p.insert(insert_at, "cdr3_aa", p["CDR3"])
    elif "cdr3_aa" in p.columns and "CDR3" in p.columns:
        p["cdr3_aa"] = p["CDR3"]
    return p


def _cleanup_stale_selection_files(folder: Path, *, keep_greedy: bool) -> None:
    patterns = ["*_leads.csv", "*_leads.csv.gz"]
    if not keep_greedy:
        patterns.append("*_greedy_*.csv")
    deleted = 0
    for pattern in patterns:
        for path in folder.glob(pattern):
            if not path.is_file():
                continue
            try:
                path.unlink()
                deleted += 1
            except Exception as exc:
                print(f"Warning: could not remove stale selection file {path.name}: {exc}")
    if deleted:
        print(f"Removed {deleted} stale lead/selection file(s) from previous runs")

def load_previous_db_bk(db_path: Path) -> tuple[set, set, set]:
    """Load previous antibodies and create lookup sets"""
    if not db_path.exists():
        print("Previous antibodies DB not found — skipping repeat flags")
        return set(), set(), set()

    df = pd.read_excel(db_path)
    
    df = df[["CDR3", "heavy", "light"]].dropna()
    df['CDR3'] = df['CDR3'].str[1:]
    df['heavy'] = df['heavy'].str[1:]
    df['light'] = df['light'].str[1:]

    df.rename(columns={"CDR3": "cdr3_aa", "heavy": "vh_scaffold", "light": "vl_scaffold"}, inplace=True)
    cdr3_set = set(df["cdr3_aa"])
    vh_set = set(df["vh_scaffold"] + "|" + df["cdr3_aa"])
    ab_set = set(df["vh_scaffold"] + "|" + df["vl_scaffold"] + "|" + df["cdr3_aa"])

    print(f"Loaded {len(df)} previous antibodies for exact repeat checking")
    return cdr3_set, vh_set, ab_set

def load_previous_db(db_path: Path) -> tuple[set, set, set, dict, dict, dict, dict]:
    """Load previous antibodies and create lookup sets + dicts for BARCODEs"""
    if not db_path.exists():
        print("Previous antibodies DB not found — skipping repeat flags")
        empty_set = set()
        empty_dict = {}
        return empty_set, empty_set, empty_set, empty_dict, empty_dict, empty_dict, empty_dict

    df = pd.read_excel(db_path)

    required = {"BARCODE", "CDR3", "heavy"}
    missing = sorted(required - set(df.columns))
    if missing:
        print(f"Previous antibodies DB missing columns {missing} — skipping repeat flags")
        empty_set = set()
        empty_dict = {}
        return empty_set, empty_set, empty_set, empty_dict, empty_dict, empty_dict, empty_dict

    optional = [col for col in ["light", "antigen"] if col in df.columns]
    df = df[["BARCODE", "CDR3", "heavy", *optional]].dropna(subset=["CDR3", "heavy"]).copy()
    if "light" not in df.columns:
        df["light"] = ""
    if "antigen" not in df.columns:
        df["antigen"] = ""

    df["CDR3"] = df["CDR3"].astype(str).apply(lambda x: x[1:] if x.startswith("C") else x)
    df["heavy"] = df["heavy"].astype(str).apply(lambda x: x[1:] if x.startswith("V") else x)
    df["light"] = df["light"].astype(str).apply(lambda x: x[1:] if x.startswith("V") else x)

    df.rename(columns={"CDR3": "cdr3_aa", "heavy": "vh_scaffold", "light": "vl_scaffold"}, inplace=True)
    cdr3_set = set(df["cdr3_aa"])
    vh_set = set(df["vh_scaffold"] + "|" + df["cdr3_aa"])
    ab_set = set(
        df.loc[df["vl_scaffold"].astype(str) != "", "vh_scaffold"]
        + "|"
        + df.loc[df["vl_scaffold"].astype(str) != "", "vl_scaffold"]
        + "|"
        + df.loc[df["vl_scaffold"].astype(str) != "", "cdr3_aa"]
    )

    # Dicts for BARCODE lists (sequence → ";" joined BARCODEs)
    cdr3_dict = df.groupby("cdr3_aa")["BARCODE"].apply(lambda x: ";".join(x.astype(str))).to_dict()
    vh_dict = df.groupby(df["vh_scaffold"] + "|" + df["cdr3_aa"])["BARCODE"].apply(lambda x: ";".join(x.astype(str))).to_dict()
    ab_dict = df.groupby(df["vh_scaffold"] + "|" + df["vl_scaffold"] + "|" + df["cdr3_aa"])["BARCODE"].apply(lambda x: ";".join(x.astype(str))).to_dict()
    ab_dict_ag = df.groupby(df["vh_scaffold"] + "|" + df["vl_scaffold"] + "|" + df["cdr3_aa"])["antigen"].apply(lambda x: ";".join(x.astype(str))).to_dict()

    print(f"Loaded {len(df)} previous antibodies for repeat checking (with BARCODEs)")
    return cdr3_set, vh_set, ab_set, cdr3_dict, vh_dict, ab_dict,ab_dict_ag

def run_combination(cfg, folder: Path):
    folder = Path(folder)
    c = cfg["combine"]
    save_greedy = bool(c.get("save_greedy_representatives", False))
    _cleanup_stale_selection_files(folder, keep_greedy=save_greedy)

    min_cdr3_len = c["min_cdr3_len"]
    max_cdr3_len = c["max_cdr3_len"]
    min_freq = c["min_freq"]
    min_count = c["min_count"]
    min_freq_sum = c["min_freq_sum"]
    min_reads = c["min_reads_per_round"]
    remove_non_func = c["remove_non_functional"]
    critical = c["critical_liabilities"]
    greedy_cutoff = c["greedy_cutoff"]  # e.g., 0.85

    # Load previous antibodies for exact repeat flags
    prev_db = cfg["general"]["previous_antibodies_db"]
    #cdr3_prev, vh_prev, ab_prev = load_previous_db(Path(prev_db))
    cdr3_prev, vh_prev, ab_prev, cdr3_barcode_dict, vh_barcode_dict, ab_barcode_dict,ab_ag_dict = load_previous_db(Path(prev_db))
    
    files = list(folder.glob("*.csv.gz"))

    if not files:
        print("No per-sample files found — Step 2 skipped.")
        return

    default_library = cfg.get("current_library", "")
    target_to_files = defaultdict(list)
    for f in files:
        target, block = _target_block_from_sample_file(f)
        sample_library = _library_from_sample_file(f, default_library)
        target_to_files[(target, block, sample_library)].append(f)

    target_block_library_counts = Counter((target, block) for target, block, _library in target_to_files)

    for (target, block, library), file_list in target_to_files.items():
        library_cfg = cfg.get("libraries", {}).get(library, {})
        library_type = library_cfg.get("library_type", "")
        pivot_cols = _combine_setting_for_library(cfg, "pivot_cols", library, c["pivot_cols"])
        cdr3_col = _combine_setting_for_library(cfg, "cdr3_col", library, c["cdr3_col"])
        pivot_cols = _normalize_grouping_columns(pivot_cols, library_type, cdr3_col)

        block_label = block or "NO_BLOCK"
        print(
            f"\n=== Combining target: {target} "
            f"(block={block_label}, {len(file_list)} samples, library={library}) ==="
        )

        target_name = _target_name_for_group(target, block)
        if block:
            print(f"   → Grouping within block: {block}")
        else:
            print(f"   → No block detected → using: {target_name}")

        if target_block_library_counts[(target, block)] > 1:
            target_name = f"{target_name}_{library}"
            print(f"   → Multiple libraries for target/block → using library-specific output: {target_name}")

        dfs = []
        for f in file_list:
            try:
                df = pd.read_csv(f)
                sample_name = f.name.replace(".csv.gz", "")
                df["sample"] = sample_name
                dfs.append(df)
            except Exception as e:
                print(f"Warning: Could not read {f}: {e}")

        if not dfs:
            print(f"No data for {target}")
            continue

        df = pd.concat(dfs, ignore_index=True)
        if "library" not in df.columns:
            df["library"] = library
        if "library_type" not in df.columns:
            df["library_type"] = library_type
        if "prevalent_targets_vh_cdr3" not in df.columns:
            df["prevalent_targets_vh_cdr3"] = 1
        if "contaminant_vh_cdr3" not in df.columns:
            df["contaminant_vh_cdr3"] = False
        df["prevalent_targets_vh_cdr3"] = pd.to_numeric(
            df["prevalent_targets_vh_cdr3"], errors="coerce"
        ).fillna(1).astype(int)
        df["contaminant_vh_cdr3"] = (
            df["contaminant_vh_cdr3"]
            .fillna(False)
            .astype(str)
            .str.lower()
            .isin({"true", "1", "yes"})
        )

        df, target_cdr3_col = _sync_cdr3_aliases(df, cdr3_col)
        if target_cdr3_col not in df.columns:
            print(f"No CDR3 column found for {target} — skipping")
            continue
        df[target_cdr3_col] = df[target_cdr3_col].astype(str).fillna("")

        available_pivot_cols = [col for col in pivot_cols if col in df.columns]
        missing_pivot_cols = [col for col in pivot_cols if col not in df.columns]
        if missing_pivot_cols:
            print(f"Skipping missing pivot columns for {target}: {missing_pivot_cols}")
        if target_cdr3_col not in available_pivot_cols:
            available_pivot_cols.insert(0, target_cdr3_col)
        pivot_cols_for_target = list(dict.fromkeys(available_pivot_cols))

        # Pandas drops groups that contain NA in any pivot index column. VHH runs
        # can have empty ANARCI CDR1/CDR2 fields, so normalize grouping fields.
        for col in pivot_cols_for_target:
            df[col] = df[col].fillna("").astype(str)
            df.loc[df[col].str.lower().isin({"nan", "none"}), col] = ""

        required_cdr_cols = [
            col for col in pivot_cols_for_target
            if col.upper().startswith("CDR")
        ]
        if required_cdr_cols:
            before_cdr_required = len(df)
            cdr_mask = pd.Series(True, index=df.index)
            for col in required_cdr_cols:
                cdr_mask &= df[col].str.len() > 0
            df = df[cdr_mask].copy()
            dropped = before_cdr_required - len(df)
            if dropped:
                print(f"Dropped {dropped} rows missing required CDR fields {required_cdr_cols} for {target}")
            if df.empty:
                print(f"No rows left after required CDR filters for {target}")
                continue

        # Filter CDR3 length
        df = df[
            (df[target_cdr3_col].str.len() >= min_cdr3_len)
            & (df[target_cdr3_col].str.len() <= max_cdr3_len)
        ]
        if df.empty:
            print(f"No rows left after CDR3 filters for {target}")
            continue

        # Remove non-functional
        if remove_non_func and "cdr3_functional" in df.columns:
            print("non-functional filtering")
            df = df[df["cdr3_functional"]]

        #OLD Pivot
        p = df.pivot_table(
            index=pivot_cols_for_target,
            columns="sample",
            values=["count", "freq"],
            aggfunc="sum",
            fill_value=0
        ).reset_index()


        # Get dominant chain per clone, now ranked by freq
        #dominant_chain = (
        #    df.sort_values("freq", ascending=False)          
        #    .groupby(pivot_cols, as_index=False)
        #    .head(1)
        #    .loc[:, pivot_cols + ["CHAIN"]]
        #    .rename(columns={"CHAIN": "dominant_CHAIN"})
        #    .drop_duplicates()
        #)
    
        # create  pivot table
        #pivot = (
        #    df.pivot_table(
        #        index=pivot_cols,
        #        columns="sample",
        #        values=["count", "freq"],
        #        aggfunc="sum",
        #        fill_value=0
        #    )
        #    .reset_index()
        #    )
        #print(pivot)
        #print(pivot.columns)

        #pivot.columns = [' '.join(col).strip() for col in pivot.columns.values]

        # Merge dominant chain in
        #p = pivot.merge(
        #    dominant_chain,
        #    on=pivot_cols,
        #    how="left"
        #)

        p.columns = [' '.join(col).strip() for col in p.columns.values]
        p = _ensure_clone_cdr3_alias(p, target_cdr3_col)

        contamination_meta = (
            df.assign(cross_target_reads=df["count"].where(df["contaminant_vh_cdr3"], 0))
            .groupby(pivot_cols_for_target, as_index=False)
            .agg(
                prevalent_targets_vh_cdr3=("prevalent_targets_vh_cdr3", "max"),
                contaminant_vh_cdr3=("contaminant_vh_cdr3", "any"),
                cross_target_reads=("cross_target_reads", "sum"),
            )
        )
        contaminant_samples = (
            df[df["contaminant_vh_cdr3"]]
            .groupby(pivot_cols_for_target)["sample"]
            .apply(lambda values: ";".join(sorted(set(values.astype(str)))))
            .reset_index(name="cross_target_samples")
        )
        p = p.merge(contamination_meta, on=pivot_cols_for_target, how="left")
        p = p.merge(contaminant_samples, on=pivot_cols_for_target, how="left")
        p["prevalent_targets_vh_cdr3"] = p["prevalent_targets_vh_cdr3"].fillna(1).astype(int)
        p["contaminant_vh_cdr3"] = p["contaminant_vh_cdr3"].fillna(False).astype(bool)
        p["cross_target_reads"] = pd.to_numeric(p["cross_target_reads"], errors="coerce").fillna(0).astype(int)
        p["cross_target_samples"] = p["cross_target_samples"].fillna("")

        metadata_cols = [
            col for col in [
                "library",
                "library_type",
                "HSEQ",
                "LSEQ",
                "CHAIN",
                "aa",
                "FR1",
                "CDR1",
                "FR2",
                "CDR2",
                "FR3",
                "FR4",
                "vl_scaffold",
                "anarci_anno",
            ]
            if col in df.columns and col not in pivot_cols_for_target
        ]
        if metadata_cols:
            dominant_meta = (
                df.sort_values("count", ascending=False)
                .groupby(pivot_cols_for_target, as_index=False)
                .first()[pivot_cols_for_target + metadata_cols]
            )
            p = p.merge(dominant_meta, on=pivot_cols_for_target, how="left")


        # Annotate liabilities
        p = annotate_liabilities(p, cdr3_col=target_cdr3_col)
        p["l_arg_value"] = p[target_cdr3_col].str.count("R")

        # Pseudo-counts
        count_cols = [c for c in p.columns if c.startswith("count ")]
        freq_cols = [c for c in p.columns if c.startswith("freq ")]

        if count_cols:
            #p[count_cols] = p[count_cols] + 1
            p[freq_cols] = p[count_cols] / p[count_cols].sum(axis=0)
        if not freq_cols:
            print(f"No frequency columns found for {target} — skipping")
            continue

        # Max freq
        p["max_freq"] = p[freq_cols].max(axis=1)
        freq_20_4 = [c for c in freq_cols if "4nM" in c or "20nM" in c]
        if freq_20_4:
            p["max_freq_20_4"] = p[freq_20_4].max(axis=1)


        clones_out = folder / f"{target_name}_clones_beforeremovelowfreq.csv"
        #p.to_csv(clones_out, index=False)

        # Filter low freq
        p = p[p["max_freq"] >= min_freq]

        ## DEBUGGING
        clones_out = folder / f"{target_name}_clones_debug_before_removelowreadcout.csv"
        #p.to_csv(clones_out, index=False)
        

        # Remove low-read rounds
        for count_col in count_cols[:]:
            if p[count_col].sum() < min_reads:
                freq_c = count_col.replace("count ", "freq ")
                p.drop(columns=[count_col, freq_c], inplace=True, errors='ignore')


        # Re-list
        freq_cols = [c for c in p.columns if c.startswith("freq ")]
        if not freq_cols or p.empty:
            print(f"No rows/frequency columns remain after low-read filtering for {target}")
            continue

        # Greedy clustering (dynamic column name)
        greedy_percent = int(round(greedy_cutoff * 100))
        cluster_col = f"greedy_cluster_{greedy_percent}"
        p[cluster_col] = greedy_clustering_by_levenshtein(p[target_cdr3_col].tolist(), greedy_cutoff)

        # Ranking
        p.sort_values(freq_cols[::-1], ascending=False, inplace=True)
        p["rank"] = range(1, len(p) + 1)

        # Critical adjustment
        p["critical"] = p[[c for c in p.columns if c in critical]].any(axis=1)
        p["rank_adjusted"] = p.apply(lambda r: 1e6 if r["critical"] else r["rank"], axis=1)
        p.sort_values(["rank_adjusted", "rank"], inplace=True)

        #=== Exact repeat flags ===
        if cdr3_prev:
            p["cdr3_repeat"] = p[target_cdr3_col].isin(cdr3_prev)
            p["vh_repeat"] = (p["vh_scaffold"] + "|" + p[target_cdr3_col]).isin(vh_prev)
            if 'vl_scaffold' in p.columns:
                p["ab_repeat"] = (p["vh_scaffold"] + "|" + p["vl_scaffold"] + "|" + p[target_cdr3_col]).isin(ab_prev)
            else:
                p["ab_repeat"] = False
        else:
            p["cdr3_repeat"] = False
            p["vh_repeat"] = False
            p["ab_repeat"] = False


        # === Exact repeat flags + BARCODEs ===
        #p["cdr3_repeat"] = p[cdr3_col].isin(cdr3_prev)
        #p["vh_repeat"] = (p["vh_scaffold"] + "|" + p[cdr3_col]).isin(vh_prev)
        #p["ab_repeat"] = (p["vh_scaffold"] + "|" + p["vl_scaffold"] + "|" + p[cdr3_col]).isin(ab_prev)
        
        # NEW: BARCODE columns (string, ";" joined if multiple, empty if no match)
        p["cdr3_repeat_barcode"] = p[target_cdr3_col].map(cdr3_barcode_dict).fillna("")
        p["vh_repeat_barcode"] = (p["vh_scaffold"] + "|" + p[target_cdr3_col]).map(vh_barcode_dict).fillna("")
        if 'vl_scaffold' in p.columns:
            p["ab_repeat_barcode"] = (p["vh_scaffold"] + "|" + p["vl_scaffold"] + "|" + p[target_cdr3_col]).map(ab_barcode_dict).fillna("")
            p["ab_repeat_antigen"] = (p["vh_scaffold"] + "|" + p["vl_scaffold"] + "|" + p[target_cdr3_col]).map(ab_ag_dict).fillna("")
        else:
            p["ab_repeat_barcode"] = ""
            p["ab_repeat_antigen"] = ""
        # Save clones
        clones_out = folder / f"{target_name}_clones.csv"
        p.to_csv(clones_out, index=False)
        print(f"Saved: {clones_out.name} ({len(p)} clones)")

        if save_greedy:
            reps = p.loc[p.groupby(cluster_col)["rank"].idxmin()]
            reps = reps.sort_values("rank")
            greedy_out = folder / f"{target_name}_greedy_{greedy_percent}.csv"
            reps.to_csv(greedy_out, index=False)
            print(f"Saved: {greedy_out.name} ({len(reps)} greedy representatives at {greedy_percent}% similarity)")

    print("\nStep 2 complete!")
