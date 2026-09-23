# pipeline/step2_combine.py
"""
Step 2: Combine per-sample tables per target
Now includes exact repeat checking from previous antibodies DB
"""

import pandas as pd
from pathlib import Path
from collections import defaultdict
from utilities.clustering import greedy_clustering_by_levenshtein
from utilities.liabilities import annotate_liabilities

def load_previous_db(db_path: Path) -> tuple[set, set, set]:
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

def run_combination(cfg, folder: Path):
    folder = Path(folder)

    c = cfg["combine"]

    pivot_cols = c["pivot_cols"]
    cdr3_col = c["cdr3_col"]

    min_cdr3_len = c["min_cdr3_len"]
    max_cdr3_len = c["max_cdr3_len"]
    min_freq = c["min_freq"]
    min_count = c["min_count"]
    min_freq_sum = c["min_freq_sum"]
    min_reads = c["min_reads_per_round"]
    remove_non_func = c["remove_non_functional"]
    critical = c["critical_liabilities"]
    greedy_cutoff = c["greedy_cutoff"]

    # Load previous antibodies for exact repeat flags
    prev_db = cfg["general"]["previous_antibodies_db"]
    cdr3_prev, vh_prev, ab_prev = load_previous_db(Path(prev_db))

    files = list(folder.glob("*.csv.gz"))

    if not files:
        print("No per-sample files found — Step 2 skipped.")
        return

    target_to_files = defaultdict(list)
    for f in files:
        target = f.stem.split("__")[0]
        target_to_files[target].append(f)

    for target, file_list in target_to_files.items():
        print(f"\n=== Combining target: {target} ({len(file_list)} samples) ===")

        dfs = []
        for f in file_list:
            try:
                df = pd.read_csv(f)
                sample_name = f.stem
                df["sample"] = sample_name
                dfs.append(df)
            except Exception as e:
                print(f"Warning: Could not read {f}: {e}")

        if not dfs:
            print(f"No data for {target}")
            continue

        df = pd.concat(dfs, ignore_index=True)

        # Filter CDR3 length
        df = df[(df[cdr3_col].str.len() >= min_cdr3_len) & (df[cdr3_col].str.len() <= max_cdr3_len)]

        # Remove non-functional
        if remove_non_func:
            df = df[df["cdr3_functional"]]


        # Pivot
        p = df.pivot_table(
            index=pivot_cols,
            columns="sample",
            values=["count", "freq"],
            aggfunc="sum",
            fill_value=0
        ).reset_index()

        p.columns = [' '.join(col).strip() for col in p.columns.values]

        # Annotate liabilities
        p = annotate_liabilities(p, cdr3_col=cdr3_col)
        p["l_arg_value"] = p[cdr3_col].str.count("R")

        # Pseudo-counts
        count_cols = [c for c in p.columns if c.startswith("count ")]
        freq_cols = [c for c in p.columns if c.startswith("freq ")]

        if count_cols:
            p[count_cols] = p[count_cols] + 1
            p[freq_cols] = p[count_cols] / p[count_cols].sum(axis=0)

        # Max freq
        p["max_freq"] = p[freq_cols].max(axis=1)
        freq_20_4 = [c for c in freq_cols if "4nM" in c or "20nM" in c]
        if freq_20_4:
            p["max_freq_20_4"] = p[freq_20_4].max(axis=1)

        # Filter low freq
        p = p[p["max_freq"] >= min_freq]

        # Remove low-read rounds
        for c in count_cols[:]:
            if p[c].sum() < min_reads:
                freq_c = c.replace("count ", "freq ")
                p.drop(columns=[c, freq_c], inplace=True, errors='ignore')

        # Re-list
        freq_cols = [c for c in p.columns if c.startswith("freq ")]

        # Greedy clustering
        p["greedy_cluster_0.85"] = greedy_clustering_by_levenshtein(p[cdr3_col].tolist(), greedy_cutoff)

        # Ranking
        p.sort_values(freq_cols[::-1], ascending=False, inplace=True)
        p["rank"] = range(1, len(p) + 1)

        # Critical adjustment
        p["critical"] = p[[c for c in p.columns if c in critical]].any(axis=1)
        p["rank_adjusted"] = p.apply(lambda r: 1e6 if r["critical"] else r["rank"], axis=1)
        p.sort_values(["rank_adjusted", "rank"], inplace=True)

        # === Exact repeat flags — AFTER pivot (fast & efficient) ===
        if cdr3_prev:
            p["cdr3_repeat"] = p[cdr3_col].isin(cdr3_prev)
            p["vh_repeat"] = (p["vh_scaffold"] + "|" + p[cdr3_col]).isin(vh_prev)
            p["ab_repeat"] = (p["vh_scaffold"] + "|" + p["vl_scaffold"] + "|" + p[cdr3_col]).isin(ab_prev)
        else:
            p["cdr3_repeat"] = False
            p["vh_repeat"] = False
            p["ab_repeat"] = False

        # Save clones
        clones_out = folder / f"{target}_clones.csv.gz"
        p.to_csv(clones_out, index=False, compression="gzip")
        print(f"Saved: {clones_out.name} ({len(p)} clones)")
        
        p["max_freq_sum"] = p[freq_cols].sum(axis=1)
        leads = p[p["max_freq_sum"] >= min_freq_sum]
        leads = leads.sort_values(["rank_adjusted", "rank"])
        leads_out = folder / f"{target}_leads.csv.gz"
        leads.to_csv(leads_out, index=False, compression="gzip")
        print(f"Saved: {leads_out.name} ({len(leads)} leads)")

    print("\nStep 2 complete! Exact repeat flags (cdr3_repeat, vh_repeat, ab_repeat) added.")