# pipeline/step5_plot_report.py
"""
Step 5: Generate QC and analysis plots + single HTML & PDF report
Saves individual PNGs in results/plots/ and creates results/report.html + report.pdf
Now with rarefaction curves from individual samples (*.csv.gz)
- All samples
- Only 100nM samples
"""

import pandas as pd
from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", str(Path(os.environ.get("TMPDIR", "/tmp")) / "matplotlib"))

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import base64
import json
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
from io import BytesIO
from datetime import datetime
import re

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_COLORS = {
    "A": "#4daf4a", "C": "#ffd92f", "D": "#e41a1c", "E": "#e41a1c",
    "F": "#377eb8", "G": "#999999", "H": "#984ea3", "I": "#4daf4a",
    "K": "#984ea3", "L": "#4daf4a", "M": "#4daf4a", "N": "#ff7f00",
    "P": "#a65628", "Q": "#ff7f00", "R": "#984ea3", "S": "#ff7f00",
    "T": "#ff7f00", "V": "#4daf4a", "W": "#377eb8", "Y": "#377eb8",
}
LOGO_FONT = FontProperties(family="DejaVu Sans", weight="bold")
DELPHI_PROFILE_SCORE_THRESHOLD = 0.5
DELPHI_FOLD_CHANGE_SCORE_THRESHOLD = 0.5
DELPHI_SCORE_LABELS = {"psr": "PSR", "sec": "SEC"}
SEQUENCE_LOGO_MIN_MAX_FREQ = 0.00005
DIVERSITY_METRICS = [
    ("shannon", "Shannon"),
    ("shannon_min2", "Shannon (>1 read clones)"),
    ("inv_simpsons", "Inverse Simpson"),
    ("inv_simpsons_min2", "Inverse Simpson (>1 read clones)"),
    ("inv_simpsons_2", "Inverse Simpson (>1 read clones)"),
]


def _format_threshold(value: float) -> str:
    return f"{value:.8f}".rstrip("0").rstrip(".")


def fig_to_base64(fig):
    """Convert matplotlib figure to base64 encoded PNG"""
    buf = BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')


def image_to_base64(path: Path) -> str:
    with open(path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode("utf-8")


def _safe_filename(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return value[:160] or "plot"


def save_fig_plot(fig, plots_dir: Path, filename: str, *, dpi: int = 200) -> Path:
    plots_dir.mkdir(exist_ok=True)
    path = plots_dir / _safe_filename(filename)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path


def plot_html_from_file(path: Path, title: str, caption: str = "") -> str:
    caption_html = f'<p class="caption">{caption}</p>' if caption else ""
    return (
        f'<div class="plot"><h2>{title}</h2>'
        f'<img src="data:image/png;base64,{image_to_base64(path)}">'
        f'{caption_html}</div>'
    )


def _num(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(0)


def _lead_label_table(df: pd.DataFrame) -> pd.DataFrame:
    table_cols = [col for col in ["vh_scaffold", "vl_scaffold", "cdr3_aa", "max_freq"] if col in df.columns]
    table_data = df[table_cols].copy()
    if "max_freq" in table_data.columns:
        table_data["max_freq"] = pd.to_numeric(table_data["max_freq"], errors="coerce").fillna(0).round(6)

    def label(row):
        parts = [str(row[col]) for col in ["vh_scaffold", "vl_scaffold", "cdr3_aa"] if col in row and str(row[col])]
        freq = f" ({row['max_freq']})" if "max_freq" in row else ""
        return ":".join(parts) + freq

    table_data["lead"] = table_data.apply(label, axis=1)
    return table_data[["lead"]]


def _value_counts_dict(df: pd.DataFrame, col: str) -> dict[str, int]:
    if col not in df.columns:
        return {}
    values = df[col].fillna("").astype(str).str.strip()
    values = values[values != ""]
    return {str(k): int(v) for k, v in values.value_counts().to_dict().items()}


def _count_final_lead_rows(final_lead_files: list[Path]) -> int:
    total = 0
    for path in final_lead_files:
        try:
            total += len(pd.read_excel(path))
        except Exception:
            continue
    return total


def load_ml_summary(folder: Path) -> pd.DataFrame:
    summary_path = folder / "ml_summary.csv"
    if summary_path.exists():
        try:
            return pd.read_csv(summary_path)
        except Exception:
            return pd.DataFrame()

    status_path = folder / "ml_prediction_status.json"
    if not status_path.exists():
        return pd.DataFrame()
    try:
        status = json.loads(status_path.read_text())
    except Exception:
        return pd.DataFrame()

    rows = []
    for item in status:
        input_path = item.get("source_input") or item.get("input") or item.get("output")
        target = Path(str(input_path)).stem if input_path else ""
        target = target.replace("_final_leads_delphi_input", "").replace("_final_leads", "")
        total_rows = int(item.get("rows") or 0)
        pass_rows = int(item.get("ml_rows") or 0)
        fail_rows = int(item.get("dropped_rows") if item.get("dropped_rows") is not None else max(total_rows - pass_rows, 0))
        rows.append(
            {
                "target": target,
                "status": item.get("status", ""),
                "backend": item.get("backend", ""),
                "total_rows": total_rows,
                "ml_pass": pass_rows,
                "ml_fail": fail_rows,
                "ml_pass_rate_pct": round(100 * pass_rows / total_rows, 2) if total_rows else 0.0,
                "prediction_columns": ";".join(map(str, item.get("merged_prediction_columns") or [])),
                "sequence_mode": item.get("sequence_mode", ""),
                "error": item.get("error", ""),
            }
        )
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary.to_csv(summary_path, index=False)
    return summary


def write_data_summary(folder: Path, qc: pd.DataFrame, prevalent_df: pd.DataFrame, ml_summary: pd.DataFrame | None = None) -> dict:
    clone_files = sorted(folder.glob("*_clones.csv"))
    final_lead_files = sorted((folder / "by_protein").glob("*_final_leads.xlsx"))
    leads_path = folder / "leads.xlsx"
    lead_count = 0
    if leads_path.exists():
        try:
            lead_count = len(pd.read_excel(leads_path))
        except Exception:
            lead_count = 0
    if lead_count == 0:
        lead_count = _count_final_lead_rows(final_lead_files)

    ml_summary = ml_summary if ml_summary is not None else pd.DataFrame()
    ml_total = int(_num(ml_summary, "total_rows").sum()) if not ml_summary.empty else 0
    ml_pass = int(_num(ml_summary, "ml_pass").sum()) if not ml_summary.empty else 0

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "results_folder": str(folder),
        "sample_count": int(len(qc)),
        "target_count": len(clone_files),
        "global_leads": int(lead_count),
        "per_target_final_lead_files": len(final_lead_files),
        "total_reads": int(_num(qc, "total").sum()),
        "merged_reads": int(_num(qc, "merged").sum()),
        "median_merged_pct": float(_num(qc, "merged_pct").median()) if "merged_pct" in qc.columns else 0.0,
        "total_unique_cdr3": int(_num(qc, "unique_cdr3").sum()),
        "libraries": _value_counts_dict(qc, "library"),
        "samples_with_cross_target_reads": int((_num(qc, "n_cross_target_reads") > 0).sum())
        if "n_cross_target_reads" in qc.columns
        else 0,
        "total_cross_target_clones": int(_num(qc, "n_cross_target_clones").sum()),
        "total_cross_target_reads": int(_num(qc, "n_cross_target_reads").sum()),
        "max_cross_target_reads_pct": float(_num(qc, "pct_cross_target_reads").max())
        if "pct_cross_target_reads" in qc.columns
        else 0.0,
        "prevalent_vh_cdr3_pairs_top50": int(len(prevalent_df)),
        "ml_targets": int(len(ml_summary)) if not ml_summary.empty else 0,
        "ml_total_rows": ml_total,
        "ml_pass": ml_pass,
        "ml_fail": max(ml_total - ml_pass, 0),
        "ml_pass_rate_pct": round(100 * ml_pass / ml_total, 2) if ml_total else 0.0,
    }

    pd.DataFrame([summary]).to_csv(folder / "data_summary.csv", index=False)
    with open(folder / "data_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    return summary

def rarefaction_curve(counts, max_reads=150000, steps=100, iterations=10):
    """Compute rarefaction curve for a single sample"""
    total_reads = counts.sum()
    if total_reads == 0:
        return [0], [0]

    read_points = np.linspace(1, min(total_reads, max_reads), steps, dtype=int)
    diversities = []

    for reads in read_points:
        iter_div = []
        for _ in range(iterations):
            sampled = np.random.choice(len(counts), size=reads, p=counts/total_reads, replace=True)
            iter_div.append(len(np.unique(sampled)))
        diversities.append(np.mean(iter_div))

    return read_points, diversities


def generate_rarefaction_plots(folder: Path):
    """Generate rarefaction curves from individual sample files (*.csv.gz)"""
    folder = Path(folder)

    # Load all individual sample files
    sample_files = list(folder.glob("*.csv.gz"))  # per-sample files

    if not sample_files:
        print("No individual sample files (*.csv.gz) found for rarefaction")
        return

    plots_dir = folder / "plots"
    plots_dir.mkdir(exist_ok=True)

    sample_data = {}
    for f in sample_files:
        # Clean sample name: remove .csv.gz and any leftover .csv
        sample_name = f.name.replace('.csv.gz', '').replace('.csv', '')
        try:
            df = pd.read_csv(f)
            if "count" not in df.columns:
                continue
            counts = df["count"].values
            sample_data[sample_name] = counts
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}")

    if not sample_data:
        print("No valid sample data for rarefaction")
        return

    # Sort sample names alphabetically for consistent numbering
    sorted_samples = sorted(sample_data.items(), key=lambda x: x[0])

    # === All samples ===
    plt.figure(figsize=(16, 10))
    for i, (sample_name, counts) in enumerate(sorted_samples, start=1):
        x, y = rarefaction_curve(counts)
        line, = plt.plot(x, y, alpha=0.8)
        color = line.get_color()

        # Legend label with number + clean name
        plt.plot([], [], label=f"{i} - {sample_name}", color=color)

        # Annotate number on curve
        if len(x) > 10:
            annot_x = x[int(0.9 * len(x))]
            annot_y = y[int(0.9 * len(y))]
            plt.text(annot_x, annot_y, str(i), fontsize=12, fontweight='bold',
                     color=color, ha='right', va='bottom',
                     bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    plt.xlabel("Subsampled Reads")
    plt.ylabel("Unique CDR3 Sequences")
    plt.title("Rarefaction Curves — All Samples\n(Numbers on curves and in legend correspond)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9, ncol=1)
    plt.tight_layout()
    plt.savefig(plots_dir / "rarefaction_all_samples.png", dpi=200)
    plt.close()

    # === Only 100nM samples ===
    nM100_samples = {k: v for k, v in sample_data.items() if "100nM" in k or "100 nM" in k}
    if nM100_samples:
        sorted_100nM = sorted(nM100_samples.items(), key=lambda x: x[0])
        plt.figure(figsize=(16, 10))
        for i, (sample_name, counts) in enumerate(sorted_100nM, start=1):
            x, y = rarefaction_curve(counts)
            line, = plt.plot(x, y, alpha=0.8)
            color = line.get_color()

            plt.plot([], [], label=f"{i} - {sample_name}", color=color)

            if len(x) > 10:
                annot_x = x[int(0.9 * len(x))]
                annot_y = y[int(0.9 * len(y))]
                plt.text(annot_x, annot_y, str(i), fontsize=12, fontweight='bold',
                         color=color, ha='right', va='bottom',
                         bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

        plt.xlabel("Subsampled Reads")
        plt.ylabel("Unique CDR3 Sequences")
        plt.title("Rarefaction Curves — 100nM Samples Only\n(Numbers on curves and in legend correspond)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9, ncol=1)
        plt.tight_layout()
        plt.savefig(plots_dir / "rarefaction_100nM_samples.png", dpi=200)
        plt.close()

    print("Rarefaction curves generated with numbered labels, curve annotations, and clean sample names (no .csv)")


def _prediction_score_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in df.columns if str(col).endswith("_score")]


def _delphi_label_score_columns(df: pd.DataFrame, label: str) -> list[str]:
    label = label.lower()
    prefix = f"{label}_"
    score_cols: list[str] = []
    for col in df.columns:
        name = str(col)
        lowered = name.lower()
        if lowered == f"mean_{label}_score":
            continue
        if not lowered.startswith(prefix) or not lowered.endswith("_score"):
            continue
        if "xgboost" in lowered or "transformer" in lowered:
            score_cols.append(name)
    return sorted(score_cols)


def _delphi_label_mean_score(df: pd.DataFrame, label: str) -> pd.Series:
    label = label.lower()
    mean_col = f"mean_{label}_score"
    if mean_col in df.columns:
        return pd.to_numeric(df[mean_col], errors="coerce")
    score_cols = _delphi_label_score_columns(df, label)
    if score_cols:
        return df[score_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    return pd.Series(np.nan, index=df.index)


def _short_delphi_model_label(col: str, label: str) -> str:
    label = label.lower()
    lowered = str(col).lower()
    if lowered == f"mean_{label}_score":
        return f"mean {DELPHI_SCORE_LABELS.get(label, label.upper())}"
    text = str(col)
    text = re.sub(rf"^{label}_", "", text, flags=re.IGNORECASE)
    text = re.sub(r"_score$", "", text, flags=re.IGNORECASE)
    text = re.sub(r"_ipi_(psr_trainset|sec_5000)$", "", text, flags=re.IGNORECASE)
    text = text.replace("transformer_lm_igbert", "transformer_lm igbert")
    text = text.replace("xgboost_biophysical", "xgboost biophysical")
    text = text.replace("xgboost_ablang", "xgboost ablang")
    return text.replace("_", " ")


def _delphi_consensus_score_columns(df: pd.DataFrame) -> list[str]:
    score_cols: list[str] = []
    for label in DELPHI_SCORE_LABELS:
        score_cols.extend(_delphi_label_score_columns(df, label))
    return sorted(score_cols, key=lambda col: (not col.startswith("psr_"), col))


def _delphi_consensus_score(df: pd.DataFrame) -> pd.Series:
    score_cols = _delphi_consensus_score_columns(df)
    if score_cols:
        return df[score_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    if "mean_delphi_score" in df.columns:
        return pd.to_numeric(df["mean_delphi_score"], errors="coerce")
    return pd.Series(np.nan, index=df.index)


def _delphi_psr_sec_global_mean_score(df: pd.DataFrame) -> pd.Series:
    """Fold-change highlight score: average of the PSR mean and SEC mean."""
    psr = _delphi_label_mean_score(df, "psr")
    sec = _delphi_label_mean_score(df, "sec")
    score = (psr + sec) / 2
    score[psr.isna() | sec.isna()] = np.nan
    return score


def _clean_profile_sequence(value) -> str:
    if pd.isna(value):
        return ""
    seq = re.sub(r"[^A-Za-z]", "", str(value).upper())
    return seq[1:] if seq.startswith("CAR") else seq


def _charge_at_ph(seq: str, ph: float) -> float:
    counts = {aa: seq.count(aa) for aa in "CDEHKRY"}
    positive = (
        1 / (1 + 10 ** (ph - 7.5))
        + counts["K"] / (1 + 10 ** (ph - 10.5))
        + counts["R"] / (1 + 10 ** (ph - 12.5))
        + counts["H"] / (1 + 10 ** (ph - 6.0))
    )
    negative = (
        1 / (1 + 10 ** (3.55 - ph))
        + counts["D"] / (1 + 10 ** (3.9 - ph))
        + counts["E"] / (1 + 10 ** (4.1 - ph))
        + counts["C"] / (1 + 10 ** (8.3 - ph))
        + counts["Y"] / (1 + 10 ** (10.1 - ph))
    )
    return positive - negative


def _sequence_isoelectric_point(seq: str) -> float:
    seq = _clean_profile_sequence(seq)
    if not seq:
        return np.nan
    low, high = 0.0, 14.0
    for _ in range(40):
        mid = (low + high) / 2
        if _charge_at_ph(seq, mid) > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2, 3)


def _sequence_net_charge(seq: str) -> float:
    seq = _clean_profile_sequence(seq)
    if not seq:
        return np.nan
    return float(seq.count("K") + seq.count("R") + seq.count("H") - seq.count("D") - seq.count("E"))


def _round_sort_value(value: str) -> tuple[int, str]:
    text = str(value)
    for pattern in [r"(?:round|rnd|r)[_\s-]*(\d+)", r"__(\d+)(?:__|$)", r"\b(\d+)\b"]:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return (int(match.group(1)), text)
    lowered = text.lower()
    if "input" in lowered:
        return (-1, text)
    if "negative" in lowered:
        return (-2, text)
    return (0, text)


def _freq_round_label(freq_col: str) -> str:
    text = str(freq_col).replace("freq ", "")
    parts = text.split("__")
    if len(parts) >= 3:
        return parts[2]
    match = re.search(r"(round[_\s-]*\d+|r[_\s-]*\d+)", text, flags=re.IGNORECASE)
    return match.group(1) if match else text


def make_delphi_score_plots(folder: Path, plots_dir: Path) -> list[dict[str, Path | str]]:
    entries: list[dict[str, Path | str]] = []
    for label, display in DELPHI_SCORE_LABELS.items():
        rows = []
        for final_file in sorted((folder / "by_protein").glob("*_final_leads.xlsx")):
            try:
                df = pd.read_excel(final_file)
            except Exception:
                continue
            score = _delphi_label_mean_score(df, label)
            if score.isna().all():
                continue
            target = final_file.stem.replace("_final_leads", "")
            target_scores = pd.DataFrame({"target": target, "score": score}).dropna(subset=["score"])
            if len(target_scores) > 5000:
                target_scores = target_scores.sample(5000, random_state=1)
            if target_scores.empty:
                continue
            rows.append(target_scores)

        if not rows:
            continue
        plot_df = pd.concat(rows, ignore_index=True)
        if plot_df.empty:
            continue

        target_order = (
            plot_df.groupby("target")["score"]
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )
        target_count = max(len(target_order), 1)
        fig_height = min(max(7, target_count * 0.38 + 2.5), 22)
        fig, ax = plt.subplots(figsize=(12, fig_height))
        sns.boxplot(data=plot_df, y="target", x="score", order=target_order, ax=ax, color="#6baed6", fliersize=1.8)
        ax.axvline(DELPHI_PROFILE_SCORE_THRESHOLD, color="#d95f02", linestyle="--", linewidth=1.2)
        ax.set_title(f"Delphi {display} Mean Score Distribution by Target")
        ax.set_xlabel(f"Mean {display} Delphi score")
        ax.set_ylabel("Target")
        ax.tick_params(axis="y", labelsize=9)
        ax.set_xlim(-0.02, 1.02)
        ax.grid(axis="x", color="#d0d0d0", linewidth=0.8, alpha=0.8)
        ax.text(
            DELPHI_PROFILE_SCORE_THRESHOLD + 0.01,
            0.99,
            f">={DELPHI_PROFILE_SCORE_THRESHOLD:g}",
            transform=ax.get_xaxis_transform(),
            color="#d95f02",
            fontsize=9,
            va="top",
        )
        path = save_fig_plot(fig, plots_dir, f"delphi_{label}_score_distribution.png")
        plt.close(fig)
        entries.append(
            {
                "path": path,
                "title": f"Delphi {display} Score Distribution",
                "caption": (
                    f"Distribution of mean {display} Delphi score per target. "
                    "Individual model score columns are not plotted here, keeping large multi-target runs readable."
                ),
            }
        )
    return entries


def make_biophysical_profile_plots(folder: Path, plots_dir: Path) -> list[dict[str, Path | str]]:
    entries: list[dict[str, Path | str]] = []
    for label, display in DELPHI_SCORE_LABELS.items():
        entry = _make_biophysical_profile_plot(folder, plots_dir, label, display)
        if entry:
            entries.append(entry)
    return entries


def _make_biophysical_profile_plot(folder: Path, plots_dir: Path, score_label: str, display: str) -> dict[str, Path | str] | None:
    rows = []
    for final_file in sorted((folder / "by_protein").glob("*_final_leads.xlsx")):
        try:
            df = pd.read_excel(final_file)
        except Exception:
            continue
        if df.empty:
            continue
        cdr3_col = "CDR3" if "CDR3" in df.columns else "cdr3_aa" if "cdr3_aa" in df.columns else ""
        if not cdr3_col:
            continue
        score = _delphi_label_mean_score(df, score_label)
        if score.isna().all():
            continue

        cdr3_seq = df[cdr3_col].apply(_clean_profile_sequence)
        profile = pd.DataFrame(
            {
                "target": final_file.stem.replace("_final_leads", ""),
                "delphi_score": score,
                "cdr3_arg_count": cdr3_seq.str.count("R").astype(float),
                "cdr3_asp_count": cdr3_seq.str.count("D").astype(float),
                "cdr3_trp_count": cdr3_seq.str.count("W").astype(float),
                "cdr3_glu_count": cdr3_seq.str.count("E").astype(float),
                "cdr3_length": cdr3_seq.str.len().astype(float),
                "cdr3_net_charge": cdr3_seq.apply(_sequence_net_charge),
                "cdr3_isoelectric_point": cdr3_seq.apply(_sequence_isoelectric_point),
            }
        )
        if "hseq_isoelectric_point" in df.columns:
            profile["hseq_isoelectric_point"] = pd.to_numeric(df["hseq_isoelectric_point"], errors="coerce")
        elif "HSEQ" in df.columns:
            profile["hseq_isoelectric_point"] = df["HSEQ"].apply(_sequence_isoelectric_point)
        profile = profile.dropna(subset=["delphi_score"])
        if len(profile) > 8000:
            profile = profile.sample(8000, random_state=7)
        rows.append(profile)

    if not rows:
        return None

    plot_df = pd.concat(rows, ignore_index=True)
    plot_df["prediction"] = np.where(
        plot_df["delphi_score"] >= DELPHI_PROFILE_SCORE_THRESHOLD,
        "High confidence",
        "Other",
    )
    pass_n = int((plot_df["prediction"] == "High confidence").sum())
    fail_n = int((plot_df["prediction"] == "Other").sum())

    features = [
        ("cdr3_arg_count", "Arginine count (HCDR3)", None),
        ("cdr3_asp_count", "Aspartic acid count (HCDR3)", None),
        ("cdr3_trp_count", "Tryptophan count (HCDR3)\n(Arg count=1)", lambda data: data[data["cdr3_arg_count"] == 1]),
        ("cdr3_glu_count", "Glutamic acid count (HCDR3)", None),
        ("cdr3_length", "HCDR3 loop length", None),
        ("cdr3_net_charge", "Net charge (HCDR3)", None),
        ("cdr3_isoelectric_point", "Isoelectric point (HCDR3)", None),
        ("hseq_isoelectric_point", "Isoelectric point\n(heavy chain)", None),
    ]
    features = [feature for feature in features if feature[0] in plot_df.columns]
    if not features:
        return None

    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes_flat = axes.ravel()
    colors = {"High confidence": "#5da5e8", "Other": "#f4a259"}
    for idx, (col, feature_label, subset_fn) in enumerate(features[:8]):
        ax = axes_flat[idx]
        feature_df = subset_fn(plot_df) if subset_fn else plot_df
        values = pd.to_numeric(feature_df[col], errors="coerce").replace([np.inf, -np.inf], np.nan)
        feature_df = feature_df.assign(_value=values).dropna(subset=["_value"])
        if feature_df.empty:
            ax.set_visible(False)
            continue

        clean_values = feature_df["_value"]
        if clean_values.nunique() <= 1:
            bins = 3
        elif np.allclose(clean_values.dropna(), clean_values.dropna().round()):
            bins = np.arange(np.floor(clean_values.min()) - 0.5, np.ceil(clean_values.max()) + 1.5, 1)
        else:
            bins = min(24, max(8, int(np.sqrt(len(feature_df)))))

        for pred_label in ["High confidence", "Other"]:
            group = feature_df.loc[feature_df["prediction"] == pred_label, "_value"]
            if group.empty:
                continue
            ax.hist(group, bins=bins, density=True, alpha=0.68, color=colors[pred_label], label=pred_label)
        ax.set_xlabel(feature_label, fontsize=11)
        ax.set_ylabel("Density", fontsize=11)
        ax.text(0.02, 0.96, chr(ord("a") + idx), transform=ax.transAxes, fontsize=15, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)

    for idx in range(len(features), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    handles = [
        plt.Line2D(
            [0],
            [0],
            color=colors["High confidence"],
            lw=8,
            alpha=0.68,
            label=f"High confidence (>={DELPHI_PROFILE_SCORE_THRESHOLD:g}), n={pass_n:,}",
        ),
        plt.Line2D(
            [0],
            [0],
            color=colors["Other"],
            lw=8,
            alpha=0.68,
            label=f"Other (<{DELPHI_PROFILE_SCORE_THRESHOLD:g}), n={fail_n:,}",
        ),
    ]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.01, 0.52), frameon=False, fontsize=11)
    fig.suptitle(f"Delphi {display} biophysical profile", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=[0.08, 0, 1, 0.94])
    path = save_fig_plot(fig, plots_dir, f"delphi_{score_label}_biophysical_profile.png")
    plt.close(fig)
    return {
        "path": path,
        "title": f"Delphi {display} Biophysical Profile",
        "caption": (
            f"Density profiles split by mean {display} Delphi score: "
            f"high confidence >= {DELPHI_PROFILE_SCORE_THRESHOLD:g}, "
            f"other < {DELPHI_PROFILE_SCORE_THRESHOLD:g}."
        ),
    }


def make_sequencing_depth_section(qc: pd.DataFrame, plots_dir: Path, plot_count: int) -> tuple[str, int]:
    plot_count += 1
    fig, ax = plt.subplots(figsize=(max(20, len(qc) * 0.5), 10))
    qc_melt = qc.melt(
        id_vars="name",
        value_vars=["total", "merged"],
        var_name="type",
        value_name="reads",
    )
    qc_melt = qc_melt.sort_values(["reads", "name"], ascending=[False, True])

    sns.barplot(data=qc_melt, x="name", y="reads", hue="type", ax=ax, palette="Set2")
    ax.set_title("Sequencing Depth and Merging Efficiency per Sample")
    ax.tick_params(axis="x", rotation=90, labelsize=18)
    ax.set_xlabel("Sample Name")
    ax.set_ylabel("Reads")
    ax.legend(title="Type")
    ax.grid(axis="y", linestyle="--", alpha=0.7)
    plot_path = save_fig_plot(fig, plots_dir, "sequencing_depth.png")
    plt.close(fig)
    return (
        plot_html_from_file(
            plot_path,
            f"{plot_count}. Sequencing Depth",
            "Total vs merged reads per sample (sorted by depth, wide view for many samples)",
        ),
        plot_count,
    )


def cleanup_stale_clone_frequency_plots(plots_dir: Path) -> None:
    for path in plots_dir.glob("clone_frequency_*.png"):
        try:
            path.unlink()
        except Exception as exc:
            print(f"Warning: could not remove stale clone-frequency plot {path.name}: {exc}")


def cleanup_stale_delphi_report_plots(plots_dir: Path) -> None:
    stale_paths = [plots_dir / "ml_score_distribution.png", plots_dir / "delphi_biophysical_profile.png"]
    stale_paths.extend(plots_dir.glob("delphi_*_score_distribution.png"))
    stale_paths.extend(plots_dir.glob("delphi_*_biophysical_profile.png"))
    for path in stale_paths:
        if not path.exists():
            continue
        try:
            path.unlink()
        except Exception as exc:
            print(f"Warning: could not remove stale Delphi report plot {path.name}: {exc}")


def _available_diversity_metrics(df: pd.DataFrame) -> list[tuple[str, str]]:
    metrics = []
    seen = set()
    for metric, label in DIVERSITY_METRICS:
        if metric in df.columns and metric not in seen:
            metrics.append((metric, label))
            seen.add(metric)
    return metrics


def _clean_target_part(value) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "na"} else text


def _target_labels_from_qc(qc: pd.DataFrame) -> pd.Series:
    if "antigen" in qc.columns:
        labels = qc["antigen"].apply(_clean_target_part)
    elif "name" in qc.columns:
        labels = qc["name"].apply(_clean_target_part)
    else:
        labels = pd.Series(["Target"] * len(qc), index=qc.index)

    if "block" in qc.columns:
        blocks = qc["block"].apply(_clean_target_part)
        labels = [
            f"{label}_{block}" if block and block not in label else label
            for label, block in zip(labels, blocks)
        ]
        labels = pd.Series(labels, index=qc.index)

    if "library" in qc.columns:
        libraries = qc["library"].apply(_clean_target_part)
        group_key = pd.DataFrame({"target": labels, "library": libraries})
        multi_library_targets = (
            group_key[group_key["library"] != ""]
            .groupby("target")["library"]
            .nunique()
        )
        multi_library_targets = set(multi_library_targets[multi_library_targets > 1].index)
        labels = [
            f"{label}_{library}" if label in multi_library_targets and library else label
            for label, library in zip(labels, libraries)
        ]
        labels = pd.Series(labels, index=qc.index)

    return labels.replace("", "Target")


def _round_order(values: pd.Series) -> list[str]:
    non_numeric_order = ["Input", "Negative", "NA"]

    def sort_key(value):
        text = str(value)
        if text in non_numeric_order:
            return (0, non_numeric_order.index(text), text)
        match = re.search(r"\d+", text)
        return (1, int(match.group()), text) if match else (2, 0, text)

    return sorted(values.dropna().astype(str).unique(), key=sort_key)


def _is_inverse_simpson_metric(metric: str) -> bool:
    return str(metric).startswith("inv_simpson")


def _diversity_plot_values(values: pd.Series, metric: str, display: str) -> tuple[pd.Series, str]:
    numeric = pd.to_numeric(values, errors="coerce")
    if _is_inverse_simpson_metric(metric):
        transformed = np.log10(numeric.clip(lower=0).fillna(0) + 1)
        transformed[numeric.isna()] = np.nan
        return transformed, f"log10(1 + {display})"
    return numeric, f"{display} diversity"


def make_target_diversity_plots(qc: pd.DataFrame, plots_dir: Path) -> list[dict[str, Path | str]]:
    for path in plots_dir.glob("diversity_by_target_*.png"):
        try:
            path.unlink()
        except Exception as exc:
            print(f"Warning: could not remove stale target-diversity plot {path.name}: {exc}")

    summary_path = plots_dir.parent / "diversity_by_target_summary.csv"
    if summary_path.exists():
        try:
            summary_path.unlink()
        except Exception as exc:
            print(f"Warning: could not remove stale target-diversity summary: {exc}")

    entries: list[dict[str, Path | str]] = []
    summary_rows: list[dict] = []
    metrics = _available_diversity_metrics(qc)
    if not metrics:
        return entries

    qc_base = qc.copy()
    qc_base["_target_label"] = _target_labels_from_qc(qc_base)
    hue_col = "round" if "round" in qc_base.columns else "condition" if "condition" in qc_base.columns else None
    hue_order = _round_order(qc_base[hue_col]) if hue_col == "round" else None

    for metric, display in metrics:
        plot_df = qc_base.copy()
        plot_df[metric] = pd.to_numeric(plot_df[metric], errors="coerce")
        plot_df["_diversity_value"], axis_label = _diversity_plot_values(plot_df[metric], metric, display)
        plot_df = plot_df.dropna(subset=[metric, "_diversity_value", "_target_label"])
        if plot_df.empty:
            continue

        target_summary = (
            plot_df.groupby("_target_label", as_index=False)
            .agg(
                mean_value=(metric, "mean"),
                median_value=(metric, "median"),
                min_value=(metric, "min"),
                max_value=(metric, "max"),
                samples=(metric, "size"),
            )
            .sort_values("median_value", ascending=False)
        )
        for _, row in target_summary.iterrows():
            summary_rows.append(
                {
                    "metric": metric,
                    "metric_label": display,
                    "target": row["_target_label"],
                    "mean": round(float(row["mean_value"]), 6),
                    "median": round(float(row["median_value"]), 6),
                    "min": round(float(row["min_value"]), 6),
                    "max": round(float(row["max_value"]), 6),
                    "samples": int(row["samples"]),
                }
            )

        target_order = target_summary["_target_label"].tolist()
        fig_height = min(max(6.5, len(target_order) * 0.42 + 2.0), 24)
        fig, ax = plt.subplots(figsize=(13.5, fig_height))
        sns.boxplot(
            data=plot_df,
            y="_target_label",
            x="_diversity_value",
            order=target_order,
            ax=ax,
            color="lightgray",
            fliersize=0,
        )
        if hue_col:
            hue_values = plot_df[hue_col].dropna().astype(str).nunique()
            palette = sns.color_palette("tab20" if hue_values > 10 else "tab10", hue_values)
            sns.stripplot(
                data=plot_df,
                y="_target_label",
                x="_diversity_value",
                hue=hue_col,
                hue_order=hue_order,
                order=target_order,
                palette=palette,
                jitter=0.22,
                size=4.5,
                alpha=0.85,
                ax=ax,
            )
            ax.legend(title=hue_col.title(), bbox_to_anchor=(1.02, 1), loc="upper left")
        else:
            sns.stripplot(
                data=plot_df,
                y="_target_label",
                x="_diversity_value",
                order=target_order,
                color="#1f77b4",
                jitter=0.22,
                size=4.5,
                alpha=0.85,
                ax=ax,
            )

        ax.set_title(f"{display} Diversity by Target")
        ax.set_xlabel(axis_label)
        ax.set_ylabel("Target")
        ax.tick_params(axis="y", labelsize=9)
        ax.grid(axis="x", linestyle="--", alpha=0.45)
        sns.despine(ax=ax, left=True, bottom=False)
        path = save_fig_plot(fig, plots_dir, f"diversity_by_target_{metric}.png")
        plt.close(fig)
        entries.append(
            {
                "path": path,
                "title": f"{display} Diversity by Target",
                "caption": (
                    "Box plots summarize sample diversity per target; points show individual samples "
                    + (f"colored by {hue_col}." if hue_col else "within each target.")
                    + (" Inverse Simpson values are shown as log10(1 + value)." if _is_inverse_simpson_metric(metric) else "")
                ),
            }
        )

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    return entries


def _clean_logo_sequence(value) -> str:
    if pd.isna(value):
        return ""
    seq = re.sub(r"[^A-Za-z]", "", str(value).upper())
    return "".join(aa for aa in seq if aa in AA_ALPHABET)


def _truthy_annotation(series: pd.Series) -> pd.Series:
    text = series.astype("string").fillna("").str.strip().str.lower()
    numeric = pd.to_numeric(series, errors="coerce")
    return text.isin({"true", "t", "yes", "y", "1", "pass", "passed"}) | numeric.eq(1)


def _logo_weights(df: pd.DataFrame) -> pd.Series:
    if "max_freq" in df.columns:
        weights = pd.to_numeric(df["max_freq"], errors="coerce")
    else:
        freq_cols = [col for col in df.columns if str(col).startswith("freq ")]
        count_cols = [col for col in df.columns if str(col).startswith("count ")]
        if freq_cols:
            weights = df[freq_cols].apply(pd.to_numeric, errors="coerce").max(axis=1)
        elif count_cols:
            weights = df[count_cols].apply(pd.to_numeric, errors="coerce").max(axis=1)
        elif "rank" in df.columns:
            rank = pd.to_numeric(df["rank"], errors="coerce").clip(lower=1)
            weights = 1 / rank
        else:
            weights = pd.Series(1.0, index=df.index)

    weights = weights.fillna(0).clip(lower=0)
    if float(weights.sum()) <= 0:
        return pd.Series(1.0, index=df.index)
    return weights


def _sequence_logo_matrix(sequences: list[str], weights: pd.Series) -> pd.DataFrame:
    if not sequences:
        return pd.DataFrame()

    length = len(sequences[0])
    aa_to_idx = {aa: idx for idx, aa in enumerate(AA_ALPHABET)}
    counts = np.zeros((length, len(AA_ALPHABET)), dtype=float)

    for seq, weight in zip(sequences, weights):
        if len(seq) != length:
            continue
        for pos, aa in enumerate(seq):
            idx = aa_to_idx.get(aa)
            if idx is not None:
                counts[pos, idx] += float(weight)

    row_sums = counts.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        probs = np.divide(counts, row_sums, out=np.zeros_like(counts), where=row_sums > 0)
        entropy = -np.nansum(np.where(probs > 0, probs * np.log2(probs), 0), axis=1)
    information = np.log2(len(AA_ALPHABET)) - entropy
    heights = probs * information[:, None]
    return pd.DataFrame(heights, columns=list(AA_ALPHABET), index=range(1, length + 1))


def _imgt_position_sort_key(label: str, region: str) -> tuple:
    match = re.match(r"^([A-Za-z]+)(\d+)([A-Za-z]*)$", str(label))
    if not match:
        return (9, str(label))

    chain_type, number_text, letter = match.groups()
    number = int(number_text)
    chain_order = {"H": 0, "K": 1, "L": 1}.get(chain_type[:1].upper(), 8)
    region = str(region).upper()

    if letter:
        letter_rank = ord(letter[0].upper())
    else:
        letter_rank = 0

    if region in {"CDR3", "CHAIN"} and chain_type[:1].upper() == "H" and number >= 112:
        insertion_rank = 0 if letter else 1
        letter_rank = -letter_rank
    else:
        insertion_rank = 1 if letter else 0

    return (chain_order, number, insertion_rank, letter_rank, str(label))


def _imgt_centered_positions(
    *,
    chain_type: str,
    start: int,
    end: int,
    length: int,
    insertion_number: int | None = None,
) -> list[str]:
    if length <= 0:
        return []

    base_len = end - start + 1
    if length <= base_len:
        left_n = (length + 1) // 2
        right_n = length // 2
        left = [f"{chain_type}{pos}" for pos in range(start, start + left_n)]
        right = [f"{chain_type}{pos}" for pos in range(end - right_n + 1, end + 1)]
        return left + right

    labels = [f"{chain_type}{pos}" for pos in range(start, end + 1)]
    insertion_number = insertion_number or (start + base_len // 2 - 1)
    insertion_count = length - base_len
    insertions = [f"{chain_type}{insertion_number}{chr(ord('A') + idx)}" for idx in range(insertion_count)]
    insert_after = labels.index(f"{chain_type}{insertion_number}") + 1
    return labels[:insert_after] + insertions + labels[insert_after:]


def _imgt_cdr3_positions(length: int, chain_type: str = "H") -> list[str]:
    if length <= 0:
        return []

    if length <= 13:
        left_n = (length + 1) // 2
        right_n = length // 2
        left = [f"{chain_type}{pos}" for pos in range(105, 105 + left_n)]
        right = [f"{chain_type}{pos}" for pos in range(118 - right_n, 118)]
        return left + right

    extra = length - 13
    left_extra = extra // 2
    right_extra = extra - left_extra
    left_insertions = [f"{chain_type}111{chr(ord('A') + idx)}" for idx in range(left_extra)]
    right_insertions = [f"{chain_type}112{chr(ord('A') + idx)}" for idx in reversed(range(right_extra))]
    return (
        [f"{chain_type}{pos}" for pos in range(105, 112)]
        + left_insertions
        + right_insertions
        + [f"{chain_type}{pos}" for pos in range(112, 118)]
    )


def _imgt_positions_for_region(region: str, length: int, chain_type: str = "H") -> list[str]:
    region = str(region).upper()
    if region == "CDR1":
        return _imgt_centered_positions(chain_type=chain_type, start=27, end=38, length=length, insertion_number=32)
    if region == "CDR2":
        return _imgt_centered_positions(chain_type=chain_type, start=56, end=65, length=length, insertion_number=60)
    if region == "CDR3":
        return _imgt_cdr3_positions(length, chain_type=chain_type)

    fixed_positions = {
        "FR1": [*range(1, 10), *range(11, 27)],
        "FR2": list(range(39, 56)),
        "FR3": [*range(66, 73), *range(74, 105)],
        "FR4": list(range(118, 129)),
    }
    full_ranges = {
        "FR1": (1, 26),
        "FR2": (39, 55),
        "FR3": (66, 104),
        "FR4": (118, 128),
    }
    base = fixed_positions.get(region)
    if base and len(base) == length:
        return [f"{chain_type}{pos}" for pos in base]

    start_end = full_ranges.get(region)
    if start_end:
        start, end = start_end
        if length == end - start + 1:
            return [f"{chain_type}{pos}" for pos in range(start, end + 1)]
        return _imgt_centered_positions(chain_type=chain_type, start=start, end=end, length=length)

    return [str(idx) for idx in range(1, length + 1)]


def _positioned_from_region_sequence(sequence: str, region: str, chain_type: str = "H") -> list[tuple[str, str]]:
    seq = _clean_logo_sequence(sequence)
    positions = _imgt_positions_for_region(region, len(seq), chain_type=chain_type)
    if len(positions) != len(seq):
        return []
    return list(zip(positions, seq))


def _region_sequence_from_row(row: pd.Series, region: str) -> str:
    region = str(region).upper()
    if region in row.index:
        seq = _clean_logo_sequence(row.get(region, ""))
        if seq:
            return seq
    if region == "CDR3" and "cdr3_aa" in row.index:
        return _clean_logo_sequence(row.get("cdr3_aa", ""))
    return ""


def _find_region_start(hseq: str, region_seq: str, start: int) -> int:
    if not hseq or not region_seq:
        return -1
    pos = hseq.find(region_seq, start)
    if pos >= 0:
        return pos
    if region_seq.startswith("C"):
        pos = hseq.find(region_seq[1:], start)
        if pos >= 0:
            return pos
    return -1


def _hseq_region_sequences(row: pd.Series) -> dict[str, str]:
    hseq = _clean_logo_sequence(row.get("HSEQ", ""))
    if not hseq:
        return {}

    cdr1 = _region_sequence_from_row(row, "CDR1")
    cdr2 = _region_sequence_from_row(row, "CDR2")
    cdr3 = _region_sequence_from_row(row, "CDR3")
    if not cdr3:
        return {}

    if cdr1 and cdr2:
        cdr1_start = _find_region_start(hseq, cdr1, 0)
        cdr2_start = _find_region_start(hseq, cdr2, cdr1_start + len(cdr1) if cdr1_start >= 0 else 0)
        cdr3_start = _find_region_start(hseq, cdr3, cdr2_start + len(cdr2) if cdr2_start >= 0 else 0)
        if min(cdr1_start, cdr2_start, cdr3_start) >= 0:
            cdr3_hseq = hseq[cdr3_start:cdr3_start + len(cdr3)]
            if cdr3_hseq != cdr3 and cdr3.startswith("C") and hseq[cdr3_start:cdr3_start + len(cdr3) - 1] == cdr3[1:]:
                cdr3 = cdr3[1:]
            return {
                "FR1": hseq[:cdr1_start],
                "CDR1": cdr1,
                "FR2": hseq[cdr1_start + len(cdr1):cdr2_start],
                "CDR2": cdr2,
                "FR3": hseq[cdr2_start + len(cdr2):cdr3_start],
                "CDR3": cdr3,
                "FR4": hseq[cdr3_start + len(cdr3):],
            }

    cdr3_start = _find_region_start(hseq, cdr3, 0)
    if cdr3_start >= 0:
        if hseq[cdr3_start:cdr3_start + len(cdr3)] != cdr3 and cdr3.startswith("C"):
            cdr3 = cdr3[1:]
        return {"CDR3": cdr3}
    return {}


def _positioned_from_hseq(row: pd.Series, sequence_col: str) -> list[tuple[str, str]]:
    region = str(sequence_col).upper()
    hseq_regions = _hseq_region_sequences(row)
    if not hseq_regions:
        return []

    if region in {"CDR1", "CDR2", "CDR3"}:
        seq = hseq_regions.get(region, "")
        return _positioned_from_region_sequence(seq, region) if seq else []

    if region in {"HSEQ", "CHAIN"}:
        region_cols = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
        if not all(col in hseq_regions and hseq_regions[col] for col in region_cols):
            return []
        positioned: list[tuple[str, str]] = []
        for col in region_cols:
            region_positioned = _positioned_from_region_sequence(hseq_regions[col], col)
            if not region_positioned:
                return []
            positioned.extend(region_positioned)
        return positioned
    return []


def _positioned_chain_from_anarci_regions(row: pd.Series, chain_type: str = "H") -> list[tuple[str, str]]:
    region_cols = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
    if not all(col in row.index for col in region_cols):
        return []

    positioned: list[tuple[str, str]] = []
    for region in region_cols:
        seq = _clean_logo_sequence(row.get(region, ""))
        if not seq:
            return []
        region_positioned = _positioned_from_region_sequence(seq, region, chain_type=chain_type)
        if not region_positioned:
            return []
        positioned.extend(region_positioned)
    return positioned


def _imgt_positioned_region(row: pd.Series, sequence_col: str, *, prefer_hseq: bool = False) -> list[tuple[str, str]]:
    if prefer_hseq:
        positioned = _positioned_from_hseq(row, sequence_col)
        if positioned:
            return positioned

    region = str(sequence_col).upper()
    if region == "CDR3" and "CDR3" not in row.index and "cdr3_aa" in row.index:
        return _positioned_from_region_sequence(row.get("cdr3_aa", ""), region)
    if region in {"CDR1", "CDR2", "CDR3"} and region in row.index:
        return _positioned_from_region_sequence(row.get(region, ""), region)
    if region in {"CHAIN", "HSEQ"}:
        return _positioned_chain_from_anarci_regions(row)
    return []


def _positioned_alignment_key(positioned: list[tuple[str, str]]) -> str:
    return ";".join(f"{position}:{aa}" for position, aa in positioned)


def _sequence_logo_matrix_from_positions(
    positioned_sequences: list[list[tuple[str, str]]],
    weights: pd.Series,
    *,
    region: str,
) -> pd.DataFrame:
    if not positioned_sequences:
        return pd.DataFrame()

    positions = sorted(
        {position for positioned in positioned_sequences for position, _aa in positioned},
        key=lambda label: _imgt_position_sort_key(label, region),
    )
    if not positions:
        return pd.DataFrame()

    position_to_idx = {position: idx for idx, position in enumerate(positions)}
    aa_to_idx = {aa: idx for idx, aa in enumerate(AA_ALPHABET)}
    counts = np.zeros((len(positions), len(AA_ALPHABET)), dtype=float)

    for positioned, weight in zip(positioned_sequences, weights):
        for position, aa in positioned:
            pos_idx = position_to_idx.get(position)
            aa_idx = aa_to_idx.get(aa)
            if pos_idx is not None and aa_idx is not None:
                counts[pos_idx, aa_idx] += float(weight)

    row_sums = counts.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        probs = np.divide(counts, row_sums, out=np.zeros_like(counts), where=row_sums > 0)
        entropy = -np.nansum(np.where(probs > 0, probs * np.log2(probs), 0), axis=1)
    information = np.log2(len(AA_ALPHABET)) - entropy
    heights = probs * information[:, None]
    return pd.DataFrame(heights, columns=list(AA_ALPHABET), index=positions)


def _draw_logo_letter(ax, letter: str, x: float, y: float, height: float) -> None:
    if height <= 0:
        return
    path = TextPath((0, 0), letter, size=1, prop=LOGO_FONT)
    bbox = path.get_extents()
    if bbox.width <= 0 or bbox.height <= 0:
        return
    width = 0.86
    transform = (
        Affine2D()
        .translate(-bbox.x0, -bbox.y0)
        .scale(width / bbox.width, height / bbox.height)
        .translate(x + (1 - width) / 2, y)
        + ax.transData
    )
    patch = PathPatch(path, transform=transform, facecolor=AA_COLORS.get(letter, "#333333"), lw=0)
    ax.add_patch(patch)


def _compact_imgt_tick_label(label) -> str:
    text = str(label)
    return re.sub(r"^[A-Za-z]+(?=\d)", "", text)


def _logo_tick_indices(length: int, alignment_mode: str) -> tuple[list[int], int]:
    if length <= 0:
        return [], 1
    if alignment_mode == "imgt":
        if length > 100:
            step = 10
        elif length > 50:
            step = 5
        elif length > 28:
            step = 3
        elif length > 14:
            step = 2
        else:
            step = 1
    else:
        if length > 80:
            step = 10
        elif length > 40:
            step = 5
        else:
            step = 1

    indices = set(range(0, length, step))
    indices.add(0)
    indices.add(length - 1)
    return sorted(indices), step


def _draw_sequence_logo(
    matrix: pd.DataFrame,
    *,
    target: str,
    sequence_label: str,
    sequence_length: int,
    sequence_count: int,
    alignment_mode: str = "dominant_length",
):
    fig_width = max(8.5, min(30, sequence_length * 0.42))
    fig, ax = plt.subplots(figsize=(fig_width, 4.8))
    max_bits = float(matrix.sum(axis=1).max()) if not matrix.empty else 0
    y_limit = max(1.0, min(np.log2(len(AA_ALPHABET)) * 1.12, max_bits * 1.15))

    for x_idx, (_position, row) in enumerate(matrix.iterrows()):
        y = 0.0
        for aa, height in row[row > 0].sort_values().items():
            _draw_logo_letter(ax, aa, x_idx, y, float(height))
            y += float(height)

    ax.set_xlim(0, len(matrix))
    ax.set_ylim(0, y_limit)
    tick_indices, tick_step = _logo_tick_indices(len(matrix), alignment_mode)
    tick_positions = np.asarray(tick_indices, dtype=float)
    tick_labels = [matrix.index.tolist()[idx] for idx in tick_indices]
    if alignment_mode == "imgt":
        tick_labels = [_compact_imgt_tick_label(label) for label in tick_labels]
    ax.set_xticks(tick_positions + 0.5)
    ax.set_xticklabels(
        tick_labels,
        fontsize=8 if len(matrix) <= 40 else 7,
        rotation=45 if alignment_mode == "imgt" and tick_step > 1 else 0,
        ha="right" if alignment_mode == "imgt" and tick_step > 1 else "center",
    )
    if alignment_mode == "imgt":
        suffix = " (labels thinned)" if tick_step > 1 else ""
        ax.set_xlabel(f"{sequence_label} IMGT-aligned position{suffix}")
    else:
        ax.set_xlabel(f"{sequence_label} position, dominant length {sequence_length}")
    ax.set_ylabel("Information (bits)")
    ax.set_title(f"{target}: {sequence_label} Sequence Logo ({sequence_count} unique {sequence_label}s)")
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    sns.despine(ax=ax)
    return fig


def _is_vhh_logo_source(df: pd.DataFrame) -> bool:
    for col in ["library_type", "library"]:
        if col in df.columns:
            values = df[col].dropna().astype(str).str.lower()
            if values.str.contains("vhh|nanobody").any():
                return True
    if "vl_scaffold" not in df.columns:
        return True
    light_values = df["vl_scaffold"].dropna().astype(str).str.strip()
    return light_values.empty or light_values.isin(["", "UNK", "nan", "None"]).all()


def _logo_sequence_series(df: pd.DataFrame, sequence_col: str, *, prefer_hseq: bool = False) -> pd.Series:
    if prefer_hseq and "HSEQ" in df.columns:
        hseq = df["HSEQ"].apply(_clean_logo_sequence)
        if sequence_col.upper() == "HSEQ":
            return hseq

    if sequence_col == "CDR3":
        source_col = "CDR3" if "CDR3" in df.columns else "cdr3_aa" if "cdr3_aa" in df.columns else None
        if source_col is None:
            return pd.Series("", index=df.index)
        return df[source_col].apply(_clean_logo_sequence)

    if sequence_col in {"CHAIN", "HSEQ"}:
        if prefer_hseq and "HSEQ" in df.columns:
            hseq = df["HSEQ"].apply(_clean_logo_sequence)
            if hseq.str.len().gt(0).any():
                return hseq

        if "CHAIN" in df.columns:
            chain = df["CHAIN"].apply(_clean_logo_sequence)
        else:
            chain = pd.Series("", index=df.index)

        region_cols = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
        if all(col in df.columns for col in region_cols):
            reconstructed = df[region_cols].fillna("").astype(str).agg("".join, axis=1).apply(_clean_logo_sequence)
            chain = chain.where(chain.str.len() > 0, reconstructed)
        return chain

    if sequence_col not in df.columns:
        return pd.Series("", index=df.index)
    return df[sequence_col].apply(_clean_logo_sequence)


def _anarci_logo_mask(df: pd.DataFrame) -> tuple[pd.Series | None, str]:
    if "anarci_anno" in df.columns:
        return _truthy_annotation(df["anarci_anno"]), "anarci_anno"

    if "CHAIN" in df.columns:
        chain_text = df["CHAIN"].astype("string").fillna("").str.strip().str.lower()
        chain_seq = df["CHAIN"].apply(_clean_logo_sequence)
        mask = (chain_seq.str.len() > 0) & (chain_text != "not fully annotated")
        if mask.any():
            return mask, "CHAIN"

    cdr_cols = ["CDR1", "CDR2", "CDR3"]
    if all(col in df.columns for col in cdr_cols):
        mask = pd.Series(True, index=df.index)
        for col in cdr_cols:
            mask &= df[col].apply(_clean_logo_sequence).str.len() > 0
        if mask.any():
            return mask, "CDR1/CDR2/CDR3"

    return None, ""


def _build_sequence_logo_entry(
    df: pd.DataFrame,
    *,
    target: str,
    source: str,
    source_label: str,
    sequence_col: str,
    sequence_label: str,
    plots_dir: Path,
    prefer_hseq: bool = False,
) -> tuple[dict | None, dict | None]:
    logo_df = df.copy()
    logo_df["_logo_seq"] = _logo_sequence_series(logo_df, sequence_col, prefer_hseq=prefer_hseq)
    if prefer_hseq:
        if "HSEQ" not in logo_df.columns:
            print(f"Warning: {source_label} {sequence_label} sequence logo skipped for {target}: missing HSEQ column")
            return None, None
        logo_df["_logo_hseq"] = logo_df["HSEQ"].apply(_clean_logo_sequence)
        logo_df = logo_df[logo_df["_logo_hseq"].str.len() > 0].copy()
    logo_df = logo_df[logo_df["_logo_seq"].str.len() > 0].copy()
    if logo_df.empty:
        return None, None
    input_sequences = int(len(logo_df))

    anarci_mask, anarci_filter_source = _anarci_logo_mask(logo_df)
    if anarci_mask is None:
        print(
            f"Warning: {source_label} {sequence_label} sequence logo skipped for {target}: "
            "missing anarci_anno column and ANARCI-derived CDR/CHAIN fields"
        )
        return None, None
    if "max_freq" not in logo_df.columns:
        print(f"Warning: {source_label} {sequence_label} sequence logo skipped for {target}: missing max_freq column")
        return None, None

    logo_df["max_freq"] = pd.to_numeric(logo_df["max_freq"], errors="coerce").fillna(0)
    logo_df = logo_df[
        anarci_mask
        & (logo_df["max_freq"] > SEQUENCE_LOGO_MIN_MAX_FREQ)
    ].copy()
    if logo_df.empty:
        print(
            f"Warning: {source_label} {sequence_label} sequence logo skipped for {target}: no sequences after "
            f"anarci_anno=true and max_freq > {_format_threshold(SEQUENCE_LOGO_MIN_MAX_FREQ)} filtering"
        )
        return None, None
    filtered_clone_rows = int(len(logo_df))

    logo_df = (
        logo_df.sort_values("max_freq", ascending=False)
        .drop_duplicates(subset=["_logo_seq"], keep="first")
        .copy()
    )
    filtered_sequences = int(len(logo_df))

    alignment_mode = "dominant_length"
    dominant_length = None
    selected = logo_df.copy()
    sequence_length = 0
    matrix = pd.DataFrame()

    if sequence_col.upper() in {"CDR1", "CDR2", "CDR3", "CHAIN"}:
        logo_df["_imgt_positioned"] = logo_df.apply(
            lambda row: _imgt_positioned_region(row, sequence_col, prefer_hseq=prefer_hseq),
            axis=1,
        )
        aligned = logo_df[logo_df["_imgt_positioned"].map(bool)].copy()
        if not aligned.empty:
            aligned["_logo_alignment_key"] = aligned["_imgt_positioned"].apply(_positioned_alignment_key)
            selected = (
                aligned.sort_values("max_freq", ascending=False)
                .drop_duplicates(subset=["_logo_alignment_key"], keep="first")
                .copy()
            )
            weights = _logo_weights(selected)
            matrix = _sequence_logo_matrix_from_positions(
                selected["_imgt_positioned"].tolist(),
                weights,
                region=sequence_col,
            )
            if not matrix.empty:
                alignment_mode = "imgt"
                sequence_length = len(matrix)

    if matrix.empty and prefer_hseq:
        print(
            f"Warning: {source_label} {sequence_label} sequence logo skipped for {target}: "
            "no HSEQ rows could be aligned to ANARCI/IMGT positions"
        )
        return None, None

    if matrix.empty:
        length_counts = logo_df["_logo_seq"].str.len().value_counts()
        dominant_length = int(length_counts.idxmax())
        selected = logo_df[logo_df["_logo_seq"].str.len() == dominant_length].copy()
        if selected.empty:
            return None, None

        sequences = selected["_logo_seq"].tolist()
        weights = _logo_weights(selected)
        matrix = _sequence_logo_matrix(sequences, weights)
        if matrix.empty:
            return None, None
        sequence_length = dominant_length

    fig = _draw_sequence_logo(
        matrix,
        target=target,
        sequence_label=sequence_label,
        sequence_length=sequence_length,
        sequence_count=len(selected),
        alignment_mode=alignment_mode,
    )
    filename_suffix = "imgt" if alignment_mode == "imgt" else f"len{dominant_length}"
    plot_path = save_fig_plot(
        fig,
        plots_dir,
        f"sequence_logo_{source}_{target}_{sequence_col}_{filename_suffix}.png",
    )
    plt.close(fig)

    fraction = round(100 * len(selected) / filtered_sequences, 2)
    entry = {
        "target": target,
        "source": source,
        "source_label": source_label,
        "sequence_col": sequence_col,
        "sequence_label": sequence_label,
        "sequence_source": "HSEQ" if prefer_hseq and alignment_mode == "imgt" else "annotation_columns",
        "path": plot_path,
        "sequence_length": sequence_length,
        "dominant_length": dominant_length if dominant_length is not None else "",
        "alignment_mode": alignment_mode,
        "sequences_used": int(len(selected)),
        "total_sequences": filtered_sequences,
        "fraction_pct": fraction,
        "input_sequences": input_sequences,
        "filtered_clone_rows": filtered_clone_rows,
        "anarci_filter": "true",
        "anarci_filter_source": anarci_filter_source,
        "min_max_freq": SEQUENCE_LOGO_MIN_MAX_FREQ,
    }
    summary_row = {key: value for key, value in entry.items() if key != "path"}
    return entry, summary_row


def _read_clone_table(folder: Path, target: str) -> pd.DataFrame:
    clone_file = folder / f"{target}_clones.csv"
    if not clone_file.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(clone_file)
    except Exception as exc:
        print(f"Warning: could not read clones for sequence logo: {clone_file.name}: {exc}")
        return pd.DataFrame()


def _lead_logo_dataframe(folder: Path, target: str, final_df: pd.DataFrame) -> pd.DataFrame:
    required = {"anarci_anno", "max_freq"}
    if required.issubset(final_df.columns):
        return final_df

    clone_df = _read_clone_table(folder, target)
    if clone_df.empty:
        return final_df
    if "LEAD" in clone_df.columns:
        lead_df = clone_df[_truthy_annotation(clone_df["LEAD"])].copy()
        if not lead_df.empty:
            print(f"Using LEAD rows from {target}_clones.csv for lead sequence logos")
            return lead_df
    return final_df


def make_sequence_logo_plots(folder: Path, plots_dir: Path) -> list[dict]:
    for path in plots_dir.glob("sequence_logo_*.png"):
        try:
            path.unlink()
        except Exception as exc:
            print(f"Warning: could not remove stale sequence-logo plot {path.name}: {exc}")
    summary_path = folder / "sequence_logo_summary.csv"
    if summary_path.exists():
        try:
            summary_path.unlink()
        except Exception as exc:
            print(f"Warning: could not remove stale sequence-logo summary: {exc}")

    logo_entries: list[dict] = []
    summary_rows: list[dict] = []
    for final_file in sorted((folder / "by_protein").glob("*_final_leads.xlsx")):
        target = final_file.stem.replace("_final_leads", "")
        try:
            df = pd.read_excel(final_file)
        except Exception as exc:
            print(f"Warning: could not read final leads for sequence logo: {final_file.name}: {exc}")
            continue

        lead_logo_df = _lead_logo_dataframe(folder, target, df)
        regions = ["CDR1", "CDR2", "CDR3"] if _is_vhh_logo_source(lead_logo_df) else ["CDR3"]
        for region in regions:
            entry, summary_row = _build_sequence_logo_entry(
                lead_logo_df,
                target=target,
                source="lead",
                source_label="Lead HSEQ",
                sequence_col=region,
                sequence_label=region,
                plots_dir=plots_dir,
                prefer_hseq=True,
            )
            if entry:
                logo_entries.append(entry)
                summary_rows.append(summary_row)

    for clone_file in sorted(folder.glob("*_clones.csv")):
        target = clone_file.stem.replace("_clones", "")
        df = _read_clone_table(folder, target)
        if df.empty:
            continue

        entry, summary_row = _build_sequence_logo_entry(
            df,
            target=target,
            source="all_clones",
            source_label="All Clone",
            sequence_col="CHAIN",
            sequence_label="Full Chain",
            plots_dir=plots_dir,
        )
        if entry:
            logo_entries.append(entry)
            summary_rows.append(summary_row)

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    return logo_entries


def make_plots(folder: Path, report_format: str = "both"):
    folder = Path(folder)
    plots_dir = folder / "plots"
    plots_dir.mkdir(exist_ok=True)

    sns.set_style("whitegrid")
    plt.rcParams["figure.figsize"] = (14, 8)

    # Load sample QC
    qc_path = folder / "sample_qc_table.csv"
    if not qc_path.exists():
        print("sample_qc_table.csv not found — skipping plots")
        return

    qc = pd.read_csv(qc_path)


    # Load vh_cdr3_prevalent for top prevalent table
    prevalent_path = folder / "vh_cdr3_prevalent.csv"
    prevalent_df = pd.DataFrame()
    if prevalent_path.exists():
        prevalent_df = pd.read_csv(prevalent_path)
        prevalent_df = prevalent_df.sort_values("prevalent_targets", ascending=False).head(50)
        print(f"Loaded {len(prevalent_df)} prevalent VH-CDR3 pairs for table")

    ml_summary = load_ml_summary(folder)
    summary = write_data_summary(folder, qc, prevalent_df, ml_summary)
    contamination_cols = [
        col for col in [
            "name",
            "library",
            "antigen",
            "round",
            "condition",
            "n_cross_target_clones",
            "n_cross_target_reads",
            "pct_cross_target_reads",
        ]
        if col in qc.columns
    ]
    contamination_table = qc[contamination_cols].copy() if contamination_cols else pd.DataFrame()
    if not contamination_table.empty:
        contamination_table = contamination_table.sort_values(
            [col for col in ["pct_cross_target_reads", "n_cross_target_reads"] if col in contamination_table.columns],
            ascending=False,
        )

    # Generate plot image files
    generate_rarefaction_plots(folder)
    cleanup_stale_clone_frequency_plots(plots_dir)
    cleanup_stale_delphi_report_plots(plots_dir)
    delphi_score_plots = make_delphi_score_plots(folder, plots_dir)
    biophysical_profile_plots = make_biophysical_profile_plots(folder, plots_dir)
    sequence_logo_entries = make_sequence_logo_plots(folder, plots_dir)
    target_diversity_entries = make_target_diversity_plots(qc, plots_dir)

    # Start HTML report
    html = f"""
    <html>
    <head>
        <title>IPI-NGS Antibody Discovery Pipeline Report</title>
        <meta charset="utf-8">
        <style>
            :root {{ --ink: #18212f; --muted: #617086; --line: #d8dee8; --panel: #ffffff; --soft: #f4f7fb; --accent: #176d7a; }}
            body {{ font-family: Arial, Helvetica, sans-serif; margin: 0; line-height: 1.55; background: #eef2f6; color: var(--ink); }}
            h1, h2, h3 {{ color: var(--ink); }}
            h1 {{ margin: 0; font-size: 2.2em; }}
            h2 {{ margin-top: 0; }}
            .container {{ max-width: 1440px; margin: 0 auto; background: var(--panel); min-height: 100vh; padding: 34px; }}
            .header {{ border-bottom: 2px solid var(--line); padding-bottom: 18px; margin-bottom: 22px; }}
            .meta {{ color: var(--muted); margin: 4px 0; }}
            .plot {{ margin: 34px 0; padding-top: 8px; page-break-inside: avoid; border-top: 1px solid var(--line); }}
            img {{ max-width: 100%; height: auto; border: 1px solid var(--line); border-radius: 6px; background: #fff; }}
            .caption {{ margin-top: 8px; color: var(--muted); font-size: 0.95em; }}
            .footer {{ text-align: center; margin-top: 50px; color: #888; font-size: 0.9em; }}
            .summary-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 10px; margin: 22px 0; }}
            .summary-item {{ border: 1px solid var(--line); border-radius: 6px; padding: 12px; background: var(--soft); }}
            .summary-value {{ display: block; font-size: 1.45em; font-weight: bold; color: var(--accent); }}
            .summary-label {{ display: block; font-size: 0.88em; color: var(--muted); }}
            .table-container {{ overflow-x: auto; margin: 18px 0; }}
            table {{ font-size: 11px; width: 100%; border-collapse: collapse; table-layout: auto; background: #fff; }}
            table th, table td {{ padding: 7px 8px; text-align: left; border: 1px solid var(--line); word-wrap: break-word; max-width: 280px; vertical-align: top; }}
            table th {{ background-color: #e9eef5; font-weight: bold; }}
            table tr:nth-child(even) {{ background-color: #f8fafc; }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>NGSAbDiscov Report</h1>
                <p class="meta"><strong>Institute for Protein Innovation</strong> antibody discovery pipeline</p>
                <p class="meta"><strong>Generated:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
                <p class="meta"><strong>Results folder:</strong> {folder}</p>
            </div>
            <div class="summary-grid">
                <div class="summary-item"><span class="summary-value">{summary["sample_count"]}</span><span class="summary-label">Samples</span></div>
                <div class="summary-item"><span class="summary-value">{summary["target_count"]}</span><span class="summary-label">Targets</span></div>
                <div class="summary-item"><span class="summary-value">{summary["merged_reads"]:,}</span><span class="summary-label">Merged reads</span></div>
                <div class="summary-item"><span class="summary-value">{summary["median_merged_pct"]:.1f}%</span><span class="summary-label">Median merge</span></div>
                <div class="summary-item"><span class="summary-value">{summary["global_leads"]:,}</span><span class="summary-label">Global leads</span></div>
                <div class="summary-item"><span class="summary-value">{summary["samples_with_cross_target_reads"]}</span><span class="summary-label">Samples with cross-target reads</span></div>
                <div class="summary-item"><span class="summary-value">{summary["total_cross_target_reads"]:,}</span><span class="summary-label">Cross-target reads</span></div>
                <div class="summary-item"><span class="summary-value">{summary["max_cross_target_reads_pct"]:.2f}%</span><span class="summary-label">Max cross-target reads</span></div>
                <div class="summary-item"><span class="summary-value">{summary["ml_pass_rate_pct"]:.1f}%</span><span class="summary-label">ML pass rate</span></div>
                <div class="summary-item"><span class="summary-value">{summary["ml_pass"]:,}/{summary["ml_total_rows"]:,}</span><span class="summary-label">ML pass rows</span></div>
            </div>
    """

    plot_count = 0
    sequencing_html, plot_count = make_sequencing_depth_section(qc, plots_dir, plot_count)
    html += sequencing_html

    if not contamination_table.empty:
        html += (
            '<div class="plot"><h2>Cross-Target Contamination Summary</h2>'
            '<p class="caption">Per-sample contamination columns from sample_qc_table.csv.</p>'
            f'{contamination_table.to_html(index=False, escape=False, border=0)}</div>'
        )
    else:
        html += (
            '<div class="plot"><h2>Cross-Target Contamination Summary</h2>'
            '<p>No contamination columns were available in sample_qc_table.csv.</p></div>'
        )

    # === Red Flag Table (after sequencing depth) ===
    plot_count += 1
    red_flags = qc[
        (qc["merged"] < 20000) |
        (qc["unique_cdr3"] < 500)
    ][["name", "antigen", "block", "round", "condition", "merged", "unique_cdr3"]]

    if not red_flags.empty:
        red_flags = red_flags.sort_values(["merged", "unique_cdr3"])
        table_html = red_flags.to_html(index=False, classes="red-flag")
        html += f'<div class="plot"><h2>{plot_count}. Red Flags: Low Depth or Diversity</h2><p class="caption">Samples with < 20,000 merged reads or < 500 unique CDR3</p>{table_html}</div>'
    else:
        html += f'<div class="plot"><h2>{plot_count}. Red Flags</h2><p>No samples below thresholds (merged<20000 or unique_cdr3<500) — all good!</p></div>'


    # === 2. Cross-Target Contamination (optimized for up to 96 samples) ===
    if "pct_cross_target_reads" in qc.columns:
        plot_count += 1
        fig, ax = plt.subplots(figsize=(max(20, len(qc) * 0.5), 10))

        # Sort by % descending
        qc_sorted = qc.sort_values("pct_cross_target_reads", ascending=False)

        sns.barplot(data=qc_sorted, x="name", y="pct_cross_target_reads", ax=ax, color="#1f77b4")
        ax.set_title("Cross-Target Contamination (% of reads)")
        ax.tick_params(axis='x', rotation=90, labelsize=14)
        ax.set_xlabel("Sample Name")
        ax.set_ylabel("% Cross-Target Reads")
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        plot_path = save_fig_plot(fig, plots_dir, "cross_target_contamination.png")
        html += plot_html_from_file(plot_path, f"{plot_count}. Cross-Target Contamination", "% of reads from cross-target VH-CDR3 (sorted by %, wide view for many samples)")
        plt.close(fig)

    # === 9. Top 50 Most Prevalent VH-CDR3 Pairs (optimized table) ===
    if not prevalent_df.empty:
        plot_count += 1

        # Custom table HTML with wrapper for horizontal scroll
        table_html = prevalent_df.to_html(index=False, escape=False, border=0)

        optimized_table = f'''
        <div class="table-container">
            <h3>Top 50 Most Prevalent VH-CDR3 Pairs</h3>
            {table_html}
        </div>
        '''

        html += f'<div class="plot"><h2>{plot_count}. Top 50 Most Prevalent VH-CDR3 Pairs</h2>{optimized_table}</div>'
    
    # === 3. Unique CDR3 by Antigen & Condition ===
    plot_count += 1
    fig, ax = plt.subplots()
    n_conditions = qc["condition"].nunique()
    palette = sns.color_palette("tab20" if n_conditions > 10 else "tab10", n_conditions)
    sns.barplot(data=qc, x="antigen", y="unique_cdr3", hue="condition", palette=palette, ax=ax)
    ax.set_title("Unique CDR3 Distribution by Antigen and Condition")
    ax.tick_params(axis='x', rotation=90)
    ax.legend(title="Condition", bbox_to_anchor=(1.05, 1), loc='upper left')
    plot_path = save_fig_plot(fig, plots_dir, "unique_cdr3_by_antigen_condition.png")
    html += plot_html_from_file(plot_path, f"{plot_count}. Unique CDR3 by Antigen & Condition")
    plt.close(fig)

    # === 4. Pre-selected Clone Distribution ===
    plot_count += 1
    clone_dir = folder
    lead_dir = folder / "by_protein"
    if clone_dir.exists():
        target_count = []
        for f in clone_dir.glob("*_clones.csv"):
            target = f.stem.replace("_clones", "")
            df = pd.read_csv(f)
            unique = df["cdr3_aa"].nunique()
            target_count.append({"antigen": target, "unique_cdr3": unique})
            
        if target_count:
            tc_df = pd.DataFrame(target_count).sort_values("unique_cdr3", ascending=False)
            fig, ax = plt.subplots()
            sns.barplot(data=tc_df, x="antigen", y="unique_cdr3", ax=ax)
            ax.set_title("Unique CDR3 by Antigen (combined)")
            ax.tick_params(axis='x', rotation=90)
            plot_path = save_fig_plot(fig, plots_dir, "preselected_clone_distribution.png")
            html += plot_html_from_file(plot_path, f"{plot_count}. Pre-selected Clones Distribution")
            plt.close(fig)

    for entry in sequence_logo_entries:
        plot_count += 1
        sequence_label = entry["sequence_label"]
        source_label = entry["source_label"]
        anarci_filter_text = (
            "anarci_anno=true"
            if entry.get("anarci_filter_source") == "anarci_anno"
            else f"ANARCI {entry.get('anarci_filter_source', 'annotation')} present"
        )
        if entry.get("alignment_mode") == "imgt":
            source_phrase = (
                "from final HSEQ constructions"
                if entry.get("sequence_source") == "HSEQ"
                else "from ANARCI-derived annotation columns"
            )
            logo_caption = (
                f"Weighted amino-acid information logo {source_phrase}, aligned by ANARCI/IMGT positions "
                f"after {anarci_filter_text} and max_freq > {_format_threshold(SEQUENCE_LOGO_MIN_MAX_FREQ)} filtering "
                f"({entry['sequence_length']} IMGT positions; {entry['sequences_used']}/{entry['total_sequences']} "
                f"unique {sequence_label}s, {entry['fraction_pct']}%)."
            )
        else:
            logo_caption = (
                f"Weighted amino-acid information logo for the dominant {sequence_label} length "
                f"after {anarci_filter_text} and max_freq > {_format_threshold(SEQUENCE_LOGO_MIN_MAX_FREQ)} filtering "
                f"({entry['sequence_length']} aa; {entry['sequences_used']}/{entry['total_sequences']} "
                f"unique {sequence_label}s, {entry['fraction_pct']}%)."
            )
        html += plot_html_from_file(
            entry["path"],
            f"{plot_count}. {entry['target']}: {source_label} {sequence_label} Sequence Logo",
            logo_caption,
        )

    # === Diversity plots by round ===
    for metric, display in _available_diversity_metrics(qc):
        if metric in qc.columns:
            plot_count += 1

            qc_plot = qc.copy()
            qc_plot[metric] = pd.to_numeric(qc_plot[metric], errors="coerce")
            qc_plot["_diversity_value"], axis_label = _diversity_plot_values(qc_plot[metric], metric, display)
            qc_plot = qc_plot.dropna(subset=["_diversity_value"])
            if qc_plot.empty:
                continue

            round_order = _round_order(qc_plot["round"]) if "round" in qc_plot.columns else None

            fig, ax = plt.subplots(figsize=(14, 8))
            n_conditions = qc_plot["condition"].nunique() if "condition" in qc_plot.columns else 1
            palette = sns.color_palette("tab20" if n_conditions > 10 else "tab10", n_conditions)

            sns.boxplot(data=qc_plot, x="round", y="_diversity_value", order=round_order, ax=ax, color="lightgray", fliersize=0)
            if "condition" in qc_plot.columns:
                sns.stripplot(data=qc_plot, x="round", y="_diversity_value", hue="condition",
                              palette=palette, jitter=True, alpha=0.9, ax=ax, order=round_order)
                ax.legend(title="Condition", bbox_to_anchor=(1.05, 1), loc='upper left')
            else:
                sns.stripplot(data=qc_plot, x="round", y="_diversity_value",
                              color="#1f77b4", jitter=True, alpha=0.9, ax=ax, order=round_order)

            ax.set_title(f"{display} Diversity by Round")
            ax.set_ylabel(axis_label)
            ax.set_xlabel("round")
            ax.tick_params(axis='x', rotation=45)

            plot_path = save_fig_plot(fig, plots_dir, f"{metric}_diversity_by_round.png")
            caption = "Box plots summarize samples by round; points show individual samples."
            if _is_inverse_simpson_metric(metric):
                caption += " Inverse Simpson values are shown as log10(1 + value)."
            html += plot_html_from_file(plot_path, f"{plot_count}. {display} Diversity by Round", caption)
            plt.close(fig)

    for entry in target_diversity_entries:
        plot_count += 1
        html += plot_html_from_file(entry["path"], f"{plot_count}. {entry['title']}", entry["caption"])

    # === Rarefaction curves ===
    rarefaction_files = sorted(plots_dir.glob("rarefaction_*.png"))
    for rare_file in rarefaction_files:
        plot_count += 1
        with open(rare_file, "rb") as img_file:
            img_base64 = base64.b64encode(img_file.read()).decode('utf-8')
        title = rare_file.stem.replace("rarefaction_", "").replace("_", " ")
        html += f'<div class="plot"><h2>{plot_count}. Rarefaction: {title}</h2><img src="data:image/png;base64,{img_base64}"></div>'

    if not ml_summary.empty:
        ml_table = ml_summary.copy()
        for col in ["ml_pass_rate_pct"]:
            if col in ml_table.columns:
                ml_table[col] = pd.to_numeric(ml_table[col], errors="coerce").round(2)
        html += (
            '<div class="plot"><h2>Machine Learning Summary</h2>'
            '<p class="caption">PASS/FAIL reflects rows with full sequences accepted into Delphi prediction input.</p>'
            f'{ml_table.to_html(index=False, escape=False, border=0)}</div>'
        )
    else:
        html += '<div class="plot"><h2>Machine Learning Summary</h2><p>No ML summary was available.</p></div>'

    for entry in delphi_score_plots:
        html += plot_html_from_file(
            entry["path"],
            entry["title"],
            entry["caption"],
        )

    for entry in biophysical_profile_plots:
        html += plot_html_from_file(
            entry["path"],
            entry["title"],
            entry["caption"],
        )


    # === Per-target fold change plots (with potential leads table below each plot) ===
    if lead_dir.exists():
        target_files = list(lead_dir.glob("*_final_leads.xlsx"))
        for target_file in target_files:
            target = target_file.stem.replace("_final_leads", "")
            df = pd.read_excel(target_file)
            if df.empty or "rank" not in df.columns:
                continue

            freq_cols = [col for col in df.columns if col.startswith("freq ")]
            freq_4 = [col for col in freq_cols if "4nM" in col]
            freq_20 = [col for col in freq_cols if "20nM" in col]
            freq_100 = [col for col in freq_cols if "100nM" in col]

            comparisons = [
                ("4nM", "20nM", freq_4, freq_20, "logFC_4_20"),
                ("20nM", "100nM", freq_20, freq_100, "logFC_20_100"),
            ]
            for numerator_label, denominator_label, numerator_cols, denominator_cols, logfc_col in comparisons:
                if not numerator_cols or not denominator_cols:
                    continue

                plot_count += 1
                ratio = (df[numerator_cols[0]] + 1e-8) / (df[denominator_cols[0]] + 1e-8)
                df[logfc_col] = np.log2(ratio)
                df["log_rank"] = np.log(pd.to_numeric(df["rank"], errors="coerce").clip(lower=1))
                df["mean_psr_sec_plot_score"] = _delphi_psr_sec_global_mean_score(df)
                df["global_delphi_score_pass"] = (
                    df["mean_psr_sec_plot_score"] > DELPHI_FOLD_CHANGE_SCORE_THRESHOLD
                )

                plot_df = (
                    df[[logfc_col, "log_rank", "global_delphi_score_pass"]]
                    .replace([np.inf, -np.inf], np.nan)
                    .dropna(subset=[logfc_col, "log_rank"])
                )
                fig, ax = plt.subplots(figsize=(16, 8.5))
                base = plot_df[~plot_df["global_delphi_score_pass"]]
                ax.scatter(
                    base[logfc_col],
                    base["log_rank"],
                    s=32,
                    color="#8da0ae",
                    alpha=0.55,
                    label=f"Delphi mean score <= {DELPHI_FOLD_CHANGE_SCORE_THRESHOLD:g} / no score",
                )
                group_df = plot_df[plot_df["global_delphi_score_pass"]]
                if not group_df.empty:
                    ax.scatter(
                        group_df[logfc_col],
                        group_df["log_rank"],
                        s=48,
                        color="#2ca25f",
                        alpha=0.92,
                        edgecolors="white",
                        linewidths=0.35,
                        label=f"Delphi mean score > {DELPHI_FOLD_CHANGE_SCORE_THRESHOLD:g}",
                    )
                ax.axvline(-0.6, color="gray", linestyle="--", linewidth=1.4)
                ax.axvline(0, color="#d0d0d0", linewidth=1.0)
                ax.axvline(0.6, color="gray", linestyle="--", linewidth=1.4)
                ax.axhline(np.log(20), color="gray", linestyle="--", linewidth=1.4)
                ax.set_title(f"{target}: Fold Change {numerator_label} vs {denominator_label}")
                ax.set_xlabel(f"log2({numerator_label} / {denominator_label})")
                ax.set_ylabel("log(Rank)")
                ax.grid(True, color="#d0d0d0", linewidth=0.8, alpha=0.9)
                ax.legend(loc="best")

                plot_path = save_fig_plot(
                    fig,
                    plots_dir,
                    f"fold_change_{target}_{numerator_label}_vs_{denominator_label}.png",
                )
                html += plot_html_from_file(
                    plot_path,
                    f"{plot_count}. {target}: Fold Change {numerator_label} vs {denominator_label}",
                    (
                        "Colored points use Delphi mean score: "
                        f"(PSR mean + SEC mean) / 2 > {DELPHI_FOLD_CHANGE_SCORE_THRESHOLD:g}."
                    ),
                )
                plt.close(fig)

                potential = df[
                    (df[logfc_col] > 0.6)
                    & (pd.to_numeric(df["rank"], errors="coerce") <= 20)
                ].sort_values(logfc_col, ascending=False).head(20)

                if not potential.empty:
                    table_html = _lead_label_table(potential).to_html(index=False, header=False, escape=False)
                    html += (
                        '<div class="table-container"><h3>Potential Leads '
                        f'({numerator_label}/{denominator_label} FC > 0.6 & rank <= 20)</h3>'
                        f'{table_html}</div>'
                    )
                else:
                    html += (
                        '<div class="table-container"><h3>Potential Leads '
                        f'({numerator_label}/{denominator_label} FC > 0.6 & rank <= 20)</h3>'
                        '<p>No potential leads in this quadrant.</p></div>'
                    )

    html += """
            <div class="footer">
                <p>Generated by NGSAbDiscov Pipeline</p>
            </div>
        </div>
    </body>
    </html>
    """

    # Save HTML report
    report_html = folder / "report.html"
    if report_format in ("html", "both"):
        with open(report_html, "w", encoding="utf-8") as f:
            f.write(html)
        print(f"\nStep 5 complete! HTML report saved: {report_html}")

    # Generate PDF
    if report_format in ("pdf", "both"):
        try:
            from weasyprint import HTML
            report_pdf = folder / "report.pdf"
            HTML(string=html).write_pdf(report_pdf)
            print(f"PDF report saved: {report_pdf}")
        except ImportError:
            print("weasyprint not installed — run: pip install weasyprint")
        except Exception as e:
            print(f"PDF generation failed: {e}")

    print("Open report.html in your browser — all plots are embedded!")
