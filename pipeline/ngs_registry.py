from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


def _json_default(value: Any):
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _safe_read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.DataFrame()


def _safe_count_excel_rows(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        return len(pd.read_excel(path))
    except Exception:
        return 0


def _safe_count_excel_rows_many(paths: list[Path]) -> int:
    return sum(_safe_count_excel_rows(path) for path in paths)


def _registry_path(cfg: dict, run_dir: Path, fastq_folder: Path) -> Path:
    configured = cfg.get("registry", {}).get("path")
    if configured:
        return Path(configured)
    ngs_root = fastq_folder.parent.parent if fastq_folder.parent.parent.exists() else run_dir.parent
    return ngs_root / "ngs_registry.jsonl"


def build_run_manifest(
    cfg: dict,
    *,
    run_dir: Path,
    results_dir: Path,
    sample_sheet: Path,
    fastq_folder: Path,
) -> dict:
    run_dir = Path(run_dir)
    results_dir = Path(results_dir)
    sample_sheet = Path(sample_sheet)
    fastq_folder = Path(fastq_folder)

    qc = _safe_read_csv(results_dir / "sample_qc_table.csv")
    leads_path = results_dir / "leads.xlsx"
    ml_summary = _safe_read_csv(results_dir / "ml_summary.csv")
    clone_files = sorted(results_dir.glob("*_clones.csv"))
    final_lead_files = sorted((results_dir / "by_protein").glob("*_final_leads.xlsx"))
    global_leads = _safe_count_excel_rows(leads_path) or _safe_count_excel_rows_many(final_lead_files)

    targets = sorted({path.stem.replace("_clones", "") for path in clone_files})
    library = cfg.get("current_library", "")
    library_cfg = cfg.get("libraries", {}).get(library, {})
    library_counts = {}
    libraries = []
    library_types = []
    if not qc.empty and "library" in qc.columns:
        library_counts = qc["library"].fillna("").astype(str).value_counts().to_dict()
        libraries = sorted(k for k in library_counts if k)
        library_types = sorted(qc.get("library_type", pd.Series(dtype=str)).fillna("").astype(str).unique())
    elif library:
        libraries = [library]
        if library_cfg.get("library_type"):
            library_types = [library_cfg.get("library_type")]

    manifest = {
        "schema_version": "1.0",
        "run_id": run_dir.name,
        "run_name": run_dir.name,
        "project": cfg.get("registry", {}).get("project", "IPI"),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "library": library,
        "library_type": library_cfg.get("library_type", ""),
        "libraries": libraries,
        "library_types": library_types,
        "library_counts": library_counts,
        "paths": {
            "run_dir": run_dir,
            "fastq_folder": fastq_folder,
            "sample_sheet": sample_sheet,
            "results_dir": results_dir,
            "report_html": results_dir / "report.html",
            "report_pdf": results_dir / "report.pdf",
        },
        "sample_count": int(len(qc)) if not qc.empty else 0,
        "target_count": len(targets),
        "targets": targets,
        "read_summary": {},
        "lead_summary": {
            "global_leads": global_leads,
            "per_target_files": len(final_lead_files),
        },
        "ml_summary": {},
        "outputs": {
            "clone_files": [str(path) for path in clone_files],
            "final_lead_files": [str(path) for path in final_lead_files],
        },
    }

    if not qc.empty:
        for col in ["total", "merged", "unique_cdr3", "n_cross_target_reads"]:
            if col in qc.columns:
                manifest["read_summary"][f"{col}_sum"] = int(pd.to_numeric(qc[col], errors="coerce").fillna(0).sum())
        if "merged_pct" in qc.columns:
            manifest["read_summary"]["merged_pct_median"] = float(
                pd.to_numeric(qc["merged_pct"], errors="coerce").dropna().median()
            )
        if "pct_cross_target_reads" in qc.columns:
            manifest["read_summary"]["cross_target_reads_pct_max"] = float(
                pd.to_numeric(qc["pct_cross_target_reads"], errors="coerce").dropna().max()
            )

    if not ml_summary.empty:
        for col in ["total_rows", "ml_pass", "ml_fail"]:
            if col in ml_summary.columns:
                manifest["ml_summary"][f"{col}_sum"] = int(
                    pd.to_numeric(ml_summary[col], errors="coerce").fillna(0).sum()
                )
        if "ml_pass_rate_pct" in ml_summary.columns:
            manifest["ml_summary"]["pass_rate_pct_mean"] = float(
                pd.to_numeric(ml_summary["ml_pass_rate_pct"], errors="coerce").dropna().mean()
            )

    return manifest


def _upsert_jsonl(path: Path, record: dict, key: str = "run_id") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []
    if path.exists():
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            try:
                old_record = json.loads(line)
            except json.JSONDecodeError:
                continue
            if old_record.get(key) != record.get(key):
                records.append(old_record)
    records.append(record)
    path.write_text("\n".join(json.dumps(item, default=_json_default) for item in records) + "\n")


def update_registry(
    cfg: dict,
    *,
    run_dir: Path,
    results_dir: Path,
    sample_sheet: Path,
    fastq_folder: Path,
) -> dict:
    manifest = build_run_manifest(
        cfg,
        run_dir=run_dir,
        results_dir=results_dir,
        sample_sheet=sample_sheet,
        fastq_folder=fastq_folder,
    )

    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = results_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, default=_json_default) + "\n")
    print(f"Run manifest saved: {manifest_path}")

    if cfg.get("registry", {}).get("enabled", True):
        registry_path = _registry_path(cfg, Path(run_dir), Path(fastq_folder))
        _upsert_jsonl(registry_path, manifest)
        print(f"NGS registry updated: {registry_path}")

    return manifest
