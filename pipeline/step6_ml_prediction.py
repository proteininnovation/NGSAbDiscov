from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


def _backend_available(name: str, cfg: dict | None = None) -> bool:
    if name in ("delph", "delphi"):
        try:
            from .predict_ml import delph_import_available

            return delph_import_available((cfg or {}).get("ml", {}).get("delph", {}))
        except Exception:
            return importlib.util.find_spec("delph") is not None or importlib.util.find_spec("delphi") is not None
    if name == "ipi_psr":
        return True
    return False


def _choose_backend(requested: str, cfg: dict) -> str:
    requested = (requested or cfg.get("ml", {}).get("backend", "auto")).lower()
    if requested == "delphi":
        requested = "delph"
    if requested != "auto":
        return requested
    if _backend_available("delph", cfg):
        return "delph"
    return "ipi_psr"


def _candidate_files(folder: Path, input_globs: list[str]) -> list[Path]:
    files: list[Path] = []
    for pattern in input_globs:
        files.extend(sorted(folder.glob(pattern)))

    seen: set[Path] = set()
    unique_files: list[Path] = []
    for path in files:
        resolved = path.resolve()
        if resolved in seen:
            continue
        if "_backup" in path.stem or "_predicted" in path.stem:
            continue
        seen.add(resolved)
        unique_files.append(path)
    return unique_files


def _input_table_is_empty(path: Path) -> bool:
    """Return True when an ML input table has a header but no data rows."""
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, nrows=1).empty
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, nrows=1).empty
    return False


def _target_name_from_ml_input(path_value: str | None) -> str:
    if not path_value:
        return ""
    name = Path(path_value).stem
    for suffix in ["_final_leads_delphi_input", "_final_leads", "_delphi_input", "_delph_predicted"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def _write_ml_summary(folder: Path, results: list[dict]) -> Path:
    rows = []
    for result in results:
        input_path = result.get("source_input") or result.get("input") or result.get("final_output") or result.get("output")
        total_rows = int(result.get("rows") or 0)
        pass_rows = int(result.get("ml_rows") or 0)
        fail_rows = int(result.get("dropped_rows") if result.get("dropped_rows") is not None else max(total_rows - pass_rows, 0))
        rows.append(
            {
                "target": _target_name_from_ml_input(str(input_path) if input_path else ""),
                "status": result.get("status", ""),
                "backend": result.get("backend", ""),
                "total_rows": total_rows,
                "ml_pass": pass_rows,
                "ml_fail": fail_rows,
                "ml_pass_rate_pct": round(100 * pass_rows / total_rows, 2) if total_rows else 0.0,
                "prediction_columns": ";".join(map(str, result.get("merged_prediction_columns") or [])),
                "sequence_mode": result.get("sequence_mode", ""),
                "error": result.get("error", ""),
            }
        )
    summary = pd.DataFrame(rows)
    summary_path = folder / "ml_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"ML summary saved: {summary_path}")
    return summary_path


def run_ml_prediction(folder: Path, cfg: dict | None = None, backend: str | None = None) -> list[dict]:
    folder = Path(folder)
    cfg = cfg or {}
    ml_cfg = cfg.get("ml", {})
    fail_on_error = bool(ml_cfg.get("fail_on_error", False))
    selected_backend = _choose_backend(backend, cfg)
    input_globs = ml_cfg.get("input_globs", ["by_protein/*_final_leads.xlsx", "*_clones.csv"])

    clone_files = _candidate_files(folder, input_globs)
    if not clone_files:
        print("No ML input files found — Step 6 skipped.")
        return []

    print(f"Found {len(clone_files)} files for ML prediction")
    print(f"ML backend: {selected_backend}")

    results: list[dict] = []
    for input_file in clone_files:
        print(f"\nProcessing ML input: {input_file}")
        try:
            if _input_table_is_empty(input_file):
                result = {
                    "input": str(input_file),
                    "status": "skipped_empty",
                    "backend": selected_backend,
                    "rows": 0,
                    "ml_rows": 0,
                    "dropped_rows": 0,
                }
                print(f"  Skipped empty ML input: {input_file.name}")
                results.append(result)
                continue

            if selected_backend == "delph":
                from .predict_ml import predict_delph_on_clone

                result = predict_delph_on_clone(input_file, delph_config=ml_cfg.get("delph", {}))
            elif selected_backend == "ipi_psr":
                from .predict_ml import predict_psr_on_clone

                result = predict_psr_on_clone(
                    input_file,
                    model_dir=ml_cfg.get("model_dir"),
                    include_embeddings=bool(ml_cfg.get("include_embeddings", False)),
                    batch_size=int(ml_cfg.get("batch_size", 128)),
                    device=str(ml_cfg.get("device", "cpu")),
                )
            else:
                raise ValueError(f"Unknown ML backend: {selected_backend}")

            result["status"] = "ok"
            result["backend"] = selected_backend
            print(f"  Prediction completed: {result.get('output', '')}")
        except Exception as exc:
            result = {
                "input": str(input_file),
                "status": "failed",
                "backend": selected_backend,
                "error": str(exc),
            }
            print(f"  Prediction failed for {input_file.name}: {exc}")
            if fail_on_error:
                raise
        results.append(result)

    status_path = folder / "ml_prediction_status.json"
    with open(status_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
        fh.write("\n")
    print(f"\nStep 6 complete. Status saved: {status_path}")
    _write_ml_summary(folder, results)
    return results
