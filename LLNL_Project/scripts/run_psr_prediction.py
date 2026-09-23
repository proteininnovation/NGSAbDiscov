#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pipeline.predict_ml import (  # noqa: E402
    CNN_PRED_PSR,
    Generate_Antiberta2_Embedding_HLSEQ,
    PredictionConsensus,
    XGBoost_PRED_PSR,
    _resolve_model_dir,
)


DEFAULT_INPUT_DIR = Path("LLNL_Project/results/final_tables")
DEFAULT_OUTPUT_DIR = Path("LLNL_Project/results/psr_predictions")
DEFAULT_REPORT_DIR = Path("LLNL_Project/results/psr_prediction_reports")
REQUIRED_COLUMNS = ["cdr3_aa", "vh_scaffold", "vl_scaffold", "HSEQ", "LSEQ", "count", "freq"]


def output_name(input_file: Path) -> str:
    return f"{input_file.stem}_psr_predicted.csv"


def valid_sequence_mask(df: pd.DataFrame) -> pd.Series:
    return (
        df["HSEQ"].fillna("").astype(str).str.strip().ne("")
        & df["LSEQ"].fillna("").astype(str).str.strip().ne("")
        & ~df["HSEQ"].fillna("").astype(str).str.contains("*", regex=False)
        & ~df["LSEQ"].fillna("").astype(str).str.contains("*", regex=False)
    )


def read_prediction_input(input_file: Path, max_rows: int | None = None) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(input_file)
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"{input_file.name} is missing required column(s): {', '.join(missing)}")

    input_rows = len(df)
    usable = valid_sequence_mask(df)
    df = df.loc[usable].copy()
    if max_rows is not None:
        df = df.head(max_rows).copy()

    return df, {
        "input_rows": input_rows,
        "usable_sequence_rows": int(usable.sum()),
        "dropped_sequence_rows": int((~usable).sum()),
    }


def score_dataframe(
    df: pd.DataFrame,
    *,
    model_dir: Path,
    batch_size: int,
    device: str,
    include_embeddings: bool,
    predictors: list[str],
) -> pd.DataFrame:
    scored = Generate_Antiberta2_Embedding_HLSEQ(df, batch_size=batch_size, device=device)

    errors: list[str] = []
    available_predictors = {
        "cnn": CNN_PRED_PSR,
        "xgboost": XGBoost_PRED_PSR,
    }
    for predictor_name in predictors:
        predictor = available_predictors[predictor_name]
        try:
            scored = predictor(scored, model_dir)
        except Exception as exc:
            errors.append(str(exc))
            print(f"Warning ({predictor_name}): {exc}")

    if "psr_cnn_pred_proba" not in scored.columns and "psr_xgboost_pred_proba" not in scored.columns:
        raise RuntimeError("All PSR predictors failed: " + " | ".join(errors))

    scored = PredictionConsensus(scored)
    if not include_embeddings:
        scored = scored.drop(columns=["antiberta2-cssp_emb", "HL"], errors="ignore")
    return scored


def write_summary(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "input_file",
        "output_file",
        "status",
        "input_rows",
        "usable_sequence_rows",
        "dropped_sequence_rows",
        "scored_rows",
        "psr_pass",
        "psr_fail",
        "psr_pass_rate_pct",
        "error",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run PSR/polyreactivity predictions for LLNL final CSV tables.")
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--model-dir", type=Path, default=Path("psr_model_ml"))
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--device", default="cpu")
    parser.add_argument(
        "--predictors",
        default="xgboost",
        help="Comma-separated PSR predictors to run: xgboost, cnn, or xgboost,cnn. Default avoids local Keras loading issues.",
    )
    parser.add_argument("--max-rows-per-file", type=int, default=None)
    parser.add_argument("--limit-files", type=int, default=None)
    parser.add_argument("--include-embeddings", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    model_dir = _resolve_model_dir(args.model_dir)
    predictors = [item.strip().lower() for item in args.predictors.split(",") if item.strip()]
    unknown_predictors = sorted(set(predictors) - {"cnn", "xgboost"})
    if unknown_predictors:
        raise ValueError(f"Unknown predictor(s): {', '.join(unknown_predictors)}")
    if not predictors:
        raise ValueError("At least one predictor is required.")

    input_files = sorted(args.input_dir.glob("*.csv"))
    if args.limit_files is not None:
        input_files = input_files[: args.limit_files]
    if not input_files:
        raise FileNotFoundError(f"No CSV inputs found in {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict] = []
    for input_file in input_files:
        output_file = args.output_dir / output_name(input_file)
        print(f"\nProcessing {input_file.name}")
        if output_file.exists() and not args.overwrite:
            print(f"Skipping existing output: {output_file}")
            continue

        row = {
            "input_file": str(input_file),
            "output_file": str(output_file),
            "status": "",
            "input_rows": 0,
            "usable_sequence_rows": 0,
            "dropped_sequence_rows": 0,
            "scored_rows": 0,
            "psr_pass": 0,
            "psr_fail": 0,
            "psr_pass_rate_pct": 0.0,
            "error": "",
        }
        try:
            df, stats = read_prediction_input(input_file, args.max_rows_per_file)
            row.update(stats)
            scored = score_dataframe(
                df,
                model_dir=model_dir,
                batch_size=args.batch_size,
                device=args.device,
                include_embeddings=args.include_embeddings,
                predictors=predictors,
            )
            scored.to_csv(output_file, index=False)
            row["status"] = "ok"
            row["scored_rows"] = len(scored)
            row["psr_pass"] = int((scored["psr_pred_mean"] == 1).sum())
            row["psr_fail"] = int((scored["psr_pred_mean"] == 0).sum())
            row["psr_pass_rate_pct"] = round(100 * row["psr_pass"] / row["scored_rows"], 2) if row["scored_rows"] else 0.0
            print(f"Saved: {output_file}")
            print(f"PASS: {row['psr_pass']} | FAIL: {row['psr_fail']}")
        except Exception as exc:
            row["status"] = "failed"
            row["error"] = str(exc)
            print(f"Failed: {exc}")
        summary_rows.append(row)

    summary_path = args.report_dir / "psr_prediction_summary.csv"
    write_summary(summary_path, summary_rows)
    print(f"\nSummary written: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
