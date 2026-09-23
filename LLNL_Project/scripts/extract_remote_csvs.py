#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import shlex
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


DEFAULT_WORKBOOK = Path("/Users/Hoan.Nguyen/Downloads/LLNL_2nd delivery_from_Deepash_2026-06-26 (1).xlsx")
DEFAULT_DEST = Path("LLNL_Project/results/extracted_selected_samples")
DEFAULT_REPORT_DIR = Path("LLNL_Project/results/selected_sample_reports")


@dataclass(frozen=True)
class AntigenRequest:
    antigen: str
    miseq: str
    block: str
    sample1: str
    sample2: str


def clean(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def load_requests(workbook: Path) -> list[AntigenRequest]:
    df = pd.read_excel(workbook, sheet_name=0)
    required = ["Antigen Name", "Discovery Block", "MiSeq", "NGS_sample1", "NGS_sample2"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Workbook is missing required column(s): {', '.join(missing)}")

    requests: list[AntigenRequest] = []
    seen: set[tuple[str, str, str, str, str]] = set()
    for _, row in df.iterrows():
        antigen = clean(row["Antigen Name"])
        miseq = clean(row["MiSeq"])
        block = clean(row["Discovery Block"])
        sample1 = clean(row["NGS_sample1"])
        sample2 = clean(row["NGS_sample2"])
        if not antigen or not miseq:
            continue
        key = (antigen, miseq, block, sample1, sample2)
        if key not in seen:
            requests.append(AntigenRequest(*key))
            seen.add(key)
    return requests


def remote_find(host: str, remote_root: str, miseq_values: list[str]) -> list[str]:
    remote_command = f"find {shlex.quote(remote_root)} -type f -name '*.csv.gz'"
    result = subprocess.run(["ssh", host, remote_command], check=True, capture_output=True, text=True)
    files = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    return [
        remote_file
        for remote_file in files
        if Path(remote_file).parent.name == "results"
        and any(f"{miseq.lower()}/results/" in remote_file.lower() for miseq in miseq_values)
    ]


def load_remote_paths_file(path: Path, miseq_values: list[str]) -> list[str]:
    files = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    return [
        remote_file
        for remote_file in files
        if Path(remote_file).parent.name == "results"
        and any(f"{miseq.lower()}/results/" in remote_file.lower() for miseq in miseq_values)
    ]


def _selected_samples(request: AntigenRequest) -> list[tuple[str, str]]:
    samples = []
    if request.sample1:
        samples.append(("NGS_sample1", request.sample1))
    if request.sample2:
        samples.append(("NGS_sample2", request.sample2))
    return samples


def _sample_match(filename: str, request: AntigenRequest, *, all_antigen_files: bool) -> tuple[str, str] | None:
    if request.antigen not in filename:
        return None
    if all_antigen_files:
        return ("all_antigen_files", "")
    if request.block and f"__{request.block}__" not in filename:
        return None
    for sample_column, sample in _selected_samples(request):
        if sample and sample in filename:
            return (sample_column, sample)
    return None


def match_files(
    requests: list[AntigenRequest],
    remote_files: list[str],
    *,
    all_antigen_files: bool = False,
) -> tuple[list[dict], list[dict]]:
    files_by_miseq: dict[str, list[str]] = defaultdict(list)
    for remote_file in remote_files:
        for request in requests:
            if f"{request.miseq.lower()}/" in remote_file.lower():
                files_by_miseq[request.miseq].append(remote_file)

    matched_rows: list[dict] = []
    missing_rows: list[dict] = []
    seen_paths: set[str] = set()

    for request in requests:
        row_matches: list[tuple[str, str, str]] = []
        for remote_file in files_by_miseq.get(request.miseq, []):
            match = _sample_match(Path(remote_file).name, request, all_antigen_files=all_antigen_files)
            if match:
                sample_column, sample = match
                row_matches.append((remote_file, sample_column, sample))

        row_matches = sorted(set(row_matches))
        if not row_matches:
            missing_rows.append(
                {
                    "antigen": request.antigen,
                    "miseq": request.miseq,
                    "block": request.block,
                    "ngs_sample1": request.sample1,
                    "ngs_sample2": request.sample2,
                }
            )
            continue

        for remote_file, sample_column, sample in row_matches:
            seen_paths.add(remote_file)
            matched_rows.append(
                {
                    "antigen": request.antigen,
                    "miseq": request.miseq,
                    "block": request.block,
                    "ngs_sample1": request.sample1,
                    "ngs_sample2": request.sample2,
                    "matched_sample_column": sample_column,
                    "matched_sample": sample,
                    "remote_path": remote_file,
                    "filename": Path(remote_file).name,
                }
            )

    matched_rows.sort(key=lambda item: (item["miseq"], item["antigen"], item["remote_path"]))
    return matched_rows, missing_rows


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_reports(report_dir: Path, matched_rows: list[dict], missing_rows: list[dict]) -> Path:
    report_dir.mkdir(parents=True, exist_ok=True)
    matched_path = report_dir / "matched_remote_csvs.csv"
    missing_path = report_dir / "missing_antigens.csv"
    paths_path = report_dir / "matched_remote_paths.txt"

    write_csv(
        matched_path,
        matched_rows,
        [
            "antigen",
            "miseq",
            "block",
            "ngs_sample1",
            "ngs_sample2",
            "matched_sample_column",
            "matched_sample",
            "remote_path",
            "filename",
        ],
    )
    write_csv(missing_path, missing_rows, ["antigen", "miseq", "block", "ngs_sample1", "ngs_sample2"])
    unique_paths = sorted({row["remote_path"] for row in matched_rows})
    paths_path.write_text("\n".join(unique_paths) + ("\n" if unique_paths else ""))
    return paths_path


def copy_matches(host: str, paths_file: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run(["rsync", "-av", "--files-from", str(paths_file), f"{host}:/", str(dest)], check=True)


def copy_matches_from_local(paths_file: Path, source_dir: Path, dest: Path) -> list[str]:
    dest.mkdir(parents=True, exist_ok=True)
    missing: list[str] = []
    paths = [line.strip() for line in paths_file.read_text().splitlines() if line.strip()]
    for remote_path in paths:
        relative_path = Path(remote_path.lstrip("/"))
        source = source_dir / relative_path
        target = dest / relative_path
        if not source.exists():
            missing.append(remote_path)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    return missing


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract LLNL csv.gz files from remote NGS result folders.")
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--remote-host", default="ipinode1")
    parser.add_argument("--remote-root", default="/alphafold/combio/NGS")
    parser.add_argument(
        "--remote-paths-file",
        type=Path,
        help="Use an existing newline-delimited remote path list instead of running ssh find.",
    )
    parser.add_argument("--dest", type=Path, default=DEFAULT_DEST)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--copy", action="store_true", help="Copy matched files after writing reports.")
    parser.add_argument(
        "--local-source-dir",
        type=Path,
        help="Copy from an existing local mirror of remote absolute paths instead of rsyncing from remote host.",
    )
    parser.add_argument(
        "--all-antigen-files",
        action="store_true",
        help="Match every direct results/*.csv.gz file for each antigen, ignoring NGS_sample1/NGS_sample2.",
    )
    args = parser.parse_args()

    requests = load_requests(args.workbook)
    miseq_values = sorted({request.miseq for request in requests})
    print(f"Loaded {len(requests)} antigen rows across MiSeq values: {', '.join(miseq_values)}")

    if args.remote_paths_file:
        remote_files = load_remote_paths_file(args.remote_paths_file, miseq_values)
    else:
        remote_files = remote_find(args.remote_host, args.remote_root, miseq_values)
    print(f"Remote files listed: {len(remote_files)}")

    match_label = "all antigen files" if args.all_antigen_files else "NGS_sample1/NGS_sample2 only"
    print(f"Match mode: {match_label}")

    matched_rows, missing_rows = match_files(requests, remote_files, all_antigen_files=args.all_antigen_files)
    paths_file = write_reports(args.report_dir, matched_rows, missing_rows)
    unique_match_count = len({row["remote_path"] for row in matched_rows})
    print(f"Matched remote files: {unique_match_count}")
    print(f"Antigen rows without matches: {len(missing_rows)}")
    print(f"Reports written to: {args.report_dir}")

    if args.copy:
        if args.local_source_dir:
            missing_local = copy_matches_from_local(paths_file, args.local_source_dir, args.dest)
            if missing_local:
                missing_path = args.report_dir / "missing_local_source_files.txt"
                missing_path.write_text("\n".join(missing_local) + "\n")
                raise FileNotFoundError(f"{len(missing_local)} matched files were missing from {args.local_source_dir}")
        else:
            copy_matches(args.remote_host, paths_file, args.dest)
        print(f"Copied matched files to: {args.dest}")
    else:
        print("Dry run only. Re-run with --copy to copy matched files.")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
