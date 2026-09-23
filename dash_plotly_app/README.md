# Dash/Plotly Clone Selection App

This folder contains the Codex-built interactive clone-selection app for
NGSAbDiscov results. It is intentionally separated from the core FASTQ/clone/ML
pipeline code.

## Purpose

The app reads an existing NGSAbDiscov results folder and helps select antibody
clones interactively. It does not rerun FASTQ processing, clone grouping, or ML
prediction.

Input data for one result folder:

```text
<results_folder>/by_protein/*_final_leads.xlsx
```

Input data for multi-MiSeq browsing:

```text
<results_root>/<MiSeq run>/<results folder>/by_protein/*_final_leads.xlsx
```

Main app file:

```text
dash_plotly_app/clone_selection_app.py
```

## Run Directly

From the NGSAbDiscov repository:

```bash
/usr/local/lib/NGS_pipeline/annaconda3/bin/python dash_plotly_app/clone_selection_app.py \
  --results-root /alphafold/combio/NGS \
  --host 127.0.0.1 \
  --port 8050 \
  --max-table-rows 500
```

## Run Through NGSAbDiscov CLI

```bash
/usr/local/lib/NGS_pipeline/annaconda3/bin/python __main__.py clone-app \
  --results-root /alphafold/combio/NGS \
  --host 127.0.0.1 \
  --port 8050 \
  --max-table-rows 500
```

## Public Team Server

Run this on `ipinode1` only when you want the app visible to other lab machines:

```bash
cd /alphafold/combio/software/NGSAbDiscov

nohup /usr/local/lib/NGS_pipeline/annaconda3/bin/python __main__.py clone-app \
  --results-root /alphafold/combio/NGS \
  --host 0.0.0.0 \
  --port 8051 \
  --max-table-rows 500 \
  > /alphafold/combio/NGS/<RUN_FOLDER>/clone_app_8051.log 2>&1 &

echo $! > /alphafold/combio/NGS/<RUN_FOLDER>/clone_app_8051.pid
```

Then open:

```text
http://ipinode1:8051
```

## Current Features

- MiSeq/result-folder selector across `/alphafold/combio/NGS`
- single-target lead-selection workflow
- max frequency, rank, CDR3, Delphi, PSR, and SEC filters
- default Delphi mean score cutoff of `0.3`
- ratio-based F_N/`2uM` control filtering: `selected_condition_freq / F_N_freq`
  with a default cutoff of `0.3`
- PSR-like background flag when `F_N_freq / selected_condition_freq > 3`
- fold-change filtering by explicit condition pairs such as `4nM / 20nM`, or manual sample columns
- condition-aware enrichment filters for available `100nM`, `20nM`, and `4nM`
  columns
- CDR3 diversity mode with one representative per greedy CDR3 cluster
- top 500 visible clone rows
- clone table above plots
- compact viewport-fitting layout with tabbed plots
- box/lasso selection on the fold-change plot updates the clone table; use
  `Show all filtered clones` to undo the plot selection
- selected-clone frequency profiles
- compact CDR3 distance tree, defaulting to top representative clusters
- sortable/filterable clone table
- CSV downloads for selected or filtered clones
- full `HSEQ` and `LSEQ` retained in downloads, while the visible table shows
  only sequence lengths for speed
