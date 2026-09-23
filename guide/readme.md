# NGSAbDiscov Run Guide

This guide is for running the NGS antibody discovery pipeline on `ipinode1`.
The current production code is:

```text
/alphafold/combio/software/NGSAbDiscov
```

Use the bundled Python environment:

```bash
/usr/local/lib/NGS_pipeline/annaconda3/bin/python
```

## Quick Start

```bash
cd /alphafold/combio/software/NGSAbDiscov

/usr/local/lib/NGS_pipeline/annaconda3/bin/python __main__.py validate-sample-sheet \
  --config config.yaml \
  --sample-sheet /path/to/SampleSheet.xlsx \
  --fastq-folder /path/to/Fastq \
  --lib mixed

/usr/local/lib/NGS_pipeline/annaconda3/bin/python __main__.py run-all \
  --config config.yaml \
  --sample-sheet /path/to/SampleSheet.xlsx \
  --fastq-folder /path/to/Fastq \
  --output-folder results_run_name \
  --lib mixed \
  --report-format html \
  --run-ml \
  --ml-backend delphi
```

For most sequencing runs, use `--lib mixed` and let the sample sheet `library`
column decide whether each sample is Fab or VHH.

## Sample Sheet

The minimum required columns are:

```text
TubeBarcode
Sample_Name
library
```

`TubeBarcode` must match the FASTQ filename prefix. FASTQ files should look like:

```text
<TubeBarcode>_S*_L001_R1_001.fastq.gz
<TubeBarcode>_S*_L001_R2_001.fastq.gz
```

`Sample_Name` must have five fields separated by double underscores:

```text
target__block__round__arm__condition
```

Example:

```text
hDKK2_175-259_nMBP_cSt_AH__Block209_Fab__Round4__F_P__20nM
```

The `library` column should usually be:

```text
fab
vhh
```

These names are resolved by `config.yaml` to the configured library definitions
such as `standard_fab` and `vhh_full`.

All samples in the same `target + block` group should use the same library. A
single MiSeq run can contain both Fab and VHH targets, but do not mix Fab and
VHH rows inside the same target/block group.

Legacy sheets with `Sample_ID` or `Description` may still be accepted, but new
sheets should use `TubeBarcode` and `Sample_Name`.

## Validation

Always validate before running:

```bash
cd /alphafold/combio/software/NGSAbDiscov

/usr/local/lib/NGS_pipeline/annaconda3/bin/python __main__.py validate-sample-sheet \
  --config config.yaml \
  --sample-sheet /alphafold/combio/NGS/<RUN_FOLDER>/SampleSheet.xlsx \
  --fastq-folder /alphafold/combio/NGS/<RUN_FOLDER>/Fastq \
  --lib mixed
```

Validation checks:

- required columns
- library names
- `Sample_Name` format
- duplicate `TubeBarcode` or `Sample_Name`
- FASTQ R1/R2 pairing
- target/block groups with mixed libraries

If validation fails, fix a copy of the sample sheet and validate again. Keep the
original sample sheet unchanged when possible.

## Full Pipeline

Run from the software directory:

```bash
cd /alphafold/combio/software/NGSAbDiscov

/usr/local/lib/NGS_pipeline/annaconda3/bin/python __main__.py run-all \
  --config config.yaml \
  --sample-sheet /alphafold/combio/NGS/<RUN_FOLDER>/SampleSheet_validated.xlsx \
  --fastq-folder /alphafold/combio/NGS/<RUN_FOLDER>/Fastq \
  --output-folder results_<RUN_NAME> \
  --lib mixed \
  --report-format html \
  --run-ml \
  --ml-backend delphi
```

The output folder is created under the MiSeq run folder, beside `Fastq`.

Example:

```text
/alphafold/combio/NGS/NGS_20250727_Miseq115/results_full_miseq115_all_20260727_192136
```

## Long Runs With nohup

Use `nohup` for full runs so the job continues after closing the local computer.

```bash
cd /alphafold/combio/software/NGSAbDiscov

RUN_DIR=/alphafold/combio/NGS/NGS_YYYYMMDD_MiseqXXX
STAMP=$(date +%Y%m%d_%H%M%S)
OUT=results_full_miseqXXX_${STAMP}
LOG=${RUN_DIR}/miseqXXX_fullrun_${STAMP}.log
PID_FILE=${RUN_DIR}/miseqXXX_fullrun_${STAMP}.pid

nohup /usr/local/lib/NGS_pipeline/annaconda3/bin/python __main__.py run-all \
  --config config.yaml \
  --sample-sheet ${RUN_DIR}/SampleSheet_validated.xlsx \
  --fastq-folder ${RUN_DIR}/Fastq \
  --output-folder ${OUT} \
  --lib mixed \
  --report-format html \
  --run-ml \
  --ml-backend delphi \
  > ${LOG} 2>&1 &

echo $! > ${PID_FILE}
echo "Output folder: ${RUN_DIR}/${OUT}"
echo "Log: ${LOG}"
echo "PID file: ${PID_FILE}"
```

Check status:

```bash
PID=$(cat /path/to/run.pid)
ps -p "$PID" -o pid,stat,pcpu,pmem,etime,args
tail -f /path/to/run.log
```

## Main Outputs

The result folder contains:

```text
sample_qc_table.csv
data_summary.csv
ml_summary.csv
vh_cdr3_prevalent.csv
*_clones.csv
by_protein/*_final_leads.xlsx
cluster_repeat/*.csv
plots/*.png
report.html
run_manifest.json
```

The pipeline also updates the NGS registry:

```text
/alphafold/combio/NGS/ngs_registry.jsonl
```

## Clone Selection App

The Dash/Plotly clone-selection app reads an existing results folder. It does
not rerun FASTQ processing or ML.

App source code lives in:

```text
dash_plotly_app/clone_selection_app.py
```

Run on `ipinode1`:

```bash
cd /alphafold/combio/software/NGSAbDiscov

/usr/local/lib/NGS_pipeline/annaconda3/bin/python __main__.py clone-app \
  --results-root /alphafold/combio/NGS \
  --host 127.0.0.1 \
  --port 8050 \
  --max-table-rows 500
```

To keep the app running after logout:

```bash
cd /alphafold/combio/software/NGSAbDiscov
RESULTS_ROOT=/alphafold/combio/NGS
LOG=/alphafold/combio/NGS/<RUN_FOLDER>/clone_selection_app.log
PIDFILE=/alphafold/combio/NGS/<RUN_FOLDER>/clone_selection_app.pid

nohup /usr/local/lib/NGS_pipeline/annaconda3/bin/python __main__.py clone-app \
  --results-root "$RESULTS_ROOT" \
  --host 127.0.0.1 \
  --port 8050 \
  --max-table-rows 500 \
  > "$LOG" 2>&1 &

echo $! > "$PIDFILE"
```

From a laptop, open an SSH tunnel in another terminal:

```bash
ssh -L 8050:127.0.0.1:8050 ipinode1
```

Then open:

```text
http://127.0.0.1:8050
```

The app supports:

- MiSeq/result-folder selection from `/alphafold/combio/NGS`
- single-target lead-selection workflow
- automatic refresh when MiSeq, target, or filters change
- `max_freq`, Delphi, PSR, SEC, rank, and CDR3 filters
- default Delphi mean score cutoff of `0.3`
- exclusion of clones present in F_N, negative, or `2uM` control samples
- enrichment/depletion filtering using explicit condition pairs such as `4nM / 20nM`, or manual sample columns
- condition-aware enrichment filtering from available `100nM`, `20nM`, and `4nM` frequency columns
- CDR3 diversity mode with one representative per CDR3 cluster
- compact viewport-fitting layout with tabbed plots
- fold-change scatter plot from condition pairs or selected frequency columns
- computed `log2_fold_change` shown in the table and CSV export
- compact CDR3 distance tree, defaulting to top representative clusters
- selected-clone frequency profile
- sortable/filterable clone table
- CSV download for selected or filtered clones, including full `HSEQ` and `LSEQ`

## Clone Grouping

Fab clone files are grouped by:

```text
CDR3
vh_scaffold
vl_scaffold
anarci_anno
```

VHH clone files are grouped by:

```text
CDR1
CDR2
CDR3
vh_scaffold
```

The combination step groups within:

```text
target
block
library
```

This prevents samples from different blocks or library types from being merged
into the same clone table.

## Machine Learning Defaults

ML uses the Delphi backend when launched with:

```bash
--run-ml --ml-backend delphi
```

Default ML input selection per target:

```text
top 500 clones by max_freq
plus any clone with max_freq >= 0.0005
```

Rows must have valid sequence input before Delphi prediction.

Fab ML requires:

```text
full HSEQ
full LSEQ
no stop codon in sequence input
```

Fab `HSEQ` and `LSEQ` are generated from the scaffold library CSV:

```text
/alphafold/combio/software/NGSAbDiscov/data/IPI_VLVH_LIB_ALL.csv
```

Fab sequence logic:

```text
HSEQ = VH scaffold AA + CDR3 + FR4
LSEQ = VL scaffold AA
```

VHH ML requires:

```text
full HSEQ
no LSEQ
no stop codon in sequence input
```

VHH uses heavy-only Delphi models.

Incomplete or invalid sequence rows remain in the final lead workbook, but they
are not sent to Delphi and will not receive ML scores.

## Expected Final Workbook Checks

For Fab final lead files:

```text
HSEQ should be present for all valid Fab rows
LSEQ should be present for all valid Fab rows
mean_psr_score / mean_sec_score / mean_delphi_score should be present for ML-scored rows
```

For VHH final lead files:

```text
HSEQ should be present for valid VHH rows
LSEQ should be empty
mean_psr_score / mean_sec_score / mean_delphi_score should be present for ML-scored rows
```

## Troubleshooting

If validation reports a library mismatch, check whether a row has `library=vhh`
inside a Fab block, or `library=fab` inside a VHH block.

If validation reports bad FASTQ pairing, check that each `TubeBarcode` has
exactly one matching R1 and one matching R2 file.

If Fab `LSEQ` is missing, check that `vl_scaffold` names match official names in
`IPI_VLVH_LIB_ALL.csv`.

If ML is slow, check the active Delphi process:

```bash
ps -eo pid,ppid,stat,pcpu,pmem,etime,args | grep delphi.py | grep -v grep
```

If a run stops before report generation, check the log first:

```bash
tail -200 /path/to/run.log
```

Then rerun only the needed step if appropriate:

```bash
/usr/local/lib/NGS_pipeline/annaconda3/bin/python __main__.py ml-prediction \
  --config config.yaml \
  --folder /path/to/results_folder \
  --backend delphi

/usr/local/lib/NGS_pipeline/annaconda3/bin/python __main__.py generate-plots \
  --folder /path/to/results_folder \
  --report-format html
```

## Notes

Use only `config.yaml` for current Fab, VHH, and mixed-library runs. The older
`config_VHH.yaml` is obsolete for normal operation.

Do not manually rename scaffold calls in output files. Scaffold names should
come from the configured barcode/library definitions and annotation logic.
