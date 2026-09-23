# LLNL_Project

Subproject workspace for LLNL antibody discovery NGS runs.

## Layout

- `Fastq/`: paired FASTQ inputs for the current run.
- `sample_sheets/`: validated Excel sample sheets.
- `results/`: pipeline outputs created beside `Fastq/`.
- `docs/`: project notes and supporting documents.
- `config.yaml`: LLNL-specific pipeline config, copied from the root config with `registry.project` set to `LLNL_Project` and `current_library` set to `mixed`.

## Run

From the repository root:

```bash
python __main__.py validate-sample-sheet \
  --config LLNL_Project/config.yaml \
  --sample-sheet LLNL_Project/sample_sheets/SampleSheet_validated.xlsx \
  --fastq-folder LLNL_Project/Fastq \
  --lib mixed
```

```bash
python __main__.py run-all \
  --config LLNL_Project/config.yaml \
  --sample-sheet LLNL_Project/sample_sheets/SampleSheet_validated.xlsx \
  --fastq-folder LLNL_Project/Fastq \
  --output-folder results \
  --lib mixed \
  --report-format html \
  --run-ml \
  --ml-backend delphi
```

The sample sheet should include `TubeBarcode`, `Sample_Name`, and `library`.
Use `fab` or `vhh` in the `library` column.
