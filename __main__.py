"""Command-line interface for the IPI NGS antibody discovery pipeline."""

from pathlib import Path

import click

from pipeline.config_loader import MIXED_LIBRARY, load_config, resolve_library_name
from pipeline.sample_sheet import SampleSheetValidationError, validate_sample_sheet


def _load_config(config: str, library: str | None = None, output_folder: str | None = None) -> dict:
    cfg = load_config(config)
    if library:
        cfg["current_library"] = resolve_library_name(library, cfg, allow_mixed=True)
        cfg["_library_cli_override"] = cfg["current_library"] != MIXED_LIBRARY
    else:
        cfg["_library_cli_override"] = False
    if output_folder:
        cfg["general"]["output_folder"] = output_folder
    return cfg


def _handle_input_error(exc: Exception):
    if isinstance(exc, (SampleSheetValidationError, ValueError)):
        raise click.ClickException(str(exc)) from exc
    raise exc


@click.group()
def cli():
    """NGSAbDiscov antibody discovery pipeline."""


@cli.command("validate-sample-sheet")
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--sample-sheet", required=True, help="Path to sample sheet Excel")
@click.option("--fastq-folder", required=True, help="Path to FASTQ folder")
@click.option("--library", "--lib", "-l", default=None, help="Override library: fab, vhh, mixed, or a concrete library name")
def validate_sample_sheet_cmd(config, sample_sheet, fastq_folder, library):
    """Validate sample sheet fields and FASTQ R1/R2 pairing."""
    try:
        cfg = _load_config(config, library)
        sheet = validate_sample_sheet(
            Path(sample_sheet),
            Path(fastq_folder),
            allowed_libraries=tuple(cfg["libraries"].keys()),
            current_library=cfg["current_library"],
            library_aliases=cfg.get("library_aliases"),
            force_library=cfg.get("_library_cli_override", False),
            require_fastq=True,
        )
    except Exception as exc:
        _handle_input_error(exc)

    library_counts = sheet["library"].value_counts().to_dict() if "library" in sheet.columns else {}
    click.echo(f"Validation passed: {len(sheet)} samples. Libraries: {library_counts}")


@cli.command()
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--sample-sheet", required=True, help="Path to sample sheet Excel")
@click.option("--fastq-folder", required=True, help="Path to FASTQ folder")
@click.option("--output-folder", default=None, help="Output folder name under the MiSeq run folder")
@click.option("--check-repeats", is_flag=True, help="Reserved for early repeat check compatibility")
@click.option("--library", "--lib", "-l", default=None, help="Override library: fab, vhh, mixed, or a concrete library name")
def process(config, sample_sheet, fastq_folder, output_folder, check_repeats, library):
    """Step 1: Process raw FASTQ files."""
    try:
        cfg = _load_config(config, library, output_folder)
        from pipeline.step1_process import run_processing

        run_processing(cfg, Path(sample_sheet), Path(fastq_folder), cfg["general"]["output_folder"])
    except Exception as exc:
        _handle_input_error(exc)


@cli.command()
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--folder", required=True, help="Results folder containing per-sample files")
@click.option("--library", "--lib", "-l", default=None, help="Override library: fab, vhh, mixed, or a concrete library name")
def combine(config, folder, library):
    """Step 2: Combine per-sample tables per target."""
    try:
        cfg = _load_config(config, library)
        from pipeline.step2_combine import run_combination

        run_combination(cfg, Path(folder))
    except Exception as exc:
        _handle_input_error(exc)


@cli.command("pick-leads")
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--folder", required=True, help="Results folder")
def pick_leads(config, folder):
    """Step 3: Global lead selection."""
    cfg = _load_config(config)
    from pipeline.step3_pick_leads import run_pick_leads

    run_pick_leads(cfg, Path(folder))


@cli.command("check-repeats")
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--folder", required=True, help="Results folder")
def check_repeats(config, folder):
    """Step 4: Repeat check."""
    cfg = _load_config(config)
    from pipeline.step4_repeat_check import run_repeat_check

    run_repeat_check(cfg, Path(folder))


@cli.command("generate-plots")
@click.option("--folder", required=True, help="Results folder")
@click.option("--report-format", default="both", type=click.Choice(["html", "pdf", "both"]))
def generate_plots(folder, report_format):
    """Step 5: Generate QC, analysis plots, and report files."""
    from pipeline.step5_plot_report import make_plots

    make_plots(Path(folder), report_format=report_format)


@cli.command("ml-prediction")
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--folder", required=True, help="Results folder")
@click.option("--backend", default=None, type=click.Choice(["auto", "ipi_psr", "delph", "delphi"]))
def ml_prediction(config, folder, backend):
    """Step 6: Optional ML prediction on clone or lead files."""
    cfg = _load_config(config)
    from pipeline.step6_ml_prediction import run_ml_prediction

    run_ml_prediction(Path(folder), cfg=cfg, backend=backend)


@cli.command("stage-b-export")
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--folder", required=True, help="Results folder containing by_protein/*_final_leads.xlsx")
@click.option("--output", default=None, help="Output CSV path. Default: <folder>/stage_b/abforge_stage_b_candidates.csv")
@click.option("--target-sequences", default=None, help="CSV/TSV/XLSX with target and antigen sequence columns")
@click.option("--top-per-target", default=None, type=int, help="Optional top-N candidates per target for Stage B")
@click.option("--max-total", default=None, type=int, help="Optional total candidate cap after cross-target sorting")
def stage_b_export(config, folder, output, target_sequences, top_per_target, max_total):
    """Export Stage A leads for AbForge/OpenDDE structural scoring."""
    cfg = _load_config(config)
    stage_b_cfg = cfg.get("stage_b", {})
    target_sequence_table = target_sequences or stage_b_cfg.get("target_sequence_table")

    from pipeline.stage_b_abforge import export_stage_b_candidates

    result = export_stage_b_candidates(
        Path(folder),
        output_csv=Path(output) if output else None,
        target_sequence_table=Path(target_sequence_table) if target_sequence_table else None,
        top_per_target=top_per_target,
        max_total=max_total,
    )
    click.echo("Stage B export complete")
    click.echo(f"  candidates: {result.output_csv}")
    click.echo(f"  ready only: {result.ready_csv}")
    click.echo(f"  summary: {result.summary_json}")
    click.echo(f"  rows: {result.rows:,}; ready: {result.ready_rows:,}; targets: {result.targets:,}")
    if result.missing_hseq_rows or result.missing_antigen_rows:
        click.echo(f"  missing HSEQ: {result.missing_hseq_rows:,}; missing antigen sequence: {result.missing_antigen_rows:,}")


@cli.command("stage-b-import")
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--folder", required=True, help="Results folder containing by_protein/*_final_leads.xlsx")
@click.option("--scores", required=True, help="AbForge/OpenDDE scored CSV")
@click.option("--output-dir", default=None, help="Output folder. Default: <folder>/stage_b")
@click.option("--target-sequences", default=None, help="CSV/TSV/XLSX with target and antigen sequence columns")
@click.option("--write-back", is_flag=True, help="Overwrite by_protein workbooks instead of writing stage_b/enriched_by_protein")
def stage_b_import(config, folder, scores, output_dir, target_sequences, write_back):
    """Import AbForge/OpenDDE structural scores back into Stage A leads."""
    cfg = _load_config(config)
    stage_b_cfg = cfg.get("stage_b", {})
    target_sequence_table = target_sequences or stage_b_cfg.get("target_sequence_table")

    from pipeline.stage_b_abforge import import_stage_b_scores

    result = import_stage_b_scores(
        Path(folder),
        Path(scores),
        output_dir=Path(output_dir) if output_dir else None,
        target_sequence_table=Path(target_sequence_table) if target_sequence_table else None,
        write_back=write_back,
    )
    click.echo("Stage B import complete")
    click.echo(f"  combined csv: {result.combined_csv}")
    click.echo(f"  combined xlsx: {result.combined_xlsx}")
    click.echo(f"  summary: {result.summary_json}")
    click.echo(f"  rows: {result.rows:,}; scored: {result.scored_rows:,}; files written: {result.files_written:,}")


@cli.command("clone-app")
@click.option("--results-folder", "--folder", default=None, help="Single results folder containing by_protein/*_final_leads.xlsx")
@click.option("--results-root", default=None, help="Root folder containing multiple MiSeq NGSAbDiscov results")
@click.option("--host", default="127.0.0.1", help="Bind host for the Dash app")
@click.option("--port", default=8050, type=int, help="Bind port for the Dash app")
@click.option("--debug", is_flag=True, help="Run Dash in debug mode")
@click.option("--max-table-rows", default=500, type=int, help="Maximum rows sent to the browser table")
def clone_app(results_folder, results_root, host, port, debug, max_table_rows):
    """Run the interactive clone-selection Dash app."""
    from dash_plotly_app.clone_selection_app import run_app

    run_app(
        Path(results_folder) if results_folder else None,
        results_root=Path(results_root) if results_root else None,
        host=host,
        port=port,
        debug=debug,
        max_table_rows=max_table_rows,
    )


@cli.command("run-all")
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--sample-sheet", required=True, help="Path to sample sheet Excel")
@click.option("--fastq-folder", required=True, help="Path to FASTQ folder")
@click.option("--output-folder", default=None, help="Output folder name under the MiSeq run folder")
@click.option("--check-repeats", is_flag=True, help="Reserved for early repeat check compatibility")
@click.option("--library", "--lib", "-l", default=None, help="Override library: fab, vhh, mixed, or a concrete library name")
@click.option("--report-format", default="both", type=click.Choice(["html", "pdf", "both"]))
@click.option("--run-ml/--skip-ml", default=None, help="Override ml.enabled from config")
@click.option("--ml-backend", default=None, type=click.Choice(["auto", "ipi_psr", "delph", "delphi"]))
def run_all(config, sample_sheet, fastq_folder, output_folder, check_repeats, library, report_format, run_ml, ml_backend):
    """Run the full pipeline end-to-end."""
    try:
        cfg = _load_config(config, library, output_folder)
        if run_ml is not None:
            cfg["ml"]["enabled"] = run_ml

        fastq_folder_path = Path(fastq_folder)
        out = fastq_folder_path.parent / cfg["general"]["output_folder"]

        from pipeline.step1_process import run_processing
        from pipeline.step2_combine import run_combination
        from pipeline.step3_pick_leads import run_pick_leads
        from pipeline.step4_repeat_check import run_repeat_check
        from pipeline.step5_plot_report import make_plots
        from pipeline.step6_ml_prediction import run_ml_prediction
        from pipeline.ngs_registry import update_registry

        click.echo("=== STEP 1: Processing FASTQ ===")
        run_processing(cfg, Path(sample_sheet), fastq_folder_path, cfg["general"]["output_folder"])

        click.echo("=== STEP 2: Combining per target ===")
        run_combination(cfg, out)

        click.echo("=== STEP 3: Picking leads ===")
        run_pick_leads(cfg, out)

        click.echo("=== STEP 4: Repeat check ===")
        run_repeat_check(cfg, out)

        if cfg.get("ml", {}).get("enabled", True):
            click.echo("=== STEP 5: ML prediction ===")
            run_ml_prediction(out, cfg=cfg, backend=ml_backend)
        else:
            click.echo("=== STEP 5: ML prediction skipped by config/CLI ===")

        click.echo("=== STEP 6: Plotting final report ===")
        make_plots(out, report_format=report_format)

        update_registry(
            cfg,
            run_dir=fastq_folder_path.parent,
            results_dir=out,
            sample_sheet=Path(sample_sheet),
            fastq_folder=fastq_folder_path,
        )
    except Exception as exc:
        _handle_input_error(exc)

    click.echo("\nFULL PIPELINE COMPLETE")
    click.echo(f"Results in: {out}")
    click.echo("  sample_qc_table.csv")
    click.echo("  *_clones.csv")
    click.echo("  by_protein/*_final_leads.xlsx")
    click.echo("  data_summary.csv / ml_summary.csv")
    click.echo("  plots/*.png")
    click.echo("  report.html / report.pdf")
    click.echo("  run_manifest.json")


if __name__ == "__main__":
    cli()
