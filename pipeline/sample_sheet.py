from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from pipeline.config_loader import MIXED_LIBRARY


DESCRIPTION_FIELDS = ["antigen", "block", "round", "arm", "condition"]
INVALID_FILE_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1F]')
FAB_TOKEN = re.compile(r"(^|[_\-\s])fab($|[_\-\s])", re.IGNORECASE)
VHH_TOKEN = re.compile(r"(^|[_\-\s])vhh($|[_\-\s])", re.IGNORECASE)


class SampleSheetValidationError(ValueError):
    """Raised when a sample sheet or its FASTQ pairs are not usable."""


def _format_errors(title: str, errors: list[str]) -> str:
    lines = [title, "=" * len(title)]
    lines.extend(f"- {error}" for error in errors)
    return "\n".join(lines)


def _normalize_sample_sheet(sheet: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    sheet = sheet.copy()
    warnings: list[str] = []

    for col in ["TubeBarcode", "Sample_ID", "Sample_Name", "Description", "library"]:
        if col in sheet.columns:
            sheet[col] = sheet[col].astype("string").fillna("").str.strip()

    if "TubeBarcode" in sheet.columns:
        if "Sample_ID" in sheet.columns:
            empty_tube_barcode = sheet["TubeBarcode"].eq("")
            from_sample_id = empty_tube_barcode & sheet["Sample_ID"].ne("")
            if from_sample_id.any():
                sheet.loc[from_sample_id, "TubeBarcode"] = sheet.loc[from_sample_id, "Sample_ID"]
                warnings.append(
                    f"Filled {int(from_sample_id.sum())} blank TubeBarcode value(s) from legacy Sample_ID."
                )
            mismatched = (
                sheet["TubeBarcode"].ne("")
                & sheet["Sample_ID"].ne("")
                & (sheet["TubeBarcode"] != sheet["Sample_ID"])
            )
            if mismatched.any():
                warnings.append(
                    f"Sample_ID is ignored because TubeBarcode is canonical; "
                    f"{int(mismatched.sum())} row(s) have different TubeBarcode/Sample_ID values."
                )
            sheet = sheet.drop(columns=["Sample_ID"])
    elif "Sample_ID" in sheet.columns:
        sheet["TubeBarcode"] = sheet["Sample_ID"]
        sheet = sheet.drop(columns=["Sample_ID"])
        warnings.append("Legacy Sample_ID column was used as TubeBarcode; please rename this column to TubeBarcode.")
    else:
        raise SampleSheetValidationError(
            "Sample sheet must include 'TubeBarcode'. TubeBarcode is the FASTQ sample identifier."
        )

    if "Sample_Name" in sheet.columns:
        if "Description" in sheet.columns:
            duplicated_id_name = (
                sheet["TubeBarcode"].ne("")
                & sheet["Description"].ne("")
                & (sheet["Sample_Name"] == sheet["TubeBarcode"])
            )
            if duplicated_id_name.any():
                sheet.loc[duplicated_id_name, "Sample_Name"] = sheet.loc[duplicated_id_name, "Description"]
                warnings.append(
                    f"Replaced {int(duplicated_id_name.sum())} legacy Sample_Name value(s) that duplicated "
                    "TubeBarcode/Sample_ID with Description."
                )

            empty_sample_name = sheet["Sample_Name"].eq("")
            from_description = empty_sample_name & sheet["Description"].ne("")
            if from_description.any():
                sheet.loc[from_description, "Sample_Name"] = sheet.loc[from_description, "Description"]
                warnings.append(
                    f"Filled {int(from_description.sum())} blank Sample_Name value(s) from legacy Description."
                )

            mismatched = (
                sheet["Sample_Name"].ne("")
                & sheet["Description"].ne("")
                & (sheet["Sample_Name"] != sheet["Description"])
            )
            if mismatched.any():
                warnings.append(
                    f"Description is ignored because Sample_Name is canonical; "
                    f"{int(mismatched.sum())} row(s) have different Sample_Name/Description values."
                )
            sheet = sheet.drop(columns=["Description"])
    elif "Description" in sheet.columns:
        sheet["Sample_Name"] = sheet["Description"]
        sheet = sheet.drop(columns=["Description"])
        warnings.append("Legacy Description column was used as Sample_Name; please rename this column to Sample_Name.")
    else:
        raise SampleSheetValidationError(
            "Sample sheet must include 'Sample_Name'. Sample_Name replaces the old Description column and must use "
            "target__block__round__arm__condition format."
        )

    return sheet, warnings


def parse_sample_name(sheet: pd.DataFrame) -> pd.DataFrame:
    split = sheet["Sample_Name"].str.split("__", n=4, expand=True)
    split.columns = DESCRIPTION_FIELDS
    split = split.apply(lambda col: col.astype("string").fillna("").str.strip())
    sheet = sheet.drop(columns=[c for c in DESCRIPTION_FIELDS if c in sheet.columns], errors="ignore")
    return pd.concat([sheet, split], axis=1)


def parse_description(sheet: pd.DataFrame) -> pd.DataFrame:
    """Backward-compatible alias for old callers."""
    return parse_sample_name(sheet)


def _library_kind(library: object) -> str:
    value = "" if pd.isna(library) else str(library).strip().lower()
    if "vhh" in value:
        return "vhh"
    if "fab" in value:
        return "fab"
    return ""


def _find_fastq_pair(tube_barcode: str, fastq_files: list[Path], read_num: int) -> tuple[Path | None, list[Path]]:
    token_matches = [
        path
        for path in fastq_files
        if path.name.split(".")[0].split("_")[0] == tube_barcode and f"_R{read_num}_" in path.name
    ]
    if token_matches:
        return (token_matches[0] if len(token_matches) == 1 else None), token_matches

    prefix_matches = [
        path
        for path in fastq_files
        if path.name.startswith(tube_barcode) and f"_R{read_num}_" in path.name
    ]
    return (prefix_matches[0] if len(prefix_matches) == 1 else None), prefix_matches


def validate_sample_sheet(
    sample_sheet_path: Path,
    fastq_folder: Path | None = None,
    *,
    allowed_libraries: list[str] | tuple[str, ...] | None = None,
    current_library: str | None = None,
    library_aliases: dict[str, str] | None = None,
    force_library: bool = False,
    require_fastq: bool = False,
    allow_missing_fastq: bool = False,
) -> pd.DataFrame:
    sample_sheet_path = Path(sample_sheet_path)
    if not sample_sheet_path.exists():
        raise SampleSheetValidationError(f"Sample sheet not found: {sample_sheet_path}")

    try:
        sheet = pd.read_excel(
            sample_sheet_path,
            dtype={"TubeBarcode": str, "Sample_ID": str, "Sample_Name": str, "Description": str, "library": str},
        )
    except Exception as exc:
        raise SampleSheetValidationError(f"Failed to read sample sheet: {exc}") from exc

    sheet, warnings = _normalize_sample_sheet(sheet)

    errors: list[str] = []
    aliases = {str(k).strip().lower(): str(v).strip() for k, v in (library_aliases or {}).items()}

    if current_library and allowed_libraries and current_library not in allowed_libraries:
        if current_library != MIXED_LIBRARY:
            errors.append(
                f"Unknown library '{current_library}'. Known libraries: {', '.join(sorted(allowed_libraries))}"
            )

    def resolve_row_library(raw_value: object, excel_row: int, tube_barcode: str) -> str:
        value = "" if pd.isna(raw_value) else str(raw_value).strip()
        if not value:
            errors.append(
                f"Row {excel_row} ({tube_barcode}): empty library; expected one of "
                f"fab, vhh, or {', '.join(sorted(allowed_libraries or []))}."
            )
            return ""
        resolved = aliases.get(value.lower(), value)
        if resolved == MIXED_LIBRARY:
            errors.append(
                f"Row {excel_row} ({tube_barcode}): sample library cannot be 'mixed'; "
                "use fab, vhh, or a concrete config library name."
            )
            return ""
        if allowed_libraries and resolved not in allowed_libraries:
            known_aliases = sorted(alias for alias, target in aliases.items() if target in allowed_libraries)
            known_values = ", ".join(dict.fromkeys([*known_aliases, *(sorted(allowed_libraries))]))
            errors.append(
                f"Row {excel_row} ({tube_barcode}): unknown library '{value}'. "
                f"Use one of: {known_values}."
            )
            return ""
        return resolved

    if not errors:
        has_library_column = "library" in sheet.columns
        if force_library and current_library and current_library != MIXED_LIBRARY:
            sheet["library"] = current_library
            if has_library_column:
                warnings.append(f"Sample sheet 'library' column ignored because CLI forced library '{current_library}'.")
            else:
                warnings.append(f"Sample sheet has no 'library' column; CLI forced library '{current_library}'.")
        elif has_library_column:
            sheet["library"] = sheet["library"].astype("string").fillna("").str.strip()
        else:
            errors.append(
                "Missing required 'library' column. Use fab, vhh, or exact config names such as standard_fab or vhh_full."
            )

    if not errors:
        seen_ids: set[str] = set()
        seen_sample_names: set[str] = set()
        for idx, row in sheet.iterrows():
            excel_row = idx + 2
            tube_barcode = str(row["TubeBarcode"]).strip()
            sample_name = str(row["Sample_Name"]).strip()
            if "library" in sheet.columns:
                sheet.at[idx, "library"] = resolve_row_library(row.get("library", ""), excel_row, tube_barcode)

            if not tube_barcode:
                errors.append(f"Row {excel_row}: empty TubeBarcode.")
            elif tube_barcode in seen_ids:
                errors.append(f"Row {excel_row}: duplicate TubeBarcode '{tube_barcode}'.")
            elif INVALID_FILE_CHARS.search(tube_barcode):
                errors.append(f"Row {excel_row}: TubeBarcode has filesystem-unsafe characters: '{tube_barcode}'.")
            seen_ids.add(tube_barcode)

            if not sample_name or sample_name.lower() == "nan":
                errors.append(f"Row {excel_row} ({tube_barcode}): Sample_Name is empty.")
                continue
            if sample_name in seen_sample_names:
                errors.append(
                    f"Row {excel_row} ({tube_barcode}): duplicate Sample_Name '{sample_name}' would overwrite output files."
                )
            if INVALID_FILE_CHARS.search(sample_name):
                errors.append(
                    f"Row {excel_row} ({tube_barcode}): Sample_Name has filesystem-unsafe characters: '{sample_name}'."
                )
            seen_sample_names.add(sample_name)

            parts = sample_name.split("__")
            if len(parts) != len(DESCRIPTION_FIELDS):
                errors.append(
                    f"Row {excel_row} ({tube_barcode}): Sample_Name has {len(parts)} parts; expected "
                    "target__block__round__arm__condition."
                )
            elif any(not part.strip() for part in parts):
                errors.append(f"Row {excel_row} ({tube_barcode}): Sample_Name contains an empty '__' field.")

            library_kind = _library_kind(row.get("library", ""))
            if library_kind == "vhh" and FAB_TOKEN.search(sample_name):
                errors.append(
                    f"Row {excel_row} ({tube_barcode}): Sample_Name looks Fab-specific but library is "
                    f"'{row.get('library', '')}'."
                )
            if library_kind == "fab" and VHH_TOKEN.search(sample_name):
                errors.append(
                    f"Row {excel_row} ({tube_barcode}): Sample_Name looks VHH-specific but library is "
                    f"'{row.get('library', '')}'."
                )

    if "library" in sheet.columns:
        valid_name_mask = sheet["Sample_Name"].astype(str).str.split("__").str.len() == len(DESCRIPTION_FIELDS)
        parsed_for_checks = parse_sample_name(sheet.loc[valid_name_mask].copy())
        group_cols = ["antigen", "block"]
        library_counts = (
            parsed_for_checks.groupby(group_cols, dropna=False)["library"]
            .nunique()
            .reset_index(name="n_libraries")
        )
        mixed_groups = library_counts[library_counts["n_libraries"] > 1]
        for group in mixed_groups.itertuples(index=False):
            mask = (
                (parsed_for_checks["antigen"] == group.antigen)
                & (parsed_for_checks["block"] == group.block)
            )
            examples = (
                parsed_for_checks.loc[mask, ["TubeBarcode", "library"]]
                .drop_duplicates()
                .head(8)
                .apply(lambda row: f"{row['TubeBarcode']}:{row['library']}", axis=1)
                .tolist()
            )
            errors.append(
                f"Target/block '{group.antigen} / {group.block}' contains multiple libraries: "
                + ", ".join(examples)
            )

    fastq_files: list[Path] = []
    if fastq_folder is not None:
        fastq_folder = Path(fastq_folder)
        if not fastq_folder.exists():
            errors.append(f"FASTQ folder not found: {fastq_folder}")
        else:
            fastq_files = sorted(fastq_folder.glob("*.fastq.gz"))
            if require_fastq and not fastq_files:
                errors.append(f"No *.fastq.gz files found in {fastq_folder}")

    if fastq_folder is not None and fastq_files and not errors:
        read1_paths: list[Path | None] = []
        read2_paths: list[Path | None] = []
        for idx, row in sheet.iterrows():
            excel_row = idx + 2
            tube_barcode = str(row["TubeBarcode"]).strip()

            read1, read1_matches = _find_fastq_pair(tube_barcode, fastq_files, 1)
            read2, read2_matches = _find_fastq_pair(tube_barcode, fastq_files, 2)

            if len(read1_matches) > 1:
                errors.append(
                    f"Row {excel_row} ({tube_barcode}): multiple R1 FASTQ files matched: "
                    + ", ".join(path.name for path in read1_matches)
                )
            if len(read2_matches) > 1:
                errors.append(
                    f"Row {excel_row} ({tube_barcode}): multiple R2 FASTQ files matched: "
                    + ", ".join(path.name for path in read2_matches)
                )
            if read1 is None or read2 is None:
                message = f"Row {excel_row} ({tube_barcode}): missing FASTQ pair"
                if allow_missing_fastq:
                    warnings.append(message)
                else:
                    errors.append(message)

            read1_paths.append(read1)
            read2_paths.append(read2)

        sheet["read1"] = read1_paths
        sheet["read2"] = read2_paths

    if errors:
        raise SampleSheetValidationError(
            _format_errors(
                "SAMPLE SHEET VALIDATION FAILED",
                errors
                + [
                    "Sample_Name format must be: target__block__round__arm__condition.",
                    "FASTQ names should begin with TubeBarcode and contain _R1_ / _R2_.",
                ],
            )
        )

    parsed = parse_sample_name(sheet)
    parsed.attrs["validation_warnings"] = warnings
    return parsed
