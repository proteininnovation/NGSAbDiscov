# pipeline/step1_process.py
"""
Step 1: Raw FASTQ processing — faithful recreation of your original fabs3.py
All scientific logic, columns, and stats preserved exactly.
Fully configurable via YAML with library switching.
"""

import os
import sys
import subprocess
import json
import gzip as gz
import pandas as pd
import shutil
from pathlib import Path
from collections import Counter, defaultdict
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from time import asctime
from pprint import pprint
from utilities.liabilities import annotate_liabilities
from utilities.diversity_index import shannon_diversity, evenness,inverse_simpson_index
from anarci import anarci
from abnumber import Chain
import regex as re
from pandarallel import pandarallel
from pipeline.sample_sheet import (
    SampleSheetValidationError,
    validate_sample_sheet as validate_ngs_sample_sheet,
)

pandarallel.initialize(25)

# Global config — set in run_processing
cfg = None

min_aa_len=40
FULL_HEAVY_MIN_AA_LEN = 105
HEAVY_FR4 = "WGQGTLVTVSS"


def _prepend_tool_dirs_to_path() -> None:
    """Make subprocess-backed annotation tools visible in non-interactive runs."""
    candidates = [Path(sys.executable).parent]
    if cfg:
        for key in ("fastp_path", "mafft_path", "hmmscan_path"):
            value = cfg.get("general", {}).get(key)
            if value:
                candidates.append(Path(value).expanduser().parent)

    existing = os.environ.get("PATH", "")
    paths = [str(path) for path in candidates if path and path.exists()]
    os.environ["PATH"] = os.pathsep.join([*paths, existing])


def _clean_sequence_part(value) -> str:
    if pd.isna(value):
        return ""
    value = str(value).strip()
    return "" if value.lower() in {"nan", "none"} else value


def _build_hseq_from_regions(row: pd.Series) -> str:
    pieces = [_clean_sequence_part(row.get(col, "")) for col in ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]]
    if all(pieces):
        assembled = "".join(pieces)
        if len(assembled) >= FULL_HEAVY_MIN_AA_LEN and HEAVY_FR4 in assembled:
            return assembled

    chain = _clean_sequence_part(row.get("CHAIN", ""))
    cdr3 = _clean_sequence_part(row.get("CDR3", ""))
    if (
        chain
        and chain != "Not Fully Annotated"
        and cdr3
        and cdr3 in chain
        and HEAVY_FR4 in chain
        and len(chain) >= FULL_HEAVY_MIN_AA_LEN
    ):
        return chain
    return ""


def _preserve_cdr3_sources(df: pd.DataFrame) -> pd.DataFrame:
    """Use ANARCI CDR3 first and retain motif extraction as an audit/fallback."""
    if "cdr3_aa_manualsearch" in df.columns:
        manual_cdr3 = df["cdr3_aa_manualsearch"].fillna("").astype(str)
    else:
        manual_cdr3 = df.get("cdr3_aa", pd.Series("", index=df.index)).fillna("").astype(str)

    anarci_cdr3 = df.get("CDR3_ANARCI", df.get("CDR3", pd.Series("", index=df.index)))
    anarci_cdr3 = anarci_cdr3.fillna("").astype(str)
    primary_cdr3 = anarci_cdr3.where(anarci_cdr3.ne(""), manual_cdr3)

    df["CDR3_ANARCI"] = anarci_cdr3
    df["cdr3_aa_manualsearch"] = manual_cdr3
    df["CDR3"] = primary_cdr3
    df["cdr3_aa"] = primary_cdr3
    df["cdr3_aa_len"] = primary_cdr3.str.len()
    df["cdr3_source"] = ""
    df.loc[manual_cdr3.ne(""), "cdr3_source"] = "motif_fallback"
    df.loc[anarci_cdr3.ne(""), "cdr3_source"] = "anarci"
    df["cdr3_concordant"] = pd.NA
    both_present = anarci_cdr3.ne("") & manual_cdr3.ne("")
    df.loc[both_present, "cdr3_concordant"] = (
        anarci_cdr3.loc[both_present] == manual_cdr3.loc[both_present]
    )
    return df

def validate_sample_sheet(sample_sheet_path: Path) -> pd.DataFrame:
    if not sample_sheet_path.exists():
        print(f"\nERROR: Sample sheet not found: {sample_sheet_path}\n")
        sys.exit(1)

    try:
        sheet = pd.read_excel(
            sample_sheet_path,
            dtype={"TubeBarcode": str, "Sample_ID": str, "Sample_Name": str, "Description": str, "library": str},
        )
    except Exception as e:
        print(f"\nERROR: Failed to read sample sheet: {e}\n")
        sys.exit(1)

    for col in ["TubeBarcode", "Sample_ID", "Sample_Name", "Description", "library"]:
        if col in sheet.columns:
            sheet[col] = sheet[col].astype("string").fillna("").str.strip()

    if "TubeBarcode" in sheet.columns:
        if "Sample_ID" in sheet.columns:
            empty_tube_barcode = sheet["TubeBarcode"].eq("")
            from_sample_id = empty_tube_barcode & sheet["Sample_ID"].ne("")
            if from_sample_id.any():
                sheet.loc[from_sample_id, "TubeBarcode"] = sheet.loc[from_sample_id, "Sample_ID"]
                print(f"WARNING: Filled {int(from_sample_id.sum())} blank TubeBarcode value(s) from legacy Sample_ID.")
            mismatched = (
                sheet["TubeBarcode"].ne("")
                & sheet["Sample_ID"].ne("")
                & (sheet["TubeBarcode"] != sheet["Sample_ID"])
            )
            if mismatched.any():
                print(
                    f"WARNING: Sample_ID is ignored because TubeBarcode is canonical; "
                    f"{int(mismatched.sum())} row(s) differ."
                )
            sheet = sheet.drop(columns=["Sample_ID"])
    elif "Sample_ID" in sheet.columns:
        sheet["TubeBarcode"] = sheet["Sample_ID"]
        sheet = sheet.drop(columns=["Sample_ID"])
        print("WARNING: Legacy Sample_ID column was used as TubeBarcode; please rename it to TubeBarcode.")
    else:
        print("\nERROR: Missing required 'TubeBarcode' column.\n")
        sys.exit(1)
    id_col = "TubeBarcode"

    if "Sample_Name" in sheet.columns:
        if "Description" in sheet.columns:
            duplicated_id_name = (
                sheet["TubeBarcode"].ne("")
                & sheet["Description"].ne("")
                & (sheet["Sample_Name"] == sheet["TubeBarcode"])
            )
            if duplicated_id_name.any():
                sheet.loc[duplicated_id_name, "Sample_Name"] = sheet.loc[duplicated_id_name, "Description"]
                print(
                    f"WARNING: Replaced {int(duplicated_id_name.sum())} legacy Sample_Name value(s) "
                    "that duplicated TubeBarcode/Sample_ID with Description."
                )

            empty_sample_name = sheet["Sample_Name"].eq("")
            from_description = empty_sample_name & sheet["Description"].ne("")
            if from_description.any():
                sheet.loc[from_description, "Sample_Name"] = sheet.loc[from_description, "Description"]
                print(f"WARNING: Filled {int(from_description.sum())} blank Sample_Name value(s) from legacy Description.")
            mismatched = (
                sheet["Sample_Name"].ne("")
                & sheet["Description"].ne("")
                & (sheet["Sample_Name"] != sheet["Description"])
            )
            if mismatched.any():
                print(
                    f"WARNING: Description is ignored because Sample_Name is canonical; "
                    f"{int(mismatched.sum())} row(s) differ."
                )
            sheet = sheet.drop(columns=["Description"])
    elif "Description" in sheet.columns:
        sheet["Sample_Name"] = sheet["Description"]
        sheet = sheet.drop(columns=["Description"])
        print("WARNING: Legacy Description column was used as Sample_Name; please rename it to Sample_Name.")
    else:
        print("\nERROR: Missing required 'Sample_Name' column.\n")
        sys.exit(1)

    if "library" not in sheet.columns:
        print("\nERROR: Missing required 'library' column.\n")
        sys.exit(1)

    errors = []

    # Check for duplicates and invalid names
    seen_names = set()
    for idx, row in sheet.iterrows():
        sample_id = str(row[id_col]).strip() if pd.notna(row[id_col]) else ""
        sample_name = str(row["Sample_Name"]).strip()
        library = str(row["library"]).strip()

        if not sample_id:
            errors.append(f"Row {idx+2}: Empty TubeBarcode")

        if sample_id in seen_names:
            errors.append(f"Row {idx+2}: Duplicate TubeBarcode '{sample_id}'")
        seen_names.add(sample_id)

        if re.search(r'[<>:"/\\|?*\x00-\x1F]', sample_id):
            errors.append(f"Row {idx+2}: Invalid characters in TubeBarcode '{sample_id}'")

        if sample_name in ["", "nan"]:
            errors.append(f"Row {idx+2} ({sample_id}): Sample_Name is empty")
            continue

        if library in ["", "nan"]:
            errors.append(f"Row {idx+2} ({sample_id}): library is empty")

        parts = sample_name.split("__")
        if len(parts) != 5:
            errors.append(f"Row {idx+2} ({sample_id}): Sample_Name has {len(parts)} parts (need 5)\n   Value: '{sample_name}'")

        if any(not p.strip() for p in parts):
            errors.append(f"Row {idx+2} ({sample_id}): Empty field in Sample_Name: '{sample_name}'")

    if errors:
        print("\n" + "="*70)
        print("SAMPLE SHEET VALIDATION FAILED - SAMPLE NAME / DESCRIPTION ISSUES")
        print("="*70)
        for e in errors:
            print(f"   • {e}")
        print("\nPlease fix the sample sheet and rerun.")
        print("Common fixes:")
        print("   - No empty or duplicate TubeBarcode values")
        print("   - No special characters in TubeBarcode or Sample_Name: <>:?*")
        print("   - Sample_Name must have exactly 5 parts: target__block__round__arm__condition")
        print("   - library must use exact config names such as standard_fab or vhh_full")
        sys.exit(1)

    print(f"Sample sheet validation passed: {len(sheet)} samples (no name/description errors)\n")
    return sheet


def parse_description(sheet: pd.DataFrame) -> pd.DataFrame:
    """
    Parse Sample_Name into antigen, block, round, arm, condition
    Avoids duplicate columns
    """
    split = sheet["Sample_Name"].str.split("__", n=5, expand=True)
    
    # Clean column names
    split.columns = ["antigen", "block", "round", "arm", "condition"]
    
    # Strip whitespace from all parsed columns
    split = split.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    
    # Drop any existing parsed columns to avoid duplicates
    cols_to_drop = ["antigen", "block", "round", "arm", "condition"]
    sheet = sheet.drop(columns=[c for c in cols_to_drop if c in sheet.columns], errors='ignore')
    
    # Concat only once
    return pd.concat([sheet, split], axis=1)

def parse_description_nk(sheet: pd.DataFrame) -> pd.DataFrame:
    split = sheet["Sample_Name"].str.split("__", n=5, expand=True)
    split.columns = ["antigen", "block", "round", "arm", "condition"]
    # Critical: strip whitespace from all parsed fields
    split = split.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    print(sheet)
    return pd.concat([sheet, split], axis=1)

def merge_fastq_pair(name: str, fastq1: Path, fastq2: Path) -> tuple[dict, Path]:
    folder = fastq1.parent
    out_fastq = folder / f"{name}_merged.fastq.gz"
    json_out = folder / f"{name}_merged.json"

    fp = cfg["fastp"]

    cmd = [
        cfg["general"]["fastp_path"],
        "--in1", str(fastq1),
        "--in2", str(fastq2),
        "--merge",
        "--merged_out", str(out_fastq),
        "--json", str(json_out),
        "--qualified_quality_phred", str(fp["qualified_quality_phred"]),
        "--unqualified_percent_limit", str(fp["unqualified_percent_limit"]),
        "--length_required", str(fp["length_required"]),
        "--n_base_limit", str(fp["n_base_limit"]),
        "--overlap_len_require", str(fp["overlap_len_require"]),
        "--overlap_diff_limit", str(fp["overlap_diff_limit"]),
        "--thread", str(fp["thread"]),
    ]

    if fp["correction"]:
        cmd.append("--correction")
    if fp["disable_adapter_trimming"]:
        cmd.append("--disable_adapter_trimming")
    if fp["disable_trim_poly_g"]:
        cmd.append("--disable_trim_poly_g")

    print(cmd)
    print(f"{asctime()} Running fastp merge for {name}...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"fastp failed for {name}")
        print(result.stderr)
        raise RuntimeError("fastp failed")

    print(result.stdout)

    data = json.load(open(json_out))
    read_count = {
        "total": data["read1_before_filtering"]["total_reads"],
        "merged": data["merged_and_filtered"]["total_reads"],
        "low_quality": int(data["filtering_result"]["low_quality_reads"] / 2),
        "too_many_N": int(data["filtering_result"]["too_many_N_reads"] / 2),
        "too_short": int(data["filtering_result"]["too_short_reads"] / 2),
        "too_long": int(data["filtering_result"]["too_long_reads"] / 2),
    }

    pprint(read_count)
    return read_count, out_fastq

def load_fastq_to_df(fastq_path: Path) -> pd.DataFrame:
    seqs = Counter()
    with gz.open(fastq_path, "rt") as fin:
        for title, sequence, quality in FastqGeneralIterator(fin):
            seqs[str(Seq(sequence).reverse_complement())] += 1

    print(asctime(), f"# reads: {sum(seqs.values())}")
    print(asctime(), f"# unique DNA: {len(seqs)}")

    df = pd.DataFrame({"nt": list(seqs.keys()), "count": list(seqs.values())})
    del seqs
    return df

def split_at_delimiters(seq: str, delimiters: list[str], split_downstream: bool = True) -> int:
    for delim in delimiters:
        idx = seq.find(delim)
        if idx != -1:
            return idx + len(delim) if split_downstream else idx
    return -1

def _translate_delimiters(delimiters: list[str]) -> tuple[str, ...]:
    translated = []
    for delimiter in delimiters:
        delimiter = str(delimiter).strip().upper()
        if delimiter and len(delimiter) % 3 == 0:
            peptide = str(Seq(delimiter).translate())
            if "*" not in peptide and peptide not in translated:
                translated.append(peptide)
    return tuple(translated)


def _oriented_nt(row: pd.Series) -> str:
    nt = _clean_sequence_part(row.get("nt", "")).upper()
    return nt if row.get("aa_strand", "forward") == "forward" else str(Seq(nt).reverse_complement())


def _nt_for_aa_interval(row: pd.Series, aa_start: int, aa_end: int) -> tuple[str, int, int]:
    frame = int(row.get("aa_frame", 0) or 0)
    start = frame + 3 * aa_start
    end = frame + 3 * aa_end
    nucleotide = _oriented_nt(row)[start:end]
    expected = _clean_sequence_part(row.get("aa", ""))[aa_start:aa_end]
    if not nucleotide or str(Seq(nucleotide).translate()) != expected:
        return "", -1, -1
    return nucleotide, start, end


def _extract_anarci_cdr3(row: pd.Series) -> tuple[str, int, int]:
    cdr3 = _clean_sequence_part(row.get("CDR3", ""))
    aa = _clean_sequence_part(row.get("aa", ""))
    if not cdr3 or not aa:
        return "", -1, -1

    fr3 = _clean_sequence_part(row.get("FR3", ""))
    fr4 = _clean_sequence_part(row.get("FR4", ""))
    context = f"{fr3}{cdr3}{fr4}"
    context_start = aa.find(context) if fr3 and fr4 else -1
    aa_start = context_start + len(fr3) if context_start >= 0 else aa.find(cdr3)
    if aa_start < 0:
        return "", -1, -1
    return _nt_for_aa_interval(row, aa_start, aa_start + len(cdr3))


def _extract_motif_cdr3(
    row: pd.Series,
    upstream_peptides: tuple[str, ...],
    downstream_peptides: tuple[str, ...],
) -> tuple[str, str, int, int]:
    aa = _clean_sequence_part(row.get("aa", ""))
    candidates: list[tuple[int, int, int]] = []
    for upstream in upstream_peptides:
        upstream_start = aa.find(upstream)
        while upstream_start >= 0:
            cdr3_start = upstream_start + len(upstream)
            for downstream in downstream_peptides:
                cdr3_end = aa.find(downstream, cdr3_start)
                if cdr3_end >= cdr3_start:
                    candidates.append((cdr3_end - cdr3_start, cdr3_start, cdr3_end))
            upstream_start = aa.find(upstream, upstream_start + 1)
    if not candidates:
        return "", "", -1, -1

    _, aa_start, aa_end = min(candidates)
    peptide = aa[aa_start:aa_end]
    nucleotide, nt_start, nt_end = _nt_for_aa_interval(row, aa_start, aa_end)
    if not nucleotide:
        return "", "", -1, -1
    return peptide, nucleotide, nt_start, nt_end


def find_hcdr3(df: pd.DataFrame, upseq: list[str], downseq: list[str], read_count: dict) -> tuple[pd.DataFrame, dict]:
    """Resolve CDR3 from ANARCI first, with amino-acid motif extraction as fallback."""
    df = df.copy()
    upstream_peptides = _translate_delimiters(upseq)
    downstream_peptides = _translate_delimiters(downseq)

    motif_results = df.apply(
        lambda row: _extract_motif_cdr3(row, upstream_peptides, downstream_peptides),
        axis=1,
        result_type="expand",
    )
    motif_results.columns = [
        "cdr3_aa_manualsearch",
        "cdr3_nt_manualsearch",
        "cdr3_beg_manualsearch",
        "cdr3_end_manualsearch",
    ]
    for column in motif_results.columns:
        df[column] = motif_results[column]

    anarci_results = df.apply(_extract_anarci_cdr3, axis=1, result_type="expand")
    anarci_results.columns = ["cdr3_nt_anarci", "cdr3_beg_anarci", "cdr3_end_anarci"]
    for column in anarci_results.columns:
        df[column] = anarci_results[column]

    df = _preserve_cdr3_sources(df)
    anarci_mask = df["CDR3_ANARCI"].ne("")
    motif_mask = df["cdr3_aa_manualsearch"].ne("")
    accepted_mask = df["cdr3_aa"].ne("")

    anarci_nt_valid = df["cdr3_nt_anarci"].fillna("").astype(str).ne("")
    manual_matches_primary = (
        df["cdr3_aa_manualsearch"].fillna("").astype(str) == df["cdr3_aa"]
    )
    use_manual_nt = (~anarci_nt_valid) & manual_matches_primary & motif_mask
    df["cdr3_nt"] = df["cdr3_nt_anarci"].fillna("").astype(str)
    df.loc[use_manual_nt, "cdr3_nt"] = df.loc[use_manual_nt, "cdr3_nt_manualsearch"]
    df["cdr3_nt_source"] = ""
    df.loc[use_manual_nt, "cdr3_nt_source"] = "motif"
    df.loc[anarci_nt_valid, "cdr3_nt_source"] = "anarci_coordinates"

    df["cdr3_beg"] = df["cdr3_beg_anarci"]
    df["cdr3_end"] = df["cdr3_end_anarci"]
    df.loc[use_manual_nt, "cdr3_beg"] = df.loc[use_manual_nt, "cdr3_beg_manualsearch"]
    df.loc[use_manual_nt, "cdr3_end"] = df.loc[use_manual_nt, "cdr3_end_manualsearch"]
    df["cdr3_mod3"] = df["cdr3_nt"].str.len() % 3
    df["cdr3_functional"] = (~df["cdr3_aa"].str.contains("\\*", na=False)) & (
        df["cdr3_nt"].eq("") | df["cdr3_mod3"].eq(0)
    )

    total_reads = int(df["count"].sum())
    total_sequences = len(df)
    motif_reads = int(df.loc[motif_mask, "count"].sum())
    anarci_reads = int(df.loc[anarci_mask, "count"].sum())
    fallback_mask = (~anarci_mask) & motif_mask
    accepted_reads = int(df.loc[accepted_mask, "count"].sum())
    concordant_mask = df["cdr3_concordant"].map(
        lambda value: bool(value) if not pd.isna(value) else False
    )
    discordant_mask = anarci_mask & motif_mask & (~concordant_mask)

    read_count["reads_no_cdr3_edges"] = total_reads - motif_reads
    read_count["seqs_no_cdr3_edges"] = total_sequences - int(motif_mask.sum())
    read_count["reads_no_cdr3_annotation"] = total_reads - accepted_reads
    read_count["seqs_no_cdr3_annotation"] = total_sequences - int(accepted_mask.sum())
    read_count["reads_cdr3_from_anarci"] = anarci_reads
    read_count["seqs_cdr3_from_anarci"] = int(anarci_mask.sum())
    read_count["reads_cdr3_from_motif_fallback"] = int(df.loc[fallback_mask, "count"].sum())
    read_count["seqs_cdr3_from_motif_fallback"] = int(fallback_mask.sum())
    read_count["reads_cdr3_discordant"] = int(df.loc[discordant_mask, "count"].sum())
    read_count["seqs_cdr3_discordant"] = int(discordant_mask.sum())

    df = df.loc[accepted_mask].copy()
    df.sort_values("count", ascending=False, inplace=True)
    return df, read_count

def get_label_from_barcode(seq: str, barcodes: dict, errors_allowed: int = 1) -> str:
    if errors_allowed > 0:
        for label, barcode in barcodes.items():
            if re.search(f"({barcode}){{e<={errors_allowed}}}", seq):
                return label
    for label, barcode in barcodes.items():
        if barcode in seq:
            return label
    return "UNK"

def _scaffold_settings_for_library(config: dict, library_name: str) -> dict:
    """Return an isolated scaffold configuration for one concrete library."""
    libraries = config.get("libraries", {})
    if library_name not in libraries:
        raise ValueError(f"Unknown concrete library '{library_name}'")

    library_cfg = libraries[library_name]
    library_type = str(library_cfg.get("library_type", "")).strip().lower()
    if library_type not in {"fab", "vhh"}:
        raise ValueError(
            f"Library '{library_name}' must declare library_type as 'fab' or 'vhh', "
            f"not '{library_type or '<empty>'}'"
        )

    vh_barcodes = dict(library_cfg.get("vh_barcodes") or {})
    if not vh_barcodes:
        raise ValueError(f"Library '{library_name}' has no VH scaffold set")

    settings = {
        "library_type": library_type,
        "vh_barcodes": vh_barcodes,
        "vh_region": tuple(library_cfg.get("vh_barcode_region", (0, -1))),
        "vl_barcodes": None,
        "vl_region": None,
    }
    if library_type == "fab":
        vl_barcodes = dict(library_cfg.get("vl_barcodes") or {})
        if not vl_barcodes:
            raise ValueError(f"Fab library '{library_name}' has no VL scaffold set")
        settings["vl_barcodes"] = vl_barcodes
        settings["vl_region"] = tuple(library_cfg.get("vl_barcode_region", (0, -1)))

    return settings


def find_vh_vl(
    df: pd.DataFrame,
    vh_barcodes: dict,
    vl_barcodes: dict = None,
    *,
    vh_barcode_region: tuple[int, int],
    vl_barcode_region: tuple[int, int] | None = None,
):
    df["vh_scaffold"] = df["nt"].apply(
        lambda x: get_label_from_barcode(
            x[vh_barcode_region[0]:vh_barcode_region[1]], vh_barcodes, 1
        )
    )
    if vl_barcodes:
        if vl_barcode_region is None:
            raise ValueError("A VL barcode region is required when a VL scaffold set is supplied")
        df["vl_scaffold"] = df["nt"].apply(
            lambda x: get_label_from_barcode(
                x[vl_barcode_region[0]:vl_barcode_region[1]], vl_barcodes, 1
            )
        )
    elif "vl_scaffold" in df.columns:
        df.drop(columns=["vl_scaffold"], inplace=True)
    return df

def get_full_anarci_anno(df: pd.DataFrame) -> pd.DataFrame:
    _prepend_tool_dirs_to_path()
    if shutil.which("hmmscan") is None:
        raise RuntimeError(
            "hmmscan was not found on PATH; ANARCI/AbNumber cannot annotate CDR1/CDR2/CDR3. "
            "Install hmmer or set general.hmmscan_path in the config."
        )

    def get_ANARCI(seq: str):
        if len(seq) <= min_aa_len:
            return "Not Fully Annotated|||||||"
        try:
            chain = Chain(seq, scheme="imgt")
            return f"{chain}|{chain.fr1_seq}|{chain.cdr1_seq}|{chain.fr2_seq}|{chain.cdr2_seq}|{chain.fr3_seq}|{chain.cdr3_seq}|{chain.fr4_seq}"
        except Exception:
            return "Not Fully Annotated|||||||"

    df["anarci"] = df["aa"].parallel_apply(get_ANARCI)
    df[["CHAIN", "FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3","FR4"]] = df["anarci"].str.split("|", expand=True)
    for col in ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]:
        df[col] = df[col].fillna("").astype(str)
    return df

def _translate_best_orf(seq: str) -> tuple[str, str, int]:
    candidates: list[tuple[str, str, int]] = []
    forward = str(seq)
    reverse = str(Seq(forward).reverse_complement())
    for strand, oriented in (("forward", forward), ("reverse", reverse)):
        for frame in range(3):
            coding_sequence = oriented[frame:]
            coding_sequence = coding_sequence[: len(coding_sequence) - (len(coding_sequence) % 3)]
            candidates.append((str(Seq(coding_sequence).translate(to_stop=True)), strand, frame))
    return max(candidates, key=lambda item: len(item[0])) if candidates else ("", "forward", 0)


def NT2AA(seq: str) -> str:
    return _translate_best_orf(seq)[0]

def consolidate(df: pd.DataFrame) -> pd.DataFrame:
    drop_cols = [
        "nt", "cdr3_beg", "cdr3_end", "cdr3_beg_anarci", "cdr3_end_anarci",
        "cdr3_beg_manualsearch", "cdr3_end_manualsearch", "cdr3_mod3", "anarci",
        "aa_strand", "aa_frame",
    ]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns])
    #df["cdr3_aa"]=df["CDR3"]
    #df["cdr3_aa_len"] = df["CDR3"].str.len()
    # FOR  ALFA TAG: group_cols = ["cdr3_aa", "anarci_anno", "vh_scaffold", "vl_scaffold", "CDR1", "CDR2", "CDR3","FR1","FR2","FR3","FR4"]
    group_cols = ["cdr3_aa", "anarci_anno", "vh_scaffold", "vl_scaffold", "CDR1", "CDR2", "CDR3"]
    group_cols = [col for col in group_cols if col in df.columns]
    agg_cols = {"count": "sum", "cdr3_aa_len": "first"}
    for col in [
        "FR1", "FR2", "FR3", "FR4", "HSEQ", "CHAIN", "aa",
        "CDR3_ANARCI", "cdr3_aa_manualsearch", "cdr3_nt", "cdr3_nt_anarci",
        "cdr3_nt_manualsearch", "cdr3_source", "cdr3_nt_source",
        "cdr3_concordant", "cdr3_functional",
    ]:
        if col in df.columns and col not in group_cols:
            agg_cols[col] = "first"
    df = df.groupby(group_cols).agg(agg_cols).reset_index()


    #df = (
    #        df.sort_values("count", ascending=False)
    #        .groupby(group_cols, as_index=False)
    #        .head(1)                                      # keep the row of the chain with highest count
    #        .assign(count = lambda x:
    #            df.groupby(group_cols)["count"].transform("sum")
    #        )
    #        # Keep only the grouping columns + CHAIN + total_count
    #        .loc[:, group_cols + ["CHAIN", "count", "cdr3_aa_len"]]
    #        .reset_index(drop=True)
    #        )

    df.sort_values("count", ascending=False, inplace=True)
    df["rank"] = range(1, len(df) + 1)
    df["freq"] = df["count"] / df["count"].sum()
    return df

def remove_crap(df: pd.DataFrame, read_count: dict, min_h3_len=1, max_h3_len=30,
                functional_only=False, keep_vh_unk=True, keep_vl_unk=True,
                min_freq=0.0, min_count=1) -> tuple[pd.DataFrame, dict]:
    n_reads = df["count"].sum()
    n_seqs = len(df)

    df = df[(df["cdr3_aa_len"] >= min_h3_len) & (df["cdr3_aa_len"] <= max_h3_len)]
    read_count["reads_with_short_cdr3"] = n_reads - df["count"].sum()
    read_count["seqs_with_short_cdr3"] = n_seqs - len(df)

    n_reads = df["count"].sum()
    n_seqs = len(df)

    if functional_only:
        df = df[df["cdr3_functional"]]
        read_count["reads_with_nonfunctional_cdr3"] = n_reads - df["count"].sum()
        read_count["seqs_with_nonfunctional_cdr3"] = n_seqs - len(df)

    n_reads = df["count"].sum()
    n_seqs = len(df)

    df = df[(df["freq"] >= min_freq) & (df["count"] >= min_count)]
    read_count["reads_with_low_frequency"] = n_reads - df["count"].sum()
    read_count["seqs_with_low_frequency"] = n_seqs - len(df)

    n_reads = df["count"].sum()
    n_seqs = len(df)

    if not keep_vh_unk:
        df = df[df["vh_scaffold"] != "UNK"]
    if not keep_vl_unk and "vl_scaffold" in df.columns:
        df = df[df["vl_scaffold"] != "UNK"]

    df = df[~df["cdr3_aa"].str.contains("X", na=False)]
    read_count["reads_ambiguous"] = n_reads - df["count"].sum()
    read_count["seqs_ambiguous"] = n_seqs - len(df)

    df.sort_values("count", ascending=False, inplace=True)
    df["rank"] = range(1, len(df) + 1)
    return df, read_count

def run_processing(cfg_in, sample_sheet: Path, fastq_folder: Path, output_folder="results"):
    global cfg
    cfg = cfg_in

    lib = cfg["current_library"]

    fastq_folder = Path(fastq_folder)
    configured_output = output_folder or cfg["general"].get("output_folder", "results")
    cfg["general"]["output_folder"] = configured_output
    output_dir = fastq_folder.parent / configured_output
    output_dir.mkdir(parents=True, exist_ok=True)

    sheet = validate_ngs_sample_sheet(
        sample_sheet,
        fastq_folder,
        allowed_libraries=tuple(cfg["libraries"].keys()),
        current_library=lib,
        library_aliases=cfg.get("library_aliases"),
        force_library=cfg.get("_library_cli_override", False),
        require_fastq=True,
        allow_missing_fastq=cfg["processing"].get("allow_missing_fastq", False),
    )
    warnings = sheet.attrs.get("validation_warnings", [])
    if warnings:
        print("\nSample sheet validation warnings:")
        for warning in warnings:
            print(f"  - {warning}")
        sheet = sheet.dropna(subset=["read1", "read2"])

    library_counts = sheet["library"].value_counts().to_dict() if "library" in sheet.columns else {lib: len(sheet)}
    print(f"Sample sheet validation passed: {len(sheet)} samples. Libraries: {library_counts}")

    sample_qc_table = pd.DataFrame(columns=[
        'name', 'library', 'library_type', 'antigen', 'block', 'round', 'arm', 'condition',
        'total', 'merged', 'low_quality', 'too_many_N', 'too_short', 'too_long',
        'reads_no_cdr3_edges', 'seqs_no_cdr3_edges',
        'reads_no_cdr3_annotation', 'seqs_no_cdr3_annotation',
        'reads_cdr3_from_anarci', 'seqs_cdr3_from_anarci',
        'reads_cdr3_from_motif_fallback', 'seqs_cdr3_from_motif_fallback',
        'reads_cdr3_discordant', 'seqs_cdr3_discordant',
        'reads_with_short_cdr3', 'seqs_with_short_cdr3',
        'reads_with_low_frequency', 'seqs_with_low_frequency',
        'reads_ambiguous', 'seqs_ambiguous',
        'unique_dna', 'unique_dna_pct', 'merged_pct',
        'anarci_anno','anarci_anno_unique_dna',
        'unique_fr1','unique_fr2','unique_fr3',
        'unique_cdr1','unique_cdr2',
        'unique_cdr3', 'total_cdr3', 'shannon', 'evenness',
        'unique_cdr3_min2', 'total_cdr3_min2', 'shannon_min2', 'evenness_min2',
        'VH_UNK', 'VL_UNK',
        'n_cross_target_clones', 'n_cross_target_reads', 'pct_cross_target_reads'
    ])

    per_sample_files = []

    for _, row in sheet.iterrows():
        name = row["Sample_Name"]
        fastq1 = row["read1"]
        fastq2 = row["read2"]
        row_library = row.get("library", lib)
        scaffold_settings = _scaffold_settings_for_library(cfg, row_library)
        library_type = scaffold_settings["library_type"]

        print(f"\n{asctime()} === Processing {name} [{row_library}] ===")

        read_count, merged_fastq = merge_fastq_pair(name, fastq1, fastq2)

        df = load_fastq_to_df(merged_fastq)

        translations = df["nt"].apply(_translate_best_orf)
        df["aa"] = translations.str[0]
        df["aa_strand"] = translations.str[1]
        df["aa_frame"] = translations.str[2]
        df["aa_len"] = df["aa"].str.len()

        df = get_full_anarci_anno(df)

        df['anarci_anno'] = (df['CHAIN'] != "Not Fully Annotated")
        unique_dna = len(df)

        df, read_count = find_hcdr3(df, cfg["processing"]["yyc_nt"], cfg["processing"]["wgq_nt"], read_count)
        df = find_vh_vl(
            df,
            scaffold_settings["vh_barcodes"],
            scaffold_settings["vl_barcodes"],
            vh_barcode_region=scaffold_settings["vh_region"],
            vl_barcode_region=scaffold_settings["vl_region"],
        )

        df["HSEQ"] = df.apply(_build_hseq_from_regions, axis=1)



        # === RAW DF WITH METADATA ===
        df_before_consolidate = df.copy()

        #Add metadata using pd.Series (broadcasts scalar safely, no index alignment)
        metadata = {
            "library": row_library,
            "library_type": library_type,
            "antigen": row["antigen"],
            "block": row["block"],
            "round": row["round"],
            "arm": row["arm"],
            "condition": row["condition"],
        }
        for col, value in metadata.items():
            df_before_consolidate[col] = pd.Series([value] * len(df_before_consolidate), index=df_before_consolidate.index)

        # === FINAL DF ===
        df = consolidate(df)
        # Add metadata to final df using same safe method
        for col, value in metadata.items():
            df[col] = pd.Series([value] * len(df), index=df.index)




        df, read_count = remove_crap(
            df, read_count,
            min_h3_len=cfg["processing"]["hcdr3_min_len"],
            max_h3_len=cfg["processing"]["hcdr3_max_len"],
            functional_only=False,
            keep_vh_unk=True,
            keep_vl_unk=True,
            min_freq=cfg["processing"]["read_freq_min"],
            min_count=cfg["processing"]["read_count_min"]
        )


        df = annotate_liabilities(df)

        # QC stats
        read_count.update({
            "name": name,
            "library": row_library,
            "library_type": library_type,
            "antigen": row["antigen"],
            "block": row["block"],
            "round": row["round"],
            "arm": row["arm"],
            "condition": row["condition"],
            "unique_dna": unique_dna,
            "unique_dna_pct": round(100 * unique_dna / read_count["total"], 2) if read_count["total"] > 0 else 0,
            "merged_pct": round(100 * read_count["merged"] / read_count["total"], 2) if read_count["total"] > 0 else 0,
            "anarci_anno":round(100*df_before_consolidate.loc[df_before_consolidate.anarci_anno]['count'].sum()/df_before_consolidate['count'].sum(),2) if df_before_consolidate['count'].sum() else 0.0,
            "anarci_anno_unique_dna":round(100*len(df_before_consolidate.loc[df_before_consolidate.anarci_anno])/len(df_before_consolidate),2) if len(df_before_consolidate) else 0.0,
            "unique_fr1":len(df_before_consolidate['FR1'].value_counts()),
            "unique_fr2":len(df_before_consolidate['FR2'].value_counts()),
            "unique_fr3":len(df_before_consolidate['FR3'].value_counts()),
            "unique_cdr1":len(df_before_consolidate['CDR1'].value_counts()),
            "unique_cdr2":len(df_before_consolidate['CDR2'].value_counts()),
            "unique_cdr3": len(df_before_consolidate['CDR3'].value_counts()),
            "total_cdr3": df_before_consolidate["count"].sum(),
            "shannon": shannon_diversity(df["count"]),
            "inv_simpsons": inverse_simpson_index(df["count"]),
            "evenness": evenness(shannon_diversity(df["count"]), len(df)),
            "unique_cdr3_min2": len(df_before_consolidate[df_before_consolidate["count"] > 1]['CDR3'].value_counts()),
            "total_cdr3_min2": df_before_consolidate[df_before_consolidate["count"] > 1]["count"].sum(),
            "shannon_min2": shannon_diversity(df[df["count"] > 1]["count"]),
            "evenness_min2": evenness(shannon_diversity(df[df["count"] > 1]["count"]), len(df[df["count"] > 1])),
            "VH_UNK": (df_before_consolidate["vh_scaffold"] == "UNK").sum(),
            "VL_UNK": (df_before_consolidate["vl_scaffold"] == "UNK").sum() if "vl_scaffold" in df_before_consolidate.columns else 0,
            "n_cross_target_clones": 0,
            "n_cross_target_reads": 0,
            "pct_cross_target_reads": 0.0
        })

        sample_qc_table = pd.concat([sample_qc_table, pd.DataFrame([read_count])], ignore_index=True)




        # Save
        out_file = output_dir / f"{name}_{row['block']}.csv.gz"
        os.makedirs(f"{output_dir}/raw", exist_ok=True)
        out_raw = f"{output_dir}/raw/{name}_{row['block']}_raw.csv.gz"

        out_file.parent.mkdir(parents=True, exist_ok=True)
        print(out_file)
        df.to_csv(out_file, index=False, compression="gzip")

        total_raw=df_before_consolidate["count"].sum()
        if total_raw > 0:
            df_before_consolidate = df_before_consolidate.drop(
                columns=[
                    'anarci', 'cdr3_beg', 'cdr3_end', 'cdr3_beg_anarci', 'cdr3_end_anarci',
                    'cdr3_beg_manualsearch', 'cdr3_end_manualsearch', 'cdr3_mod3',
                    'antigen', 'block', 'round', 'arm', 'condition', 'aa', 'aa_strand', 'aa_frame'
                ],
                errors='ignore',
            )
            df_before_consolidate["freq"] =  df_before_consolidate["count"] / total_raw
            df_before_consolidate.to_csv(out_raw, index=False, compression="gzip")

        per_sample_files.append(out_file)


    # === NEW: Cross-target VH-CDR3 contamination (updated logic) ===
    print("\nCalculating cross-target VH-CDR3 contamination (excluding top-frequency target)...")
    
    data = []
    for f in per_sample_files:
        df = pd.read_csv(f, compression="gzip")
        if df.empty or "antigen" not in df.columns or "count" not in df.columns:
            continue
        sample_stem = f.name.replace(".csv.gz", "")
        antigen = df["antigen"].iloc[0]
        total_reads = df["count"].sum()
        if total_reads == 0:
            continue
        cdr3_col = "cdr3_aa" if "cdr3_aa" in df.columns else "CDR3"
        if cdr3_col not in df.columns or "vh_scaffold" not in df.columns:
            continue
        sample_data = df[["vh_scaffold", cdr3_col, "count"]].copy()
        if cdr3_col != "cdr3_aa":
            sample_data = sample_data.rename(columns={cdr3_col: "cdr3_aa"})
        sample_data = (
            sample_data
            .groupby(["vh_scaffold", "cdr3_aa"], dropna=False, as_index=False)["count"]
            .sum()
        )
        sample_data["sample"] = sample_stem
        sample_data["antigen"] = antigen
        sample_data["freq"] = sample_data["count"] / total_reads
        data.append(sample_data)

    pair_to_prevalent = {}
    pair_to_top_antigens = {}

    if data:
        contam_df = pd.concat(data, ignore_index=True)

        pair_cols = ["vh_scaffold", "cdr3_aa"]
        antigen_max = (
            contam_df
            .groupby([*pair_cols, "antigen"], dropna=False)["freq"]
            .max()
            .rename("antigen_max_freq")
            .reset_index()
        )
        pair_summary = (
            antigen_max
            .groupby(pair_cols, dropna=False)
            .agg(
                prevalent_targets=("antigen", "nunique"),
                global_max_freq=("antigen_max_freq", "max"),
                targets_list=("antigen", lambda s: ";".join(sorted(map(str, pd.unique(s)))))
            )
            .reset_index()
        )

        contam_df = contam_df.merge(antigen_max, on=[*pair_cols, "antigen"], how="left")
        contam_df = contam_df.merge(pair_summary[pair_cols + ["prevalent_targets", "global_max_freq"]], on=pair_cols, how="left")
        contam_df["contaminant_vh_cdr3"] = (
            (contam_df["prevalent_targets"] > 1)
            & (contam_df["antigen_max_freq"] < contam_df["global_max_freq"])
        )

        top_antigen_rows = antigen_max.merge(pair_summary[pair_cols + ["global_max_freq"]], on=pair_cols, how="left")
        top_antigen_rows = top_antigen_rows[top_antigen_rows["antigen_max_freq"] == top_antigen_rows["global_max_freq"]]
        top_antigens = (
            top_antigen_rows
            .groupby(pair_cols, dropna=False)["antigen"]
            .agg(lambda s: ";".join(sorted(map(str, pd.unique(s)))))
            .rename("top_frequency_antigen")
            .reset_index()
        )

        top_sample_rows = contam_df[contam_df["freq"] == contam_df["global_max_freq"]]
        top_samples = (
            top_sample_rows
            .groupby(pair_cols, dropna=False)
            .agg(
                top_frequency_sample=("sample", lambda s: ";".join(sorted(map(str, pd.unique(s))))),
                top_frequency_sample_reads=("count", "sum"),
            )
            .reset_index()
        )

        top_total_reads = (
            contam_df[contam_df["antigen_max_freq"] == contam_df["global_max_freq"]]
            .groupby(pair_cols, dropna=False)["count"]
            .sum()
            .rename("top_frequency_total_reads")
            .reset_index()
        )

        pair_info = (
            pair_summary
            .merge(top_antigens, on=pair_cols, how="left")
            .merge(top_samples, on=pair_cols, how="left")
            .merge(top_total_reads, on=pair_cols, how="left")
        )
        pair_info["top_frequency_antigen"] = pair_info["top_frequency_antigen"].fillna("N/A")
        pair_info["top_frequency_sample"] = pair_info["top_frequency_sample"].fillna("N/A")
        pair_info["top_frequency_sample_reads"] = pair_info["top_frequency_sample_reads"].fillna(0).astype(int)
        pair_info["top_frequency_total_reads"] = pair_info["top_frequency_total_reads"].fillna(0).astype(int)

        pair_to_prevalent = {
            (row.vh_scaffold, row.cdr3_aa): int(row.prevalent_targets)
            for row in pair_info.itertuples(index=False)
        }
        pair_to_top_antigens = {
            (row.vh_scaffold, row.cdr3_aa): set(str(row.top_frequency_antigen).split(";"))
            for row in pair_info.itertuples(index=False)
        }

        # Per-sample contamination stats (using per-occurrence contaminant flag)
        contam_group = contam_df[contam_df["contaminant_vh_cdr3"]]
        if not contam_group.empty:
            contam_stats = contam_group.groupby("sample").apply(
                lambda g: pd.Series({
                    "n_cross_target_clones": g[["vh_scaffold", "cdr3_aa"]].drop_duplicates().shape[0],
                    "n_cross_target_reads": g["count"].sum()
                })
            ).reset_index()
    
            sample_totals = contam_df.groupby("sample")["count"].sum().reset_index(name="total_reads")
    
            contam_stats = contam_stats.merge(sample_totals, on="sample", how="left")
            contam_stats["total_reads"] = contam_stats["total_reads"].fillna(1)
            contam_stats["pct_cross_target_reads"] = round(100 * contam_stats["n_cross_target_reads"] / contam_stats["total_reads"], 2)
    
            for _, stat in contam_stats.iterrows():
                sample_name = stat["sample"].rsplit("_", 1)[0]
                mask = sample_qc_table["name"] == sample_name
                if mask.any():
                    sample_qc_table.loc[mask, ["n_cross_target_clones", "n_cross_target_reads", "pct_cross_target_reads"]] = [
                        stat["n_cross_target_clones"],
                        stat["n_cross_target_reads"],
                        stat["pct_cross_target_reads"]
                    ]
    
        empty_cols = [
            "vh_cdr3", "vh_scaffold", "cdr3_aa", "prevalent_targets",
            "top_frequency_antigen", "top_frequency_sample",
            "top_frequency_sample_reads", "top_frequency_total_reads",
            "global_max_freq", "targets_list"
        ]
        cross_df = pair_info[pair_info["prevalent_targets"] > 1].copy()
        cross_df["vh_cdr3"] = cross_df["vh_scaffold"].astype(str) + "-" + cross_df["cdr3_aa"].astype(str)
        cross_df = cross_df[empty_cols]
        if not cross_df.empty:
            cross_df = cross_df.sort_values("prevalent_targets", ascending=False)
        cross_df = cross_df[empty_cols]
    
        cross_df.to_csv(output_dir / "vh_cdr3_prevalent.csv", index=False)
        print(f"Saved global cross-target summary: {output_dir / 'vh_cdr3_prevalent.csv'} ({len(cross_df)} cross-target VH-CDR3 pairs)")
    else:
        print("No valid data for contamination calculation — assuming no cross-target contamination.")
        empty_cols = ["vh_cdr3", "vh_scaffold", "cdr3_aa", "prevalent_targets",
                      "top_frequency_antigen", "top_frequency_sample",
                      "top_frequency_sample_reads", "top_frequency_total_reads",
                      "global_max_freq", "targets_list"]
        pd.DataFrame(columns=empty_cols).to_csv(output_dir / "vh_cdr3_prevalent.csv", index=False)
    
    # === Update per-sample files (NO MERGE - avoids duplication bug) ===
    print("Updating per-sample files with prevalence and contamination flags...")
    for f in per_sample_files:
        df = pd.read_csv(f, compression="gzip")
    
        if df.empty:
            # Safe defaults if empty
            df["prevalent_targets_vh_cdr3"] = 1
            df["contaminant_vh_cdr3"] = False
            df.to_csv(f, index=False, compression="gzip")
            continue
    
        antigen = df["antigen"].iloc[0] if "antigen" in df.columns and not df["antigen"].empty else "UNKNOWN"
    
        # Create pair keys
        pairs = list(zip(df["vh_scaffold"], df["cdr3_aa"]))
    
        # Add prevalent_targets_vh_cdr3 (global)
        df["prevalent_targets_vh_cdr3"] = [pair_to_prevalent.get(p, 1) for p in pairs]
    
        # Add contaminant_vh_cdr3 (sample-specific: True only if cross-target AND this antigen is NOT top)
        df["contaminant_vh_cdr3"] = [
            (prev > 1) and (antigen not in pair_to_top_antigens.get(p, set()))
            for p, prev in zip(pairs, df["prevalent_targets_vh_cdr3"])
        ]
    
        df["prevalent_targets_vh_cdr3"] = df["prevalent_targets_vh_cdr3"].astype(int)
        df["contaminant_vh_cdr3"] = df["contaminant_vh_cdr3"].astype(bool)
    
        df.to_csv(f, index=False, compression="gzip")
    
    sample_qc_table.to_csv(output_dir / "sample_qc_table.csv", index=False)
    print("\nStep 1 complete!")
