import pandas as pd

from pipeline.config_loader import load_config, resolve_library_name
from pipeline.step1_process import (
    _translate_best_orf,
    find_hcdr3,
    _preserve_cdr3_sources,
    _scaffold_settings_for_library,
    find_vh_vl,
)


def test_mixed_library_scaffold_sets_are_isolated_per_sample():
    config = {
        "libraries": {
            "standard_fab": {
                "library_type": "fab",
                "vh_barcodes": {"FAB_VH": "AAAACCCCGG"},
                "vl_barcodes": {"FAB_VL": "TTTTGGGGCC"},
                "vh_barcode_region": [0, 40],
                "vl_barcode_region": [0, 40],
            },
            "fab4": {
                "library_type": "fab",
                "vh_barcodes": {"FAB4_VH": "GGGGAAAACC"},
                "vl_barcodes": {"FAB4_VL": "AAAATTTTCC"},
                "vh_barcode_region": [0, 40],
                "vl_barcode_region": [0, 40],
            },
            "vhh_full": {
                "library_type": "vhh",
                "vh_barcodes": {"VHH": "CCCCAAAATT"},
                # Even if a VHH config contains legacy VL entries, they must be ignored.
                "vl_barcodes": {"WRONG_VHH_VL": "GGGGTTTTAA"},
                "vh_barcode_region": [0, 40],
                "vl_barcode_region": [0, 40],
            },
        }
    }

    fab_settings = _scaffold_settings_for_library(config, "standard_fab")
    fab_df = find_vh_vl(
        pd.DataFrame({"nt": ["CCCCAAAATTAAAACCCCGGTTTTGGGGCC"]}),
        fab_settings["vh_barcodes"],
        fab_settings["vl_barcodes"],
        vh_barcode_region=fab_settings["vh_region"],
        vl_barcode_region=fab_settings["vl_region"],
    )
    assert fab_df.loc[0, "vh_scaffold"] == "FAB_VH"
    assert fab_df.loc[0, "vl_scaffold"] == "FAB_VL"

    fab4_settings = _scaffold_settings_for_library(config, "fab4")
    fab4_df = find_vh_vl(
        pd.DataFrame({"nt": ["AAAACCCCGGGGGGAAAACCAAAATTTTCC"]}),
        fab4_settings["vh_barcodes"],
        fab4_settings["vl_barcodes"],
        vh_barcode_region=fab4_settings["vh_region"],
        vl_barcode_region=fab4_settings["vl_region"],
    )
    assert fab4_df.loc[0, "vh_scaffold"] == "FAB4_VH"
    assert fab4_df.loc[0, "vl_scaffold"] == "FAB4_VL"

    vhh_settings = _scaffold_settings_for_library(config, "vhh_full")
    assert vhh_settings["vl_barcodes"] is None
    vhh_df = find_vh_vl(
        pd.DataFrame(
            {
                "nt": ["AAAACCCCGGTTTTGGGGCCCCCAAAATT"],
                "vl_scaffold": ["STALE_FAB_VALUE"],
            }
        ),
        vhh_settings["vh_barcodes"],
        vhh_settings["vl_barcodes"],
        vh_barcode_region=vhh_settings["vh_region"],
        vl_barcode_region=vhh_settings["vl_region"],
    )
    assert vhh_df.loc[0, "vh_scaffold"] == "VHH"
    assert "vl_scaffold" not in vhh_df.columns


def test_fab_aliases_use_standard_fab_and_fab4_remains_separate():
    config = load_config("config.yaml")

    assert resolve_library_name("fab", config) == "standard_fab"
    assert resolve_library_name("fab_standard", config) == "standard_fab"
    assert resolve_library_name("Fab4", config) == "fab4"
    assert config["libraries"]["fab4"]["library_type"] == "fab"


def test_anarci_cdr3_is_primary_and_manual_is_retained_for_audit():
    df = pd.DataFrame(
        {
            "cdr3_aa_manualsearch": ["MANUAL_ONE", "MANUAL_TWO"],
            "CDR3": ["ANARCI_ONE", None],
        }
    )

    result = _preserve_cdr3_sources(df)

    assert result["CDR3_ANARCI"].tolist() == ["ANARCI_ONE", ""]
    assert result["CDR3"].tolist() == ["ANARCI_ONE", "MANUAL_TWO"]
    assert result["cdr3_aa"].tolist() == ["ANARCI_ONE", "MANUAL_TWO"]
    assert result["cdr3_aa_manualsearch"].tolist() == ["MANUAL_ONE", "MANUAL_TWO"]
    assert result["cdr3_source"].tolist() == ["anarci", "motif_fallback"]


def _cdr3_test_frame(anarci_cdr3: str, *, with_motifs: bool = True) -> pd.DataFrame:
    if with_motifs:
        nt = "TACTACTGCGCTCGTTTTTGGGGACAG"
        aa = "YYCARFWGQ"
        fr3 = "YYC"
        fr4 = "WGQ"
    else:
        nt = "CAACAACAATGTGCTCGTTTTTGGGCTGCT"
        aa = "QQQCARFWAA"
        fr3 = "QQQC"
        fr4 = "WAA"
    return pd.DataFrame(
        {
            "nt": [nt],
            "count": [10],
            "aa": [aa],
            "aa_strand": ["forward"],
            "aa_frame": [0],
            "FR3": [fr3],
            "CDR3": [anarci_cdr3],
            "FR4": [fr4],
        }
    )


def test_anarci_coordinates_recover_cdr3_nt_and_synonymous_wgq_crosscheck():
    result, counts = find_hcdr3(
        _cdr3_test_frame("ARF"),
        ["TACTACTGC", "TATTACTGC"],
        ["TGGGGACAA", "TGGGGACAG"],
        {},
    )

    assert result.loc[0, "cdr3_aa"] == "ARF"
    assert result.loc[0, "cdr3_nt"] == "GCTCGTTTT"
    assert result.loc[0, "cdr3_nt_source"] == "anarci_coordinates"
    assert result.loc[0, "cdr3_aa_manualsearch"] == "ARF"
    assert result.loc[0, "cdr3_concordant"] == True
    assert counts["reads_cdr3_from_anarci"] == 10
    assert counts["reads_no_cdr3_annotation"] == 0


def test_anarci_valid_read_is_kept_without_motif_edges():
    result, counts = find_hcdr3(
        _cdr3_test_frame("ARF", with_motifs=False),
        ["TACTACTGC", "TATTACTGC"],
        ["TGGGGACAA", "TGGGGACAG"],
        {},
    )

    assert len(result) == 1
    assert result.loc[0, "cdr3_source"] == "anarci"
    assert result.loc[0, "cdr3_aa_manualsearch"] == ""
    assert counts["reads_no_cdr3_edges"] == 10
    assert counts["reads_no_cdr3_annotation"] == 0


def test_motif_annotation_is_used_only_when_anarci_cdr3_is_missing():
    result, counts = find_hcdr3(
        _cdr3_test_frame(""),
        ["TACTACTGC", "TATTACTGC"],
        ["TGGGGACAA", "TGGGGACAG"],
        {},
    )

    assert result.loc[0, "cdr3_aa"] == "ARF"
    assert result.loc[0, "CDR3"] == "ARF"
    assert result.loc[0, "cdr3_source"] == "motif_fallback"
    assert result.loc[0, "cdr3_nt"] == "GCTCGTTTT"
    assert counts["reads_cdr3_from_motif_fallback"] == 10


def test_best_orf_returns_translation_strand_and_frame():
    peptide, strand, frame = _translate_best_orf("ATGGCC")

    assert peptide == "MA"
    assert strand == "forward"
    assert frame == 0
