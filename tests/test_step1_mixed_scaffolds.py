import pandas as pd

from pipeline.config_loader import load_config, resolve_library_name
from pipeline.step1_process import (
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


def test_manual_cdr3_never_backfills_anarci_cdr3():
    df = pd.DataFrame(
        {
            "cdr3_aa": ["MANUAL_ONE", "MANUAL_TWO"],
            "CDR3": ["ANARCI_ONE", None],
        }
    )

    result = _preserve_cdr3_sources(df)

    assert result["CDR3"].tolist() == ["ANARCI_ONE", ""]
    assert result["cdr3_aa"].tolist() == ["MANUAL_ONE", "MANUAL_TWO"]
    assert result["cdr3_aa_manualsearch"].tolist() == ["MANUAL_ONE", "MANUAL_TWO"]
    assert result["cdr3_aa_len"].tolist() == [10, 10]
