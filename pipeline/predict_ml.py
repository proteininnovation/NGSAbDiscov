"""Machine-learning prediction helpers for NGSAbDiscov.

The Delphi backend is label-agnostic: configured commands may target
``psr_filter``, ``sec_filter``, ``spr_filter``, or any other Delphi label
available in the installed model/database. The older IPI in-house backend is
still PSR-specific and remains exposed as ``predict_psr_on_clone``.
"""

from __future__ import annotations

import importlib
import math
import re
import sys
import shlex
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL_DIR = REPO_ROOT / "psr_model_ml"
DEFAULT_ANTIBERTA2_MODEL = "alchemab/antiberta2-cssp"
FR4 = "WGQGTLVTVSS"
FAB_LIGHT_CONSTANT = (
    "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
)
FULL_HEAVY_MIN_AA_LEN = 105
SEQUENCE_BIOPHYSICAL_COLUMNS = [
    "cdr3_arg_count",
    "cdr3_asp_count",
    "cdr3_trp_count",
    "cdr3_glu_count",
    "cdr3_length",
    "cdr3_net_charge",
    "cdr3_isoelectric_point",
    "hseq_isoelectric_point",
]
DEFAULT_FINAL_EXCLUDE_COLUMNS = [
    "prevalent_targets_vh_cdr3",
    "contaminant_vh_cdr3",
    "cross_target_reads",
    "cross_target_samples",
    "library",
    "library_type",
    "CHAIN",
    "aa",
    "FR1",
    "FR2",
    "FR3",
    "FR4",
    "ML_SEQUENCE_OK",
    "HSEQ_SOURCE",
    "LSEQ_SOURCE",
]


def _read_table(input_file: str | Path) -> pd.DataFrame:
    input_file = str(input_file)
    if input_file.lower().endswith(".csv"):
        return pd.read_csv(input_file)
    if input_file.lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(input_file)
    raise ValueError("Input must be .csv, .xlsx, or .xls")


def _write_table(df: pd.DataFrame, output_file: str | Path) -> None:
    output_file = str(output_file)
    if output_file.lower().endswith(".csv"):
        df.to_csv(output_file, index=False)
    else:
        df.to_excel(output_file, index=False)


def _default_output_path(input_file: str | Path, suffix: str) -> Path:
    path = Path(input_file)
    return path.with_name(f"{path.stem}_{suffix}.xlsx")


def _resolve_model_dir(model_dir: str | Path | None = None) -> Path:
    candidates = []
    if model_dir:
        candidates.append(Path(model_dir))
    candidates.extend(
        [
            DEFAULT_MODEL_DIR,
            Path("/alphafold/combio/software/NGSAbDiscov/psr_model_ml"),
            Path("/alphafold/combio/software/IPIAbDiscov/psr_model_ml"),
        ]
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return Path(model_dir) if model_dir else DEFAULT_MODEL_DIR


def _load_ipi_library(model_dir: Path, library_csv: str | Path | None = None) -> pd.DataFrame:
    candidates = []
    if library_csv:
        candidates.append(Path(library_csv))
    candidates.extend([model_dir / "IPI_VLVH_LIB_ALL.csv", REPO_ROOT / "data/IPI_VLVH_LIB_ALL.csv"])
    for candidate in candidates:
        if candidate.exists():
            return pd.read_csv(candidate)
    return pd.DataFrame(columns=["Name", "AA"])


def _clean_sequence(value: Any) -> str:
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value.lower() in {"nan", "none"}:
        return ""
    return re.sub(r"[^A-Za-z]", "", value).upper()


def _has_stop_codon(value: Any) -> bool:
    if pd.isna(value):
        return False
    return "*" in str(value)


def _clean_cdr3_sequence(value: Any) -> str:
    seq = _clean_sequence(value)
    return seq[1:] if seq.startswith("CAR") else seq


def _nonempty_lower_values(df: pd.DataFrame, column: str) -> set[str]:
    if column not in df.columns:
        return set()
    values = set()
    for value in df[column].dropna():
        text = str(value).strip().lower()
        if text and text not in {"nan", "none", "null", "n/a", "na"}:
            values.add(text)
    return values


def _infer_library_type_from_df(df: pd.DataFrame) -> str | None:
    library_types = _nonempty_lower_values(df, "library_type")
    if len(library_types) == 1:
        return next(iter(library_types))
    if len(library_types) > 1:
        return "mixed"

    libraries = _nonempty_lower_values(df, "library")
    if libraries and all("vhh" in value for value in libraries):
        return "vhh"
    if libraries and all("fab" in value for value in libraries):
        return "fab"
    return None


def _column_has_values(df: pd.DataFrame, column: str) -> bool:
    if column not in df.columns:
        return False
    values = df[column].fillna("").astype(str).str.strip().str.lower()
    return bool((values != "").loc[~values.isin({"nan", "none", "null", "n/a", "na"})].any())


def _matching_clone_table(input_file: Path) -> Path | None:
    if input_file.parent.name != "by_protein":
        return None
    stem = input_file.stem.replace("_final_leads", "")
    clone_path = input_file.parent.parent / f"{stem}_clones.csv"
    return clone_path if clone_path.exists() else None


def _clone_merge_keys(df: pd.DataFrame, clone_df: pd.DataFrame) -> list[str]:
    candidates = [
        ["BARCODE"],
        ["CDR1", "CDR2", "CDR3", "vh_scaffold"],
        ["CDR3", "vh_scaffold", "vl_scaffold", "anarci_anno"],
        ["CDR3", "vh_scaffold", "vl_scaffold"],
        ["cdr3_aa", "vh_scaffold", "vl_scaffold", "anarci_anno"],
        ["cdr3_aa", "vh_scaffold", "vl_scaffold"],
        ["CDR3", "vh_scaffold"],
        ["cdr3_aa", "vh_scaffold"],
    ]
    for keys in candidates:
        if all(key in df.columns and key in clone_df.columns for key in keys):
            return keys
    return []


def _enrich_from_clone_table(input_file: Path, df: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    clone_path = _matching_clone_table(input_file)
    if clone_path is None:
        return df, False
    try:
        clone_df = pd.read_csv(clone_path)
    except Exception as exc:
        print(f"  Warning: could not read clone table for ML enrichment: {clone_path.name}: {exc}")
        return df, False

    keys = _clone_merge_keys(df, clone_df)
    if not keys:
        return df, False

    helper_cols = [
        "library",
        "library_type",
        "HSEQ",
        "LSEQ",
        "CHAIN",
        "aa",
        "FR1",
        "FR2",
        "FR3",
        "FR4",
        "CDR1",
        "CDR2",
        "CDR3",
        "cdr3_aa",
        "vh_scaffold",
        "vl_scaffold",
        "anarci_anno",
    ]
    cols_to_add = [
        col
        for col in helper_cols
        if col not in keys
        and col in clone_df.columns
        and (col not in df.columns or not _column_has_values(df, col))
    ]
    if not cols_to_add:
        return df, False

    left = df.copy()
    right = clone_df[keys + cols_to_add].drop_duplicates(keys).copy()
    for key in keys:
        left[key] = left[key].fillna("").astype(str)
        right[key] = right[key].fillna("").astype(str)

    merged = left.merge(right, on=keys, how="left", suffixes=("", "__clone"))
    for col in cols_to_add:
        clone_col = f"{col}__clone"
        if clone_col not in merged.columns:
            continue
        if col not in left.columns:
            merged[col] = merged[clone_col]
        else:
            missing = ~merged[col].apply(lambda value: bool(_clean_sequence(value)))
            merged.loc[missing, col] = merged.loc[missing, clone_col]
        merged.drop(columns=[clone_col], inplace=True)

    print(f"  Enriched Delphi input from {clone_path.name} using keys: {', '.join(keys)}")
    return merged, True


def _merge_library_type_delph_config(delph_config: dict, library_type: str | None) -> dict:
    effective = dict(delph_config or {})
    if not library_type:
        return effective

    command_by_type = effective.get("command_by_library_type") or effective.get("commands_by_library_type") or {}
    if library_type in command_by_type:
        effective["command"] = command_by_type[library_type]

    overrides = effective.get("by_library_type") or effective.get("library_type_overrides") or {}
    if library_type in overrides:
        override = overrides[library_type] or {}
        effective.update(override)

    if library_type == "vhh":
        effective.setdefault("force_no_lseq", True)
        effective.setdefault("require_lseq", False)
    return effective


def _normalize_command_templates(command_template: str | list[str] | tuple[str, ...] | None) -> list[str]:
    if command_template is None:
        return []
    if isinstance(command_template, (list, tuple)):
        return [str(command) for command in command_template if str(command).strip()]
    return [str(command_template)]


def _path_is_within(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def _extract_delph_generated_paths(stdout: str, parent: Path) -> list[Path]:
    paths: list[Path] = []
    markers = ["Saved predictions to:", "[log] Writing to:", "->", "→"]
    for line in stdout.splitlines():
        for marker in markers:
            if marker not in line:
                continue
            candidate = line.split(marker, 1)[1].strip()
            candidate = candidate.split(" | ", 1)[0].strip()
            candidate = candidate.rstrip(".,;")
            if not candidate.startswith("/"):
                continue
            path = Path(candidate)
            if _path_is_within(path, parent):
                paths.append(path)
            break
    return list(dict.fromkeys(paths))


def _library_sequence_lookup(ipi_library: pd.DataFrame) -> dict[str, str]:
    if "Name" not in ipi_library.columns:
        return {}
    seq_col = None
    for candidate in ["AA", "Amino acid seq", "Amino acid sequence"]:
        if candidate in ipi_library.columns:
            seq_col = candidate
            break
    if seq_col is None:
        return {}
    lookup = {}
    for _, row in ipi_library.iterrows():
        name = str(row.get("Name", "")).strip()
        if not name or _has_stop_codon(row[seq_col]):
            continue
        lookup[name] = _clean_sequence(row[seq_col])
    return lookup


def _is_full_heavy_sequence(seq: str, cdr3_value: str) -> bool:
    if _has_stop_codon(seq) or _has_stop_codon(cdr3_value):
        return False
    seq = _clean_sequence(seq)
    cdr3_value = _clean_cdr3_sequence(cdr3_value)
    if not seq or not cdr3_value:
        return False
    return len(seq) >= FULL_HEAVY_MIN_AA_LEN and FR4 in seq and cdr3_value in seq


def _is_full_light_sequence(seq: str) -> bool:
    if _has_stop_codon(seq):
        return False
    return len(_clean_sequence(seq)) >= 80


def _assemble_hseq_from_regions(row: pd.Series, cdr3_value: str) -> str:
    pieces: list[str] = []
    for col in ["FR1", "CDR1", "FR2", "CDR2", "FR3"]:
        if _has_stop_codon(row.get(col, "")):
            return ""
        piece = _clean_sequence(row.get(col, ""))
        if not piece:
            return ""
        pieces.append(piece)
    if _has_stop_codon(row.get("CDR3", "")) or _has_stop_codon(cdr3_value) or _has_stop_codon(row.get("FR4", "")):
        return ""
    cdr3_piece = _clean_sequence(row.get("CDR3", "")) or _clean_sequence(cdr3_value)
    fr4 = _clean_sequence(row.get("FR4", ""))
    if not cdr3_piece or not fr4:
        return ""
    assembled = "".join([*pieces, cdr3_piece, fr4])
    return assembled if _is_full_heavy_sequence(assembled, cdr3_piece) else ""


def _normalize_scaffold(value: Any, prefix: str = "V") -> str:
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if not value or value == "UNK":
        return value
    return value if value.startswith(prefix) else prefix + value


def _scaffold_lookup_candidates(value: Any) -> list[str]:
    if pd.isna(value):
        return []
    raw = str(value).strip()
    if not raw or raw == "UNK":
        return []

    candidates = [raw]
    if raw.startswith("V") and len(raw) > 1:
        candidates.append(raw[1:])
    else:
        candidates.append(f"V{raw}")

    # Current library CSVs may use H1-69/K3-20 while older files used
    # VH1-69/VK3-20. Try both forms without changing the reported scaffold.
    if raw.startswith(("H", "K", "L")):
        candidates.append(f"V{raw}")
    if raw.startswith(("VH", "VK", "VL")) and len(raw) > 2:
        candidates.append(raw[1:])
    return list(dict.fromkeys(candidates))


def _lookup_scaffold_sequence(vh_lookup: dict[str, str], value: Any) -> tuple[str, str]:
    for candidate in _scaffold_lookup_candidates(value):
        seq = vh_lookup.get(candidate, "")
        if seq:
            return seq, candidate
    return "", ""


def _row_uses_heavy_only_sequence(row: pd.Series, force_no_lseq: bool) -> bool:
    row_library_type = str(row.get("library_type", "") or "").strip().lower()
    row_library_name = str(row.get("library", "") or "").strip().lower()
    return (
        force_no_lseq
        or row_library_type in {"vhh", "vhh_full"}
        or row_library_name in {"vhh", "vhh_full"}
    )


def Generate_FullVHVL(
    ngs_leads: pd.DataFrame,
    *,
    model_dir: str | Path | None = None,
    library_csv: str | Path | None = None,
    heavy_scaffold: str = "vh_scaffold",
    light_scaffold: str = "vl_scaffold",
    cdr3: str = "cdr3_aa",
    allow_cdr3_fallback: bool = False,
    force_no_lseq: bool = False,
    append_light_constant: bool = False,
) -> pd.DataFrame:
    """Generate Fab/VHH HSEQ and LSEQ columns from IPI scaffold calls.

    Fab rows intentionally follow the legacy IPI report logic:
    HSEQ = VH scaffold AA + CDR3 without the leading CAR cysteine + FR4,
    LSEQ = VL scaffold AA. ANARCI chain/FR-region assembly is only used for
    heavy-only VHH rows where no light chain is expected.
    """
    model_dir = _resolve_model_dir(model_dir)
    ipi_library = _load_ipi_library(model_dir, library_csv)
    vh_lookup = _library_sequence_lookup(ipi_library)

    df = ngs_leads.copy().reset_index(drop=True)
    if cdr3 not in df.columns:
        raise ValueError(f"Missing required CDR3 column: {cdr3}")
    if heavy_scaffold not in df.columns:
        df[heavy_scaffold] = ""
    if light_scaffold not in df.columns:
        df[light_scaffold] = ""

    hseq: list[str] = []
    lseq: list[str] = []
    hseq_source: list[str] = []
    lseq_source: list[str] = []
    for _, row in df.iterrows():
        cdr3_has_stop = (
            _has_stop_codon(row.get(cdr3, ""))
            or _has_stop_codon(row.get("CDR3", ""))
            or _has_stop_codon(row.get("cdr3_aa", ""))
        )
        cdr3_value = "" if cdr3_has_stop else _clean_cdr3_sequence(row[cdr3])

        heavy_base, heavy_key = _lookup_scaffold_sequence(vh_lookup, row[heavy_scaffold])
        light_base, light_key = _lookup_scaffold_sequence(vh_lookup, row[light_scaffold])
        row_force_no_lseq = _row_uses_heavy_only_sequence(row, force_no_lseq)

        if not cdr3_value:
            hseq.append("")
            hseq_source.append("stop_codon_cdr3" if cdr3_has_stop else "missing_cdr3")
            if row_force_no_lseq:
                lseq.append("")
                lseq_source.append("not_applicable")
            else:
                light_seq = light_base
                if light_seq and append_light_constant:
                    light_seq = f"{light_seq}{FAB_LIGHT_CONSTANT}"
                lseq.append(light_seq)
                if light_seq:
                    source = f"scaffold_library:{light_key}"
                    lseq_source.append(f"{source}_plus_constant" if append_light_constant else source)
                else:
                    lseq_source.append("missing")
            continue

        if row_force_no_lseq:
            existing_hseq = _clean_sequence(row.get("HSEQ", ""))
            chain_hseq = _clean_sequence(row.get("CHAIN", ""))
            if _is_full_heavy_sequence(existing_hseq, cdr3_value):
                heavy_seq = existing_hseq
                heavy_source = "existing_hseq"
            elif _is_full_heavy_sequence(chain_hseq, cdr3_value):
                heavy_seq = chain_hseq
                heavy_source = "anarci_chain"
            else:
                heavy_seq = _assemble_hseq_from_regions(row, cdr3_value)
                heavy_source = "fr_cdr_regions" if heavy_seq else ""
            if not heavy_seq:
                if heavy_base:
                    heavy_seq = f"{heavy_base}{cdr3_value}{FR4}"
                    heavy_source = "scaffold_library"
                elif allow_cdr3_fallback:
                    heavy_seq = cdr3_value
                    heavy_source = "cdr3_only"
                else:
                    heavy_seq = ""
                    heavy_source = "incomplete"
            lseq.append("")
            lseq_source.append("not_applicable")
        else:
            if heavy_base:
                heavy_seq = f"{heavy_base}{cdr3_value}{FR4}"
                heavy_source = f"scaffold_library:{heavy_key}"
            elif allow_cdr3_fallback:
                heavy_seq = cdr3_value
                heavy_source = "cdr3_only"
            else:
                heavy_seq = ""
                heavy_source = "incomplete"

            light_seq = light_base
            if light_seq and append_light_constant:
                light_seq = f"{light_seq}{FAB_LIGHT_CONSTANT}"
            lseq.append(light_seq)
            if light_seq:
                source = f"scaffold_library:{light_key}"
                lseq_source.append(f"{source}_plus_constant" if append_light_constant else source)
            else:
                lseq_source.append("missing")
        hseq.append(heavy_seq)
        hseq_source.append(heavy_source)

    df["HSEQ"] = hseq
    df["LSEQ"] = lseq
    df["HSEQ_SOURCE"] = hseq_source
    df["LSEQ_SOURCE"] = lseq_source
    if "BARCODE" not in df.columns:
        df["BARCODE"] = df.index
    return df


def _batch_loader(data, batch_size):
    for start in range(0, len(data), batch_size):
        end = min(start + batch_size, len(data))
        yield start, end, data[start:end]


def _insert_space_every_other_except_cls(input_string):
    parts = input_string.split("[CLS]")
    modified_parts = ["".join([char + " " for char in part]).strip() for part in parts]
    return " [CLS] ".join(modified_parts)


@lru_cache(maxsize=2)
def _antiberta2_components(model_name: str = DEFAULT_ANTIBERTA2_MODEL):
    from transformers import RoFormerForMaskedLM, RoFormerTokenizer

    tokenizer = RoFormerTokenizer.from_pretrained(model_name)
    model = RoFormerForMaskedLM.from_pretrained(model_name)
    model.eval()
    return tokenizer, model


def Generate_Antiberta2_Embedding_HLSEQ(
    df: pd.DataFrame,
    *,
    model_name: str = DEFAULT_ANTIBERTA2_MODEL,
    batch_size: int = 128,
    device: str = "cpu",
) -> pd.DataFrame:
    import torch

    dat = df.copy()
    linker = "[CLS][CLS]"
    dat["HL"] = np.where(dat["LSEQ"].astype(str).str.len() > 0, dat["HSEQ"] + linker + dat["LSEQ"], dat["HSEQ"])
    sequences = (
        dat["HL"]
        .astype(str)
        .apply(_insert_space_every_other_except_cls)
        .str.replace("  ", " ", regex=False)
        .values
    )
    tokenizer, model = _antiberta2_components(model_name)
    torch_device = torch.device(device)
    model.to(torch_device)

    max_length = 256
    embeddings = torch.empty((len(sequences), 1024), device="cpu")
    n_batches = math.ceil(len(sequences) / batch_size) if len(sequences) else 0
    for batch_idx, (start, end, batch) in enumerate(_batch_loader(sequences, batch_size), start=1):
        print(f"AntiBERTa2 batch {batch_idx}/{n_batches}")
        encoded = [
            tokenizer.encode(
                seq,
                padding="max_length",
                truncation=True,
                max_length=max_length,
                return_special_tokens_mask=True,
            )
            for seq in batch
        ]
        x = torch.tensor(encoded, device=torch_device)
        attention_mask = (x != tokenizer.pad_token_id).float()
        with torch.no_grad():
            outputs = model(x, attention_mask=attention_mask, output_hidden_states=True).hidden_states[-1]
        pooled = []
        for seq_output, mask in zip(outputs, attention_mask):
            pooled.append(seq_output[mask == 1, :].mean(0).detach().cpu())
        embeddings[start:end] = torch.stack(pooled)

    dat["antiberta2-cssp_emb"] = embeddings.tolist()
    return dat


def _load_model_file(path: Path):
    try:
        import joblib

        return joblib.load(path)
    except Exception as joblib_exc:
        try:
            from tensorflow import keras

            return keras.models.load_model(path)
        except Exception as keras_exc:
            raise RuntimeError(f"Could not load model {path}: joblib={joblib_exc}; keras={keras_exc}") from keras_exc


def _model_predict(model, values):
    try:
        return model.predict(values, verbose=0)
    except TypeError:
        return model.predict(values)


def CNN_PRED_PSR(dat: pd.DataFrame, model_dir: Path) -> pd.DataFrame:
    emb = pd.DataFrame(dat["antiberta2-cssp_emb"].to_list(), index=dat.index)
    model_path = model_dir / "model_june_2025/CNN_MODEL_ANTIBERTA2-CSSP_IPI_ELISA_NGS.keras"
    if not model_path.exists():
        raise FileNotFoundError(f"CNN PSR model not found: {model_path}")
    model = _load_model_file(model_path)
    pred_psr = _model_predict(model, emb.values)
    dat["psr_cnn_pred_proba"] = np.asarray(pred_psr).reshape(-1)
    dat["psr_cnn_pred"] = (dat["psr_cnn_pred_proba"] > 0.5).astype(int)
    return dat


def XGBoost_PRED_PSR(dat: pd.DataFrame, model_dir: Path) -> pd.DataFrame:
    emb = pd.DataFrame(dat["antiberta2-cssp_emb"].to_list(), index=dat.index)
    model_path = model_dir / "model_june_2025/XGBOOST_MODEL_ANTIBERTA2-CSSP_IPI_ELISA_NGS.keras"
    if not model_path.exists():
        raise FileNotFoundError(f"XGBoost PSR model not found: {model_path}")
    model = _load_model_file(model_path)
    pred_psr = model.predict(emb.values)
    pred_psr_prob = model.predict_proba(emb.values)[:, 1]
    dat["psr_xgboost_pred"] = pred_psr
    dat["psr_xgboost_pred_proba"] = pred_psr_prob
    return dat


def PredictionConsensus(dat: pd.DataFrame) -> pd.DataFrame:
    proba_cols = [col for col in ["psr_cnn_pred_proba", "psr_xgboost_pred_proba"] if col in dat.columns]
    if not proba_cols:
        raise RuntimeError("No PSR probability columns were produced.")
    dat["psr_pred_mean_proba"] = dat[proba_cols].mean(axis=1)
    dat["psr_pred_mean"] = (dat["psr_pred_mean_proba"] >= 0.5).astype(int)
    return dat


def predict_psr_on_clone(
    input_file: str | Path,
    *,
    model_dir: str | Path | None = None,
    output_file: str | Path | None = None,
    include_embeddings: bool = False,
    batch_size: int = 128,
    device: str = "cpu",
) -> dict:
    input_file = Path(input_file)
    output_file = Path(output_file) if output_file else _default_output_path(input_file, "psr_predicted")
    model_dir = _resolve_model_dir(model_dir)

    print(f"Processing PSR prediction: {input_file}")
    df = _read_table(input_file)
    if df.empty:
        _write_table(df, output_file)
        return {"input": str(input_file), "output": str(output_file), "rows": 0, "pass": 0, "fail": 0}

    df = Generate_FullVHVL(df, model_dir=model_dir)
    df = Generate_Antiberta2_Embedding_HLSEQ(df, batch_size=batch_size, device=device)

    errors = []
    for predictor in [CNN_PRED_PSR, XGBoost_PRED_PSR]:
        try:
            df = predictor(df, model_dir)
        except Exception as exc:
            errors.append(str(exc))
            print(f"Warning: {exc}")

    if "psr_cnn_pred_proba" not in df.columns and "psr_xgboost_pred_proba" not in df.columns:
        raise RuntimeError("All PSR predictors failed: " + " | ".join(errors))

    df = PredictionConsensus(df)
    if not include_embeddings:
        df = df.drop(columns=["antiberta2-cssp_emb", "HL"], errors="ignore")

    _write_table(df, output_file)
    passed = int((df["psr_pred_mean"] == 1).sum())
    failed = int((df["psr_pred_mean"] == 0).sum())
    print(f"Predictions saved to: {output_file}")
    print(f"PASS: {passed} | FAIL: {failed}")
    return {"input": str(input_file), "output": str(output_file), "rows": len(df), "pass": passed, "fail": failed}


def _delph_pythonpath_candidates(delph_config: dict | None = None) -> list[Path]:
    delph_config = delph_config or {}
    configured_paths = []
    for key in ["package_path", "pythonpath"]:
        value = delph_config.get(key)
        if not value:
            continue
        if isinstance(value, (list, tuple)):
            configured_paths.extend(value)
        else:
            configured_paths.append(value)

    candidates: list[Path] = []
    for value in configured_paths:
        path = Path(value)
        candidates.extend([path, path / "src", path / "python"])
        if path.name in {"delphi", "delph"}:
            candidates.append(path.parent)
    return candidates


def configure_delph_pythonpath(delph_config: dict | None = None) -> list[str]:
    added: list[str] = []
    for path in _delph_pythonpath_candidates(delph_config):
        if not path.exists():
            continue
        path_str = str(path)
        if path_str not in sys.path:
            sys.path.insert(0, path_str)
            added.append(path_str)
    return added


def delph_import_available(delph_config: dict | None = None) -> bool:
    delph_config = delph_config or {}
    if delph_config.get("command"):
        return True
    configure_delph_pythonpath(delph_config)
    module_name = delph_config.get("module")
    module_names = [module_name] if module_name else []
    module_names.extend(["delphi", "delph"])
    return any(name and importlib.util.find_spec(name) is not None for name in module_names)


def _import_delph_module(delph_config: dict | None = None):
    delph_config = delph_config or {}
    configure_delph_pythonpath(delph_config)
    module_name = delph_config.get("module")
    candidates = [module_name] if module_name else []
    candidates.extend(["delphi", "delph"])
    for name in candidates:
        if not name:
            continue
        try:
            return importlib.import_module(name)
        except ImportError:
            continue
    searched = ", ".join(str(path) for path in _delph_pythonpath_candidates(delph_config))
    raise ImportError(
        "Could not import 'delphi' or 'delph'. "
        f"Check ml.delph.package_path/module in config.yaml. Searched paths: {searched or '(sys.path only)'}"
    )


def _run_delph_command(input_file: Path, output_file: Path, command_template: str) -> dict:
    command = command_template.format(
        input=input_file,
        output=output_file,
        input_q=shlex.quote(str(input_file)),
        output_q=shlex.quote(str(output_file)),
    )
    result = subprocess.run(shlex.split(command), capture_output=True, text=True)
    if result.returncode != 0:
        exc = RuntimeError(f"Delphi command failed: {result.stderr.strip()}")
        exc.stdout = result.stdout.strip()
        exc.stderr = result.stderr.strip()
        raise exc
    stdout = result.stdout.strip()
    saved_match = re.search(r"Saved predictions to:\s*(.+)", stdout)
    actual_output = saved_match.group(1).strip() if saved_match else str(output_file)
    return {
        "input": str(input_file),
        "output": actual_output,
        "stdout": stdout,
        "prediction_prefix": _prediction_prefix_from_command(command_template),
    }


def _prediction_prefix_from_command(command_template: str) -> str:
    match = re.search(r"--target\s+([A-Za-z0-9_-]+)", command_template)
    if not match:
        return ""
    target = match.group(1).strip().lower()
    return target.removesuffix("_filter")


def _run_delph_commands(
    input_file: Path,
    output_file: Path,
    command_templates: list[str],
    *,
    continue_on_error: bool = False,
) -> dict:
    results = []
    failures = []
    for idx, command_template in enumerate(command_templates, start=1):
        command_output = output_file
        if len(command_templates) > 1:
            command_output = output_file.with_name(f"{output_file.stem}_{idx}{output_file.suffix}")
        try:
            command_result = _run_delph_command(input_file, command_output, command_template)
            command_result["command"] = command_template
            results.append(command_result)
        except Exception as exc:
            failure = {
                "command": command_template,
                "error": str(exc),
                "stdout": getattr(exc, "stdout", ""),
                "stderr": getattr(exc, "stderr", ""),
            }
            failures.append(failure)
            if not continue_on_error:
                raise
            print(f"  Delphi command skipped/failed: {exc}")

    if not results:
        error_text = " | ".join(f"{failure['command']}: {failure['error']}" for failure in failures)
        raise RuntimeError(f"All Delphi commands failed. {error_text}")

    if len(results) == 1:
        result = results[0]
        result["outputs"] = [result.get("output", "")]
        if failures:
            result["failed_commands"] = failures
        return result
    return {
        "input": str(input_file),
        "output": ", ".join(result.get("output", "") for result in results),
        "outputs": [result.get("output", "") for result in results],
        "commands": results,
        "failed_commands": failures,
    }


def _delph_prediction_column(column: str) -> bool:
    model = r"(?:rf|xgboost|cnn|transformer(?:_lm|_onehot)?|delph)"
    return bool(re.match(rf"^(?:[a-z0-9]+_)*{model}_.+_(score|label|pred|pred_proba)$", column))


def _delph_score_columns(
    df: pd.DataFrame,
    columns: list[str] | None = None,
    *,
    label_prefix: str | None = None,
) -> list[str]:
    candidates = columns or list(df.columns)
    score_cols: list[str] = []
    for col in candidates:
        col = str(col)
        lowered = col.lower()
        if lowered in {"mean_delphi_score", "mean_psr_score", "mean_sec_score"}:
            continue
        if label_prefix and not lowered.startswith(f"{label_prefix.lower()}_"):
            continue
        if not lowered.endswith("_score"):
            continue
        if not _delph_prediction_column(col):
            continue
        if "xgboost" in lowered or "transformer" in lowered:
            score_cols.append(col)
    return [col for col in dict.fromkeys(score_cols) if col in df.columns]


def _mean_delphi_score(
    df: pd.DataFrame,
    columns: list[str] | None = None,
    *,
    label_prefix: str | None = None,
) -> pd.Series:
    score_cols = _delph_score_columns(df, columns, label_prefix=label_prefix)
    if not score_cols:
        return pd.Series(np.nan, index=df.index)
    scores = df[score_cols].apply(pd.to_numeric, errors="coerce")
    return scores.mean(axis=1)


def _charge_at_ph(seq: str, ph: float) -> float:
    counts = {aa: seq.count(aa) for aa in "CDEHKRY"}
    positive = (
        1 / (1 + 10 ** (ph - 7.5))
        + counts["K"] / (1 + 10 ** (ph - 10.5))
        + counts["R"] / (1 + 10 ** (ph - 12.5))
        + counts["H"] / (1 + 10 ** (ph - 6.0))
    )
    negative = (
        1 / (1 + 10 ** (3.55 - ph))
        + counts["D"] / (1 + 10 ** (3.9 - ph))
        + counts["E"] / (1 + 10 ** (4.1 - ph))
        + counts["C"] / (1 + 10 ** (8.3 - ph))
        + counts["Y"] / (1 + 10 ** (10.1 - ph))
    )
    return positive - negative


def _sequence_isoelectric_point(value: Any) -> float:
    seq = _clean_sequence(value)
    if not seq:
        return float("nan")
    low, high = 0.0, 14.0
    for _ in range(40):
        mid = (low + high) / 2
        if _charge_at_ph(seq, mid) > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2, 3)


def _sequence_net_charge(value: Any) -> float:
    seq = _clean_sequence(value)
    if not seq:
        return float("nan")
    return float(seq.count("K") + seq.count("R") + seq.count("H") - seq.count("D") - seq.count("E"))


def _add_sequence_biophysical_columns(df: pd.DataFrame) -> pd.DataFrame:
    cdr3_col = "CDR3" if "CDR3" in df.columns else "cdr3_aa" if "cdr3_aa" in df.columns else ""
    if cdr3_col:
        cdr3_seq = df[cdr3_col].apply(_clean_cdr3_sequence)
        df["cdr3_arg_count"] = cdr3_seq.str.count("R").astype(float)
        df["cdr3_asp_count"] = cdr3_seq.str.count("D").astype(float)
        df["cdr3_trp_count"] = cdr3_seq.str.count("W").astype(float)
        df["cdr3_glu_count"] = cdr3_seq.str.count("E").astype(float)
        df["cdr3_length"] = cdr3_seq.str.len().astype(float)
        df["cdr3_net_charge"] = cdr3_seq.apply(_sequence_net_charge)
        df["cdr3_isoelectric_point"] = cdr3_seq.apply(_sequence_isoelectric_point)
    if "HSEQ" in df.columns:
        df["hseq_isoelectric_point"] = df["HSEQ"].apply(_sequence_isoelectric_point)
    return df


def _refresh_final_workbook_with_predictions(
    source_input: Path,
    prepared_input: Path,
    prediction_outputs: list[Path],
    *,
    prediction_prefixes: dict[Path, str] | None = None,
    exclude_columns: list[str] | None = None,
) -> dict:
    final_df = _read_table(source_input)
    prepared_df = _read_table(prepared_input)

    if "BARCODE" not in final_df.columns:
        final_df["BARCODE"] = final_df.index.astype(str)
    final_df["BARCODE"] = final_df["BARCODE"].astype(str).str.strip()

    if "BARCODE" not in prepared_df.columns:
        prepared_df["BARCODE"] = prepared_df.index.astype(str)
    prepared_df["BARCODE"] = prepared_df["BARCODE"].astype(str).str.strip()

    old_prediction_cols = [col for col in final_df.columns if _delph_prediction_column(str(col))]
    old_prediction_cols.extend(
        [
            col
            for col in [
                "ML_SEQUENCE_OK",
                "mean_psr_score",
                "mean_sec_score",
                "mean_delphi_score",
                *SEQUENCE_BIOPHYSICAL_COLUMNS,
            ]
            if col in final_df.columns
        ]
    )
    if old_prediction_cols:
        final_df = final_df.drop(columns=list(dict.fromkeys(old_prediction_cols)))

    prepared_indexed = prepared_df.drop_duplicates("BARCODE").set_index("BARCODE")
    final_df["ML_SEQUENCE_OK"] = final_df["BARCODE"].isin(prepared_indexed.index)
    for col in ["HSEQ", "LSEQ", "HSEQ_SOURCE", "LSEQ_SOURCE"]:
        if col not in prepared_indexed.columns:
            continue
        mapped = final_df["BARCODE"].map(prepared_indexed[col])
        if col not in final_df.columns:
            final_df[col] = ""
        final_df.loc[mapped.notna(), col] = mapped[mapped.notna()]

    merged_columns: list[str] = []
    prediction_prefixes = prediction_prefixes or {}
    for output_path in prediction_outputs:
        if not output_path.exists():
            continue
        pred_df = _read_table(output_path)
        if "BARCODE" not in pred_df.columns and len(pred_df) == len(prepared_df):
            pred_df["BARCODE"] = prepared_df["BARCODE"].values
        if "BARCODE" not in pred_df.columns:
            print(f"  Warning: prediction output has no BARCODE and cannot be merged: {output_path}")
            continue
        pred_df["BARCODE"] = pred_df["BARCODE"].astype(str).str.strip()
        pred_indexed = pred_df.drop_duplicates("BARCODE").set_index("BARCODE")
        pred_cols = [col for col in pred_indexed.columns if _delph_prediction_column(str(col))]
        prefix = prediction_prefixes.get(output_path.resolve(), "").strip()
        for col in pred_cols:
            out_col = str(col)
            if prefix and not out_col.startswith(f"{prefix}_"):
                out_col = f"{prefix}_{out_col}"
            final_df[out_col] = final_df["BARCODE"].map(pred_indexed[col])
            merged_columns.append(out_col)

    final_df["mean_psr_score"] = _mean_delphi_score(final_df, merged_columns, label_prefix="psr")
    final_df["mean_sec_score"] = _mean_delphi_score(final_df, merged_columns, label_prefix="sec")
    final_df["mean_delphi_score"] = _mean_delphi_score(final_df, merged_columns)
    final_df = _add_sequence_biophysical_columns(final_df)

    rows = len(final_df)
    ml_rows = int(final_df["ML_SEQUENCE_OK"].sum())
    drop_cols = [
        col
        for col in (DEFAULT_FINAL_EXCLUDE_COLUMNS if exclude_columns is None else exclude_columns)
        if col in final_df.columns
    ]
    if drop_cols:
        final_df = final_df.drop(columns=drop_cols, errors="ignore")
    _write_table(final_df, source_input)
    return {
        "final_output": str(source_input),
        "rows": rows,
        "ml_rows": ml_rows,
        "dropped_rows": max(rows - ml_rows, 0),
        "merged_prediction_columns": list(dict.fromkeys(merged_columns)),
    }


def _cleanup_delph_artifacts(source_input: Path, prepared_input: Path, result: dict) -> dict:
    parent = source_input.parent
    keep = {source_input.resolve()}
    paths: list[Path] = []

    if prepared_input.resolve() not in keep:
        paths.append(prepared_input)
    for output in [*result.get("prediction_outputs", []), *result.get("outputs", [])]:
        if output:
            paths.append(Path(output))
    if result.get("output"):
        for value in str(result["output"]).split(","):
            value = value.strip()
            if value:
                paths.append(Path(value))

    command_results = result.get("commands") or [result]
    for command_result in command_results:
        stdout = str(command_result.get("stdout", ""))
        paths.extend(_extract_delph_generated_paths(stdout, parent))
    for failure in result.get("failed_commands", []):
        paths.extend(_extract_delph_generated_paths(str(failure.get("stdout", "")), parent))
        paths.extend(_extract_delph_generated_paths(str(failure.get("stderr", "")), parent))

    paths.extend(parent.glob(f"{prepared_input.name}.*.emb.csv"))
    unique_paths = []
    for path in paths:
        if not _path_is_within(path, parent):
            continue
        if path.resolve() in keep:
            continue
        if path not in unique_paths:
            unique_paths.append(path)

    deleted = []
    errors = []
    for path in sorted(unique_paths, key=lambda p: len(str(p)), reverse=True):
        try:
            if path.is_dir():
                shutil.rmtree(path)
                deleted.append(str(path))
            elif path.exists():
                path.unlink()
                deleted.append(str(path))
        except Exception as exc:
            errors.append({"path": str(path), "error": str(exc)})
    return {"deleted": deleted, "errors": errors}


def cleanup_existing_delph_artifacts(source_input: Path) -> dict:
    parent = source_input.parent
    stem = source_input.stem
    paths: list[Path] = []
    patterns = [
        f"{stem}_delphi_input*",
        f"{stem}_delphi_dropped_sequences*",
        f"{stem}_delph_predicted*",
        f"predict_{stem}_delphi_input*",
    ]
    for pattern in patterns:
        paths.extend(parent.glob(pattern))

    keep = {source_input.resolve()}
    deleted = []
    errors = []
    for path in sorted(set(paths), key=lambda p: len(str(p)), reverse=True):
        try:
            if not _path_is_within(path, parent) or path.resolve() in keep:
                continue
            if path.is_dir():
                shutil.rmtree(path)
                deleted.append(str(path))
            elif path.exists():
                path.unlink()
                deleted.append(str(path))
        except Exception as exc:
            errors.append({"path": str(path), "error": str(exc)})
    return {"deleted": deleted, "errors": errors}


def _positive_int_config(config: dict, key: str) -> int | None:
    value = config.get(key)
    if value in (None, ""):
        return None
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        print(f"  Warning: ignoring invalid Delphi {key}={value!r}; expected integer")
        return None
    return parsed if parsed > 0 else None


def _float_config(config: dict, key: str) -> float | None:
    value = config.get(key)
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        print(f"  Warning: ignoring invalid Delphi {key}={value!r}; expected number")
        return None


def _limit_delph_prediction_rows(df: pd.DataFrame, effective_config: dict) -> tuple[pd.DataFrame, bool]:
    max_rows = _positive_int_config(effective_config, "max_prediction_rows")
    min_max_freq = _float_config(effective_config, "min_prediction_max_freq")
    if max_rows is None and min_max_freq is None:
        return df, False

    if df.empty:
        return df, False

    rank_col = str(effective_config.get("prediction_rank_column") or "max_freq")
    keep = pd.Series(False, index=df.index)
    criteria: list[str] = []
    rank_values = None

    if rank_col in df.columns:
        rank_values = pd.to_numeric(df[rank_col], errors="coerce")
        if min_max_freq is not None:
            keep |= rank_values.fillna(float("-inf")) >= min_max_freq
            criteria.append(f"{rank_col} >= {min_max_freq:g}")
        if max_rows is not None:
            top_index = rank_values.fillna(float("-inf")).sort_values(ascending=False).head(max_rows).index
            keep.loc[top_index] = True
            criteria.append(f"top {max_rows} by {rank_col}")
    else:
        if min_max_freq is not None:
            print(f"  Warning: cannot apply min_prediction_max_freq; missing {rank_col} column")
        if max_rows is not None:
            keep.loc[df.head(max_rows).index] = True
            criteria.append(f"first {max_rows} rows")

    limited = df.loc[keep].copy()
    if limited.empty:
        fallback_rows = max_rows or len(df)
        limited = df.head(fallback_rows).copy()
        criteria.append(f"fallback first {len(limited)} rows")

    if len(limited) == len(df):
        return df, False

    print(
        f"Limited Delphi prediction input: {len(df)} -> {len(limited)} rows "
        f"({'; '.join(criteria)})"
    )
    return limited, True


def _write_sequence_columns_to_input(source_input: Path, sequence_df: pd.DataFrame) -> None:
    sequence_cols = [col for col in ["HSEQ", "LSEQ"] if col in sequence_df.columns]
    if not sequence_cols:
        return

    try:
        final_df = _read_table(source_input)
    except Exception as exc:
        print(f"  Warning: could not write sequence columns to {source_input.name}: {exc}")
        return

    if "BARCODE" not in final_df.columns:
        final_df["BARCODE"] = final_df.index.astype(str)
    final_df["BARCODE"] = final_df["BARCODE"].astype(str).str.strip()

    seq_df = sequence_df.copy()
    if "BARCODE" not in seq_df.columns:
        seq_df["BARCODE"] = seq_df.index.astype(str)
    seq_df["BARCODE"] = seq_df["BARCODE"].astype(str).str.strip()
    seq_indexed = seq_df.drop_duplicates("BARCODE").set_index("BARCODE")

    changed = False
    for col in sequence_cols:
        mapped = final_df["BARCODE"].map(seq_indexed[col])
        if col not in final_df.columns:
            final_df[col] = ""
            changed = True
        values = mapped.fillna("")
        update_mask = values.astype(str).str.strip().ne("")
        if update_mask.any():
            final_df.loc[update_mask, col] = mapped[update_mask]
            changed = True

    if changed:
        _write_table(final_df, source_input)
        print(f"  Wrote sequence columns to final workbook before Delphi prediction: {source_input.name}")


def _prepare_delph_command_input(input_file: Path, delph_config: dict) -> tuple[Path, dict]:
    base_config = delph_config or {}
    if base_config.get("prepare_input", True) is False:
        df = _read_table(input_file)
        df, _ = _enrich_from_clone_table(input_file, df)
        library_type = _infer_library_type_from_df(df)
        return input_file, _merge_library_type_delph_config(base_config, library_type)

    df = _read_table(input_file)
    df, enriched_from_clone = _enrich_from_clone_table(input_file, df)
    library_type = _infer_library_type_from_df(df)
    effective_config = _merge_library_type_delph_config(base_config, library_type)
    changed = enriched_from_clone
    old_prediction_cols = [col for col in df.columns if _delph_prediction_column(str(col))]
    if old_prediction_cols:
        df = df.drop(columns=old_prediction_cols)
        changed = True
    cdr3_column = effective_config.get("cdr3_column") or ("CDR3" if "CDR3" in df.columns else "cdr3_aa")
    if cdr3_column not in df.columns and cdr3_column == "CDR3" and "cdr3_aa" in df.columns:
        df["CDR3"] = df["cdr3_aa"]
        changed = True
    if cdr3_column not in df.columns and cdr3_column == "cdr3_aa" and "CDR3" in df.columns:
        df["cdr3_aa"] = df["CDR3"]
        changed = True
    force_no_lseq = bool(effective_config.get("force_no_lseq", library_type == "vhh"))

    if cdr3_column in df.columns:
        df = Generate_FullVHVL(
            df,
            model_dir=effective_config.get("model_dir"),
            library_csv=effective_config.get("library_csv"),
            heavy_scaffold=effective_config.get("heavy_scaffold", "vh_scaffold"),
            light_scaffold=effective_config.get("light_scaffold", "vl_scaffold"),
            cdr3=cdr3_column,
            allow_cdr3_fallback=effective_config.get("allow_cdr3_fallback", False),
            force_no_lseq=force_no_lseq,
            append_light_constant=effective_config.get("append_light_constant", False),
        )
        changed = True
        if effective_config.get("write_sequence_columns_to_input", True):
            _write_sequence_columns_to_input(input_file, df)
    elif force_no_lseq:
        df["LSEQ"] = ""
        df["LSEQ_SOURCE"] = "not_applicable"
        changed = True

    require_full_hseq = effective_config.get("require_full_hseq", True)
    require_lseq = effective_config.get("require_lseq", True)
    drop_incomplete = effective_config.get("drop_incomplete_sequences", True)
    write_dropped_audit = bool(effective_config.get("write_dropped_audit", True))
    if cdr3_column in df.columns and (require_full_hseq or require_lseq):
        hseq_ok = (
            df.apply(lambda row: _is_full_heavy_sequence(row.get("HSEQ", ""), row.get(cdr3_column, "")), axis=1)
            if require_full_hseq
            else pd.Series(True, index=df.index)
        )
        lseq_ok = (
            df["LSEQ"].apply(_is_full_light_sequence)
            if require_lseq and "LSEQ" in df.columns
            else pd.Series(True, index=df.index)
        )
        df["_ml_sequence_ok"] = hseq_ok & lseq_ok
        if drop_incomplete and not bool(df["_ml_sequence_ok"].all()):
            dropped = df.loc[~df["_ml_sequence_ok"]].copy()
            dropped_path = _default_output_path(input_file, "delphi_dropped_sequences")
            if write_dropped_audit:
                _write_table(dropped, dropped_path)
            df = df.loc[df["_ml_sequence_ok"]].copy()
            requirements = ["full HSEQ"] if require_full_hseq else []
            if require_lseq:
                requirements.append("full LSEQ")
            requirement_text = " and ".join(requirements) or "required sequences"
            audit_text = f" Audit file: {dropped_path}" if write_dropped_audit else ""
            print(f"Dropped {len(dropped)} rows from Delphi input because {requirement_text} could not be generated.{audit_text}")
            if df.empty:
                audit_text = f" See dropped-sequence audit file: {dropped_path}" if write_dropped_audit else ""
                raise ValueError(
                    f"No rows remain for Delphi prediction after {requirement_text} validation. "
                    f"{audit_text}"
                )
            changed = True

    if "CDR3" not in df.columns and "cdr3_aa" in df.columns:
        df["CDR3"] = df["cdr3_aa"].astype(str).str.removeprefix("C")
        changed = True

    if "BARCODE" not in df.columns:
        df["BARCODE"] = df.index.astype(str)
        changed = True

    df, limited = _limit_delph_prediction_rows(df, effective_config)
    if limited:
        changed = True

    if not changed:
        return input_file, effective_config

    prepared_input = _default_output_path(input_file, "delphi_input")
    _write_table(df, prepared_input)
    return prepared_input, effective_config


def predict_delph_on_clone(
    input_file: str | Path,
    *,
    output_file: str | Path | None = None,
    delph_config: dict | None = None,
) -> dict:
    input_file = Path(input_file)
    output_file = Path(output_file) if output_file else _default_output_path(input_file, "delph_predicted")
    delph_config = delph_config or {}

    initial_df = _read_table(input_file)
    initial_library_type = _infer_library_type_from_df(initial_df)
    initial_effective_config = _merge_library_type_delph_config(delph_config, initial_library_type)
    if initial_effective_config.get("cleanup_existing_outputs", False):
        cleanup_result = cleanup_existing_delph_artifacts(input_file)
        if cleanup_result["deleted"]:
            print(f"  Removed {len(cleanup_result['deleted'])} stale Delphi artifact(s) before prediction")
        if cleanup_result["errors"]:
            print(f"  Warning: stale cleanup had {len(cleanup_result['errors'])} error(s)")

    prepared_input, effective_config = _prepare_delph_command_input(input_file, delph_config)
    command_templates = _normalize_command_templates(effective_config.get("command"))
    if command_templates:
        result = _run_delph_commands(
            prepared_input,
            output_file,
            command_templates,
            continue_on_error=bool(effective_config.get("continue_on_command_error", False)),
        )
        result["source_input"] = str(input_file)
        result["prepared_input"] = str(prepared_input)
        prediction_outputs = [Path(path) for path in result.get("outputs", []) if path]
        if result.get("output") and not prediction_outputs:
            prediction_outputs = [Path(str(result["output"]))]
        result["prediction_outputs"] = [str(path) for path in prediction_outputs]
        if effective_config.get("write_back_to_input", False):
            command_results = result.get("commands") or [result]
            prediction_prefixes = {
                Path(str(command_result.get("output"))).resolve(): str(command_result.get("prediction_prefix", ""))
                for command_result in command_results
                if command_result.get("output")
            }
            merge_result = _refresh_final_workbook_with_predictions(
                input_file,
                prepared_input,
                prediction_outputs,
                prediction_prefixes=prediction_prefixes,
                exclude_columns=effective_config.get("final_exclude_columns"),
            )
            result.update(merge_result)
            result["output"] = merge_result["final_output"]
            print(
                f"  Prediction columns merged into final workbook: {merge_result['final_output']} "
                f"({len(merge_result['merged_prediction_columns'])} columns)"
            )
        if effective_config.get("cleanup_intermediate_outputs", False):
            cleanup_result = _cleanup_delph_artifacts(input_file, prepared_input, result)
            result["cleanup"] = cleanup_result
            if cleanup_result["deleted"]:
                print(f"  Cleaned {len(cleanup_result['deleted'])} Delphi intermediate artifact(s)")
            if cleanup_result["errors"]:
                print(f"  Warning: cleanup had {len(cleanup_result['errors'])} error(s)")
        if effective_config.get("force_no_lseq"):
            result["sequence_mode"] = "heavy_only"
        return result

    df = _read_table(prepared_input)
    module = _import_delph_module(effective_config)
    callable_names = effective_config.get(
        "callables",
        ["predict_dataframe", "predict_df", "score_dataframe", "score_df", "predict", "score"],
    )

    last_error = None
    for callable_name in callable_names:
        func = getattr(module, callable_name, None)
        if not callable(func):
            continue
        try:
            result = func(df.copy())
        except TypeError as exc:
            last_error = exc
            sequence_col = effective_config.get("sequence_column")
            if sequence_col is None:
                sequence_col = "HSEQ" if "HSEQ" in df.columns else "cdr3_aa"
            try:
                result = func(df[sequence_col].astype(str).tolist())
            except Exception as seq_exc:
                last_error = seq_exc
                continue
        except Exception as exc:
            last_error = exc
            continue

        if isinstance(result, pd.DataFrame):
            out_df = result
        else:
            out_df = df.copy()
            values = np.asarray(result)
            if values.ndim == 1 and len(values) == len(out_df):
                out_df["delph_score"] = values
            else:
                out_df["delph_result"] = str(result)
        _write_table(out_df, output_file)
        return {"input": str(input_file), "output": str(output_file), "rows": len(out_df), "backend": "delph"}

    guidance = (
        "Delphi package imported, but no supported callable was found. Configure ml.delph.command "
        "or ml.delph.module/ml.delph.callables in config.yaml."
    )
    if last_error:
        guidance += f" Last error: {last_error}"
    raise RuntimeError(guidance)
