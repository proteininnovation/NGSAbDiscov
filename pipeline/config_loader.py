from pathlib import Path

import yaml

MIXED_LIBRARY = "mixed"


def _as_path(value, base_dir: Path | None = None) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value)
    if not path.is_absolute() and base_dir is not None:
        path = base_dir / path
    return path


def _infer_library_type(name: str, library_cfg: dict) -> str:
    explicit = library_cfg.get("type") or library_cfg.get("library_type")
    if explicit:
        return str(explicit).lower()
    if "vhh" in name.lower():
        return "vhh"
    if not library_cfg.get("vl_barcodes"):
        return "vhh"
    return "fab"


def _first_library_by_type(libraries: dict, library_type: str) -> str | None:
    for name, library_cfg in libraries.items():
        if _infer_library_type(name, library_cfg) == library_type:
            return name
    return None


def _normalize_library_aliases(cfg: dict) -> dict[str, str]:
    libraries = cfg.get("libraries", {})
    aliases = {str(k).lower(): str(v) for k, v in (cfg.get("library_aliases") or {}).items()}

    fab_default = "standard_fab" if "standard_fab" in libraries else aliases.get("fab") or _first_library_by_type(libraries, "fab")
    vhh_default = aliases.get("vhh") or ("vhh_full" if "vhh_full" in libraries else _first_library_by_type(libraries, "vhh"))
    if fab_default:
        aliases["fab"] = fab_default
        aliases["fab_standard"] = fab_default
    if vhh_default:
        aliases.setdefault("vhh", vhh_default)
    aliases.setdefault(MIXED_LIBRARY, MIXED_LIBRARY)

    for name in libraries:
        aliases.setdefault(str(name).lower(), name)

    cfg["library_aliases"] = aliases
    return aliases


def resolve_library_name(value: str | None, cfg: dict, *, allow_mixed: bool = False) -> str | None:
    if value is None:
        return None
    requested = str(value).strip()
    if not requested:
        return None

    aliases = cfg.get("library_aliases") or _normalize_library_aliases(cfg)
    resolved = aliases.get(requested.lower(), requested)
    if resolved == MIXED_LIBRARY:
        if allow_mixed:
            return MIXED_LIBRARY
        known = ", ".join(sorted(cfg.get("libraries", {})))
        raise ValueError(f"Library 'mixed' is only valid for run selection. Known concrete libraries: {known}")
    if resolved not in cfg.get("libraries", {}):
        known_aliases = ", ".join(sorted(aliases))
        known_libraries = ", ".join(sorted(cfg.get("libraries", {})))
        raise ValueError(
            f"Unknown library '{requested}'. Known libraries: {known_libraries}. "
            f"Known aliases: {known_aliases}"
        )
    return resolved


def load_config(config_path="config.yaml"):
    with open(config_path) as f:
        cfg = yaml.safe_load(f) or {}

    cfg.setdefault("general", {})
    cfg.setdefault("processing", {})
    cfg.setdefault("combine", {})
    cfg.setdefault("pick_leads", {})
    cfg.setdefault("repeat_check", {})
    cfg.setdefault("fastp", {})
    cfg.setdefault("libraries", {})
    _normalize_library_aliases(cfg)

    base_dir = _as_path(cfg["general"].get("base_dir", ".")) or Path(".")
    cfg["general"]["base_dir"] = base_dir
    cfg["general"].setdefault("output_folder", "results")

    previous_db = cfg["general"].get("previous_antibodies_db", "data/All_mAb_20251106_FACS_BLI.xlsx")
    cfg["general"]["previous_antibodies_db"] = _as_path(previous_db, base_dir)

    current_library = cfg.get("current_library")
    if not current_library:
        current_library = next(iter(cfg["libraries"]), None)
        cfg["current_library"] = current_library
    if current_library:
        cfg["current_library"] = resolve_library_name(current_library, cfg, allow_mixed=True)

    for library_name, library_cfg in cfg["libraries"].items():
        library_cfg.setdefault("library_type", _infer_library_type(library_name, library_cfg))
        library_cfg.setdefault("vh_barcodes", {})
        if library_cfg["library_type"] == "fab":
            library_cfg.setdefault("vl_barcodes", {})

    cfg.setdefault("ml", {})
    cfg["ml"].setdefault("enabled", True)
    cfg["ml"].setdefault("backend", "auto")
    cfg["ml"].setdefault("fail_on_error", False)
    cfg["ml"].setdefault("model_dir", "psr_model_ml")
    cfg["ml"].setdefault("input_globs", ["by_protein/*_final_leads.xlsx", "*_clones.csv"])
    cfg["ml"]["model_dir"] = _as_path(cfg["ml"].get("model_dir"), base_dir)
    cfg["ml"].setdefault("delph", {})
    cfg["ml"]["delph"].setdefault("package_path", "/alphafold/combio/software/delphi")
    cfg["ml"]["delph"].setdefault("model_dir", cfg["ml"].get("model_dir"))
    if cfg["ml"]["delph"].get("model_dir"):
        cfg["ml"]["delph"]["model_dir"] = _as_path(cfg["ml"]["delph"]["model_dir"], base_dir)
    if cfg["ml"]["delph"].get("package_path"):
        cfg["ml"]["delph"]["package_path"] = _as_path(cfg["ml"]["delph"]["package_path"], base_dir)
    pythonpath = cfg["ml"]["delph"].get("pythonpath")
    if pythonpath:
        if not isinstance(pythonpath, list):
            pythonpath = [pythonpath]
        cfg["ml"]["delph"]["pythonpath"] = [_as_path(path, base_dir) for path in pythonpath]

    cfg.setdefault("registry", {})
    cfg["registry"].setdefault("enabled", True)
    cfg["registry"].setdefault("path", None)
    cfg["registry"].setdefault("project", "IPI")
    if cfg["registry"].get("path"):
        cfg["registry"]["path"] = _as_path(cfg["registry"]["path"], base_dir)

    cfg.setdefault("stage_b", {})
    cfg["stage_b"].setdefault("target_sequence_table", None)
    if cfg["stage_b"].get("target_sequence_table"):
        cfg["stage_b"]["target_sequence_table"] = _as_path(cfg["stage_b"]["target_sequence_table"], base_dir)

    return cfg
