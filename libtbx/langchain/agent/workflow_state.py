"""
Workflow State Utilities for PHENIX AI Agent.

This module provides:
- File categorization by type and purpose
- History analysis (what programs have been run)
- Experiment type detection (X-ray vs Cryo-EM)
- Workflow state detection (delegating to WorkflowEngine)
- Prompt formatting for LLM

The actual workflow logic is defined in:
- knowledge/workflows.yaml (state machine)
- knowledge/file_categories.yaml (file categorization rules)
- agent/workflow_engine.py (YAML interpreter)
"""

from __future__ import absolute_import, division, print_function
import os
import re
import fnmatch


# =============================================================================
# FILE CATEGORIZATION
# =============================================================================

def _load_category_rules():
    """Load file category rules from YAML."""
    from libtbx.langchain.knowledge.yaml_loader import load_file_categories
    return load_file_categories()


def _match_pattern(filename, pattern):
    """Check if filename matches a pattern (supports * wildcards)."""
    # Convert pattern to regex-friendly fnmatch pattern
    return fnmatch.fnmatch(filename.lower(), pattern.lower())


def _categorize_files(available_files):
    """
    Categorize files by type and purpose.

    Uses rules from knowledge/file_categories.yaml.

    Returns dict with keys for BOTH:
    - Subcategories: refined, phaser_output, predicted, etc.
    - Parent categories: model, search_model, ligand, map, mtz, sequence

    Files in subcategories are automatically "bubbled up" to their parent
    semantic categories. For example:
        - refined -> also in model
        - phaser_output -> also in model
        - predicted -> also in search_model
        - processed_predicted -> also in search_model
    """
    # Try to load YAML rules
    category_rules = _load_category_rules()

    if category_rules:
        files = _categorize_files_yaml(available_files, category_rules)
    else:
        # Fallback to hardcoded rules if YAML not available
        files = _categorize_files_hardcoded(available_files)

    # Bubble up subcategories to their parent semantic categories
    files = _bubble_up_to_parents(files, category_rules)

    # Post-processing: If we have exactly one half-map and no full maps,
    # treat it as a full map. Half-maps only make sense in pairs for FSC.
    # A user providing a single map (even if named like a half-map) wants to use it.
    if "half_map" in files and "full_map" in files:
        if len(files["half_map"]) == 1 and len(files["full_map"]) == 0:
            files["full_map"].append(files["half_map"][0])
            files["half_map"] = []

    return files


# Mapping from subcategory to parent semantic category
# This is the source of truth for bubbling up
SUBCATEGORY_TO_PARENT = {
    # Model subcategories (positioned, ready for refinement)
    "refined": "model",
    "rsr_output": "model",
    "phaser_output": "model",
    "autobuild_output": "model",
    "docked": "model",
    "with_ligand": "model",
    "ligand_fit_output": "ligand",  # Ligand fragment from ligandfit, not a full model
    "model_cif": "model",
    "unclassified_pdb": "model",  # Generic PDB files bubble to model category
                                  # but _has_placed_model checks subcategories to determine placement

    # Search model subcategories (templates, NOT positioned)
    "predicted": "search_model",
    "processed_predicted": "search_model",
    "pdb_template": "search_model",

    # Ligand subcategories
    "ligand_pdb": "ligand",
    "ligand_cif": "ligand",

    # Map subcategories
    "full_map": "map",
    "half_map": "map",
    "optimized_full_map": "map",
    "sharpened": "map",

    # Data MTZ subcategories (measured Fobs, R-free)
    "original_data_mtz": "data_mtz",
    "phased_data_mtz": "data_mtz",

    # Map coefficients MTZ subcategories (calculated phases)
    "refine_map_coeffs": "map_coeffs_mtz",
    "denmod_map_coeffs": "map_coeffs_mtz",
    "predict_build_map_coeffs": "map_coeffs_mtz",

    # Intermediate - these should NOT be bubbled up or tracked
    # Set to "intermediate" parent so they're excluded from model/search_model
    "intermediate_mr": "intermediate",
    "autobuild_temp": "intermediate",
    "carryover_temp": "intermediate",
}


def _bubble_up_to_parents(files, category_rules=None):
    """
    Ensure files in subcategories also appear in their parent semantic categories.

    This enables programs to request by parent category (e.g., "model")
    while files are categorized into specific subcategories (e.g., "refined").

    Args:
        files: Dict of category -> list of files
        category_rules: Optional YAML rules (for dynamic parent lookup)

    Returns:
        Updated files dict with parent categories populated
    """
    # Ensure parent categories exist
    for parent in ["model", "search_model", "ligand", "intermediate", "map", "data_mtz", "map_coeffs_mtz", "sequence"]:
        if parent not in files:
            files[parent] = []

    # Build parent mapping from YAML if available
    parent_map = dict(SUBCATEGORY_TO_PARENT)  # Start with hardcoded

    if category_rules:
        for cat_name, cat_def in category_rules.items():
            parent = cat_def.get("parent_category")
            if parent:
                parent_map[cat_name] = parent

    # Bubble up each subcategory to its parent
    for subcat, parent in parent_map.items():
        if parent is None:
            continue  # Don't bubble up if no parent
        if subcat not in files:
            continue
        if parent not in files:
            files[parent] = []

        for f in files[subcat]:
            if f not in files[parent]:
                files[parent].append(f)

    # Also ensure backward compatibility: "pdb" category contains all model and search_model
    if "pdb" not in files:
        files["pdb"] = []
    for f in files.get("model", []):
        if f not in files["pdb"]:
            files["pdb"].append(f)
    for f in files.get("search_model", []):
        if f not in files["pdb"]:
            files["pdb"].append(f)

    return files


def _categorize_files_yaml(available_files, rules):
    """
    Categorize files using YAML-defined rules.

    This handles two types of categories:
    1. Extension-based primary categories (data_mtz, map_coeffs_mtz, map, sequence)
    2. Pattern-based subcategories with semantic parents (refined->model, predicted->search_model)

    Files are first matched to subcategories by patterns, then bubbled up to parent categories.
    """
    # Initialize all categories from YAML
    files = {cat: [] for cat in rules.keys()}

    # Ensure semantic parent categories exist
    for parent in ["model", "search_model", "ligand", "map", "data_mtz", "map_coeffs_mtz", "sequence", "intermediate"]:
        if parent not in files:
            files[parent] = []

    # Also ensure pdb exists for backward compatibility
    if "pdb" not in files:
        files["pdb"] = []

    # Group categories by extension for primary matching
    # Only include categories that are PURELY extension-based (no patterns)
    # Categories with patterns are handled in Step 2
    ext_to_categories = {}
    ext_to_excludes = {}  # Track excludes for each category
    for cat_name, cat_def in rules.items():
        # Skip semantic parent categories - they don't match by extension
        if cat_def.get("is_semantic_parent"):
            continue
        # Skip categories that have patterns - they need pattern matching in Step 2
        if cat_def.get("patterns"):
            continue
        for ext in cat_def.get("extensions", []):
            if ext not in ext_to_categories:
                ext_to_categories[ext] = []
            ext_to_categories[ext].append(cat_name)
            # Store excludes for this category
            if cat_def.get("excludes"):
                ext_to_excludes[cat_name] = cat_def.get("excludes", [])

    for f in available_files:
        f_lower = f.lower()
        basename = os.path.basename(f_lower)
        _, ext = os.path.splitext(f_lower)

        # Step 1: Primary categorization by extension (for non-PDB files)
        primary_categories = ext_to_categories.get(ext, [])
        for cat in primary_categories:
            # Check excludes before adding
            excludes = ext_to_excludes.get(cat, [])
            excluded = False
            for exc_pattern in excludes:
                if _match_pattern(basename, exc_pattern):
                    excluded = True
                    break
            if not excluded and f not in files[cat]:
                files[cat].append(f)

        # Step 2: Subcategorization by patterns
        # This is where we match PDB files to specific subcategories like "refined", "predicted"
        for cat_name, cat_def in rules.items():
            # Skip semantic parent categories
            if cat_def.get("is_semantic_parent"):
                continue

            # Skip deprecated categories
            if cat_def.get("is_deprecated"):
                continue

            # Check if this category uses subcategory_of (old style) or parent_category (new semantic style)
            old_parent = cat_def.get("subcategory_of")  # e.g., refined -> pdb
            semantic_parent = cat_def.get("parent_category")  # e.g., refined -> model

            # For old-style subcategories, check if file is in parent
            if old_parent and f not in files.get(old_parent, []):
                continue

            # For new-style semantic subcategories with extension requirements
            if semantic_parent and "extensions" in cat_def:
                cat_extensions = cat_def.get("extensions", [])
                if not any(f_lower.endswith(e) for e in cat_extensions):
                    continue

            # Check patterns
            patterns = cat_def.get("patterns", [])
            excludes = cat_def.get("excludes", [])
            max_len = cat_def.get("max_basename_length")

            # Check excludes first
            excluded = False
            for exc_pattern in excludes:
                if _match_pattern(basename, exc_pattern):
                    excluded = True
                    break

            if excluded:
                continue

            # Check max basename length
            if max_len and len(os.path.splitext(basename)[0]) > max_len:
                continue

            # Check patterns
            matched = False
            if not patterns and semantic_parent:
                # No patterns = use extension matching only (already checked above)
                matched = any(f_lower.endswith(e) for e in cat_def.get("extensions", []))
            else:
                for pattern in patterns:
                    if pattern == "*":
                        # Wildcard matches all (used with excludes)
                        matched = True
                        break
                    if _match_pattern(basename, pattern):
                        matched = True
                        break

            if matched and f not in files[cat_name]:
                files[cat_name].append(f)

                # Also add to "also_in" categories
                for also_cat in cat_def.get("also_in", []):
                    if also_cat in files and f not in files[also_cat]:
                        files[also_cat].append(f)

    return files


def _categorize_files_hardcoded(available_files):
    """
    Categorize files by type and purpose.

    Returns dict with keys:
        data_mtz, map_coeffs_mtz, pdb, sequence, map, ligand_cif, ligand_pdb,
        phaser_output, refined, rsr_output, with_ligand, ligand_fit, predicted,
        full_map, half_map (for cryo-EM half-maps)
    """
    files = {
        "data_mtz": [],  # Reflection data with Fobs, R-free (for refinement)
        "map_coeffs_mtz": [],  # Map coefficients with phases (for ligand fitting)
        "pdb": [],
        "sequence": [],
        "map": [],  # All map files (for backward compatibility)
        "full_map": [],  # Full cryo-EM maps
        "half_map": [],  # Half maps (usually come in pairs)
        "ligand_cif": [],
        "ligand_pdb": [],
        "phaser_output": [],
        "refined": [],  # X-ray refinement output
        "rsr_output": [],  # real_space_refine output (cryo-EM)
        "with_ligand": [],
        "ligand_fit_output": [],
        "predicted": [],
        "processed_predicted": [],
        "autobuild_output": [],
        "docked": [],  # dock_in_map output
        "intermediate_mr": [],  # Intermediate MR files (never use for refinement)
    }

    # X-ray data file extensions (PHENIX can read these)
    xray_data_extensions = ('.mtz', '.sca', '.hkl', '.sdf')

    def is_half_map(basename):
        """Detect if a map file is a half-map based on naming conventions."""
        name_lower = basename.lower()
        name_no_ext = os.path.splitext(name_lower)[0]

        if 'half' in name_lower:
            return True
        if re.search(r'[_\-]?[12]$', name_no_ext):
            return True
        if re.search(r'[_\-][ab]$', name_no_ext):
            return True

        return False

    # Import shared MTZ classification
    try:
        from libtbx.langchain.agent.file_utils import classify_mtz_type
    except ImportError:
        from agent.file_utils import classify_mtz_type

    for f in available_files:
        f_lower = f.lower()
        basename = os.path.basename(f_lower)

        # Primary type categorization
        if f_lower.endswith(xray_data_extensions):
            # Classify into data_mtz or map_coeffs_mtz
            mtz_type = classify_mtz_type(f)  # Pass full path
            files[mtz_type].append(f)
        elif f_lower.endswith('.pdb'):
            files["pdb"].append(f)

            # Subcategorize PDBs by their origin/purpose
            if 'phaser' in basename or basename.startswith('phaser'):
                files["phaser_output"].append(f)

            if 'refine' in basename and 'real_space' not in basename and 'rsr' not in basename:
                files["refined"].append(f)

            # RSR output detection - real_space_refine outputs contain 'real_space_refined'
            # e.g., model_real_space_refined_000.pdb
            if 'real_space_refined' in basename or 'rsr_' in basename or '_rsr' in basename:
                files["rsr_output"].append(f)

            if 'with_ligand' in basename:
                files["with_ligand"].append(f)
            if 'ligand_fit' in basename or 'ligandfit' in basename:
                files["ligand_fit_output"].append(f)

            if 'predict' in basename or 'alphafold' in basename or 'colabfold' in basename:
                files["predicted"].append(f)

            if 'processed' in basename:
                files["processed_predicted"].append(f)

            is_autobuild = (
                ('autobuild' in basename or 'auto_build' in basename) or
                ('AutoBuild' in os.path.basename(f)) or
                ('overall_best' in basename) or
                ('build' in basename or 'built' in basename) or
                'buccaneer' in basename or
                'arp_warp' in basename or
                ('shelxe' in basename and 'trace' in basename)
            )
            if is_autobuild and 'predict' not in basename:
                files["autobuild_output"].append(f)

            if 'dock' in basename and 'map' in basename:
                files["docked"].append(f)
            # Also match placed_model* from dock_in_map output
            if basename.startswith('placed_model') or '_placed' in basename:
                files["docked"].append(f)

            # Intermediate MR files from dock_in_map - never use for refinement
            if basename.startswith('run_mr') or fnmatch.fnmatch(basename, '*mr.[0-9]*'):
                files["intermediate_mr"].append(f)

            if (basename.startswith('lig') and len(basename) < 20) or 'ligand' in basename:
                if not any(x in basename for x in ['ligand_fit', 'ligandfit', 'with_ligand']):
                    files["ligand_pdb"].append(f)

        elif f_lower.endswith(('.fa', '.fasta', '.seq', '.dat')):
            files["sequence"].append(f)
        elif f_lower.endswith(('.mrc', '.ccp4', '.map')):
            files["map"].append(f)
            if is_half_map(basename):
                files["half_map"].append(f)
            else:
                files["full_map"].append(f)
        elif f_lower.endswith('.cif'):
            if 'refine' in basename:
                files["pdb"].append(f)
                files["refined"].append(f)
            else:
                files["ligand_cif"].append(f)

    return files


def _detect_experiment_type(files, history_info=None):
    """
    Determine if this is X-ray crystallography or Cryo-EM.

    Logic:
    - If mtriage has been run → cryo-EM (definitive)
    - Has map (full or half) but no MTZ → cryo-EM
    - Has MTZ → X-ray (even if map also present)
    - Neither → unknown (default to X-ray)
    """
    if history_info and history_info.get("mtriage_done"):
        return "cryoem"

    if history_info and history_info.get("rsr_done"):
        return "cryoem"

    has_data_mtz = bool(files.get("data_mtz")) or bool(files.get("map_coeffs_mtz"))
    has_map = bool(files.get("map")) or bool(files.get("full_map")) or bool(files.get("half_map"))

    if has_map and not has_data_mtz:
        return "cryoem"
    else:
        return "xray"


# =============================================================================
# HISTORY ANALYSIS
# =============================================================================

def _is_failed_result(result):
    """
    Check if a result string indicates failure.

    Uses specific patterns to avoid false positives like "No ERROR detected".

    Args:
        result: Result string from history entry

    Returns:
        bool: True if result indicates failure
    """
    if not result:
        return False

    result_upper = result.upper()

    # Specific failure patterns (avoid matching "No ERROR detected" etc.)
    failure_patterns = [
        'FAILED',           # Common failure indicator
        'SORRY:',           # Phenix error prefix
        'SORRY ',           # Phenix error prefix with space
        'ERROR:',           # Error with colon
        'ERROR ',           # Error as prefix
        ': ERROR',          # Error after colon
        'TRACEBACK',        # Python exception
        'EXCEPTION',        # Exception indicator
    ]

    return any(pattern in result_upper for pattern in failure_patterns)


# Cache for done tracking config loaded from YAML
_DONE_TRACKING_CACHE = None

# Allowed count_field values — rejects typos at load time
ALLOWED_COUNT_FIELDS = {"refine_count", "rsr_count", "phaser_count"}


def _load_done_tracking_configs():
    """Load done_tracking configuration from programs.yaml.

    Returns list of dicts for all programs with history_detection, each with:
      - flag: done flag name (e.g., 'dock_done')
      - strategy: 'set_flag' (default), 'run_once', or 'count'
      - count_field: counter name for strategy='count' (e.g., 'refine_count')
      - markers: list of strings to match in combined text (OR logic,
                 substring matching)
      - exclude_markers: list of strings that reject a match (checked FIRST)
      - alt_markers: optional list of alternate marker strings
      - alt_requires: optional list of strings that must ALL be present
                      alongside alt_markers (AND logic)
      - success_flag: optional additional flag set on success

    Validates count_field against ALLOWED_COUNT_FIELDS at load time to
    prevent typos from silently creating garbage attributes.
    """
    global _DONE_TRACKING_CACHE
    if _DONE_TRACKING_CACHE is not None:
        return _DONE_TRACKING_CACHE

    configs = []
    try:
        try:
            from libtbx.langchain.knowledge.yaml_loader import load_programs
        except ImportError:
            from knowledge.yaml_loader import load_programs
        programs = load_programs()
        for name, defn in programs.items():
            if not isinstance(defn, dict):
                continue
            tracking = defn.get("done_tracking", {})
            detection = tracking.get("history_detection")
            if not detection or not isinstance(detection, dict):
                continue

            strategy = tracking.get("strategy", "set_flag")
            count_field = tracking.get("count_field")

            # Validate count_field at load time
            if strategy == "count":
                if not count_field:
                    print("  [workflow_state] Warning: %s has strategy='count' "
                          "but no count_field" % name)
                    continue
                if count_field not in ALLOWED_COUNT_FIELDS:
                    raise ValueError(
                        "Unknown count_field %r in done_tracking for %s. "
                        "Allowed: %s" % (count_field, name, ALLOWED_COUNT_FIELDS))

            configs.append({
                "program": name,
                "flag": tracking.get("flag"),
                "strategy": strategy,
                "count_field": count_field,
                "markers": detection.get("markers", []),
                "exclude_markers": detection.get("exclude_markers", []),
                "alt_markers": detection.get("alt_markers", []),
                "alt_requires": detection.get("alt_requires", []),
                "success_flag": detection.get("success_flag"),
            })
    except Exception as e:
        print("  [workflow_state] Warning: could not load done_tracking "
              "from programs.yaml: %s" % str(e))

    _DONE_TRACKING_CACHE = configs
    return configs


def _set_done_flags(info, combined, result):
    """Set done flags from YAML history_detection for all strategies.

    Handles: set_flag, run_once, and count strategies.
    The run_once filtering itself happens in program_registration;
    here we just set the flag and increment counts.

    Args:
        info: The info dict being built by _analyze_history
        combined: Lowercase string of program + command
        result: Result string from history entry
    """
    if _is_failed_result(result):
        return  # All strategies require success

    configs = _load_done_tracking_configs()

    for config in configs:
        flag = config["flag"]
        if not flag:
            continue

        # Exclude markers take precedence — checked FIRST
        if any(m in combined for m in config["exclude_markers"]):
            continue

        # Check primary markers (OR logic, substring matching)
        matched = any(m in combined for m in config["markers"])

        # Check alt_markers with alt_requires (AND logic)
        if not matched and config["alt_markers"] and config["alt_requires"]:
            if (any(m in combined for m in config["alt_markers"]) and
                    all(r in combined for r in config["alt_requires"])):
                matched = True

        if matched:
            info[flag] = True

            # Count strategy: increment count field
            if config["strategy"] == "count" and config["count_field"]:
                info[config["count_field"]] = info.get(config["count_field"], 0) + 1

            # Optional success flag
            if config["success_flag"]:
                info[config["success_flag"]] = True


def _analyze_history(history):
    """
    Extract information about what has been done from history.

    Returns dict with program completion flags and counts.

    Note: Simple done flags for programs with run_once: true are auto-generated
    from programs.yaml via program_registration. Complex flags (counts, success
    conditions) are still handled manually below.
    """
    # =========================================================================
    # Initialize flags from YAML done_tracking configs
    # =========================================================================
    info = {
        "programs_run": set(),
        # predict_and_build flags — no history_detection (Python-only cascade)
        "predict_done": False,
        "predict_full_done": False,
        # Post-ligandfit refinement tracking
        "needs_post_ligandfit_refine": False,
        # Metrics extracted from history
        "last_program": None,
        "last_r_free": None,
        "last_map_cc": None,
        "last_clashscore": None,
        "last_tfz": None,
        "resolution": None,
        "anomalous_resolution": None,
        "anomalous_measurability": None,
        "has_anomalous": False,
        "strong_anomalous": False,
        "has_twinning": False,
        "twin_law": None,
        "twin_fraction": None,
        "has_ncs": False,  # NCS detected in data
    }

    # Initialize done flags and count fields from YAML (single source of truth)
    configs = _load_done_tracking_configs()
    for config in configs:
        if config["flag"]:
            info[config["flag"]] = False
        if config["strategy"] == "count" and config["count_field"]:
            info[config["count_field"]] = 0
        if config.get("success_flag"):
            info[config["success_flag"]] = False

    if not history:
        return info

    # =========================================================================
    # Process history for flags and metrics
    # =========================================================================
    for entry in history:
        prog = ""
        cmd = ""

        if isinstance(entry, str):
            prog = entry.lower()
            cmd = entry.lower()
        elif isinstance(entry, dict):
            prog = (entry.get("program") or "").lower()
            cmd = entry.get("command", "").lower()

            # Extract metrics
            # NOTE: history from session has 'analysis' key, but after transport has 'metrics'
            analysis = entry.get("analysis", entry.get("metrics", {}))
            if isinstance(analysis, dict):
                if analysis.get("r_free"):
                    info["last_r_free"] = analysis["r_free"]
                if analysis.get("map_cc"):
                    info["last_map_cc"] = analysis["map_cc"]
                if analysis.get("clashscore"):
                    info["last_clashscore"] = analysis["clashscore"]
                if analysis.get("tfz"):
                    info["last_tfz"] = analysis["tfz"]
                if analysis.get("resolution"):
                    info["resolution"] = analysis["resolution"]
                if analysis.get("anomalous_resolution"):
                    info["anomalous_resolution"] = analysis["anomalous_resolution"]
                    info["has_anomalous"] = True
                elif analysis.get("has_anomalous"):
                    info["has_anomalous"] = analysis["has_anomalous"]
                # Store anomalous measurability for decision making
                if analysis.get("anomalous_measurability"):
                    info["anomalous_measurability"] = analysis["anomalous_measurability"]
                    # Strong anomalous signal if measurability > 0.10
                    if analysis["anomalous_measurability"] > 0.10:
                        info["has_anomalous"] = True
                        info["strong_anomalous"] = True
                # Twinning threshold (0.20) from workflows.yaml shared section
                if analysis.get("twin_law") and analysis.get("twin_fraction"):
                    twin_frac = analysis["twin_fraction"]
                    if twin_frac > 0.20 and not analysis.get("no_twinning_suspected"):
                        info["has_twinning"] = True
                        info["twin_law"] = analysis["twin_law"]
                        info["twin_fraction"] = twin_frac
                # NCS detection (from map_symmetry or similar)
                if analysis.get("ncs_found") or analysis.get("has_ncs"):
                    info["has_ncs"] = True
        else:
            continue

        combined = prog + " " + cmd

        info["programs_run"].add(prog)
        info["last_program"] = prog

        # =====================================================================
        # Program done flags — YAML-driven detection for all strategies
        # =====================================================================
        result = entry.get("result", "") if isinstance(entry, dict) else ""

        # Handles: set_flag, run_once, and count strategies.
        # Covers all programs with history_detection in programs.yaml.
        _set_done_flags(info, combined, result)

        # The ONE remaining Python-only case: predict_and_build cascade.
        # When predict runs fully, it also sets refine_done + refine_count.
        if "predict_and_build" in combined:
            if not _is_failed_result(result):
                info["predict_done"] = True
                if "stop_after_predict=true" not in combined and "stop_after_predict=True" not in combined:
                    info["predict_full_done"] = True
                    info["refine_done"] = True
                    info["refine_count"] += 1

        # Track whether refinement is needed after ligandfit.
        # After ligandfit adds a ligand, the complex always needs re-refinement.
        # This flag is True when ligandfit succeeded but no refine happened after.
        if not _is_failed_result(result):
            if "ligandfit" in combined:
                info["needs_post_ligandfit_refine"] = True
            elif "refine" in combined and "real_space" not in combined:
                info["needs_post_ligandfit_refine"] = False

    return info


# =============================================================================
# WORKFLOW STATE DETECTION
# =============================================================================


def detect_workflow_state(history, available_files, analysis=None, maximum_automation=True,
                         use_yaml_engine=True, directives=None):
    """
    Determine current workflow state based on history and files.

    This function delegates to the YAML-driven WorkflowEngine for state detection.

    Args:
        history: List of cycle records from client
        available_files: List of available file paths
        analysis: Current log analysis dict (optional)
        maximum_automation: If True, use fully automated cryo-EM path
        use_yaml_engine: If True, use YAML-driven WorkflowEngine (default: True)
        directives: Optional user directives dict

    Returns:
        dict: {
            state: str,              # State name
            experiment_type: str,    # "xray" or "cryoem"
            valid_programs: list,    # Programs allowed in this state
            reason: str,             # Human-readable explanation
            conditions: dict,        # Conditional program availability
            automation_path: str,    # "stepwise" or "automated" (cryo-EM only)
            categorized_files: dict, # Pre-categorized files (full_map, half_map, etc.)
        }
    """
    # Categorize files
    files = _categorize_files(available_files)

    # Analyze history
    history_info = _analyze_history(history)

    # Determine experiment type
    experiment_type = _detect_experiment_type(files, history_info)

    # Use YAML-driven workflow engine
    if use_yaml_engine:
        try:
            # Lazy import to avoid circular dependencies
            from libtbx.langchain.agent.workflow_engine import WorkflowEngine

            engine = WorkflowEngine()
            state = engine.get_workflow_state(experiment_type, files, history_info, analysis,
                                             directives, maximum_automation)

            state["categorized_files"] = files
            # Set automation_path for both experiment types
            state["automation_path"] = "automated" if maximum_automation else "stepwise"
            return state
        except Exception as e:
            import sys
            print("Warning: YAML workflow engine failed: %s" % e, file=sys.stderr)

    # Fallback: return minimal state (should not happen if YAML is properly configured)
    return {
        "state": "unknown",
        "experiment_type": experiment_type,
        "valid_programs": ["STOP"],
        "reason": "Workflow engine unavailable",
        "conditions": {},
        "automation_path": "automated" if maximum_automation else "stepwise",
        "categorized_files": files,
    }


# =============================================================================
# VALIDATION
# =============================================================================

def validate_program_choice(chosen_program, workflow_state):
    """
    Validate that a program choice is allowed in the current state.

    Args:
        chosen_program: Program the LLM chose
        workflow_state: Dict from detect_workflow_state()

    Returns:
        tuple: (is_valid: bool, error_message: str or None)
    """
    if chosen_program is None:
        return True, None

    if chosen_program == "STOP":
        return True, None

    valid = workflow_state["valid_programs"]

    if chosen_program in valid:
        return True, None

    try:
        from libtbx.langchain.knowledge.yaml_loader import get_all_programs
        all_known_programs = get_all_programs()
    except Exception:
        all_known_programs = []  # Graceful degradation

    if chosen_program in all_known_programs:
        error = (
            "Program '%s' is not valid in state '%s'. "
            "Valid programs: %s. Reason: %s"
        ) % (
            chosen_program,
            workflow_state["state"],
            ", ".join(valid),
            workflow_state["reason"]
        )
    else:
        error = "Unknown program '%s'. Valid programs: %s" % (chosen_program, ", ".join(valid))

    return False, error


# =============================================================================
# PROMPT FORMATTING
# =============================================================================

def format_workflow_for_prompt(workflow_state):
    """
    Format workflow state for inclusion in LLM prompt.

    Args:
        workflow_state: Output from detect_workflow_state()

    Returns:
        str: Formatted text for prompt
    """
    lines = []

    lines.append("### WORKFLOW STATE: %s" % workflow_state["state"])
    lines.append("Experiment type: %s" % workflow_state["experiment_type"])

    if workflow_state.get("automation_path"):
        lines.append("Automation path: %s" % workflow_state["automation_path"])

    lines.append("")
    lines.append(workflow_state["reason"])
    lines.append("")
    lines.append("**VALID PROGRAMS FOR THIS STATE:**")
    lines.append(", ".join(workflow_state["valid_programs"]))
    lines.append("")
    lines.append("⚠️ You MUST choose a program from the list above, or set \"stop\": true.")
    lines.append("Choosing an invalid program will cause a validation error.")

    if workflow_state.get("conditions"):
        lines.append("")
        lines.append("Conditional availability:")
        for prog, condition in workflow_state["conditions"].items():
            lines.append("  - %s: requires %s" % (prog, condition))

    # Show stepwise mode hint for both cryo-EM and X-ray
    if workflow_state.get("automation_path") == "stepwise":
        stepwise_states = ["cryoem_analyzed", "xray_initial", "xray_placed"]
        if workflow_state["state"] in stepwise_states:
            lines.append("")
            lines.append("NOTE (Stepwise mode): predict_and_build will use stop_after_predict=true")

    return "\n".join(lines)
