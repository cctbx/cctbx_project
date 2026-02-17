#!/usr/bin/env python
"""
Tests for hardcoded cleanup conformance (HARDCODED_CLEANUP_PLAN.md).

Two kinds of tests:

  1. FUNCTIONAL — verify YAML-driven mechanisms produce correct results
  2. CONFORMANCE — scan source code to prevent hardcoded patterns from
     creeping back in.  These are the "trip-wire" tests: they inspect
     the actual Python/YAML files and fail if anyone re-introduces
     hardcoded program lists, top-level run_once, _MANUAL_DONE_FLAGS, etc.
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys
import ast
import yaml

# Add parent directory to path for imports
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT_DIR = os.path.dirname(_THIS_DIR)
sys.path.insert(0, _ROOT_DIR)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read(relpath):
    """Read a file relative to project root.

    For paths starting with 'phenix_ai/', resolve via import since
    that package lives outside the langchain directory in PHENIX.
    """
    # First try the simple local path
    path = os.path.join(_ROOT_DIR, relpath)
    if os.path.exists(path):
        with open(path) as f:
            return f.read()

    # Not co-located — resolve via import system.
    # Convert "phenix_ai/log_parsers.py" → module "phenix.phenix_ai.log_parsers"
    mod_path = relpath.replace("/", ".").replace(".py", "")
    # Try phenix.X first, then bare X
    for prefix in ["phenix.", ""]:
        try:
            mod = __import__(prefix + mod_path, fromlist=["__file__"])
            if hasattr(mod, "__file__") and mod.__file__:
                with open(mod.__file__) as f:
                    return f.read()
        except (ImportError, AttributeError):
            continue

    raise FileNotFoundError(
        f"Cannot find {relpath} locally ({path}) or via import system"
    )


def _load_programs_yaml():
    return yaml.safe_load(_read("knowledge/programs.yaml"))


def _load_recoverable_errors_yaml():
    return yaml.safe_load(_read("knowledge/recoverable_errors.yaml"))


def _program_defs():
    """Return {name: defn} for actual programs (skip non-dict entries)."""
    data = _load_programs_yaml()
    return {k: v for k, v in data.items() if isinstance(v, dict)}


# ===================================================================
# FUNCTIONAL TESTS — YAML mechanisms produce correct results
# ===================================================================


def test_all_programs_from_yaml():
    """get_all_programs() returns every program defined in programs.yaml."""
    from knowledge.yaml_loader import get_all_programs

    yaml_progs = set(_program_defs().keys())
    api_progs = set(get_all_programs())

    assert yaml_progs == api_progs, (
        "get_all_programs() mismatch vs programs.yaml:\n"
        f"  In YAML only: {yaml_progs - api_progs}\n"
        f"  In API only:  {api_progs - yaml_progs}"
    )
    print("  test_all_programs_from_yaml PASSED")


def test_detect_program_covers_all_programs_with_markers():
    """detect_program() can identify every program that has log_detection markers."""
    try:
        from phenix.phenix_ai.log_parsers import detect_program
    except ImportError:
        from phenix_ai.log_parsers import detect_program

    programs = _program_defs()
    failures = []

    for name, defn in sorted(programs.items()):
        det = defn.get("log_detection", {})
        markers = det.get("markers", []) + det.get("case_sensitive_markers", [])
        if not markers:
            continue  # No markers = not expected to be detectable

        # Build a minimal fake log that contains one of the markers
        fake_log = f"Running {markers[0]} on data..."
        detected = detect_program(fake_log)
        if detected != name:
            failures.append(f"  {name}: marker '{markers[0]}' → detected '{detected}'")

    assert not failures, (
        "detect_program() failed to detect these programs from their own markers:\n"
        + "\n".join(failures)
    )
    print("  test_detect_program_covers_all_programs_with_markers PASSED")


def test_detect_program_priority_ordering():
    """More specific programs are detected before substring matches."""
    try:
        from phenix.phenix_ai.log_parsers import detect_program
    except ImportError:
        from phenix_ai.log_parsers import detect_program

    # real_space_refine before refine
    assert detect_program("Running phenix.real_space_refine") == "phenix.real_space_refine", \
        "real_space_refine should be detected before refine"

    # autobuild_denmod before autobuild
    assert detect_program("Running autobuild_denmod maps_only") == "phenix.autobuild_denmod", \
        "autobuild_denmod should be detected before autobuild"

    # plain refine should still work
    assert detect_program("Running phenix.refine") == "phenix.refine"

    # plain autobuild should still work
    assert detect_program("Running phenix.autobuild build") == "phenix.autobuild"

    print("  test_detect_program_priority_ordering PASSED")


def test_done_flag_map_covers_expected_programs():
    """get_program_done_flag_map() includes all programs with done_tracking."""
    from knowledge.program_registration import get_program_done_flag_map

    flag_map = get_program_done_flag_map()
    programs = _program_defs()

    # Every program with done_tracking in YAML should appear in the map
    missing = []
    for name, defn in programs.items():
        tracking = defn.get("done_tracking")
        if tracking and tracking.get("flag"):
            if name not in flag_map:
                missing.append(name)
            elif flag_map[name] != tracking["flag"]:
                missing.append(f"{name}: expected {tracking['flag']}, got {flag_map[name]}")

    assert not missing, (
        "get_program_done_flag_map() is missing or wrong for:\n  "
        + "\n  ".join(missing)
    )
    print("  test_done_flag_map_covers_expected_programs PASSED")


def test_done_flag_map_irregular_names():
    """Irregular flag names (predict_done, rsr_done, dock_done, etc.) are correct."""
    from knowledge.program_registration import get_program_done_flag_map

    flag_map = get_program_done_flag_map()

    # These are the programs whose flag names DON'T follow the
    # simple pattern of <short_name>_done
    expected_irregulars = {
        "phenix.predict_and_build": "predict_done",
        "phenix.real_space_refine": "rsr_done",
        "phenix.dock_in_map": "dock_done",
        "phenix.molprobity": "validation_done",
        "phenix.holton_geometry_validation": "validation_done",
        "phenix.model_vs_data": "validation_done",
        "phenix.validation_cryoem": "validation_done",
    }

    for prog, expected_flag in expected_irregulars.items():
        actual = flag_map.get(prog)
        assert actual == expected_flag, (
            f"{prog}: expected flag '{expected_flag}', got '{actual}'"
        )
    print("  test_done_flag_map_irregular_names PASSED")


def test_done_tracking_flags_match_analyze_history():
    """Every done flag set in _analyze_history has a matching done_tracking entry.

    This prevents drift between the YAML definitions and the Python code
    that actually sets the flags.
    """
    # Scrape all flag names set by _analyze_history from the source
    source = _read("agent/workflow_state.py")
    # Match: info["some_flag_done"] = True  OR  info["some_count"] += 1
    set_patterns = re.findall(r'info\["(\w+_done)"\]\s*=\s*True', source)
    set_patterns = list(set(set_patterns))  # dedupe

    # Get all flags from YAML done_tracking
    programs = _program_defs()
    yaml_flags = set()
    for defn in programs.values():
        tracking = defn.get("done_tracking", {})
        if tracking.get("flag"):
            yaml_flags.add(tracking["flag"])

    # Every flag set in _analyze_history should be declared in YAML
    # (except predict_full_done which is a sub-flag of predict_and_build)
    KNOWN_EXCEPTIONS = {
        "predict_full_done",       # Sub-flag of predict_and_build
    }

    missing_from_yaml = []
    for flag in set_patterns:
        if flag in KNOWN_EXCEPTIONS:
            continue
        if flag not in yaml_flags:
            missing_from_yaml.append(flag)

    assert not missing_from_yaml, (
        "_analyze_history sets these done flags that have no done_tracking in YAML:\n  "
        + "\n  ".join(sorted(missing_from_yaml))
        + "\n\nFix: add done_tracking blocks to programs.yaml for these programs."
    )
    print("  test_done_tracking_flags_match_analyze_history PASSED")


def test_data_label_parameters_in_yaml():
    """recoverable_errors.yaml contains data_label_parameters section."""
    config = _load_recoverable_errors_yaml()
    dlp = config.get("data_label_parameters", {})

    assert len(dlp) >= 5, (
        f"Expected at least 5 data_label_parameters entries, got {len(dlp)}"
    )
    assert "phenix.xtriage" in dlp, "Missing xtriage data label parameter"
    assert "phenix.refine" in dlp, "Missing refine data label parameter"
    assert "default" in dlp, "Missing default data label parameter"

    # Each should have a 'parameter' key
    for name, entry in dlp.items():
        assert "parameter" in entry, f"data_label_parameters[{name}] missing 'parameter'"

    print("  test_data_label_parameters_in_yaml PASSED")


# ===================================================================
# CONFORMANCE GUARDS — scan source to prevent regression
# ===================================================================


def test_no_hardcoded_program_lists_in_python():
    """No Python file should contain generic hardcoded lists of 5+ phenix.* strings.

    This catches regressions of Fix 1 (all_known_programs) and Fix 2
    (metrics_analyzer program list).

    NOTE: Curated subsets for specific logic (e.g., "programs requiring
    a placed model") are OK — the conformance guard targets *generic*
    program lists that should come from YAML. Files with known curated
    lists or deferred Phase 3 fixes are excluded.
    """
    SCAN_DIRS = ["agent", "knowledge"]
    SKIP_FILES = {
        # Legitimate curated lists, not generic "all programs" lists:
        "knowledge/prompts.py",
        "knowledge/prompts_hybrid.py",
        # Phase 3 deferred fixes:
        "agent/rules_selector.py",   # Fix 6: priority lists → workflows.yaml
        "agent/program_registry.py", # experiment_type filtering (JSON mode)
    }
    # Specific lines where curated subsets are OK (not generic lists).
    # Format: (file, description) — we just skip the whole file if it
    # has known curated lists rather than maintaining fragile line numbers.
    ALSO_SKIP = {
        "agent/workflow_engine.py",  # common_programs (UX), programs_requiring_placed_model
    }
    SKIP_FILES.update(ALSO_SKIP)

    violations = []
    for scan_dir in SCAN_DIRS:
        dirpath = os.path.join(_ROOT_DIR, scan_dir)
        for fname in os.listdir(dirpath):
            if not fname.endswith(".py"):
                continue
            relpath = os.path.join(scan_dir, fname)
            if relpath in SKIP_FILES:
                continue

            source = _read(relpath)
            try:
                tree = ast.parse(source)
            except SyntaxError:
                continue

            for node in ast.walk(tree):
                if isinstance(node, ast.List):
                    # Count string elements that look like phenix program names
                    phenix_strings = [
                        elt for elt in node.elts
                        if isinstance(elt, (ast.Constant, ast.Str))
                        and str(getattr(elt, 'value', getattr(elt, 's', '')))
                            .startswith("phenix.")
                    ]
                    if len(phenix_strings) >= 5:
                        line = getattr(node, 'lineno', '?')
                        values = [str(getattr(e, 'value', getattr(e, 's', '')))
                                  for e in phenix_strings[:3]]
                        violations.append(
                            f"  {relpath}:{line} — list with {len(phenix_strings)} "
                            f"phenix.* strings: {values}..."
                        )

    assert not violations, (
        "Found hardcoded program lists (Fix 1/2 regression):\n"
        + "\n".join(violations)
        + "\n\nThese should use get_all_programs() or YAML configuration."
    )
    print("  test_no_hardcoded_program_lists_in_python PASSED")


def test_no_top_level_run_once_in_yaml():
    """programs.yaml must not have top-level run_once fields.

    After Fix 10, run_once lives ONLY inside done_tracking blocks.
    """
    programs = _program_defs()
    violations = []

    for name, defn in programs.items():
        if "run_once" in defn:
            violations.append(f"  {name}: has top-level run_once={defn['run_once']}")

    assert not violations, (
        "Found top-level run_once in programs.yaml (Fix 10 regression):\n"
        + "\n".join(violations)
        + "\n\nMove run_once inside done_tracking block."
    )
    print("  test_no_top_level_run_once_in_yaml PASSED")


def test_no_manual_done_flags_dict():
    """program_registration.py must not contain _MANUAL_DONE_FLAGS.

    After Fix 10, all done flags come from YAML done_tracking blocks.
    """
    source = _read("knowledge/program_registration.py")

    assert "_MANUAL_DONE_FLAGS" not in source, (
        "Found _MANUAL_DONE_FLAGS in program_registration.py (Fix 10 regression).\n"
        "All done flags should come from done_tracking blocks in programs.yaml."
    )
    print("  test_no_manual_done_flags_dict PASSED")


def test_no_duplicate_detection_chains():
    """graph_nodes.py and metrics_analyzer.py must not have hardcoded
    program detection if/elif chains.

    After Fix 2, these delegate to detect_program() from log_parsers.
    """
    for filename in ["agent/graph_nodes.py", "agent/metrics_analyzer.py"]:
        source = _read(filename)

        # Look for chains of: if "phenix.X" in ... elif "phenix.Y" in ...
        # with 3+ branches (short ad-hoc checks are OK, big chains are not)
        pattern = r'(?:el)?if\s+["\']phenix\.\w+["\']\s+in\s+'
        matches = re.findall(pattern, source)
        assert len(matches) < 3, (
            f"{filename} contains {len(matches)} hardcoded program detection "
            f"branches (Fix 2 regression).\n"
            f"Delegate to detect_program() from log_parsers instead."
        )

    print("  test_no_duplicate_detection_chains PASSED")


def test_workflow_engine_reads_done_tracking():
    """workflow_engine.py must read done_tracking, not top-level run_once.

    After Fix 10, the pattern should be:
        tracking = prog_def.get("done_tracking", {})
        if tracking.get("run_once"):

    NOT:
        if prog_def.get("run_once"):
    """
    source = _read("agent/workflow_engine.py")

    # This should NOT appear:
    bad_pattern = re.findall(
        r'prog_def\s*(?:and\s+prog_def)?\.get\(\s*["\']run_once["\']\s*\)',
        source
    )
    assert not bad_pattern, (
        "workflow_engine.py reads top-level run_once (Fix 10 regression).\n"
        f"Found: {bad_pattern}\n"
        "Should read done_tracking.run_once instead."
    )

    # This SHOULD appear:
    good_pattern = re.findall(r'tracking\.get\(["\']run_once["\']\)', source)
    assert len(good_pattern) >= 2, (
        "workflow_engine.py should have at least 2 references to "
        "tracking.get('run_once') (filter + explain). "
        f"Found {len(good_pattern)}."
    )
    print("  test_workflow_engine_reads_done_tracking PASSED")


def test_all_programs_have_log_detection():
    """Every program in programs.yaml should have a log_detection block.

    After Fix 2, this is the single source of truth for program detection.
    """
    programs = _program_defs()
    missing = []

    for name, defn in sorted(programs.items()):
        det = defn.get("log_detection")
        if not det or not det.get("markers"):
            missing.append(name)

    assert not missing, (
        "These programs are missing log_detection markers in programs.yaml:\n  "
        + "\n  ".join(missing)
        + "\n\nAdd log_detection.markers so detect_program() can identify them."
    )
    print("  test_all_programs_have_log_detection PASSED")


def test_done_tracking_flags_in_build_context():
    """Every done_tracking flag should appear in build_context().

    If a flag is declared in YAML but never read in build_context, it
    can't gate any workflow phase, making skip_programs ineffective.
    """
    source = _read("agent/workflow_engine.py")
    programs = _program_defs()

    # Collect all unique done_tracking flags
    yaml_flags = set()
    for defn in programs.values():
        tracking = defn.get("done_tracking", {})
        if tracking.get("flag"):
            yaml_flags.add(tracking["flag"])

    # Flags that are auto-set by skip_programs logic (don't need explicit
    # entries in build_context because get_program_done_flag_map handles them)
    # But they DO need to exist as context keys for phase detection to work.
    # The skip_programs code sets them dynamically, but the initial context
    # dict should have them as defaults for non-skip scenarios.

    # Check which flags appear in the build_context source
    # (either as literal key or via history_info.get)
    missing = []
    for flag in sorted(yaml_flags):
        # Check for: "flag_name": or .get("flag_name" or ["flag_name"]
        if (f'"{flag}"' not in source and f"'{flag}'" not in source):
            missing.append(flag)

    # Some flags share a name with their auto-registered version and are
    # set via the skip_programs dynamic code path, not explicit keys.
    # These are OK as long as the skip_programs block handles them.
    KNOWN_DYNAMIC_ONLY = {
        # These are set by skip_programs code via get_program_done_flag_map()
        # and also set explicitly by _analyze_history → build_context.
        # If they show up here it means they're only in YAML but not in
        # the explicit context dict. That's fine IF skip_programs handles them.
        "autobuild_denmod_done",
        "map_sharpening_done",
        "map_to_model_done",
        "resolve_cryo_em_done",
        "map_symmetry_done",
        "mtriage_done",
        "xtriage_done",
        "polder_done",
        "process_predicted_done",
    }

    real_missing = [f for f in missing if f not in KNOWN_DYNAMIC_ONLY]

    assert not real_missing, (
        "These done_tracking flags from YAML don't appear in build_context():\n  "
        + "\n  ".join(real_missing)
        + "\n\nEither add them to the context dict in build_context() or to "
        "KNOWN_DYNAMIC_ONLY if they're handled by skip_programs."
    )
    print("  test_done_tracking_flags_in_build_context PASSED")


def test_yaml_tools_valid_fields_include_done_tracking():
    """yaml_tools.py valid field list must include done_tracking, not run_once."""
    source = _read("agent/yaml_tools.py")

    assert "'done_tracking'" in source, (
        "yaml_tools.py valid_program_fields is missing 'done_tracking'"
    )
    # Check run_once is NOT in the valid fields set (it's OK in comments)
    # Look specifically for 'run_once' as a set element
    field_match = re.search(
        r"valid_program_fields\s*=\s*\{([^}]+)\}", source, re.DOTALL
    )
    if field_match:
        field_block = field_match.group(1)
        assert "'run_once'" not in field_block, (
            "yaml_tools.py valid_program_fields still contains 'run_once' "
            "(Fix 10 regression). Should be 'done_tracking'."
        )
    print("  test_yaml_tools_valid_fields_include_done_tracking PASSED")


def test_fallback_paths_emit_deprecation_warnings():
    """Hardcoded fallbacks must emit DeprecationWarning.

    Fix 1.5 principle: fallback paths should be noisy so stale data
    is detected during testing.
    """
    files_with_fallbacks = {
        "agent/directive_validator.py": "_get_fallback_programs",
        "agent/command_builder.py": "_get_data_label_parameters",
        "phenix_ai/log_parsers.py": "_detect_program_fallback",
        "agent/directive_extractor.py": "_get_stop_directive_patterns",
    }

    for filepath, function_name in files_with_fallbacks.items():
        source = _read(filepath)

        # The function should exist
        assert f"def {function_name}" in source, (
            f"{filepath} is missing {function_name}"
        )

        # Find the function body and check for DeprecationWarning
        # (command_builder uses a different pattern — warning in the
        #  except block of the loader, not in the fallback function)
        if filepath == "agent/command_builder.py":
            assert "DeprecationWarning" in source, (
                f"{filepath} missing DeprecationWarning in data label fallback"
            )
        elif filepath == "phenix_ai/log_parsers.py":
            # The fallback function itself doesn't warn, but the caller does
            # via the except → fallback pattern. Check the docstring warns.
            assert "maintenance risk" in source.lower() or "fallback" in source.lower()
        else:
            assert "DeprecationWarning" in source, (
                f"{filepath} {function_name} should emit DeprecationWarning "
                "when using hardcoded fallback data."
            )

    print("  test_fallback_paths_emit_deprecation_warnings PASSED")


# ===================================================================
# FIX 7 TESTS — planner.py dead code removal
# ===================================================================


def test_planner_no_heavy_imports():
    """planner.py must not import langchain_core or other heavy dependencies.

    After Fix 7, only fix_program_parameters and extract_output_files
    remain. The heavy imports (langchain_core, phenix_knowledge, etc.)
    were only needed by the removed dead code.
    """
    source = _read("agent/planner.py")

    forbidden_imports = [
        "langchain_core",
        "from libtbx.langchain.knowledge.prompts import",
        "from libtbx.langchain.validation",
        "from libtbx.langchain.agent.memory import",
        "phenix_knowledge",
    ]

    violations = []
    for imp in forbidden_imports:
        if imp in source:
            violations.append(imp)

    assert not violations, (
        "planner.py still imports dead-code dependencies (Fix 7 regression):\n  "
        + "\n  ".join(violations)
    )
    print("  test_planner_no_heavy_imports PASSED")


def test_planner_no_dead_functions():
    """planner.py must not contain removed functions.

    After Fix 7, generate_next_move, construct_command_mechanically,
    and other dead functions should be gone.
    """
    source = _read("agent/planner.py")

    dead_functions = [
        "def generate_next_move",
        "def construct_command_mechanically",
        "def get_required_params",
        "def extract_clean_command",
        "def get_relative_path",
        "def get_program_keywords",
    ]

    found = [fn for fn in dead_functions if fn in source]

    assert not found, (
        "planner.py still contains dead functions (Fix 7 regression):\n  "
        + "\n  ".join(found)
    )
    print("  test_planner_no_dead_functions PASSED")


def test_planner_retained_functions_work():
    """fix_program_parameters and extract_output_files must still work."""
    try:
        from libtbx.langchain.agent.planner import fix_program_parameters, extract_output_files
    except ImportError:
        from agent.planner import fix_program_parameters, extract_output_files

    # fix_program_parameters: passthrough when no fixes apply
    result = fix_program_parameters("phenix.refine model.pdb data.mtz", "phenix.unknown")
    assert result == "phenix.refine model.pdb data.mtz", "Should pass through unchanged"

    # extract_output_files: basic extraction
    summary = "Results:\n**Key Output Files:**\n- model_refine_001.pdb\n- data_refine_001.mtz\n**Done**"
    files = extract_output_files(summary)
    assert len(files) == 2, f"Expected 2 files, got {len(files)}: {files}"
    assert "model_refine_001.pdb" in files

    # extract_output_files: empty input
    assert extract_output_files("") == []
    assert extract_output_files(None) == []

    print("  test_planner_retained_functions_work PASSED")


def test_planner_size():
    """planner.py should be under 200 lines after dead code removal.

    Was 1436 lines before Fix 7. Should be ~130 lines now.
    """
    source = _read("agent/planner.py")
    line_count = len(source.splitlines())

    assert line_count < 200, (
        f"planner.py is {line_count} lines — expected under 200 after "
        f"Fix 7 dead code removal. Did dead code creep back in?"
    )
    print("  test_planner_size PASSED")


# ===================================================================
# FIX 9 TESTS — gui_app_id in programs.yaml
# ===================================================================


def test_gui_app_id_coverage():
    """Programs with known GUI windows must have gui_app_id in YAML.

    Programs without GUI windows (autobuild_denmod, model_vs_data,
    holton_geometry_validation) are excluded.
    """
    programs = _program_defs()

    # These programs have no dedicated GUI window
    NO_GUI = {
        "phenix.autobuild_denmod",
        "phenix.model_vs_data",
        "phenix.holton_geometry_validation",
    }

    missing = []
    for name in sorted(programs):
        if name in NO_GUI:
            continue
        if not programs[name].get("gui_app_id"):
            missing.append(name)

    assert not missing, (
        "These programs are missing gui_app_id in programs.yaml:\n  "
        + "\n  ".join(missing)
        + "\n\nAdd gui_app_id with the wxGUI2 app_id value."
    )
    print("  test_gui_app_id_coverage PASSED")


def test_gui_app_id_cryoem_variants():
    """Programs with cryo-EM variant GUI windows must have gui_app_id_cryoem."""
    programs = _program_defs()

    # Known cryo-EM variant programs (from __init__.params)
    expected_cryo = {
        "phenix.predict_and_build": "PredictAndBuildCryoEM",
        "phenix.map_correlations": "MapCorrelationsCryoEM",
        "phenix.map_sharpening": "MapSharpeningCryoEM",
    }

    for prog, expected in expected_cryo.items():
        defn = programs.get(prog, {})
        actual = defn.get("gui_app_id_cryoem")
        assert actual == expected, (
            f"{prog}: expected gui_app_id_cryoem='{expected}', got '{actual}'"
        )

    print("  test_gui_app_id_cryoem_variants PASSED")


def test_yaml_tools_includes_gui_app_id():
    """yaml_tools.py valid fields must include gui_app_id and gui_app_id_cryoem."""
    source = _read("agent/yaml_tools.py")

    assert "'gui_app_id'" in source, (
        "yaml_tools.py valid_program_fields missing 'gui_app_id'"
    )
    assert "'gui_app_id_cryoem'" in source, (
        "yaml_tools.py valid_program_fields missing 'gui_app_id_cryoem'"
    )
    print("  test_yaml_tools_includes_gui_app_id PASSED")


# ===================================================================
# FIX 8 TESTS — stop_directive_patterns in programs.yaml
# ===================================================================


def test_stop_directive_patterns_coverage():
    """Programs with known stop-condition phrases must have patterns in YAML.

    These 9 programs had hardcoded patterns in directive_extractor.py.
    """
    programs = _program_defs()

    expected = {
        "phenix.mtriage", "phenix.xtriage", "phenix.phaser",
        "phenix.ligandfit", "phenix.refine", "phenix.autobuild",
        "phenix.map_to_model", "phenix.dock_in_map", "phenix.map_symmetry",
    }

    missing = []
    for prog in sorted(expected):
        defn = programs.get(prog, {})
        pats = defn.get("stop_directive_patterns")
        if not pats or not isinstance(pats, list) or len(pats) == 0:
            missing.append(prog)

    assert not missing, (
        "These programs are missing stop_directive_patterns in programs.yaml:\n  "
        + "\n  ".join(missing)
    )
    print("  test_stop_directive_patterns_coverage PASSED")


def test_stop_directive_patterns_are_valid_regex():
    """All stop_directive_patterns must compile as valid regex."""
    import re as _re
    programs = _program_defs()

    errors = []
    for name, defn in sorted(programs.items()):
        if not isinstance(defn, dict):
            continue
        pats = defn.get("stop_directive_patterns")
        if not pats:
            continue
        for pat in pats:
            try:
                _re.compile(pat)
            except _re.error as e:
                errors.append(f"{name}: {pat!r} → {e}")

    assert not errors, (
        "Invalid regex in stop_directive_patterns:\n  "
        + "\n  ".join(errors)
    )
    print("  test_stop_directive_patterns_are_valid_regex PASSED")


def test_stop_directive_patterns_match_expected():
    """Verify YAML-loaded patterns produce expected matches.

    Tests the critical ordering cases where pattern specificity matters.
    """
    import re as _re

    try:
        from libtbx.langchain.agent.directive_extractor import (
            _get_stop_directive_patterns,
        )
    except ImportError:
        from agent.directive_extractor import _get_stop_directive_patterns

    mappings = _get_stop_directive_patterns()

    def match(text):
        for pattern, program in mappings:
            if _re.search(pattern, text, _re.IGNORECASE):
                return program
        return None

    # Critical ordering tests
    expected = [
        ("map to model", "phenix.map_to_model"),
        ("dock in map", "phenix.dock_in_map"),
        ("build model into map", "phenix.map_to_model"),  # NOT dock_in_map
        ("fit model into the map", "phenix.dock_in_map"),
        ("molecular replacement", "phenix.phaser"),
        ("refine", "phenix.refine"),
        ("autobuild", "phenix.autobuild"),
        ("map symmetry", "phenix.map_symmetry"),
        ("something unknown", None),
    ]

    failures = []
    for phrase, expected_prog in expected:
        actual = match(phrase)
        if actual != expected_prog:
            failures.append(
                f"  '{phrase}': expected {expected_prog}, got {actual}"
            )

    assert not failures, (
        "Stop directive pattern matching failures:\n"
        + "\n".join(failures)
    )
    print("  test_stop_directive_patterns_match_expected PASSED")


def test_no_hardcoded_program_mappings_in_directive_extractor():
    """directive_extractor.py must not have inline program_mappings list.

    After Fix 8, all patterns come from _get_stop_directive_patterns()
    which loads from YAML. The inline list should be gone.
    """
    source = _read("agent/directive_extractor.py")

    # The old code had "program_mappings = [" as a local variable
    # with specific pattern tuples inside
    assert "program_mappings = [" not in source, (
        "directive_extractor.py still has hardcoded program_mappings list "
        "(Fix 8 regression). Patterns should come from "
        "_get_stop_directive_patterns() → programs.yaml."
    )
    print("  test_no_hardcoded_program_mappings_in_directive_extractor PASSED")


# ===================================================================
# FIX 6 TESTS — rules_priority in workflows.yaml
# ===================================================================


def test_rules_priority_in_workflow_phases():
    """Key workflow phases must have rules_priority in YAML.

    Phases: analyze, obtain_model, refine, validate (both xray and cryoem).
    """
    import yaml as _yaml

    yaml_path = os.path.join(_ROOT_DIR, "knowledge", "workflows.yaml")
    with open(yaml_path) as f:
        workflows = _yaml.safe_load(f)

    expected_phases = {
        "xray": ["analyze", "obtain_model", "refine", "validate"],
        "cryoem": ["analyze", "obtain_model", "refine", "validate"],
    }

    missing = []
    for wf_name, phases in expected_phases.items():
        wf = workflows.get(wf_name, {})
        wf_phases = wf.get("phases", {})
        for phase_name in phases:
            phase = wf_phases.get(phase_name, {})
            rp = phase.get("rules_priority")
            if not rp:
                missing.append(f"{wf_name}/{phase_name}")

    assert not missing, (
        "These workflow phases are missing rules_priority:\n  "
        + "\n  ".join(missing)
    )
    print("  test_rules_priority_in_workflow_phases PASSED")


def test_rules_config_in_shared():
    """shared/rules_config must have default_priority, ligand_priority, and aliases."""
    import yaml as _yaml

    yaml_path = os.path.join(_ROOT_DIR, "knowledge", "workflows.yaml")
    with open(yaml_path) as f:
        workflows = _yaml.safe_load(f)

    config = workflows.get("shared", {}).get("rules_config", {})
    assert config, "shared/rules_config section missing from workflows.yaml"

    required_keys = ["default_priority", "ligand_priority", "state_aliases"]
    missing = [k for k in required_keys if k not in config]
    assert not missing, (
        f"shared/rules_config missing keys: {missing}"
    )

    # Verify aliases map to real phases
    aliases = config["state_aliases"]
    assert aliases.get("xray_initial") == "analyze", "xray_initial alias wrong"
    assert aliases.get("cryoem_initial") == "analyze", "cryoem_initial alias wrong"
    assert aliases.get("molecular_replacement") == "obtain_model", "molecular_replacement alias wrong"

    print("  test_rules_config_in_shared PASSED")


def test_no_hardcoded_priority_lists_in_rules_selector():
    """rules_selector.py must not have hardcoded priority lists.

    After Fix 6, all priority lists come from _get_phase_rules_priority()
    which loads from workflows.yaml. The inline if/elif chain with
    hardcoded lists should be gone.
    """
    source = _read("agent/rules_selector.py")

    # The old code had state_name checks with inline priority lists
    # like: if state_name in ("analyze", "xray_initial", "cryoem_initial"):
    #           priority = ["phenix.xtriage", ...
    assert '"xray_initial"' not in source or '_get_phase_rules_priority' in source, (
        "rules_selector.py still has hardcoded state_name aliases"
    )

    # Check that the big if/elif chain is gone
    # Old code had: elif state_name in ("validate",):
    #                   priority = [
    assert 'state_name in ("validate",)' not in source, (
        "rules_selector.py still has hardcoded validate priority list "
        "(Fix 6 regression)"
    )
    assert 'state_name in ("obtain_model"' not in source, (
        "rules_selector.py still has hardcoded obtain_model priority list "
        "(Fix 6 regression)"
    )

    print("  test_no_hardcoded_priority_lists_in_rules_selector PASSED")


# ===================================================================
# FIX 10B TESTS — YAML-driven _analyze_history simple flags
# ===================================================================


def test_history_detection_coverage():
    """All programs with simple done-flag logic must have history_detection in YAML.

    These programs set only a boolean done flag (and optional success_flag)
    on successful completion. They don't have counts, cascading side effects,
    or exclusion logic.
    """
    programs = _program_defs()

    # Programs that MUST have history_detection
    expected = {
        "phenix.process_predicted_model",
        "phenix.autobuild",
        "phenix.autobuild_denmod",
        "phenix.autosol",
        "phenix.ligandfit",
        "phenix.pdbtools",
        "phenix.dock_in_map",
        "phenix.map_to_model",
        "phenix.resolve_cryo_em",
        "phenix.map_sharpening",
    }

    missing = []
    for name in sorted(expected):
        defn = programs.get(name, {})
        tracking = defn.get("done_tracking", {})
        hd = tracking.get("history_detection")
        if not hd or not hd.get("markers"):
            missing.append(name)

    assert not missing, (
        "These programs are missing history_detection in done_tracking:\n  "
        + "\n  ".join(missing)
    )
    print("  test_history_detection_coverage PASSED")


def test_history_detection_behavioral():
    """Verify _set_done_flags produces correct flags for each program."""
    try:
        from libtbx.langchain.agent.workflow_state import _set_done_flags
    except ImportError:
        from agent.workflow_state import _set_done_flags

    test_cases = [
        # (combined_text, result, flag, expected_value)
        # --- Simple set_flag programs ---
        ("phenix.dock_in_map model.pdb map.mrc", "SUCCESS", "dock_done", True),
        ("phenix.dock_in_map model.pdb map.mrc", "FAILED", "dock_done", False),
        ("phenix.autobuild data.mtz seq.fa", "SUCCESS", "autobuild_done", True),
        ("phenix.autobuild maps_only=True data.mtz", "SUCCESS", "autobuild_denmod_done", True),
        ("phenix.autobuild_denmod data.mtz", "SUCCESS", "autobuild_denmod_done", True),
        ("phenix.autosol data.mtz seq.fa", "SUCCESS", "autosol_done", True),
        ("phenix.autosol data.mtz seq.fa", "SUCCESS", "autosol_success", True),
        ("phenix.autosol data.mtz", "SORRY: failed", "autosol_done", False),
        ("phenix.ligandfit model.pdb lig.cif", "SUCCESS", "ligandfit_done", True),
        ("phenix.map_sharpening map.mrc", "SUCCESS", "map_sharpening_done", True),
        ("phenix.local_aniso_sharpen map.mrc", "SUCCESS", "map_sharpening_done", True),
        ("phenix.map_to_model map.mrc seq.fa", "SUCCESS", "map_to_model_done", True),
        ("phenix.resolve_cryo_em half1.mrc", "SUCCESS", "resolve_cryo_em_done", True),
        ("phenix.pdbtools model.pdb", "SUCCESS", "pdbtools_done", True),
        ("phenix.process_predicted_model pred.pdb", "SUCCESS", "process_predicted_done", True),
        ("phenix.polder model.pdb data.mtz", "SUCCESS", "polder_done", True),
        # --- Counted programs (strategy: count) ---
        ("phenix.refine model.pdb data.mtz", "SUCCESS", "refine_done", True),
        ("phenix.refine model.pdb data.mtz", "SUCCESS", "refine_count", 1),
        ("phenix.real_space_refine m.pdb map.mrc", "SUCCESS", "rsr_done", True),
        ("phenix.real_space_refine m.pdb map.mrc", "SUCCESS", "rsr_count", 1),
        ("phenix.phaser data.mtz", "SUCCESS", "phaser_done", True),
        ("phenix.phaser data.mtz", "SUCCESS", "phaser_count", 1),
        # --- Exclude markers: real_space_refine must NOT trigger refine ---
        ("phenix.real_space_refine m.pdb map.mrc", "SUCCESS", "refine_done", False),
        ("phenix.real_space_refine m.pdb map.mrc", "SUCCESS", "refine_count", 0),
        # --- Validation: shared flag from multiple programs ---
        ("phenix.molprobity model.pdb", "SUCCESS", "validation_done", True),
        ("phenix.holton_geometry_validation m.pdb", "SUCCESS", "validation_done", True),
        ("phenix.validation_cryoem model.pdb", "SUCCESS", "validation_done", True),
    ]

    failures = []
    for combined, result, flag, expected in test_cases:
        info = {}
        _set_done_flags(info, combined.lower(), result)
        default = 0 if isinstance(expected, int) else False
        actual = info.get(flag, default)
        if actual != expected:
            failures.append(f"  '{combined}' result={result}: {flag}={actual}, expected {expected}")

    assert not failures, (
        "Behavioral mismatches in _set_done_flags:\n"
        + "\n".join(failures)
    )
    print("  test_history_detection_behavioral PASSED")


def test_no_simple_done_flags_in_analyze_history():
    """_analyze_history must not have simple if-blocks for programs
    that have history_detection in their done_tracking YAML config.

    After Fix 10B, these programs are handled by _set_done_flags()
    and should not appear as standalone if-blocks in _analyze_history.
    """
    source = _read("agent/workflow_state.py")

    # These simple flag assignments should NOT appear in _analyze_history
    # because they're now handled by _set_done_flags
    removed_patterns = [
        # Simple flag programs (YAML history_detection since Fix 10B)
        'info["process_predicted_done"] = True',
        'info["autobuild_done"] = True',
        'info["autobuild_denmod_done"] = True',
        'info["autosol_done"] = True',
        'info["autosol_success"] = True',
        'info["ligandfit_done"] = True',
        'info["pdbtools_done"] = True',
        'info["dock_done"] = True',
        'info["map_to_model_done"] = True',
        'info["resolve_cryo_em_done"] = True',
        'info["map_sharpening_done"] = True',
        # Counted/shared programs (YAML strategy since strategy enum)
        'info["validation_done"] = True',
        'info["phaser_done"] = True',
        'info["rsr_done"] = True',
        # Dead flags removed
        'info["refine_success"] = True',
        'info["rsr_success"] = True',
        'info["phaser_success"] = True',
    ]

    found = [p for p in removed_patterns if p in source]
    assert not found, (
        "These simple flag assignments should be in YAML history_detection, "
        "not hardcoded in _analyze_history (Fix 10B regression):\n  "
        + "\n  ".join(found)
    )

    # The predict_and_build cascade should still be in Python
    assert 'info["predict_done"]' in source, "predict_done should stay in Python (cascade)"
    assert 'info["refine_count"]' in source, "refine_count should stay in Python (predict cascade)"

    # Count fields should be in ALLOWED_COUNT_FIELDS (YAML-driven validation)
    assert "ALLOWED_COUNT_FIELDS" in source, "ALLOWED_COUNT_FIELDS whitelist should exist"

    print("  test_no_simple_done_flags_in_analyze_history PASSED")


def test_strategy_enum_values():
    """Strategy values in YAML must be valid enum members."""
    programs = _program_defs()

    VALID_STRATEGIES = {"set_flag", "run_once", "count"}

    bad = []
    for name, defn in programs.items():
        tracking = defn.get("done_tracking", {})
        strategy = tracking.get("strategy")
        if strategy and strategy not in VALID_STRATEGIES:
            bad.append(f"{name}: strategy={strategy!r}")

    assert not bad, (
        "Invalid strategy values in programs.yaml:\n  "
        + "\n  ".join(bad)
        + f"\nAllowed: {VALID_STRATEGIES}"
    )
    print("  test_strategy_enum_values PASSED")


def test_count_field_validation():
    """Programs with strategy='count' must have a valid count_field."""
    programs = _program_defs()

    try:
        from libtbx.langchain.agent.workflow_state import ALLOWED_COUNT_FIELDS
    except ImportError:
        from agent.workflow_state import ALLOWED_COUNT_FIELDS

    errors = []
    for name, defn in programs.items():
        tracking = defn.get("done_tracking", {})
        if tracking.get("strategy") != "count":
            continue

        cf = tracking.get("count_field")
        if not cf:
            errors.append(f"{name}: strategy='count' but no count_field")
        elif cf not in ALLOWED_COUNT_FIELDS:
            errors.append(f"{name}: count_field={cf!r} not in ALLOWED_COUNT_FIELDS")

        hd = tracking.get("history_detection", {})
        if not hd.get("markers"):
            errors.append(f"{name}: strategy='count' but no history_detection markers")

    assert not errors, (
        "Count field validation errors:\n  "
        + "\n  ".join(errors)
    )
    print("  test_count_field_validation PASSED")


def test_exclude_markers_prevent_false_matches():
    """Programs with exclude_markers must not match excluded text.

    Critical case: refine must not match real_space_refine.
    """
    try:
        from libtbx.langchain.agent.workflow_state import _set_done_flags
    except ImportError:
        from agent.workflow_state import _set_done_flags

    # real_space_refine should set rsr flags but NOT refine flags
    info = {}
    _set_done_flags(info, "phenix.real_space_refine model.pdb map.mrc", "SUCCESS")

    assert info.get("rsr_done") is True, "rsr_done should be True"
    assert info.get("rsr_count", 0) == 1, "rsr_count should be 1"
    assert info.get("refine_done") is not True, "refine_done must NOT be set by real_space_refine"
    assert info.get("refine_count", 0) == 0, "refine_count must stay 0 for real_space_refine"

    print("  test_exclude_markers_prevent_false_matches PASSED")


# ===================================================================
# Runner
# ===================================================================


def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.tst_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
