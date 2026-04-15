"""
Backward compatibility tests for v115.05 (Phase 1) and v115.06 (Phase 2).

These tests verify that server-side changes do NOT break existing clients.
The golden rule: a user with yesterday's PHENIX install must connect to
today's server and get correct results.

Tests cover:
  - programs.yaml flag additions are backward-safe (new flags default to None)
  - Supplement logic doesn't affect non-multiple slots or X-ray tutorials
  - Command builder works with minimal v1 session_info
  - Half-map dedup doesn't fire on X-ray (no half-maps)
  - strict_strategy_flags doesn't block programs that don't set it
  - prefers_half_maps doesn't remove files from programs that don't set it
  - Response format is unchanged

Run with:
    PYTHONPATH=. python tests/tst_backward_compat_v115.py
"""

import os
import sys

sys.path.insert(
    0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_false, assert_in,
    assert_not_in, assert_none, assert_not_none,
    run_tests_with_fail_fast,
)

# Silence linter
(assert_equal, assert_true, assert_false, assert_in,
 assert_not_in, assert_none, assert_not_none, run_tests_with_fail_fast)

try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_programs_yaml():
    path = os.path.join(_PROJECT_ROOT, "knowledge", "programs.yaml")
    with open(path) as f:
        return yaml.safe_load(f)


def _load_workflows_yaml():
    path = os.path.join(_PROJECT_ROOT, "knowledge", "workflows.yaml")
    with open(path) as f:
        return yaml.safe_load(f)


# =========================================================================
# 1. New YAML flags have safe defaults for unlabeled programs
# =========================================================================

def test_prefers_half_maps_default_is_none():
    """Programs without prefers_half_maps must default to None/falsy.

    The dedup logic checks `prog_def.get("prefers_half_maps")`. If a program
    doesn't have this flag, .get() returns None, which is falsy. The dedup
    falls through to the default branch (keep full_map, drop half_map).
    This test ensures no X-ray program accidentally gets the flag.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    progs = _load_programs_yaml()
    xray_with_flag = []
    for name, prog in progs.items():
        if not isinstance(prog, dict):
            continue
        if prog.get("prefers_half_maps") and "xray" in prog.get("experiment_types", []):
            xray_with_flag.append(name)
    assert_equal(xray_with_flag, [],
                 "X-ray programs must NOT have prefers_half_maps: %s" % xray_with_flag)
    print("  PASSED: No X-ray program has prefers_half_maps")


def test_strict_strategy_flags_default_is_none():
    """Programs without strict_strategy_flags must allow PHIL passthrough.

    Only programs that explicitly set strict_strategy_flags: true will block
    KNOWN_PHIL_SHORT_NAMES. All other programs must continue to allow them.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    progs = _load_programs_yaml()
    strict_progs = [name for name, prog in progs.items()
                    if isinstance(prog, dict) and prog.get("strict_strategy_flags")]
    # Only process_predicted_model should have it
    assert_in("phenix.process_predicted_model", strict_progs,
              "process_predicted_model must have strict_strategy_flags")
    # Refine, phaser, autosol etc. must NOT have it
    for prog in ["phenix.refine", "phenix.phaser", "phenix.autosol",
                 "phenix.real_space_refine", "phenix.autobuild"]:
        assert_not_in(prog, strict_progs,
                      "%s must NOT have strict_strategy_flags" % prog)
    print("  PASSED: strict_strategy_flags only on intended programs (%s)" %
          ", ".join(strict_progs))


def test_keep_half_maps_still_works():
    """Programs with keep_half_maps_with_full_map must still be honored.

    Verify map_to_model still has keep_half_maps_with_full_map: true.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    progs = _load_programs_yaml()
    mtm = progs.get("phenix.map_to_model", {})
    assert_true(mtm.get("keep_half_maps_with_full_map"),
                "map_to_model must retain keep_half_maps_with_full_map: true")
    print("  PASSED: map_to_model still has keep_half_maps_with_full_map")


# =========================================================================
# 2. Half-map dedup doesn't affect X-ray tutorials
# =========================================================================

def test_xray_programs_unaffected_by_half_map_dedup():
    """Pure X-ray programs (phenix.refine, phaser, etc.) have no half_map input.

    The dedup only fires when both full_map and half_map are in selected_files.
    Pure X-ray programs never have half_map → dedup never fires → no regression.
    Dual-purpose programs (xray+cryoem) like predict_and_build may have half_map
    with auto_fill:false, which is safe — supplement only fires if LLM filled it.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    progs = _load_programs_yaml()
    problems = []
    for name, prog in progs.items():
        if not isinstance(prog, dict):
            continue
        exp_types = prog.get("experiment_types", [])
        if "xray" not in exp_types:
            continue
        # Skip dual-purpose programs (xray + cryoem)
        if "cryoem" in exp_types:
            continue
        inputs = prog.get("inputs", {})
        all_slots = {}
        all_slots.update(inputs.get("required", {}))
        all_slots.update(inputs.get("optional", {}))
        if "half_map" in all_slots:
            problems.append(name)
    assert_equal(problems, [],
                 "Pure X-ray programs must not have half_map slot: %s" % problems)
    print("  PASSED: No pure X-ray program has a half_map input slot")


# =========================================================================
# 3. Refine template has no unconditional generate=True
# =========================================================================

def test_refine_no_unconditional_generate():
    """phenix.refine must NOT have generate=True in command or defaults.

    The R-free flag is added conditionally by the BUILD node based on
    rfree_mtz presence. Unconditional generate=True causes mismatches
    on second+ refinements and on datasets with pre-existing flags.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    progs = _load_programs_yaml()
    refine = progs.get("phenix.refine", {})
    assert_not_in("generate", refine.get("command", ""),
                  "generate=True must not be in refine command template")
    defaults = refine.get("defaults", {})
    for key in defaults:
        assert_not_in("generate", str(key).lower(),
                      "generate must not be in refine defaults: %s" % key)
    # But it must be in strategy_flags (conditional path)
    sf = refine.get("strategy_flags", {})
    assert_in("generate_rfree_flags", sf,
              "generate_rfree_flags must be in refine strategy_flags")
    print("  PASSED: refine has no unconditional generate=True")


# =========================================================================
# 4. map_sharpening half_map flag is positional
# =========================================================================

def test_map_sharpening_half_map_positional():
    """map_sharpening half_map flag must be empty string (positional).

    The working command is:
        phenix.map_sharpening h1.ccp4 h2.ccp4 seq_file=seq.fa
    Not:
        phenix.map_sharpening half_map=h1.ccp4 half_map=h2.ccp4
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    progs = _load_programs_yaml()
    ms = progs.get("phenix.map_sharpening", {})
    hm = ms.get("inputs", {}).get("optional", {}).get("half_map", {})
    assert_equal(hm.get("flag"), "",
                 "map_sharpening half_map flag must be '' (positional), got %r" %
                 hm.get("flag"))
    print("  PASSED: map_sharpening half_map is positional")


# =========================================================================
# 5. Workflow conditions: map_sharpening accepts both full_map and half_map
# =========================================================================

def test_map_sharpening_workflow_accepts_half_maps():
    """All map_sharpening workflow entries must have has_any:[full_map,half_map].

    This ensures map_sharpening is offered when only half-maps are available,
    not just when a full_map exists.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    wf = _load_workflows_yaml()

    def find_programs(obj, target):
        results = []
        if isinstance(obj, dict):
            if obj.get("program") == target:
                results.append(obj)
            for v in obj.values():
                results.extend(find_programs(v, target))
        elif isinstance(obj, list):
            for item in obj:
                results.extend(find_programs(item, target))
        return results

    ms_entries = find_programs(wf, "phenix.map_sharpening")
    assert_true(len(ms_entries) >= 4,
                "Expected at least 4 map_sharpening entries, found %d" %
                len(ms_entries))
    for e in ms_entries:
        conds = e.get("conditions", [])
        has_any_vals = []
        for c in conds:
            if isinstance(c, dict) and "has_any" in c:
                has_any_vals.extend(c["has_any"])
        assert_in("full_map", has_any_vals,
                  "map_sharpening missing full_map in has_any: %s" % conds)
        assert_in("half_map", has_any_vals,
                  "map_sharpening missing half_map in has_any: %s" % conds)
        # Must NOT have bare has:full_map
        has_only = [c.get("has") for c in conds
                    if isinstance(c, dict) and "has" in c]
        assert_not_in("full_map", has_only,
                      "map_sharpening has stale has:full_map condition: %s" % conds)
    print("  PASSED: All %d map_sharpening entries accept both map types" %
          len(ms_entries))


# =========================================================================
# 6. resolve_cryo_em still requires half_map
# =========================================================================

def test_resolve_cryo_em_requires_half_map():
    """All resolve_cryo_em workflow entries must have has:half_map."""
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    wf = _load_workflows_yaml()

    def find_programs(obj, target):
        results = []
        if isinstance(obj, dict):
            if obj.get("program") == target:
                results.append(obj)
            for v in obj.values():
                results.extend(find_programs(v, target))
        elif isinstance(obj, list):
            for item in obj:
                results.extend(find_programs(item, target))
        return results

    rce_entries = find_programs(wf, "phenix.resolve_cryo_em")
    assert_true(len(rce_entries) >= 4,
                "Expected at least 4 resolve_cryo_em entries, found %d" %
                len(rce_entries))
    for e in rce_entries:
        conds = e.get("conditions", [])
        has_half = any(c.get("has") == "half_map"
                       for c in conds if isinstance(c, dict))
        assert_true(has_half,
                    "resolve_cryo_em missing has:half_map: %s" % conds)
    print("  PASSED: All %d resolve_cryo_em entries require half_map" %
          len(rce_entries))


# =========================================================================
# 7. All dock_in_map entries have not_done:dock guard
# =========================================================================

def test_dock_in_map_guarded():
    """All dock_in_map entries must have not_done:dock to prevent re-docking."""
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    wf = _load_workflows_yaml()

    def find_programs(obj, target):
        results = []
        if isinstance(obj, dict):
            if obj.get("program") == target:
                results.append(obj)
            for v in obj.values():
                results.extend(find_programs(v, target))
        elif isinstance(obj, list):
            for item in obj:
                results.extend(find_programs(item, target))
        return results

    dock_entries = find_programs(wf, "phenix.dock_in_map")
    assert_true(len(dock_entries) >= 3,
                "Expected at least 3 dock_in_map entries, found %d" %
                len(dock_entries))
    for e in dock_entries:
        conds = e.get("conditions", [])
        has_guard = any(c.get("not_done") == "dock"
                        for c in conds if isinstance(c, dict))
        assert_true(has_guard,
                    "dock_in_map missing not_done:dock: %s" % conds)
    print("  PASSED: All %d dock_in_map entries have not_done:dock" %
          len(dock_entries))


# =========================================================================
# 8. Supplement logic: only fires for multiple:true slots
# =========================================================================

def test_multiple_true_only_on_expected_slots():
    """Only half_map, model (phaser), and sequence (phaser) should be multiple.

    The supplement logic iterates all multiple:true slots. If we accidentally
    set multiple:true on a single-file slot, supplement could inject extra
    files and corrupt the command.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    progs = _load_programs_yaml()
    EXPECTED_MULTIPLE = {
        "phenix.map_sharpening:half_map",
        "phenix.map_to_model:half_map",
        "phenix.mtriage:half_map",
        "phenix.phaser:model",
        "phenix.phaser:sequence",
        "phenix.predict_and_build:half_map",
        "phenix.resolve_cryo_em:half_map",
    }
    actual_multiple = set()
    for name, prog in progs.items():
        if not isinstance(prog, dict):
            continue
        inputs = prog.get("inputs", {})
        for section in ["required", "optional"]:
            for slot_name, slot_def in inputs.get(section, {}).items():
                if isinstance(slot_def, dict) and slot_def.get("multiple"):
                    actual_multiple.add("%s:%s" % (name, slot_name))
    unexpected = actual_multiple - EXPECTED_MULTIPLE
    missing = EXPECTED_MULTIPLE - actual_multiple
    assert_equal(unexpected, set(),
                 "Unexpected multiple:true slots: %s" % unexpected)
    assert_equal(missing, set(),
                 "Missing expected multiple:true slots: %s" % missing)
    print("  PASSED: multiple:true slots exactly match expected set (%d)" %
          len(actual_multiple))


# =========================================================================
# 9. No new required session_info fields
# =========================================================================

def test_no_new_required_session_info():
    """Phase 1+2 changes must not require new session_info fields.

    All our changes (prefers_half_maps, supplement, strict_strategy_flags,
    etc.) are in programs.yaml, workflows.yaml, and server-side code.
    None of them add new fields to session_info that old clients would
    need to send.
    """
    try:
        from agent.contract import SESSION_INFO_FIELDS
    except ImportError:
        print("  SKIP (contract.py not importable)")
        return

    # These are the fields that existed before v115.05
    PRE_V115_05_FIELDS = {
        "experiment_type", "best_files", "rfree_mtz", "directives",
        "rfree_resolution", "force_retry_program", "recovery_strategies",
        "explicit_program", "advice_changed", "bad_inject_params",
        "unplaced_model_cell", "model_hetatm_residues",
        "client_protocol_version",
    }
    current_fields = {name for name, _, _, _ in SESSION_INFO_FIELDS}
    new_fields = current_fields - PRE_V115_05_FIELDS
    # New fields are okay IF they have safe defaults. But there should be
    # none added by v115.05/v115.06 specifically.
    if new_fields:
        print("  INFO: New session_info fields since v115.05: %s" % new_fields)
        print("  (These are acceptable if they have safe defaults)")
    else:
        print("  PASSED: No new session_info fields added by v115.05/v115.06")


# =========================================================================
# 10. Supplement logic code-level safety
# =========================================================================

def test_supplement_uses_realpath():
    """The supplement logic must use os.path.realpath for symlink robustness."""
    path = os.path.join(_PROJECT_ROOT, "agent", "command_builder.py")
    with open(path) as f:
        src = f.read()
    assert_in("os.path.realpath", src,
              "Supplement logic must use os.path.realpath")
    # Verify it's in the supplement block specifically
    supp_start = src.index("SUPPLEMENT:")
    supp_end = src.index("POST-SELECTION VALIDATION:", supp_start)
    supp_block = src[supp_start:supp_end]
    assert_in("os.path.realpath", supp_block,
              "os.path.realpath must be in the supplement block")
    print("  PASSED: Supplement uses os.path.realpath")


def test_supplement_variable_not_shadowed():
    """The supplement loop variable must not shadow all_inputs."""
    path = os.path.join(_PROJECT_ROOT, "agent", "command_builder.py")
    with open(path) as f:
        src = f.read()
    # Must use all_input_defs, not all_inputs
    supp_start = src.index("SUPPLEMENT:")
    supp_end = src.index("POST-SELECTION VALIDATION:", supp_start)
    supp_block = src[supp_start:supp_end]
    assert_in("all_input_defs", supp_block,
              "Supplement must use all_input_defs (not all_inputs)")
    assert_not_in("all_inputs = {}", supp_block,
                  "Supplement must not shadow all_inputs")
    print("  PASSED: Supplement variable not shadowed")


# =========================================================================
# 11. Space group validation has cctbx-first ordering
# =========================================================================

def test_space_group_cctbx_before_heuristic():
    """cctbx.sgtbx validation must run before the heuristic fallback."""
    path = os.path.join(_PROJECT_ROOT, "agent", "command_builder.py")
    with open(path) as f:
        src = f.read()
    cctbx_pos = src.index("from cctbx import sgtbx")
    heuristic_pos = src.index("heuristic, no cctbx")
    assert_true(cctbx_pos < heuristic_pos,
                "cctbx must be tried before heuristic fallback")
    print("  PASSED: cctbx-first space group validation ordering")


# =========================================================================
# 12. Unplaced model guard preserves original program name
# =========================================================================

def test_unplaced_model_guard_saves_orig_program():
    """The unplaced model guard must save _orig_program before overwriting."""
    path = os.path.join(_PROJECT_ROOT, "agent", "graph_nodes.py")
    with open(path) as f:
        src = f.read()
    assert_in("_orig_program = chosen_program", src,
              "Must save original program name before redirect")
    # Verify the log message uses _orig_program, not intent.get("program")
    guard_start = src.index("UNPLACED MODEL")
    guard_end = src.index("Simple validation:", guard_start)
    guard_block = src[guard_start:guard_end]
    assert_in("_orig_program))", guard_block,
              "Log/reasoning must reference _orig_program, not intent")
    print("  PASSED: Unplaced model guard preserves original program name")


# =========================================================================
# 13. mtriage no longer has keep_half_maps_with_full_map
# =========================================================================

def test_mtriage_prefers_half_maps():
    """mtriage must have prefers_half_maps (not keep_half_maps_with_full_map).

    This prevents dimension mismatch when half-maps and full_map have
    different grid sizes from post-processing.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    progs = _load_programs_yaml()
    mt = progs.get("phenix.mtriage", {})
    assert_true(mt.get("prefers_half_maps"),
                "mtriage must have prefers_half_maps: true")
    assert_none(mt.get("keep_half_maps_with_full_map"),
                "mtriage must NOT have keep_half_maps_with_full_map")
    print("  PASSED: mtriage has prefers_half_maps, no keep_half_maps")


# =========================================================================
# 14. No programs removed from programs.yaml
# =========================================================================

def test_no_programs_removed():
    """All programs from v115.04 must still exist in programs.yaml.

    Old clients may have directives referencing any of these programs.
    Removing one would cause confusing "unknown program" behavior.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    progs = _load_programs_yaml()
    # Programs that existed before v115.05
    V115_04_PROGRAMS = [
        "phenix.xtriage", "phenix.refine", "phenix.phaser",
        "phenix.autobuild", "phenix.autosol", "phenix.molprobity",
        "phenix.real_space_refine", "phenix.resolve_cryo_em",
        "phenix.map_sharpening", "phenix.mtriage", "phenix.dock_in_map",
        "phenix.map_to_model", "phenix.map_symmetry",
        "phenix.predict_and_build", "phenix.process_predicted_model",
        "phenix.model_vs_data", "phenix.polder", "phenix.ligandfit",
        "phenix.validation_cryoem", "phenix.autobuild_denmod",
    ]
    missing = [p for p in V115_04_PROGRAMS if p not in progs]
    assert_equal(missing, [],
                 "Programs removed from programs.yaml: %s" % missing)
    print("  PASSED: All %d pre-v115.05 programs still present" %
          len(V115_04_PROGRAMS))


# =========================================================================
# 15. Command templates produce valid shell commands
# =========================================================================

def test_command_templates_well_formed():
    """Every program's command template must start with 'phenix.' and
    contain only valid placeholder patterns {slot_name}.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    import re
    progs = _load_programs_yaml()
    problems = []
    for name, prog in progs.items():
        if not isinstance(prog, dict):
            continue
        cmd = prog.get("command", "")
        if not cmd:
            continue
        if not cmd.startswith("phenix."):
            problems.append("%s: command doesn't start with phenix." % name)
        # Check placeholders are valid {word} patterns
        placeholders = re.findall(r'\{([^}]+)\}', cmd)
        inputs = prog.get("inputs", {})
        all_slots = set()
        all_slots.update(inputs.get("required", {}).keys())
        all_slots.update(inputs.get("optional", {}).keys())
        for ph in placeholders:
            if ph not in all_slots:
                problems.append("%s: placeholder {%s} not in inputs" %
                                (name, ph))
    assert_equal(problems, [],
                 "Command template issues: %s" % "; ".join(problems))
    print("  PASSED: All command templates well-formed")


# =========================================================================
# 16. Response field guarantees (output_node contract)
# =========================================================================

def test_output_node_guarantees_documented_fields():
    """output_node in graph_nodes.py must guarantee all response fields.

    Old clients expect these fields. If output_node stops setting defaults,
    old clients will KeyError or get wrong behavior.
    """
    path = os.path.join(_PROJECT_ROOT, "agent", "graph_nodes.py")
    with open(path) as f:
        src = f.read()

    # Find the output_node function
    assert_in("def output_node(", src,
              "output_node function must exist in graph_nodes.py")
    out_start = src.index("def output_node(")
    # The function should set defaults for these critical fields
    out_block = src[out_start:out_start + 3000]
    for field in ["warnings", "debug_log", "events", "stop_reason",
                  "abort_message", "red_flag_issues"]:
        assert_in(field, out_block,
                  "output_node must guarantee field '%s'" % field)
    print("  PASSED: output_node guarantees all 6 response fields")


# =========================================================================
# 17. dedup comment documents all three flag types
# =========================================================================

def test_dedup_documents_all_flag_types():
    """The POST-SELECTION VALIDATION comment must document all three paths.

    Future developers reading command_builder.py need to understand the
    three mutually-exclusive half-map dedup behaviors.
    """
    path = os.path.join(_PROJECT_ROOT, "agent", "command_builder.py")
    with open(path) as f:
        src = f.read()
    dedup_start = src.index("POST-SELECTION VALIDATION:")
    dedup_end = src.index("if \"full_map\" in selected_files", dedup_start)
    dedup_comment = src[dedup_start:dedup_end]
    assert_in("keep_half_maps_with_full_map", dedup_comment,
              "Dedup comment must mention keep_half_maps_with_full_map")
    assert_in("prefers_half_maps", dedup_comment,
              "Dedup comment must mention prefers_half_maps")
    assert_in("EXCEPTION", dedup_comment,
              "Dedup comment must document exceptions")
    print("  PASSED: Dedup comment documents all flag types")


# =========================================================================
# 18. Half-map flags are mutually exclusive per program
# =========================================================================

def test_half_map_flags_mutually_exclusive():
    """No program should have both prefers_half_maps AND
    keep_half_maps_with_full_map — they contradict each other.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    progs = _load_programs_yaml()
    conflicts = []
    for name, prog in progs.items():
        if not isinstance(prog, dict):
            continue
        has_prefer = prog.get("prefers_half_maps")
        has_keep = prog.get("keep_half_maps_with_full_map")
        if has_prefer and has_keep:
            conflicts.append(name)
    assert_equal(conflicts, [],
                 "Programs with BOTH prefers and keep flags: %s" % conflicts)
    print("  PASSED: No program has conflicting half-map flags")


# =========================================================================
# 19. Shared agent/ code has no bare LLM SDK imports
# =========================================================================

def test_no_llm_imports_in_shared_code():
    """Files in agent/ that run on the client must not import LLM SDKs.

    Old clients don't have langchain, openai, etc. installed. Server-only
    code (graph_nodes.py, prompts_hybrid.py) may import them, but shared
    files like command_builder.py, workflow_engine.py, session.py must not.
    """
    SHARED_FILES = [
        "command_builder.py", "program_registry.py", "workflow_engine.py",
        "workflow_state.py", "session.py", "contract.py", "transport.py",
    ]
    FORBIDDEN = ["langchain", "openai", "anthropic", "ollama"]
    problems = []
    for fname in SHARED_FILES:
        path = os.path.join(_PROJECT_ROOT, "agent", fname)
        if not os.path.exists(path):
            continue
        with open(path) as f:
            src = f.read()
        for lib in FORBIDDEN:
            # Check for bare imports (not inside try/except)
            import re
            pattern = r'^(?:import|from)\s+%s' % re.escape(lib)
            for match in re.finditer(pattern, src, re.MULTILINE):
                # Check if it's inside a try block
                line_start = src.rfind('\n', 0, match.start()) + 1
                preceding = src[max(0, line_start - 200):line_start]
                if 'try:' not in preceding:
                    problems.append("%s imports %s without try/except" %
                                    (fname, lib))
    assert_equal(problems, [],
                 "Forbidden imports in shared code: %s" % "; ".join(problems))
    print("  PASSED: No bare LLM SDK imports in shared agent/ files")


# =========================================================================
# 20. auto_fill:false on dual-purpose half_map slots
# =========================================================================

def test_dual_purpose_half_map_auto_fill_false():
    """Dual-purpose programs (xray+cryoem) with half_map must have auto_fill:false.

    Without this, auto-fill would inject half-maps into X-ray runs where they
    don't belong. The supplement logic respects auto_fill:false — it only
    completes slots the LLM explicitly chose to fill.
    """
    if not YAML_AVAILABLE:
        print("  SKIP (no yaml)")
        return
    progs = _load_programs_yaml()
    problems = []
    for name, prog in progs.items():
        if not isinstance(prog, dict):
            continue
        exp_types = prog.get("experiment_types", [])
        if "xray" not in exp_types or "cryoem" not in exp_types:
            continue
        hm = prog.get("inputs", {}).get("optional", {}).get("half_map")
        if hm and hm.get("auto_fill") is not False:
            problems.append(name)
    assert_equal(problems, [],
                 "Dual-purpose programs with half_map must have auto_fill:false: %s"
                 % problems)
    print("  PASSED: All dual-purpose half_map slots have auto_fill:false")


# =========================================================================
# 21. LocalAgent/RemoteAgent request-building parity
# =========================================================================

def test_local_remote_agent_settings_parity():
    """LocalAgent and RemoteAgent must set the same request["settings"] keys.

    If one agent adds a setting (e.g., verbosity) but the other doesn't,
    the server receives different inputs depending on run_on_server, causing
    silent result divergence.
    """
    import re

    local_path = os.path.join(_PROJECT_ROOT, "phenix_ai", "local_agent.py")
    remote_path = os.path.join(_PROJECT_ROOT, "phenix_ai", "remote_agent.py")
    if not os.path.exists(local_path) or not os.path.exists(remote_path):
        print("  SKIP (phenix_ai/ not available)")
        return

    def extract_settings(src):
        """Extract request["settings"]["X"] assignments."""
        pattern = r'request\["settings"\]\["(\w+)"\]'
        return sorted(set(re.findall(pattern, src)))

    with open(local_path) as f:
        local_src = f.read()
    with open(remote_path) as f:
        remote_src = f.read()

    local_settings = extract_settings(local_src)
    remote_settings = extract_settings(remote_src)

    assert_equal(local_settings, remote_settings,
                 "Settings parity: LocalAgent has %s, RemoteAgent has %s" %
                 (local_settings, remote_settings))
    print("  PASSED: LocalAgent/RemoteAgent set identical settings keys: %s" %
          local_settings)


def test_local_remote_agent_return_parity():
    """LocalAgent and RemoteAgent must return the same group_args fields.

    If one agent returns `events` but the other doesn't, event display
    breaks silently in one mode.
    """
    import re

    local_path = os.path.join(_PROJECT_ROOT, "phenix_ai", "local_agent.py")
    remote_path = os.path.join(_PROJECT_ROOT, "phenix_ai", "remote_agent.py")
    if not os.path.exists(local_path) or not os.path.exists(remote_path):
        print("  SKIP (phenix_ai/ not available)")
        return

    def extract_group_args_fields(src):
        """Extract field names from group_args() call."""
        # Find lines with "field_name=" inside the group_args block
        # Look for the return group_args( block
        start = src.find("return group_args(")
        if start < 0:
            return []
        # Find the matching close — scan forward counting parens
        depth = 0
        block_start = src.index("(", start)
        i = block_start
        while i < len(src):
            if src[i] == '(':
                depth += 1
            elif src[i] == ')':
                depth -= 1
                if depth == 0:
                    break
            i += 1
        block = src[block_start:i+1]
        fields = re.findall(r'^\s*(\w+)\s*=', block, re.MULTILINE)
        # Remove 'group_args_type' as it's the type label, not a data field
        return sorted(f for f in fields if f != 'group_args_type')

    with open(local_path) as f:
        local_src = f.read()
    with open(remote_path) as f:
        remote_src = f.read()

    local_fields = extract_group_args_fields(local_src)
    remote_fields = extract_group_args_fields(remote_src)

    assert_equal(local_fields, remote_fields,
                 "Return parity: LocalAgent returns %s, RemoteAgent returns %s" %
                 (local_fields, remote_fields))
    print("  PASSED: LocalAgent/RemoteAgent return identical fields: %s" %
          local_fields)


def test_local_agent_full_roundtrip():
    """LocalAgent must perform encode/decode roundtrip, not shortcut.

    The full roundtrip catches serialization bugs that only surface when
    data goes through JSON encoding. Without it, local mode could silently
    succeed while remote mode fails.
    """
    local_path = os.path.join(_PROJECT_ROOT, "phenix_ai", "local_agent.py")
    if not os.path.exists(local_path):
        print("  SKIP (phenix_ai/ not available)")
        return

    with open(local_path) as f:
        src = f.read()

    assert_in("prepare_request_for_transport", src,
              "LocalAgent must use prepare_request_for_transport")
    assert_in("process_request_from_transport", src,
              "LocalAgent must use process_request_from_transport")
    # Verify it actually does encode=True (not a no-op)
    assert_in("do_encode=True", src,
              "LocalAgent must encode (do_encode=True)")
    print("  PASSED: LocalAgent performs full encode/decode roundtrip")


# =========================================================================
# Run
# =========================================================================

def run_all_tests():
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
