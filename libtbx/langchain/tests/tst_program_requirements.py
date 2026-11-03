"""Tests for v116.10 Tier 2.1 — declarative `requirements:` schema.

The Tier 2.1 cycle added an optional `requirements:` block to
entries in `programs.yaml`.  When present, the block is evaluated
by `WorkflowEngine._check_requirements()` and applied as an
additional filter pass inside `_filter_programs_missing_data_inputs`.

The schema is boolean-only (no file I/O), uses the existing
condition vocabulary (`has`/`has_any`/`not_has`/`not_done`) plus
two new keywords (`done`/`not_done`), and is fully backward-
compatible — programs without a `requirements:` block are
unaffected by the new code path.

This test file covers:

  A. Parser unit tests
     - Each clause type evaluated in isolation
     - Unknown clause keyword handling (soft warning + skip)
     - Malformed input (defensive: return True)
     - Empty / missing `requires:` list

  B. Filter integration tests
     - autobuild scenarios (the initial declaration's target):
       * phased MTZ alone               → kept
       * data MTZ + model               → kept
       * data MTZ + placed model        → kept
       * data MTZ alone (no model)      → filtered
       * model alone (no data)          → filtered
       * empty session                  → filtered

  C. Backward compatibility
     - A program without `requirements:` is unaffected
     - The existing xtriage/mtriage filters still work
     - Order of operations: requirements pass runs AFTER existing
       xtriage/mtriage checks

These tests use the libtbx-stub pattern so they run without a real
PHENIX install.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import types


# =============================================================================
# Path setup
# =============================================================================

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.normpath(os.path.join(_HERE, ".."))

for d in (_PROJECT_ROOT, os.path.join(_PROJECT_ROOT, "agent")):
    if d not in sys.path:
        sys.path.insert(0, d)


# =============================================================================
# libtbx stubs
# =============================================================================

# The current program registry used by the test fixtures.  Tests
# override this to install per-test program definitions.
_TEST_PROGRAM_REGISTRY = {}


def _install_stubs():
    class _Stub:
        def __init__(self):
            pass

    def _ensure(modname):
        if modname not in sys.modules:
            sys.modules[modname] = types.ModuleType(modname)
        return sys.modules[modname]

    _ensure("libtbx")
    _ensure("libtbx.langchain")
    _ensure("libtbx.langchain.agent")
    _ensure("libtbx.langchain.knowledge")

    yl = _ensure("libtbx.langchain.knowledge.yaml_loader")

    yl.get_workflow_steps = lambda et: {
        "analyze": {"description": "Analyze data"},
        "obtain_model": {"description": "Obtain model"},
        "molecular_replacement": {"description": "MR"},
        "experimental_phasing": {"description": "Experimental phasing"},
        "build_from_phases": {"description": "Build from phases"},
        "refine": {"description": "Refine"},
        "combine_ligand": {"description": "Combine ligand"},
        "validate": {"description": "Validate"},
        "complete": {"description": "Complete", "stop": True},
    }
    yl.get_workflow_targets = lambda et, m: None
    yl.get_metric_threshold = lambda *args, **kwargs: None
    yl.get_program = lambda p: _TEST_PROGRAM_REGISTRY.get(p)

    pr = _ensure("libtbx.langchain.agent.program_registry")
    pr.ProgramRegistry = _Stub


try:
    from libtbx.langchain.agent.workflow_engine import WorkflowEngine
except ImportError:
    _install_stubs()
    try:
        from libtbx.langchain.agent.workflow_engine import WorkflowEngine
    except ImportError:
        from workflow_engine import WorkflowEngine


def _engine():
    return WorkflowEngine()


def _set_registry(programs):
    """Install a per-test program registry.  Tests pass a dict of
    program_name → program_def (where program_def may include a
    `requirements:` block).
    """
    global _TEST_PROGRAM_REGISTRY
    _TEST_PROGRAM_REGISTRY = dict(programs)


def _clear_registry():
    global _TEST_PROGRAM_REGISTRY
    _TEST_PROGRAM_REGISTRY = {}


# =============================================================================
# A. Parser unit tests — _check_requirements()
# =============================================================================

def test_has_clause_passes_when_flag_set():
    print("Test: has_clause_passes_when_flag_set")
    eng = _engine()
    block = {"requires": [{"has": "model"}]}
    ctx = {"has_model": True}
    assert eng._check_requirements(block, ctx) is True
    print("  PASS")


def test_has_clause_fails_when_flag_missing():
    print("Test: has_clause_fails_when_flag_missing")
    eng = _engine()
    block = {"requires": [{"has": "model"}]}
    ctx = {}
    assert eng._check_requirements(block, ctx) is False
    print("  PASS")


def test_has_clause_fails_when_flag_explicitly_false():
    print("Test: has_clause_fails_when_flag_explicitly_false")
    eng = _engine()
    block = {"requires": [{"has": "model"}]}
    ctx = {"has_model": False}
    assert eng._check_requirements(block, ctx) is False
    print("  PASS")


def test_has_any_passes_when_one_set():
    print("Test: has_any_passes_when_one_set")
    eng = _engine()
    block = {"requires": [{"has_any": ["data_mtz", "phased_data_mtz"]}]}
    ctx = {"has_data_mtz": True}
    assert eng._check_requirements(block, ctx) is True
    print("  PASS")


def test_has_any_passes_when_second_set():
    print("Test: has_any_passes_when_second_set")
    eng = _engine()
    block = {"requires": [{"has_any": ["data_mtz", "phased_data_mtz"]}]}
    ctx = {"has_phased_data_mtz": True}
    assert eng._check_requirements(block, ctx) is True
    print("  PASS")


def test_has_any_fails_when_none_set():
    print("Test: has_any_fails_when_none_set")
    eng = _engine()
    block = {"requires": [{"has_any": ["data_mtz", "phased_data_mtz"]}]}
    ctx = {}
    assert eng._check_requirements(block, ctx) is False
    print("  PASS")


def test_not_has_passes_when_flag_absent():
    print("Test: not_has_passes_when_flag_absent")
    eng = _engine()
    block = {"requires": [{"not_has": "model"}]}
    ctx = {}
    assert eng._check_requirements(block, ctx) is True
    print("  PASS")


def test_not_has_fails_when_flag_set():
    print("Test: not_has_fails_when_flag_set")
    eng = _engine()
    block = {"requires": [{"not_has": "model"}]}
    ctx = {"has_model": True}
    assert eng._check_requirements(block, ctx) is False
    print("  PASS")


def test_done_passes_when_flag_set():
    print("Test: done_passes_when_flag_set")
    eng = _engine()
    block = {"requires": [{"done": "xtriage"}]}
    ctx = {"xtriage_done": True}
    assert eng._check_requirements(block, ctx) is True
    print("  PASS")


def test_done_fails_when_flag_absent():
    print("Test: done_fails_when_flag_absent")
    eng = _engine()
    block = {"requires": [{"done": "xtriage"}]}
    ctx = {}
    assert eng._check_requirements(block, ctx) is False
    print("  PASS")


def test_not_done_passes_when_flag_absent():
    print("Test: not_done_passes_when_flag_absent")
    eng = _engine()
    block = {"requires": [{"not_done": "autobuild"}]}
    ctx = {}
    assert eng._check_requirements(block, ctx) is True
    print("  PASS")


def test_not_done_fails_when_flag_set():
    print("Test: not_done_fails_when_flag_set")
    eng = _engine()
    block = {"requires": [{"not_done": "autobuild"}]}
    ctx = {"autobuild_done": True}
    assert eng._check_requirements(block, ctx) is False
    print("  PASS")


def test_and_of_clauses_all_must_pass():
    print("Test: and_of_clauses_all_must_pass")
    eng = _engine()
    block = {"requires": [
        {"has": "model"},
        {"has_any": ["data_mtz", "phased_data_mtz"]},
    ]}
    # Both pass
    ctx = {"has_model": True, "has_data_mtz": True}
    assert eng._check_requirements(block, ctx) is True
    # First fails, second passes — overall fail
    ctx = {"has_data_mtz": True}
    assert eng._check_requirements(block, ctx) is False
    # First passes, second fails — overall fail
    ctx = {"has_model": True}
    assert eng._check_requirements(block, ctx) is False
    print("  PASS")


def test_empty_requires_list_passes():
    print("Test: empty_requires_list_passes")
    eng = _engine()
    block = {"requires": []}
    ctx = {}
    assert eng._check_requirements(block, ctx) is True
    print("  PASS")


def test_missing_requires_key_passes():
    """A `requirements:` block without `requires:` is treated as
    no constraint.  Defensive — never filter on malformed input."""
    print("Test: missing_requires_key_passes")
    eng = _engine()
    block = {}
    ctx = {}
    assert eng._check_requirements(block, ctx) is True
    print("  PASS")


def test_non_dict_block_passes():
    """A non-dict input is treated as no constraint."""
    print("Test: non_dict_block_passes")
    eng = _engine()
    ctx = {}
    assert eng._check_requirements(None, ctx) is True
    assert eng._check_requirements("malformed", ctx) is True
    assert eng._check_requirements([], ctx) is True
    print("  PASS")


def test_unknown_clause_keyword_skipped_with_warning():
    """Unknown keyword: print warning, treat clause as satisfied,
    continue evaluating remaining clauses.  No exception unless
    PHENIX_AGENT_STRICT is set."""
    print("Test: unknown_clause_keyword_skipped_with_warning")
    eng = _engine()
    block = {"requires": [
        {"frobnicate": "model"},   # unknown — skipped
        {"has": "data_mtz"},        # known — must pass
    ]}
    # First clause unknown (skipped → satisfied);
    # second clause requires has_data_mtz which is set → True
    ctx = {"has_data_mtz": True}
    assert eng._check_requirements(block, ctx) is True
    # Without has_data_mtz, the second clause fails → False
    ctx = {}
    assert eng._check_requirements(block, ctx) is False
    print("  PASS")


def test_strict_mode_raises_on_unknown_keyword():
    print("Test: strict_mode_raises_on_unknown_keyword")
    eng = _engine()
    block = {"requires": [{"frobnicate": "model"}]}
    ctx = {}
    os.environ["PHENIX_AGENT_STRICT"] = "1"
    try:
        try:
            eng._check_requirements(block, ctx)
            raised = False
        except ValueError:
            raised = True
        assert raised, "Strict mode should raise ValueError on unknown keyword"
    finally:
        del os.environ["PHENIX_AGENT_STRICT"]
    print("  PASS")


def test_any_of_passes_when_one_sub_clause_passes():
    """any_of with mixed clause types: passes if any sub-clause is true."""
    print("Test: any_of_passes_when_one_sub_clause_passes")
    eng = _engine()
    block = {"requires": [{"any_of": [
        {"has": "model"},
        {"done": "phaser"},
    ]}]}
    # has_model True → passes
    assert eng._check_requirements(block, {"has_model": True}) is True
    # phaser_done True → passes (different flag family)
    assert eng._check_requirements(block, {"phaser_done": True}) is True
    print("  PASS")


def test_any_of_fails_when_no_sub_clauses_pass():
    print("Test: any_of_fails_when_no_sub_clauses_pass")
    eng = _engine()
    block = {"requires": [{"any_of": [
        {"has": "model"},
        {"done": "phaser"},
    ]}]}
    # Neither has_model nor phaser_done → fails
    assert eng._check_requirements(block, {}) is False
    # Unrelated flag set → still fails
    assert eng._check_requirements(
        block, {"has_data_mtz": True}) is False
    print("  PASS")


def test_any_of_with_nested_has_any():
    """any_of sub-clauses can themselves be has_any clauses."""
    print("Test: any_of_with_nested_has_any")
    eng = _engine()
    block = {"requires": [{"any_of": [
        {"has_any": ["model", "placed_model"]},
        {"done": "autosol"},
    ]}]}
    # has_placed_model True → first branch passes (via has_any)
    assert eng._check_requirements(
        block, {"has_placed_model": True}) is True
    # autosol_done True → second branch passes (via done)
    assert eng._check_requirements(
        block, {"autosol_done": True}) is True
    # Neither → fails
    assert eng._check_requirements(block, {}) is False
    print("  PASS")


def test_any_of_with_empty_sub_clauses_fails():
    """An any_of with no sub-clauses can satisfy nothing, so fails."""
    print("Test: any_of_with_empty_sub_clauses_fails")
    eng = _engine()
    block = {"requires": [{"any_of": []}]}
    # No sub-clauses → no way to satisfy → fail
    assert eng._check_requirements(block, {"has_anything": True}) is False
    print("  PASS")


def test_any_of_with_non_list_skipped_defensively():
    """Malformed any_of (non-list) is skipped, not treated as fail."""
    print("Test: any_of_with_non_list_skipped_defensively")
    eng = _engine()
    block = {"requires": [{"any_of": "not a list"}]}
    # Malformed → defensively passes (no fail-closed on bad YAML)
    assert eng._check_requirements(block, {}) is True
    print("  PASS")


def test_and_of_has_any_plus_any_of():
    """Top-level AND of a has_any and an any_of: both must pass."""
    print("Test: and_of_has_any_plus_any_of")
    eng = _engine()
    block = {"requires": [
        {"has_any": ["data_mtz", "phased_data_mtz"]},
        {"any_of": [
            {"has": "model"},
            {"done": "autosol"},
        ]},
    ]}
    # data + model: both pass
    assert eng._check_requirements(
        block, {"has_data_mtz": True, "has_model": True}) is True
    # data + autosol_done: both pass (model not present but autosol_done is)
    assert eng._check_requirements(
        block, {"has_data_mtz": True, "autosol_done": True}) is True
    # data alone: first passes, second fails
    assert eng._check_requirements(
        block, {"has_data_mtz": True}) is False
    # model alone: first fails
    assert eng._check_requirements(
        block, {"has_model": True}) is False
    print("  PASS")


# =============================================================================
# B. Filter integration tests — _filter_programs_missing_data_inputs
#    with autobuild's declaration
# =============================================================================

# The autobuild requirements: block as it appears in programs.yaml.
# Tests install this into _TEST_PROGRAM_REGISTRY.
#
# v116.10 Tier 2.1 update: second clause uses any_of to accept either
# file-based (has_*) OR history-based (_done) flags.  This matches the
# existing explanation pattern in explain_unavailable_program.
AUTOBUILD_DEF = {
    "experiment_types": ["xray"],
    "requirements": {
        "requires": [
            {"has_any": ["data_mtz", "phased_data_mtz"]},
            {"any_of": [
                {"has": "phased_data_mtz"},
                {"has": "model"},
                {"has": "placed_model"},
                {"has": "placed_model_from_history"},
                {"done": "phaser"},
                {"done": "autosol"},
            ]},
        ]
    }
}


def _install_autobuild_only():
    _set_registry({"phenix.autobuild": AUTOBUILD_DEF})


def test_autobuild_kept_with_phased_data_alone():
    """A phased MTZ satisfies both clauses (it has both data and
    phase information embedded), so autobuild should be kept even
    without a separate model."""
    print("Test: autobuild_kept_with_phased_data_alone")
    _install_autobuild_only()
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild"],
            {"has_phased_data_mtz": True},
        )
        assert result == ["phenix.autobuild"], (
            "Expected autobuild kept; got %r" % result)
    finally:
        _clear_registry()
    print("  PASS")


def test_autobuild_kept_with_data_and_model():
    print("Test: autobuild_kept_with_data_and_model")
    _install_autobuild_only()
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild"],
            {"has_data_mtz": True, "has_model": True},
        )
        assert result == ["phenix.autobuild"]
    finally:
        _clear_registry()
    print("  PASS")


def test_autobuild_kept_with_data_and_placed_model():
    """After phaser MR, has_placed_model is set; autobuild should
    accept this as the model side of the requirement."""
    print("Test: autobuild_kept_with_data_and_placed_model")
    _install_autobuild_only()
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild"],
            {"has_data_mtz": True, "has_placed_model": True},
        )
        assert result == ["phenix.autobuild"]
    finally:
        _clear_registry()
    print("  PASS")


def test_autobuild_filtered_with_data_alone_no_model_no_phases():
    """Raw MTZ data without model or phases: autobuild cannot run.
    This is the original bug case — autobuild would crash at runtime
    with 'MTZ lacks phase/FOM columns'."""
    print("Test: autobuild_filtered_with_data_alone_no_model_no_phases")
    _install_autobuild_only()
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild"],
            {"has_data_mtz": True},  # data only — no model, no phases
        )
        assert result == [], (
            "Expected autobuild filtered; got %r" % result)
    finally:
        _clear_registry()
    print("  PASS")


def test_autobuild_filtered_with_model_alone_no_data():
    """A model without diffraction data is not enough for autobuild."""
    print("Test: autobuild_filtered_with_model_alone_no_data")
    _install_autobuild_only()
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild"],
            {"has_model": True},
        )
        assert result == []
    finally:
        _clear_registry()
    print("  PASS")


def test_autobuild_filtered_with_empty_session():
    print("Test: autobuild_filtered_with_empty_session")
    _install_autobuild_only()
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild"],
            {},
        )
        assert result == []
    finally:
        _clear_registry()
    print("  PASS")


def test_autobuild_kept_after_autosol_done():
    """S5A regression: autosol succeeded → autobuild must be available.

    The simulator (and some real session paths) set autosol_done=True
    after autosol completes, but may not set has_phased_data_mtz on
    the same cycle.  The rule must accept autosol_done as satisfying
    the phase-information requirement.
    """
    print("Test: autobuild_kept_after_autosol_done")
    _install_autobuild_only()
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild"],
            {"has_data_mtz": True, "autosol_done": True},
        )
        assert result == ["phenix.autobuild"], (
            "S5A: autobuild must survive when autosol_done=True; "
            "got %r" % result)
    finally:
        _clear_registry()
    print("  PASS")


def test_autobuild_kept_after_phaser_done():
    """Analogous to S5A but for the MR path: phaser_done satisfies
    the phase/model requirement."""
    print("Test: autobuild_kept_after_phaser_done")
    _install_autobuild_only()
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild"],
            {"has_data_mtz": True, "phaser_done": True},
        )
        assert result == ["phenix.autobuild"]
    finally:
        _clear_registry()
    print("  PASS")


def test_autobuild_kept_with_placed_model_from_history():
    """A model placed in earlier session history satisfies the model
    requirement.  Matches existing explanation pattern."""
    print("Test: autobuild_kept_with_placed_model_from_history")
    _install_autobuild_only()
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild"],
            {"has_data_mtz": True, "has_placed_model_from_history": True},
        )
        assert result == ["phenix.autobuild"]
    finally:
        _clear_registry()
    print("  PASS")


def test_autobuild_still_filtered_when_only_data_no_history():
    """Negative case unchanged: raw data alone, no phase info, no
    model, no completed phasing program — autobuild filtered."""
    print("Test: autobuild_still_filtered_when_only_data_no_history")
    _install_autobuild_only()
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild"],
            {"has_data_mtz": True},
        )
        assert result == [], (
            "Raw data alone (no model, no phases, no done flags) "
            "should still filter autobuild; got %r" % result)
    finally:
        _clear_registry()
    print("  PASS")


# =============================================================================
# C. Backward compatibility — programs WITHOUT requirements: blocks
# =============================================================================

def test_program_without_requirements_is_unaffected():
    """Foundational backward-compat test.  A program with no
    `requirements:` block in its programs.yaml entry must not be
    touched by the new filter pass."""
    print("Test: program_without_requirements_is_unaffected")
    # Install a program that has NO requirements: block
    _set_registry({
        "phenix.some_program": {
            "experiment_types": ["xray"],
            # No requirements: block here
        }
    })
    try:
        eng = _engine()
        # Empty context — would fail any requirements check that existed
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.some_program"], {})
        assert result == ["phenix.some_program"], (
            "Program without requirements: should pass; got %r" % result)
    finally:
        _clear_registry()
    print("  PASS")


def test_program_not_in_registry_is_unaffected():
    """If get_program returns None (program not in programs.yaml),
    the program is kept as-is.  Behavior should match a program
    that's known but has no `requirements:` block."""
    print("Test: program_not_in_registry_is_unaffected")
    _set_registry({})  # empty registry
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.unknown"], {})
        assert result == ["phenix.unknown"]
    finally:
        _clear_registry()
    print("  PASS")


def test_existing_xtriage_filter_still_works():
    """The existing Phase 4b filter for xtriage runs BEFORE the new
    declarative pass.  Confirm it still operates regardless of any
    requirements: block."""
    print("Test: existing_xtriage_filter_still_works")
    _set_registry({})  # no requirements: blocks installed
    try:
        eng = _engine()
        # No mtz data — xtriage should still be filtered by the
        # existing Phase 4b check
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.xtriage"], {})
        assert result == [], (
            "xtriage should still be filtered by Phase 4b; got %r"
            % result)
    finally:
        _clear_registry()
    print("  PASS")


def test_existing_mtriage_filter_still_works():
    print("Test: existing_mtriage_filter_still_works")
    _set_registry({})
    try:
        eng = _engine()
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.mtriage"], {})
        assert result == []
    finally:
        _clear_registry()
    print("  PASS")


def test_order_of_operations_xtriage_before_requirements():
    """If both the existing xtriage filter AND a requirements:
    block could remove xtriage, the existing filter removes it
    first.  This verifies the order documented in PATH_Y_DESIGN.md.

    Note: xtriage is unlikely to ever get a requirements: block in
    practice, but this test exercises the order regardless.
    """
    print("Test: order_of_operations_xtriage_before_requirements")
    _set_registry({
        "phenix.xtriage": {
            "experiment_types": ["xray"],
            "requirements": {
                "requires": [{"has": "data_mtz"}]
            }
        }
    })
    try:
        eng = _engine()
        # No data_mtz: existing filter will remove first
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.xtriage"], {})
        assert result == []
    finally:
        _clear_registry()
    print("  PASS")


def test_mixed_program_list_with_and_without_requirements():
    """A mixed list of programs (some with requirements:, some
    without) should be filtered correctly per program."""
    print("Test: mixed_program_list_with_and_without_requirements")
    _set_registry({
        "phenix.autobuild": AUTOBUILD_DEF,
        "phenix.some_program": {"experiment_types": ["xray"]},
    })
    try:
        eng = _engine()

        # Session has neither data nor model: autobuild filtered;
        # some_program survives (no requirements)
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild", "phenix.some_program"], {})
        assert "phenix.autobuild" not in result
        assert "phenix.some_program" in result

        # Phased data: autobuild kept; some_program still kept
        result = eng._filter_programs_missing_data_inputs(
            ["phenix.autobuild", "phenix.some_program"],
            {"has_phased_data_mtz": True})
        assert "phenix.autobuild" in result
        assert "phenix.some_program" in result
    finally:
        _clear_registry()
    print("  PASS")


# =============================================================================
# Test runner (matches the project's tst_*.py convention)
# =============================================================================

def run_all_tests():
    try:
        from libtbx.langchain.tests.tst_utils import (
            run_tests_with_fail_fast)
    except ImportError:
        try:
            from tests.tst_utils import run_tests_with_fail_fast
        except ImportError:
            _standalone_runner()
            return
    run_tests_with_fail_fast()


def _standalone_runner():
    test_fns = [(k, v) for k, v in sorted(globals().items())
                if k.startswith("test_") and callable(v)]
    passed = 0
    failed = 0
    for name, fn in test_fns:
        try:
            fn()
            passed += 1
        except Exception as e:
            print("  FAIL: %s — %s" % (name, e))
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    _standalone_runner()
