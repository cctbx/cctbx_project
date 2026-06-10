"""
Sandbox tests for v118 Section F: BUILD experiment_type threading +
R-free auto-fill fix.

Section F fixes a cycle-1 BUILD bug surfaced in production run
AIAgent_245 (OpenAI + "refine model and stop"): refine command
emitted without `generate_rfree_flags=True`, causing phenix.refine
to fail on missing R-free flags.

Root cause: experiment_type is locked on the session only AFTER the
first program returns history_record.  For cycle 1 before any
program runs, session.get_experiment_type() returns "".  This empty
value reached BUILD via CommandContext.from_state, and the R-free
auto-fill in _apply_invariants was guarded by
`experiment_type == "xray"`, which failed.

Two fixes:

  F1 — CommandContext.from_state falls back to
       workflow_state["experiment_type"] when session_info's is
       empty.  Mirrors the existing resolution fallback pattern.

  F2 — _apply_invariants removes the redundant
       `experiment_type == "xray"` guard from the phenix.refine
       R-free block.  phenix.refine is intrinsically xray; the
       program name guard is sufficient.

Tests K_F1-K_F3 cover F1.
Tests K_F4-K_F8 cover F2.
K_F9 is an integration test invoking the full pipeline.
"""

import os
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

COMMAND_BUILDER_PATH = os.path.join(ROOT, "agent", "command_builder.py")


# --------------------------------------------------------------------
# Stub setup for command_builder.py imports
# --------------------------------------------------------------------
# command_builder.py has phenix-side dependencies (file_utils,
# pattern_manager, program_registry) that aren't available in the
# sandbox.  Stub them out before import.
def _install_stubs():
    """Install minimal stub modules so command_builder.py can import."""
    if "libtbx.langchain.agent.file_utils" in sys.modules:
        return  # already installed

    # libtbx.langchain.agent.file_utils — provides matches_exclude_pattern
    stub_file_utils = types.ModuleType(
        "libtbx.langchain.agent.file_utils")
    stub_file_utils.matches_exclude_pattern = (
        lambda *a, **kw: False)

    # libtbx.langchain.agent.pattern_manager — extract_*
    stub_pattern_manager = types.ModuleType(
        "libtbx.langchain.agent.pattern_manager")
    stub_pattern_manager.extract_cycle_number = (
        lambda *a, **kw: None)
    stub_pattern_manager.extract_all_numbers = (
        lambda *a, **kw: [])

    # libtbx.langchain.agent.program_registry — ProgramRegistry
    stub_program_registry = types.ModuleType(
        "libtbx.langchain.agent.program_registry")

    class _StubProgramRegistry(object):
        def __init__(self):
            pass
        def get_invariants(self, prog):
            return []
        def get_program(self, prog):
            return None

    stub_program_registry.ProgramRegistry = _StubProgramRegistry

    # Parent module stubs
    sys.modules.setdefault("libtbx", types.ModuleType("libtbx"))
    sys.modules.setdefault(
        "libtbx.langchain",
        types.ModuleType("libtbx.langchain"))
    sys.modules.setdefault(
        "libtbx.langchain.agent",
        types.ModuleType("libtbx.langchain.agent"))

    sys.modules[
        "libtbx.langchain.agent.file_utils"] = stub_file_utils
    sys.modules[
        "libtbx.langchain.agent.pattern_manager"] = stub_pattern_manager
    sys.modules[
        "libtbx.langchain.agent.program_registry"] = stub_program_registry


_install_stubs()


# --------------------------------------------------------------------
# Module import helper
# --------------------------------------------------------------------

def _import_command_builder():
    """Import command_builder.py module under test."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "command_builder_under_test", COMMAND_BUILDER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------
# Site 1 (F1): CommandContext.from_state
# --------------------------------------------------------------------

def test_k_f1_from_state_falls_back_to_workflow_state():
    """K_F1: from_state falls back to workflow_state.experiment_type
    when session_info.experiment_type is empty.

    This is the AIAgent_245 cycle-1 condition: session hasn't locked
    experiment_type yet (no prior program returned), but the workflow
    engine has already inferred "xray" from input file extensions.
    """
    cb = _import_command_builder()

    state = {
        "cycle_number": 1,
        "session_info": {
            "experiment_type": "",       # empty — not yet locked
            "best_files": {},
            "rfree_mtz": None,
        },
        "workflow_state": {
            "state": "xray_initial",
            "experiment_type": "xray",   # inferred by workflow engine
        },
    }

    context = cb.CommandContext.from_state(state)

    assert context.experiment_type == "xray", (
        "K_F1 failed: should fall back to workflow_state value 'xray', "
        "got %r" % context.experiment_type)

    print("  PASS: K_F1 (from_state falls back to workflow_state)")


def test_k_f2_from_state_prefers_session_when_set():
    """K_F2: from_state prefers session_info.experiment_type when set.

    Backward compatibility: once the session locks an experiment_type,
    that value is authoritative.  Workflow_state is only the fallback.
    """
    cb = _import_command_builder()

    state = {
        "cycle_number": 2,
        "session_info": {
            "experiment_type": "xray",   # session-locked
            "best_files": {},
            "rfree_mtz": None,
        },
        "workflow_state": {
            "state": "xray_initial",
            "experiment_type": "cryoem",  # hypothetical conflict
        },
    }

    context = cb.CommandContext.from_state(state)

    assert context.experiment_type == "xray", (
        "K_F2 failed: should prefer session_info value 'xray', got %r"
        % context.experiment_type)

    print("  PASS: K_F2 (from_state prefers session over workflow_state)")


def test_k_f3_from_state_returns_empty_when_both_absent():
    """K_F3: from_state returns empty string when both are absent.

    No regression: if neither source has experiment_type, the default
    "" is preserved.  Downstream code should handle this case.
    """
    cb = _import_command_builder()

    state = {
        "cycle_number": 1,
        "session_info": {
            "best_files": {},
            "rfree_mtz": None,
        },
        "workflow_state": {
            "state": "unknown",
        },
    }

    context = cb.CommandContext.from_state(state)

    assert context.experiment_type == "", (
        "K_F3 failed: should return '' when both absent, got %r"
        % context.experiment_type)

    print("  PASS: K_F3 (from_state returns '' when both absent)")


# --------------------------------------------------------------------
# Site 2 (F2): R-free auto-fill in _apply_invariants
# --------------------------------------------------------------------

class _StubRegistry(object):
    """Minimal registry stub for testing _apply_invariants in isolation.

    _apply_invariants iterates over invariants from the registry (the
    data-driven invariants list).  We return an empty list so the
    test focuses on the phenix.refine-specific R-free block at the
    end of _apply_invariants.
    """
    def get_invariants(self, program):
        return []

    def get_program(self, program):
        return None


def _make_builder_with_stub_registry():
    """Build a CommandBuilder with a stub registry for unit testing."""
    cb = _import_command_builder()
    builder = cb.CommandBuilder()
    builder._registry = _StubRegistry()
    return builder, cb


def _make_context(cb, experiment_type="", rfree_mtz=None):
    """Build a CommandContext for unit testing _apply_invariants."""
    return cb.CommandContext(
        experiment_type=experiment_type,
        rfree_mtz=rfree_mtz,
    )


def test_k_f4_refine_empty_experiment_type_no_rfree_mtz_auto_fills():
    """K_F4: phenix.refine + empty experiment_type + no rfree_mtz
    → generate_rfree_flags=True is auto-filled.

    This is the AIAgent_245 reproducer.  Pre-F2, the empty
    experiment_type would have caused the auto-fill block to be
    skipped, and the command would have been emitted without
    R-free generation.
    """
    builder, cb = _make_builder_with_stub_registry()
    context = _make_context(cb, experiment_type="", rfree_mtz=None)

    files = {"model": "test.pdb", "data_mtz": "test.mtz"}
    strategy = {}

    files_out, strategy_out = builder._apply_invariants(
        "phenix.refine", files, strategy, context)

    assert strategy_out.get("generate_rfree_flags") is True, (
        "K_F4 failed: R-free auto-fill should fire with empty "
        "experiment_type and no rfree_mtz. strategy=%r" % strategy_out)

    print("  PASS: K_F4 (refine + empty exp_type + no rfree_mtz "
          "→ auto-filled)")


def test_k_f5_refine_xray_experiment_type_auto_fills():
    """K_F5: phenix.refine + experiment_type=xray + no rfree_mtz
    → generate_rfree_flags=True (unchanged behavior).
    """
    builder, cb = _make_builder_with_stub_registry()
    context = _make_context(cb, experiment_type="xray", rfree_mtz=None)

    strategy = {}
    files_out, strategy_out = builder._apply_invariants(
        "phenix.refine", {}, strategy, context)

    assert strategy_out.get("generate_rfree_flags") is True, (
        "K_F5 failed: R-free auto-fill should fire with xray "
        "experiment_type. strategy=%r" % strategy_out)

    print("  PASS: K_F5 (refine + xray exp_type → auto-filled)")


def test_k_f6_refine_with_rfree_mtz_locked_strips_generate():
    """K_F6: phenix.refine + rfree_mtz locked → strip generate.

    The else-branch must still fire correctly: when rfree_mtz is
    locked (from a prior successful refine), any existing
    generate_rfree_flags should be removed.
    """
    builder, cb = _make_builder_with_stub_registry()
    context = _make_context(
        cb, experiment_type="xray",
        rfree_mtz="/path/to/locked.mtz")

    strategy = {"generate_rfree_flags": True}  # LLM tried to set it
    files_out, strategy_out = builder._apply_invariants(
        "phenix.refine", {}, strategy, context)

    assert "generate_rfree_flags" not in strategy_out, (
        "K_F6 failed: generate_rfree_flags should be stripped when "
        "rfree_mtz is locked. strategy=%r" % strategy_out)

    print("  PASS: K_F6 (refine + rfree_mtz locked → strip generate)")


def test_k_f7_refine_with_existing_generate_in_strategy_idempotent():
    """K_F7: phenix.refine with generate_rfree_flags=True already
    in strategy → not overwritten (idempotent).

    If the LLM (or earlier code) already set the flag, F2 should not
    overwrite it.  This case fires the `if "generate_rfree_flags" not
    in strategy:` check, which returns False, so no action is taken.
    """
    builder, cb = _make_builder_with_stub_registry()
    context = _make_context(cb, experiment_type="", rfree_mtz=None)

    strategy = {"generate_rfree_flags": True}  # already set
    files_out, strategy_out = builder._apply_invariants(
        "phenix.refine", {}, strategy, context)

    assert strategy_out.get("generate_rfree_flags") is True, (
        "K_F7 failed: existing True should be preserved, got %r"
        % strategy_out.get("generate_rfree_flags"))

    # Also: don't unexpectedly add other keys
    assert set(strategy_out.keys()) == {"generate_rfree_flags"}, (
        "K_F7 failed: should not add other keys, got %r"
        % list(strategy_out.keys()))

    print("  PASS: K_F7 (refine with existing generate → idempotent)")


def test_k_f8_non_refine_program_block_does_not_fire():
    """K_F8: non-refine programs don't trigger the R-free block.

    Verify that the program-name guard (now the only guard, post-F2)
    correctly excludes other programs from this block.
    """
    builder, cb = _make_builder_with_stub_registry()
    context = _make_context(cb, experiment_type="xray", rfree_mtz=None)

    for prog in ["phenix.xtriage", "phenix.phaser",
                 "phenix.autosol", "phenix.ligandfit",
                 "phenix.polder", "phenix.real_space_refine"]:
        strategy = {}
        files_out, strategy_out = builder._apply_invariants(
            prog, {}, strategy, context)

        assert "generate_rfree_flags" not in strategy_out, (
            "K_F8 failed: R-free block should not fire for %s, "
            "strategy=%r" % (prog, strategy_out))

    print("  PASS: K_F8 (non-refine programs → R-free block skipped)")


# --------------------------------------------------------------------
# Integration test (K_F9)
# --------------------------------------------------------------------

def test_k_f9_aiagent_245_integration_reproducer():
    """K_F9: end-to-end integration test reproducing AIAgent_245.

    Simulates the full state dict that BUILD would see for cycle 1
    of an OpenAI "refine and stop" run:
      - session_info: empty experiment_type, no rfree_mtz, no
        best_files
      - workflow_state: experiment_type="xray" (inferred), state=
        "xray_initial"
      - directives: stop_conditions only (typical for bare advice)
      - cycle_number: 1
      - strategy: empty (planner LLM didn't add generate_rfree_flags
        for bare advice)

    Verifies that the resulting CommandContext + _apply_invariants
    produces a strategy with generate_rfree_flags=True.

    Note: this is a unit-level integration that exercises the patched
    code paths (CommandContext.from_state + _apply_invariants) end-to-
    end.  It does NOT invoke CommandBuilder.build() because that
    requires the full data-driven registry and file selection
    pipeline, which is heavier than needed to verify F1+F2.
    """
    builder, cb = _make_builder_with_stub_registry()

    # State dict mirroring AIAgent_245 cycle 1
    state = {
        "cycle_number": 1,
        "session_info": {
            "experiment_type": "",       # not yet locked (cycle 1)
            "best_files": {},
            "rfree_mtz": None,
            "rfree_resolution": None,
            "recovery_strategies": {},
        },
        "workflow_state": {
            "state": "xray_initial",
            "experiment_type": "xray",   # inferred from input files
            "categorized_files": {
                "data_mtz": ["/path/nsf-d2.mtz"],
                "model": ["/path/nsf-d2_noligand.pdb"],
            },
        },
        "directives": {
            "stop_conditions": {
                "after_program": "phenix.refine",
                "stop_after_requested": True,
                "skip_validation": True,
            }
        },
        "history": [],
        "strategy": {},  # empty — bare advice, planner didn't hint
    }

    # F1: build context from state
    context = cb.CommandContext.from_state(state)

    # F1 verification
    assert context.experiment_type == "xray", (
        "K_F9 failed (F1): expected experiment_type='xray' via "
        "workflow_state fallback, got %r" % context.experiment_type)
    assert context.rfree_mtz is None, (
        "K_F9 failed: rfree_mtz should be None for cycle 1, got %r"
        % context.rfree_mtz)

    # F2: simulate _apply_invariants on phenix.refine
    files = {"model": "/path/nsf-d2_noligand.pdb",
             "data_mtz": "/path/nsf-d2.mtz"}
    strategy = {}  # empty strategy (bare advice case)

    files_out, strategy_out = builder._apply_invariants(
        "phenix.refine", files, strategy, context)

    # F2 verification: R-free auto-fill must have fired
    assert strategy_out.get("generate_rfree_flags") is True, (
        "K_F9 failed (F2): expected generate_rfree_flags=True from "
        "auto-fill, got strategy=%r" % strategy_out)

    print("  PASS: K_F9 (AIAgent_245 end-to-end reproducer)")


# --------------------------------------------------------------------
# Test runner
# --------------------------------------------------------------------

def run_all_tests():
    tests = [
        # F1: from_state experiment_type threading
        ("K_F1_from_state_falls_back_to_workflow_state",
         test_k_f1_from_state_falls_back_to_workflow_state),
        ("K_F2_from_state_prefers_session_when_set",
         test_k_f2_from_state_prefers_session_when_set),
        ("K_F3_from_state_returns_empty_when_both_absent",
         test_k_f3_from_state_returns_empty_when_both_absent),
        # F2: R-free auto-fill
        ("K_F4_refine_empty_experiment_type_no_rfree_mtz_auto_fills",
         test_k_f4_refine_empty_experiment_type_no_rfree_mtz_auto_fills),
        ("K_F5_refine_xray_experiment_type_auto_fills",
         test_k_f5_refine_xray_experiment_type_auto_fills),
        ("K_F6_refine_with_rfree_mtz_locked_strips_generate",
         test_k_f6_refine_with_rfree_mtz_locked_strips_generate),
        ("K_F7_refine_with_existing_generate_in_strategy_idempotent",
         test_k_f7_refine_with_existing_generate_in_strategy_idempotent),
        ("K_F8_non_refine_program_block_does_not_fire",
         test_k_f8_non_refine_program_block_does_not_fire),
        # Integration
        ("K_F9_aiagent_245_integration_reproducer",
         test_k_f9_aiagent_245_integration_reproducer),
    ]
    passed = 0
    failed = 0
    for name, fn in tests:
        print("Test: %s" % name)
        try:
            fn()
            passed += 1
        except AssertionError as e:
            print("  FAIL: %s" % e)
            failed += 1
        except Exception as e:
            import traceback
            print("  ERROR: %s" % e)
            print(traceback.format_exc())
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
