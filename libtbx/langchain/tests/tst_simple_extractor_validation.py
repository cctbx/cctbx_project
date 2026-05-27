"""K_H14_1: extract_directives_simple validation closure (v119.H14.1).

H14 added the sentinel + Hermann-Mauguin shape check in
``validate_directives``, but Tom's ollama xtriage tutorial run on
2026-05-26 surfaced a gap: when ollama's LLM fails to return parseable
JSON, ``extract_directives`` falls back to ``extract_directives_simple``,
which has its OWN space_group regex AND returns the dict directly
without going through ``validate_directives``.  The H14 sentinel check
was bypassed.

The bug surfaced as `space_group=Not explicitly mentio` in Tom's
extracted-directives display — the truncated value was captured by
the simple extractor's regex (the truncation comes from the regex's
own ``{1,20}`` quantifier, NOT from any LLM length cap) and reached
the agent unchanged.

H14.1 fixes two things:

  1. Latent VALID_STOP_CONDITIONS gap: ``start_with_program`` was set
     by ``_resolve_after_program`` and consumed downstream by
     ``workflow_engine.py``, ``ai_agent.py``, and ``ai_analysis.py``,
     but wasn't in ``VALID_STOP_CONDITIONS``.  ``validate_directives``
     would log "Unknown stop condition start_with_program" and drop
     it.  Pre-H14.1 this didn't matter because validate_directives
     ran BEFORE ``_apply_workflow_intent_fallback`` (which set the
     key), so the key was added post-validation.  H14.1 makes
     validate_directives idempotent for this key.

  2. The actual bypass: H14.1 adds ``validate_directives(...)`` as
     the final step of ``extract_directives_simple``.  Now both the
     LLM-extracted-JSON path AND the simple-pattern-match path apply
     the same sanity check.

This file pins:

  - Tom's exact xtriage case: simple-extractor + "Not explicitly
    mentio" → space_group dropped, log line emitted
  - Item 1's multi-action+stop case (from rules-only path):
    start_with_program=phenix.phaser PRESERVED through
    validate_directives
  - Other simple-extractor outputs (resolution, max_refine_cycles,
    skip_programs, etc.) survive validation
  - The log parameter to extract_directives_simple is optional
    (backward compat with external callers)
  - validate_directives is idempotent on already-validated dicts
    (extract_directives_simple now ends with it, but external code
    may also call validate_directives on the result; running twice
    must be a no-op)

Total: 12 tests.
"""
from __future__ import absolute_import, division, print_function

import os
import sys


def _try_import_de():
    """Import the directive_extractor module."""
    try:
        from libtbx.langchain.agent import directive_extractor as de
        return de, None
    except ImportError:
        pass
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        if parent not in sys.path:
            sys.path.insert(0, parent)
        from agent import directive_extractor as de
        return de, None
    except ImportError as e:
        return None, str(e)


# Tom's preprocessed advice from the 2026-05-26 ollama xtriage run.
# Exact text from session summary, lines 192-203 of ollama_xtriage.log.
# Original log lines had Markdown trailing double-spaces (hard-line-break
# syntax) on each line; stripped here since they don't affect the
# extractor's regex match for "Space group: Not explicitly mentioned"
# and they trip libtbx.find_clutter.
TOMS_ADVICE = """1. **Input Files Found**: porin.mtz, p9_hires.mtz, sec17.sca
2. **Experiment Type**: X-ray crystallography (analysis only)
3. **Primary Goal**: Run Xtriage to analyze data quality, including checking for twinning, space group identification, ice rings, and Wilson plot evaluation.
4. **Key Parameters**:
   - Resolution limit: 1.7 Å (p9_hires.mtz), 2.3 Å (porin.mtz), 3.3 Å (sec17.sca)
   - Space group: Not explicitly mentioned
5. **Program Parameters**: None
6. **Special Instructions**:
   - Analyze porin.mtz for merohedral twinning.
   - Analyze p9_hires.mtz for high-resolution data quality (minimal unusual features).
   - Analyze sec17.sca for potential ice rings.
7. **Stop Condition**: None"""


# =====================================================================
# §A: Tom's exact ollama xtriage case (Gate C closure)
# =====================================================================

def test_tom_xtriage_space_group_dropped():
    """The exact preprocessed advice from Tom's ollama xtriage tutorial
    run must now produce a directives dict with space_group DROPPED.
    Pre-H14.1: space_group='Not explicitly mentio' (gate C failed).
    Post-H14.1: space_group absent + log line emitted.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    log = []
    result = de.extract_directives_simple(TOMS_ADVICE, log=log.append)

    default = result.get("program_settings", {}).get("default", {})
    assert "space_group" not in default, (
        "Tom's 'Not explicitly mentio' value should be DROPPED. "
        "Got default=%r" % default)

    # And the canonical log line should be emitted
    drop_lines = [
        l for l in log
        if "Dropping invalid space_group" in l
        and "Not explicitly mentio" in l
    ]
    assert drop_lines, (
        "Expected 'Dropping invalid space_group value: ...Not "
        "explicitly mentio...' log line. Got: %r" % log)
    print("  PASS: test_tom_xtriage_space_group_dropped")


def test_tom_xtriage_resolution_preserved():
    """Same case as above — resolution (1.7) must still come through
    intact.  Validates that H14.1 doesn't over-correct."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    result = de.extract_directives_simple(TOMS_ADVICE)
    default = result.get("program_settings", {}).get("default", {})
    assert default.get("resolution") == 1.7, (
        "Resolution 1.7 must be preserved (only space_group "
        "should be dropped). Got default=%r" % default)
    print("  PASS: test_tom_xtriage_resolution_preserved")


# =====================================================================
# §B: Item 1 path through simple extractor (start_with_program
# preservation is critical here — Item 1's fix lives in this path)
# =====================================================================

def test_item1_multi_action_preserves_start_with_program():
    """When extract_directives_simple sees explicit MR + stop, it sets
    start_with_program=phenix.phaser via _resolve_after_program.
    validate_directives must NOT strip this key — it's a valid stop
    condition that downstream consumers in workflow_engine.py and
    ai_agent.py read.

    This is the regression-guard for the latent VALID_STOP_CONDITIONS
    bug that H14.1 fixed.  Pre-H14.1, applying validate_directives at
    the end of extract_directives_simple would have stripped the key,
    breaking Item 1's fix.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    result = de.extract_directives_simple(
        "do molecular replacement then refine and stop")
    sc = result.get("stop_conditions", {})
    assert sc.get("start_with_program") == "phenix.phaser", (
        "Explicit MR + stop must set start_with_program=phenix.phaser. "
        "Got: %r" % sc)
    assert sc.get("after_program") == "phenix.refine", (
        "Multi-action + stop should set after_program=last_action. "
        "Got: %r" % sc)
    print("  PASS: test_item1_multi_action_preserves_start_with_program")


def test_item1_1029B_sad_pattern_unaffected():
    """The 1029B-sad-style pattern (goal phrase + stop, no explicit
    MR keyword) should NOT trigger start_with_program=phaser even with
    validate_directives applied.

    Post-H14, "solve the structure" was removed from
    _ACTION_TABLE["solve"]["keywords"], so this advice produces only
    the "refine" action.  Item 1 fix preserved.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    result = de.extract_directives_simple(
        "solve the structure by standard procedures. "
        "anomalous scatterer is selenium. stop after refinement")
    sc = result.get("stop_conditions", {})
    # after_program=phenix.refine (the only detected action) + stop signal
    assert sc.get("after_program") == "phenix.refine", (
        "Goal phrase + 'stop after refinement' should set "
        "after_program=phenix.refine. Got: %r" % sc)
    # And no start_with_program (n==1 with stop branch, not n>1)
    assert "start_with_program" not in sc, (
        "n=1 branch should NOT set start_with_program. Got: %r" % sc)
    print("  PASS: test_item1_1029B_sad_pattern_unaffected")


# =====================================================================
# §C: VALID_STOP_CONDITIONS now includes start_with_program
# =====================================================================

def test_valid_stop_conditions_includes_start_with_program():
    """The latent bug fix: start_with_program must be in
    VALID_STOP_CONDITIONS.  Otherwise validate_directives logs
    'Unknown stop condition start_with_program' and drops the key.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    assert "start_with_program" in de.VALID_STOP_CONDITIONS, (
        "VALID_STOP_CONDITIONS must include start_with_program "
        "(H14.1 latent bug fix).  Got keys: %r"
        % list(de.VALID_STOP_CONDITIONS.keys()))
    assert de.VALID_STOP_CONDITIONS["start_with_program"] is str, (
        "start_with_program should be typed as str. Got: %r"
        % de.VALID_STOP_CONDITIONS["start_with_program"])
    print("  PASS: test_valid_stop_conditions_includes_start_with_program")


def test_validate_directives_preserves_start_with_program():
    """Direct test of the VALID_STOP_CONDITIONS fix: a dict with
    start_with_program in stop_conditions should survive
    validate_directives unchanged."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    directives = {
        "stop_conditions": {
            "after_program": "phenix.refine",
            "start_with_program": "phenix.phaser",
            "skip_validation": True,
            "stop_after_requested": True,
        }
    }
    log = []
    result = de.validate_directives(directives, log=log.append)
    sc = result.get("stop_conditions", {})
    assert sc.get("start_with_program") == "phenix.phaser", (
        "start_with_program must survive validate_directives. "
        "Got: %r" % sc)
    # And NO "Unknown stop condition" log line for start_with_program
    unknown_lines = [
        l for l in log
        if "Unknown stop condition" in l
        and "start_with_program" in l
    ]
    assert not unknown_lines, (
        "validate_directives should NOT log 'Unknown stop "
        "condition start_with_program'.  Got: %r" % unknown_lines)
    print("  PASS: test_validate_directives_preserves_start_with_program")


# =====================================================================
# §D: Other simple-extractor outputs survive validate_directives
# =====================================================================

def test_simple_atom_type_survives_validation():
    """phenix.autosol.atom_type set by simple extractor must survive."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    result = de.extract_directives_simple(
        "solve with phenix.autosol using selenium as the "
        "anomalous scatterer")
    autosol = result.get("program_settings", {}).get(
        "phenix.autosol", {})
    assert autosol.get("atom_type") == "Se", (
        "atom_type='Se' must survive H14.1 validation. "
        "Got autosol=%r" % autosol)
    print("  PASS: test_simple_atom_type_survives_validation")


def test_simple_max_refine_cycles_survives_validation():
    """max_refine_cycles=1 from 'stop after first refine' must survive.

    The simple extractor sets max_refine_cycles=1 only when both a
    stop signal ("stop after") and an explicit count word ("first",
    "one", "1", "single") modify a "refine" verb.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    result = de.extract_directives_simple(
        "stop after the first refine")
    sc = result.get("stop_conditions", {})
    assert sc.get("max_refine_cycles") == 1, (
        "max_refine_cycles=1 must survive H14.1 validation. "
        "Got: %r" % sc)
    print("  PASS: test_simple_max_refine_cycles_survives_validation")


def test_simple_skip_programs_survives_validation():
    """workflow_preferences.skip_programs must survive validation
    when extract_directives_simple sets it."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    result = de.extract_directives_simple("do not run molprobity")
    wp = result.get("workflow_preferences", {})
    skip = wp.get("skip_programs", [])
    assert "phenix.molprobity" in skip, (
        "skip_programs must survive H14.1 validation. "
        "Got wp=%r" % wp)
    print("  PASS: test_simple_skip_programs_survives_validation")


def test_simple_valid_space_group_survives_validation():
    """A real space-group symbol from the advice must survive."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    result = de.extract_directives_simple(
        "space group P 1 21 1, refine the structure")
    default = result.get("program_settings", {}).get("default", {})
    assert default.get("space_group") == "P 1 21 1", (
        "Valid space group 'P 1 21 1' must survive H14.1 validation. "
        "Got default=%r" % default)
    print("  PASS: test_simple_valid_space_group_survives_validation")


# =====================================================================
# §E: Backward compatibility — log param is optional
# =====================================================================

def test_extract_directives_simple_log_param_optional():
    """The four external callers in ai_agent.py and run_ai_analysis.py
    call extract_directives_simple WITHOUT a log argument.  The
    H14.1 signature change must not break them."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    # Should work without log param (backward compat)
    result = de.extract_directives_simple("refine the structure")
    assert isinstance(result, dict), (
        "extract_directives_simple without log param must still "
        "return a dict. Got: %r" % type(result))
    # And explicitly passing None for log should also work
    result2 = de.extract_directives_simple(
        "refine the structure", log=None)
    assert result == result2, (
        "extract_directives_simple(advice) and "
        "extract_directives_simple(advice, log=None) must produce "
        "the same result")
    print("  PASS: test_extract_directives_simple_log_param_optional")


# =====================================================================
# §F: validate_directives idempotency on simple-extractor output
# =====================================================================

def test_validate_directives_idempotent_after_simple():
    """Now that extract_directives_simple ends with validate_directives,
    calling validate_directives a SECOND time on the result must be a
    no-op (external code may also call validate_directives — if it's
    not idempotent we'd see weird interactions).
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    # A nontrivial advice that exercises multiple branches
    advice = (
        "space group P 1 21 1, anomalous scatterer is selenium, "
        "do one refinement and stop")
    result_once = de.extract_directives_simple(advice)
    result_twice = de.validate_directives(result_once)
    assert result_once == result_twice, (
        "validate_directives must be idempotent. "
        "First pass: %r\nSecond pass: %r"
        % (result_once, result_twice))
    print("  PASS: test_validate_directives_idempotent_after_simple")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §A: Tom's gate C closure (2)
    test_tom_xtriage_space_group_dropped()
    test_tom_xtriage_resolution_preserved()
    # §B: Item 1 preservation through validation (2)
    test_item1_multi_action_preserves_start_with_program()
    test_item1_1029B_sad_pattern_unaffected()
    # §C: latent VALID_STOP_CONDITIONS fix (2)
    test_valid_stop_conditions_includes_start_with_program()
    test_validate_directives_preserves_start_with_program()
    # §D: other simple-extractor outputs (4)
    test_simple_atom_type_survives_validation()
    test_simple_max_refine_cycles_survives_validation()
    test_simple_skip_programs_survives_validation()
    test_simple_valid_space_group_survives_validation()
    # §E: backward compat (1)
    test_extract_directives_simple_log_param_optional()
    # §F: idempotency (1)
    test_validate_directives_idempotent_after_simple()


if __name__ == "__main__":
    print("K_H14_1: extract_directives_simple validation closure (v119.H14.1)")
    print("=" * 70)
    run_all_tests()
    print("=" * 70)
    print("K_H14_1 complete.")
