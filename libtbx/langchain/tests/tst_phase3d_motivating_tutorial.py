"""Tests for v116.10 Phase 3d motivating tutorial — dock-and-stop
with sequence + map inputs.

Phase 3d reclassified phenix.dock_in_map and phenix.map_symmetry
as standalone programs.  The motivating case is the "dock and
stop" workflow with sequence + map (no model):

  Pre-Phase-3d behavior:
    1. _initialize_plan_inner hits v116_10_elif
    2. Generates a multi-stage plan (xtriage → predict → dock)
    3. skip_to_program("phenix.dock_in_map") skips the predict
       stage
    4. dock_in_map runs without a predicted model → FAILS at
       runtime ("no model file found")

  Post-Phase-3d behavior:
    1. _initialize_plan_inner hits single_program_skip
    2. No plan generated
    3. workflow_engine takes over via the cryo-EM state machine
    4. State analyzes → routes to obtain_model → predict_and_build
       generates a model → dock_in_map gets the model it needs

This file fills the gap flagged in the v116.10 review: Phase 3d
was verified only by decision-tree traces and unit tests.  No
tutorial currently exercises the dock-and-stop path with
sequence + map inputs.

This file provides:

  1. A `tutorial_expectations.yaml` entry describing what the
     tutorial expects to happen at each step (Tom adds this to
     the tutorial corpus' expectations file).

  2. A small test fixture that asserts the expected control flow
     when run against a real session (the test is skip-aware if
     no session is provided, so it doesn't fail CI).

  3. A README explaining what data files the tutorial needs.

The test integrates with the existing tst_dock_and_stop.py
decision-tree tests, which already verify the classification
logic.  This file verifies the END-TO-END behavior: given
sequence + map inputs and "dock and stop" advice, does the
workflow correctly invoke predict_and_build before dock_in_map?
"""

from __future__ import absolute_import, division, print_function

import os
import json
import sys


# =============================================================================
# Tutorial spec: what data files are needed and what's expected
# =============================================================================

TUTORIAL_NAME = "dock_and_stop_with_sequence_and_map"

TUTORIAL_SPEC = {
    "name": TUTORIAL_NAME,
    "description": (
        "Phase 3d motivating case: user provides sequence + map "
        "and asks 'dock and stop'. Tests that the workflow "
        "correctly chains predict_and_build → dock_in_map "
        "instead of skipping the prediction prerequisite."),
    "inputs": {
        "sequence": "seq.fa",
        "full_map": "map.ccp4",
    },
    "user_advice": "Dock the model into the map and stop",
    "expected_directives": {
        "intent": "task",
        "stop_conditions": {
            "after_program": "phenix.dock_in_map",
        },
    },
    "expected_experiment_type": "cryoem",
    "expected_workflow_path": [
        # Phase 3d post-fix: workflow_engine drives the sequence,
        # generating a model before docking.
        {
            "cycle": 1,
            "program": "phenix.mtriage",
            "rationale": "Analyze map quality first",
        },
        {
            "cycle": 2,
            "program": "phenix.predict_and_build",
            "rationale": (
                "Sequence present, no model — predict model "
                "before docking"),
        },
        {
            "cycle": 3,
            "program": "phenix.dock_in_map",
            "rationale": (
                "Predicted model now available — dock into map "
                "(user-requested stop target)"),
        },
    ],
    "expected_stop_after_cycle": 3,
    "anti_pattern": {
        # The pre-Phase-3d failure mode that this tutorial guards
        # against
        "description": (
            "dock_in_map runs in cycle 1 without a predicted "
            "model, fails with 'no model file found'"),
        "indicator_in_log": "no model file found",
    },
}


def get_expectations_yaml_entry():
    """Return the YAML snippet to add to tutorial_expectations.yaml.

    This is what Tom adds to his tutorial corpus' expectations
    file.  Format mirrors the convention seen in the other
    v116.10 entries in tutorial_expectations.yaml.
    """
    return """
# Phase 3d motivating case (v116.10).
# Verifies that dock-and-stop with sequence + map correctly
# chains predict_and_build → dock_in_map instead of attempting
# dock_in_map directly (the pre-Phase-3d failure mode).
dock_and_stop_with_sequence_and_map:
  description: |
    User provides sequence + map and asks 'dock and stop'.
    Workflow should chain predict → dock, not skip prediction.
  inputs:
    sequence: seq.fa
    full_map: map.ccp4
  user_advice: "Dock the model into the map and stop"
  expected:
    intent: task
    experiment_type: cryoem
    stop_after: phenix.dock_in_map
    programs_run_in_order:
      - phenix.mtriage         # cycle 1: analyze map
      - phenix.predict_and_build  # cycle 2: produce model
      - phenix.dock_in_map     # cycle 3: dock (stop target)
    stop_after_cycle: 3
    must_not_appear_in_log:
      # Pre-Phase-3d failure mode
      - "no model file found"
"""


# =============================================================================
# Decision-tree verification (no session needed)
# =============================================================================

def test_phase_3d_motivating_case_decision_tree():
    """For the motivating-case directives (intent=task,
    stop_after=phenix.dock_in_map), verify the decision tree
    routes to single_program_skip (Phase 3d post-fix), not
    v116_10_elif (Phase 3d pre-fix).

    This is the unit-test layer.  Tutorial-level verification is
    test_phase_3d_motivating_case_against_session (below).
    """
    print("Test: phase_3d_motivating_case_decision_tree")

    # Import the smoke-test infrastructure from the sibling
    # test file
    try:
        from tests.tst_initialize_plan_smoke import (
            _run_decision_tree, _SessionStub, _STANDALONE,
            _require_program_sets)
    except ImportError:
        try:
            from tst_initialize_plan_smoke import (
                _run_decision_tree, _SessionStub, _STANDALONE,
                _require_program_sets)
        except ImportError:
            print("  SKIP — tst_initialize_plan_smoke not importable")
            return

    if not _require_program_sets():
        return

    if "phenix.dock_in_map" not in _STANDALONE:
        print("  SKIP — Phase 3d not deployed")
        return

    session = _SessionStub(directives={
        "intent": "task",
        "stop_conditions": {
            "after_program": "phenix.dock_in_map",
        },
    })
    result = _run_decision_tree(
        session,
        processed_advice="Dock the model into the map and stop")

    assert result["outcome"] == "single_program_skip", (
        "Phase 3d post-fix expects single_program_skip "
        "(workflow_engine drives), got %s" % result["outcome"])
    assert not result["proceeds_to_plan_generation"], (
        "Should not generate a plan — workflow_engine takes over")
    print("  PASS — decision tree routes correctly")


# =============================================================================
# Session-based verification (skip-aware)
# =============================================================================

def _find_session_fixture():
    """Look for a real run's session.json that exercised this
    tutorial.  Returns None if not found (tests skip cleanly).

    Looking in conventional locations the tutorial test harness
    might use; Tom can add a fixture path here.
    """
    candidates = [
        # Test fixture directory
        os.path.join(os.path.dirname(__file__),
                     "fixtures",
                     "dock_and_stop_session.json"),
        # Environment variable override
        os.environ.get("DOCK_AND_STOP_FIXTURE", ""),
    ]
    for c in candidates:
        if c and os.path.exists(c):
            return c
    return None


def _extract_program_history(session_data):
    """Pull the ordered list of program names from session.data.cycles."""
    cycles = session_data.get("cycles", [])
    history = []
    for c in cycles:
        if not isinstance(c, dict):
            continue
        prog = c.get("program") or ""
        result = (c.get("result") or "")
        history.append({
            "program": prog,
            "result": result,
            "success": "SUCCESS" in str(result).upper(),
        })
    return history


def test_phase_3d_motivating_case_against_session():
    """End-to-end: given a session.json from a real run, verify
    the workflow chained predict_and_build → dock_in_map.

    Skips cleanly if no session fixture is available — this test
    is intended to run against a real tutorial result, which Tom
    wires up.
    """
    print("Test: phase_3d_motivating_case_against_session")

    fixture = _find_session_fixture()
    if not fixture:
        print("  SKIP — no session fixture available "
              "(set DOCK_AND_STOP_FIXTURE env var to enable)")
        return

    with open(fixture, encoding='utf-8') as f:
        session = json.load(f)

    history = _extract_program_history(session)
    programs_run = [h["program"] for h in history if h["success"]]

    # Verify the expected sequence appears in order.  Allow
    # extra programs between the expected ones (e.g. dry_run
    # cycles, validation steps).
    expected = [
        step["program"]
        for step in TUTORIAL_SPEC["expected_workflow_path"]
    ]

    # Find each expected program in the order it appears
    j = 0
    for prog in programs_run:
        if j < len(expected) and prog == expected[j]:
            j += 1

    missing = expected[j:]
    assert not missing, (
        "Expected programs in order %r, but missing tail: %r\n"
        "Actual programs run: %r"
        % (expected, missing, programs_run))

    # Verify the anti-pattern doesn't appear
    log = session.get("log", "")
    indicator = TUTORIAL_SPEC["anti_pattern"]["indicator_in_log"]
    if indicator and indicator in log:
        raise AssertionError(
            "Pre-Phase-3d failure indicator found in log: %r\n"
            "This means dock_in_map ran before predict_and_build."
            % indicator)

    print("  PASS — predict→dock chain verified, no failure indicator")


def test_phase_3d_motivating_case_expectations_yaml_format():
    """The YAML entry for tutorial_expectations.yaml is
    syntactically valid.

    Sanity check — the returned string parses as YAML.
    """
    print("Test: phase_3d_motivating_case_expectations_yaml_format")
    try:
        import yaml
    except ImportError:
        print("  SKIP — yaml not available")
        return
    entry = get_expectations_yaml_entry()
    data = yaml.safe_load(entry)
    assert TUTORIAL_NAME in data, (
        "Tutorial name not found in parsed YAML")
    spec = data[TUTORIAL_NAME]
    assert "description" in spec
    assert "inputs" in spec
    assert "expected" in spec
    assert spec["expected"]["stop_after"] == "phenix.dock_in_map"
    programs = spec["expected"]["programs_run_in_order"]
    assert "phenix.predict_and_build" in programs, (
        "predict_and_build must be in expected programs (Phase 3d "
        "fix; pre-fix would skip it)")
    assert "phenix.dock_in_map" in programs, (
        "dock_in_map must be in expected programs (it's the "
        "user's stop target)")
    print("  PASS — YAML entry well-formed")


# =============================================================================
# Test runner (standard format)
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
    test_fns = [v for k, v in sorted(globals().items())
                if k.startswith("test_") and callable(v)]
    passed = 0
    failed = 0
    for fn in test_fns:
        try:
            fn()
            passed += 1
        except Exception as e:
            print("  FAIL: %s" % e)
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    _standalone_runner()
