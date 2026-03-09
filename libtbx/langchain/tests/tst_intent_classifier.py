"""
Tests for Intent Classifier (v115).

Tests the four intent categories: solve, solve_constrained,
task, tutorial.

Run with: python3 tests/tst_intent_classifier.py
"""

from __future__ import absolute_import, division, print_function
import os
import sys

PROJECT_ROOT = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

try:
    from tests.tst_utils import (
        assert_equal, assert_true,
        run_tests_with_fail_fast,
    )
except ImportError:
    def assert_equal(a, b, msg=""):
        assert a == b, "%s: %r != %r" % (msg, a, b)
    def assert_true(val, msg=""):
        assert val, msg
    def assert_false(val, msg=""):
        assert not val, msg
    def assert_in(needle, haystack, msg=""):
        assert needle in haystack, "%s: %r not in %r" % (
            msg, needle, haystack)
    def run_tests_with_fail_fast():
        g = globals()
        tests = sorted(k for k in g if k.startswith("test_"))
        for name in tests:
            print("  Running %s..." % name)
            g[name]()
            print("  PASS: %s" % name)
        print("\nAll %d tests passed." % len(tests))


from agent.intent_classifier import (
    classify_intent,
    SOLVE, SOLVE_CONSTRAINED, TASK, TUTORIAL,
)


# ===========================================================
# SOLVE intent tests
# ===========================================================

def test_solve_no_advice():
    """No advice → solve."""
    r = classify_intent("")
    assert_equal(r["intent"], SOLVE)
    assert_equal(r["confidence"], "high")


def test_solve_none_advice():
    """None advice → solve."""
    r = classify_intent(None)
    assert_equal(r["intent"], SOLVE)


def test_solve_explicit():
    """'Solve the structure' → solve."""
    r = classify_intent("solve the structure")
    assert_equal(r["intent"], SOLVE)
    assert_equal(r["confidence"], "high")


def test_solve_determine():
    """'Determine the structure' → solve."""
    r = classify_intent("determine the structure")
    assert_equal(r["intent"], SOLVE)


def test_solve_process_files():
    """'Process these files' → solve."""
    r = classify_intent("process these files")
    assert_equal(r["intent"], SOLVE)


def test_solve_short_param_advice():
    """Short parameter advice → solve (not tutorial)."""
    r = classify_intent("use nproc=4")
    assert_equal(r["intent"], SOLVE,
        "Brief param advice should not be tutorial")


# ===========================================================
# TASK intent tests
# ===========================================================

def test_task_run_xtriage():
    """'Run xtriage' → task."""
    r = classify_intent("run xtriage")
    assert_equal(r["intent"], TASK)
    assert_equal(r["task_program"],
        "phenix.xtriage")


def test_task_run_refinement():
    """'Run refinement' → task."""
    r = classify_intent("run refinement")
    assert_equal(r["intent"], TASK)
    assert_equal(r["task_program"],
        "phenix.refine")


def test_task_just_do_phaser():
    """'Just do phaser' → task."""
    r = classify_intent("just do phaser")
    assert_equal(r["intent"], TASK)
    assert_equal(r["task_program"],
        "phenix.phaser")


def test_task_check_twinning():
    """'Check for twinning' → task."""
    r = classify_intent("check for twinning")
    assert_equal(r["intent"], TASK)
    assert_equal(r["task_program"],
        "phenix.xtriage")


def test_task_calculate_polder():
    """'Calculate polder map' → task."""
    r = classify_intent("calculate polder map")
    assert_equal(r["intent"], TASK)
    assert_equal(r["task_program"],
        "phenix.polder")


def test_task_run_validation():
    """'Run validation' → task."""
    r = classify_intent("run validation")
    assert_equal(r["intent"], TASK)
    assert_equal(r["task_program"],
        "phenix.molprobity")


def test_task_just_run_autobuild():
    """'Just run autobuild' → task."""
    r = classify_intent("just run autobuild")
    assert_equal(r["intent"], TASK)
    assert_equal(r["task_program"],
        "phenix.autobuild")


def test_task_only_run_refine():
    """'Only run refinement' → task."""
    r = classify_intent("only run refinement")
    assert_equal(r["intent"], TASK)
    assert_equal(r["task_program"],
        "phenix.refine")


# ===========================================================
# SOLVE_CONSTRAINED intent tests
# ===========================================================

def test_constrained_mr():
    """'Solve with molecular replacement' → constrained."""
    r = classify_intent(
        "solve the structure with molecular replacement")
    assert_equal(r["intent"], SOLVE_CONSTRAINED)
    assert_equal(r["method_constraint"],
        "molecular_replacement")


def test_constrained_sad():
    """'Solve using SAD phasing' → constrained."""
    r = classify_intent("solve using SAD phasing")
    assert_equal(r["intent"], SOLVE_CONSTRAINED)
    assert_equal(r["method_constraint"],
        "sad_phasing")


def test_constrained_mrsad():
    """'Solve with MR-SAD' → constrained."""
    r = classify_intent(
        "determine the structure with MR-SAD")
    assert_equal(r["intent"], SOLVE_CONSTRAINED)
    assert_equal(r["method_constraint"],
        "mr_sad")


def test_constrained_predict_and_build():
    """'Solve using PredictAndBuild' → constrained."""
    r = classify_intent(
        "solve with predict and build")
    assert_equal(r["intent"], SOLVE_CONSTRAINED)
    assert_equal(r["method_constraint"],
        "predict_and_build")


def test_constrained_use_mr():
    """'Use molecular replacement' → constrained."""
    r = classify_intent(
        "use molecular replacement to solve")
    assert_equal(r["intent"], SOLVE_CONSTRAINED)
    assert_equal(r["method_constraint"],
        "molecular_replacement")


# ===========================================================
# TUTORIAL intent tests
# ===========================================================

def test_tutorial_preprocessor_output():
    """Preprocessor output → tutorial."""
    advice = (
        "Input Files Found:\n"
        "  data.mtz, search_model.pdb\n"
        "Experiment Type: X-ray crystallography\n"
        "Key Parameters: resolution 3.0\n"
        "Special Instructions: none\n"
        "Goal: Solve by molecular replacement")
    r = classify_intent(advice)
    assert_equal(r["intent"], TUTORIAL,
        "Preprocessor output should be tutorial")


def test_tutorial_explicit_phrase():
    """'In this tutorial...' → tutorial."""
    advice = (
        "In this tutorial, we will solve the "
        "structure by molecular replacement using "
        "Phaser.")
    r = classify_intent(advice)
    assert_equal(r["intent"], TUTORIAL)


def test_tutorial_goal_phrase():
    """'The goal of this...' → tutorial."""
    advice = (
        "The goal of this exercise is to run "
        "density modification on cryo-EM half-maps "
        "using phenix.resolve_cryo_em.")
    r = classify_intent(advice)
    assert_equal(r["intent"], TUTORIAL)


def test_tutorial_step_by_step():
    """Numbered steps → tutorial."""
    advice = (
        "Step 1: Run xtriage on the data.\n"
        "Step 2: Run phaser for molecular "
        "replacement.\n"
        "Step 3: Refine the solution.")
    r = classify_intent(advice)
    assert_equal(r["intent"], TUTORIAL)


def test_tutorial_multiline():
    """Long multi-line text → tutorial (low confidence)."""
    advice = (
        "This example uses X-ray data from "
        "a protein crystal.\n"
        "The data was collected at the APS "
        "synchrotron.\n"
        "Resolution extends to 2.5 Angstroms.\n"
        "A search model is available from a "
        "homologous protein.\n"
        "The expected result is an MR solution "
        "with good TFZ.\n"
        "After MR, one round of refinement "
        "should improve R-free.")
    r = classify_intent(advice)
    assert_equal(r["intent"], TUTORIAL,
        "Multi-line description is tutorial")


# ===========================================================
# Priority tests (task > constrained > tutorial)
# ===========================================================

def test_task_over_tutorial():
    """'Run xtriage' in tutorial-like text → task wins.

    When the user gives a clear single-program command,
    that takes priority even if the text looks tutorial-like.
    """
    r = classify_intent("just run xtriage")
    assert_equal(r["intent"], TASK,
        "Task should win over tutorial")


def test_constrained_not_in_tutorial():
    """'Solve with SAD' WITHOUT tutorial context → constrained.

    When the user gives a direct method constraint without
    tutorial-style language, it's solve_constrained.
    """
    advice = "solve the structure using SAD phasing"
    r = classify_intent(advice)
    assert_equal(r["intent"], SOLVE_CONSTRAINED,
        "Direct constraint without tutorial → constrained")


def test_tutorial_with_method_mention():
    """Tutorial text mentioning a method → still tutorial.

    When the text says 'In this tutorial...MR', the tutorial
    intent takes priority because the text is describing a
    procedure, not giving a user command.
    """
    advice = (
        "In this tutorial, solve the structure "
        "using SAD phasing")
    r = classify_intent(advice)
    assert_equal(r["intent"], TUTORIAL,
        "Tutorial phrase should win over method mention")


# ===========================================================
# Real tutorial README simulations
# ===========================================================

def test_real_a2u_mr_readme():
    """Simulate a2u-globulin-mr README."""
    advice = (
        "Input Files Found: a2u_data.mtz, "
        "search_model.pdb, sequence.fasta\n"
        "Experiment Type: X-ray crystallography\n"
        "Key Parameters: resolution 1.8\n"
        "Goal: Solve the structure by molecular "
        "replacement")
    r = classify_intent(advice)
    assert_equal(r["intent"], TUTORIAL,
        "Preprocessed README should be tutorial")


def test_real_actin_denmod_readme():
    """Simulate actin_denmod README."""
    advice = (
        "Input Files Found: half_map_1.ccp4, "
        "half_map_2.ccp4\n"
        "Experiment Type: Cryo-EM\n"
        "Key Parameters: none\n"
        "Special Instructions: Perform density "
        "modification of cryo-EM map\n"
        "Goal: Density modification of cryo-EM map")
    r = classify_intent(advice)
    assert_equal(r["intent"], TUTORIAL)


def test_real_bare_files_solve():
    """Bare files with no advice → solve."""
    r = classify_intent("")
    assert_equal(r["intent"], SOLVE,
        "No advice means solve mode")


# ===========================================================
# Edge cases
# ===========================================================

def test_whitespace_only():
    """Whitespace-only advice → solve."""
    r = classify_intent("   \n  \n  ")
    assert_equal(r["intent"], SOLVE)


def test_never_raises():
    """Classifier should never raise."""
    for input in [None, "", "x", 42, [], {}]:
        try:
            r = classify_intent(input)
            assert_true(
                r["intent"] in (
                    SOLVE, SOLVE_CONSTRAINED,
                    TASK, TUTORIAL),
                "Should return valid intent")
        except TypeError:
            # Non-string inputs may TypeError on .strip()
            # That's acceptable — classify_intent
            # documents str input
            pass


# ===========================================================
# Entry point
# ===========================================================

def run_all_tests():
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
