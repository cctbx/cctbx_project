"""
Fabricated Stop Condition lines stripped before directive extraction (v115).

LLM-fabricated "Stop Condition:" lines from the advice preprocessor are
stripped before directive extraction while real user stop commands are kept.
"""

from __future__ import absolute_import, division, print_function
import os
import sys

# Ensure project root is on path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

try:
  from tests.tst_utils import (
    assert_equal, assert_true, assert_false,
    assert_in, run_tests_with_fail_fast,
  )
except ImportError:
  # Minimal fallback when tst_utils not available
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

from agent.directive_extractor import (
  _strip_preprocessor_stop_condition,
  extract_directives_simple,
)


# =====================================================================
# Tests for _strip_preprocessor_stop_condition()
# =====================================================================

def test_strip_fabricated_mr_stop():
  """Strip 'Stop after molecular replacement' from preprocessed advice."""
  advice = (
    "Input Files Found: data.mtz, model.pdb\n"
    "Primary Goal: Solve structure by MR\n"
    "Stop Condition: Stop after molecular replacement "
    "(i.e., once Phaser returns an MR solution).\n"
    "Key Parameters: resolution 2.5"
  )
  result = _strip_preprocessor_stop_condition(advice)
  assert_false(
    "Stop Condition" in result,
    "Fabricated stop condition should be stripped")
  assert_false(
    "Primary Goal" in result,
    "Primary Goal header should be stripped")
  assert_in("Input Files Found", result,
    "Other content should be preserved")
  assert_in("Key Parameters", result,
    "Content after stop line should be preserved")


def test_strip_fabricated_autobuild_stop():
  """Strip 'Stop after AutoBuild' from preprocessed advice."""
  advice = (
    "Input Files Found: data.mtz, model.pdb\n"
    "Experiment Type: X-ray crystallography\n"
    "Goal: Rebuild model\n"
    "Stop Condition: Stop after AutoBuild\n"
    "Key Parameters: resolution 3.0"
  )
  result = _strip_preprocessor_stop_condition(advice)
  assert_false(
    "Stop Condition" in result,
    "Fabricated AutoBuild stop should be stripped")
  assert_false(
    "Goal:" in result,
    "Goal header should be stripped in preprocessor context")
  assert_in("Key Parameters", result,
    "Key Parameters should be preserved")


def test_strip_goal_with_denmod():
  """Strip 'Goal: Density modification' that triggers denmod patterns."""
  advice = (
    "Input Files Found: half_map1.mrc, half_map2.mrc\n"
    "Experiment Type: cryo-EM\n"
    "Goal: Density modification of cryo-EM map\n"
    "Key Parameters: None\n"
    "Special Instructions: None"
  )
  result = _strip_preprocessor_stop_condition(advice)
  assert_false(
    "Goal:" in result,
    "Goal header should be stripped in preprocessor context")
  assert_in("Input Files Found", result,
    "Input Files Found should be preserved")


def test_preserve_goal_in_prose():
  """'goal' in a sentence is NOT stripped (not a header)."""
  advice = "The goal is to refine this model to completion."
  result = _strip_preprocessor_stop_condition(advice)
  assert_in("goal is to refine", result,
    "Prose 'goal' in sentence must be preserved")


def test_preserve_raw_user_goal_header():
  """'Goal: R-free below 0.25' preserved when no preprocessor context."""
  advice = "Goal: get R-free below 0.25 and fit the ligand"
  result = _strip_preprocessor_stop_condition(advice)
  assert_in("R-free below 0.25", result,
    "Raw user 'Goal:' must be preserved without "
    "preprocessor context")
  assert_in("fit the ligand", result,
    "Full user text must be preserved")


def test_strip_goal_with_preprocessor_context():
  """'Goal:' stripped when preprocessor signatures are present."""
  advice = (
    "Experiment Type: X-ray crystallography\n"
    "Goal: Solve by molecular replacement\n"
    "Key Parameters: resolution 2.5"
  )
  result = _strip_preprocessor_stop_condition(advice)
  assert_false(
    "Goal:" in result,
    "Goal should be stripped in preprocessor context")
  assert_in("Experiment Type", result,
    "Experiment Type preserved")
  assert_in("Key Parameters", result,
    "Key Parameters preserved")


def test_strip_none_stop():
  """Strip 'Stop Condition: None' (also preprocessor output)."""
  advice = (
    "Goal: Rebuild model\n"
    "Stop Condition: None\n"
    "Parameters: resolution 3.0"
  )
  result = _strip_preprocessor_stop_condition(advice)
  assert_false(
    "Stop Condition" in result,
    "'None' stop should also be stripped")
  assert_in("Parameters: resolution 3.0", result,
    "Content should be preserved")


def test_strip_case_insensitive():
  """Strip regardless of case."""
  advice = (
    "STOP CONDITION: Stop after autobuild completes.\n"
    "Other stuff here."
  )
  result = _strip_preprocessor_stop_condition(advice)
  assert_false(
    "STOP CONDITION" in result,
    "Case-insensitive stripping should work")
  assert_in("Other stuff here.", result,
    "Other content preserved")


def test_preserve_real_user_stop():
  """Real user prose 'stop after refinement' is NOT stripped."""
  advice = "stop after refinement and validate the model"
  result = _strip_preprocessor_stop_condition(advice)
  assert_equal(
    result.strip(), advice.strip(),
    "Real user stop prose must be preserved")


def test_preserve_real_user_stop_in_sentence():
  """'stop' in a sentence is not a 'Stop Condition:' line."""
  advice = (
    "Please stop after the first refinement cycle.\n"
    "Then I will inspect the model."
  )
  result = _strip_preprocessor_stop_condition(advice)
  assert_in("stop after the first", result,
    "Prose stop in sentence must be preserved")


def test_empty_input():
  """Empty and None inputs are handled safely."""
  assert_equal(
    _strip_preprocessor_stop_condition(""), "",
    "Empty string should return empty")
  assert_equal(
    _strip_preprocessor_stop_condition(None), None,
    "None should return None")


def test_no_blank_line_explosion():
  """Stripping doesn't create excessive blank lines."""
  advice = (
    "Line one\n\n"
    "Stop Condition: Stop after MR\n\n"
    "Line two"
  )
  result = _strip_preprocessor_stop_condition(advice)
  # Should not have more than 2 consecutive newlines
  assert_false(
    "\n\n\n" in result,
    "Should not have 3+ consecutive newlines")


# =====================================================================
# Tests for extract_directives_simple() with fabricated stops
# =====================================================================

def test_simple_no_fabricated_mr_stop():
  """Simple extraction ignores fabricated MR stop from preprocessor."""
  preprocessed = (
    "Input Files Found: data.mtz, model.pdb, seq.fasta\n"
    "Experiment Type: X-ray, molecular replacement\n"
    "Primary Goal: Solve by molecular replacement\n"
    "Stop Condition: Stop after molecular replacement "
    "(i.e., once Phaser returns an MR solution).\n"
    "Key Parameters: resolution 2.5"
  )
  result = extract_directives_simple(preprocessed)
  after_prog = result.get(
    "stop_conditions", {}).get("after_program", "")
  assert_false(
    after_prog,
    "Fabricated MR stop should not produce after_program "
    "(got %s)" % after_prog)


def test_simple_no_fabricated_denmod_stop():
  """Simple extraction ignores fabricated denmod stop."""
  preprocessed = (
    "Input Files Found: half_map1.mrc, half_map2.mrc\n"
    "Experiment Type: cryo-EM\n"
    "Goal: Density modification of cryo-EM map\n"
    "Stop Condition: Stop after running resolve_cryo_em "
    "once\n"
    "Key Parameters: None"
  )
  result = extract_directives_simple(preprocessed)
  after_prog = result.get(
    "stop_conditions", {}).get("after_program", "")
  assert_false(
    after_prog,
    "Fabricated denmod stop should not produce "
    "after_program (got %s)" % after_prog)


def test_simple_all_plan_fabricated_stops():
  """All 7 fabricated stops from PLAN.md are blocked."""
  cases = [
    ("a2u-globulin-mr",
     "Stop Condition: Stop after molecular replacement"),
    ("ICAM-mr",
     "Stop Condition: Stop after MR placement and "
     "initial refinement"),
    ("AF_exoV_PredictAndBuild",
     "Stop Condition: Stop after PredictAndBuild "
     "completes"),
    ("7rpq_AF_reference_model",
     "Stop Condition: Stop after two refinement runs"),
    ("AF_exoV_MRSAD",
     "Stop Condition: Stop after MR and MR-SAD"),
    ("a2u-globulin-rebuild",
     "Stop Condition: Stop after AutoBuild"),
    ("actin_denmod",
     "Stop Condition: Stop after running "
     "resolve_cryo_em once"),
  ]
  for tutorial, stop_line in cases:
    advice = (
      "Input Files Found: data.mtz, model.pdb\n"
      "Experiment Type: X-ray crystallography\n"
      "Goal: Solve structure\n"
      "%s\n"
      "Key Parameters: None"
    ) % stop_line
    result = extract_directives_simple(advice)
    after_prog = result.get(
      "stop_conditions", {}).get("after_program", "")
    assert_false(
      after_prog,
      "%s: fabricated stop extracted as after_program=%s"
      % (tutorial, after_prog))


# =====================================================================
# Tests that real stop commands STILL WORK
# =====================================================================

def test_simple_real_stop_after_refine():
  """'stop after refinement' still produces after_program."""
  result = extract_directives_simple(
    "stop after refinement, skip validation")
  stop = result.get("stop_conditions", {})
  assert_equal(
    stop.get("after_program"), "phenix.refine",
    "Real 'stop after refinement' must be preserved")
  assert_true(
    stop.get("skip_validation"),
    "skip_validation should be set")


def test_simple_real_run_xtriage():
  """'run xtriage to check for twinning' still works."""
  result = extract_directives_simple(
    "run xtriage to check for twinning")
  stop = result.get("stop_conditions", {})
  assert_equal(
    stop.get("after_program"), "phenix.xtriage",
    "Explicit 'run xtriage' command must be preserved")


def test_simple_real_stop_after_cycle():
  """'stop after cycle 3' still works."""
  result = extract_directives_simple(
    "stop after cycle 3")
  stop = result.get("stop_conditions", {})
  assert_equal(
    stop.get("after_cycle"), 3,
    "Explicit 'stop after cycle N' must be preserved")


def test_simple_real_run_phaser():
  """'just run phaser' still works."""
  result = extract_directives_simple(
    "just run phaser to test MR solution")
  stop = result.get("stop_conditions", {})
  assert_equal(
    stop.get("after_program"), "phenix.phaser",
    "Explicit 'run phaser' command must be preserved")


def test_simple_real_polder():
  """'calculate polder map' still works."""
  result = extract_directives_simple(
    "calculate a polder map for chain B residue 100")
  stop = result.get("stop_conditions", {})
  assert_equal(
    stop.get("after_program"), "phenix.polder",
    "Explicit polder command must be preserved")


# =====================================================================
# Entry point
# =====================================================================

def run_all_tests():
  run_tests_with_fail_fast()


if __name__ == "__main__":
  run_all_tests()
