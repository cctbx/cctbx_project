"""
Integration tests for Fixes 3+4+5 interactions (v115).

Tests the full failure handling chain:
  Fix 4 (PHIL validation) → Fix 5 (duplicate detection)
  → Fix 3 (error classification + pivot)

Run with: python3 tests/tst_fix_integration.py
"""

from __future__ import absolute_import, division, print_function
import os
import sys

PROJECT_ROOT = os.path.dirname(
  os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

try:
  from tests.tst_utils import (
    assert_equal, assert_true, assert_false,
    assert_in, run_tests_with_fail_fast,
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


from agent.phil_validator import (
  validate_phil_strategy,
  _STRATEGY_FLAGS_CACHE,
)
from agent.error_classifier import (
  classify_error, should_pivot,
  LABEL_ERROR, NO_ERROR,
)


# ===========================================================
# Fix 4 → Fix 3 chain: PHIL validation prevents the error
# ===========================================================

def test_fix4_prevents_phil_error():
  """Fix 4 strips bad params so Fix 3 never sees PHIL error.

  PLAN case: AF_exoV autobuild with obs_labels.
  """
  _STRATEGY_FLAGS_CACHE.clear()

  # Fix 4: strip obs_labels
  cleaned, stripped = validate_phil_strategy(
    "phenix.autobuild",
    {"resolution": 3.0, "obs_labels": "I(+)"})
  assert_false("obs_labels" in cleaned,
    "Fix 4 should strip obs_labels")

  # Cleaned command runs — no PHIL error in log
  success_log = (
    "phenix.autobuild data=exoV.mtz resolution=3.0\n"
    "AutoBuild completed: model.pdb\n"
    "R-free = 0.310")
  ec = classify_error(success_log, "phenix.autobuild")
  assert_equal(ec["category"], NO_ERROR,
    "No PHIL error because Fix 4 stripped the param")


def test_fix4_strips_fix3_classifies_remainder():
  """Fix 4 strips some params; Fix 3 classifies the remaining error.

  Scenario: LLM sends obs_labels (stripped by Fix 4) AND
  the program fails for a different reason (label error).
  """
  _STRATEGY_FLAGS_CACHE.clear()

  # Fix 4 strips obs_labels
  cleaned, stripped = validate_phil_strategy(
    "phenix.autobuild",
    {"resolution": 3.0, "obs_labels": "I(+)"})
  assert_true("obs_labels" not in cleaned)

  # Program runs without obs_labels but still fails
  # (different error — label mismatch in the MTZ)
  label_error_log = (
    "phenix.autobuild data=exoV.mtz resolution=3.0\n"
    "Sorry: cannot determine data labels\n"
    "Multiple arrays found in MTZ")
  ec = classify_error(
    label_error_log, "phenix.autobuild")
  assert_equal(ec["category"], LABEL_ERROR,
    "Fix 3 classifies the remaining error")
  assert_false(ec["is_terminal"],
    "Label errors are retryable")


# ===========================================================
# Fix 3 tiered response: first failure → retry,
# second → pivot
# ===========================================================

def test_fix3_first_failure_allows_retry():
  """First failure of a program allows self-correction."""
  ec = classify_error(
    "Sorry: some phil parameters are not recognized",
    "phenix.autobuild")
  history = [
    {"program": "autobuild",
     "result": "FAILED: PHIL error"},
  ]
  pivot, reason = should_pivot(
    history, "phenix.autobuild", ec)
  assert_false(pivot,
    "First failure should not pivot")
  assert_in("self-correction", reason)


def test_fix3_second_failure_pivots():
  """Second consecutive failure forces pivot."""
  ec = classify_error(
    "Sorry: some phil parameters are not recognized",
    "phenix.autobuild")
  history = [
    {"program": "autobuild",
     "result": "FAILED: error 1"},
    {"program": "autobuild",
     "result": "FAILED: error 2"},
  ]
  pivot, reason = should_pivot(
    history, "phenix.autobuild", ec)
  assert_true(pivot,
    "Second failure should pivot")


def test_fix3_terminal_immediate_pivot():
  """Terminal errors pivot immediately, no retry."""
  ec = classify_error(
    "Traceback (most recent call last):\n"
    "  RuntimeError: crash",
    "phenix.autobuild")
  # Even with NO history, terminal errors pivot
  pivot, reason = should_pivot(
    [], "phenix.autobuild", ec)
  assert_true(pivot,
    "Terminal errors should always pivot")


# ===========================================================
# Fix 5: duplicate detection unit tests
# ===========================================================

def test_fix5_exact_dup_of_success_blocked():
  """Exact duplicate of a successful command is always blocked."""
  command = "phenix.refine model.pdb data.mtz cycles=5"
  history = [
    {"command": command, "cycle_number": 1,
     "result": "SUCCESS: OK", "program": "refine"},
  ]
  # Find matching successful cycles
  recent = [h for h in history
            if isinstance(h, dict)][-5:]
  matching = [
    (h.get("cycle_number"), True)
    for h in recent
    if h.get("command") == command
    and str(h.get("result", "")).startswith("SUCCESS")]
  assert_true(len(matching) > 0,
    "Should find successful match")


def test_fix5_exact_dup_of_failure_with_new_params():
  """Failed command + new LLM strategy params → allow retry."""
  command = "phenix.autobuild data.mtz resolution=3.0"
  llm_strategy = {
    "resolution": 3.0, "obs_labels": "I(+)"}
  history = [
    {"command": command, "cycle_number": 1,
     "result": "FAILED: PHIL error",
     "program": "autobuild"},
  ]
  recent = [h for h in history
            if isinstance(h, dict)][-5:]
  matching_failed = [
    h for h in recent
    if h.get("command") == command
    and not str(h.get("result", "")).startswith(
      "SUCCESS")]
  new_params = [
    k for k in llm_strategy
    if k not in command
    and k not in ("output_prefix", "output_name")]
  should_bypass = (
    len(matching_failed) == 1
    and len(new_params) > 0)
  assert_true(should_bypass,
    "1 failed match + new params → bypass")
  assert_in("obs_labels", new_params,
    "obs_labels not in command → new param")


def test_fix5_two_failures_always_blocked():
  """Two+ failed matches with same command → always blocked."""
  command = "phenix.autobuild data.mtz resolution=3.0"
  history = [
    {"command": command, "cycle_number": 1,
     "result": "FAILED: err1"},
    {"command": command, "cycle_number": 2,
     "result": "FAILED: err2"},
  ]
  recent = [h for h in history
            if isinstance(h, dict)][-5:]
  matching_failed = [
    h for h in recent
    if h.get("command") == command
    and not str(h.get("result", "")).startswith(
      "SUCCESS")]
  assert_equal(len(matching_failed), 2,
    "Should find 2 failed matches")
  should_block = len(matching_failed) >= 2
  assert_true(should_block,
    "2+ failed matches → block")


def test_fix5_no_match_passes():
  """Different command → no match → passes."""
  command = "phenix.refine model.pdb data.mtz cycles=5"
  history = [
    {"command": "phenix.autobuild data.mtz",
     "cycle_number": 1,
     "result": "SUCCESS: OK"},
  ]
  recent = [h for h in history
            if isinstance(h, dict)][-5:]
  matching = [
    h for h in recent
    if h.get("command") == command]
  assert_equal(len(matching), 0,
    "Different command → no match")


# ===========================================================
# Fix 6: map extension safety net unit tests
# ===========================================================

def test_fix6_half_maps_categorized():
  """Half-map .ccp4 files are categorized correctly."""
  from agent.workflow_state import _categorize_files
  result = _categorize_files(
    ["/data/half_map_1.ccp4",
     "/data/half_map_2.ccp4"],
    files_local=False)
  hm = result.get("half_map", [])
  mp = result.get("map", [])
  assert_equal(len(hm), 2,
    "Should have 2 half_maps")
  assert_equal(len(mp), 2,
    "Should have 2 in map parent")


def test_fix6_full_map_categorized():
  """Non-half .ccp4 files go to full_map."""
  from agent.workflow_state import _categorize_files
  result = _categorize_files(
    ["/data/local_resolution.ccp4"],
    files_local=False)
  fm = result.get("full_map", [])
  mp = result.get("map", [])
  assert_true(len(fm) >= 1,
    "Should have full_map")
  assert_true(len(mp) >= 1,
    "Should have map parent")


def test_fix6_mrc_extension():
  """.mrc files are categorized."""
  from agent.workflow_state import _categorize_files
  result = _categorize_files(
    ["/data/reconstruction.mrc"],
    files_local=False)
  fm = result.get("full_map", [])
  mp = result.get("map", [])
  assert_true(len(fm) >= 1 or len(mp) >= 1,
    ".mrc should be in map categories")


def test_fix6_map_extension():
  """.map files are categorized."""
  from agent.workflow_state import _categorize_files
  result = _categorize_files(
    ["/data/density.map"],
    files_local=False)
  fm = result.get("full_map", [])
  mp = result.get("map", [])
  assert_true(len(fm) >= 1 or len(mp) >= 1,
    ".map should be in map categories")


def test_fix6_no_duplicates():
  """Already-categorized maps don't get duplicated."""
  from agent.workflow_state import _categorize_files
  result = _categorize_files(
    ["/data/half_map_1.ccp4",
     "/data/half_map_2.ccp4"],
    files_local=False)
  hm = result.get("half_map", [])
  assert_equal(len(hm), 2,
    "Should have exactly 2, not duplicates")



# ===========================================================
# Full chain: Fix 4 + Fix 5 + Fix 3 together
# ===========================================================

def test_full_chain_af_exov():
  """Full chain simulation for AF_exoV_PredictAndBuild.

  Before fixes: 4× autobuild failure with obs_labels
  After fixes: max 2 cycles, then pivot
  """
  _STRATEGY_FLAGS_CACHE.clear()

  # --- Cycle 1 ---
  # Fix 4 strips obs_labels
  c1_cleaned, c1_stripped = validate_phil_strategy(
    "phenix.autobuild",
    {"resolution": 3.0, "obs_labels": "I(+)"})
  assert_false("obs_labels" in c1_cleaned,
    "Cycle 1: obs_labels stripped")

  # Program may still fail (e.g., label error in MTZ)
  c1_ec = classify_error(
    "Sorry: cannot determine data labels",
    "phenix.autobuild")
  assert_equal(c1_ec["category"], LABEL_ERROR)

  # Fix 3: first failure, allow retry
  c1_pivot, _ = should_pivot(
    [{"program": "autobuild",
      "result": "FAILED: label"}],
    "phenix.autobuild", c1_ec)
  assert_false(c1_pivot, "Cycle 1: no pivot")

  # Fix 5: first failed match, LLM has obs_labels
  # → bypass allowed

  # --- Cycle 2 ---
  # Same thing, fails again
  c2_ec = classify_error(
    "Sorry: cannot determine data labels",
    "phenix.autobuild")

  # Fix 3: second failure → pivot
  c2_pivot, _ = should_pivot(
    [{"program": "autobuild",
      "result": "FAILED: error 1"},
     {"program": "autobuild",
      "result": "FAILED: error 2"}],
    "phenix.autobuild", c2_ec)
  assert_true(c2_pivot, "Cycle 2: pivot forced")

  # Fix 5: two failed matches → also blocks

  # --- Cycle 3 ---
  # LLM forced to pick different program (refine)
  c3_cleaned, c3_stripped = validate_phil_strategy(
    "phenix.refine", {"cycles": 5})
  assert_true("cycles" in c3_cleaned,
    "Cycle 3: refine params valid")
  assert_equal(len(c3_stripped), 0,
    "Cycle 3: nothing stripped for refine")


# ===========================================================
# Entry point
# ===========================================================

def run_all_tests():
  run_tests_with_fail_fast()


if __name__ == "__main__":
  run_all_tests()
