"""
Tests for Fix 4: PHIL parameter validation (v115).

Verifies that LLM-injected PHIL parameters are validated against
each program's strategy_flags from programs.yaml, with
unrecognized params stripped before command assembly.

Run with: python3 tests/tst_phil_validation.py
"""

from __future__ import absolute_import, division, print_function
import os
import sys

# Ensure project root is on path
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


# Import after path setup
from agent.phil_validator import (
  validate_phil_strategy,
  _STRATEGY_FLAGS_CACHE,
)

# Clear cache before tests (may have None from failed loads)
_STRATEGY_FLAGS_CACHE.clear()


# =============================================================
# Tests for validate_phil_strategy()
# =============================================================

def test_valid_autobuild_params():
  """Valid autobuild params pass through unchanged."""
  strategy = {
    "rebuild_in_place": False,
    "resolution": 2.5,
    "nproc": 4,
  }
  cleaned, stripped = validate_phil_strategy(
    "phenix.autobuild", strategy)
  assert_equal(cleaned, strategy,
    "Valid params should pass through unchanged")


def test_strip_invalid_autobuild_params():
  """Invalid params for autobuild are stripped."""
  strategy = {
    "rebuild_in_place": False,     # valid
    "obs_labels": "I(+)",          # INVALID for autobuild
    "resolution": 2.5,             # valid
    "wavelength": 0.9792,          # INVALID for autobuild
  }
  cleaned, stripped = validate_phil_strategy(
    "phenix.autobuild", strategy)
  assert_true("rebuild_in_place" in cleaned,
    "Valid param should be kept")
  assert_true("resolution" in cleaned,
    "Valid param should be kept")
  assert_false("obs_labels" in cleaned,
    "obs_labels is not a valid autobuild param")
  assert_false("wavelength" in cleaned,
    "wavelength is not a valid autobuild param")
  # Check stripped list
  stripped_keys = [k for k, v in stripped]
  assert_in("obs_labels", stripped_keys,
    "obs_labels should be in stripped list")
  assert_in("wavelength", stripped_keys,
    "wavelength should be in stripped list")


def test_strip_invalid_process_predicted_model():
  """Invalid params for process_predicted_model are stripped."""
  strategy = {
    "maximum_rmsd": 1.5,           # valid
    "wavelength": 0.9792,          # INVALID
    "anomalous": True,             # INVALID
  }
  cleaned, stripped = validate_phil_strategy(
    "phenix.process_predicted_model", strategy)
  assert_true("maximum_rmsd" in cleaned,
    "maximum_rmsd is valid for process_predicted_model")
  assert_false("wavelength" in cleaned,
    "wavelength is invalid for process_predicted_model")
  assert_false("anomalous" in cleaned,
    "anomalous is invalid for process_predicted_model")


def test_empty_strategy():
  """Empty strategy passes through."""
  cleaned, stripped = validate_phil_strategy(
    "phenix.refine", {})
  assert_equal(cleaned, {},
    "Empty strategy should return empty")


def test_none_strategy():
  """None strategy passes through."""
  cleaned, stripped = validate_phil_strategy(
    "phenix.refine", None)
  assert_equal(cleaned, None,
    "None strategy should return None")


def test_no_strategy_flags_program():
  """Programs with no strategy_flags pass through."""
  strategy = {"some_param": True}
  cleaned, stripped = validate_phil_strategy(
    "phenix.molprobity", strategy)
  # molprobity has no strategy_flags, so validation can't
  # strip anything — params pass through
  assert_equal(cleaned, strategy,
    "Programs without strategy_flags should pass through")


def test_preserve_build_pipeline_keys():
  """Build pipeline keys like 'ligand' are never stripped."""
  strategy = {
    "ligand": "ATP",               # pipeline key
    "output_prefix": "cycle5",     # pipeline key
    "nproc": 4,                    # valid PHIL
    "obs_labels": "I(+)",          # INVALID
  }
  cleaned, stripped = validate_phil_strategy(
    "phenix.ligandfit", strategy)
  assert_true("ligand" in cleaned,
    "ligand is a pipeline key, must be preserved")
  assert_true("output_prefix" in cleaned,
    "output_prefix is a pipeline key, must be preserved")
  assert_true("nproc" in cleaned,
    "nproc is valid for ligandfit")
  assert_false("obs_labels" in cleaned,
    "obs_labels is invalid for ligandfit")


def test_refine_valid_params():
  """Valid refine params pass through."""
  strategy = {
    "cycles": 5,
    "simulated_annealing": True,
    "ordered_solvent": True,
    "output_prefix": "refine_001",
  }
  cleaned, stripped = validate_phil_strategy(
    "phenix.refine", strategy)
  assert_equal(cleaned, strategy,
    "All valid refine params should pass through")


def test_refine_strip_invalid():
  """Invalid params for refine are stripped."""
  strategy = {
    "cycles": 5,                   # valid
    "atom_type": "Se",             # INVALID (autosol only)
    "stop_after_predict": True,    # INVALID (P&B only)
  }
  cleaned, stripped = validate_phil_strategy(
    "phenix.refine", strategy)
  assert_true("cycles" in cleaned,
    "cycles is valid for refine")
  assert_false("atom_type" in cleaned,
    "atom_type is invalid for refine")
  assert_false("stop_after_predict" in cleaned,
    "stop_after_predict is invalid for refine")


def test_unknown_program():
  """Unknown programs pass through (YAML lookup returns None)."""
  strategy = {"some_param": True}
  cleaned, stripped = validate_phil_strategy(
    "phenix.unknown_program", strategy)
  assert_equal(cleaned, strategy,
    "Unknown programs should pass through")


def test_plan_case_obs_labels_autobuild():
  """PLAN case: obs_labels injected into autobuild (AF_exoV)."""
  # From PLAN: AF_exoV_PredictAndBuild/llm failed 4x because
  # LLM injected autobuild.input.xray_data.obs_labels=I(+)
  strategy = {
    "resolution": 3.0,
    "obs_labels": "I(+)",
  }
  cleaned, stripped = validate_phil_strategy(
    "phenix.autobuild", strategy)
  assert_false("obs_labels" in cleaned,
    "obs_labels must be stripped from autobuild "
    "(caused 4 failures in AF_exoV)")
  assert_true("resolution" in cleaned,
    "resolution is valid for autobuild")


def test_plan_case_wavelength_process_predicted_model():
  """PLAN case: wavelength passed to process_predicted_model."""
  # From PLAN: AF_exoV_MRSAD/expert failed 3x because
  # wavelength and anomalous were passed to PPM
  strategy = {
    "maximum_rmsd": 1.5,
    "wavelength": 0.9792,
    "anomalous": True,
  }
  cleaned, stripped = validate_phil_strategy(
    "phenix.process_predicted_model", strategy)
  assert_true("maximum_rmsd" in cleaned,
    "maximum_rmsd is the correct PPM param")
  assert_false("wavelength" in cleaned,
    "wavelength must be stripped from PPM")
  assert_false("anomalous" in cleaned,
    "anomalous must be stripped from PPM")


def test_rewrite_refine_reference_model():
  """Params from _STRATEGY_REWRITES must survive validation.

  graph_nodes.py rewrites 'reference_model.enabled' →
  'reference_model_enabled' for phenix.refine.  The validator
  must not strip the rewritten key.
  """
  _STRATEGY_FLAGS_CACHE.clear()
  strategy = {
    "cycles": 5,
    "reference_model_enabled": True,
    "reference_model_use_starting": True,
  }
  cleaned, stripped = validate_phil_strategy(
    "phenix.refine", strategy)
  assert_true("reference_model_enabled" in cleaned,
    "reference_model_enabled is a valid rewritten param")
  assert_true("reference_model_use_starting" in cleaned,
    "reference_model_use_starting is a valid rewritten "
    "param")
  assert_equal(len(stripped), 0,
    "No params should be stripped")


def test_rewrite_resolve_cryo_em_mask_atoms():
  """mask_atoms must be BLOCKED for phenix.resolve_cryo_em.

  Phase 3 Bug 3: mask_atoms=True was interpreted by PHENIX as
  strategy.mask_atoms_atom_radius="True" (numeric field) causing
  RuntimeError.  mask_atoms was removed from strategy_flags and
  added to _BLOCKED_PARAMS.  The validator must strip it.
  """
  _STRATEGY_FLAGS_CACHE.clear()
  strategy = {
    "resolution": 3.0,
    "mask_atoms": True,
  }
  cleaned, stripped = validate_phil_strategy(
    "phenix.resolve_cryo_em", strategy)
  assert_false("mask_atoms" in cleaned,
    "mask_atoms must be BLOCKED (causes RuntimeError)")
  assert_true("resolution" in cleaned,
    "resolution is a valid param")
  stripped_keys = [k for k, v in stripped]
  assert_true("mask_atoms" in stripped_keys,
    "mask_atoms must appear in stripped list")


def test_plan_case_rmsd_ambiguity():
  """PLAN case: rmsd is ambiguous for process_predicted_model.

  From PLAN: AF_exoV_MRSAD/llm_think had 3x failures
  because 'rmsd=1.0' matches multiple PHIL scopes in PPM.
  The correct param is 'maximum_rmsd'.  Our validator
  should strip 'rmsd' (not in strategy_flags) and keep
  'maximum_rmsd' (is in strategy_flags).
  """
  _STRATEGY_FLAGS_CACHE.clear()
  # Case 1: both present
  strategy = {"rmsd": 1.0, "maximum_rmsd": 1.5}
  cleaned, stripped = validate_phil_strategy(
    "phenix.process_predicted_model", strategy)
  assert_true("maximum_rmsd" in cleaned,
    "maximum_rmsd is the correct PPM param")
  assert_false("rmsd" in cleaned,
    "Ambiguous rmsd must be stripped")

  # Case 2: only rmsd (LLM used wrong name)
  _STRATEGY_FLAGS_CACHE.clear()
  strategy2 = {"rmsd": 1.0}
  cleaned2, stripped2 = validate_phil_strategy(
    "phenix.process_predicted_model", strategy2)
  assert_false("rmsd" in cleaned2,
    "Ambiguous rmsd must be stripped even alone")
  assert_equal(len(stripped2), 1,
    "rmsd should appear in stripped list")


# =============================================================
# Entry point
# =============================================================

def run_all_tests():
  run_tests_with_fail_fast()


if __name__ == "__main__":
  run_all_tests()
