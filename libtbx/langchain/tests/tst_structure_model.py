"""
Unit tests for the Structure Model (Phase 1, Step 1.4).

Run standalone:
  python tests/tst_structure_model.py

Run from PHENIX:
  libtbx.python tests/tst_structure_model.py

No PHENIX dependencies — pure Python with dict inputs.

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import json
import os
import sys
import traceback

# Ensure project root on sys.path
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_THIS_DIR)
if _PROJECT_ROOT not in sys.path:
  sys.path.insert(0, _PROJECT_ROOT)

from agent.structure_model import (
  StructureModel, Hypothesis, _safe_float, _deep_merge,
)


def run_tests():
  """Run all Structure Model tests."""
  passed = 0
  failed = 0
  skipped = 0

  def test(name, fn):
    nonlocal passed, failed, skipped
    try:
      result = fn()
      if result == "SKIP":
        print("  SKIP: %s" % name)
        skipped += 1
      else:
        print("  PASS: %s" % name)
        passed += 1
    except Exception as e:
      print("  FAIL: %s — %s" % (name, e))
      traceback.print_exc()
      failed += 1

  print("=" * 60)
  print("Structure Model Unit Tests")
  print("=" * 60)
  print()

  # --- Hypothesis ---
  print("Hypothesis data structure")
  test("hypothesis_init_defaults",
    test_hypothesis_init_defaults)
  test("hypothesis_status_validation",
    test_hypothesis_status_validation)
  test("hypothesis_is_active",
    test_hypothesis_is_active)
  test("hypothesis_is_resolved",
    test_hypothesis_is_resolved)
  test("hypothesis_roundtrip",
    test_hypothesis_roundtrip)
  test("hypothesis_from_dict_tolerant",
    test_hypothesis_from_dict_tolerant)
  print()

  # --- StructureModel init and basics ---
  print("StructureModel init and basics")
  test("init_defaults",
    test_init_defaults)
  test("safe_float",
    test_safe_float)
  test("deep_merge",
    test_deep_merge)
  print()

  # --- update_from_xtriage ---
  print("update_from_xtriage")
  test("xtriage_basic",
    test_xtriage_basic)
  test("xtriage_twinning",
    test_xtriage_twinning)
  test("xtriage_anomalous",
    test_xtriage_anomalous)
  test("xtriage_never_raises",
    test_xtriage_never_raises)
  test("xtriage_twinning_inferred",
    test_xtriage_twinning_inferred)
  test("xtriage_twinning_explicit_overrides",
    test_xtriage_twinning_explicit_overrides)
  test("xtriage_anomalous_inferred",
    test_xtriage_anomalous_inferred)
  print()

  # --- update_from_phaser ---
  print("update_from_phaser")
  test("phaser_basic",
    test_phaser_basic)
  test("phaser_space_group_update",
    test_phaser_space_group_update)
  test("phaser_never_raises",
    test_phaser_never_raises)
  print()

  # --- update_from_validation ---
  print("update_from_validation")
  test("validation_r_factors",
    test_validation_r_factors)
  test("validation_model_contents",
    test_validation_model_contents)
  test("validation_geometry",
    test_validation_geometry)
  test("validation_diff_peaks",
    test_validation_diff_peaks)
  test("validation_ligand_rscc",
    test_validation_ligand_rscc)
  test("validation_progress_tracking",
    test_validation_progress_tracking)
  test("validation_never_raises",
    test_validation_never_raises)
  print()

  # --- Problem detection ---
  print("Problem detection")
  test("problem_rfree_gap",
    test_problem_rfree_gap)
  test("problem_clashscore",
    test_problem_clashscore)
  test("problem_rama_outliers",
    test_problem_rama_outliers)
  test("problem_poor_ligand",
    test_problem_poor_ligand)
  test("problem_unmodeled_density",
    test_problem_unmodeled_density)
  test("problem_severity_ordering",
    test_problem_severity_ordering)
  test("problems_rebuild_each_cycle",
    test_problems_rebuild_each_cycle)
  test("problem_rfree_stalled",
    test_problem_rfree_stalled)
  test("problem_rfree_regression",
    test_problem_rfree_regression)
  test("problem_twinning",
    test_problem_twinning)
  test("problem_negative_peaks",
    test_problem_negative_peaks)
  print()

  # --- Multi-cycle accumulation ---
  print("Multi-cycle accumulation")
  test("state_persists_across_cycles",
    test_state_persists_across_cycles)
  print()

  # --- Runtime cache ---
  print("Runtime cache")
  test("kdtree_cache_slot",
    test_kdtree_cache_slot)
  print()

  # --- Chain completeness ---
  print("Chain completeness")
  test("update_chain_completeness",
    test_update_chain_completeness)
  test("chain_completeness_add_new",
    test_chain_completeness_add_new)
  test("chain_completeness_via_validation",
    test_chain_completeness_via_validation)
  print()

  # --- Strategy blacklist ---
  print("Strategy blacklist")
  test("blacklist_add",
    test_blacklist_add)
  test("blacklist_is_blacklisted",
    test_blacklist_is_blacklisted)
  test("blacklist_reason",
    test_blacklist_reason)
  test("blacklist_cap",
    test_blacklist_cap)
  print()

  # --- Hypothesis management ---
  print("Hypothesis management")
  test("add_hypothesis_basic",
    test_add_hypothesis_basic)
  test("single_active_budget",
    test_single_active_budget)
  test("add_resolved_no_conflict",
    test_add_resolved_no_conflict)
  test("get_active_hypothesis",
    test_get_active_hypothesis)
  test("get_hypothesis_by_id",
    test_get_hypothesis_by_id)
  test("hypothesis_cap",
    test_hypothesis_cap)
  print()

  # --- Enantiomorph tracking ---
  print("Enantiomorph tracking")
  test("rfree_trend_per_spacegroup",
    test_rfree_trend_per_spacegroup)
  test("rfree_trend_all",
    test_rfree_trend_all)
  print()

  # --- get_metric ---
  print("get_metric")
  test("get_metric_rfree",
    test_get_metric_rfree)
  test("get_metric_geometry",
    test_get_metric_geometry)
  test("get_metric_data_chars",
    test_get_metric_data_chars)
  test("get_metric_ligand_cc",
    test_get_metric_ligand_cc)
  test("get_metric_ligand_rscc_alias",
    test_get_metric_ligand_rscc_alias)
  test("get_metric_ligand_z_rscc",
    test_get_metric_ligand_z_rscc)
  test("get_metric_r_free_gap",
    test_get_metric_r_free_gap)
  test("get_metric_counts",
    test_get_metric_counts)
  test("get_metric_waters_is_int",
    test_get_metric_waters_is_int)
  test("get_metric_missing",
    test_get_metric_missing)
  print()

  # --- Summaries ---
  print("Summaries")
  test("summary_brief_empty",
    test_summary_brief_empty)
  test("summary_brief_populated",
    test_summary_brief_populated)
  test("summary_normal_populated",
    test_summary_normal_populated)
  test("summary_detailed_populated",
    test_summary_detailed_populated)
  test("summary_never_raises",
    test_summary_never_raises)
  print()

  # --- Serialization ---
  print("Serialization")
  test("roundtrip_empty",
    test_roundtrip_empty)
  test("roundtrip_populated",
    test_roundtrip_populated)
  test("from_dict_tolerant",
    test_from_dict_tolerant)
  test("from_dict_none",
    test_from_dict_none)
  test("json_serializable",
    test_json_serializable)
  test("roundtrip_preserves_hypotheses",
    test_roundtrip_preserves_hypotheses)
  test("fingerprint_changes",
    test_fingerprint_changes)
  print()

  # --- Annotate progress ---
  print("Annotate progress")
  test("annotate_last_progress",
    test_annotate_last_progress)
  print()

  # Summary
  print("=" * 60)
  total = passed + failed + skipped
  print("Results: %d/%d passed, %d failed, %d skipped"
    % (passed, total, failed, skipped))
  print("=" * 60)

  if failed > 0:
    sys.exit(1)


# ── Fixtures ──────────────────────────────────────────

def _make_validation_result():
  """Realistic validation_result dict."""
  return {
    "model_contents": {
      "chains": ["A", "B"],
      "residue_count": 480,
      "ligands": [
        {"name": "ATP", "chain": "A", "resid": 301,
         "n_atoms": 31},
      ],
      "waters": 187,
      "ions": [
        {"name": "MG", "chain": "A", "resid": 401},
      ],
      "has_hetatm": True,
    },
    "geometry": {
      "rama_favored": 0.968,
      "rama_outliers": 0.005,
      "rama_outlier_list": ["A/Arg78"],
      "rotamer_outliers": 0.012,
      "rotamer_outlier_list": ["A/Leu45"],
      "clashscore": 4.2,
      "bonds_rmsd": 0.007,
      "angles_rmsd": 1.1,
    },
    "data_model": {
      "ligand_rscc": [
        {"name": "ATP", "chain": "A", "resid": 301,
         "rscc": 0.87, "z_rscc": 1.2, "b_mean": 22.5,
         "occupancy": 1.0},
      ],
    },
    "diff_peaks": {
      "positive": [
        {"height": 5.2, "near_residue": "A/His47",
         "near_chain": "A", "distance": 2.1,
         "xyz": [10.0, 20.0, 30.0]},
      ],
      "negative": [
        {"height": -4.5, "near_residue": "B/Gly100",
         "near_chain": "B", "distance": 1.8,
         "xyz": [15.0, 25.0, 35.0]},
      ],
      "peak_count": 3,
    },
  }


def _make_log_metrics():
  """Realistic log_metrics dict."""
  return {
    "r_work": 0.210,
    "r_free": 0.248,
    "resolution": 2.1,
    "program": "phenix.refine",
  }


def _make_xtriage_results():
  """Realistic xtriage results dict."""
  return {
    "resolution": 2.1,
    "completeness": 0.98,
    "redundancy": 6.2,
    "space_group": "P 21 21 21",
    "unit_cell": [56.1, 72.3, 89.4, 90, 90, 90],
    "is_twinned": False,
    "twin_law": None,
    "twin_fraction": None,
    "has_anomalous_data": True,
    "anomalous_d_min": 3.0,
    "anomalous_measurability": 0.05,
    "i_over_sigma": 15.2,
    "data_quality": "good",
  }


def _make_phaser_results():
  """Realistic phaser results dict."""
  return {
    "tfz": 14.2,
    "llg": 890.0,
    "space_group": "P 21 21 21",
    "unit_cell": [56.1, 72.3, 89.4, 90, 90, 90],
    "n_copies": 2,
  }


def _make_populated_model():
  """Build a fully populated StructureModel."""
  sm = StructureModel()
  sm.update_from_xtriage(_make_xtriage_results())
  sm.update_from_phaser(_make_phaser_results())
  sm.update_from_validation(
    _make_validation_result(),
    _make_log_metrics(),
    cycle_number=4,
    program_name="phenix.refine",
  )
  return sm


# ── Hypothesis tests ──────────────────────────────────

def test_hypothesis_init_defaults():
  h = Hypothesis(id="h1", statement="test")
  assert h.id == "h1"
  assert h.statement == "test"
  assert h.status == "proposed"
  assert h.test_program == ""
  assert h.test_parameters == {}
  assert h.proposed_at_cycle == 0
  assert h.resolved_at_cycle is None
  assert h.test_cycles_remaining == 1
  assert h.revalidation_reason == ""


def test_hypothesis_status_validation():
  h = Hypothesis(id="h1", statement="test",
                 status="invalid_status")
  assert h.status == "proposed"  # falls back


def test_hypothesis_is_active():
  h1 = Hypothesis(id="h1", statement="test",
                  status="testing")
  assert h1.is_active is True
  h2 = Hypothesis(id="h2", statement="test",
                  status="pending")
  assert h2.is_active is True
  h3 = Hypothesis(id="h3", statement="test",
                  status="proposed")
  assert h3.is_active is False
  h4 = Hypothesis(id="h4", statement="test",
                  status="confirmed")
  assert h4.is_active is False


def test_hypothesis_is_resolved():
  for status in ("confirmed", "refuted", "abandoned"):
    h = Hypothesis(id="h1", statement="t",
                   status=status)
    assert h.is_resolved is True, (
      "status=%s should be resolved" % status
    )
  for status in ("proposed", "testing", "pending"):
    h = Hypothesis(id="h1", statement="t",
                   status=status)
    assert h.is_resolved is False, (
      "status=%s should not be resolved" % status
    )


def test_hypothesis_roundtrip():
  h = Hypothesis(
    id="h1",
    statement="Peak near His47 is Zn",
    test_program="phenix.refine",
    test_parameters={"add_ion": "Zn"},
    confirm_if="octahedral coordination",
    refute_if="B-factor > 50",
    status="testing",
    proposed_at_cycle=5,
    resolved_at_cycle=None,
    test_cycles_remaining=2,
    revalidation_reason="",
  )
  d = h.to_dict()
  h2 = Hypothesis.from_dict(d)
  assert h2.id == h.id
  assert h2.statement == h.statement
  assert h2.test_program == h.test_program
  assert h2.test_parameters == h.test_parameters
  assert h2.confirm_if == h.confirm_if
  assert h2.refute_if == h.refute_if
  assert h2.status == h.status
  assert h2.proposed_at_cycle == h.proposed_at_cycle
  assert h2.resolved_at_cycle == h.resolved_at_cycle
  assert h2.test_cycles_remaining == \
    h.test_cycles_remaining


def test_hypothesis_from_dict_tolerant():
  # Minimal dict — should not raise
  h = Hypothesis.from_dict({"id": "x"})
  assert h.id == "x"
  assert h.statement == "unknown"
  # Non-dict — should not raise
  h = Hypothesis.from_dict("garbage")
  assert h.id == "unknown"
  # Empty dict
  h = Hypothesis.from_dict({})
  assert h.id == "unknown"


# ── StructureModel init ──────────────────────────────

def test_init_defaults():
  sm = StructureModel()
  assert sm.data_characteristics["resolution"] is None
  assert sm.model_state["waters"] == 0
  assert sm.model_state["ligands"] == []
  assert sm.progress == []
  assert sm.strategy_blacklist == []
  assert sm.hypotheses == []
  geom = sm.model_state["geometry"]
  assert geom["clashscore"] is None
  dp = sm.model_state["diff_peaks"]
  assert dp["positive"] == []
  assert dp["peak_count"] == 0


def test_safe_float():
  assert _safe_float(3.14) == 3.14
  assert _safe_float("2.5") == 2.5
  assert _safe_float(None) is None
  assert _safe_float("bad") is None
  assert _safe_float([1]) is None


def test_deep_merge():
  target = {
    "a": 1,
    "b": {"x": 10, "y": 20},
    "c": 3,
  }
  source = {
    "a": 100,
    "b": {"x": 99, "z": 30},
    "d": 4,
  }
  _deep_merge(target, source)
  assert target["a"] == 100
  assert target["b"]["x"] == 99
  assert target["b"]["y"] == 20  # preserved
  assert target["b"]["z"] == 30  # added
  assert target["c"] == 3        # untouched
  assert target["d"] == 4        # added


# ── update_from_xtriage ──────────────────────────────

def test_xtriage_basic():
  sm = StructureModel()
  sm.update_from_xtriage(_make_xtriage_results())
  dc = sm.data_characteristics
  assert dc["resolution"] == 2.1
  assert dc["completeness"] == 0.98
  assert dc["redundancy"] == 6.2
  assert dc["space_group"] == "P 21 21 21"
  assert dc["unit_cell"] == (
    56.1, 72.3, 89.4, 90, 90, 90
  )
  assert dc["i_over_sigma"] == 15.2
  assert dc["data_quality"] == "good"


def test_xtriage_twinning():
  sm = StructureModel()
  sm.update_from_xtriage({
    "is_twinned": True,
    "twin_law": "h,-k,-l",
    "twin_fraction": 0.45,
  })
  tw = sm.data_characteristics["twinning"]
  assert tw["is_twinned"] is True
  assert tw["twin_law"] == "h,-k,-l"
  assert tw["twin_fraction"] == 0.45


def test_xtriage_anomalous():
  sm = StructureModel()
  sm.update_from_xtriage({
    "has_anomalous_data": False,
    "anomalous_d_min": None,
    "anomalous_measurability": 0.0,
  })
  anom = sm.data_characteristics["anomalous"]
  assert anom["has_anomalous_data"] is False
  assert anom["anomalous_d_min"] is None
  assert anom["anomalous_measurability"] == 0.0


def test_xtriage_never_raises():
  sm = StructureModel()
  # None input
  sm.update_from_xtriage(None)
  # String input
  sm.update_from_xtriage("garbage")
  # Bad value types
  sm.update_from_xtriage({"resolution": "bad"})
  assert sm.data_characteristics["resolution"] is None


def test_xtriage_twinning_inferred():
  """is_twinned inferred from twin_fraction when not
  explicitly provided (covers _extract_xtriage path)."""
  sm = StructureModel()
  # Only twin_fraction, no is_twinned key
  sm.update_from_xtriage({
    "twin_fraction": 0.45,
    "twin_law": "h,-k,-l",
  })
  tw = sm.data_characteristics["twinning"]
  assert tw["is_twinned"] is True
  assert tw["twin_fraction"] == 0.45
  # Low fraction should NOT flag as twinned
  sm2 = StructureModel()
  sm2.update_from_xtriage({"twin_fraction": 0.02})
  assert sm2.data_characteristics[
    "twinning"
  ]["is_twinned"] is False


def test_xtriage_twinning_explicit_overrides():
  """Explicit is_twinned=False should NOT be overridden
  by twin_fraction inference."""
  sm = StructureModel()
  sm.update_from_xtriage({
    "is_twinned": False,
    "twin_fraction": 0.45,
  })
  tw = sm.data_characteristics["twinning"]
  # Explicit False takes precedence
  assert tw["is_twinned"] is False


def test_xtriage_anomalous_inferred():
  """has_anomalous_data inferred from measurability."""
  sm = StructureModel()
  sm.update_from_xtriage({
    "anomalous_measurability": 0.08,
  })
  anom = sm.data_characteristics["anomalous"]
  assert anom["has_anomalous_data"] is True
  # Low measurability → no signal
  sm2 = StructureModel()
  sm2.update_from_xtriage({
    "anomalous_measurability": 0.02,
  })
  assert sm2.data_characteristics[
    "anomalous"
  ]["has_anomalous_data"] is False


# ── update_from_phaser ───────────────────────────────

def test_phaser_basic():
  sm = StructureModel()
  sm.update_from_phaser(_make_phaser_results())
  dc = sm.data_characteristics
  assert dc["mr_tfz"] == 14.2
  assert dc["mr_llg"] == 890.0
  assert dc["n_copies_asu"] == 2


def test_phaser_space_group_update():
  sm = StructureModel()
  sm.update_from_xtriage({
    "space_group": "P 3 2 1"
  })
  # Phaser confirms a different space group
  sm.update_from_phaser({
    "space_group": "P 32 2 1"
  })
  assert sm.data_characteristics["space_group"] == (
    "P 32 2 1"
  )


def test_phaser_never_raises():
  sm = StructureModel()
  sm.update_from_phaser(None)
  sm.update_from_phaser("garbage")
  sm.update_from_phaser({"tfz": "not_a_number"})


# ── update_from_validation ───────────────────────────

def test_validation_r_factors():
  sm = StructureModel()
  sm.update_from_validation(
    None, {"r_work": 0.21, "r_free": 0.25},
    cycle_number=1, program_name="phenix.refine",
  )
  assert sm.model_state["r_work"] == 0.21
  assert sm.model_state["r_free"] == 0.25


def test_validation_model_contents():
  sm = StructureModel()
  sm.update_from_validation(
    _make_validation_result(),
    _make_log_metrics(),
    cycle_number=1,
    program_name="phenix.refine",
  )
  ms = sm.model_state
  assert len(ms["chains"]) == 2
  assert ms["chains"][0]["chain_id"] == "A"
  assert ms["waters"] == 187
  assert len(ms["ligands"]) == 1
  assert ms["ligands"][0]["name"] == "ATP"
  assert len(ms["ions"]) == 1
  assert ms["ions"][0]["name"] == "MG"


def test_validation_geometry():
  sm = StructureModel()
  sm.update_from_validation(
    _make_validation_result(),
    _make_log_metrics(),
    cycle_number=1,
    program_name="phenix.refine",
  )
  geom = sm.model_state["geometry"]
  assert geom["rama_favored"] == 0.968
  assert geom["clashscore"] == 4.2
  assert geom["bonds_rmsd"] == 0.007


def test_validation_diff_peaks():
  sm = StructureModel()
  sm.update_from_validation(
    _make_validation_result(),
    _make_log_metrics(),
    cycle_number=1,
    program_name="phenix.refine",
  )
  dp = sm.model_state["diff_peaks"]
  assert len(dp["positive"]) == 1
  assert dp["positive"][0]["height"] == 5.2
  assert len(dp["negative"]) == 1
  assert dp["peak_count"] == 3


def test_validation_ligand_rscc():
  sm = StructureModel()
  sm.update_from_validation(
    _make_validation_result(),
    _make_log_metrics(),
    cycle_number=1,
    program_name="phenix.refine",
  )
  ligs = sm.model_state["ligands"]
  assert len(ligs) == 1
  assert ligs[0]["rscc"] == 0.87
  assert ligs[0]["z_rscc"] == 1.2
  assert ligs[0]["b_mean"] == 22.5


def test_validation_progress_tracking():
  sm = StructureModel()
  sm.update_from_xtriage({"space_group": "P 21 21 21"})
  sm.update_from_validation(
    _make_validation_result(),
    {"r_work": 0.21, "r_free": 0.25},
    cycle_number=3,
    program_name="phenix.refine",
  )
  assert len(sm.progress) == 1
  entry = sm.progress[0]
  assert entry["cycle"] == 3
  assert entry["program"] == "phenix.refine"
  assert entry["r_free"] == 0.25
  assert entry["space_group"] == "P 21 21 21"


def test_validation_never_raises():
  sm = StructureModel()
  # All None inputs
  sm.update_from_validation(None, None, 1, "test")
  # Bad types
  sm.update_from_validation("bad", "bad", "bad", "bad")
  # Partial validation result
  sm.update_from_validation(
    {"model_contents": None},
    {},
    cycle_number=1,
    program_name="test",
  )


# ── Problem detection ────────────────────────────────

def test_problem_rfree_gap():
  sm = StructureModel()
  sm.update_from_validation(
    None,
    {"r_work": 0.20, "r_free": 0.35},
    cycle_number=1,
    program_name="phenix.refine",
  )
  problems = sm.model_state["problems"]
  gap_probs = [
    p for p in problems
    if p["type"] == "r_free_gap"
  ]
  assert len(gap_probs) == 1
  assert "overfitting" in gap_probs[0]["description"]
  assert gap_probs[0]["severity"] == "high"


def test_problem_clashscore():
  sm = StructureModel()
  vr = {"geometry": {"clashscore": 25.0}}
  sm.update_from_validation(
    vr, {}, cycle_number=1, program_name="test",
  )
  problems = sm.model_state["problems"]
  clash_probs = [
    p for p in problems
    if p["type"] == "clashscore"
  ]
  assert len(clash_probs) == 1
  assert clash_probs[0]["severity"] == "high"


def test_problem_rama_outliers():
  sm = StructureModel()
  vr = {
    "geometry": {
      "rama_outlier_list": [
        "A/Arg78", "A/Leu45", "B/Pro100", "B/Gly55",
      ],
    },
  }
  sm.update_from_validation(
    vr, {}, cycle_number=1, program_name="test",
  )
  problems = sm.model_state["problems"]
  rama_probs = [
    p for p in problems
    if p["type"] == "rama_outlier"
  ]
  assert len(rama_probs) == 1
  assert rama_probs[0]["severity"] == "high"
  assert "4" in rama_probs[0]["description"]


def test_problem_poor_ligand():
  sm = StructureModel()
  vr = {
    "model_contents": {
      "chains": ["A"],
      "ligands": [
        {"name": "ATP", "chain": "A",
         "resid": 301, "n_atoms": 31},
      ],
      "waters": 0,
      "ions": [],
    },
    "data_model": {
      "ligand_rscc": [
        {"name": "ATP", "chain": "A", "resid": 301,
         "rscc": 0.45},
      ],
    },
  }
  sm.update_from_validation(
    vr, {}, cycle_number=1, program_name="test",
  )
  problems = sm.model_state["problems"]
  lig_probs = [
    p for p in problems
    if p["type"] == "poor_ligand_fit"
  ]
  assert len(lig_probs) == 1
  assert lig_probs[0]["severity"] == "high"  # <0.5


def test_problem_unmodeled_density():
  sm = StructureModel()
  vr = {
    "diff_peaks": {
      "positive": [
        {"height": 8.5, "near_residue": "A/His47"},
      ],
      "negative": [],
      "peak_count": 1,
    },
  }
  sm.update_from_validation(
    vr, {}, cycle_number=1, program_name="test",
  )
  problems = sm.model_state["problems"]
  density_probs = [
    p for p in problems
    if p["type"] == "unmodeled_density"
  ]
  assert len(density_probs) == 1
  assert density_probs[0]["severity"] == "high"
  assert "8.5" in density_probs[0]["description"]


def test_problem_severity_ordering():
  sm = StructureModel()
  vr = {
    "geometry": {
      "clashscore": 12.0,  # medium
      "rama_outlier_list": ["A/Arg78"],  # medium (1)
    },
  }
  sm.update_from_validation(
    vr,
    {"r_work": 0.20, "r_free": 0.35},  # high gap
    cycle_number=1,
    program_name="test",
  )
  problems = sm.model_state["problems"]
  assert len(problems) >= 2
  # High severity should come first
  assert problems[0]["severity"] == "high"


def test_problems_rebuild_each_cycle():
  """Problems are rebuilt from scratch, not accumulated."""
  sm = StructureModel()
  # Cycle 1: high clashscore
  vr1 = {"geometry": {"clashscore": 25.0}}
  sm.update_from_validation(
    vr1, {}, cycle_number=1, program_name="test",
  )
  assert any(
    p["type"] == "clashscore"
    for p in sm.model_state["problems"]
  )
  # Cycle 2: clashscore fixed
  vr2 = {"geometry": {"clashscore": 3.0}}
  sm.update_from_validation(
    vr2, {}, cycle_number=2, program_name="test",
  )
  # Clashscore problem should be gone
  assert not any(
    p["type"] == "clashscore"
    for p in sm.model_state["problems"]
  )


# ── Chain completeness ───────────────────────────────

def test_update_chain_completeness():
  sm = StructureModel()
  # First set up chains from validation
  vr = {
    "model_contents": {
      "chains": ["A", "B"],
      "ligands": [],
      "waters": 0,
      "ions": [],
    },
  }
  sm.update_from_validation(
    vr, {}, cycle_number=1, program_name="test",
  )
  # Now add completeness data
  sm.update_chain_completeness([
    {"chain_id": "A", "residues_built": 270,
     "residues_expected": 280,
     "completeness_fraction": 0.964,
     "gaps": [{"start": 145, "end": 160}],
     "avg_b_factor": 25.0,
     "disordered_fraction": 0.05},
    {"chain_id": "B", "residues_built": 240,
     "residues_expected": 280,
     "completeness_fraction": 0.857,
     "gaps": [],
     "avg_b_factor": 35.0,
     "disordered_fraction": 0.15},
  ])
  chains = sm.model_state["chains"]
  assert chains[0]["completeness"] == 0.964
  assert chains[0]["residues_built"] == 270
  assert len(chains[0]["gaps"]) == 1
  assert chains[1]["completeness"] == 0.857


def test_chain_completeness_add_new():
  sm = StructureModel()
  # No existing chains
  sm.update_chain_completeness([
    {"chain_id": "C", "residues_built": 100,
     "residues_expected": 100,
     "completeness_fraction": 1.0,
     "gaps": [],
     "avg_b_factor": 20.0,
     "disordered_fraction": 0.0},
  ])
  chains = sm.model_state["chains"]
  assert len(chains) == 1
  assert chains[0]["chain_id"] == "C"
  assert chains[0]["completeness"] == 1.0


def test_chain_completeness_via_validation():
  """chain_completeness in validation_result is
  auto-consumed by update_from_validation."""
  sm = StructureModel()
  vr = {
    "model_contents": {
      "chains": ["A", "B"],
      "ligands": [],
      "waters": 50,
      "ions": [],
    },
    "chain_completeness": [
      {"chain_id": "A", "residues_built": 270,
       "residues_expected": 280,
       "completeness_fraction": 0.964,
       "gaps": [{"start": 145, "end": 160}],
       "avg_b_factor": 25.0,
       "disordered_fraction": 0.05},
      {"chain_id": "B", "residues_built": 200,
       "residues_expected": 280,
       "completeness_fraction": 0.714,
       "gaps": [],
       "avg_b_factor": 40.0,
       "disordered_fraction": 0.20},
    ],
  }
  sm.update_from_validation(
    vr, {"r_free": 0.30},
    cycle_number=1, program_name="phenix.refine",
  )
  chains = sm.model_state["chains"]
  assert len(chains) == 2
  assert chains[0]["chain_id"] == "A"
  assert chains[0]["completeness"] == 0.964
  assert chains[0]["residues_built"] == 270
  assert len(chains[0]["gaps"]) == 1
  assert chains[1]["chain_id"] == "B"
  assert chains[1]["completeness"] == 0.714
  assert chains[1]["disordered_fraction"] == 0.20


# ── Strategy blacklist ───────────────────────────────

def test_blacklist_add():
  sm = StructureModel()
  sm.blacklist_strategy(
    "MR_beta.pdb", "R-free stuck at 0.48",
    {"r_free": 0.48, "cycle": 5},
  )
  assert len(sm.strategy_blacklist) == 1
  entry = sm.strategy_blacklist[0]
  assert entry["strategy_id"] == "MR_beta.pdb"
  assert "stuck" in entry["reason"]
  assert entry["metrics_at_retreat"]["r_free"] == 0.48


def test_blacklist_is_blacklisted():
  sm = StructureModel()
  sm.blacklist_strategy("MR_beta.pdb", "failed")
  assert sm.is_blacklisted("MR_beta.pdb") is True
  assert sm.is_blacklisted("MR_alpha.pdb") is False


def test_blacklist_reason():
  sm = StructureModel()
  sm.blacklist_strategy("MR_beta.pdb", "stuck at 0.48")
  reason = sm.get_blacklist_reason("MR_beta.pdb")
  assert "stuck" in reason
  assert sm.get_blacklist_reason("other") is None


def test_blacklist_cap():
  sm = StructureModel()
  for i in range(25):
    sm.blacklist_strategy("strat_%d" % i, "reason")
  assert len(sm.strategy_blacklist) <= 20


# ── Hypothesis management ────────────────────────────

def test_add_hypothesis_basic():
  sm = StructureModel()
  h = Hypothesis(id="h1", statement="test",
                 status="proposed")
  assert sm.add_hypothesis(h) is True
  assert len(sm.hypotheses) == 1


def test_single_active_budget():
  sm = StructureModel()
  h1 = Hypothesis(id="h1", statement="first",
                  status="testing")
  assert sm.add_hypothesis(h1) is True
  # Second active hypothesis should be rejected
  h2 = Hypothesis(id="h2", statement="second",
                  status="testing")
  assert sm.add_hypothesis(h2) is False
  assert len(sm.hypotheses) == 1
  # Pending also counts as active
  h3 = Hypothesis(id="h3", statement="third",
                  status="pending")
  assert sm.add_hypothesis(h3) is False


def test_add_resolved_no_conflict():
  sm = StructureModel()
  # Active hypothesis
  h1 = Hypothesis(id="h1", statement="active",
                  status="testing")
  sm.add_hypothesis(h1)
  # Resolved hypothesis should be accepted
  h2 = Hypothesis(id="h2", statement="done",
                  status="confirmed")
  assert sm.add_hypothesis(h2) is True
  assert len(sm.hypotheses) == 2


def test_get_active_hypothesis():
  sm = StructureModel()
  assert sm.get_active_hypothesis() is None
  h = Hypothesis(id="h1", statement="test",
                 status="testing")
  sm.add_hypothesis(h)
  active = sm.get_active_hypothesis()
  assert active is not None
  assert active.id == "h1"


def test_get_hypothesis_by_id():
  sm = StructureModel()
  h = Hypothesis(id="h1", statement="test")
  sm.add_hypothesis(h)
  found = sm.get_hypothesis("h1")
  assert found is not None
  assert found.statement == "test"
  assert sm.get_hypothesis("h999") is None


def test_hypothesis_cap():
  sm = StructureModel()
  for i in range(35):
    h = Hypothesis(
      id="h%d" % i, statement="test %d" % i,
      status="confirmed",
      resolved_at_cycle=i,
    )
    sm.add_hypothesis(h)
  assert len(sm.hypotheses) <= 30


# ── Enantiomorph tracking ────────────────────────────

def test_rfree_trend_per_spacegroup():
  sm = StructureModel()
  # Run 3 cycles in P3121
  sm.data_characteristics["space_group"] = "P 31 21"
  for i in range(1, 4):
    sm.update_from_validation(
      None,
      {"r_free": 0.50 - i * 0.01},
      cycle_number=i,
      program_name="phenix.refine",
    )
  # Switch to P3221
  sm.data_characteristics["space_group"] = "P 32 21"
  for i in range(4, 7):
    sm.update_from_validation(
      None,
      {"r_free": 0.45 - (i - 3) * 0.05},
      cycle_number=i,
      program_name="phenix.refine",
    )
  # Check per-SG trends
  trend_31 = sm.get_rfree_trend(
    space_group="P 31 21"
  )
  assert len(trend_31) == 3
  assert trend_31[0] == (1, 0.49)
  assert trend_31[2] == (3, 0.47)

  trend_32 = sm.get_rfree_trend(
    space_group="P 32 21"
  )
  assert len(trend_32) == 3
  assert trend_32[0] == (4, 0.40)
  assert trend_32[2] == (6, 0.30)


def test_rfree_trend_all():
  sm = StructureModel()
  sm.update_from_validation(
    None, {"r_free": 0.48},
    cycle_number=1, program_name="test",
  )
  sm.update_from_validation(
    None, {"r_free": 0.35},
    cycle_number=2, program_name="test",
  )
  trend = sm.get_rfree_trend()  # no filter
  assert len(trend) == 2
  assert trend[0] == (1, 0.48)
  assert trend[1] == (2, 0.35)


# ── get_metric ───────────────────────────────────────

def test_get_metric_rfree():
  sm = _make_populated_model()
  assert sm.get_metric("r_free") == 0.248
  assert sm.get_metric("r_work") == 0.210


def test_get_metric_geometry():
  sm = _make_populated_model()
  assert sm.get_metric("clashscore") == 4.2
  assert abs(
    sm.get_metric("rama_favored") - 0.968
  ) < 0.001


def test_get_metric_data_chars():
  sm = _make_populated_model()
  assert sm.get_metric("resolution") == 2.1
  assert sm.get_metric("tfz") == 14.2
  assert sm.get_metric("llg") == 890.0


def test_get_metric_ligand_cc():
  sm = _make_populated_model()
  assert sm.get_metric("ligand_cc") == 0.87


def test_get_metric_missing():
  sm = StructureModel()
  assert sm.get_metric("r_free") is None
  assert sm.get_metric("nonexistent") is None


# ── Summaries ────────────────────────────────────────

def test_summary_brief_empty():
  sm = StructureModel()
  assert sm.get_summary("brief") == ""


def test_summary_brief_populated():
  sm = _make_populated_model()
  brief = sm.get_summary("brief")
  assert "chain" in brief.lower() or "Chain" in brief
  assert "ATP" in brief
  assert "187 waters" in brief


def test_summary_normal_populated():
  sm = _make_populated_model()
  normal = sm.get_summary("normal")
  assert "R-work" in normal or "R-free" in normal
  assert "2.1" in normal  # resolution
  assert "ATP" in normal
  assert "Clash" in normal or "clash" in normal


def test_summary_detailed_populated():
  sm = _make_populated_model()
  detailed = sm.get_summary("detailed")
  assert "P 21 21 21" in detailed
  assert "Chain A" in detailed
  assert "Ligand: ATP" in detailed or "ATP" in detailed
  assert "Waters:" in detailed
  assert "Clashscore" in detailed


def test_summary_never_raises():
  sm = StructureModel()
  # Corrupt internal state
  sm.model_state = "garbage"
  result = sm.get_summary("normal")
  assert isinstance(result, str)


# ── Serialization ────────────────────────────────────

def test_roundtrip_empty():
  sm = StructureModel()
  d = sm.to_dict()
  sm2 = StructureModel.from_dict(d)
  assert sm2.data_characteristics["resolution"] is None
  assert sm2.model_state["waters"] == 0
  assert sm2.progress == []
  assert sm2.strategy_blacklist == []
  assert sm2.hypotheses == []


def test_roundtrip_populated():
  sm = _make_populated_model()
  sm.blacklist_strategy("strat1", "failed")
  sm.annotate_last_progress("added 187 waters")

  d = sm.to_dict()
  sm2 = StructureModel.from_dict(d)

  # Data characteristics
  assert sm2.data_characteristics["resolution"] == 2.1
  assert sm2.data_characteristics["space_group"] == (
    "P 21 21 21"
  )
  # Model state
  assert sm2.model_state["r_free"] == 0.248
  assert len(sm2.model_state["chains"]) == 2
  assert len(sm2.model_state["ligands"]) == 1
  assert sm2.model_state["waters"] == 187
  # Progress
  assert len(sm2.progress) == 1
  assert sm2.progress[0]["annotation"] == (
    "added 187 waters"
  )
  # Blacklist
  assert len(sm2.strategy_blacklist) == 1
  assert sm2.is_blacklisted("strat1")


def test_from_dict_tolerant():
  """Should handle partial/extra keys."""
  d = {
    "data_characteristics": {
      "resolution": 3.0,
      "extra_key": "ignored_by_merge",
    },
    "model_state": {
      "r_free": 0.30,
    },
    # Missing progress, blacklist, hypotheses
  }
  sm = StructureModel.from_dict(d)
  assert sm.data_characteristics["resolution"] == 3.0
  assert sm.model_state["r_free"] == 0.30
  assert sm.progress == []
  assert sm.hypotheses == []


def test_from_dict_none():
  sm = StructureModel.from_dict(None)
  assert sm.data_characteristics["resolution"] is None


def test_json_serializable():
  """to_dict() output must be fully JSON-serializable."""
  sm = _make_populated_model()
  sm.blacklist_strategy("strat", "reason")
  h = Hypothesis(id="h1", statement="test",
                 status="confirmed",
                 resolved_at_cycle=5)
  sm.add_hypothesis(h)
  d = sm.to_dict()
  # This will raise if not JSON-serializable
  json_str = json.dumps(d)
  assert len(json_str) > 0
  # Verify round-trip through JSON
  d2 = json.loads(json_str)
  sm2 = StructureModel.from_dict(d2)
  assert sm2.model_state["r_free"] == 0.248


def test_roundtrip_preserves_hypotheses():
  sm = StructureModel()
  h1 = Hypothesis(
    id="h1", statement="Peak is Zn",
    status="confirmed", resolved_at_cycle=8,
    test_program="phenix.refine",
    test_parameters={"add_ion": "Zn"},
  )
  h2 = Hypothesis(
    id="h2", statement="Alt conf at Arg78",
    status="refuted", resolved_at_cycle=10,
  )
  sm.add_hypothesis(h1)
  sm.add_hypothesis(h2)
  d = sm.to_dict()
  # Through JSON
  d2 = json.loads(json.dumps(d))
  sm2 = StructureModel.from_dict(d2)
  assert len(sm2.hypotheses) == 2
  assert sm2.hypotheses[0].id == "h1"
  assert sm2.hypotheses[0].status == "confirmed"
  assert sm2.hypotheses[0].test_parameters == (
    {"add_ion": "Zn"}
  )
  assert sm2.hypotheses[1].id == "h2"
  assert sm2.hypotheses[1].status == "refuted"


def test_fingerprint_changes():
  sm = StructureModel()
  fp1 = sm.compute_fingerprint()
  assert isinstance(fp1, str)
  assert len(fp1) > 0
  # Update R-free — fingerprint should change
  sm.model_state["r_free"] = 0.30
  fp2 = sm.compute_fingerprint()
  assert fp2 != fp1
  # Blacklist a strategy — fingerprint should change
  sm.blacklist_strategy("strat", "reason")
  fp3 = sm.compute_fingerprint()
  assert fp3 != fp2


# ── Annotate progress ────────────────────────────────

def test_annotate_last_progress():
  sm = StructureModel()
  sm.update_from_validation(
    None, {"r_free": 0.30},
    cycle_number=1, program_name="test",
  )
  sm.annotate_last_progress("added 100 waters")
  assert sm.progress[-1]["annotation"] == (
    "added 100 waters"
  )
  # Annotate on empty progress should not raise
  sm2 = StructureModel()
  sm2.annotate_last_progress("no-op")


# ── New problem detection tests ─────────────────────

def test_problem_rfree_stalled():
  """R-free stalled for 3+ cycles is detected."""
  sm = StructureModel()
  # 3 cycles with nearly identical R-free
  for i in range(1, 4):
    sm.update_from_validation(
      None,
      {"r_free": 0.350 + i * 0.0005},
      cycle_number=i,
      program_name="phenix.refine",
    )
  problems = sm.model_state["problems"]
  stall = [
    p for p in problems
    if p["type"] == "r_free_stalled"
  ]
  assert len(stall) == 1, (
    "Expected r_free_stalled problem, got: %s"
    % [p["type"] for p in problems]
  )
  assert stall[0]["severity"] == "high"


def test_problem_rfree_regression():
  """R-free increase > 0.01 is flagged."""
  sm = StructureModel()
  sm.update_from_validation(
    None, {"r_free": 0.30},
    cycle_number=1, program_name="test",
  )
  sm.update_from_validation(
    None, {"r_free": 0.32},
    cycle_number=2, program_name="test",
  )
  problems = sm.model_state["problems"]
  reg = [
    p for p in problems
    if p["type"] == "r_free_regression"
  ]
  assert len(reg) == 1, (
    "Expected r_free_regression problem, got: %s"
    % [p["type"] for p in problems]
  )
  assert "0.300" in reg[0]["description"]
  assert "0.320" in reg[0]["description"]


def test_problem_twinning():
  """Twinned data is flagged as a problem."""
  sm = StructureModel()
  sm.update_from_xtriage({
    "is_twinned": True,
    "twin_fraction": 0.42,
  })
  # Trigger problem detection via validation
  sm.update_from_validation(
    None, {"r_free": 0.40},
    cycle_number=1, program_name="test",
  )
  problems = sm.model_state["problems"]
  tw = [
    p for p in problems
    if p["type"] == "twinning"
  ]
  assert len(tw) == 1, (
    "Expected twinning problem, got: %s"
    % [p["type"] for p in problems]
  )
  assert tw[0]["severity"] == "high"
  assert "0.42" in tw[0]["description"]


def test_problem_negative_peaks():
  """Deep negative holes flagged as misplaced atoms."""
  sm = StructureModel()
  vr = {
    "diff_peaks": {
      "positive": [],
      "negative": [
        {"height": -6.5, "near_residue": "B/Gly100",
         "near_chain": "B", "distance": 1.8},
      ],
      "peak_count": 0,
    },
  }
  sm.update_from_validation(
    vr, {}, cycle_number=1, program_name="test",
  )
  problems = sm.model_state["problems"]
  misplaced = [
    p for p in problems
    if p["type"] == "misplaced_atoms"
  ]
  assert len(misplaced) == 1, (
    "Expected misplaced_atoms problem, got: %s"
    % [p["type"] for p in problems]
  )
  assert misplaced[0]["severity"] == "high"
  assert "B/Gly100" in misplaced[0]["description"]


def test_state_persists_across_cycles():
  """Chain/ligand data from cycle N survives when
  cycle N+1 validation has only R-factors."""
  sm = StructureModel()
  # Cycle 1: full validation with model contents
  vr1 = _make_validation_result()
  sm.update_from_validation(
    vr1, {"r_work": 0.25, "r_free": 0.30},
    cycle_number=1, program_name="phenix.refine",
  )
  assert len(sm.model_state["chains"]) == 2
  assert len(sm.model_state["ligands"]) == 1
  assert sm.model_state["ligands"][0]["rscc"] == 0.87
  assert sm.model_state["waters"] == 187
  # Cycle 2: only R-factors (validation returned None)
  sm.update_from_validation(
    None, {"r_work": 0.22, "r_free": 0.27},
    cycle_number=2, program_name="phenix.refine",
  )
  # R-factors updated
  assert sm.model_state["r_free"] == 0.27
  assert sm.model_state["r_work"] == 0.22
  # But chains/ligands/waters PERSIST from cycle 1
  assert len(sm.model_state["chains"]) == 2
  assert len(sm.model_state["ligands"]) == 1
  assert sm.model_state["ligands"][0]["rscc"] == 0.87
  assert sm.model_state["waters"] == 187
  # And progress has both cycles
  assert len(sm.progress) == 2
  assert sm.progress[0]["r_free"] == 0.30
  assert sm.progress[1]["r_free"] == 0.27


def test_kdtree_cache_slot():
  """KD-tree cache attributes exist and are not
  serialized."""
  sm = StructureModel()
  assert sm._kdtree_cache is None
  assert sm._kdtree_cache_key is None
  # Set a fake cache value
  sm._kdtree_cache = "fake_tree"
  sm._kdtree_cache_key = ("model.pdb", 123, 456)
  # Verify not serialized
  d = sm.to_dict()
  assert "_kdtree_cache" not in d
  assert "_kdtree_cache_key" not in d
  # Verify fresh after round-trip
  sm2 = StructureModel.from_dict(d)
  assert sm2._kdtree_cache is None
  assert sm2._kdtree_cache_key is None


# ── New get_metric tests ────────────────────────────

def test_get_metric_ligand_rscc_alias():
  """ligand_rscc is an alias for ligand_cc."""
  sm = _make_populated_model()
  cc = sm.get_metric("ligand_cc")
  rscc = sm.get_metric("ligand_rscc")
  assert cc == rscc
  assert cc == 0.87


def test_get_metric_ligand_z_rscc():
  """ligand_z_rscc returns worst Z-score."""
  sm = _make_populated_model()
  z = sm.get_metric("ligand_z_rscc")
  assert z == 1.2
  # No z_rscc data -> None
  sm2 = StructureModel()
  assert sm2.get_metric("ligand_z_rscc") is None


def test_get_metric_r_free_gap():
  """r_free_gap returns r_free - r_work."""
  sm = _make_populated_model()
  gap = sm.get_metric("r_free_gap")
  assert gap is not None
  expected = 0.248 - 0.210
  assert abs(gap - expected) < 0.001
  # No R-factors -> None
  sm2 = StructureModel()
  assert sm2.get_metric("r_free_gap") is None


def test_get_metric_counts():
  """Count metrics return correct values."""
  sm = _make_populated_model()
  assert sm.get_metric("ligand_count") == 1
  assert sm.get_metric("ion_count") == 1
  assert sm.get_metric("water_count") == 187
  assert sm.get_metric("positive_peak_count") == 1
  assert sm.get_metric("negative_peak_count") == 1
  # Empty model
  sm2 = StructureModel()
  assert sm2.get_metric("ligand_count") == 0
  assert sm2.get_metric("positive_peak_count") == 0


def test_get_metric_waters_is_int():
  """waters and water_count return int, not float."""
  sm = _make_populated_model()
  w = sm.get_metric("waters")
  assert isinstance(w, int), (
    "Expected int, got %s" % type(w).__name__
  )
  assert w == 187
  wc = sm.get_metric("water_count")
  assert isinstance(wc, int)
  assert wc == 187


# ── Entry point ──────────────────────────────────────

if __name__ == "__main__":
  run_tests()
