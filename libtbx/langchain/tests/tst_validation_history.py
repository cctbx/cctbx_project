"""
Unit tests for ValidationHistory (Phase 1, Step 1.6).

Run standalone:
  python tests/tst_validation_history.py

Run from PHENIX:
  libtbx.python tests/tst_validation_history.py

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

from agent.validation_history import (
  ValidationHistory,
  _extract_metric,
  _extract_all_metrics,
  _metric_direction,
  _safe_float,
)


def run_tests():
  """Run all ValidationHistory tests."""
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
  print("ValidationHistory Unit Tests")
  print("=" * 60)
  print()

  # --- Basics ---
  print("Basics")
  test("init_empty",
    test_init_empty)
  test("safe_float_helper",
    test_safe_float_helper)
  test("metric_direction",
    test_metric_direction)
  print()

  # --- Recording ---
  print("Recording")
  test("record_basic",
    test_record_basic)
  test("record_none_validation",
    test_record_none_validation)
  test("record_idempotent",
    test_record_idempotent)
  test("record_cap",
    test_record_cap)
  test("record_never_raises",
    test_record_never_raises)
  print()

  # --- Retrieval ---
  print("Retrieval")
  test("get_at_cycle",
    test_get_at_cycle)
  test("get_at_cycle_missing",
    test_get_at_cycle_missing)
  test("get_latest",
    test_get_latest)
  test("get_latest_empty",
    test_get_latest_empty)
  test("get_previous",
    test_get_previous)
  test("get_previous_by_cycle",
    test_get_previous_by_cycle)
  test("get_previous_insufficient",
    test_get_previous_insufficient)
  print()

  # --- Metric extraction ---
  print("Metric extraction")
  test("extract_r_factors",
    test_extract_r_factors)
  test("extract_geometry",
    test_extract_geometry)
  test("extract_ligand_cc",
    test_extract_ligand_cc)
  test("extract_peak_count",
    test_extract_peak_count)
  test("extract_water_count",
    test_extract_water_count)
  test("extract_residue_count",
    test_extract_residue_count)
  test("extract_missing",
    test_extract_missing)
  test("extract_all_metrics",
    test_extract_all_metrics)
  print()

  # --- Metric series ---
  print("Metric series")
  test("series_rfree",
    test_series_rfree)
  test("series_clashscore",
    test_series_clashscore)
  test("series_missing_gaps",
    test_series_missing_gaps)
  test("series_space_group_filter",
    test_series_space_group_filter)
  test("series_never_raises",
    test_series_never_raises)
  print()

  # --- Phase start metrics ---
  print("Phase start metrics")
  test("phase_start_exact",
    test_phase_start_exact)
  test("phase_start_closest",
    test_phase_start_closest)
  test("phase_start_missing",
    test_phase_start_missing)
  print()

  # --- Flat metrics ---
  print("Flat metrics")
  test("flat_metrics_populated",
    test_flat_metrics_populated)
  test("flat_metrics_missing",
    test_flat_metrics_missing)
  print()

  # --- Delta computation ---
  print("Delta computation")
  test("metrics_delta",
    test_metrics_delta)
  test("metrics_delta_missing_cycle",
    test_metrics_delta_missing_cycle)
  test("metrics_delta_partial_overlap",
    test_metrics_delta_partial_overlap)
  print()

  # --- is_improving ---
  print("is_improving")
  test("improving_rfree",
    test_improving_rfree)
  test("improving_rama_favored",
    test_improving_rama_favored)
  test("improving_insufficient",
    test_improving_insufficient)
  test("not_improving_stalled",
    test_not_improving_stalled)
  test("improving_with_space_group",
    test_improving_with_space_group)
  test("recent_values_basic",
    test_recent_values_basic)
  test("recent_values_empty",
    test_recent_values_empty)
  test("metric_at_cycle",
    test_metric_at_cycle)
  print()

  # --- Gate evaluator usage pattern ---
  print("Gate evaluator usage pattern")
  test("phase_start_comparison",
    test_phase_start_comparison)
  test("threshold_after_n_cycles",
    test_threshold_after_n_cycles)
  print()

  # --- Explanation engine usage pattern ---
  print("Explanation engine usage pattern")
  test("consecutive_cycle_delta",
    test_consecutive_cycle_delta)
  print()

  # --- Serialization ---
  print("Serialization")
  test("roundtrip_empty",
    test_roundtrip_empty)
  test("roundtrip_populated",
    test_roundtrip_populated)
  test("roundtrip_through_json",
    test_roundtrip_through_json)
  test("from_dict_tolerant",
    test_from_dict_tolerant)
  test("from_dict_none",
    test_from_dict_none)
  test("json_serializable",
    test_json_serializable)
  test("sanitize_fallback",
    test_sanitize_fallback)
  print()

  # --- Multi-cycle scenario ---
  print("Multi-cycle scenario")
  test("five_cycle_scenario",
    test_five_cycle_scenario)
  test("enantiomorph_scenario",
    test_enantiomorph_scenario)
  print()

  # Summary
  print("=" * 60)
  total = passed + failed + skipped
  print(
    "Results: %d/%d passed, %d failed, %d skipped"
    % (passed, total, failed, skipped)
  )
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
         "rscc": 0.87, "rscc_z": 1.2, "b_mean": 22.5,
         "occupancy": 1.0},
      ],
    },
    "diff_peaks": {
      "positive": [
        {"height": 5.2, "near_residue": "A/His47"},
      ],
      "negative": [],
      "peak_count": 3,
    },
  }


def _make_log_metrics(**overrides):
  """Realistic log_metrics dict with overrides."""
  m = {
    "r_work": 0.210,
    "r_free": 0.248,
    "resolution": 2.1,
    "program": "phenix.refine",
  }
  m.update(overrides)
  return m


def _build_five_cycle_history():
  """Build a ValidationHistory with 5 cycles."""
  vh = ValidationHistory()
  # Cycle 1: xtriage (no R-factors, no validation)
  vh.record(1, "phenix.xtriage", None,
            {"resolution": 2.1},
            space_group="P 21 21 21")
  # Cycle 2: phaser
  vh.record(2, "phenix.phaser", None,
            {"tfz": 14.2, "llg": 890.0},
            space_group="P 21 21 21")
  # Cycle 3: refine round 1
  vr3 = {
    "geometry": {
      "clashscore": 12.0,
      "rama_favored": 0.92,
    },
  }
  vh.record(3, "phenix.refine", vr3,
            {"r_work": 0.30, "r_free": 0.38},
            space_group="P 21 21 21")
  # Cycle 4: refine round 2
  vr4 = {
    "geometry": {
      "clashscore": 6.5,
      "rama_favored": 0.95,
    },
    "model_contents": {"waters": 120},
  }
  vh.record(4, "phenix.refine", vr4,
            {"r_work": 0.25, "r_free": 0.31},
            space_group="P 21 21 21")
  # Cycle 5: refine round 3
  vh.record(5, "phenix.refine",
            _make_validation_result(),
            _make_log_metrics(),
            space_group="P 21 21 21")
  return vh


# ── Basics ────────────────────────────────────────────

def test_init_empty():
  vh = ValidationHistory()
  assert vh.snapshots == []
  assert vh.cycle_count() == 0


def test_safe_float_helper():
  assert _safe_float(3.14) == 3.14
  assert _safe_float("2.5") == 2.5
  assert _safe_float(None) is None
  assert _safe_float("bad") is None


def test_metric_direction():
  assert _metric_direction("r_free") == (
    "lower_is_better"
  )
  assert _metric_direction("clashscore") == (
    "lower_is_better"
  )
  assert _metric_direction("rama_favored") == (
    "higher_is_better"
  )
  assert _metric_direction("ligand_cc") == (
    "higher_is_better"
  )
  assert _metric_direction("unknown") == "neutral"


# ── Recording ─────────────────────────────────────────

def test_record_basic():
  vh = ValidationHistory()
  vh.record(
    1, "phenix.refine",
    _make_validation_result(),
    _make_log_metrics(),
    space_group="P 21 21 21",
  )
  assert vh.cycle_count() == 1
  s = vh.snapshots[0]
  assert s["cycle_number"] == 1
  assert s["program_name"] == "phenix.refine"
  assert s["space_group"] == "P 21 21 21"
  assert s["validation_result"] is not None
  assert s["log_metrics"]["r_free"] == 0.248
  assert "timestamp" in s


def test_record_none_validation():
  vh = ValidationHistory()
  vh.record(1, "phenix.xtriage", None,
            {"resolution": 2.1})
  assert vh.cycle_count() == 1
  s = vh.snapshots[0]
  assert s["validation_result"] is None
  assert s["log_metrics"]["resolution"] == 2.1


def test_record_idempotent():
  vh = ValidationHistory()
  vh.record(1, "phenix.refine", None,
            {"r_free": 0.30})
  vh.record(1, "phenix.refine", None,
            {"r_free": 0.28})
  # Should replace, not duplicate
  assert vh.cycle_count() == 1
  assert vh.snapshots[0]["log_metrics"]["r_free"] == (
    0.28
  )


def test_record_cap():
  vh = ValidationHistory()
  for i in range(110):
    vh.record(i, "test", None, {"r_free": 0.5})
  assert vh.cycle_count() <= 100


def test_record_never_raises():
  vh = ValidationHistory()
  # Bad types
  vh.record("bad", None, None, None)
  vh.record(1, None, "not_a_dict", 42)
  # Should not crash
  assert True


# ── Retrieval ─────────────────────────────────────────

def test_get_at_cycle():
  vh = _build_five_cycle_history()
  s = vh.get_at_cycle(3)
  assert s is not None
  assert s["cycle_number"] == 3
  assert s["program_name"] == "phenix.refine"


def test_get_at_cycle_missing():
  vh = _build_five_cycle_history()
  assert vh.get_at_cycle(99) is None


def test_get_latest():
  vh = _build_five_cycle_history()
  s = vh.get_latest()
  assert s is not None
  assert s["cycle_number"] == 5


def test_get_latest_empty():
  vh = ValidationHistory()
  assert vh.get_latest() is None


def test_get_previous():
  vh = _build_five_cycle_history()
  s = vh.get_previous()  # second-to-last
  assert s is not None
  assert s["cycle_number"] == 4


def test_get_previous_by_cycle():
  vh = _build_five_cycle_history()
  s = vh.get_previous(cycle_number=4)
  assert s is not None
  assert s["cycle_number"] == 3


def test_get_previous_insufficient():
  vh = ValidationHistory()
  vh.record(1, "test", None, {"r_free": 0.5})
  assert vh.get_previous() is None
  assert vh.get_previous(cycle_number=1) is None


# ── Metric extraction ────────────────────────────────

def _make_full_snapshot():
  """A snapshot with all metric types populated."""
  return {
    "cycle_number": 5,
    "program_name": "phenix.refine",
    "log_metrics": {
      "r_work": 0.210,
      "r_free": 0.248,
      "resolution": 2.1,
    },
    "validation_result": _make_validation_result(),
    "space_group": "P 21 21 21",
  }


def test_extract_r_factors():
  s = _make_full_snapshot()
  assert _extract_metric(s, "r_free") == 0.248
  assert _extract_metric(s, "r_work") == 0.210
  assert _extract_metric(s, "resolution") == 2.1


def test_extract_geometry():
  s = _make_full_snapshot()
  assert _extract_metric(s, "clashscore") == 4.2
  assert abs(
    _extract_metric(s, "rama_favored") - 0.968
  ) < 0.001
  assert _extract_metric(s, "bonds_rmsd") == 0.007


def test_extract_ligand_cc():
  s = _make_full_snapshot()
  assert _extract_metric(s, "ligand_cc") == 0.87


def test_extract_peak_count():
  s = _make_full_snapshot()
  assert _extract_metric(s, "peak_count") == 3


def test_extract_water_count():
  s = _make_full_snapshot()
  assert _extract_metric(s, "water_count") == 187


def test_extract_residue_count():
  s = _make_full_snapshot()
  assert _extract_metric(s, "residue_count") == 480


def test_extract_missing():
  s = {"log_metrics": {}, "validation_result": None}
  assert _extract_metric(s, "r_free") is None
  assert _extract_metric(s, "clashscore") is None
  assert _extract_metric(s, "nonexistent") is None


def test_extract_all_metrics():
  s = _make_full_snapshot()
  flat = _extract_all_metrics(s)
  assert flat["r_free"] == 0.248
  assert flat["r_work"] == 0.210
  assert flat["clashscore"] == 4.2
  assert flat["ligand_cc"] == 0.87
  assert flat["water_count"] == 187
  assert flat["residue_count"] == 480
  assert flat["peak_count"] == 3
  assert flat["resolution"] == 2.1
  # Should not have None values
  for v in flat.values():
    assert v is not None


# ── Metric series ─────────────────────────────────────

def test_series_rfree():
  vh = _build_five_cycle_history()
  series = vh.get_metric_series("r_free")
  # Cycles 3, 4, 5 have r_free
  assert len(series) == 3
  assert series[0] == (3, 0.38)
  assert series[1] == (4, 0.31)
  assert series[2] == (5, 0.248)
  # All decreasing
  values = [v for _, v in series]
  assert values == sorted(values, reverse=True)


def test_series_clashscore():
  vh = _build_five_cycle_history()
  series = vh.get_metric_series("clashscore")
  # Cycles 3, 4, 5
  assert len(series) == 3
  assert series[0] == (3, 12.0)
  assert series[1] == (4, 6.5)
  assert series[2] == (5, 4.2)


def test_series_missing_gaps():
  """Cycles without a metric are skipped."""
  vh = _build_five_cycle_history()
  series = vh.get_metric_series("tfz")
  # Only cycle 2 has tfz
  assert len(series) == 1
  assert series[0] == (2, 14.2)


def test_series_space_group_filter():
  """Space group filtering for enantiomorph tracking."""
  vh = ValidationHistory()
  # 3 cycles in P3121
  for i in range(1, 4):
    vh.record(i, "phenix.refine", None,
              {"r_free": 0.50 - i * 0.02},
              space_group="P 31 21")
  # 3 cycles in P3221
  for i in range(4, 7):
    vh.record(i, "phenix.refine", None,
              {"r_free": 0.45 - (i - 3) * 0.05},
              space_group="P 32 21")

  # All r_free values
  all_series = vh.get_metric_series("r_free")
  assert len(all_series) == 6

  # P3121 only
  p31 = vh.get_metric_series(
    "r_free", space_group="P 31 21"
  )
  assert len(p31) == 3
  assert p31[0] == (1, 0.48)
  assert p31[2] == (3, 0.44)

  # P3221 only
  p32 = vh.get_metric_series(
    "r_free", space_group="P 32 21"
  )
  assert len(p32) == 3
  assert p32[0] == (4, 0.40)
  assert p32[2] == (6, 0.30)


def test_series_never_raises():
  vh = ValidationHistory()
  # Empty history
  assert vh.get_metric_series("r_free") == []
  # Bad metric name
  assert vh.get_metric_series("nonexistent") == []


# ── Phase start metrics ───────────────────────────────

def test_phase_start_exact():
  vh = _build_five_cycle_history()
  s = vh.get_phase_start_metrics(3)
  assert s is not None
  assert s["cycle_number"] == 3


def test_phase_start_closest():
  """If exact cycle not recorded, find closest after."""
  vh = ValidationHistory()
  vh.record(2, "test", None, {"r_free": 0.40})
  vh.record(5, "test", None, {"r_free": 0.30})
  # Ask for cycle 3 — should get cycle 5
  s = vh.get_phase_start_metrics(3)
  assert s is not None
  assert s["cycle_number"] == 5


def test_phase_start_missing():
  vh = ValidationHistory()
  assert vh.get_phase_start_metrics(1) is None


# ── Flat metrics ──────────────────────────────────────

def test_flat_metrics_populated():
  vh = _build_five_cycle_history()
  flat = vh.get_flat_metrics(5)
  assert "r_free" in flat
  assert flat["r_free"] == 0.248
  assert "clashscore" in flat
  assert flat["clashscore"] == 4.2
  assert "water_count" in flat
  assert flat["water_count"] == 187


def test_flat_metrics_missing():
  vh = _build_five_cycle_history()
  flat = vh.get_flat_metrics(99)
  assert flat == {}


# ── Delta computation ─────────────────────────────────

def test_metrics_delta():
  vh = _build_five_cycle_history()
  delta = vh.get_metrics_delta(3, 5)
  assert "r_free" in delta
  d = delta["r_free"]
  assert d["before"] == 0.38
  assert d["after"] == 0.248
  assert abs(d["delta"] - (0.248 - 0.38)) < 0.001
  # Clashscore should also be present
  assert "clashscore" in delta
  cd = delta["clashscore"]
  assert cd["before"] == 12.0
  assert cd["after"] == 4.2


def test_metrics_delta_missing_cycle():
  vh = _build_five_cycle_history()
  delta = vh.get_metrics_delta(3, 99)
  assert delta == {}


def test_metrics_delta_partial_overlap():
  """Only metrics present in both cycles are included."""
  vh = ValidationHistory()
  vh.record(1, "xtriage", None,
            {"resolution": 2.1})
  vh.record(2, "refine", None,
            {"r_free": 0.30, "resolution": 2.1})
  delta = vh.get_metrics_delta(1, 2)
  # resolution is in both
  assert "resolution" in delta
  # r_free only in cycle 2 — NOT in delta
  assert "r_free" not in delta


# ── is_improving ──────────────────────────────────────

def test_improving_rfree():
  vh = _build_five_cycle_history()
  # r_free: 0.38 → 0.31 → 0.248 (decreasing = good)
  result = vh.is_improving("r_free", n_recent=3)
  assert result is True


def test_improving_rama_favored():
  vh = _build_five_cycle_history()
  # rama_favored: 0.92 → 0.95 → 0.968 (increasing = good)
  result = vh.is_improving(
    "rama_favored", n_recent=3
  )
  assert result is True


def test_improving_insufficient():
  vh = ValidationHistory()
  vh.record(1, "test", None, {"r_free": 0.30})
  result = vh.is_improving("r_free", n_recent=3)
  assert result is None  # insufficient data


def test_not_improving_stalled():
  vh = ValidationHistory()
  for i in range(1, 5):
    vh.record(i, "test", None,
              {"r_free": 0.35})  # constant
  result = vh.is_improving("r_free", n_recent=3)
  assert result is False  # not strictly less


def test_improving_with_space_group():
  """is_improving respects space group filter."""
  vh = ValidationHistory()
  # P3121: stalled
  for i in range(1, 4):
    vh.record(i, "refine", None,
              {"r_free": 0.45},
              space_group="P 31 21")
  # P3221: improving
  for i in range(4, 7):
    vh.record(i, "refine", None,
              {"r_free": 0.45 - (i - 3) * 0.05},
              space_group="P 32 21")
  # P3121 is NOT improving
  assert vh.is_improving(
    "r_free", n_recent=3,
    space_group="P 31 21"
  ) is False
  # P3221 IS improving
  assert vh.is_improving(
    "r_free", n_recent=3,
    space_group="P 32 21"
  ) is True


def test_recent_values_basic():
  vh = _build_five_cycle_history()
  # Last 3 R-free values
  vals = vh.get_recent_values("r_free", n=3)
  assert len(vals) == 3
  assert vals == [0.38, 0.31, 0.248]


def test_recent_values_empty():
  vh = ValidationHistory()
  assert vh.get_recent_values("r_free") == []


def test_metric_at_cycle():
  vh = _build_five_cycle_history()
  assert vh.get_metric_at_cycle("r_free", 5) == 0.248
  assert vh.get_metric_at_cycle("r_free", 1) is None
  assert vh.get_metric_at_cycle("r_free", 99) is None
  assert vh.get_metric_at_cycle(
    "clashscore", 3
  ) == 12.0


# ── Gate evaluator usage pattern ──────────────────────

def test_phase_start_comparison():
  """Gate evaluator pattern: compare stage-start metrics
  against current metrics for monotonic progress gate."""
  vh = _build_five_cycle_history()
  # Phase started at cycle 3 (first refine)
  start = vh.get_phase_start_metrics(3)
  assert start is not None
  start_rf = _extract_metric(start, "r_free")
  assert start_rf == 0.38
  start_cs = _extract_metric(start, "clashscore")
  assert start_cs == 12.0
  # Current is cycle 5
  current_rf = vh.get_metric_at_cycle("r_free", 5)
  assert current_rf == 0.248
  current_cs = vh.get_metric_at_cycle(
    "clashscore", 5
  )
  assert current_cs == 4.2
  # Metrics improved (lower is better for both)
  assert current_rf < start_rf
  assert current_cs < start_cs


def test_threshold_after_n_cycles():
  """Gate evaluator pattern: check if r_free > threshold
  after N cycles in current step."""
  vh = ValidationHistory()
  # Phase starts at cycle 5 with R-free stuck
  for i in range(5, 8):
    vh.record(i, "phenix.refine", None,
              {"r_free": 0.46})
  # Get last 2 values
  recent = vh.get_recent_values("r_free", n=2)
  assert len(recent) == 2
  # All above 0.45 threshold → retreat trigger
  assert all(v > 0.45 for v in recent)


# ── Explanation engine usage pattern ──────────────────

def test_consecutive_cycle_delta():
  """Explanation engine pattern: compute metrics_before
  and metrics_after for per-cycle commentary."""
  vh = _build_five_cycle_history()
  # Commentary for cycle 5 needs cycle 4 vs cycle 5
  delta = vh.get_metrics_delta(4, 5)
  assert "r_free" in delta
  rf = delta["r_free"]
  assert rf["before"] == 0.31
  assert rf["after"] == 0.248
  assert abs(rf["delta"] - (0.248 - 0.31)) < 0.001
  # Clashscore also available
  assert "clashscore" in delta
  cs = delta["clashscore"]
  assert cs["before"] == 6.5
  assert cs["after"] == 4.2


# ── Serialization ─────────────────────────────────────

def test_roundtrip_empty():
  vh = ValidationHistory()
  d = vh.to_dict()
  vh2 = ValidationHistory.from_dict(d)
  assert vh2.cycle_count() == 0


def test_roundtrip_populated():
  vh = _build_five_cycle_history()
  d = vh.to_dict()
  vh2 = ValidationHistory.from_dict(d)
  assert vh2.cycle_count() == 5
  # Check a specific snapshot
  s = vh2.get_at_cycle(5)
  assert s is not None
  assert s["program_name"] == "phenix.refine"
  assert s["log_metrics"]["r_free"] == 0.248
  # Metric series should work on restored history
  series = vh2.get_metric_series("r_free")
  assert len(series) == 3
  assert series[-1] == (5, 0.248)


def test_roundtrip_through_json():
  """Full round-trip: to_dict → JSON → from_dict."""
  vh = _build_five_cycle_history()
  d = vh.to_dict()
  json_str = json.dumps(d)
  d2 = json.loads(json_str)
  vh2 = ValidationHistory.from_dict(d2)
  assert vh2.cycle_count() == 5
  # Verify data survived JSON
  flat = vh2.get_flat_metrics(5)
  assert flat["r_free"] == 0.248
  assert flat["clashscore"] == 4.2


def test_from_dict_tolerant():
  # Missing snapshots key
  vh = ValidationHistory.from_dict({"extra": True})
  assert vh.cycle_count() == 0
  # Snapshots with non-dict entries (should skip)
  vh = ValidationHistory.from_dict({
    "snapshots": [
      {"cycle_number": 1, "program_name": "test",
       "log_metrics": {}, "validation_result": None,
       "space_group": None, "timestamp": 0},
      "garbage",
      42,
    ],
  })
  assert vh.cycle_count() == 1


def test_from_dict_none():
  vh = ValidationHistory.from_dict(None)
  assert vh.cycle_count() == 0


def test_json_serializable():
  """to_dict output must be fully JSON-serializable."""
  vh = _build_five_cycle_history()
  d = vh.to_dict()
  json_str = json.dumps(d)
  assert len(json_str) > 0


def test_sanitize_fallback():
  """_sanitize_snapshot fallback when deepcopy fails.

  In production, validation_result may contain cctbx
  objects that aren't deepcopy-able. The fallback
  should preserve scalar fields and drop the rest.
  """
  from agent.validation_history import (
    _sanitize_snapshot,
  )

  class NotCopyable(object):
    def __deepcopy__(self, memo):
      raise TypeError("cannot deepcopy")
    def __copy__(self):
      raise TypeError("cannot copy")

  snapshot = {
    "cycle_number": 5,
    "program_name": "phenix.refine",
    "validation_result": {
      "geometry": NotCopyable(),
    },
    "log_metrics": {"r_free": 0.25},
    "space_group": "P 21 21 21",
    "timestamp": 1234567890.0,
  }
  result = _sanitize_snapshot(snapshot)
  # Fallback should keep scalars
  assert result["cycle_number"] == 5
  assert result["program_name"] == "phenix.refine"
  assert result["space_group"] == "P 21 21 21"
  # Should be JSON-safe
  json.dumps(result)


# ── Multi-cycle scenarios ─────────────────────────────

def test_five_cycle_scenario():
  """Spec test: Run 5 cycles, verify get_metric_series
  returns correct R-free series."""
  vh = _build_five_cycle_history()

  # R-free series
  rf_series = vh.get_metric_series("r_free")
  assert len(rf_series) == 3
  cycles = [c for c, _ in rf_series]
  values = [v for _, v in rf_series]
  assert cycles == [3, 4, 5]
  assert abs(values[0] - 0.38) < 0.001
  assert abs(values[1] - 0.31) < 0.001
  assert abs(values[2] - 0.248) < 0.001

  # Clashscore series
  cs_series = vh.get_metric_series("clashscore")
  assert len(cs_series) == 3
  cs_values = [v for _, v in cs_series]
  assert cs_values == [12.0, 6.5, 4.2]

  # Phase start at cycle 3
  ps = vh.get_phase_start_metrics(3)
  assert ps is not None
  assert ps["log_metrics"]["r_free"] == 0.38

  # Delta from stage start to current
  delta = vh.get_metrics_delta(3, 5)
  assert abs(
    delta["r_free"]["delta"] - (0.248 - 0.38)
  ) < 0.001

  # is_improving
  assert vh.is_improving("r_free") is True
  assert vh.is_improving("clashscore") is True


def test_enantiomorph_scenario():
  """Spec test: track R-free for both enantiomorphs,
  confirm trajectories don't conflate."""
  vh = ValidationHistory()

  # P3121: 3 cycles, R-free stuck at 0.48
  for i in range(1, 4):
    vh.record(i, "phenix.refine", None,
              {"r_free": 0.48},
              space_group="P 31 21")

  # P3221: 3 cycles, R-free improving
  rfrees = [0.40, 0.35, 0.30]
  for i, rf in enumerate(rfrees, start=4):
    vh.record(i, "phenix.refine", None,
              {"r_free": rf},
              space_group="P 32 21")

  # All 6 cycles recorded
  all_rf = vh.get_metric_series("r_free")
  assert len(all_rf) == 6

  # P3121: stalled at 0.48
  p31 = vh.get_metric_series(
    "r_free", space_group="P 31 21"
  )
  assert len(p31) == 3
  assert all(v == 0.48 for _, v in p31)

  # P3221: improving
  p32 = vh.get_metric_series(
    "r_free", space_group="P 32 21"
  )
  assert len(p32) == 3
  assert p32[0] == (4, 0.40)
  assert p32[2] == (6, 0.30)

  # P3121 is NOT improving
  assert vh.is_improving(
    "r_free", n_recent=3,
    space_group="P 31 21"
  ) is False

  # P3221 IS improving
  assert vh.is_improving(
    "r_free", n_recent=3,
    space_group="P 32 21"
  ) is True

  # Survives serialization
  d = vh.to_dict()
  vh2 = ValidationHistory.from_dict(
    json.loads(json.dumps(d))
  )
  p32_restored = vh2.get_metric_series(
    "r_free", space_group="P 32 21"
  )
  assert len(p32_restored) == 3
  assert p32_restored[2] == (6, 0.30)


# ── Entry point ──────────────────────────────────────

if __name__ == "__main__":
  run_tests()
