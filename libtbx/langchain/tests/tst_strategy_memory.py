"""Tests for agent/strategy_memory.py (v113 Thinking Agent)."""

from __future__ import absolute_import, division, print_function
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agent.strategy_memory import StrategyMemory


def test_empty_construction():
  print("Test: empty_construction")
  m = StrategyMemory()
  assert m.data_quality == ""
  assert m.phasing_strategy == ""
  assert m.concerns == []
  assert m.decisions == []
  assert m.r_free_history == []
  assert m.programs_run == []
  assert m.summary() == "(empty)"
  print("  PASSED")


def test_to_dict_from_dict_roundtrip():
  print("Test: to_dict_from_dict_roundtrip")
  m = StrategyMemory()
  m.data_quality = "twinned"
  m.phasing_strategy = "SAD"
  m.concerns = ["weak anomalous"]
  m.r_free_history = [0.30, 0.28, 0.27]
  m.programs_run = ["phenix.xtriage", "phenix.autosol"]
  m.decisions = [[1, "Try SAD first"]]

  d = m.to_dict()
  j = json.dumps(d)
  m2 = StrategyMemory.from_dict(json.loads(j))
  assert m2.to_dict() == d
  print("  PASSED")


def test_from_dict_tolerates_bad_input():
  print("Test: from_dict_tolerates_bad_input")
  m1 = StrategyMemory.from_dict(None)
  assert m1.data_quality == ""
  m2 = StrategyMemory.from_dict("garbage")
  assert m2.phasing_strategy == ""
  m3 = StrategyMemory.from_dict({"extra_key": "ignored"})
  assert m3.concerns == []
  m4 = StrategyMemory.from_dict({"data_quality": "good"})
  assert m4.data_quality == "good"
  print("  PASSED")


def test_update():
  print("Test: update")
  m = StrategyMemory()
  m.update({
    "data_quality": "anisotropic",
    "phasing_strategy": "MR",
    "concerns": ["ice rings", "weak diffraction"],
  })
  assert m.data_quality == "anisotropic"
  assert m.phasing_strategy == "MR"
  assert len(m.concerns) == 2
  assert "ice rings" in m.concerns

  # Second update should not duplicate concerns
  m.update({"concerns": ["ice rings", "twinning"]})
  assert m.concerns.count("ice rings") == 1
  assert "twinning" in m.concerns
  print("  PASSED")


def test_update_bad_input():
  print("Test: update_bad_input")
  m = StrategyMemory()
  m.update(None)
  m.update("garbage")
  m.update({"concerns": "not a list"})
  assert m.data_quality == ""
  print("  PASSED")


def test_concerns_cap():
  print("Test: concerns_cap")
  m = StrategyMemory()
  m.update({"concerns": ["c%d" % i for i in range(15)]})
  assert len(m.concerns) <= 10
  print("  PASSED")


def test_record_outcome():
  print("Test: record_outcome")
  m = StrategyMemory()
  m.record_outcome("phenix.refine", "SUCCESS", {"r_free": 0.30})
  m.record_outcome("phenix.refine", "SUCCESS", {"r_free": 0.28})
  assert m.programs_run == ["phenix.refine", "phenix.refine"]
  assert m.r_free_history == [0.30, 0.28]

  # Without metrics
  m.record_outcome("phenix.xtriage", "SUCCESS")
  assert len(m.programs_run) == 3
  assert len(m.r_free_history) == 2

  # Bad r_free value
  m.record_outcome("phenix.refine", "SUCCESS", {"r_free": "bad"})
  assert len(m.r_free_history) == 2
  print("  PASSED")


def test_programs_run_cap():
  print("Test: programs_run_cap")
  m = StrategyMemory()
  for i in range(40):
    m.record_outcome("prog_%d" % i, "SUCCESS")
  assert len(m.programs_run) <= 30
  print("  PASSED")


def test_record_decision():
  print("Test: record_decision")
  m = StrategyMemory()
  m.record_decision(1, "Start with MR")
  m.record_decision(3, "Switch to SAD")
  assert len(m.decisions) == 2
  assert m.decisions[0] == [1, "Start with MR"]
  print("  PASSED")


def test_decisions_cap():
  print("Test: decisions_cap")
  m = StrategyMemory()
  for i in range(15):
    m.record_decision(i, "decision %d" % i)
  assert len(m.decisions) <= 10
  print("  PASSED")


def test_metrics_stalled_not_enough_data():
  print("Test: metrics_stalled_not_enough_data")
  m = StrategyMemory()
  assert m.metrics_stalled() == False
  m.r_free_history = [0.30]
  assert m.metrics_stalled() == False
  m.r_free_history = [0.30, 0.29]
  assert m.metrics_stalled() == False
  print("  PASSED")


def test_metrics_stalled_improving():
  print("Test: metrics_stalled_improving")
  m = StrategyMemory()
  m.r_free_history = [0.30, 0.28, 0.26]
  assert m.metrics_stalled() == False
  print("  PASSED")


def test_metrics_stalled_plateau():
  print("Test: metrics_stalled_plateau")
  m = StrategyMemory()
  m.r_free_history = [0.28, 0.280, 0.281]
  assert m.metrics_stalled() == True
  # Getting worse also counts as stalled
  m.r_free_history = [0.28, 0.29, 0.30]
  assert m.metrics_stalled() == True
  print("  PASSED")


def test_summary():
  print("Test: summary")
  m = StrategyMemory()
  assert m.summary() == "(empty)"
  m.data_quality = "good"
  m.r_free_history = [0.25]
  assert "quality=good" in m.summary()
  assert "r_free=0.250" in m.summary()
  print("  PASSED")


# =========================================================================

def run_tests_with_fail_fast():
  tests = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
  passed = 0
  failed = 0
  for t in tests:
    name = t.__name__
    try:
      t()
      passed += 1
      print("  PASS: %s" % name)
    except Exception as e:
      failed += 1
      print("  FAIL: %s — %s" % (name, e))
      break
  print("\n%d passed, %d failed" % (passed, failed))
  return failed == 0


def run_all_tests():
  run_tests_with_fail_fast()


if __name__ == "__main__":
  success = run_tests_with_fail_fast()
  sys.exit(0 if success else 1)
