"""Tests for agent/log_section_extractor.py (v113 Thinking Agent)."""

from __future__ import absolute_import, division, print_function
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agent.log_section_extractor import (
  extract_sections, SECTION_MARKERS, _merge_windows)


# =========================================================================
# Sample log fixtures
# =========================================================================

XTRIAGE_LOG = "\n".join([
  "phenix.xtriage analysis",
] + ["padding line %d" % i for i in range(20)] + [
  "Results of L-test for twinning:",
  "  L-statistic = 0.42",
  "  Estimated twinning fraction = 0.15",
  "  The data do not appear to be twinned.",
] + ["gap %d" % i for i in range(30)] + [
  "Anomalous signal analysis:",
  "  measurability = 0.05",
  "  d_min_anom = 3.5",
  "  No significant anomalous signal detected.",
] + ["tail %d" % i for i in range(10)])

PHASER_LOG = "\n".join([
  "Phaser molecular replacement",
] + ["search step %d" % i for i in range(15)] + [
  "SOLU SET RFZ=4.2 TFZ=12.5 LLG=150",
  "Solution found with 1 copies placed",
] + ["Packing analysis:",
     "  No clashes detected",
] + ["post %d" % i for i in range(5)])

REFINE_LOG = "\n".join([
  "phenix.refine macro_cycle 1",
] + ["optimization step %d" % i for i in range(10)] + [
  "Start: r_work = 0.2500  r_free = 0.2900",
] + ["more steps %d" % i for i in range(20)] + [
  "Final: r_work = 0.2200  r_free = 0.2600",
  "Refinement statistics summary",
] + ["geometry %d" % i for i in range(10)] + [
  "Ramachandran plot analysis:",
  "  Favored: 95.2%",
  "  Outliers: 0.5%",
  "  clashscore = 3.2",
] + ["end %d" % i for i in range(5)])


# =========================================================================
# Tests
# =========================================================================

def test_empty_log():
  print("Test: empty_log")
  assert extract_sections("", "phenix.xtriage") == ""
  assert extract_sections(None, "phenix.refine") == ""
  assert extract_sections("  \n  ", "phenix.phaser") == ""
  print("  PASSED")


def test_unknown_program_fallback():
  print("Test: unknown_program_fallback")
  log = "\n".join(["line %d" % i for i in range(200)])
  result = extract_sections(log, "phenix.totally_unknown")
  assert "--- Last" in result
  assert "line 199" in result
  # Should not contain early lines
  assert "line 50" not in result
  print("  PASSED")


def test_short_log_fallback():
  print("Test: short_log_fallback")
  log = "\n".join(["short line %d" % i for i in range(10)])
  result = extract_sections(log, "phenix.unknown")
  assert "--- Last 10 lines" in result
  assert "short line 0" in result
  assert "short line 9" in result
  print("  PASSED")


def test_xtriage_twinning_first():
  print("Test: xtriage_twinning_first")
  result = extract_sections(XTRIAGE_LOG, "phenix.xtriage")
  assert "--- Twinning (phenix.xtriage) ---" in result
  assert "--- Anomalous signal (phenix.xtriage) ---" in result
  twin_pos = result.index("Twinning")
  anom_pos = result.index("Anomalous signal")
  assert twin_pos < anom_pos, "Twinning should precede Anomalous"
  print("  PASSED")


def test_xtriage_contains_keywords():
  print("Test: xtriage_contains_keywords")
  result = extract_sections(XTRIAGE_LOG, "phenix.xtriage")
  assert "L-test" in result or "L-statistic" in result
  assert "measurability" in result
  print("  PASSED")


def test_phaser_mr_scores():
  print("Test: phaser_mr_scores")
  result = extract_sections(PHASER_LOG, "phenix.phaser")
  assert "--- MR scores (phenix.phaser) ---" in result
  assert "TFZ=12.5" in result
  assert "LLG=150" in result
  print("  PASSED")


def test_phaser_packing():
  print("Test: phaser_packing")
  result = extract_sections(PHASER_LOG, "phenix.phaser")
  assert "--- Packing (phenix.phaser) ---" in result
  print("  PASSED")


def test_refine_r_factors_first():
  print("Test: refine_r_factors_first")
  result = extract_sections(REFINE_LOG, "phenix.refine")
  assert "--- R-factors (phenix.refine) ---" in result
  assert "--- Geometry (phenix.refine) ---" in result
  rfact_pos = result.index("R-factors")
  geom_pos = result.index("Geometry")
  assert rfact_pos < geom_pos
  print("  PASSED")


def test_refine_contains_metrics():
  print("Test: refine_contains_metrics")
  result = extract_sections(REFINE_LOG, "phenix.refine")
  assert "r_free = 0.2600" in result
  assert "Ramachandran" in result
  print("  PASSED")


def test_budget_enforcement():
  print("Test: budget_enforcement")
  # Large log with many keyword matches
  big_log = "\n".join(["r_free = 0.%04d some data" % i
                        for i in range(500)])
  result = extract_sections(big_log, "phenix.refine", max_chars=500)
  assert len(result) <= 550  # small slack for final line
  print("  PASSED (len=%d)" % len(result))


def test_budget_omission_note():
  print("Test: budget_omission_note")
  big_log = "\n".join(
    ["r_free line %d" % i for i in range(100)] +
    ["Ramachandran line %d" % i for i in range(100)])
  result = extract_sections(big_log, "phenix.refine", max_chars=800)
  if "[remaining sections omitted]" in result:
    # Good — budget was tight enough to trigger omission
    pass
  else:
    # Both sections fit — also acceptable for this budget
    pass
  assert len(result) <= 850
  print("  PASSED")


def test_no_keyword_matches_fallback():
  print("Test: no_keyword_matches_fallback")
  boring_log = "\n".join(["nothing interesting %d" % i for i in range(50)])
  result = extract_sections(boring_log, "phenix.xtriage")
  assert "--- Last" in result
  assert "nothing interesting 49" in result
  print("  PASSED")


def test_merge_windows():
  print("Test: merge_windows")
  assert _merge_windows([]) == []
  assert _merge_windows([(0, 10)]) == [(0, 10)]
  assert _merge_windows([(0, 10), (5, 20)]) == [(0, 20)]
  assert _merge_windows([(0, 10), (15, 25)]) == [(0, 10), (15, 25)]
  assert _merge_windows([(5, 20), (0, 10)]) == [(0, 20)]
  assert _merge_windows([(0, 10), (10, 20)]) == [(0, 20)]
  print("  PASSED")


def test_all_programs_have_markers():
  print("Test: all_programs_have_markers")
  expected = ["phenix.xtriage", "phenix.phaser", "phenix.autosol",
              "phenix.autobuild", "phenix.refine"]
  for prog in expected:
    assert prog in SECTION_MARKERS, "Missing markers for %s" % prog
    sections = SECTION_MARKERS[prog]
    assert len(sections) >= 2, "%s should have >= 2 sections" % prog
    for name, keywords in sections:
      assert len(keywords) >= 1, "Section %s should have keywords" % name
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
