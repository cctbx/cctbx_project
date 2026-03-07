"""
Unit tests for explanation_prompts.py (Phase 5).

Run standalone:
  python tests/tst_explanation_prompts.py

No PHENIX dependencies — pure Python.

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import traceback

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_THIS_DIR)
if _PROJECT_ROOT not in sys.path:
  sys.path.insert(0, _PROJECT_ROOT)

from knowledge.explanation_prompts import (
  generate_cycle_commentary,
  generate_stage_summary,
  generate_final_report,
  generate_stopped_report,
)
from agent.structure_model import (
  StructureModel, Hypothesis,
)
from knowledge.plan_schema import (
  StageDef, StructurePlan,
)


def run_tests():
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
  print("Explanation Prompts Unit Tests")
  print("=" * 60)
  print()

  # --- Per-cycle commentary ---
  print("Per-cycle commentary")
  test("commentary_rfree_improved",
    test_commentary_rfree_improved)
  test("commentary_rfree_worsened",
    test_commentary_rfree_worsened)
  test("commentary_with_model_contents",
    test_commentary_with_model_contents)
  test("commentary_with_problems",
    test_commentary_with_problems)
  test("commentary_cryoem",
    test_commentary_cryoem)
  test("commentary_no_data",
    test_commentary_no_data)
  test("commentary_none_inputs",
    test_commentary_none_inputs)
  print()

  # --- Stage summary ---
  print("Stage summary")
  test("phase_complete",
    test_phase_complete)
  test("phase_skipped",
    test_stage_skipped)
  test("phase_failed",
    test_stage_failed)
  test("phase_with_next",
    test_stage_with_next)
  test("phase_none_input",
    test_stage_none_input)
  print()

  # --- Final report ---
  print("Final report")
  test("final_report_full",
    test_final_report_full)
  test("final_report_minimal",
    test_final_report_minimal)
  test("final_report_with_hypotheses",
    test_final_report_with_hypotheses)
  test("final_report_none_inputs",
    test_final_report_none_inputs)
  print()

  # --- Stopped report ---
  print("Stopped report")
  test("stopped_report_with_blacklist",
    test_stopped_report_with_blacklist)
  test("stopped_report_recommendations",
    test_stopped_report_recommendations)
  test("stopped_report_none_inputs",
    test_stopped_report_none_inputs)
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

def _make_sm():
  """Build a populated StructureModel."""
  sm = StructureModel()
  sm.data_characteristics["resolution"] = 2.1
  sm.data_characteristics["space_group"] = (
    "P 21 21 21"
  )
  sm.data_characteristics["experiment_type"] = "xray"
  sm.model_state["r_free"] = 0.28
  sm.model_state["r_work"] = 0.22
  sm.model_state["waters"] = 187
  sm.model_state["chains"] = [
    {"chain_id": "A", "completeness": 0.96},
    {"chain_id": "B", "completeness": 0.88},
  ]
  sm.model_state["ligands"] = [
    {"name": "ATP", "chain": "A",
     "resid": 301, "rscc": 0.87},
  ]
  sm.model_state["ions"] = [
    {"name": "MG", "chain": "A", "resid": 500},
  ]
  sm.model_state["geometry"] = {
    "clashscore": 4.2,
    "rama_favored": 0.968,
    "rama_outliers": 0.005,
    "rotamer_outliers": 0.02,
    "bonds_rmsd": 0.007,
    "angles_rmsd": 1.1,
  }
  return sm


def _make_plan():
  """Build a completed plan."""
  stages = [
    StageDef(
      id="data_assessment",
      programs=["phenix.xtriage"],
      max_cycles=1,
      description="Assess data quality",
    ),
    StageDef(
      id="molecular_replacement",
      programs=["phenix.phaser"],
      max_cycles=1,
      success_criteria={"tfz": ">8"},
      description="Find MR solution",
    ),
    StageDef(
      id="initial_refinement",
      programs=["phenix.refine"],
      max_cycles=3,
      success_criteria={"r_free": "<0.35"},
      description="Initial refinement",
    ),
    StageDef(
      id="final_refinement",
      programs=["phenix.refine"],
      max_cycles=3,
      success_criteria={"r_free": "<0.25"},
      description="Final refinement",
    ),
  ]
  # Mark all complete with metrics
  stages[0].status = "complete"
  stages[0].cycles_used = 1
  stages[1].status = "complete"
  stages[1].cycles_used = 1
  stages[1].result_metrics = {"tfz": 14.2}
  stages[2].status = "complete"
  stages[2].cycles_used = 3
  stages[2].result_metrics = {"r_free": 0.32}
  stages[3].status = "complete"
  stages[3].cycles_used = 3
  stages[3].result_metrics = {"r_free": 0.24}

  return StructurePlan(
    goal="Solve kinase at 2.1A",
    stages=stages,
    template_id="mr_refine",
  )


# ── Per-cycle commentary ──────────────────────────────

def test_commentary_rfree_improved():
  sm = _make_sm()
  before = {"r_free": 0.30, "r_work": 0.24}
  after = {"r_free": 0.28, "r_work": 0.22}
  text = generate_cycle_commentary(
    sm, 5, "phenix.refine", before, after
  )
  assert "Cycle 5" in text
  assert "phenix.refine" in text
  assert "improved" in text
  assert "0.300" in text  # before
  assert "0.280" in text  # after


def test_commentary_rfree_worsened():
  sm = _make_sm()
  before = {"r_free": 0.28}
  after = {"r_free": 0.32}
  text = generate_cycle_commentary(
    sm, 6, "phenix.refine", before, after
  )
  assert "worsened" in text


def test_commentary_with_model_contents():
  sm = _make_sm()
  after = {"r_free": 0.28}
  text = generate_cycle_commentary(
    sm, 5, "phenix.refine", {}, after
  )
  assert "187 waters" in text
  assert "ATP" in text
  assert "RSCC" in text


def test_commentary_with_problems():
  sm = _make_sm()
  # Add a problem
  sm.model_state["geometry"]["clashscore"] = 25.0
  sm._detect_problems({})
  after = {"r_free": 0.28}
  text = generate_cycle_commentary(
    sm, 5, "phenix.refine", {}, after
  )
  assert "Issue" in text or "clashscore" in text.lower()


def test_commentary_cryoem():
  sm = StructureModel()
  sm.model_state["model_map_cc"] = 0.72
  before = {"model_map_cc": 0.65}
  after = {"model_map_cc": 0.72}
  text = generate_cycle_commentary(
    sm, 3, "phenix.real_space_refine",
    before, after,
  )
  assert "model-map CC" in text
  assert "improved" in text


def test_commentary_no_data():
  text = generate_cycle_commentary(
    None, 1, "phenix.xtriage", {}, {}
  )
  assert text == ""


def test_commentary_none_inputs():
  text = generate_cycle_commentary(
    None, 1, None, None, None
  )
  assert isinstance(text, str)


# ── Stage summary ─────────────────────────────────────

def test_phase_complete():
  sm = _make_sm()
  stage = StageDef(
    id="initial_refinement",
    description="Initial refinement",
    max_cycles=3,
  )
  stage.status = "complete"
  stage.cycles_used = 3
  stage.result_metrics = {"r_free": 0.32}

  text = generate_stage_summary(
    sm, None, stage, None
  )
  assert "completed" in text
  assert "3 cycle" in text
  assert "0.320" in text


def test_stage_skipped():
  stage = StageDef(id="model_rebuilding")
  stage.status = "skipped"
  text = generate_stage_summary(
    None, None, stage, None
  )
  assert "skipped" in text


def test_stage_failed():
  stage = StageDef(id="initial_refinement")
  stage.status = "failed"
  stage.cycles_used = 2
  text = generate_stage_summary(
    None, None, stage, None
  )
  assert "failed" in text
  assert "2 cycle" in text


def test_stage_with_next():
  sm = _make_sm()
  completed = StageDef(
    id="initial_refinement",
    description="Initial refinement",
  )
  completed.status = "complete"
  completed.cycles_used = 3
  next_p = StageDef(
    id="model_rebuilding",
    description="Rebuild model",
  )
  text = generate_stage_summary(
    sm, None, completed, next_p
  )
  # Summary should mention the completed stage
  assert "initial_refinement" in text
  assert "3 cycle" in text
  # "Next:" line was removed (shown in GUI gate
  # transition block instead), so next_stage info
  # should NOT appear in the summary text.
  assert "model_rebuilding" not in text


def test_stage_none_input():
  text = generate_stage_summary(
    None, None, None, None
  )
  assert text == ""


# ── Final report ──────────────────────────────────────

def test_final_report_full():
  sm = _make_sm()
  plan = _make_plan()
  report = generate_final_report(sm, plan)
  assert "STRUCTURE DETERMINATION REPORT" in report
  assert "2.10" in report  # resolution
  assert "P 21 21 21" in report  # space group
  assert "R-free: 0.2800" in report
  assert "187" in report  # waters
  assert "ATP" in report  # ligand
  assert "STAGE TIMELINE" in report
  assert "[done]" in report
  assert "data_assessment" in report
  # Check stage metrics
  assert "14.200" in report  # TFZ


def test_final_report_minimal():
  sm = StructureModel()
  report = generate_final_report(sm, None)
  assert "STRUCTURE DETERMINATION REPORT" in report


def test_final_report_with_hypotheses():
  sm = _make_sm()
  h1 = Hypothesis(
    id="h1",
    statement="Density near His47 is Zn ion",
    status="confirmed",
    resolved_at_cycle=8,
  )
  h2 = Hypothesis(
    id="h2",
    statement="ATP is in wrong orientation",
    status="refuted",
    resolved_at_cycle=10,
  )
  sm.hypotheses = [h1, h2]
  report = generate_final_report(sm, None)
  assert "HYPOTHESES" in report
  assert "Zn ion" in report
  assert "confirmed" in report
  assert "refuted" in report


def test_final_report_none_inputs():
  report = generate_final_report(None, None)
  assert isinstance(report, str)
  assert len(report) > 0


# ── Stopped report ────────────────────────────────────

def test_stopped_report_with_blacklist():
  sm = _make_sm()
  sm.model_state["r_free"] = 0.48
  sm.blacklist_strategy(
    "MR_with_beta",
    "R-free stuck at 0.48 after 3 cycles",
    {"r_free": 0.48},
  )
  plan = _make_plan()
  plan.stages[2].status = "failed"
  plan.stages[3].status = "pending"

  report = generate_stopped_report(sm, plan)
  assert "STOPPED" in report
  assert "R-free: 0.4800" in report
  assert "STRATEGIES THAT FAILED" in report
  assert "stuck at 0.48" in report
  assert "RECOMMENDATIONS" in report


def test_stopped_report_recommendations():
  sm = _make_sm()
  sm.model_state["r_free"] = 0.46
  sm.model_state["geometry"]["clashscore"] = 20.0
  sm.model_state["ligands"][0]["rscc"] = 0.45

  report = generate_stopped_report(sm, None)
  assert "RECOMMENDATIONS" in report
  # Should recommend checking search model
  assert "search model" in report.lower() or \
    "data quality" in report.lower()
  # Should recommend checking clashscore
  assert "clashscore" in report.lower() or \
    "side chain" in report.lower()
  # Should recommend checking ligand
  assert "ligand" in report.lower()


def test_stopped_report_none_inputs():
  report = generate_stopped_report(None, None)
  assert "STOPPED" in report
  assert "No structure model" in report


# ── Entry point ──────────────────────────────────────

if __name__ == "__main__":
  run_tests()
