"""
Unit tests for hypothesis_evaluator.py (Phase 4).

Run standalone:
  python tests/tst_hypothesis_evaluator.py

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

from agent.hypothesis_evaluator import (
  evaluate_hypotheses,
  revalidate_confirmed,
  build_hypothesis_prompt,
  _check_criteria,
)
from agent.structure_model import (
  StructureModel, Hypothesis,
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
  print("Hypothesis Evaluator Unit Tests")
  print("=" * 60)
  print()

  # --- Criteria checking ---
  print("Criteria checking")
  test("single_criterion_met",
    test_single_criterion_met)
  test("single_criterion_not_met",
    test_single_criterion_not_met)
  test("compound_and_all_met",
    test_compound_and_all_met)
  test("compound_and_partial",
    test_compound_and_partial)
  test("criterion_missing_metric",
    test_criterion_missing_metric)
  test("empty_criteria",
    test_empty_criteria)
  print()

  # --- Lifecycle ---
  print("Lifecycle")
  test("testing_to_pending",
    test_testing_to_pending)
  test("pending_countdown",
    test_pending_countdown)
  test("confirmed",
    test_confirmed)
  test("refuted",
    test_refuted)
  test("abandoned_timeout",
    test_abandoned_timeout)
  test("waiting_inconclusive",
    test_waiting_inconclusive)
  test("no_active_hypotheses",
    test_no_active_hypotheses)
  test("none_structure_model",
    test_none_structure_model)
  print()

  # --- Revalidation ---
  print("Revalidation")
  test("revalidate_still_valid",
    test_revalidate_still_valid)
  test("revalidate_ligand_degraded",
    test_revalidate_ligand_degraded)
  test("revalidate_rfree_degraded",
    test_revalidate_rfree_degraded)
  test("revalidate_no_confirmed",
    test_revalidate_no_confirmed)
  print()

  # --- Prompt building ---
  print("Prompt building")
  test("prompt_active_hypothesis",
    test_prompt_active_hypothesis)
  test("prompt_revalidation_reason",
    test_prompt_revalidation_reason)
  test("prompt_opportunity",
    test_prompt_opportunity)
  test("prompt_nothing_to_hypothesize",
    test_prompt_nothing_to_hypothesize)
  test("prompt_none_sm",
    test_prompt_none_sm)
  print()

  # --- Full scenario ---
  print("Full scenario")
  test("ligand_hypothesis_lifecycle",
    test_ligand_hypothesis_lifecycle)
  test("metal_ion_with_latency",
    test_metal_ion_with_latency)
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

def _make_sm(**kwargs):
  """Build a StructureModel with given metrics."""
  sm = StructureModel()
  for k, v in kwargs.items():
    if k in ("r_free", "r_work", "model_map_cc"):
      sm.model_state[k] = v
    elif k == "clashscore":
      sm.model_state["geometry"]["clashscore"] = v
    elif k == "ligand_cc":
      # Create a ligand with this RSCC
      sm.model_state["ligands"] = [
        {"name": "ATP", "chain": "A",
         "resid": 301, "rscc": v},
      ]
    elif k == "diff_peaks":
      sm.model_state["diff_peaks"] = v
  return sm


def _make_hypothesis(**kwargs):
  """Build a hypothesis with defaults."""
  defaults = {
    "id": "h1",
    "statement": "Test hypothesis",
    "test_program": "phenix.refine",
    "confirm_if": "r_free < 0.30",
    "refute_if": "r_free > 0.45",
    "status": "proposed",
    "proposed_at_cycle": 1,
    "test_cycles_remaining": 1,
  }
  defaults.update(kwargs)
  return Hypothesis(**defaults)


# ── Criteria checking ─────────────────────────────────

def test_single_criterion_met():
  sm = _make_sm(r_free=0.25)
  met, details = _check_criteria(
    "r_free < 0.30", sm
  )
  assert met is True
  assert len(details) == 1
  assert details[0]["met"] is True


def test_single_criterion_not_met():
  sm = _make_sm(r_free=0.35)
  met, details = _check_criteria(
    "r_free < 0.30", sm
  )
  assert met is False


def test_compound_and_all_met():
  sm = _make_sm(r_free=0.25, clashscore=3.0)
  met, details = _check_criteria(
    "r_free < 0.30 AND clashscore < 5", sm
  )
  assert met is True
  assert len(details) == 2
  assert all(d["met"] for d in details)


def test_compound_and_partial():
  sm = _make_sm(r_free=0.25, clashscore=8.0)
  met, details = _check_criteria(
    "r_free < 0.30 AND clashscore < 5", sm
  )
  assert met is False
  met_count = sum(1 for d in details if d["met"])
  assert met_count == 1


def test_criterion_missing_metric():
  sm = StructureModel()  # no metrics at all
  met, details = _check_criteria(
    "r_free < 0.30", sm
  )
  assert met is False


def test_empty_criteria():
  sm = _make_sm(r_free=0.25)
  met, details = _check_criteria("", sm)
  assert met is False


# ── Lifecycle ─────────────────────────────────────────

def test_testing_to_pending():
  """testing → pending, decrement countdown."""
  sm = _make_sm(r_free=0.35)
  h = _make_hypothesis(
    status="testing",
    test_cycles_remaining=2,
  )
  sm.hypotheses = [h]
  results = evaluate_hypotheses(sm, 3)
  assert len(results) == 1
  assert results[0].action == "countdown"
  assert h.status == "pending"
  assert h.test_cycles_remaining == 1


def test_pending_countdown():
  """pending with remaining cycles → decrement."""
  sm = _make_sm(r_free=0.35)
  h = _make_hypothesis(
    status="pending",
    test_cycles_remaining=1,
  )
  sm.hypotheses = [h]
  results = evaluate_hypotheses(sm, 4)
  assert len(results) == 1
  assert results[0].action == "countdown"
  assert h.test_cycles_remaining == 0


def test_confirmed():
  """pending with 0 remaining → evaluate → confirmed."""
  sm = _make_sm(r_free=0.25)
  h = _make_hypothesis(
    status="pending",
    test_cycles_remaining=0,
    confirm_if="r_free < 0.30",
    refute_if="r_free > 0.45",
  )
  sm.hypotheses = [h]
  results = evaluate_hypotheses(sm, 5)
  assert len(results) == 1
  assert results[0].action == "confirmed"
  assert h.status == "confirmed"
  assert h.resolved_at_cycle == 5


def test_refuted():
  """pending with 0 remaining → evaluate → refuted."""
  sm = _make_sm(r_free=0.48)
  h = _make_hypothesis(
    status="pending",
    test_cycles_remaining=0,
    confirm_if="r_free < 0.30",
    refute_if="r_free > 0.45",
  )
  sm.hypotheses = [h]
  results = evaluate_hypotheses(sm, 5)
  assert len(results) == 1
  assert results[0].action == "refuted"
  assert h.status == "refuted"
  assert h.resolved_at_cycle == 5


def test_abandoned_timeout():
  """Inconclusive after timeout → abandoned."""
  sm = _make_sm(r_free=0.35)  # between thresholds
  h = _make_hypothesis(
    status="pending",
    test_cycles_remaining=0,
    confirm_if="r_free < 0.30",
    refute_if="r_free > 0.45",
    proposed_at_cycle=1,
  )
  sm.hypotheses = [h]
  # Cycle 10: well past timeout
  # Timeout = proposed(1) + remaining(1) +
  #   abandon(2) + 2 = 6, so cycle 10 > 6
  results = evaluate_hypotheses(sm, 10)
  assert len(results) == 1
  assert results[0].action == "abandoned"
  assert h.status == "abandoned"


def test_waiting_inconclusive():
  """Inconclusive but not yet timed out."""
  sm = _make_sm(r_free=0.35)
  h = _make_hypothesis(
    status="pending",
    test_cycles_remaining=0,
    confirm_if="r_free < 0.30",
    refute_if="r_free > 0.45",
    proposed_at_cycle=3,
  )
  sm.hypotheses = [h]
  # Cycle 5: only 2 cycles since proposed — not
  # past timeout yet
  results = evaluate_hypotheses(sm, 5)
  assert len(results) == 1
  assert results[0].action == "waiting"
  assert h.status == "pending"


def test_no_active_hypotheses():
  sm = _make_sm(r_free=0.30)
  # Only resolved hypotheses
  h = _make_hypothesis(status="confirmed")
  sm.hypotheses = [h]
  results = evaluate_hypotheses(sm, 5)
  assert len(results) == 0


def test_none_structure_model():
  results = evaluate_hypotheses(None, 5)
  assert results == []


# ── Revalidation ──────────────────────────────────────

def test_revalidate_still_valid():
  sm = _make_sm(r_free=0.25, ligand_cc=0.85)
  h = _make_hypothesis(
    status="confirmed",
    statement="ATP ligand is correctly placed",
    confirm_if="ligand_cc > 0.7",
  )
  sm.hypotheses = [h]
  results = revalidate_confirmed(sm)
  assert len(results) == 1
  assert results[0].action == "still_valid"
  assert h.status == "confirmed"


def test_revalidate_ligand_degraded():
  sm = _make_sm(r_free=0.25, ligand_cc=0.40)
  h = _make_hypothesis(
    status="confirmed",
    statement="ligand is correctly placed",
    confirm_if="ligand_cc > 0.7",
    resolved_at_cycle=5,
  )
  sm.hypotheses = [h]
  results = revalidate_confirmed(sm)
  assert len(results) == 1
  assert results[0].action == "revalidate"
  assert h.status == "proposed"
  assert "RSCC" in h.revalidation_reason


def test_revalidate_rfree_degraded():
  sm = _make_sm(r_free=0.50)
  h = _make_hypothesis(
    status="confirmed",
    statement="model is improved",
    confirm_if="r_free < 0.30",
  )
  sm.hypotheses = [h]
  results = revalidate_confirmed(sm)
  assert len(results) == 1
  assert results[0].action == "revalidate"
  assert h.status == "proposed"


def test_revalidate_no_confirmed():
  sm = _make_sm(r_free=0.25)
  results = revalidate_confirmed(sm)
  assert results == []


# ── Prompt building ───────────────────────────────────

def test_prompt_active_hypothesis():
  sm = _make_sm(r_free=0.35)
  h = _make_hypothesis(
    status="testing",
    statement="Difference density is Zn2+",
    test_program="phenix.refine",
    test_parameters={"add_ion": "Zn"},
    confirm_if="coordination > 4",
    refute_if="b_factor > 80",
    test_cycles_remaining=1,
  )
  sm.hypotheses = [h]
  prompt = build_hypothesis_prompt(sm)
  assert "ACTIVE HYPOTHESIS" in prompt
  assert "Zn2+" in prompt
  assert "DO NOT propose" in prompt
  assert "1 cycle(s)" in prompt


def test_prompt_revalidation_reason():
  sm = _make_sm(r_free=0.35)
  h = _make_hypothesis(
    status="testing",
    statement="Metal ion placement",
    revalidation_reason="B-factor rose to 90",
  )
  sm.hypotheses = [h]
  prompt = build_hypothesis_prompt(sm)
  assert "RE-EVALUATING" in prompt
  assert "B-factor" in prompt


def test_prompt_opportunity():
  sm = _make_sm(r_free=0.35)
  sm.model_state["diff_peaks"] = {
    "positive": [
      {"height": 6.5,
       "near_residue": "A/His47"},
      {"height": 5.2,
       "near_residue": "B/Glu100"},
    ],
    "negative": [],
  }
  prompt = build_hypothesis_prompt(sm)
  assert "HYPOTHESIS OPPORTUNITY" in prompt
  assert "His47" in prompt
  assert "test_program" in prompt


def test_prompt_nothing_to_hypothesize():
  sm = _make_sm(r_free=0.25)
  # No problems, no peaks
  prompt = build_hypothesis_prompt(sm)
  assert prompt == ""


def test_prompt_none_sm():
  prompt = build_hypothesis_prompt(None)
  assert prompt == ""


# ── Full scenarios ────────────────────────────────────

def test_ligand_hypothesis_lifecycle():
  """Full lifecycle: propose → test → confirm.

  Simulates: "ATP ligand should improve R-free"
  proposed at cycle 3, tested by ligandfit at 4,
  evaluated at cycle 5.
  """
  sm = _make_sm(r_free=0.32)

  # Propose
  h = Hypothesis(
    id="h_atp",
    statement="ATP ligand fits into density",
    test_program="phenix.ligandfit",
    test_parameters={"ligand_code": "ATP"},
    confirm_if="ligand_cc > 0.7",
    refute_if="ligand_cc < 0.4",
    proposed_at_cycle=3,
    test_cycles_remaining=1,
  )
  h.status = "testing"
  assert sm.add_hypothesis(h) is True

  # Cycle 4: test running (testing → pending)
  r1 = evaluate_hypotheses(sm, 4)
  assert len(r1) == 1
  assert r1[0].action == "countdown"
  assert h.status == "pending"
  assert h.test_cycles_remaining == 0

  # Cycle 5: evaluate — ligand placed, CC good
  sm.model_state["ligands"] = [
    {"name": "ATP", "chain": "A",
     "resid": 301, "rscc": 0.85},
  ]
  r2 = evaluate_hypotheses(sm, 5)
  assert len(r2) == 1
  assert r2[0].action == "confirmed"
  assert h.status == "confirmed"
  assert h.resolved_at_cycle == 5

  # Verify active slot is now free
  assert not sm.has_active_hypothesis()


def test_metal_ion_with_latency():
  """Hypothesis with 2-cycle latency.

  Metal ion placement needs 1 cycle to refine + 1
  cycle to stabilize B-factors before evaluation.
  """
  sm = _make_sm(r_free=0.30)

  h = Hypothesis(
    id="h_zn",
    statement="Difference peak near His47 is Zn",
    test_program="phenix.refine",
    test_parameters={"add_ion": "Zn"},
    confirm_if="clashscore < 5",
    refute_if="clashscore > 20",
    proposed_at_cycle=5,
    test_cycles_remaining=2,
  )
  h.status = "testing"
  sm.hypotheses = [h]

  # Cycle 6: test running → pending, countdown 2→1
  r1 = evaluate_hypotheses(sm, 6)
  assert r1[0].action == "countdown"
  assert h.test_cycles_remaining == 1
  assert h.status == "pending"

  # Cycle 7: still stabilizing → countdown 1→0
  r2 = evaluate_hypotheses(sm, 7)
  assert r2[0].action == "countdown"
  assert h.test_cycles_remaining == 0

  # Cycle 8: NOW evaluate — clashscore is good
  sm.model_state["geometry"]["clashscore"] = 3.5
  r3 = evaluate_hypotheses(sm, 8)
  assert r3[0].action == "confirmed"
  assert h.status == "confirmed"


# ── Entry point ──────────────────────────────────────

if __name__ == "__main__":
  run_tests()
