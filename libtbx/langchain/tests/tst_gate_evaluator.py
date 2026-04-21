"""
Unit tests for gate_evaluator.py (Stage 3, Steps 3.1-3.3).

Run standalone:
  python tests/tst_gate_evaluator.py

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

from agent.gate_evaluator import (
  GateEvaluator,
  parse_criterion, evaluate_criterion,
  apply_hysteresis,
  parse_gate_condition, evaluate_gate_condition,
)
from knowledge.plan_schema import (
  StageDef, StructurePlan,
  STAGE_ACTIVE, STAGE_COMPLETE,
)
from agent.structure_model import StructureModel


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
  print("Gate Evaluator Unit Tests")
  print("=" * 60)
  print()

  # --- Criterion parsing ---
  print("Criterion parsing")
  test("parse_less_than", test_parse_less_than)
  test("parse_greater_than", test_parse_greater_than)
  test("parse_less_equal", test_parse_less_equal)
  test("parse_boolean_true", test_parse_boolean_true)
  test("parse_boolean_false",
    test_parse_boolean_false)
  test("parse_empty", test_parse_empty)
  test("parse_garbage", test_parse_garbage)
  print()

  # --- Criterion evaluation ---
  print("Criterion evaluation")
  test("eval_less_than", test_eval_less_than)
  test("eval_greater_than", test_eval_greater_than)
  test("eval_less_equal", test_eval_less_equal)
  test("eval_equal", test_eval_equal)
  test("eval_boolean", test_eval_boolean)
  test("eval_none_actual", test_eval_none_actual)
  print()

  # --- Hysteresis ---
  print("Hysteresis")
  test("hysteresis_less_than",
    test_hysteresis_less_than)
  test("hysteresis_greater_than",
    test_hysteresis_greater_than)
  test("hysteresis_boolean",
    test_hysteresis_boolean)
  test("hysteresis_boundary",
    test_hysteresis_boundary)
  print()

  # --- Gate condition parsing ---
  print("Gate condition parsing")
  test("gate_parse_simple",
    test_gate_parse_simple)
  test("gate_parse_after_cycles",
    test_gate_parse_after_cycles)
  test("gate_parse_range",
    test_gate_parse_range)
  test("gate_parse_empty",
    test_gate_parse_empty)
  print()

  # --- Gate condition evaluation ---
  print("Gate condition evaluation")
  test("gate_eval_fires",
    test_gate_eval_fires)
  test("gate_eval_after_cycles_blocks",
    test_gate_eval_after_cycles_blocks)
  test("gate_eval_after_cycles_fires",
    test_gate_eval_after_cycles_fires)
  test("gate_eval_range",
    test_gate_eval_range)
  print()

  # --- Success checking ---
  print("Success checking")
  test("success_all_met",
    test_success_all_met)
  test("success_partial",
    test_success_partial)
  test("success_none_met",
    test_success_none_met)
  test("success_hysteresis_blocks",
    test_success_hysteresis_blocks)
  test("success_boolean",
    test_success_boolean)
  print()

  # --- Skip conditions ---
  print("Skip conditions")
  test("skip_met", test_skip_met)
  test("skip_not_met", test_skip_not_met)
  test("skip_no_data", test_skip_no_data)
  print()

  # --- Stage exhaustion ---
  print("Stage exhaustion")
  test("exhausted_partial_success",
    test_exhausted_partial_success)
  test("exhausted_no_success",
    test_exhausted_no_success)
  print()

  # --- Retreat ---
  print("Retreat")
  test("retreat_basic",
    test_retreat_basic)
  test("retreat_blocked_by_blacklist",
    test_retreat_blocked_by_blacklist)
  test("retreat_blocked_by_counter",
    test_retreat_blocked_by_counter)
  test("retreat_blocked_by_cooldown",
    test_retreat_blocked_by_cooldown)
  test("retreat_depth_limit",
    test_retreat_depth_limit)
  test("retreat_explicit_deep_target",
    test_retreat_explicit_deep_target)
  test("retreat_blocked_by_progress",
    test_retreat_blocked_by_progress)
  test("retreat_allowed_when_worsening",
    test_retreat_allowed_when_worsening)
  print()

  # --- Full evaluation scenarios ---
  print("Full evaluation scenarios")
  test("scenario_continue",
    test_scenario_continue)
  test("scenario_advance",
    test_scenario_advance)
  test("scenario_advance_skip",
    test_scenario_advance_skip)
  test("scenario_skip_preserves_status",
    test_scenario_skip_preserves_status)
  test("scenario_retreat",
    test_scenario_retreat)
  test("scenario_stop",
    test_scenario_stop)
  test("scenario_no_plan",
    test_scenario_no_plan)
  test("scenario_plan_complete",
    test_scenario_plan_complete)
  test("scenario_no_structure_model",
    test_scenario_no_structure_model)
  test("scenario_dict_structure_model",
    test_scenario_dict_structure_model)
  print()

  # --- Anti-oscillation scenario ---
  print("Anti-oscillation")
  test("oscillation_prevented",
    test_oscillation_prevented)
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

def _make_plan():
  """Build a realistic MR+refine plan."""
  stages = [
    StageDef(
      id="data_assessment",
      programs=["phenix.xtriage"],
      max_cycles=1,
      success_criteria={
        "xtriage_completed": "true",
      },
    ),
    StageDef(
      id="molecular_replacement",
      programs=["phenix.phaser"],
      max_cycles=1,
      success_criteria={"tfz": ">8", "llg": ">100"},
      gate_conditions=[
        {"if": "tfz < 5", "action": "stop_report"},
      ],
    ),
    StageDef(
      id="initial_refinement",
      programs=["phenix.refine"],
      max_cycles=3,
      success_criteria={"r_free": "<0.35"},
      gate_conditions=[{
        "if": "r_free > 0.45 after 2 cycles",
        "action":
          "retreat_to molecular_replacement",
      }],
    ),
    StageDef(
      id="model_rebuilding",
      programs=["phenix.autobuild"],
      max_cycles=1,
      success_criteria={"r_free": "<0.30"},
      skip_if="r_free < 0.28",
    ),
    StageDef(
      id="final_refinement",
      programs=["phenix.refine"],
      max_cycles=3,
      success_criteria={"r_free": "<0.25"},
    ),
  ]
  plan = StructurePlan(
    goal="Test plan", stages=stages,
    template_id="mr_refine",
  )
  return plan


def _make_sm(**kwargs):
  """Build a StructureModel with given metrics."""
  sm = StructureModel()
  for k, v in kwargs.items():
    if k in ("r_free", "r_work", "model_map_cc"):
      sm.model_state[k] = v
    elif k == "clashscore":
      sm.model_state["geometry"]["clashscore"] = v
    elif k == "rama_favored":
      sm.model_state["geometry"]["rama_favored"] = v
    elif k in ("tfz", "mr_tfz"):
      sm.data_characteristics["mr_tfz"] = v
    elif k in ("llg", "mr_llg"):
      sm.data_characteristics["mr_llg"] = v
    elif k == "resolution":
      sm.data_characteristics["resolution"] = v
  return sm


# ── Criterion parsing ─────────────────────────────────

def test_parse_less_than():
  op, val = parse_criterion("<0.35")
  assert op == "<"
  assert abs(val - 0.35) < 1e-9


def test_parse_greater_than():
  op, val = parse_criterion(">8")
  assert op == ">"
  assert abs(val - 8.0) < 1e-9


def test_parse_less_equal():
  op, val = parse_criterion("<=0.30")
  assert op == "<="
  assert abs(val - 0.30) < 1e-9


def test_parse_boolean_true():
  op, val = parse_criterion("true")
  assert op == "=="
  assert val is True


def test_parse_boolean_false():
  op, val = parse_criterion("false")
  assert op == "=="
  assert val is False


def test_parse_empty():
  op, val = parse_criterion("")
  assert op is None


def test_parse_garbage():
  op, val = parse_criterion("not_a_number")
  assert op is None


# ── Criterion evaluation ──────────────────────────────

def test_eval_less_than():
  assert evaluate_criterion(0.30, "<", 0.35) is True
  assert evaluate_criterion(0.35, "<", 0.35) is False
  assert evaluate_criterion(0.40, "<", 0.35) is False


def test_eval_greater_than():
  assert evaluate_criterion(10.0, ">", 8.0) is True
  assert evaluate_criterion(8.0, ">", 8.0) is False
  assert evaluate_criterion(5.0, ">", 8.0) is False


def test_eval_less_equal():
  assert evaluate_criterion(0.30, "<=", 0.30) is True
  assert evaluate_criterion(0.25, "<=", 0.30) is True
  assert evaluate_criterion(0.35, "<=", 0.30) is False


def test_eval_equal():
  assert evaluate_criterion(8.0, "==", 8.0) is True
  assert evaluate_criterion(7.9, "==", 8.0) is False


def test_eval_boolean():
  assert evaluate_criterion(True, "==", True) is True
  assert evaluate_criterion(False, "==", True) is False
  assert evaluate_criterion(1, "==", True) is True


def test_eval_none_actual():
  assert evaluate_criterion(None, "<", 0.35) is False


# ── Hysteresis ────────────────────────────────────────

def test_hysteresis_less_than():
  # <0.35 with 1.5% buffer → ~0.34475
  adj = apply_hysteresis(0.35, "<")
  assert adj < 0.35
  assert abs(adj - 0.34475) < 0.001


def test_hysteresis_greater_than():
  # >8 with 1.5% buffer → ~8.12
  adj = apply_hysteresis(8.0, ">")
  assert adj > 8.0
  assert abs(adj - 8.12) < 0.01


def test_hysteresis_boolean():
  # Boolean — no adjustment
  assert apply_hysteresis(True, "==") is True


def test_hysteresis_boundary():
  """A value that crosses the raw threshold but not
  the hysteresis threshold should NOT trigger."""
  # Raw threshold: <0.35, adjusted: ~0.34475
  # Value: 0.348 — below 0.35 but above 0.34475
  adj = apply_hysteresis(0.35, "<")
  assert evaluate_criterion(0.348, "<", adj) is False
  # Value: 0.340 — below both
  assert evaluate_criterion(0.340, "<", adj) is True


# ── Gate condition parsing ────────────────────────────

def test_gate_parse_simple():
  p = parse_gate_condition("tfz < 5")
  assert p is not None
  assert p["metric"] == "tfz"
  assert p["operator"] == "<"
  assert p["value"] == 5.0
  assert p["after_cycles"] is None


def test_gate_parse_after_cycles():
  p = parse_gate_condition(
    "r_free > 0.45 after 2 cycles"
  )
  assert p is not None
  assert p["metric"] == "r_free"
  assert p["operator"] == ">"
  assert p["value"] == 0.45
  assert p["after_cycles"] == 2


def test_gate_parse_range():
  p = parse_gate_condition(
    "r_free 0.35-0.45 after 3 cycles"
  )
  assert p is not None
  assert p["operator"] == "range"
  assert p["value"] == (0.35, 0.45)
  assert p["after_cycles"] == 3


def test_gate_parse_empty():
  assert parse_gate_condition("") is None
  assert parse_gate_condition(None) is None


# ── Gate condition evaluation ─────────────────────────

def test_gate_eval_fires():
  p = parse_gate_condition("tfz < 5")
  assert evaluate_gate_condition(p, 3.0, 1) is True
  assert evaluate_gate_condition(p, 8.0, 1) is False


def test_gate_eval_after_cycles_blocks():
  p = parse_gate_condition(
    "r_free > 0.45 after 2 cycles"
  )
  # Only 1 cycle — should NOT fire
  assert evaluate_gate_condition(
    p, 0.50, 1
  ) is False


def test_gate_eval_after_cycles_fires():
  p = parse_gate_condition(
    "r_free > 0.45 after 2 cycles"
  )
  # 2 cycles — should fire
  assert evaluate_gate_condition(
    p, 0.50, 2
  ) is True
  # Value below threshold — no fire
  assert evaluate_gate_condition(
    p, 0.40, 2
  ) is False


def test_gate_eval_range():
  p = parse_gate_condition("r_free 0.35-0.45")
  assert evaluate_gate_condition(
    p, 0.40, 1
  ) is True
  assert evaluate_gate_condition(
    p, 0.30, 1
  ) is False
  assert evaluate_gate_condition(
    p, 0.50, 1
  ) is False


# ── Success checking ──────────────────────────────────

def test_success_all_met():
  ev = GateEvaluator()
  sm = _make_sm(r_free=0.30)
  all_met, details = ev._check_success(
    {"r_free": "<0.35"}, sm
  )
  assert all_met is True
  assert details[0]["met"] is True


def test_success_partial():
  ev = GateEvaluator()
  sm = _make_sm(tfz=10.0, llg=80.0)
  all_met, details = ev._check_success(
    {"tfz": ">8", "llg": ">100"}, sm
  )
  assert all_met is False
  met_names = [d["criterion"] for d in details
               if d["met"]]
  assert "tfz" in met_names
  assert "llg" not in met_names


def test_success_none_met():
  ev = GateEvaluator()
  sm = _make_sm(r_free=0.45)
  all_met, details = ev._check_success(
    {"r_free": "<0.35"}, sm
  )
  assert all_met is False


def test_success_hysteresis_blocks():
  """Value crosses raw threshold but not hysteresis."""
  ev = GateEvaluator()
  # Threshold <0.35, hysteresis ~0.34475
  # Value 0.348 — below 0.35 but above 0.34475
  sm = _make_sm(r_free=0.348)
  all_met, details = ev._check_success(
    {"r_free": "<0.35"}, sm
  )
  assert all_met is False
  # Value 0.340 — below hysteresis
  sm2 = _make_sm(r_free=0.340)
  all_met2, _ = ev._check_success(
    {"r_free": "<0.35"}, sm2
  )
  assert all_met2 is True


def test_success_boolean():
  ev = GateEvaluator()
  # Boolean criteria with dict "structure model"
  all_met, _ = ev._check_success(
    {"xtriage_completed": "true"},
    {"xtriage_completed": True},
  )
  assert all_met is True


# ── Skip conditions ───────────────────────────────────

def test_skip_met():
  ev = GateEvaluator()
  stage = StageDef(
    id="test", skip_if="r_free < 0.28",
  )
  sm = _make_sm(r_free=0.25)
  result = ev._check_skip(stage, sm)
  assert result is not None
  assert result.action == "skip"


def test_skip_not_met():
  ev = GateEvaluator()
  stage = StageDef(
    id="test", skip_if="r_free < 0.28",
  )
  sm = _make_sm(r_free=0.35)
  result = ev._check_skip(stage, sm)
  assert result is None


def test_skip_no_data():
  ev = GateEvaluator()
  stage = StageDef(
    id="test", skip_if="r_free < 0.28",
  )
  sm = StructureModel()  # no r_free
  result = ev._check_skip(stage, sm)
  assert result is None


# ── Stage exhaustion ──────────────────────────────────

def test_exhausted_partial_success():
  ev = GateEvaluator()
  plan = _make_plan()
  # Set up initial_refinement as active + exhausted
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.stages[2].cycles_used = 3
  sm = _make_sm(r_free=0.36)  # close but not <0.35
  # But r_free itself is evaluable — partially met
  # Actually 0.36 does NOT meet <0.35, so 0/1 met
  # Let me use multi-criteria where one is met
  plan.stages[2].success_criteria = {
    "r_free": "<0.40",  # met
    "clashscore": "<5",  # not met (no data)
  }
  sm2 = _make_sm(r_free=0.36)
  result = ev._handle_exhaustion(
    plan, plan.stages[2], sm2, 5,
  )
  assert result.action == "advance"
  assert "1/2 criteria met" in result.reason


def test_exhausted_no_success():
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.stages[2].cycles_used = 3
  sm = _make_sm(r_free=0.45)
  result = ev._handle_exhaustion(
    plan, plan.stages[2], sm, 5,
  )
  assert result.action == "advance"
  assert "All steps completed" in result.reason


# ── Retreat ───────────────────────────────────────────

def test_retreat_basic():
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  result = ev._evaluate_retreat(
    plan, "molecular_replacement",
    "r_free > 0.45", 0.48, 5, None,
  )
  assert result.action == "retreat"
  assert result.new_phase_id == (
    "molecular_replacement"
  )
  assert result.blacklist_entry is not None


def test_retreat_blocked_by_blacklist():
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  sm = _make_sm(r_free=0.48)
  sm.blacklist_strategy(
    "molecular_replacement",
    "previous attempt failed",
  )
  result = ev._evaluate_retreat(
    plan, "molecular_replacement",
    "r_free > 0.45", 0.48, 5, sm,
  )
  assert result.action == "continue"
  assert "blacklisted" in result.reason


def test_retreat_blocked_by_counter():
  ev = GateEvaluator()
  plan = _make_plan()
  plan.retreat_count = 2  # max reached
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  result = ev._evaluate_retreat(
    plan, "molecular_replacement",
    "r_free > 0.45", 0.48, 10, None,
  )
  assert result.action == "continue"
  assert "max retreats" in result.reason


def test_retreat_blocked_by_cooldown():
  ev = GateEvaluator()
  plan = _make_plan()
  plan.retreat_count = 1
  plan.last_retreat_cycle = 4
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  # Cycle 5 — only 1 cycle since last retreat
  result = ev._evaluate_retreat(
    plan, "molecular_replacement",
    "r_free > 0.45", 0.48, 5, None,
  )
  assert result.action == "continue"
  assert "cooldown" in result.reason


def test_retreat_depth_limit():
  """Without explicit_target, retreat is limited
  to one stage back."""
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2  # initial_refinement
  plan.stages[2].status = STAGE_ACTIVE
  # Try to retreat to data_assessment (2 stages back)
  result = ev._evaluate_retreat(
    plan, "data_assessment",
    "r_free > 0.50", 0.52, 5, None,
    explicit_target=False,
  )
  assert result.action == "retreat"
  # Should be clamped to one-back
  assert result.new_phase_id == (
    "molecular_replacement"
  ), (
    "Expected molecular_replacement, got %s"
    % result.new_phase_id
  )


def test_retreat_explicit_deep_target():
  """With explicit_target, depth limit is bypassed."""
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  result = ev._evaluate_retreat(
    plan, "data_assessment",
    "twinning discovered", None, 5, None,
    explicit_target=True,
  )
  assert result.action == "retreat"
  assert result.new_phase_id == "data_assessment"


def test_retreat_blocked_by_progress():
  """Retreat blocked when metric is improving since
  stage start (monotonic progress gate)."""
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.stages[2].start_cycle = 3
  plan.stages[2].cycles_used = 2

  sm = _make_sm(r_free=0.46)
  # Simulate progress: stage started at r_free=0.50
  sm.progress = [
    {"cycle": 3, "r_free": 0.50,
     "r_work": None, "model_map_cc": None},
    {"cycle": 4, "r_free": 0.46,
     "r_work": None, "model_map_cc": None},
  ]

  result = ev._evaluate_retreat(
    plan, "molecular_replacement",
    "r_free > 0.45 after 2 cycles", 0.46,
    5, sm,
    explicit_target=True,
  )
  assert result.action == "continue", (
    "Expected continue (improving), got %s: %s"
    % (result.action, result.reason)
  )
  assert "improving" in result.reason


def test_retreat_allowed_when_worsening():
  """Retreat allowed when metric is worse than at
  stage start."""
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.stages[2].start_cycle = 3
  plan.stages[2].cycles_used = 2

  sm = _make_sm(r_free=0.48)
  # Stage started at r_free=0.45 — WORSENED
  sm.progress = [
    {"cycle": 3, "r_free": 0.45,
     "r_work": None, "model_map_cc": None},
    {"cycle": 4, "r_free": 0.48,
     "r_work": None, "model_map_cc": None},
  ]

  result = ev._evaluate_retreat(
    plan, "molecular_replacement",
    "r_free > 0.45 after 2 cycles", 0.48,
    5, sm,
    explicit_target=True,
  )
  assert result.action == "retreat", (
    "Expected retreat (worsening), got %s: %s"
    % (result.action, result.reason)
  )


# ── Full evaluation scenarios ─────────────────────────

def test_scenario_continue():
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.stages[2].cycles_used = 1
  sm = _make_sm(r_free=0.40)  # improving but not there
  result = ev.evaluate(plan, sm, None, 3)
  assert result.action == "continue"


def test_scenario_advance():
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.stages[2].cycles_used = 2
  sm = _make_sm(r_free=0.30)  # below <0.35 threshold
  result = ev.evaluate(plan, sm, None, 4)
  assert result.action == "advance"
  assert "success criteria met" in result.reason


def test_scenario_advance_skip():
  """Stage with skip_if met should skip."""
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 3  # model_rebuilding
  plan.stages[3].status = STAGE_ACTIVE
  sm = _make_sm(r_free=0.25)  # < 0.28, skip_if met
  result = ev.evaluate(plan, sm, None, 6)
  assert result.action == "skip"
  assert "skip_if" in result.reason


def test_scenario_skip_preserves_status():
  """After skip + advance, stage should be SKIPPED
  not COMPLETE."""
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 3  # model_rebuilding
  plan.stages[3].status = STAGE_ACTIVE
  sm = _make_sm(r_free=0.25)
  result = ev.evaluate(plan, sm, None, 6)
  assert result.action == "skip"
  # Simulate the caller's handling
  plan.skip_stage(plan.stages[3].id)
  plan.advance()
  # Stage 3 should be SKIPPED, not COMPLETE
  assert plan.stages[3].status == "skipped", (
    "Expected skipped, got %s"
    % plan.stages[3].status
  )
  # Should have advanced to final_refinement
  curr = plan.current_stage()
  assert curr is not None
  assert curr.id == "final_refinement"


def test_scenario_retreat():
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.stages[2].cycles_used = 2
  sm = _make_sm(r_free=0.48)  # > 0.45 after 2 cycles
  result = ev.evaluate(plan, sm, None, 4)
  assert result.action == "retreat"
  assert result.new_phase_id == (
    "molecular_replacement"
  )


def test_scenario_stop():
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 1  # MR stage
  plan.stages[1].status = STAGE_ACTIVE
  plan.stages[1].cycles_used = 1
  sm = _make_sm(tfz=3.0)  # < 5 → stop_report
  result = ev.evaluate(plan, sm, None, 2)
  assert result.action == "stop"
  assert "tfz" in result.reason.lower()


def test_scenario_no_plan():
  ev = GateEvaluator()
  result = ev.evaluate(None, None, None, 1)
  assert result.action == "continue"
  assert "no plan" in result.reason


def test_scenario_plan_complete():
  ev = GateEvaluator()
  plan = _make_plan()
  for p in plan.stages:
    p.status = STAGE_COMPLETE
  result = ev.evaluate(
    plan, _make_sm(r_free=0.22), None, 10
  )
  assert result.action == "stop"
  assert "complete" in result.reason


def test_scenario_no_structure_model():
  """Without a structure model, gate evaluator
  should continue (cannot check criteria)."""
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.stages[2].cycles_used = 1
  result = ev.evaluate(plan, None, None, 3)
  assert result.action == "continue"


def test_scenario_dict_structure_model():
  """Gate evaluator works with plain dict."""
  ev = GateEvaluator()
  plan = _make_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  sm_dict = {"r_free": 0.30}
  result = ev.evaluate(plan, sm_dict, None, 4)
  assert result.action == "advance"


# ── Anti-oscillation scenario ─────────────────────────

def test_oscillation_prevented():
  """Simulate retreat → retry → retreat cycle.
  Second retreat should be blocked."""
  ev = GateEvaluator()
  plan = _make_plan()

  # First attempt: refinement stalls
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.stages[2].cycles_used = 2
  sm = _make_sm(r_free=0.48)

  r1 = ev.evaluate(plan, sm, None, 4)
  assert r1.action == "retreat"

  # Execute retreat
  plan.retreat_to(
    r1.new_phase_id, cycle_number=4
  )

  # Second attempt: redo MR and refine, stalls again
  plan.advance()  # back to refinement
  plan.stages[2].status = STAGE_ACTIVE
  plan.stages[2].cycles_used = 2

  r2 = ev.evaluate(plan, sm, None, 7)
  assert r2.action == "retreat"

  # Execute second retreat
  plan.retreat_to(
    r2.new_phase_id, cycle_number=7
  )

  # Third attempt: should be blocked
  plan.advance()
  plan.stages[2].status = STAGE_ACTIVE
  plan.stages[2].cycles_used = 2

  r3 = ev.evaluate(plan, sm, None, 10)
  # Max retreats reached (2) — should continue
  # instead of retreating
  assert r3.action != "retreat", (
    "Expected non-retreat, got %s: %s"
    % (r3.action, r3.reason)
  )


# ── Entry point ──────────────────────────────────────

if __name__ == "__main__":
  run_tests()
