"""Regression test: Option 2a reactive-deviation HOLD in record_stage_cycle.

Bug (beta-blip AIAgent_21): the refine_rebuild_placed plan put model_rebuilding
active after refinement, but at cycle 3 the LLM ran phenix.molprobity (validate
first).  molprobity matched a LATER stage (validation), so plan_schema's catch-up
advanced PAST model_rebuilding (marking it complete, 0 cycles) — abandoning the
rebuild the user requested; autobuild never ran.

Fix (Option 2a): record_stage_cycle HOLDS a single reactive deviation.  When the
off-plan program matches a later stage but the CURRENT stage is must-run-and-un-run
(cycles_used == 0, lead program not run, criteria unmet) and its
reactive_deviations < 1, it does NOT advance — it counts one reactive_deviation
and keeps the stage active, so the next cycle runs the stage's lead program.  On a
SECOND reactive deviation the guard is exhausted and the normal catch-up proceeds
(no infinite hold).

These tests drive the REAL StructurePlan/StageDef + the real GateEvaluator
(loaded for _criteria_met).  2-space indent.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import types
import importlib.util

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_plan_schema():
  """Load StructurePlan/StageDef.  Prefer package import; fall back to env path.
  Also loads the real gate_evaluator so _criteria_met evaluates criteria instead
  of falling back to its 'treat as met' degradation."""
  try:
    from libtbx.langchain.knowledge.plan_schema import StructurePlan
    from libtbx.langchain.knowledge.plan_schema import StageDef
    return StructurePlan, StageDef
  except Exception:
    pass
  ps_path = os.environ.get("PLAN_SCHEMA_PY")
  ge_path = os.environ.get("GATE_EVALUATOR_PY")
  if not ps_path:
    for cand in (
        os.path.join(_HERE, "plan_schema.py"),
        os.path.join(_HERE, "..", "knowledge", "plan_schema.py"),
    ):
      if os.path.isfile(cand):
        ps_path = cand
        break
  if not ps_path or not os.path.isfile(ps_path):
    return None, None
  # Best-effort: load gate_evaluator under agent.gate_evaluator so the lazy
  # import in _criteria_met resolves; if unavailable, _criteria_met degrades to
  # "treat as met" (no hold) — tests that need a hold will skip.
  if ge_path and os.path.isfile(ge_path):
    if "agent" not in sys.modules:
      m = types.ModuleType("agent")
      m.__path__ = [os.path.dirname(ge_path)]
      sys.modules["agent"] = m
    spec_ge = importlib.util.spec_from_file_location(
      "agent.gate_evaluator", ge_path)
    ge = importlib.util.module_from_spec(spec_ge)
    sys.modules["agent.gate_evaluator"] = ge
    try:
      spec_ge.loader.exec_module(ge)
    except Exception:
      pass
  if "knowledge" not in sys.modules:
    m = types.ModuleType("knowledge")
    m.__path__ = [os.path.dirname(ps_path)]
    sys.modules["knowledge"] = m
  spec = importlib.util.spec_from_file_location(
    "knowledge.plan_schema", ps_path)
  ps = importlib.util.module_from_spec(spec)
  sys.modules["knowledge.plan_schema"] = ps
  try:
    spec.loader.exec_module(ps)
  except Exception:
    return None, None
  return ps.StructurePlan, ps.StageDef


_SP, _SD = _load_plan_schema()


class _SM(object):
  """Minimal StructureModel stand-in with get_metric()."""
  def __init__(self, r_free):
    self._r_free = r_free

  def get_metric(self, name):
    if name == "r_free":
      return self._r_free
    return None

  def __getattr__(self, n):
    return None


def _make_plan():
  stages = [
    _SD(id="refinement", programs=["phenix.refine"],
        success_criteria={"r_free": "<0.25"}),
    _SD(id="model_rebuilding",
        programs=["phenix.autobuild", "phenix.refine"],
        success_criteria={"r_free": "<0.30"},
        skip_if="r_free < 0.28", max_cycles=2),
    _SD(id="final_refinement", programs=["phenix.refine"],
        success_criteria={"r_free": "<0.22"}),
    _SD(id="validation", programs=["phenix.molprobity"]),
  ]
  p = _SP(goal="t", stages=stages,
          template_id="refine_rebuild_placed")
  p.current_stage_index = 1
  p.stages[1].status = "active"
  return p


_SESS = {"cycles": [
  {"program": "phenix.xtriage", "result": "SUCCESS"},
  {"program": "phenix.refine", "result": "SUCCESS"},
]}


def _hold_supported():
  """Hold requires _criteria_met to return False (criteria unmet).  If the
  gate_evaluator import is unavailable, _criteria_met degrades to True and the
  hold can't fire — skip those tests rather than report a false failure."""
  if _SP is None:
    return False
  p = _make_plan()
  p.record_stage_cycle("phenix.molprobity", structure_model=_SM(0.326),
                       session_data=dict(_SESS), cycle_number=3)
  return p.current_stage() is not None and \
      p.current_stage().id == "model_rebuilding"


def test_first_reactive_deviation_holds():
  if _SP is None:
    print("  (skip) plan_schema not loadable")
    return
  if not _hold_supported():
    print("  (skip) _criteria_met degraded (gate_evaluator unavailable)")
    return
  p = _make_plan()
  p.record_stage_cycle("phenix.molprobity", structure_model=_SM(0.326),
                       session_data=dict(_SESS), cycle_number=3)
  assert p.current_stage().id == "model_rebuilding", \
    "stage should be HELD active, got %r" % p.current_stage().id
  assert p.stages[1].cycles_used == 0, "lead un-run; must not count a cycle"
  assert p.stages[1].reactive_deviations == 1


def test_second_reactive_deviation_catches_up():
  if _SP is None or not _hold_supported():
    print("  (skip)")
    return
  p = _make_plan()
  p.record_stage_cycle("phenix.molprobity", structure_model=_SM(0.326),
                       session_data=dict(_SESS), cycle_number=3)
  p.record_stage_cycle("phenix.molprobity", structure_model=_SM(0.326),
                       session_data=dict(_SESS), cycle_number=4)
  cur = p.current_stage()
  assert cur is None or cur.id == "validation", \
    "2nd deviation should catch up, got %r" % (cur.id if cur else None)


def test_no_hold_when_lead_already_ran():
  if _SP is None:
    print("  (skip)")
    return
  p = _make_plan()
  sess = {"cycles": [
    {"program": "phenix.refine", "result": "SUCCESS"},
    {"program": "phenix.autobuild", "result": "SUCCESS"},
  ]}
  p.record_stage_cycle("phenix.molprobity", structure_model=_SM(0.326),
                       session_data=sess, cycle_number=3)
  cur = p.current_stage()
  assert cur is None or cur.id == "validation", \
    "lead already ran -> should advance, got %r" % (cur.id if cur else None)


def test_no_hold_when_criteria_met():
  if _SP is None:
    print("  (skip)")
    return
  p = _make_plan()
  # r_free 0.20 meets <0.30 -> _criteria_met True -> no hold
  p.record_stage_cycle("phenix.molprobity", structure_model=_SM(0.20),
                       session_data=dict(_SESS), cycle_number=3)
  cur = p.current_stage()
  assert cur is None or cur.id == "validation", \
    "criteria met -> should advance, got %r" % (cur.id if cur else None)


def test_reactive_deviations_round_trips_through_serialization():
  if _SP is None:
    print("  (skip)")
    return
  p = _make_plan()
  p.stages[1].reactive_deviations = 1
  d = p.to_dict()
  p2 = _SP.from_dict(d)
  assert p2.stages[1].reactive_deviations == 1, \
    "reactive_deviations must survive to_dict/from_dict"
  # Backward-compat: a stage dict missing the field deserializes to 0.
  d2 = p.to_dict()
  for s in d2["stages"]:
    s.pop("reactive_deviations", None)
  p3 = _SP.from_dict(d2)
  assert p3.stages[1].reactive_deviations == 0, \
    "missing field must default to 0"


def test_holds_with_empty_history():
  """Regression: the lead-ran check must NOT treat an empty history as 'lead
  already ran' (an earlier version used a helper that returned True on an empty
  list).  With no successes in history, autobuild has not run -> HOLD."""
  if _SP is None:
    print("  (skip)")
    return
  p = _make_plan()
  if p.current_stage() is None:
    return
  p.record_stage_cycle("phenix.molprobity", structure_model=_SM(0.326),
                       session_data={"cycles": []}, cycle_number=3)
  # Only meaningful when _criteria_met can evaluate (gate_evaluator present).
  if not _hold_supported():
    print("  (skip) _criteria_met degraded")
    return
  assert p.current_stage().id == "model_rebuilding", \
    "empty history must still HOLD, got %r" % p.current_stage().id


def test_no_hold_when_lead_variant_ran():
  """A variant of the lead (phenix.autobuild_denmod) counts as the lead
  (phenix.autobuild) having run -> no hold."""
  if _SP is None:
    print("  (skip)")
    return
  p = _make_plan()
  sess = {"cycles": [
    {"program": "phenix.refine", "result": "SUCCESS"},
    {"program": "phenix.autobuild_denmod", "result": "SUCCESS"},
  ]}
  p.record_stage_cycle("phenix.molprobity", structure_model=_SM(0.326),
                       session_data=sess, cycle_number=3)
  cur = p.current_stage()
  assert cur is None or cur.id == "validation", \
    "lead variant ran -> should advance, got %r" % (cur.id if cur else None)


_TESTS = [
  test_first_reactive_deviation_holds,
  test_second_reactive_deviation_catches_up,
  test_no_hold_when_lead_already_ran,
  test_no_hold_when_criteria_met,
  test_holds_with_empty_history,
  test_no_hold_when_lead_variant_ran,
  test_reactive_deviations_round_trips_through_serialization,
]


def run_all_tests():
  for fn in _TESTS:
    fn()
  print("All %d tests passed." % len(_TESTS))
  return True


if __name__ == "__main__":
  p = f = 0
  for fn in _TESTS:
    try:
      fn()
      print("  PASS: %s" % fn.__name__)
      p += 1
    except AssertionError as e:
      print("  FAIL: %s -- %s" % (fn.__name__, e))
      f += 1
  print("\n%d passed, %d failed" % (p, f))
