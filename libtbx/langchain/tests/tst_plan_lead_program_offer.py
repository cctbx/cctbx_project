"""Regression test: Option 2a parts 2+3 — surface and offer the current active
stage's un-run lead program.

After plan_schema HOLDS model_rebuilding active (part 1), the agent must actually
be OFFERED phenix.autobuild next cycle, or it will just pick the workflow engine's
default program again.  Two pieces:

- ai_agent._get_plan_current_unrun_lead_program: returns the current ACTIVE
  stage's lead program when it has not run (status active, cycles_used 0, lead not
  in history).  `_get_plan_next_stage_programs` misses this because it only looks
  at status=="pending".

- graph_nodes PERCEIVE: injects that lead program into valid_programs even when
  other non-STOP programs are already valid (the prior injection only fired when
  valid_programs was STOP-only).

These tests exercise the decision LOGIC of both pieces in isolation (the real
functions pull in heavy graph/session machinery).  The logic mirrors the deployed
code and is kept in sync.  2-space indent.
"""

from __future__ import absolute_import, division, print_function


def _current_unrun_lead(plan_data, cycles):
  """Mirror of ai_agent._get_plan_current_unrun_lead_program."""
  if not plan_data or not isinstance(plan_data, dict):
    return ""
  stages = (plan_data.get("stages")
            or plan_data.get("phases") or [])
  idx = plan_data.get("current_stage_index")
  if not isinstance(idx, int) or idx < 0 or idx >= len(stages):
    return ""
  stage = stages[idx]
  if not isinstance(stage, dict):
    return ""
  if stage.get("status") != "active":
    return ""
  if int(stage.get("cycles_used", 0) or 0) != 0:
    return ""
  progs = stage.get("programs", []) or []
  if not progs:
    return ""
  lead = progs[0]
  for h in (cycles or []):
    if isinstance(h, dict):
      prog = h.get("program") or ""
      result = h.get("result") or ""
      if prog and "SUCCESS" in str(result).upper():
        if prog == lead or prog.startswith(lead + "_"):
          return ""
  return lead


def _perceive_inject(valid_programs, lead):
  """Mirror of the graph_nodes PERCEIVE lead-program injection (part 2)."""
  valid = list(valid_programs)
  if lead and lead not in valid:
    valid.append(lead)
  return valid


_PLAN = {
  "stages": [
    {"id": "refinement", "programs": ["phenix.refine"],
     "status": "complete", "cycles_used": 1},
    {"id": "model_rebuilding",
     "programs": ["phenix.autobuild", "phenix.refine"],
     "status": "active", "cycles_used": 0},
    {"id": "validation", "programs": ["phenix.molprobity"],
     "status": "pending", "cycles_used": 0},
  ],
  "current_stage_index": 1,
}
_CYC = [
  {"program": "phenix.xtriage", "result": "SUCCESS"},
  {"program": "phenix.refine", "result": "SUCCESS"},
]


def test_lead_is_autobuild_when_active_and_unrun():
  assert _current_unrun_lead(_PLAN, _CYC) == "phenix.autobuild"


def test_no_lead_when_autobuild_already_ran():
  cyc = _CYC + [{"program": "phenix.autobuild", "result": "SUCCESS"}]
  assert _current_unrun_lead(_PLAN, cyc) == ""


def test_no_lead_when_variant_already_ran():
  cyc = _CYC + [{"program": "phenix.autobuild_denmod", "result": "SUCCESS"}]
  assert _current_unrun_lead(_PLAN, cyc) == ""


def test_no_lead_when_stage_pending_not_active():
  plan = {"stages": [dict(s) for s in _PLAN["stages"]],
          "current_stage_index": 1}
  plan["stages"][1] = dict(plan["stages"][1])
  plan["stages"][1]["status"] = "pending"
  assert _current_unrun_lead(plan, _CYC) == ""


def test_no_lead_when_cycles_used_nonzero():
  plan = {"stages": [dict(s) for s in _PLAN["stages"]],
          "current_stage_index": 1}
  plan["stages"][1] = dict(plan["stages"][1])
  plan["stages"][1]["cycles_used"] = 1
  assert _current_unrun_lead(plan, _CYC) == ""


def test_perceive_injects_lead_alongside_other_programs():
  """The key fix: autobuild is added even when molprobity is already valid."""
  lead = _current_unrun_lead(_PLAN, _CYC)
  valid = _perceive_inject(["phenix.molprobity", "STOP"], lead)
  assert "phenix.autobuild" in valid, "lead must be offered: %r" % valid
  assert "phenix.molprobity" in valid, "existing programs preserved"


def test_perceive_noop_when_no_lead():
  valid = _perceive_inject(["phenix.molprobity", "STOP"], "")
  assert valid == ["phenix.molprobity", "STOP"]


def test_perceive_no_duplicate():
  valid = _perceive_inject(["phenix.autobuild", "STOP"], "phenix.autobuild")
  assert valid.count("phenix.autobuild") == 1


_TESTS = [
  test_lead_is_autobuild_when_active_and_unrun,
  test_no_lead_when_autobuild_already_ran,
  test_no_lead_when_variant_already_ran,
  test_no_lead_when_stage_pending_not_active,
  test_no_lead_when_cycles_used_nonzero,
  test_perceive_injects_lead_alongside_other_programs,
  test_perceive_noop_when_no_lead,
  test_perceive_no_duplicate,
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
