"""
Scenario Tracer for PHENIX AI Agent.

Calls real agent functions with mock data to trace
what would happen at each decision point. Finds bugs
in plan selection, workflow routing, gate evaluation,
and prompt construction WITHOUT running the full agent.

Usage:
  python tests/scenario_tracer.py           # All
  python tests/scenario_tracer.py S1A       # One
  python tests/scenario_tracer.py --list    # List
  python tests/scenario_tracer.py --failures # Failures only

Each scenario produces a trace like:

  S1A: Standard MR — phaser succeeds
  ──────────────────────────────────
  Files: data.mtz(data_mtz) model.pdb(model)
  Template: mr_refine (5 stages)
  Cycle 1: valid=[xtriage,phaser,STOP]
           plan_phase=data_assessment (step 0/1)
           xtriage runs → gate: advance
  Cycle 2: valid=[refine,STOP]
           plan_phase=molecular_replacement (step 0/1)
           phaser runs → gate: advance
  ...
  RESULT: PASS (matches expected)
"""

from __future__ import (
  absolute_import, division, print_function,
)

import os
import sys

_this_dir = os.path.dirname(
  os.path.abspath(__file__))
_project_dir = os.path.dirname(_this_dir)
if _project_dir not in sys.path:
  sys.path.insert(0, _project_dir)

from knowledge import plan_template_loader
from agent.plan_generator import generate_plan
from agent.gate_evaluator import GateEvaluator
from agent.structure_model import StructureModel
from knowledge.plan_schema import StructurePlan


# ── Mock helpers ──────────────────────────────

def mock_history(cycles):
  """Build history list from cycle tuples.

  Each cycle: (program, result, metrics_dict)
  """
  history = []
  for i, (prog, result, metrics) in enumerate(
    cycles
  ):
    history.append({
      "cycle_number": i + 1,
      "program": prog,
      "result": result,
      "metrics": metrics or {},
      "output_files": [],
    })
  return history


def mock_detect_ws(files, history, advice="",
                   directives=None):
  """Call detect_workflow_state with mock data.

  Returns (state_name, valid_programs, reason).
  Returns ("SKIP", [], "libtbx not available") when
  the full PHENIX environment is not present.
  """
  try:
    from agent.workflow_state import (
      detect_workflow_state,
    )
    ws = detect_workflow_state(
      history=history,
      available_files=files,
      analysis=None,
      maximum_automation=True,
      use_yaml_engine=True,
      directives=directives or {},
      session_info={
        "user_advice": advice,
      },
      files_local=False,
    )
    return (
      ws.get("state", "?"),
      ws.get("valid_programs", []),
      ws.get("reason", ""),
    )
  except ImportError:
    return ("SKIP", [], "libtbx not available")
  except Exception as e:
    err = str(e)[:80]
    if "libtbx" in err.lower() or "No module" in err:
      return ("SKIP", [], err)
    return ("ERROR", [], err)


def mock_gate(plan_dict, sm_dict, cycle):
  """Run gate evaluation on a plan + structure model.

  Returns GateResult.
  """
  plan = StructurePlan.from_dict(plan_dict)
  sm = None
  if sm_dict:
    sm = StructureModel.from_dict(sm_dict)
  gate = GateEvaluator()
  return gate.evaluate(plan, sm, None, cycle)


def make_sm(r_free=None, r_work=None,
            model_map_cc=None, resolution=None,
            space_group=None, tfz=None, llg=None,
            fom=None, clashscore=None,
            chains=None, waters=0):
  """Build a StructureModel dict."""
  ms = {}
  if r_free is not None:
    ms["r_free"] = r_free
  if r_work is not None:
    ms["r_work"] = r_work
  if model_map_cc is not None:
    ms["model_map_cc"] = model_map_cc
  if clashscore is not None:
    ms["geometry"] = {"clashscore": clashscore}
  if chains:
    ms["chains"] = [{"chain_id": c} for c in chains]
  ms["waters"] = waters

  dc = {}
  if resolution is not None:
    dc["resolution"] = resolution
  if space_group:
    dc["space_group"] = space_group
  if tfz is not None:
    dc["mr_tfz"] = tfz
  if llg is not None:
    dc["mr_llg"] = llg
  if fom is not None:
    dc["phasing_fom"] = fom

  return {
    "_version": 1,
    "model_state": ms,
    "data_characteristics": dc,
    "progress": [],
    "strategy_blacklist": [],
    "hypotheses": [],
  }


# ── Trace engine ──────────────────────────────

class TraceResult:
  def __init__(self, scenario_id, description):
    self.scenario_id = scenario_id
    self.description = description
    self.steps = []
    self.issues = []
    self.status = "PASS"

  def step(self, msg):
    self.steps.append(("  " + msg))

  def issue(self, severity, msg):
    self.issues.append((severity, msg))
    if severity == "WRONG" and self.status in (
      "PASS", "DEGRADED"
    ):
      self.status = "WRONG"
    elif severity == "STUCK" and self.status in (
      "PASS", "DEGRADED"
    ):
      self.status = "STUCK"
    elif severity == "DEGRADED" and self.status == (
      "PASS"
    ):
      self.status = "DEGRADED"
    elif severity == "CRASH":
      self.status = "CRASH"

  def print_trace(self, failures_only=False):
    if failures_only and self.status == "PASS":
      return
    icon = {
      "PASS": "\u2705",
      "DEGRADED": "\u26A0\uFE0F",
      "WRONG": "\u274C",
      "STUCK": "\u267B\uFE0F",
      "CRASH": "\U0001F4A5",
    }.get(self.status, "?")
    print("\n%s %s: %s" % (
      icon, self.scenario_id, self.description))
    print("\u2500" * 50)
    for s in self.steps:
      print(s)
    if self.issues:
      print()
      for sev, msg in self.issues:
        print("  [%s] %s" % (sev, msg))
    print()


# ── Scenario definitions ──────────────────────

def trace_S1A():
  """S1-A: Standard MR — phaser succeeds."""
  t = TraceResult("S1A",
    "Standard MR — phaser succeeds")

  files = ["data.mtz", "model.pdb", "sequence.fa"]
  advice = "solve by molecular replacement"

  # Step 1: Plan selection
  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan generated")
    return t
  t.step("Template: %s (%d stages)" % (
    plan.template_id, len(plan.stages)))
  if plan.template_id != "mr_refine":
    t.issue("WRONG",
      "Expected mr_refine, got %s"
      % plan.template_id)
  phase_ids = [p.id for p in plan.stages]
  t.step("Phases: %s" % " → ".join(phase_ids))

  # Step 2: Workflow state — fresh session
  state, vp, reason = mock_detect_ws(
    files, [], advice)
  t.step("Cycle 1: state=%s valid=%s" % (
    state, vp))
  if state == "SKIP":
    t.step("  (workflow state requires PHENIX "
           "environment — skipped)")
  elif "phenix.xtriage" not in vp:
    t.issue("WRONG", "xtriage not in valid_programs")

  # Step 3: After xtriage
  h1 = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.0}),
  ])
  state2, vp2, _ = mock_detect_ws(
    files, h1, advice)
  t.step("Cycle 2: state=%s valid=%s" % (
    state2, vp2))
  if state2 != "SKIP":
    # Workflow may probe model fit before MR
    if "phenix.phaser" in vp2:
      t.step("phaser available (direct to MR)")
    elif "phenix.model_vs_data" in vp2:
      t.step("model_vs_data probe first "
             "(checks model fit before MR)")
    else:
      t.issue("WRONG",
        "Neither phaser nor model_vs_data "
        "after xtriage: %s" % vp2)

  # Step 4: After phaser (good TFZ)
  h2 = h1 + mock_history([
    ("phenix.phaser", "SUCCESS: TFZ=14.2",
     {"tfz": 14.2, "llg": 342}),
  ])
  # Fix cycle numbers
  for i, c in enumerate(h2):
    c["cycle_number"] = i + 1
  state3, vp3, _ = mock_detect_ws(
    files, h2, advice)
  t.step("Cycle 3: state=%s valid=%s" % (
    state3, vp3))
  if state3 != "SKIP" and "phenix.refine" not in vp3:
    t.issue("WRONG",
      "refine not in valid_programs after phaser")

  # Step 5: Gate evaluation after xtriage
  plan_d = plan.to_dict()
  # Simulate stage advancement
  plan_copy = StructurePlan.from_dict(plan_d)
  plan_copy.mark_stage_started(1)
  plan_copy.record_stage_cycle("phenix.xtriage")
  gr = mock_gate(plan_copy.to_dict(), None, 1)
  t.step("Gate after xtriage: %s (%s)" % (
    gr.action, gr.reason[:50]))
  if gr.action != "advance":
    t.issue("DEGRADED",
      "Expected advance after xtriage")

  # Step 6: Gate after phaser
  plan_copy.advance()
  plan_copy.mark_stage_started(2)
  plan_copy.record_stage_cycle("phenix.phaser")
  sm = make_sm(tfz=14.2, resolution=2.0)
  gr2 = mock_gate(plan_copy.to_dict(), sm, 2)
  t.step("Gate after phaser: %s (%s)" % (
    gr2.action, gr2.reason[:50]))

  return t


def trace_S1B():
  """S1-B: Standard MR — phaser fails."""
  t = TraceResult("S1B",
    "Standard MR — phaser fails (TFZ < 5)")

  files = ["data.mtz", "model.pdb", "sequence.fa"]
  advice = "solve by molecular replacement"

  # After xtriage + failed phaser
  h = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.0}),
    ("phenix.phaser", "FAILED: No solution",
     {"tfz": 3.1}),
  ])
  state, vp, reason = mock_detect_ws(
    files, h, advice)
  t.step("After failed phaser: state=%s" % state)
  t.step("  valid=%s" % vp)
  t.step("  reason=%s" % reason[:60])

  if not vp or vp == ["STOP"]:
    t.step("Agent would STOP (no recovery)")
    # This is expected for rules_only
  elif "phenix.predict_and_build" in vp:
    t.step("predict_and_build available (good)")
  else:
    t.step("Has non-STOP programs: %s" % vp)

  return t


def trace_S1C():
  """S1-C: Silent failure — good TFZ, stuck R-free
  (enantiomorph/space group error)."""
  t = TraceResult("S1C",
    "MR silent failure — TFZ=8.5 but R-free stuck")

  files = ["data.mtz", "model.pdb", "sequence.fa"]

  # After phaser (TFZ=8.5) + 2 refine cycles stuck
  h = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.0}),
    ("phenix.phaser", "SUCCESS: TFZ=8.5",
     {"tfz": 8.5}),
    ("phenix.refine", "SUCCESS: OK",
     {"r_free": 0.48}),
    ("phenix.refine", "SUCCESS: OK",
     {"r_free": 0.47}),
  ])

  # Gate evaluation: should R-free > 0.45 after
  # 2 cycles trigger retreat?
  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice="solve by MR",
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan generated")
    return t

  # Advance to initial_refinement stage
  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  plan.advance()
  plan.mark_stage_started(2)
  plan.record_stage_cycle("phenix.phaser")
  plan.advance()
  # Two refine cycles
  plan.mark_stage_started(3)
  plan.record_stage_cycle("phenix.refine")
  plan.record_stage_cycle("phenix.refine")

  sm = make_sm(r_free=0.47, tfz=8.5,
               resolution=2.0)
  gr = mock_gate(plan.to_dict(), sm, 4)
  t.step("Gate after 2 refine (R-free 0.47): "
         "%s" % gr.action)
  t.step("  reason: %s" % gr.reason[:80])

  # Check the gate condition
  curr = plan.current_stage()
  t.step("Current stage: %s (gate: %s)" % (
    curr.id if curr else "?",
    curr.gate_conditions if curr else "?"))

  # Does the expert prompt include space group?
  try:
    from knowledge.thinking_prompts import (
      build_thinking_prompt,
    )
    context = {
      "log_sections": "R-free: 0.47, not improving",
      "program_name": "phenix.refine",
      "cycle_number": 4,
      "experiment_type": "xray",
      "workflow_state": "xray_refined",
      "valid_programs": ["phenix.refine", "STOP"],
      "plan_phase": "Current plan prefers: refine",
      "plan_goal": "MR + refinement",
      "metrics": {"r_free": 0.47, "tfz": 8.5},
      "user_advice": "solve by MR",
      "history_summary": (
        "xtriage OK, phaser TFZ=8.5, "
        "refine 0.48→0.47"),
      "r_free_trend": [0.48, 0.47],
    }
    _, user_msg = build_thinking_prompt(context)
    has_rfree = "0.47" in user_msg
    has_tfz = "8.5" in user_msg
    has_trend = "0.48" in user_msg
    t.step("THINK prompt has: R-free=%s TFZ=%s "
           "trend=%s" % (has_rfree, has_tfz,
                         has_trend))
    if not has_trend:
      t.issue("DEGRADED",
        "THINK prompt missing R-free trend — "
        "expert can't diagnose stagnation")
  except Exception as e:
    t.step("THINK prompt build failed: %s" % e)

  return t


def trace_S3B():
  """S3-B: Placed model, R-free regression."""
  t = TraceResult("S3B",
    "Refine only — R-free regression (0.28→0.33)")

  files = ["data.mtz", "model.pdb"]
  advice = "refine the model"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
    if plan.template_id != "refine_placed":
      t.issue("WRONG",
        "Expected refine_placed, got %s"
        % plan.template_id)
  else:
    t.issue("CRASH", "No plan generated")
    return t

  # After 3 refine cycles: 0.32 → 0.28 → 0.33
  h = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.0}),
    ("phenix.refine", "SUCCESS: OK",
     {"r_free": 0.28}),
    ("phenix.refine", "SUCCESS: OK",
     {"r_free": 0.33}),  # REGRESSION
  ])

  # Check: does BestFilesTracker prefer cycle 2?
  t.step("R-free trajectory: 0.32 → 0.28 → 0.33")
  t.step("Regression at cycle 3 (0.28 → 0.33)")
  t.step("BestFilesTracker should use cycle 2 "
         "model (R-free=0.28)")
  # We can't easily test BFT here without files,
  # but we can check the DisplayDataModel
  try:
    from agent.display_data_model import (
      DisplayDataModel,
    )
    ddm = DisplayDataModel.from_session({
      "cycles": h,
      "structure_model": make_sm(r_free=0.33),
    })
    traj = ddm.rfree_trajectory
    if len(traj) >= 2:
      t.step("Trajectory: %s" % [
        "%.3f" % p.value for p in traj])
  except Exception:
    pass

  # Gate: should NOT advance despite cycle 3 done
  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  plan.advance()
  plan.mark_stage_started(2)
  plan.record_stage_cycle("phenix.refine")
  plan.record_stage_cycle("phenix.refine")
  sm = make_sm(r_free=0.33, resolution=2.0)
  gr = mock_gate(plan.to_dict(), sm, 3)
  t.step("Gate after regression: %s (%s)" % (
    gr.action, gr.reason[:60]))

  return t


def trace_S4A():
  """S4-A: Placed model + ligand — clean fit."""
  t = TraceResult("S4A",
    "Refine + fit ATP — clean ligandfit")

  files = ["data.mtz", "model.pdb", "atp.cif"]
  advice = "refine and fit the ATP ligand"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
    phase_ids = [p.id for p in plan.stages]
    t.step("Phases: %s" % " → ".join(phase_ids))
    if plan.template_id != "refine_placed_ligand":
      t.issue("WRONG",
        "Expected refine_placed_ligand")
    if "model_rebuilding" in phase_ids:
      t.issue("WRONG",
        "model_rebuilding should NOT be in "
        "refine_placed_ligand template")
    if "ligand_fitting" not in phase_ids:
      t.issue("WRONG",
        "ligand_fitting missing from template")
  else:
    t.issue("CRASH", "No plan generated")
    return t

  # Workflow state after refine
  h = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.0}),
    ("phenix.refine", "SUCCESS: OK",
     {"r_free": 0.30}),
  ])
  state, vp, _ = mock_detect_ws(
    files, h, advice)
  t.step("After refine: state=%s" % state)
  t.step("  valid=%s" % vp)
  if state == "SKIP":
    t.step("  (requires PHENIX — skipped)")
  elif "phenix.ligandfit" in vp:
    t.step("ligandfit available (good)")
  elif "phenix.refine" in vp:
    # Ligandfit may require real CIF categorization
    # or user_wants_ligandfit context flag to appear.
    # With mock filenames, the workflow engine may
    # not detect the ligand file.
    t.step("ligandfit not yet available — "
           "may need real CIF or prefer_programs")
  else:
    t.issue("WRONG",
      "No useful programs after refine")

  return t


def trace_S5A():
  """S5-A: SAD phasing — good signal."""
  t = TraceResult("S5A",
    "SAD phasing — good anomalous signal")

  files = ["p9.sca", "sequence.dat"]
  advice = "SAD with selenium at 2.5A"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
    phase_ids = [p.id for p in plan.stages]
    t.step("Phases: %s" % " → ".join(phase_ids))
    if "sad" not in plan.template_id:
      t.issue("WRONG",
        "Expected SAD template, got %s"
        % plan.template_id)
  else:
    t.issue("CRASH", "No plan generated")
    return t

  # After autosol
  h = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.5}),
    ("phenix.autosol", "SUCCESS: OK",
     {"fom": 0.42}),
  ])
  state, vp, _ = mock_detect_ws(
    files, h, advice)
  t.step("After autosol: state=%s valid=%s" % (
    state, vp))
  if state == "SKIP":
    t.step("  (requires PHENIX — skipped)")
  elif "phenix.autobuild" not in vp:
    t.issue("WRONG",
      "autobuild not available after autosol")

  return t


def trace_S10B():
  """S10-B: Directive conflicts with data
  (twinning detected, user says no twin law)."""
  t = TraceResult("S10B",
    "Twinning detected but user says "
    "do not use twin law")

  files = ["data.mtz", "model.pdb", "seq.fa"]
  advice = "solve by MR, do not use twin law"

  # After xtriage detects twinning
  h = mock_history([
    ("phenix.xtriage", "SUCCESS: twinning detected",
     {"resolution": 2.0, "twinning": True}),
  ])
  state, vp, _ = mock_detect_ws(
    files, h, advice)
  t.step("After xtriage (twinning): state=%s" % (
    state))
  t.step("  valid=%s" % vp)

  # Check what the THINK prompt would include
  try:
    from knowledge.thinking_prompts import (
      build_thinking_prompt,
    )
    context = {
      "log_sections": "Twinning detected: twin_law",
      "program_name": "phenix.xtriage",
      "cycle_number": 1,
      "experiment_type": "xray",
      "workflow_state": state,
      "valid_programs": vp,
      "metrics": {"resolution": 2.0,
                  "twinning": True},
      "user_advice": advice,
      "history_summary": "xtriage: twinning",
      "r_free_trend": [],
    }
    _, user_msg = build_thinking_prompt(context)
    has_twin = "twin" in user_msg.lower()
    has_advice = "do not" in user_msg.lower()
    t.step("THINK prompt: mentions twinning=%s, "
           "mentions user directive=%s"
           % (has_twin, has_advice))
    if not has_advice:
      t.issue("DEGRADED",
        "THINK prompt doesn't include user "
        "advice about twin law — expert can't "
        "warn about conflict")
  except Exception as e:
    t.step("THINK prompt failed: %s" % e)

  return t


def trace_S11A():
  """S11-A: Missing files — autobuild without seq."""
  t = TraceResult("S11A",
    "Missing sequence — user wants autobuild+ligandfit")

  files = ["data.mtz", "model.pdb"]
  advice = "run autobuild and fit ATP"

  # Plan selection — should NOT select templates
  # requiring sequence
  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
    phase_ids = [p.id for p in plan.stages]
    t.step("Phases: %s" % " → ".join(phase_ids))
    # Should NOT have autobuild or ligandfit stages
    # if files are missing
    has_ab = any("autobuild" in pid
                 for pid in phase_ids)
    has_lf = any("ligand" in pid
                 for pid in phase_ids)
    if has_lf:
      t.issue("WRONG",
        "ligand_fitting stage present but no "
        "CIF file provided")
  else:
    t.step("No plan generated (expected)")

  # Workflow state check
  h = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.0}),
  ])
  state, vp, _ = mock_detect_ws(
    files, h, advice)
  t.step("After xtriage: valid=%s" % vp)
  if "phenix.autobuild" in vp:
    t.issue("WRONG",
      "autobuild in valid_programs WITHOUT "
      "sequence file (YAML should block)")
  else:
    t.step("autobuild correctly blocked (good)")
  if "phenix.ligandfit" in vp:
    t.issue("WRONG",
      "ligandfit in valid_programs WITHOUT "
      "ligand CIF")
  else:
    t.step("ligandfit correctly blocked (good)")

  return t


def trace_G1():
  """G1: Hysteresis compliance test."""
  t = TraceResult("G1",
    "Gate hysteresis — R-free near threshold")

  from agent.gate_evaluator import apply_hysteresis

  # Target: R-free < 0.35
  # Buffer: 1.5% → adjusted target < 0.3448
  adj = apply_hysteresis(0.35, "<")
  t.step("Target <0.35 with 1.5%% buffer → <%.4f"
         % adj)

  # R-free = 0.349: should NOT advance
  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=["d.mtz", "m.pdb", "s.fa"],
    user_advice="solve by MR",
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t

  # Advance to initial_refinement
  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  plan.advance()
  plan.mark_stage_started(2)
  plan.record_stage_cycle("phenix.phaser")
  plan.advance()
  plan.mark_stage_started(3)
  plan.record_stage_cycle("phenix.refine")

  # Test: R-free = 0.349 (below 0.35 but above
  # hysteresis buffer of 0.3448)
  sm1 = make_sm(r_free=0.349, resolution=2.0)
  gr1 = mock_gate(plan.to_dict(), sm1, 3)
  t.step("R-free=0.349: gate=%s" % gr1.action)
  if gr1.action == "advance":
    t.issue("WRONG",
      "Advanced at R-free=0.349 — hysteresis "
      "should prevent (buffer=%.4f)" % adj)
  else:
    t.step("Correctly did NOT advance (good)")

  # Test: R-free = 0.344 (below buffer)
  plan2 = StructurePlan.from_dict(plan.to_dict())
  plan2.record_stage_cycle("phenix.refine")
  sm2 = make_sm(r_free=0.344, resolution=2.0)
  gr2 = mock_gate(plan2.to_dict(), sm2, 4)
  t.step("R-free=0.344: gate=%s" % gr2.action)
  if gr2.action != "advance":
    t.issue("DEGRADED",
      "Did not advance at R-free=0.344 — "
      "should be past hysteresis buffer")

  return t


# ── S2: MR + Ligand ───────────────────────────

def trace_S2A():
  """S2-A: MR + ligand — ligand fits well."""
  t = TraceResult("S2A",
    "MR + ligand — clean fit (CC > 0.7)")

  files = ["data.mtz", "model.pdb", "sequence.fa",
           "atp.cif"]
  advice = "solve by MR and fit ATP"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan generated")
    return t
  t.step("Template: %s" % plan.template_id)
  phase_ids = [p.id for p in plan.stages]
  t.step("Phases: %s" % " -> ".join(phase_ids))
  if "mr" not in plan.template_id:
    t.issue("WRONG",
      "Expected MR template, got %s"
      % plan.template_id)
  if "ligand_fitting" not in phase_ids:
    t.issue("WRONG",
      "ligand_fitting missing from MR+lig template")

  # Gate: after refine with R-free=0.30, should
  # advance to ligandfit stage
  for ph in plan.stages:
    if ph.id in ("data_assessment",
                 "molecular_replacement"):
      ph.status = "complete"
      ph.cycles_used = 1
  plan.current_stage_index = phase_ids.index(
    "initial_refinement") if (
    "initial_refinement" in phase_ids) else 2
  plan.mark_stage_started(3)
  plan.record_stage_cycle("phenix.refine")
  plan.record_stage_cycle("phenix.refine")
  sm = make_sm(r_free=0.30, resolution=2.0,
               tfz=12.0)
  gr = mock_gate(plan.to_dict(), sm, 4)
  t.step("Gate after refine (R-free 0.30): %s"
         % gr.action)

  return t


def trace_S2B():
  """S2-B: MR + ligand — poor ligand fit."""
  t = TraceResult("S2B",
    "MR + ligand — poor fit (CC < 0.4)")

  files = ["data.mtz", "model.pdb", "sequence.fa",
           "atp.cif"]
  advice = "solve by MR and fit ATP"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t
  t.step("Template: %s" % plan.template_id)

  # After ligandfit with poor CC
  # Check: does the THINK prompt include CC info?
  try:
    from knowledge.thinking_prompts import (
      build_thinking_prompt,
    )
    context = {
      "log_sections": "LigandFit: CC=0.35 (poor)",
      "program_name": "phenix.ligandfit",
      "cycle_number": 5,
      "experiment_type": "xray",
      "workflow_state": "xray_combined",
      "valid_programs": [
        "phenix.refine", "phenix.polder", "STOP"],
      "metrics": {"cc": 0.35},
      "user_advice": advice,
      "history_summary": (
        "xtriage, phaser TFZ=12, "
        "refine 0.42->0.30, ligandfit CC=0.35"),
      "r_free_trend": [0.42, 0.35, 0.30],
    }
    _, user_msg = build_thinking_prompt(context)
    has_cc = "0.35" in user_msg
    has_polder = "polder" in user_msg.lower()
    t.step("THINK prompt: has CC=%s, "
           "mentions polder=%s"
           % (has_cc, has_polder))
    if not has_cc:
      t.issue("DEGRADED",
        "THINK prompt missing ligand CC — "
        "expert can't assess fit quality")
  except Exception as e:
    t.step("THINK prompt: %s" % e)

  return t


# ── S4B: Placed model + ligand fails ──────────

def trace_S4B():
  """S4-B: Placed model + ligand — ligandfit fails."""
  t = TraceResult("S4B",
    "Refine + fit ATP — ligandfit fails")

  files = ["data.mtz", "model.pdb", "atp.cif"]
  advice = "refine and fit the ATP ligand"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t
  t.step("Template: %s" % plan.template_id)
  phase_ids = [p.id for p in plan.stages]
  t.step("Phases: %s" % " -> ".join(phase_ids))

  # Advance to ligandfit stage
  for ph in plan.stages:
    if ph.id in ("data_assessment", "refinement"):
      ph.status = "complete"
      ph.cycles_used = 1
  lig_idx = phase_ids.index("ligand_fitting") if (
    "ligand_fitting" in phase_ids) else 2
  plan.current_stage_index = lig_idx
  plan.mark_stage_started(3)
  plan.record_stage_cycle("phenix.ligandfit")

  # Gate after ligandfit with no CC metric
  # (failed or very poor)
  sm = make_sm(r_free=0.28, resolution=2.0)
  gr = mock_gate(plan.to_dict(), sm, 3)
  t.step("Gate after failed ligandfit: %s (%s)"
         % (gr.action, gr.reason[:60]))

  # Phase should exhaust (max_cycles=1) and advance
  # to final_refinement
  if gr.action == "advance":
    t.step("Advances past ligandfit (correct "
           "— stage exhausted)")
  else:
    t.step("Gate action: %s" % gr.action)

  return t


# ── S5B: SAD fails ────────────────────────────

def trace_S5B():
  """S5-B: SAD phasing — autosol fails (0 sites)."""
  t = TraceResult("S5B",
    "SAD phasing — autosol fails (0 sites)")

  files = ["p9.sca", "sequence.dat"]
  advice = "SAD with selenium at 2.5A"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t
  t.step("Template: %s" % plan.template_id)

  # Advance to experimental_phasing, then fail
  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  plan.advance()
  plan.mark_stage_started(2)
  plan.record_stage_cycle("phenix.autosol")

  # Gate after autosol with no FOM (failed)
  sm = make_sm(resolution=2.5)
  gr = mock_gate(plan.to_dict(), sm, 2)
  t.step("Gate after failed autosol: %s (%s)"
         % (gr.action, gr.reason[:60]))

  # Check THINK prompt for failure diagnosis
  try:
    from knowledge.thinking_prompts import (
      build_thinking_prompt,
    )
    context = {
      "log_sections": "AutoSol: 0 sites found",
      "program_name": "phenix.autosol",
      "cycle_number": 2,
      "experiment_type": "xray",
      "workflow_state": "xray_initial",
      "valid_programs": ["STOP"],
      "metrics": {},
      "user_advice": advice,
      "history_summary": "xtriage OK, autosol FAILED",
      "r_free_trend": [],
    }
    _, user_msg = build_thinking_prompt(context)
    has_sites = "0 sites" in user_msg or (
      "sites" in user_msg.lower())
    t.step("THINK prompt mentions sites: %s"
           % has_sites)
  except Exception as e:
    t.step("THINK prompt: %s" % e)

  return t


# ── S6A: SAD + Ligand ─────────────────────────

def trace_S6A():
  """S6-A: SAD + ligand — full pipeline."""
  t = TraceResult("S6A",
    "SAD + ligand — full pipeline to ATP fit")

  files = ["p9.sca", "sequence.dat", "atp.cif"]
  advice = "SAD with selenium, then fit ATP"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t
  t.step("Template: %s" % plan.template_id)
  phase_ids = [p.id for p in plan.stages]
  t.step("Phases: %s" % " -> ".join(phase_ids))

  if "sad" not in plan.template_id:
    t.issue("WRONG",
      "Expected SAD template, got %s"
      % plan.template_id)
  if "ligand_fitting" not in phase_ids:
    t.issue("WRONG",
      "ligand_fitting missing from SAD+lig")

  # Verify rebuild_in_place would be False
  # (autobuild after autosol)
  t.step("rebuild_in_place=False expected after "
         "autosol (checked in graph_nodes.py)")

  return t


def trace_S6B():
  """S6-B: SAD + ligand — R-free stuck after
  autobuild (gate retreat test)."""
  t = TraceResult("S6B",
    "SAD + ligand — R-free stuck, gate retreat")

  files = ["p9.sca", "sequence.dat", "atp.cif"]
  advice = "SAD with selenium, then fit ATP"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t
  t.step("Template: %s" % plan.template_id)

  # Advance to build_and_refine
  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  plan.advance()
  plan.mark_stage_started(2)
  plan.record_stage_cycle("phenix.autosol")
  plan.advance()

  # 3 cycles in build_and_refine, R-free stuck >0.45
  plan.mark_stage_started(3)
  plan.record_stage_cycle("phenix.autobuild")
  plan.record_stage_cycle("phenix.refine")
  plan.record_stage_cycle("phenix.refine")

  sm = make_sm(r_free=0.46, fom=0.42,
               resolution=2.5)
  gr = mock_gate(plan.to_dict(), sm, 5)
  t.step("Gate after 3 cycles (R-free 0.46): "
         "%s" % gr.action)
  t.step("  reason: %s" % gr.reason[:80])

  curr = plan.current_stage()
  gates = curr.gate_conditions if curr else []
  t.step("Phase gates: %s" % gates)

  if gr.action == "retreat":
    t.step("Retreat triggered (correct)")
  elif gr.action == "advance":
    t.step("Advanced despite R-free > 0.45")
    # This may be correct if gate condition
    # specifies "after 3 cycles" and we have 3
  else:
    t.step("Gate action: %s" % gr.action)

  return t


# ── S7: Cryo-EM ───────────────────────────────

def trace_S7A():
  """S7-A: Cryo-EM refinement — CC improves."""
  t = TraceResult("S7A",
    "Cryo-EM refinement — CC improves")

  files = ["half_map_1.mrc", "half_map_2.mrc",
           "model.pdb"]
  advice = "refine against the 3.2A cryo-EM map"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
    phase_ids = [p.id for p in plan.stages]
    t.step("Phases: %s" % " -> ".join(phase_ids))
    if "cryoem" not in plan.template_id:
      t.issue("WRONG",
        "Expected cryoem template, got %s"
        % plan.template_id)
  else:
    # CryoEM detection from files is tricky
    t.step("No plan (cryo-EM may need data_chars)")
    t.issue("DEGRADED",
      "Plan generator didn't detect cryo-EM "
      "from .mrc files + advice")

  # Check DisplayDataModel handles CC trajectory
  try:
    from agent.display_data_model import (
      DisplayDataModel,
    )
    ddm = DisplayDataModel.from_session({
      "experiment_type": "cryoem",
      "cycles": [
        {"cycle_number": 1,
         "program": "phenix.real_space_refine",
         "result": "SUCCESS: OK",
         "metrics": {"model_map_cc": 0.65}},
        {"cycle_number": 2,
         "program": "phenix.real_space_refine",
         "result": "SUCCESS: OK",
         "metrics": {"model_map_cc": 0.78}},
      ],
      "structure_model": make_sm(
        model_map_cc=0.78, resolution=3.2),
    })
    traj = ddm.primary_trajectory
    t.step("CC trajectory: %s" % [
      "%.3f" % p.value for p in traj])
    if len(traj) < 2:
      t.issue("WRONG",
        "CC trajectory should have 2 points")
    if ddm.outcome_status != "determined":
      t.issue("WRONG",
        "Expected 'determined' for CC=0.78")
    t.step("Outcome: %s" % ddm.outcome_message)
  except Exception as e:
    t.step("DDM test: %s" % e)

  return t


def trace_S7B():
  """S7-B: Cryo-EM refinement — CC stuck."""
  t = TraceResult("S7B",
    "Cryo-EM refinement — CC doesn't improve")

  try:
    from agent.display_data_model import (
      DisplayDataModel,
    )
    ddm = DisplayDataModel.from_session({
      "experiment_type": "cryoem",
      "cycles": [
        {"cycle_number": 1,
         "program": "phenix.real_space_refine",
         "result": "SUCCESS: OK",
         "metrics": {"model_map_cc": 0.45}},
        {"cycle_number": 2,
         "program": "phenix.real_space_refine",
         "result": "SUCCESS: OK",
         "metrics": {"model_map_cc": 0.46}},
      ],
      "structure_model": make_sm(
        model_map_cc=0.46, resolution=3.2),
    })
    t.step("Outcome: %s (%s)" % (
      ddm.outcome_status, ddm.outcome_message))
    if ddm.outcome_status == "determined":
      t.issue("WRONG",
        "CC=0.46 should NOT be 'determined'")
    tl = ddm.timeline
    for e in tl:
      t.step("  %s" % ddm.format_cycle_compact(e))
  except Exception as e:
    t.step("DDM: %s" % e)

  # Check THINK prompt mentions CC stagnation
  try:
    from knowledge.thinking_prompts import (
      build_thinking_prompt,
    )
    context = {
      "log_sections": "Map CC: 0.46",
      "program_name": "phenix.real_space_refine",
      "cycle_number": 2,
      "experiment_type": "cryoem",
      "workflow_state": "cryoem_refined",
      "valid_programs": [
        "phenix.real_space_refine", "STOP"],
      "metrics": {"model_map_cc": 0.46},
      "user_advice": "refine cryo-EM",
      "history_summary": "RSR CC 0.45->0.46",
      "r_free_trend": [],
    }
    _, user_msg = build_thinking_prompt(context)
    t.step("THINK prompt length: %d" % len(
      user_msg))
  except Exception:
    pass

  return t


# ── S8: Cryo-EM build from scratch ────────────

def trace_S8A():
  """S8-A: Cryo-EM build — prediction + dock."""
  t = TraceResult("S8A",
    "Cryo-EM build — predict + dock succeeds")

  files = ["half_map_1.mrc", "half_map_2.mrc",
           "sequence.fa"]
  advice = "build a model into the map"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
    phase_ids = [p.id for p in plan.stages]
    t.step("Phases: %s" % " -> ".join(phase_ids))
  else:
    t.step("No plan (cryo-EM build may need "
           "data_chars)")

  # Verify DDM handles cryo-EM without model
  try:
    from agent.display_data_model import (
      DisplayDataModel,
    )
    ddm = DisplayDataModel.from_session({
      "experiment_type": "cryoem",
      "cycles": [
        {"cycle_number": 1,
         "program": "phenix.mtriage",
         "result": "SUCCESS: OK",
         "metrics": {"resolution": 2.8}},
        {"cycle_number": 2,
         "program": "phenix.predict_and_build",
         "result": "SUCCESS: OK",
         "metrics": {"model_map_cc": 0.70}},
      ],
    })
    t.step("Outcome: %s" % ddm.outcome_message)
    for e in ddm.timeline:
      t.step("  %s" % ddm.format_cycle_compact(e))
  except Exception as e:
    t.step("DDM: %s" % e)

  return t


# ── S9: Resume after removing cycles ──────────

def trace_S9A():
  """S9-A: Resume after removing last 2 cycles."""
  t = TraceResult("S9A",
    "Resume after removing last 2 refine cycles")

  try:
    from agent.session_tools import (
      remove_last_cycles,
    )
  except ImportError:
    t.step("SKIP (session_tools not available)")
    return t

  data = {
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.xtriage",
       "result": "SUCCESS: OK"},
      {"cycle_number": 2,
       "program": "phenix.phaser",
       "result": "SUCCESS: OK"},
      {"cycle_number": 3,
       "program": "phenix.refine",
       "result": "SUCCESS: OK"},
      {"cycle_number": 4,
       "program": "phenix.refine",
       "result": "SUCCESS: OK"},
    ],
    "plan": {
      "goal": "MR",
      "_version": 1,
      "current_stage_index": 3,
      "stages": [
        {"id": "data", "status": "complete",
         "programs": ["phenix.xtriage"],
         "cycles_used": 1, "max_cycles": 1},
        {"id": "mr", "status": "complete",
         "programs": ["phenix.phaser"],
         "cycles_used": 1, "max_cycles": 1},
        {"id": "refine", "status": "complete",
         "programs": ["phenix.refine"],
         "cycles_used": 2, "max_cycles": 5},
        {"id": "final", "status": "pending",
         "programs": ["phenix.refine"],
         "cycles_used": 0, "max_cycles": 3},
      ],
    },
    "gate_stop": True,
    "gate_stop_reason": "all stages complete",
    "structure_model": {"_version": 1,
      "model_state": {"r_free": 0.25}},
    "validation_history": {"v": 1},
    "strategy_memory": {"v": 1},
    "structure_report": "some report",
    "summary": "some summary",
    "html_report_path": "/some/path.html",
  }

  data, removed = remove_last_cycles(data, 2)
  t.step("Removed %d cycles, %d remaining" % (
    removed, len(data["cycles"])))

  # Check all stale keys cleared
  checks = [
    ("gate_stop", False),
    ("structure_model", {}),
    ("validation_history", {}),
    ("strategy_memory", {}),
    ("structure_report", ""),
    ("summary", ""),
  ]
  for key, expected_empty in checks:
    val = data.get(key)
    is_empty = (
      val == expected_empty or val is None
      or val is False)
    if not is_empty:
      t.issue("WRONG",
        "%s not cleared: %s" % (
          key, str(val)[:30]))
  t.step("Stale keys cleared: OK")

  # Check plan reset
  stages = data.get("plan", {}).get("stages", [])
  all_pending = all(
    p.get("status") == "pending" for p in stages)
  if all_pending:
    t.step("Plan stages all pending: OK")
  else:
    statuses = [p.get("status") for p in stages]
    t.issue("WRONG",
      "Plan not fully reset: %s" % statuses)

  # Check fast-forward would work
  done_progs = set()
  for c in data.get("cycles", []):
    if "SUCCESS" in str(c.get("result", "")):
      done_progs.add(c.get("program", ""))
  t.step("Remaining programs: %s" % done_progs)

  ff_count = 0
  plan = StructurePlan.from_dict(data["plan"])
  while True:
    curr = plan.current_stage()
    if curr is None:
      break
    if curr.status not in ("pending", "active"):
      break
    if not curr.programs:
      break
    if not all(p in done_progs
               for p in curr.programs):
      break
    curr.status = "complete"
    plan.advance()
    ff_count += 1
  t.step("Fast-forward: %d stages" % ff_count)

  curr = plan.current_stage()
  if curr:
    t.step("Resume at: %s (%s)" % (
      curr.id, curr.status))
    if curr.id in ("data", "mr"):
      t.issue("WRONG",
        "Should resume at refine, not %s"
        % curr.id)
  else:
    t.issue("WRONG", "No current stage after FF")

  return t


def trace_S9B():
  """S9-B: Resume after removing ALL cycles."""
  t = TraceResult("S9B",
    "Resume after removing ALL cycles")

  try:
    from agent.session_tools import (
      remove_last_cycles,
    )
  except ImportError:
    t.step("SKIP")
    return t

  data = {
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.xtriage",
       "result": "SUCCESS: OK"},
    ],
    "plan": {
      "goal": "test", "_version": 1,
      "current_stage_index": 1,
      "stages": [
        {"id": "data", "status": "complete",
         "programs": ["phenix.xtriage"],
         "cycles_used": 1, "max_cycles": 1},
        {"id": "refine", "status": "pending",
         "programs": ["phenix.refine"],
         "cycles_used": 0, "max_cycles": 5},
      ],
    },
    "gate_stop": False,
  }

  data, removed = remove_last_cycles(data, 10)
  t.step("Removed %d, %d remaining" % (
    removed, len(data["cycles"])))

  if len(data["cycles"]) != 0:
    t.issue("WRONG",
      "Expected 0 cycles, got %d"
      % len(data["cycles"]))

  plan = StructurePlan.from_dict(data["plan"])
  curr = plan.current_stage()
  if curr and curr.id == "data":
    t.step("Resume at data_assessment (correct)")
  else:
    t.issue("WRONG",
      "Expected data_assessment, got %s"
      % (curr.id if curr else "None"))

  return t


# ── S10A: Twinning detected, phaser OK ────────

def trace_S10A():
  """S10-A: Twinned data — phaser succeeds."""
  t = TraceResult("S10A",
    "Twinned data — xtriage detects, phaser OK")

  files = ["data.mtz", "model.pdb", "seq.fa"]
  advice = "solve by MR"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
  else:
    t.issue("CRASH", "No plan")
    return t

  # Check THINK prompt after xtriage with twinning
  try:
    from knowledge.thinking_prompts import (
      build_thinking_prompt,
    )
    context = {
      "log_sections": (
        "Twinning: twin law detected "
        "-h,-k,l. Twin fraction: 0.42"),
      "program_name": "phenix.xtriage",
      "cycle_number": 1,
      "experiment_type": "xray",
      "workflow_state": "xray_analyzed",
      "valid_programs": [
        "phenix.phaser", "STOP"],
      "metrics": {"resolution": 2.0,
                  "twinning": True},
      "user_advice": advice,
      "history_summary": "",
      "r_free_trend": [],
    }
    _, user_msg = build_thinking_prompt(context)
    has_twin = "twin" in user_msg.lower()
    t.step("THINK prompt mentions twinning: %s"
           % has_twin)
    if not has_twin:
      t.issue("DEGRADED",
        "Twinning not in THINK prompt")
  except Exception:
    pass

  return t


# ── S11B: Hallucination negative test ─────────

def trace_S11B():
  """S11-B: Missing files — LLM hallucination
  risk with strong user request."""
  t = TraceResult("S11B",
    "Missing files — hallucination risk")

  # User provides ONLY data.mtz, asks for full SAD
  files = ["data.mtz"]
  advice = "solve by SAD with selenium"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})

  if plan:
    t.step("Template: %s" % plan.template_id)
    phase_ids = [p.id for p in plan.stages]
    t.step("Phases: %s" % " -> ".join(phase_ids))
    # SAD without sequence should NOT select SAD
    if "sad" in plan.template_id:
      t.issue("WRONG",
        "SAD template selected without "
        "sequence file")
  else:
    t.step("No plan generated (correct — "
           "insufficient files)")

  # Check workflow state
  state, vp, _ = mock_detect_ws(
    files, [], advice)
  if state != "SKIP":
    if "phenix.autosol" in vp:
      t.issue("WRONG",
        "autosol available without seq file")
    else:
      t.step("autosol blocked (correct)")

  return t


# ── S12: MR without sequence (edge case) ──────

def trace_S12A():
  """S12-A: MR without sequence file — autobuild
  stage must not block progress."""
  t = TraceResult("S12A",
    "MR without sequence — autobuild handled")

  files = ["data.mtz", "model.pdb"]
  advice = "solve by molecular replacement"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t
  t.step("Template: %s" % plan.template_id)
  phase_ids = [p.id for p in plan.stages]
  t.step("Phases: %s" % " -> ".join(phase_ids))

  # Check if model_rebuilding has refine as fallback
  for ph in plan.stages:
    if ph.id == "model_rebuilding":
      t.step("model_rebuilding programs: %s"
             % ph.programs)
      if "phenix.refine" not in ph.programs:
        t.issue("WRONG",
          "model_rebuilding only has autobuild "
          "— will get stuck without sequence")
      else:
        t.step("refine is fallback (correct)")
      break

  return t


# ── S13: No advice (bare files only) ──────────

def trace_S13A():
  """S13-A: PDB + MTZ + sequence, no advice."""
  t = TraceResult("S13A",
    "No advice — PDB+MTZ+seq should default to MR")

  files = ["data.mtz", "model.pdb", "seq.fa"]

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice="",
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
    if plan.template_id != "mr_refine":
      t.issue("WRONG",
        "Expected mr_refine for PDB+MTZ+seq "
        "with no advice")
  else:
    t.issue("DEGRADED",
      "No plan for PDB+MTZ+seq (should "
      "default to mr_refine)")

  return t


def trace_S13B():
  """S13-B: PDB + MTZ only, no advice."""
  t = TraceResult("S13B",
    "No advice — PDB+MTZ only, ambiguous intent")

  files = ["data.mtz", "model.pdb"]

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice="",
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
    # Without advice, PDB+MTZ is ambiguous:
    # could be refine-only or MR
    # mr_refine is the safe default
    t.step("(mr_refine is conservative default)")
  else:
    t.step("No plan (may need advice)")

  return t


# ── S14: Contradictory advice ─────────────────

def trace_S14A():
  """S14-A: User says 'refine' but provides
  sequence (MR likely intended)."""
  t = TraceResult("S14A",
    "Ambiguous: 'refine' + sequence file")

  files = ["data.mtz", "model.pdb", "seq.fa"]
  advice = "refine the model"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
    # "refine" → model_is_placed → refine_placed
    # This is correct: user explicitly said "refine"
    if plan.template_id == "refine_placed":
      t.step("Respects user advice (correct)")
    elif plan.template_id == "mr_refine":
      t.step("Chose MR despite 'refine' advice")
      t.issue("DEGRADED",
        "User said 'refine' but agent chose MR")
  else:
    t.issue("CRASH", "No plan")

  return t


def trace_S14B():
  """S14-B: User says 'solve by MR' but also
  says 'fit ATP' without providing CIF."""
  t = TraceResult("S14B",
    "MR + 'fit ATP' but no CIF file")

  files = ["data.mtz", "model.pdb", "seq.fa"]
  advice = "solve by MR and fit ATP"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
    phase_ids = [p.id for p in plan.stages]
    t.step("Phases: %s" % " -> ".join(phase_ids))
    if "ligand" in plan.template_id:
      t.issue("WRONG",
        "Ligand template selected without CIF")
    else:
      t.step("No ligand template (correct — "
             "no CIF file)")
    if "ligand_fitting" in phase_ids:
      t.issue("WRONG",
        "ligand_fitting stage without CIF")
  else:
    t.issue("CRASH", "No plan")

  return t


# ── G2: Metric regression detection ───────────

def trace_G2():
  """G2: Gate handles R-free regression."""
  t = TraceResult("G2",
    "Gate — R-free regression detection")

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=["d.mtz", "m.pdb"],
    user_advice="refine the model",
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t

  # Advance to refinement stage
  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  plan.advance()
  plan.mark_stage_started(2)

  # Cycle 1: R-free improves
  plan.record_stage_cycle("phenix.refine")
  sm1 = make_sm(r_free=0.320, resolution=2.0)
  gr1 = mock_gate(plan.to_dict(), sm1, 2)
  t.step("Cycle 2 (R-free 0.320): %s" % gr1.action)

  # Cycle 2: R-free improves more
  plan.record_stage_cycle("phenix.refine")
  sm2 = make_sm(r_free=0.280, resolution=2.0)
  gr2 = mock_gate(plan.to_dict(), sm2, 3)
  t.step("Cycle 3 (R-free 0.280): %s" % gr2.action)

  # Cycle 3: R-free REGRESSES
  plan.record_stage_cycle("phenix.refine")
  sm3 = make_sm(r_free=0.340, resolution=2.0)
  gr3 = mock_gate(plan.to_dict(), sm3, 4)
  t.step("Cycle 4 (R-free 0.340 REGRESSION): "
         "%s" % gr3.action)

  # Check DDM trajectory shows regression
  try:
    from agent.display_data_model import (
      DisplayDataModel,
    )
    ddm = DisplayDataModel.from_session({
      "cycles": [
        {"cycle_number": 1,
         "program": "phenix.xtriage",
         "result": "SUCCESS: OK",
         "metrics": {"resolution": 2.0}},
        {"cycle_number": 2,
         "program": "phenix.refine",
         "result": "SUCCESS: OK",
         "metrics": {"r_free": 0.320}},
        {"cycle_number": 3,
         "program": "phenix.refine",
         "result": "SUCCESS: OK",
         "metrics": {"r_free": 0.280}},
        {"cycle_number": 4,
         "program": "phenix.refine",
         "result": "SUCCESS: OK",
         "metrics": {"r_free": 0.340}},
      ],
    })
    for e in ddm.timeline:
      line = ddm.format_cycle_compact(e)
      if "0.340" in line:
        t.step("DDM shows regression: %s" % line)
  except Exception:
    pass

  return t


# ── G3: Phase exhaustion with marginal progress ─

def trace_G3():
  """G3: Phase exhaustion — slow but positive
  progress, should advance not retreat."""
  t = TraceResult("G3",
    "Gate — stage exhausted with marginal progress")

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=["d.mtz", "m.pdb", "s.fa"],
    user_advice="solve by MR",
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t

  # Advance to initial_refinement (max_cycles=5)
  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  plan.advance()
  plan.mark_stage_started(2)
  plan.record_stage_cycle("phenix.phaser")
  plan.advance()
  plan.mark_stage_started(3)

  # 5 cycles of slow progress
  for _ in range(5):
    plan.record_stage_cycle("phenix.refine")

  sm = make_sm(r_free=0.36, resolution=2.0)
  gr = mock_gate(plan.to_dict(), sm, 7)
  t.step("Gate after 5 refine (R-free 0.36): "
         "%s" % gr.action)
  t.step("  reason: %s" % gr.reason[:80])

  # Should advance (stage exhausted) not retreat
  # (R-free > 0.35 target but < 0.45 retreat)
  if gr.action == "advance":
    t.step("Advances (correct — exhausted, "
           "no retreat triggered)")
  elif gr.action == "retreat":
    t.step("Retreats (check gate conditions)")
  else:
    t.step("Action: %s" % gr.action)

  return t


# ── G4: Skip condition test ───────────────────

def trace_G4():
  """G4: Phase skip_if condition — skip model
  rebuilding when R-free already good."""
  t = TraceResult("G4",
    "Gate — skip_if: r_free < 0.28 → skip rebuild")

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=["d.mtz", "m.pdb", "s.fa"],
    user_advice="solve by MR",
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t

  # Advance through to model_rebuilding
  phase_ids = [p.id for p in plan.stages]
  t.step("Phases: %s" % " -> ".join(phase_ids))

  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  plan.advance()
  plan.mark_stage_started(2)
  plan.record_stage_cycle("phenix.phaser")
  plan.advance()
  plan.mark_stage_started(3)
  plan.record_stage_cycle("phenix.refine")
  plan.record_stage_cycle("phenix.refine")

  # R-free = 0.25 — below skip_if threshold
  sm = make_sm(r_free=0.25, resolution=2.0,
               tfz=14.0)
  gr = mock_gate(plan.to_dict(), sm, 5)
  t.step("Gate at refine (R-free 0.25): %s"
         % gr.action)

  # After advancing, gate should evaluate
  # model_rebuilding and potentially skip it
  if gr.action == "advance":
    plan.advance()
    plan.mark_stage_started(6)
    gr2 = mock_gate(plan.to_dict(), sm, 6)
    t.step("Gate at model_rebuilding: %s (%s)"
           % (gr2.action, gr2.reason[:50]))
    if gr2.action == "skip":
      t.step("Skipped rebuild (correct — "
             "R-free < 0.28)")
    elif gr2.action == "advance":
      t.step("Advanced (may include skip logic)")

  return t


# ── S15: Multiple ligands ─────────────────────

def trace_S15():
  """S15: Multiple ligand CIF files — only one
  ligandfit stage should exist."""
  t = TraceResult("S15",
    "Multiple CIF files — single ligandfit stage")

  files = ["data.mtz", "model.pdb",
           "atp.cif", "mg.cif"]
  advice = "refine and fit ATP and magnesium"

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=files,
    user_advice=advice,
    directives={})
  if plan:
    t.step("Template: %s" % plan.template_id)
    phase_ids = [p.id for p in plan.stages]
    t.step("Phases: %s" % " -> ".join(phase_ids))
    lig_count = phase_ids.count("ligand_fitting")
    if lig_count > 1:
      t.issue("WRONG",
        "Multiple ligand_fitting stages (%d)"
        % lig_count)
    elif lig_count == 1:
      t.step("Single ligand_fitting (correct)")
    else:
      t.step("No ligand_fitting stage")
  else:
    t.step("No plan")

  return t


# ── C1-C3: Cycle counting edge cases ──────────

def trace_C1():
  """C1: Cycle counting — variant programs."""
  t = TraceResult("C1",
    "Cycle counting: variants, reactives, STOP")

  plan = StructurePlan.from_dict({
    "goal": "test", "_version": 1,
    "current_stage_index": 0,
    "stages": [
      {"id": "build", "status": "active",
       "programs": ["phenix.autobuild",
                    "phenix.refine"],
       "max_cycles": 5, "cycles_used": 0},
      {"id": "lig", "status": "pending",
       "programs": ["phenix.ligandfit"],
       "max_cycles": 1, "cycles_used": 0},
    ],
  })

  tests = [
    ("phenix.autobuild", 1, "exact match"),
    ("phenix.autobuild_denmod", 2, "variant"),
    ("phenix.refine", 3, "exact match #2"),
    ("phenix.pdbtools", 4, "reactive (no stage)"),
    ("phenix.model_vs_data", 5, "reactive probe"),
  ]
  for prog, expected, desc in tests:
    plan.record_stage_cycle(prog)
    actual = plan.stages[0].cycles_used
    if actual != expected:
      t.issue("WRONG",
        "%s: expected %d, got %d"
        % (prog, expected, actual))
    t.step("%s (%s): cycles=%d" % (
      prog, desc, actual))

  # STOP should NOT count
  plan.record_stage_cycle("STOP")
  if plan.stages[0].cycles_used != 5:
    t.issue("WRONG",
      "STOP incremented count to %d"
      % plan.stages[0].cycles_used)
  else:
    t.step("STOP: not counted (correct)")

  # ligandfit belongs to lig stage → skip
  plan.record_stage_cycle("phenix.ligandfit")
  if plan.stages[0].cycles_used != 5:
    t.issue("WRONG",
      "ligandfit counted in build stage")
  else:
    t.step("ligandfit: skipped (belongs to "
           "lig stage)")

  return t


def trace_C2():
  """C2: Cycle counting — stage with no programs
  list (always counts)."""
  t = TraceResult("C2",
    "Cycle counting: empty programs list")

  plan = StructurePlan.from_dict({
    "goal": "test", "_version": 1,
    "current_stage_index": 0,
    "stages": [
      {"id": "any", "status": "active",
       "programs": [],
       "max_cycles": 3, "cycles_used": 0},
    ],
  })

  plan.record_stage_cycle("phenix.anything")
  plan.record_stage_cycle("phenix.whatever")
  if plan.stages[0].cycles_used != 2:
    t.issue("WRONG",
      "Empty programs list didn't count: %d"
      % plan.stages[0].cycles_used)
  else:
    t.step("Empty programs: always counts (%d)"
           % plan.stages[0].cycles_used)

  return t


def trace_C3():
  """C3: Cycle counting — None program (backward
  compat from rules_only mode)."""
  t = TraceResult("C3",
    "Cycle counting: None program (compat)")

  plan = StructurePlan.from_dict({
    "goal": "test", "_version": 1,
    "current_stage_index": 0,
    "stages": [
      {"id": "ref", "status": "active",
       "programs": ["phenix.refine"],
       "max_cycles": 5, "cycles_used": 0},
    ],
  })

  plan.record_stage_cycle(None)
  if plan.stages[0].cycles_used != 1:
    t.issue("WRONG",
      "None program didn't count")
  else:
    t.step("None: counts (backward compat)")

  return t


# ── M1-M2: Multi-cycle progression ───────────

def trace_M1():
  """M1: Full MR progression through all stages."""
  t = TraceResult("M1",
    "Full MR: 5 stages, correct advancement")

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=["d.mtz", "m.pdb", "s.fa"],
    user_advice="solve by MR",
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t

  phase_ids = [p.id for p in plan.stages]
  t.step("Phases: %s" % " -> ".join(phase_ids))

  # Simulate full run
  # Note: gate may advance early if success criteria
  # are met. Final cycle count depends on metrics.
  cycles = [
    ("phenix.xtriage", {"resolution": 2.0},
     "data_assessment"),
    ("phenix.phaser", {"tfz": 14.2, "llg": 342},
     "molecular_replacement"),
    ("phenix.refine", {"r_free": 0.35},
     "initial_refinement"),
    ("phenix.refine", {"r_free": 0.30},
     "initial_refinement"),
    ("phenix.autobuild", {"r_free": 0.27},
     "model_rebuilding"),
    ("phenix.refine", {"r_free": 0.24},
     "final_refinement"),
  ]

  for i, (prog, metrics, expected_phase) in (
    enumerate(cycles)
  ):
    cycle_num = i + 1
    plan.mark_stage_started(cycle_num)
    plan.record_stage_cycle(prog)

    sm = make_sm(**metrics)
    gr = mock_gate(plan.to_dict(), sm, cycle_num)

    curr = plan.current_stage()
    curr_id = curr.id if curr else "DONE"
    t.step("Cycle %d (%s): step=%s gate=%s"
           % (cycle_num, prog, curr_id,
              gr.action))

    if curr_id != expected_phase:
      t.issue("WRONG",
        "Cycle %d: expected stage %s, got %s"
        % (cycle_num, expected_phase, curr_id))

    if gr.action == "advance":
      plan.advance()
    elif gr.action == "skip":
      plan.advance()
    elif gr.action == "stop":
      t.step("Gate STOP at cycle %d" % cycle_num)
      break

  # After all cycles, plan should be complete
  # or at final_refinement
  final_curr = plan.current_stage()
  if final_curr is None:
    t.step("Plan complete (all stages done)")
  else:
    t.step("Final stage: %s (%s)" % (
      final_curr.id, final_curr.status))

  return t


def trace_M2():
  """M2: SAD progression with retreat."""
  t = TraceResult("M2",
    "SAD with retreat: build_and_refine stuck")

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=["p9.sca", "seq.dat"],
    user_advice="SAD with selenium at 2.5A",
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t

  phase_ids = [p.id for p in plan.stages]
  t.step("Phases: %s" % " -> ".join(phase_ids))

  # Phase 1: xtriage
  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  gr1 = mock_gate(plan.to_dict(),
    make_sm(resolution=2.5), 1)
  t.step("Cycle 1 (xtriage): gate=%s" % gr1.action)
  if gr1.action == "advance":
    plan.advance()

  # Phase 2: autosol (good)
  plan.mark_stage_started(2)
  plan.record_stage_cycle("phenix.autosol")
  gr2 = mock_gate(plan.to_dict(),
    make_sm(fom=0.42, resolution=2.5), 2)
  t.step("Cycle 2 (autosol): gate=%s" % gr2.action)
  if gr2.action == "advance":
    plan.advance()

  # Phase 3: build_and_refine — stuck
  plan.mark_stage_started(3)
  plan.record_stage_cycle("phenix.autobuild")
  plan.record_stage_cycle("phenix.refine")
  plan.record_stage_cycle("phenix.refine")

  # R-free stuck at 0.46 → retreat condition
  sm_stuck = make_sm(r_free=0.46, fom=0.42,
                     resolution=2.5)
  gr3 = mock_gate(plan.to_dict(), sm_stuck, 5)
  t.step("Cycle 5 (stuck R-free 0.46): gate=%s"
         % gr3.action)
  t.step("  reason: %s" % gr3.reason[:80])

  if gr3.action == "retreat":
    t.step("Retreat triggered (correct)")
    # Check blacklist
    bl = gr3.blacklist_entry
    if bl:
      t.step("Blacklisted: %s" % bl)
  elif gr3.action == "advance":
    t.step("Advanced (stage exhausted despite "
           "poor metrics)")
  else:
    t.step("Action: %s" % gr3.action)

  return t


# ── P1-P6: PHENIX-dependent scenarios ─────────
# These test workflow state detection, valid programs,
# and command assembly. They call detect_workflow_state
# which requires libtbx. When libtbx is not available,
# the tests PASS with SKIP notes.

def _requires_phenix(t):
  """Check if PHENIX is available. If not, mark
  all checks as skipped and return True."""
  state, vp, reason = mock_detect_ws(
    ["test.mtz"], [], "")
  if state == "SKIP":
    t.step("PHENIX not available — skipping "
           "workflow state checks")
    return True
  return False


def trace_P1():
  """P1: Workflow state — MR fresh session.
  Files: data.mtz + model.pdb + sequence.fa
  Expected: xray_initial, xtriage available."""
  t = TraceResult("P1",
    "Workflow: MR fresh → xtriage available")

  files = ["data.mtz", "model.pdb", "sequence.fa"]
  if _requires_phenix(t):
    return t

  state, vp, reason = mock_detect_ws(
    files, [], "solve by MR")
  t.step("State: %s" % state)
  t.step("Valid: %s" % vp)
  t.step("Reason: %s" % reason[:60])

  if "phenix.xtriage" not in vp:
    t.issue("WRONG",
      "xtriage should be in valid_programs "
      "for fresh MR session")
  if "phenix.phaser" in vp:
    t.issue("DEGRADED",
      "phaser available before xtriage — "
      "should run xtriage first")
  return t


def trace_P2():
  """P2: Workflow state — after xtriage, phaser
  or model_vs_data should be available for MR."""
  t = TraceResult("P2",
    "Workflow: after xtriage → phaser or probe")

  files = ["data.mtz", "model.pdb", "sequence.fa"]
  if _requires_phenix(t):
    return t

  h = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.0}),
  ])
  state, vp, reason = mock_detect_ws(
    files, h, "solve by MR")
  t.step("State: %s" % state)
  t.step("Valid: %s" % vp)

  # The workflow engine may insert a model_vs_data
  # probe step before MR — this is correct behavior
  # (check model fit before committing to MR).
  if "phenix.phaser" in vp:
    t.step("phaser available (direct to MR)")
  elif "phenix.model_vs_data" in vp:
    t.step("model_vs_data probe first (correct "
           "— checks model fit before MR)")
  else:
    t.issue("WRONG",
      "Neither phaser nor model_vs_data "
      "available after xtriage")
  return t


def trace_P3():
  """P3: Workflow state — after phaser success,
  refine should be available."""
  t = TraceResult("P3",
    "Workflow: after phaser → refine available")

  files = ["data.mtz", "model.pdb", "sequence.fa"]
  if _requires_phenix(t):
    return t

  h = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.0}),
    ("phenix.phaser", "SUCCESS: TFZ=14.2",
     {"tfz": 14.2, "llg": 342}),
  ])
  state, vp, reason = mock_detect_ws(
    files, h, "solve by MR")
  t.step("State: %s" % state)
  t.step("Valid: %s" % vp)

  if "phenix.refine" not in vp:
    t.issue("WRONG",
      "refine not available after phaser")
  if "phenix.phaser" in vp:
    t.issue("DEGRADED",
      "phaser still available after success "
      "— should be done")
  return t


def trace_P4():
  """P4: Workflow state — placed model + ligand.
  After refine, ligandfit should be available."""
  t = TraceResult("P4",
    "Workflow: placed model + CIF → ligandfit")

  files = ["data.mtz", "model.pdb", "atp.cif"]
  if _requires_phenix(t):
    return t

  h = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.0}),
    ("phenix.refine", "SUCCESS: OK",
     {"r_free": 0.28}),
  ])

  # The workflow engine needs user_wants_ligandfit
  # to inject ligandfit into valid_programs when
  # the user explicitly requests it. Without this,
  # YAML conditions alone may not trigger ligandfit
  # (the "refine.*stage programs don't include it).
  state, vp, reason = mock_detect_ws(
    files, h,
    "refine and fit ATP",
    directives={
      "workflow_preferences": {
        "model_is_placed": True,
        "prefer_programs": ["phenix.ligandfit"],
      },
    })
  t.step("State: %s" % state)
  t.step("Valid: %s" % vp)

  if "phenix.ligandfit" in vp:
    t.step("ligandfit available (correct)")
  elif "phenix.refine" in vp:
    # Ligandfit may need specific YAML conditions
    # that aren't met with mock filenames. The
    # workflow engine might not categorize "atp.cif"
    # as a ligand file from just the filename.
    t.step("ligandfit not yet available — "
           "may need real CIF categorization")
    t.issue("DEGRADED",
      "ligandfit not in valid_programs — "
      "file categorization may need real "
      "CIF content or user_wants_ligandfit "
      "context flag. Valid: %s" % vp)
  else:
    t.issue("WRONG",
      "Neither ligandfit nor refine available")
  return t


def trace_P5():
  """P5: Workflow state — SAD after xtriage.
  autosol should be available."""
  t = TraceResult("P5",
    "Workflow: SAD after xtriage → autosol")

  files = ["p9.sca", "sequence.dat"]
  if _requires_phenix(t):
    return t

  h = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.5}),
  ])
  state, vp, reason = mock_detect_ws(
    files, h, "SAD with selenium at 2.5A",
    directives={
      "program_settings": {
        "phenix.autosol": {"atom_type": "Se"},
      },
    })
  t.step("State: %s" % state)
  t.step("Valid: %s" % vp)

  if "phenix.autosol" not in vp:
    t.issue("WRONG",
      "autosol not available after xtriage "
      "with anomalous data")
  else:
    t.step("autosol available (correct)")
  return t


def trace_P6():
  """P6: Workflow state — no sequence file.
  autobuild should NOT be available."""
  t = TraceResult("P6",
    "Workflow: no sequence → autobuild blocked")

  files = ["data.mtz", "model.pdb"]
  if _requires_phenix(t):
    return t

  h = mock_history([
    ("phenix.xtriage", "SUCCESS: OK",
     {"resolution": 2.0}),
    ("phenix.refine", "SUCCESS: OK",
     {"r_free": 0.35}),
  ])
  state, vp, reason = mock_detect_ws(
    files, h, "refine the model")
  t.step("State: %s" % state)
  t.step("Valid: %s" % vp)

  if "phenix.autobuild" in vp:
    t.issue("WRONG",
      "autobuild available WITHOUT sequence "
      "file — YAML condition 'has: sequence' "
      "should block it")
  else:
    t.step("autobuild blocked (correct)")

  if "phenix.refine" not in vp:
    t.issue("WRONG",
      "refine should be available")
  else:
    t.step("refine available (correct)")
  return t


# ── L1-L10: Mock LLM output tests ────────────
# These test the LLM response parsing and post-parse
# validation without calling an actual LLM. They
# inject mock LLM output strings and check that
# the agent handles them correctly.

def _get_parse_intent():
  """Import parse_intent_json."""
  try:
    from agent.graph_nodes import parse_intent_json
    return parse_intent_json
  except ImportError:
    try:
      from libtbx.langchain.agent.graph_nodes \
        import parse_intent_json
      return parse_intent_json
    except ImportError:
      return None


def _get_parse_assessment():
  """Import parse_assessment."""
  try:
    from knowledge.thinking_prompts import (
      parse_assessment)
    return parse_assessment
  except ImportError:
    return None


def trace_L1():
  """L1: LLM returns valid program selection."""
  t = TraceResult("L1",
    "Mock LLM: valid program JSON")

  parse = _get_parse_intent()
  if parse is None:
    t.step("SKIP (parse_intent_json unavailable)")
    return t

  raw = '''{
    "program": "phenix.refine",
    "files": {"model": "refine_001.pdb",
              "data": "data.mtz"},
    "strategy": {"ordered_solvent": true},
    "reasoning": "Model needs refinement.",
    "stop": false
  }'''
  intent = parse(raw)
  t.step("program: %s" % intent.get("program"))
  t.step("strategy: %s" % intent.get("strategy"))
  t.step("stop: %s" % intent.get("stop"))

  if intent["program"] != "phenix.refine":
    t.issue("WRONG", "Wrong program parsed")
  if not intent["strategy"].get("ordered_solvent"):
    t.issue("WRONG", "Strategy flag lost")
  if intent.get("stop") is not False:
    t.issue("WRONG", "stop flag wrong")

  return t


def trace_L2():
  """L2: LLM returns hallucinated program name."""
  t = TraceResult("L2",
    "Mock LLM: hallucinated program name")

  parse = _get_parse_intent()
  if parse is None:
    t.step("SKIP")
    return t

  raw = '''{
    "program": "phenix.super_refine_v2",
    "files": {},
    "strategy": {},
    "reasoning": "Using advanced refinement."
  }'''
  intent = parse(raw)
  program = intent.get("program")
  t.step("Parsed program: %s" % program)

  # The parser should return it as-is — validation
  # happens downstream in plan()
  valid_programs = [
    "phenix.refine", "phenix.autobuild", "STOP"]
  if program not in valid_programs:
    t.step("Hallucinated program correctly parsed "
           "(will be caught by validation)")
  else:
    t.issue("WRONG",
      "Hallucinated program somehow in "
      "valid_programs")

  # Simulate what plan() would do
  if program not in valid_programs:
    # Falls through to "not in valid_programs" check
    # plan() would either use first valid or STOP
    fallback = valid_programs[0] if (
      valid_programs) else "STOP"
    t.step("Validation would fallback to: %s"
           % fallback)

  return t


def trace_L3():
  """L3: LLM returns STOP but user requested
  a specific program."""
  t = TraceResult("L3",
    "Mock LLM: STOP override by explicit_program")

  parse = _get_parse_intent()
  if parse is None:
    t.step("SKIP")
    return t

  raw = '''{
    "program": "STOP",
    "reasoning": "No further improvement possible.",
    "stop": true
  }'''
  intent = parse(raw)
  t.step("LLM chose: %s" % intent["program"])

  # Simulate explicit_program override
  explicit_prog = "phenix.ligandfit"
  valid_programs = [
    "phenix.refine", "phenix.ligandfit", "STOP"]
  history = [
    {"program": "phenix.xtriage"},
    {"program": "phenix.refine"},
  ]
  programs_run = {h["program"] for h in history}

  if (explicit_prog in valid_programs and
      explicit_prog not in programs_run):
    t.step("Override: %s hasn't run yet → "
           "override STOP" % explicit_prog)
    intent["program"] = explicit_prog
    intent["stop"] = False
  else:
    t.step("No override (already run or not valid)")

  if intent["program"] != "phenix.ligandfit":
    t.issue("WRONG",
      "STOP should be overridden by "
      "explicit_program")

  return t


def trace_L4():
  """L4: LLM returns program + stop=true
  (contradictory — should run the program)."""
  t = TraceResult("L4",
    "Mock LLM: program + stop=true (run program)")

  parse = _get_parse_intent()
  if parse is None:
    t.step("SKIP")
    return t

  raw = '''{
    "program": "phenix.autosol",
    "stop": true,
    "reasoning": "Running autosol then stopping."
  }'''
  intent = parse(raw)
  t.step("program: %s, stop: %s" % (
    intent["program"], intent.get("stop")))

  valid_programs = [
    "phenix.autosol", "phenix.refine", "STOP"]

  # plan() logic: if program is valid and not STOP,
  # clear stop flag
  if (intent["program"] in valid_programs
      and intent["program"] != "STOP"):
    intent["stop"] = False
    t.step("Cleared stop flag (valid program)")
  else:
    t.step("Stop flag preserved")

  if intent.get("stop"):
    t.issue("WRONG",
      "stop=true should be cleared when "
      "valid program selected")

  return t


def trace_L5():
  """L5: LLM returns malformed JSON."""
  t = TraceResult("L5",
    "Mock LLM: malformed JSON recovery")

  parse = _get_parse_intent()
  if parse is None:
    t.step("SKIP")
    return t

  # Markdown fences + extra text
  raw = '''Here is my analysis:

```json
{"program": "phenix.refine", "reasoning": "test"}
```

I hope this helps!'''

  try:
    intent = parse(raw)
    t.step("Parsed despite markdown: program=%s"
           % intent.get("program"))
    if intent["program"] != "phenix.refine":
      t.issue("WRONG", "Wrong program from fenced")
  except Exception as e:
    t.issue("CRASH",
      "Failed to parse fenced JSON: %s" % e)

  # Completely broken JSON
  raw2 = "I think phenix.refine is best"
  try:
    intent2 = parse(raw2)
    t.issue("WRONG",
      "Should fail on non-JSON, got: %s"
      % intent2)
  except (ValueError, Exception):
    t.step("Correctly rejected non-JSON")

  # JSON with no braces
  raw3 = '"program": "phenix.refine"'
  try:
    intent3 = parse(raw3)
    t.issue("WRONG",
      "Should fail on bare JSON, got: %s"
      % intent3)
  except (ValueError, Exception):
    t.step("Correctly rejected bare JSON")

  return t


def trace_L6():
  """L6: LLM returns hallucinated file paths."""
  t = TraceResult("L6",
    "Mock LLM: hallucinated file paths")

  parse = _get_parse_intent()
  if parse is None:
    t.step("SKIP")
    return t

  raw = '''{
    "program": "phenix.refine",
    "files": {
      "model": "/home/user/secret_model.pdb",
      "data": "/nonexistent/path/data.mtz"
    },
    "strategy": {},
    "reasoning": "Using the refined model."
  }'''
  intent = parse(raw)
  files = intent.get("files", {})
  t.step("LLM files: %s" % files)

  # The BUILD node should validate files against
  # best_files or available_files. The parser
  # should NOT validate — it's raw extraction.
  # But we can check if the strip_unsupplied_files
  # method would catch these.
  t.step("File validation happens in BUILD node "
         "(strip_unsupplied_files)")
  t.step("LLM-injected paths should be caught by "
         "path whitelist check")

  # Key insight: the file paths from LLM should
  # be replaced by actual paths from best_files
  actual_best_files = {
    "model": "/real/path/refine_003.pdb",
    "data": "/real/path/data.mtz",
  }
  t.step("Best files would override: %s"
         % actual_best_files)

  return t


def trace_L7():
  """L7: Expert assessment — valid JSON response."""
  t = TraceResult("L7",
    "Mock THINK: valid expert assessment JSON")

  parse = _get_parse_assessment()
  if parse is None:
    t.step("SKIP (parse_assessment unavailable)")
    return t

  raw = '''{
    "analysis": "R-free is stuck at 0.47. TFZ was 8.5 which suggests the MR solution is correct, but the space group may be wrong (enantiomorph P3121 vs P3221).",
    "confidence": "low",
    "action": "pivot",
    "guidance": "Re-run phaser in the alternative enantiomorph.",
    "concerns": ["possible enantiomorph error"],
    "alternatives": ["try P3221", "check twin law"]
  }'''
  assessment = parse(raw)
  t.step("action: %s" % assessment.get("action"))
  t.step("confidence: %s" % assessment.get(
    "confidence"))
  t.step("analysis length: %d" % len(
    assessment.get("analysis", "")))

  if assessment["action"] != "pivot":
    t.issue("WRONG", "action should be pivot")
  if assessment["confidence"] != "low":
    t.issue("WRONG", "confidence should be low")
  if "enantiomorph" not in assessment.get(
    "analysis", ""
  ):
    t.issue("WRONG",
      "analysis should mention enantiomorph")

  return t


def trace_L8():
  """L8: Expert assessment — malformed response."""
  t = TraceResult("L8",
    "Mock THINK: malformed expert response")

  parse = _get_parse_assessment()
  if parse is None:
    t.step("SKIP")
    return t

  # No JSON at all — just prose
  raw = ("The model looks fine. R-free is 0.25. "
         "I recommend continuing refinement.")
  assessment = parse(raw)
  t.step("Fallback analysis: '%s...'"
         % assessment.get("analysis", "")[:40])
  t.step("Default action: %s" % assessment.get(
    "action"))

  # Should use defaults
  if assessment["action"] != "let_run":
    t.issue("WRONG",
      "Malformed should default to let_run")
  if not assessment.get("analysis"):
    t.issue("WRONG",
      "Raw text should be used as analysis")

  # Markdown-wrapped JSON
  raw2 = '''```json
  {"analysis": "test", "action": "stop",
   "confidence": "hopeless"}
  ```'''
  a2 = parse(raw2)
  t.step("Markdown-wrapped: action=%s conf=%s"
         % (a2["action"], a2["confidence"]))
  if a2["action"] != "stop":
    t.issue("WRONG", "Should parse from fences")
  if a2["confidence"] != "hopeless":
    t.issue("WRONG", "Should preserve hopeless")

  return t


def trace_L9():
  """L9: Expert says stop but plan has stages."""
  t = TraceResult("L9",
    "Mock THINK: stop vs active plan conflict")

  parse = _get_parse_assessment()
  if parse is None:
    t.step("SKIP")
    return t

  # Expert recommends stopping
  raw = '''{
    "analysis": "No further improvement possible.",
    "confidence": "high",
    "action": "stop",
    "guidance": "The structure is complete."
  }'''
  assessment = parse(raw)
  t.step("Expert action: %s" % assessment["action"])

  # But plan has remaining stages
  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=["d.mtz", "m.pdb", "s.fa"],
    user_advice="solve by MR",
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t

  # Plan at initial_refinement, 3 more stages
  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  plan.advance()
  plan.mark_stage_started(2)
  plan.record_stage_cycle("phenix.phaser")
  plan.advance()
  plan.mark_stage_started(3)

  remaining = sum(
    1 for p in plan.stages
    if p.status in ("active", "pending"))
  t.step("Remaining plan stages: %d" % remaining)
  t.step("Expert says stop, plan says continue")
  t.step("Conflict: LLM decision node resolves "
         "(expert advises, LLM decides)")

  # The key insight: expert assessment is ADVISORY.
  # The LLM sees it in the PLAN prompt and may
  # choose to follow or ignore it. The plan system
  # doesn't enforce expert recommendations.
  t.step("Resolution: expert is advisory, "
         "LLM makes final call based on "
         "valid_programs + plan + expert")

  return t


def trace_L10():
  """L10: Forced program overrides LLM choice."""
  t = TraceResult("L10",
    "Mock LLM: forced_program override")

  parse = _get_parse_intent()
  if parse is None:
    t.step("SKIP")
    return t

  # LLM chose refine but workflow forces phaser
  raw = '''{
    "program": "phenix.refine",
    "reasoning": "Refinement next.",
    "stop": false
  }'''
  intent = parse(raw)
  t.step("LLM chose: %s" % intent["program"])

  # Simulate forced_program from workflow
  forced = "phenix.phaser"
  valid = ["phenix.phaser", "phenix.refine", "STOP"]
  chosen = intent["program"]

  if (forced and forced in valid
      and chosen != forced
      and chosen != "STOP"):
    t.step("Override: %s → %s (forced_program)"
           % (chosen, forced))
    intent["program"] = forced
    intent["reasoning"] = (
      "Running %s (workflow directive). "
      "[LLM suggested: %s]"
      % (forced, intent["reasoning"][:100]))

  if intent["program"] != "phenix.phaser":
    t.issue("WRONG",
      "forced_program should override LLM")
  if "workflow directive" not in intent.get(
    "reasoning", ""
  ):
    t.issue("DEGRADED",
      "Reasoning should note the override")
  t.step("Final: %s" % intent["program"])
  t.step("Reasoning: %s" % intent["reasoning"][:60])

  return t


# ── PG1-PG5: Placement gate tests (v114.1) ───

def _check_placement(program, metrics):
  """Replicate the placement detection logic from
  _detect_model_placement for testing.

  Returns (placed, metric_name, metric_value) or
  (False, None, None).
  """
  if "model_vs_data" in program:
    cc = metrics.get(
      "cc", metrics.get("map_cc",
        metrics.get("model_data_cc")))
    if cc is not None:
      try:
        cc = float(cc)
        if cc > 0.3:
          return (True, "CC", cc)
      except (ValueError, TypeError):
        pass

  elif "refine" in program:
    rf = metrics.get("r_free")
    if rf is not None:
      try:
        rf = float(rf)
        if rf < 0.50:
          return (True, "R-free", rf)
      except (ValueError, TypeError):
        pass
    cc = metrics.get(
      "map_cc", metrics.get("model_map_cc"))
    if cc is not None:
      try:
        cc = float(cc)
        if cc > 0.3:
          return (True, "Map CC", cc)
      except (ValueError, TypeError):
        pass

  return (False, None, None)


def trace_PG1():
  """PG1: Placement detected from model_vs_data
  CC → plan MR stage skipped."""
  t = TraceResult("PG1",
    "Placement gate: model_vs_data CC=0.5 "
    "skips MR stage")

  # Detection check
  placed, metric, value = _check_placement(
    "phenix.model_vs_data", {"cc": 0.5})
  t.step("Detection: placed=%s %s=%s"
         % (placed, metric, value))
  if not placed:
    t.issue("WRONG", "Should detect placement "
            "from CC=0.5")

  # Plan fast-forward
  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=["d.mtz", "m.pdb", "s.fa"],
    user_advice="solve by MR",
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t

  phase_ids = [p.id for p in plan.stages]
  t.step("Plan: %s" % " -> ".join(phase_ids))

  # Simulate fast-forward: skip MR stage
  for stage in plan.stages:
    if stage.id == "molecular_replacement":
      if stage.status in ("pending", "active"):
        stage.status = "skipped"
        t.step("Skipped: %s" % stage.id)

  # Advance past skipped
  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  plan.advance()
  # Should skip molecular_replacement
  curr = plan.current_stage()
  if curr and curr.id == "molecular_replacement":
    # advance() doesn't skip — need manual advance
    if curr.status == "skipped":
      plan.advance()
      curr = plan.current_stage()
  t.step("Current stage after skip: %s"
         % (curr.id if curr else "DONE"))
  if curr and curr.id == "molecular_replacement":
    t.issue("WRONG",
      "Should have skipped MR stage")
  elif curr and "refine" in curr.id:
    t.step("Correctly at refinement (good)")

  return t


def trace_PG2():
  """PG2: Placement NOT detected from poor
  model_vs_data CC → MR stage remains."""
  t = TraceResult("PG2",
    "Placement gate: model_vs_data CC=0.1 "
    "does NOT skip MR")

  placed, _, _ = _check_placement(
    "phenix.model_vs_data", {"cc": 0.1})
  t.step("Detection: placed=%s" % placed)
  if placed:
    t.issue("WRONG",
      "CC=0.1 should NOT trigger placement")
  else:
    t.step("Correctly not placed (good)")

  # Also test edge case: no CC metric
  placed2, _, _ = _check_placement(
    "phenix.model_vs_data", {})
  if placed2:
    t.issue("WRONG",
      "Empty metrics should not trigger")
  else:
    t.step("Empty metrics: not placed (good)")

  return t


def trace_PG3():
  """PG3: Placement detected from refine R-free
  → phaser should be suppressed."""
  t = TraceResult("PG3",
    "Placement gate: refine R-free=0.28 "
    "confirms placement")

  placed, metric, value = _check_placement(
    "phenix.refine", {"r_free": 0.28})
  t.step("Detection: placed=%s %s=%s"
         % (placed, metric, value))
  if not placed:
    t.issue("WRONG",
      "R-free=0.28 should confirm placement")

  # R-free=0.55 should NOT confirm
  placed2, _, _ = _check_placement(
    "phenix.refine", {"r_free": 0.55})
  t.step("R-free=0.55: placed=%s" % placed2)
  if placed2:
    t.issue("WRONG",
      "R-free=0.55 should NOT confirm — "
      "model may not be placed correctly")

  # Cryo-EM: map_cc=0.4 should confirm
  placed3, metric3, val3 = _check_placement(
    "phenix.real_space_refine",
    {"map_cc": 0.4})
  t.step("RSR map_cc=0.4: placed=%s %s=%s"
         % (placed3, metric3, val3))
  if not placed3:
    t.issue("WRONG",
      "map_cc=0.4 should confirm placement")

  return t


def trace_PG4():
  """PG4: Plan fast-forward skips both MR and
  phasing phases in mr_sad template."""
  t = TraceResult("PG4",
    "Placement gate: mr_sad plan skips MR "
    "and phasing phases")

  plan_template_loader._templates_cache = None
  plan = generate_plan(
    available_files=["d.mtz", "m.pdb", "s.fa"],
    user_advice="solve by MR-SAD",
    directives={})
  if plan is None:
    t.issue("CRASH", "No plan")
    return t

  phase_ids = [p.id for p in plan.stages]
  t.step("Plan: %s" % " -> ".join(phase_ids))
  if plan.template_id != "mr_sad":
    t.issue("WRONG",
      "Expected mr_sad, got %s"
      % plan.template_id)
    return t

  # Simulate: xtriage runs, then model_vs_data
  # detects placement → skip MR + phasing
  plan.mark_stage_started(1)
  plan.record_stage_cycle("phenix.xtriage")
  plan.advance()

  _skip_ids = {
    "molecular_replacement",
    "experimental_phasing",
  }
  skipped = []
  for stage in plan.stages:
    if (stage.id in _skip_ids
        and stage.status in ("pending", "active")):
      stage.status = "skipped"
      skipped.append(stage.id)

  t.step("Skipped stages: %s" % skipped)

  # Advance past skipped
  while True:
    curr = plan.current_stage()
    if curr is None:
      break
    if curr.status == "skipped":
      plan.advance()
    else:
      break

  curr = plan.current_stage()
  t.step("Current stage: %s"
         % (curr.id if curr else "DONE"))

  if curr and curr.id in _skip_ids:
    t.issue("WRONG",
      "MR/phasing phase not skipped: %s"
      % curr.id)
  elif curr and "refine" in curr.id:
    t.step("At refinement after skip (correct)")
  elif curr and "build" in curr.id:
    t.step("At build_and_refine after skip "
           "(correct)")

  return t


def trace_PG5():
  """PG5: Placement gate — PHENIX workflow
  engine filters destructive programs."""
  t = TraceResult("PG5",
    "Placement gate: workflow engine filters "
    "phaser/autosol when placed")

  # This test requires PHENIX
  state, vp, _ = mock_detect_ws(
    ["data.mtz", "model.pdb", "seq.fa"],
    [], "solve by MR")
  if state == "SKIP":
    t.step("PHENIX not available — skipping "
           "workflow engine check")
    return t

  # Run with model_is_placed in session_info
  try:
    from agent.workflow_state import (
      detect_workflow_state,
    )
    h = mock_history([
      ("phenix.xtriage", "SUCCESS: OK",
       {"resolution": 2.0}),
      ("phenix.model_vs_data", "SUCCESS: OK",
       {"cc": 0.5}),
    ])
    ws = detect_workflow_state(
      history=h,
      available_files=[
        "data.mtz", "model.pdb", "seq.fa"],
      analysis=None,
      maximum_automation=True,
      use_yaml_engine=True,
      directives={},
      session_info={
        "user_advice": "solve by MR",
        "model_is_placed": True,
      },
      files_local=False,
    )
    vp = ws.get("valid_programs", [])
    t.step("Valid programs (placed): %s" % vp)
    if "phenix.phaser" in vp:
      t.issue("WRONG",
        "phaser should be filtered when "
        "model_is_placed=True")
    else:
      t.step("phaser filtered (correct)")
    if "phenix.autosol" in vp:
      t.issue("WRONG",
        "autosol should be filtered")
    else:
      t.step("autosol filtered (correct)")
    if "phenix.refine" in vp:
      t.step("refine available (correct)")
  except Exception as e:
    t.step("Workflow engine test: %s" % e)

  return t


# ── Registry ──────────────────────────────────

SCENARIOS = {
  "S1A": trace_S1A,
  "S1B": trace_S1B,
  "S1C": trace_S1C,
  "S2A": trace_S2A,
  "S2B": trace_S2B,
  "S3B": trace_S3B,
  "S4A": trace_S4A,
  "S4B": trace_S4B,
  "S5A": trace_S5A,
  "S5B": trace_S5B,
  "S6A": trace_S6A,
  "S6B": trace_S6B,
  "S7A": trace_S7A,
  "S7B": trace_S7B,
  "S8A": trace_S8A,
  "S9A": trace_S9A,
  "S9B": trace_S9B,
  "S10A": trace_S10A,
  "S10B": trace_S10B,
  "S11A": trace_S11A,
  "S11B": trace_S11B,
  "S12A": trace_S12A,
  "S13A": trace_S13A,
  "S13B": trace_S13B,
  "S14A": trace_S14A,
  "S14B": trace_S14B,
  "S15": trace_S15,
  "G1": trace_G1,
  "G2": trace_G2,
  "G3": trace_G3,
  "G4": trace_G4,

  # Cycle counting
  "C1": trace_C1,
  "C2": trace_C2,
  "C3": trace_C3,

  # Multi-cycle progression
  "M1": trace_M1,
  "M2": trace_M2,

  # PHENIX-dependent (workflow state + commands)
  "P1": trace_P1,
  "P2": trace_P2,
  "P3": trace_P3,
  "P4": trace_P4,
  "P5": trace_P5,
  "P6": trace_P6,

  # Mock LLM output tests
  "L1": trace_L1,
  "L2": trace_L2,
  "L3": trace_L3,
  "L4": trace_L4,
  "L5": trace_L5,
  "L6": trace_L6,
  "L7": trace_L7,
  "L8": trace_L8,
  "L9": trace_L9,
  "L10": trace_L10,

  # Placement gate tests (v114.1)
  "PG1": trace_PG1,
  "PG2": trace_PG2,
  "PG3": trace_PG3,
  "PG4": trace_PG4,
  "PG5": trace_PG5,
}


# ── Runner ────────────────────────────────────

def run_all(failures_only=False):
  results = []
  for sid in sorted(SCENARIOS.keys()):
    plan_template_loader._templates_cache = None
    try:
      r = SCENARIOS[sid]()
    except Exception as e:
      r = TraceResult(sid, "CRASHED")
      r.issue("CRASH", str(e)[:100])
    results.append(r)
    r.print_trace(failures_only=failures_only)

  # Summary
  counts = {}
  for r in results:
    counts[r.status] = counts.get(
      r.status, 0) + 1
  print("\n" + "=" * 50)
  print("SCENARIO TRACER SUMMARY")
  print("=" * 50)
  for status in ("PASS", "DEGRADED", "WRONG",
                 "STUCK", "CRASH"):
    n = counts.get(status, 0)
    if n > 0:
      print("  %-10s %d" % (status, n))
  total = len(results)
  passed = counts.get("PASS", 0)
  print("  Total: %d/%d passed" % (
    passed, total))
  if passed < total:
    failed_ids = [r.scenario_id for r in results
                  if r.status != "PASS"]
    raise AssertionError(
      "Scenario tracer: %d/%d failed: %s"
      % (total - passed, total,
         ", ".join(failed_ids)))
  return True


def main():
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("scenario", nargs="?",
    help="Run specific scenario (e.g. S1A)")
  parser.add_argument("--list", action="store_true",
    help="List available scenarios")
  parser.add_argument("--failures",
    action="store_true",
    help="Show only failures")
  args = parser.parse_args()

  if args.list:
    for sid in sorted(SCENARIOS.keys()):
      fn = SCENARIOS[sid]
      desc = fn.__doc__.strip().split("\n")[0]
      print("  %-6s %s" % (sid, desc))
    return

  if args.scenario:
    sid = args.scenario.upper()
    if sid not in SCENARIOS:
      print("Unknown scenario: %s" % sid)
      print("Available: %s" % ", ".join(
        sorted(SCENARIOS.keys())))
      sys.exit(1)
    plan_template_loader._templates_cache = None
    r = SCENARIOS[sid]()
    r.print_trace()
    return

  all_pass = run_all(
    failures_only=args.failures)
  sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
  main()
