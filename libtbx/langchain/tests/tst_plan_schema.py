"""
Unit tests for plan_schema.py (Stage 2, Step 2.1).

Run standalone:
  python tests/tst_plan_schema.py

Run from PHENIX:
  libtbx.python tests/tst_plan_schema.py

No PHENIX dependencies — pure Python.

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import json
import os
import sys
import traceback

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_THIS_DIR)
if _PROJECT_ROOT not in sys.path:
  sys.path.insert(0, _PROJECT_ROOT)

from knowledge.plan_schema import (
  StageDef, StructurePlan, merge_directives,
  STAGE_PENDING, STAGE_ACTIVE, STAGE_COMPLETE,
  STAGE_SKIPPED,
)


def run_tests():
  """Run all plan schema tests."""
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
  print("Plan Schema Unit Tests")
  print("=" * 60)
  print()

  # --- StageDef ---
  print("StageDef")
  test("phase_init_defaults",
    test_phase_init_defaults)
  test("phase_init_full",
    test_phase_init_full)
  test("phase_roundtrip",
    test_phase_roundtrip)
  test("phase_from_dict_tolerant",
    test_phase_from_dict_tolerant)
  test("phase_runtime_state_persists",
    test_phase_runtime_state_persists)
  print()

  # --- StructurePlan basics ---
  print("StructurePlan basics")
  test("plan_init_empty",
    test_plan_init_empty)
  test("plan_current_stage",
    test_plan_current_stage)
  test("plan_current_stage_empty",
    test_plan_current_stage_empty)
  print()

  # --- Navigation ---
  print("Navigation")
  test("advance_basic",
    test_advance_basic)
  test("advance_to_end",
    test_advance_to_end)
  test("advance_skips_completed",
    test_advance_skips_completed)
  test("retreat_basic",
    test_retreat_basic)
  test("retreat_not_found",
    test_retreat_not_found)
  test("retreat_resets_current_to_pending",
    test_retreat_resets_current_to_pending)
  test("retreat_count_tracks",
    test_retreat_count_tracks)
  test("retreat_resets_downstream",
    test_retreat_resets_downstream)
  test("retreat_preserves_skipped",
    test_retreat_preserves_skipped)
  test("skip_stage",
    test_skip_stage)
  test("mark_stage_started",
    test_mark_stage_started)
  test("record_stage_cycle",
    test_record_stage_cycle)
  test("mark_phase_complete",
    test_mark_phase_complete)
  test("is_complete",
    test_is_complete)
  test("get_stage_by_id",
    test_get_stage_by_id)
  test("get_previous_phase",
    test_get_previous_phase)
  test("is_exhausted",
    test_is_exhausted)
  test("can_retreat_allowed",
    test_can_retreat_allowed)
  test("can_retreat_max_reached",
    test_can_retreat_max_reached)
  test("can_retreat_cooldown",
    test_can_retreat_cooldown)
  test("retreat_records_cycle",
    test_retreat_records_cycle)
  print()

  # --- Directives ---
  print("Directives")
  test("to_directives_single_program",
    test_to_directives_single_program)
  test("to_directives_multi_program",
    test_to_directives_multi_program)
  test("to_directives_with_strategy",
    test_to_directives_with_strategy)
  test("to_directives_with_phase_directives",
    test_to_directives_with_phase_directives)
  test("to_directives_empty_plan",
    test_to_directives_empty_plan)
  print()

  # --- Hash ---
  print("Hash")
  test("compute_hash_stable",
    test_compute_hash_stable)
  test("compute_hash_changes_on_advance",
    test_compute_hash_changes_on_advance)
  test("compute_hash_changes_on_retreat",
    test_compute_hash_changes_on_retreat)
  print()

  # --- Merge directives ---
  print("Merge directives")
  test("merge_no_conflict",
    test_merge_no_conflict)
  test("merge_user_skip_conflicts",
    test_merge_user_skip_conflicts)
  test("merge_user_overrides_after_program",
    test_merge_user_overrides_after_program)
  test("merge_user_program_settings",
    test_merge_user_program_settings)
  test("merge_plan_empty",
    test_merge_plan_empty)
  test("merge_user_empty",
    test_merge_user_empty)
  test("merge_preserves_plan_prefer",
    test_merge_preserves_plan_prefer)
  print()

  # --- Display ---
  print("Display")
  test("get_display_phases",
    test_get_display_phases)
  test("format_plan_header",
    test_format_plan_header)
  print()

  # --- Serialization ---
  print("Serialization")
  test("plan_roundtrip_empty",
    test_plan_roundtrip_empty)
  test("plan_roundtrip_populated",
    test_plan_roundtrip_populated)
  test("plan_roundtrip_through_json",
    test_plan_roundtrip_through_json)
  test("plan_from_dict_tolerant",
    test_plan_from_dict_tolerant)
  test("plan_from_dict_none",
    test_plan_from_dict_none)
  print()

  # --- Real workflow scenario ---
  print("Real workflow scenario")
  test("mr_refine_workflow",
    test_mr_refine_workflow)
  test("retreat_and_retry_workflow",
    test_retreat_and_retry_workflow)
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

def _make_mr_refine_plan():
  """Build a realistic MR+refine plan."""
  stages = [
    StageDef(
      id="data_assessment",
      programs=["phenix.xtriage"],
      max_cycles=1,
      success_criteria={"xtriage_completed": "true"},
      description="Analyze data quality",
      provides=["resolution", "twinning_status"],
    ),
    StageDef(
      id="molecular_replacement",
      programs=["phenix.phaser"],
      max_cycles=1,
      success_criteria={"tfz": ">8", "llg": ">100"},
      gate_conditions=[
        {"if": "tfz < 5",
         "action": "stop_report"},
      ],
      fallbacks=[
        {"if": "tfz < 6",
         "action": "try_alphafold"},
      ],
      description="Find MR solution",
    ),
    StageDef(
      id="initial_refinement",
      programs=["phenix.refine"],
      max_cycles=3,
      success_criteria={"r_free": "<0.35"},
      gate_conditions=[
        {"if": "r_free > 0.45 after 2 cycles",
         "action": "retreat_to molecular_replacement"},
      ],
      description="Initial refinement",
    ),
    StageDef(
      id="model_rebuilding",
      programs=["phenix.autobuild"],
      max_cycles=1,
      success_criteria={"r_free": "<0.30"},
      skip_if="r_free < 0.28",
      strategy={"rebuild_in_place": False},
      description="Rebuild model",
    ),
    StageDef(
      id="final_refinement",
      programs=["phenix.refine"],
      max_cycles=3,
      success_criteria={"r_free": "<0.25"},
      strategy={"ordered_solvent": True},
      description="Final refinement",
    ),
  ]
  return StructurePlan(
    goal="Solve kinase by MR at 2.1A",
    stages=stages,
    template_id="mr_refine",
    created_at_cycle=0,
  )


# ── StageDef tests ────────────────────────────────────

def test_phase_init_defaults():
  p = StageDef(id="test")
  assert p.id == "test"
  assert p.programs == []
  assert p.max_cycles == 5
  assert p.success_criteria == {}
  assert p.gate_conditions == []
  assert p.fallbacks == []
  assert p.skip_if == ""
  assert p.directives == {}
  assert p.strategy == {}
  assert p.status == STAGE_PENDING
  assert p.cycles_used == 0
  assert p.start_cycle is None


def test_phase_init_full():
  p = StageDef(
    id="refine",
    programs=["phenix.refine"],
    max_cycles=3,
    success_criteria={"r_free": "<0.35"},
    gate_conditions=[
      {"if": "stalled", "action": "retreat"}
    ],
    fallbacks=[{"if": "bad", "action": "rebuild"}],
    skip_if="r_free < 0.28",
    strategy={"ordered_solvent": True},
    description="Refine the model",
  )
  assert p.id == "refine"
  assert p.programs == ["phenix.refine"]
  assert p.max_cycles == 3
  assert p.success_criteria == {"r_free": "<0.35"}
  assert len(p.gate_conditions) == 1
  assert len(p.fallbacks) == 1
  assert p.skip_if == "r_free < 0.28"
  assert p.strategy["ordered_solvent"] is True


def test_phase_roundtrip():
  p = StageDef(
    id="mr",
    programs=["phenix.phaser"],
    max_cycles=2,
    success_criteria={"tfz": ">8"},
    gate_conditions=[
      {"if": "tfz < 5", "action": "stop"}
    ],
    strategy={"ensemble": True},
    description="Find MR solution",
    provides=["mr_solution"],
    if_skipped={"mr_solution": "unavailable"},
  )
  p.status = STAGE_COMPLETE
  p.cycles_used = 1
  p.start_cycle = 2
  p.end_cycle = 2
  p.result_metrics = {"tfz": 14.2}

  d = p.to_dict()
  p2 = StageDef.from_dict(d)
  assert p2.id == "mr"
  assert p2.programs == ["phenix.phaser"]
  assert p2.success_criteria == {"tfz": ">8"}
  assert p2.status == STAGE_COMPLETE
  assert p2.cycles_used == 1
  assert p2.start_cycle == 2
  assert p2.end_cycle == 2
  assert p2.result_metrics == {"tfz": 14.2}
  assert p2.provides == ["mr_solution"]
  assert p2.if_skipped == {
    "mr_solution": "unavailable"
  }


def test_phase_from_dict_tolerant():
  p = StageDef.from_dict({})
  assert p.id == "unknown"
  assert p.programs == []
  p = StageDef.from_dict("garbage")
  assert p.id == "unknown"


def test_phase_runtime_state_persists():
  """Runtime state (status, cycles) survives
  serialization."""
  p = StageDef(id="test", max_cycles=3)
  p.status = STAGE_ACTIVE
  p.cycles_used = 2
  p.start_cycle = 5
  d = p.to_dict()
  p2 = StageDef.from_dict(d)
  assert p2.status == STAGE_ACTIVE
  assert p2.cycles_used == 2
  assert p2.start_cycle == 5


# ── StructurePlan basics ─────────────────────────────

def test_plan_init_empty():
  plan = StructurePlan()
  assert plan.goal == ""
  assert plan.stages == []
  assert plan.current_stage_index == 0
  assert plan.retreat_count == 0


def test_plan_current_stage():
  plan = _make_mr_refine_plan()
  curr = plan.current_stage()
  assert curr is not None
  assert curr.id == "data_assessment"


def test_plan_current_stage_empty():
  plan = StructurePlan()
  assert plan.current_stage() is None


# ── Navigation ────────────────────────────────────────

def test_advance_basic():
  plan = _make_mr_refine_plan()
  plan.stages[0].status = STAGE_ACTIVE
  # Advance from data_assessment to MR
  ok = plan.advance()
  assert ok is True
  assert plan.current_stage().id == (
    "molecular_replacement"
  )
  assert plan.stages[0].status == STAGE_COMPLETE
  assert plan.stages[1].status == STAGE_ACTIVE


def test_advance_to_end():
  plan = _make_mr_refine_plan()
  # Advance through all stages
  for stage in plan.stages:
    stage.status = STAGE_ACTIVE
    plan.current_stage_index = (
      plan.stages.index(stage)
    )
    ok = plan.advance()
    if stage is plan.stages[-1]:
      assert ok is False, (
        "Should return False at end"
      )
    else:
      assert ok is True


def test_advance_skips_completed():
  """Advance skips over already-completed stages."""
  plan = _make_mr_refine_plan()
  plan.stages[0].status = STAGE_ACTIVE
  plan.stages[1].status = STAGE_COMPLETE  # MR done
  ok = plan.advance()
  assert ok is True
  # Should skip MR, land on initial_refinement
  assert plan.current_stage().id == (
    "initial_refinement"
  )


def test_retreat_basic():
  plan = _make_mr_refine_plan()
  # Start at initial_refinement (index 2)
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  ok = plan.retreat_to("molecular_replacement")
  assert ok is True
  assert plan.current_stage().id == (
    "molecular_replacement"
  )
  assert plan.stages[1].status == STAGE_ACTIVE
  # Target stage should be reset
  assert plan.stages[1].cycles_used == 0


def test_retreat_not_found():
  plan = _make_mr_refine_plan()
  ok = plan.retreat_to("nonexistent")
  assert ok is False


def test_retreat_resets_current_to_pending():
  """After retreat, the previously-active stage is
  reset to PENDING (retriable with new strategy)."""
  plan = _make_mr_refine_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.retreat_to("molecular_replacement")
  # Stage 2 is downstream of retreat target (1),
  # so it's reset to PENDING for retry.
  assert plan.stages[2].status == STAGE_PENDING


def test_retreat_count_tracks():
  plan = _make_mr_refine_plan()
  assert plan.retreat_count == 0
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.retreat_to("molecular_replacement")
  assert plan.retreat_count == 1
  # Second retreat
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.retreat_to("data_assessment")
  assert plan.retreat_count == 2


def test_retreat_resets_downstream():
  """Retreat resets all stages after the target to
  PENDING so they can re-run with the new strategy."""
  plan = _make_mr_refine_plan()
  # Simulate: stages 0-2 done, stage 3 active
  plan.stages[0].status = STAGE_COMPLETE
  plan.stages[0].cycles_used = 1
  plan.stages[1].status = STAGE_COMPLETE
  plan.stages[1].cycles_used = 1
  plan.stages[1].result_metrics = {"tfz": 7.0}
  plan.stages[2].status = STAGE_COMPLETE
  plan.stages[2].cycles_used = 3
  plan.stages[3].status = STAGE_ACTIVE
  plan.stages[3].cycles_used = 1
  plan.current_stage_index = 3
  # Retreat to MR (index 1)
  plan.retreat_to("molecular_replacement")
  # Target (index 1) is ACTIVE and reset
  assert plan.stages[1].status == STAGE_ACTIVE
  assert plan.stages[1].cycles_used == 0
  assert plan.stages[1].result_metrics == {}
  # Stages 2,3 are reset to PENDING
  assert plan.stages[2].status == STAGE_PENDING
  assert plan.stages[2].cycles_used == 0
  assert plan.stages[3].status == STAGE_PENDING
  assert plan.stages[3].cycles_used == 0
  # Stage 0 (before target) is untouched
  assert plan.stages[0].status == STAGE_COMPLETE
  assert plan.stages[0].cycles_used == 1
  # Stage 4 (never started) is still PENDING
  assert plan.stages[4].status == STAGE_PENDING


def test_retreat_preserves_skipped():
  """Retreat should NOT reset stages that were
  explicitly skipped (they were skipped for a reason
  independent of the failed strategy)."""
  plan = _make_mr_refine_plan()
  plan.stages[0].status = STAGE_COMPLETE
  plan.stages[1].status = STAGE_COMPLETE
  plan.stages[2].status = STAGE_COMPLETE
  plan.stages[3].status = STAGE_SKIPPED  # explicitly
  plan.stages[4].status = STAGE_ACTIVE
  plan.current_stage_index = 4
  plan.retreat_to("initial_refinement")
  # Skipped stage should stay skipped
  assert plan.stages[3].status == STAGE_SKIPPED
  # But the active/complete ones after target
  # should be reset
  assert plan.stages[4].status == STAGE_PENDING


def test_skip_stage():
  plan = _make_mr_refine_plan()
  ok = plan.skip_stage("model_rebuilding",
                       reason="r_free already <0.28")
  assert ok is True
  stage = plan.get_stage_by_id("model_rebuilding")
  assert stage.status == STAGE_SKIPPED


def test_mark_stage_started():
  plan = _make_mr_refine_plan()
  plan.mark_stage_started(1)
  curr = plan.current_stage()
  assert curr.status == STAGE_ACTIVE
  assert curr.start_cycle == 1
  # Second call doesn't overwrite start_cycle
  plan.mark_stage_started(2)
  assert curr.start_cycle == 1


def test_record_stage_cycle():
  plan = _make_mr_refine_plan()
  plan.mark_stage_started(1)
  plan.record_stage_cycle()
  plan.record_stage_cycle()
  assert plan.current_stage().cycles_used == 2


def test_mark_phase_complete():
  plan = _make_mr_refine_plan()
  plan.mark_stage_started(1)
  plan.mark_phase_complete(3,
                           metrics={"r_free": 0.34})
  curr = plan.current_stage()
  assert curr.status == STAGE_COMPLETE
  assert curr.end_cycle == 3
  assert curr.result_metrics == {"r_free": 0.34}


def test_is_complete():
  plan = _make_mr_refine_plan()
  assert plan.is_complete() is False
  # Mark all as complete
  for stage in plan.stages:
    stage.status = STAGE_COMPLETE
  assert plan.is_complete() is True
  # With skipped
  plan.stages[3].status = STAGE_SKIPPED
  assert plan.is_complete() is True
  # With one pending
  plan.stages[4].status = STAGE_PENDING
  assert plan.is_complete() is False


def test_get_stage_by_id():
  plan = _make_mr_refine_plan()
  p = plan.get_stage_by_id("initial_refinement")
  assert p is not None
  assert p.id == "initial_refinement"
  assert plan.get_stage_by_id("nope") is None


def test_get_previous_phase():
  plan = _make_mr_refine_plan()
  plan.current_stage_index = 2
  prev = plan.get_previous_phase()
  assert prev is not None
  assert prev.id == "molecular_replacement"
  # At index 0, no previous
  plan.current_stage_index = 0
  assert plan.get_previous_phase() is None


def test_is_exhausted():
  p = StageDef(id="test", max_cycles=3)
  assert p.is_exhausted() is False
  p.cycles_used = 2
  assert p.is_exhausted() is False
  p.cycles_used = 3
  assert p.is_exhausted() is True
  p.cycles_used = 5
  assert p.is_exhausted() is True


def test_can_retreat_allowed():
  plan = _make_mr_refine_plan()
  allowed, reason = plan.can_retreat(
    cycle_number=5
  )
  assert allowed is True
  assert reason == ""


def test_can_retreat_max_reached():
  plan = _make_mr_refine_plan()
  plan.retreat_count = 2
  allowed, reason = plan.can_retreat(
    cycle_number=10, max_retreats=2
  )
  assert allowed is False
  assert "max retreats" in reason


def test_can_retreat_cooldown():
  plan = _make_mr_refine_plan()
  plan.last_retreat_cycle = 5
  allowed, reason = plan.can_retreat(
    cycle_number=6, cooldown=2
  )
  assert allowed is False
  assert "cooldown" in reason
  # After enough cycles, retreat is allowed
  allowed2, _ = plan.can_retreat(
    cycle_number=7, cooldown=2
  )
  assert allowed2 is True


def test_retreat_records_cycle():
  plan = _make_mr_refine_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  plan.retreat_to(
    "molecular_replacement", cycle_number=8
  )
  assert plan.last_retreat_cycle == 8
  # Persists through round-trip
  d = plan.to_dict()
  plan2 = StructurePlan.from_dict(d)
  assert plan2.last_retreat_cycle == 8


# ── Directives ────────────────────────────────────────

def test_to_directives_single_program():
  plan = _make_mr_refine_plan()
  # data_assessment: single program phenix.xtriage
  directives = plan.to_directives()
  wf = directives.get("workflow_preferences", {})
  assert "phenix.xtriage" in wf.get(
    "prefer_programs", []
  )
  sc = directives.get("stop_conditions", {})
  assert sc.get("after_program") == (
    "phenix.xtriage"
  )


def test_to_directives_multi_program():
  """Multi-program stage should not set after_program."""
  plan = StructurePlan(
    goal="test",
    stages=[StageDef(
      id="build",
      programs=[
        "phenix.predict_and_build",
        "phenix.map_to_model",
      ],
      max_cycles=2,
    )],
  )
  directives = plan.to_directives()
  sc = directives.get("stop_conditions", {})
  assert "after_program" not in sc


def test_to_directives_with_strategy():
  plan = _make_mr_refine_plan()
  # model_rebuilding (index 3) has strategy
  plan.current_stage_index = 3
  plan.stages[3].status = STAGE_ACTIVE
  directives = plan.to_directives()
  ps = directives.get("program_settings", {})
  assert "phenix.autobuild" in ps
  assert ps["phenix.autobuild"][
    "rebuild_in_place"
  ] is False


def test_to_directives_with_phase_directives():
  """Stage-level directives merge into output."""
  stage = StageDef(
    id="test",
    programs=["phenix.refine"],
    directives={
      "workflow_preferences": {
        "start_with_program": "phenix.refine",
      },
      "program_settings": {
        "phenix.refine": {"nproc": 4},
      },
    },
  )
  plan = StructurePlan(
    goal="test", stages=[stage],
  )
  d = plan.to_directives()
  wf = d.get("workflow_preferences", {})
  assert wf.get("start_with_program") == (
    "phenix.refine"
  )
  ps = d.get("program_settings", {})
  assert ps.get("phenix.refine", {}).get(
    "nproc"
  ) == 4


def test_to_directives_empty_plan():
  plan = StructurePlan()
  assert plan.to_directives() == {}


# ── Hash ──────────────────────────────────────────────

def test_compute_hash_stable():
  plan = _make_mr_refine_plan()
  h1 = plan.compute_hash()
  h2 = plan.compute_hash()
  assert h1 == h2
  assert len(h1) == 12


def test_compute_hash_changes_on_advance():
  plan = _make_mr_refine_plan()
  plan.stages[0].status = STAGE_ACTIVE
  h1 = plan.compute_hash()
  plan.advance()
  h2 = plan.compute_hash()
  assert h1 != h2


def test_compute_hash_changes_on_retreat():
  plan = _make_mr_refine_plan()
  plan.current_stage_index = 2
  plan.stages[2].status = STAGE_ACTIVE
  h1 = plan.compute_hash()
  plan.retreat_to("molecular_replacement")
  h2 = plan.compute_hash()
  assert h1 != h2


# ── Merge directives ─────────────────────────────────

def test_merge_no_conflict():
  plan_d = {
    "workflow_preferences": {
      "prefer_programs": ["phenix.refine"],
    },
    "stop_conditions": {
      "after_program": "phenix.refine",
    },
  }
  user_d = {
    "program_settings": {
      "phenix.refine": {"nproc": 8},
    },
  }
  merged, warnings = merge_directives(
    plan_d, user_d
  )
  assert "phenix.refine" in merged[
    "workflow_preferences"
  ]["prefer_programs"]
  assert merged["program_settings"][
    "phenix.refine"
  ]["nproc"] == 8
  assert len(warnings) == 0


def test_merge_user_skip_conflicts():
  plan_d = {
    "workflow_preferences": {
      "prefer_programs": ["phenix.xtriage"],
    },
  }
  user_d = {
    "workflow_preferences": {
      "skip_programs": ["phenix.xtriage"],
    },
  }
  merged, warnings = merge_directives(
    plan_d, user_d
  )
  assert len(warnings) == 1
  assert "xtriage" in warnings[0]


def test_merge_user_overrides_after_program():
  plan_d = {
    "stop_conditions": {
      "after_program": "phenix.refine",
    },
  }
  user_d = {
    "stop_conditions": {
      "after_program": "phenix.autobuild",
    },
  }
  merged, warnings = merge_directives(
    plan_d, user_d
  )
  assert merged["stop_conditions"][
    "after_program"
  ] == "phenix.autobuild"
  assert len(warnings) == 1


def test_merge_user_program_settings():
  plan_d = {
    "program_settings": {
      "phenix.refine": {
        "ordered_solvent": True,
      },
    },
  }
  user_d = {
    "program_settings": {
      "phenix.refine": {"nproc": 4},
    },
  }
  merged, warnings = merge_directives(
    plan_d, user_d
  )
  ps = merged["program_settings"]["phenix.refine"]
  # Both settings present
  assert ps["ordered_solvent"] is True
  assert ps["nproc"] == 4


def test_merge_plan_empty():
  merged, warnings = merge_directives(
    {}, {"stop_conditions": {"after_cycle": 5}}
  )
  assert merged["stop_conditions"][
    "after_cycle"
  ] == 5
  assert len(warnings) == 0


def test_merge_user_empty():
  plan_d = {
    "workflow_preferences": {
      "prefer_programs": ["phenix.refine"],
    },
  }
  merged, warnings = merge_directives(plan_d, {})
  assert "phenix.refine" in merged[
    "workflow_preferences"
  ]["prefer_programs"]


def test_merge_preserves_plan_prefer():
  """User workflow_preferences without prefer_programs
  should not clobber plan's prefer_programs."""
  plan_d = {
    "workflow_preferences": {
      "prefer_programs": ["phenix.phaser"],
    },
  }
  user_d = {
    "workflow_preferences": {
      "skip_programs": ["phenix.xtriage"],
    },
  }
  merged, warnings = merge_directives(
    plan_d, user_d
  )
  wf = merged["workflow_preferences"]
  # Plan's prefer_programs preserved
  assert "phenix.phaser" in wf.get(
    "prefer_programs", []
  )
  # User's skip_programs also present
  assert "phenix.xtriage" in wf.get(
    "skip_programs", []
  )


# ── Display ───────────────────────────────────────────

def test_get_display_phases():
  plan = _make_mr_refine_plan()
  plan.stages[0].status = STAGE_COMPLETE
  plan.stages[1].status = STAGE_ACTIVE
  display = plan.get_display_phases()
  assert len(display) == 5
  assert display[0]["status"] == STAGE_COMPLETE
  assert display[1]["status"] == STAGE_ACTIVE
  assert display[2]["status"] == STAGE_PENDING
  assert display[0]["id"] == "data_assessment"
  assert "phenix.xtriage" in display[0]["programs"]


def test_format_plan_header():
  plan = _make_mr_refine_plan()
  plan.stages[0].status = STAGE_COMPLETE
  plan.stages[1].status = STAGE_ACTIVE
  header = plan.format_plan_header()
  assert "STRATEGY PLAN" in header
  assert "kinase" in header
  # Check status indicators
  lines = header.split("\n")
  phase_lines = [
    l for l in lines if "Stage " in l
  ]
  assert len(phase_lines) == 5
  # First stage should have ✓
  assert "✓" in phase_lines[0]
  # Second should have ●
  assert "●" in phase_lines[1]
  # Rest should have ○
  assert "○" in phase_lines[2]


# ── Serialization ─────────────────────────────────────

def test_plan_roundtrip_empty():
  plan = StructurePlan()
  d = plan.to_dict()
  plan2 = StructurePlan.from_dict(d)
  assert plan2.goal == ""
  assert plan2.stages == []


def test_plan_roundtrip_populated():
  plan = _make_mr_refine_plan()
  plan.stages[0].status = STAGE_COMPLETE
  plan.stages[0].cycles_used = 1
  plan.stages[0].start_cycle = 1
  plan.stages[0].end_cycle = 1
  plan.current_stage_index = 1
  plan.stages[1].status = STAGE_ACTIVE
  plan.retreat_count = 1
  plan.last_retreat_cycle = 5

  d = plan.to_dict()
  plan2 = StructurePlan.from_dict(d)

  assert plan2.goal == "Solve kinase by MR at 2.1A"
  assert len(plan2.stages) == 5
  assert plan2.stages[0].status == STAGE_COMPLETE
  assert plan2.stages[0].cycles_used == 1
  assert plan2.current_stage_index == 1
  assert plan2.stages[1].status == STAGE_ACTIVE
  assert plan2.template_id == "mr_refine"
  assert plan2.retreat_count == 1
  assert plan2.last_retreat_cycle == 5
  # Verify stage content survived
  p2_mr = plan2.stages[1]
  assert p2_mr.id == "molecular_replacement"
  assert "phenix.phaser" in p2_mr.programs
  assert p2_mr.success_criteria == {
    "tfz": ">8", "llg": ">100"
  }


def test_plan_roundtrip_through_json():
  plan = _make_mr_refine_plan()
  plan.stages[0].status = STAGE_COMPLETE
  plan.stages[1].status = STAGE_ACTIVE
  plan.current_stage_index = 1

  d = plan.to_dict()
  json_str = json.dumps(d)
  d2 = json.loads(json_str)
  plan2 = StructurePlan.from_dict(d2)

  assert plan2.goal == "Solve kinase by MR at 2.1A"
  assert len(plan2.stages) == 5
  assert plan2.current_stage_index == 1
  assert plan2.stages[0].status == STAGE_COMPLETE
  # Verify to_directives works after round-trip
  directives = plan2.to_directives()
  assert "phenix.phaser" in directives.get(
    "workflow_preferences", {}
  ).get("prefer_programs", [])


def test_plan_from_dict_tolerant():
  plan = StructurePlan.from_dict({
    "goal": "test",
    "extra_key": "ignored",
  })
  assert plan.goal == "test"
  assert plan.stages == []


def test_plan_from_dict_none():
  plan = StructurePlan.from_dict(None)
  assert plan.goal == ""
  assert plan.stages == []


# ── Real workflow scenarios ───────────────────────────

def test_mr_refine_workflow():
  """Simulate a successful MR + refine session.

  Walk through the plan as the gate evaluator would:
  start each stage, record cycles, mark complete,
  advance. Verify directives change at each stage.
  """
  plan = _make_mr_refine_plan()

  # Stage 1: data_assessment
  plan.mark_stage_started(1)
  plan.record_stage_cycle()
  plan.mark_phase_complete(1)
  d1 = plan.to_directives()
  assert "phenix.xtriage" in d1.get(
    "workflow_preferences", {}
  ).get("prefer_programs", [])

  # Advance to MR
  plan.advance()
  plan.mark_stage_started(2)
  plan.record_stage_cycle()
  d2 = plan.to_directives()
  assert "phenix.phaser" in d2.get(
    "workflow_preferences", {}
  ).get("prefer_programs", [])
  plan.mark_phase_complete(
    2, metrics={"tfz": 14.2}
  )

  # Advance to initial_refinement
  plan.advance()
  plan.mark_stage_started(3)
  d3 = plan.to_directives()
  assert "phenix.refine" in d3.get(
    "workflow_preferences", {}
  ).get("prefer_programs", [])
  for _ in range(3):
    plan.record_stage_cycle()
  plan.mark_phase_complete(
    5, metrics={"r_free": 0.34}
  )

  # Advance to model_rebuilding
  plan.advance()
  plan.mark_stage_started(6)
  d4 = plan.to_directives()
  ps = d4.get("program_settings", {})
  assert ps.get("phenix.autobuild", {}).get(
    "rebuild_in_place"
  ) is False
  plan.record_stage_cycle()
  plan.mark_phase_complete(
    6, metrics={"r_free": 0.28}
  )

  # Advance to final_refinement
  plan.advance()
  plan.mark_stage_started(7)
  d5 = plan.to_directives()
  ps5 = d5.get("program_settings", {})
  assert ps5.get("phenix.refine", {}).get(
    "ordered_solvent"
  ) is True
  for _ in range(3):
    plan.record_stage_cycle()
  plan.mark_phase_complete(
    9, metrics={"r_free": 0.24}
  )

  # No more stages
  ok = plan.advance()
  assert ok is False
  assert plan.is_complete() is True

  # Verify stage history
  assert plan.stages[0].status == STAGE_COMPLETE
  assert plan.stages[1].result_metrics[
    "tfz"
  ] == 14.2
  assert plan.stages[4].result_metrics[
    "r_free"
  ] == 0.24

  # Survives full round-trip
  d = plan.to_dict()
  plan2 = StructurePlan.from_dict(
    json.loads(json.dumps(d))
  )
  assert plan2.is_complete() is True
  assert plan2.stages[4].result_metrics[
    "r_free"
  ] == 0.24


def test_retreat_and_retry_workflow():
  """Simulate a retreat scenario.

  MR succeeds, but refinement stalls at R-free 0.48.
  Gate evaluator retreats to MR. Second MR attempt
  gives better solution. Refinement succeeds.
  """
  plan = _make_mr_refine_plan()

  # Stage 1: data_assessment
  plan.mark_stage_started(1)
  plan.record_stage_cycle()
  plan.mark_phase_complete(1)
  plan.advance()

  # Stage 2: MR (first attempt)
  plan.mark_stage_started(2)
  plan.record_stage_cycle()
  plan.mark_phase_complete(
    2, metrics={"tfz": 7.0}
  )
  plan.advance()
  h1 = plan.compute_hash()

  # Stage 3: Refinement stalls
  plan.mark_stage_started(3)
  plan.record_stage_cycle()
  plan.record_stage_cycle()
  # Gate evaluator decides to retreat
  assert plan.current_stage().cycles_used == 2

  # Retreat to MR
  plan.retreat_to("molecular_replacement")
  h2 = plan.compute_hash()
  assert h2 != h1  # hash changed
  assert plan.retreat_count == 1
  assert plan.current_stage().id == (
    "molecular_replacement"
  )
  assert plan.stages[2].status == STAGE_PENDING

  # MR second attempt
  plan.mark_stage_started(5)
  plan.record_stage_cycle()
  plan.mark_phase_complete(
    5, metrics={"tfz": 14.0}
  )
  plan.advance()

  # Refinement succeeds this time
  plan.mark_stage_started(6)
  for _ in range(3):
    plan.record_stage_cycle()
  plan.mark_phase_complete(
    8, metrics={"r_free": 0.32}
  )
  assert plan.stages[2].status == STAGE_COMPLETE

  # Plan continues normally
  plan.advance()
  assert plan.current_stage().id == (
    "model_rebuilding"
  )


# ── Entry point ──────────────────────────────────────

if __name__ == "__main__":
  run_tests()
