"""
Plan Schema for the Goal-Directed Agent.

Defines the data structures for multi-stage strategy
plans: StageDef (one stage in a plan) and StructurePlan
(the full plan with navigation and serialization).

The plan is generated once at session start (Stage 2.3),
consumed by the gate evaluator each cycle (Stage 3),
and displayed by the explanation engine (Stage 5).

Entry points: StageDef, StructurePlan classes.

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import copy
import hashlib
import json
import logging

logger = logging.getLogger(__name__)


# ── Stage status constants ──────────────────────────

STAGE_PENDING = "pending"
STAGE_ACTIVE = "active"
STAGE_COMPLETE = "complete"
STAGE_SKIPPED = "skipped"
STAGE_FAILED = "failed"

_VALID_PHASE_STATUSES = frozenset([
  STAGE_PENDING, STAGE_ACTIVE, STAGE_COMPLETE,
  STAGE_SKIPPED, STAGE_FAILED,
])


# ── StageDef ────────────────────────────────────────

class StageDef(object):
  """One stage in a structure determination plan.

  Defines what programs to run, when to advance, when
  to retreat, and what success looks like. Each stage
  maps to a segment of the crystallographic workflow
  (e.g. "molecular replacement", "initial refinement").

  The reactive agent sees only the directives derived
  from the current stage — it does not know about the
  plan structure. The goal-directed layer communicates
  through directives only.
  """

  def __init__(self, id, programs=None, max_cycles=5,
               success_criteria=None,
               gate_conditions=None,
               fallbacks=None, skip_if="",
               directives=None, strategy=None,
               description="",
               provides=None, if_skipped=None):
    self.id = str(id)
    self.programs = (
      list(programs) if programs else []
    )
    self.max_cycles = max(1, int(max_cycles))
    # e.g. {"r_free": "<0.35", "tfz": ">8"}
    self.success_criteria = (
      dict(success_criteria)
      if success_criteria else {}
    )
    # e.g. [{"if": "r_free > 0.45 after 2 cycles",
    #         "action": "retreat_to molecular_replacement"}]
    self.gate_conditions = (
      list(gate_conditions)
      if gate_conditions else []
    )
    # e.g. [{"if": "tfz < 6", "action": "try_alphafold"}]
    self.fallbacks = (
      list(fallbacks) if fallbacks else []
    )
    # e.g. "r_free < 0.28"
    self.skip_if = str(skip_if) if skip_if else ""
    # Translated to reactive agent directives
    # e.g. {"prefer_programs": [...]}
    self.directives = (
      dict(directives) if directives else {}
    )
    # Program-specific settings
    # e.g. {"ordered_solvent": True,
    #        "rebuild_in_place": False}
    self.strategy = (
      dict(strategy) if strategy else {}
    )
    # Human-readable description
    self.description = str(description)
    # What this stage provides to downstream stages
    # e.g. ["resolution", "twinning_status"]
    self.provides = (
      list(provides) if provides else []
    )
    # Repair rules if this stage is skipped
    # e.g. {"resolution": "extract from refine log"}
    self.if_skipped = (
      dict(if_skipped) if if_skipped else {}
    )

    # --- Runtime state (updated by gate evaluator) ---
    self.status = STAGE_PENDING
    self.cycles_used = 0
    self.start_cycle = None    # cycle number when entered
    self.end_cycle = None      # cycle number when exited
    self.result_metrics = {}   # metrics at stage completion
    # v119.H15 Item 1: deviation metadata.  Populated when a
    # stage transitions due to a plan deviation (LLM ran a
    # program belonging to a later stage out of plan order).
    # See StructurePlan.advance().  These fields are optional
    # observability data — None means "no deviation recorded
    # for this stage."  Captured here (not in a side channel)
    # so they survive serialization and resume.
    self.failure_reason = None  # e.g. "abandoned_by_deviation"
    self.entered_via = None     # e.g. "deviation" | "advance" | None

  def to_dict(self):
    """Serialize to JSON-safe dict."""
    return {
      "id": self.id,
      "programs": list(self.programs),
      "max_cycles": self.max_cycles,
      "success_criteria": dict(self.success_criteria),
      "gate_conditions": [
        dict(g) for g in self.gate_conditions
      ],
      "fallbacks": [
        dict(f) for f in self.fallbacks
      ],
      "skip_if": self.skip_if,
      "directives": dict(self.directives),
      "strategy": dict(self.strategy),
      "description": self.description,
      "provides": list(self.provides),
      "if_skipped": dict(self.if_skipped),
      # Runtime state
      "status": self.status,
      "cycles_used": self.cycles_used,
      "start_cycle": self.start_cycle,
      "end_cycle": self.end_cycle,
      "result_metrics": dict(self.result_metrics),
      # v119.H15 Item 1: deviation metadata
      "failure_reason": self.failure_reason,
      "entered_via": self.entered_via,
    }

  @classmethod
  def from_dict(cls, d):
    """Deserialize from dict. Tolerant of missing keys."""
    if not isinstance(d, dict):
      return cls(id="unknown")
    stage = cls(
      id=d.get("id", "unknown"),
      programs=d.get("programs"),
      max_cycles=d.get("max_cycles", 5),
      success_criteria=d.get("success_criteria"),
      gate_conditions=d.get("gate_conditions"),
      fallbacks=d.get("fallbacks"),
      skip_if=d.get("skip_if", ""),
      directives=d.get("directives"),
      strategy=d.get("strategy"),
      description=d.get("description", ""),
      provides=d.get("provides"),
      if_skipped=d.get("if_skipped"),
    )
    # Restore runtime state
    status = d.get("status", STAGE_PENDING)
    if status in _VALID_PHASE_STATUSES:
      stage.status = status
    stage.cycles_used = int(
      d.get("cycles_used", 0)
    )
    stage.start_cycle = d.get("start_cycle")
    stage.end_cycle = d.get("end_cycle")
    stage.result_metrics = dict(
      d.get("result_metrics", {})
    )
    # v119.H15 Item 1: deviation metadata (optional —
    # absent in pre-H15 session JSONs).
    stage.failure_reason = d.get("failure_reason")
    stage.entered_via = d.get("entered_via")
    return stage

  def is_exhausted(self):
    """Check if this stage has used all its cycles.

    Returns:
      True if cycles_used >= max_cycles.
    """
    return self.cycles_used >= self.max_cycles

  def __repr__(self):
    return (
      "StageDef(id=%r, status=%r, programs=%r)"
      % (self.id, self.status, self.programs)
    )


# ── Program matching helper ────────────────────────

def _program_matches_phase(program_name, phase_programs):
  """Check if a program name matches any stage program.

  Exact match first, then variant match: a program
  is a variant if a stage program is a prefix of the
  program name followed by _ (e.g. phenix.autobuild
  matches phenix.autobuild_denmod).

  Does NOT match across program families:
  phenix.refine does NOT match phenix.real_space_refine.

  Args:
    program_name: str, e.g. "phenix.autobuild_denmod"
    phase_programs: list of str, e.g.
      ["phenix.autobuild", "phenix.refine"]

  Returns:
    bool.
  """
  if not program_name or not phase_programs:
    return True  # No filter
  # Exact match
  if program_name in phase_programs:
    return True
  # Variant match: stage program + "_" is a prefix
  for pp in phase_programs:
    if program_name.startswith(pp + "_"):
      return True
  return False


# ── Success-criteria checker (v119.H15 Item 1) ──────

def _criteria_met(stage, structure_model):
  """Check whether a stage's success_criteria are met.

  Used by StructurePlan.advance(force=False) to decide
  between STAGE_COMPLETE and STAGE_FAILED when advancing
  away from a stage prematurely (LLM ran a later-stage
  program and triggered catch-up).

  Returns True (treat as met / let stage advance to
  COMPLETE) when:
    - stage has no success_criteria (e.g. validation
      stages with empty criteria — they always complete
      cleanly when their program runs)
    - structure_model is None (can't evaluate; preserve
      legacy "advance anyway" behavior)
    - all criteria evaluate to True against the model

  Returns False (treat as failed / mark STAGE_FAILED)
  when at least one criterion evaluates to False.

  Implementation: delegates to GateEvaluator._check_success
  via a lazy import (avoids a hard module-level
  circular dependency between plan_schema and
  gate_evaluator).  If the import fails for any reason,
  returns True (graceful degradation — error toward
  legacy behavior, never panic).
  """
  if not stage.success_criteria:
    return True
  if structure_model is None:
    return True
  try:
    try:
      from libtbx.langchain.agent.gate_evaluator \
        import GateEvaluator
    except ImportError:
      from agent.gate_evaluator import GateEvaluator
    all_met, _details = GateEvaluator()._check_success(
      stage.success_criteria, structure_model
    )
    return all_met
  except Exception:
    # Defensive: never raise from this helper.  If
    # criteria evaluation is broken for any reason,
    # fall back to legacy behavior (treat as met).
    logger.debug(
      "_criteria_met fallback: criteria evaluation "
      "failed for stage %r",
      stage.id,
      exc_info=True,
    )
    return True


# ── StructurePlan ───────────────────────────────────

class StructurePlan(object):
  """Multi-stage strategy plan for structure determination.

  Generated once at session start from templates
  (Step 2.3). Consumed by the gate evaluator after
  each cycle to decide advance/continue/retreat.
  Persisted to session JSON for resume support.

  Navigation:
    current_stage() — get the active StageDef
    advance() — move to next stage (returns False at end)
    retreat_to(id) — go back to a named stage
    mark_phase_complete() — mark current as done

  Directive translation:
    to_directives() — convert current stage to the
      reactive agent's directive format

  Display:
    get_display_phases() — list of {id, status, ...}
      for the GUI progress indicator
  """

  def __init__(self, goal="", stages=None,
               current_stage_index=0,
               strategy_hash="",
               created_at_cycle=0,
               revised_at_cycle=None,
               revision_reason="",
               template_id=""):
    self.goal = str(goal)
    self.stages = (
      list(stages) if stages else []
    )
    self.current_stage_index = max(
      0, int(current_stage_index)
    )
    self.strategy_hash = str(strategy_hash)
    self.created_at_cycle = int(created_at_cycle)
    self.revised_at_cycle = (
      int(revised_at_cycle)
      if revised_at_cycle is not None else None
    )
    self.revision_reason = str(revision_reason)
    self.template_id = str(template_id)
    # Retreat tracking (anti-oscillation)
    self.retreat_count = 0
    self.last_retreat_cycle = None

  # ── Navigation ──────────────────────────────────

  def current_stage(self):
    """Get the active StageDef.

    Returns:
      StageDef or None if index out of range or
      plan has no stages.
    """
    if not self.stages:
      return None
    idx = self.current_stage_index
    if 0 <= idx < len(self.stages):
      return self.stages[idx]
    return None

  def advance(self, force=True, structure_model=None,
              via=None):
    """Move to next stage.

    Marks the current stage and activates the next
    pending stage.  How the current stage is marked
    depends on `force` and `structure_model`:

    - force=True (default, legacy behavior): current
      stage is marked STAGE_COMPLETE if it was
      STAGE_ACTIVE.  No criteria check.  This is what
      every pre-H15 caller of advance() expected.

    - force=False: current stage's success_criteria are
      checked against structure_model.  If criteria are
      met (or empty, or model unavailable), marked
      STAGE_COMPLETE.  If criteria are NOT met, marked
      STAGE_FAILED with failure_reason populated.

    The `via` argument optionally tags the newly-active
    stage's `entered_via` field for observability —
    typically set by the catch-up call site to "deviation"
    so downstream code (and humans reading session JSON)
    can tell which transitions were plan-driven vs.
    LLM-driven.

    Args:
      force: bool, default True.  When True, current
        stage always marked COMPLETE on advance.  When
        False, criteria are evaluated.
      structure_model: StructureModel or dict-like with
        get_metric().  Only consulted when force=False.
      via: optional str.  Tagged onto the next active
        stage's entered_via field.

    Returns:
      True if advanced to a new stage.
      False if already at the last stage (plan done).
    """
    if not self.stages:
      return False
    # Mark current stage's exit status
    curr = self.current_stage()
    if curr and curr.status == STAGE_ACTIVE:
      if force:
        # Legacy behavior — pre-H15 callers preserved
        curr.status = STAGE_COMPLETE
      else:
        # Strict semantics — check success criteria
        if _criteria_met(curr, structure_model):
          curr.status = STAGE_COMPLETE
        else:
          # v119.H15 Item 1: explicit failure rather
          # than silently marking complete.  The
          # failure_reason field is observability for
          # session inspection and any future plan-
          # repair logic.  result_metrics are PRESERVED
          # (per H15 plan question 1: observation vs
          # status are orthogonal).
          curr.status = STAGE_FAILED
          curr.failure_reason = "abandoned_by_deviation"
    # Find next pending stage
    for i in range(
      self.current_stage_index + 1, len(self.stages)
    ):
      if self.stages[i].status == STAGE_PENDING:
        self.current_stage_index = i
        self.stages[i].status = STAGE_ACTIVE
        if via:
          self.stages[i].entered_via = via
        return True
    # No more pending stages — plan is done
    self.current_stage_index = len(self.stages)
    return False

  def record_plan_deviation(self, session_data,
                            from_stage_id, to_stage_id,
                            program_name, cycle_number):
    """Record a structured plan-deviation event.

    Called by record_stage_cycle when the LLM ran a
    program that matched a later stage's programs list,
    triggering a catch-up.  Stored in session_data
    under "plan_deviations" as a list of dicts.

    This data captures the LLM-vs-plan divergence in
    a form that future architectural work (e.g., a
    proper deviation-aware repair pass) can drive
    decisions from.  In H15, it's observability-only.

    Args:
      session_data: dict to append the event into.
        Modified in place.  No-op if not a dict.
      from_stage_id: str, the stage the agent was IN.
      to_stage_id: str, the stage the agent CATCHES UP to.
      program_name: str, the LLM's chosen program.
      cycle_number: int, the cycle of the deviation.
    """
    if not isinstance(session_data, dict):
      return
    events = session_data.setdefault(
      "plan_deviations", []
    )
    events.append({
      "_v": 1,
      "cycle": int(cycle_number) if cycle_number
              is not None else None,
      "from_stage": str(from_stage_id),
      "to_stage": str(to_stage_id),
      "program": str(program_name),
    })

  def retreat_to(self, stage_id, cycle_number=None):
    """Go back to a named stage.

    Resets the target stage to active, and resets all
    stages AFTER the target to pending (they need to
    re-run with the new strategy). Skipped stages are
    preserved.

    The gate evaluator handles blacklisting the failed
    strategy before calling retreat_to.

    Args:
      stage_id: str, the id of the target stage.
      cycle_number: int or None. Current cycle, used
        for retreat cooldown tracking.

    Returns:
      True if retreat succeeded.
      False if stage_id not found.
    """
    target_idx = None
    for i, stage in enumerate(self.stages):
      if stage.id == stage_id:
        target_idx = i
        break
    if target_idx is None:
      return False

    # Reset target stage
    target = self.stages[target_idx]
    target.status = STAGE_ACTIVE
    target.cycles_used = 0
    target.start_cycle = None
    target.end_cycle = None
    target.result_metrics = {}

    # Reset all stages AFTER the target to pending.
    # They may need to re-run with the new strategy.
    # Skipped stages are preserved — they were skipped
    # for structural reasons independent of the failed
    # strategy.
    for i in range(target_idx + 1, len(self.stages)):
      p = self.stages[i]
      if p.status in (
        STAGE_COMPLETE, STAGE_FAILED, STAGE_ACTIVE
      ):
        p.status = STAGE_PENDING
        p.cycles_used = 0
        p.start_cycle = None
        p.end_cycle = None
        p.result_metrics = {}

    self.current_stage_index = target_idx
    self.retreat_count += 1
    if cycle_number is not None:
      self.last_retreat_cycle = int(cycle_number)
    return True

  def skip_stage(self, stage_id, reason=""):
    """Skip a stage (e.g. skip_if condition met).

    Args:
      stage_id: str, stage to skip.
      reason: str, why skipped.

    Returns:
      True if skipped, False if not found.
    """
    for stage in self.stages:
      if stage.id == stage_id:
        stage.status = STAGE_SKIPPED
        return True
    return False

  def skip_to_program(self, program_name):
    """Skip stages until reaching one that contains program.

    Marks all earlier stages as SKIPPED and sets
    current_stage_index to the target stage.

    Called at plan initialization time when the user
    has set start_with_program or after_program to a
    target that isn't in the first plan stage (e.g.
    "predict and stop" — target is predict_and_build
    but Stage 1 is xtriage).  Without this, the agent
    would try to run xtriage first and stop with an
    "after_program_not_available" error.

    Only skips FORWARD — if the target program is in
    the current stage or an earlier (already-passed)
    stage, no skipping occurs and 0 is returned.

    Args:
      program_name: str, e.g. "phenix.predict_and_build"

    Returns:
      int: number of stages skipped, or 0 if program
        not found or already in current/earlier stage.
    """
    # Find first stage containing this program
    target_idx = None
    for i, stage in enumerate(self.stages):
      if program_name in stage.programs:
        target_idx = i
        break

    if target_idx is None:
      return 0  # Program not in any stage

    if target_idx <= self.current_stage_index:
      return 0  # Already at or past target

    # Skip everything before the target
    skipped = 0
    for i in range(self.current_stage_index, target_idx):
      if self.stages[i].status in (
              STAGE_PENDING, STAGE_ACTIVE):
        self.stages[i].status = STAGE_SKIPPED
        skipped += 1
    self.current_stage_index = target_idx
    return skipped

  def mark_stage_started(self, cycle_number):
    """Mark the current stage as started at a cycle.

    Called by the gate evaluator when a stage begins
    execution.

    Args:
      cycle_number: int.
    """
    curr = self.current_stage()
    if curr:
      curr.status = STAGE_ACTIVE
      if curr.start_cycle is None:
        curr.start_cycle = int(cycle_number)

  def record_stage_cycle(self, program_name=None,
                         structure_model=None,
                         session_data=None,
                         cycle_number=None):
    """Increment the cycle counter for current stage.

    Counting rules:
    - program_name is None → always count (compat)
    - program is STOP → never count
    - program matches current stage → count
    - program matches NO stage → count (reactive
      program like pdbtools running during stage)
    - program matches a LATER stage → catch up to that
      stage (PLAN DEVIATION) and count there
    - program matches an EARLIER stage → skip
      (e.g. xtriage re-run during refinement)

    Matching includes variant programs: e.g.
    phenix.autobuild_denmod counts as autobuild.

    v119.H15 Item 1: when the agent runs a program
    matching a later stage (the catch-up path), the
    intermediate-stage advances now use STRICT semantics:
    each intermediate stage's success_criteria are
    checked.  Stages that meet criteria are marked
    STAGE_COMPLETE; stages that don't are marked
    STAGE_FAILED with `failure_reason="abandoned_by_deviation"`.
    This replaces the pre-H15 behavior of marking every
    intermediate stage COMPLETE regardless of criteria
    (which corrupted Tom's bromodomain session — see
    H15 plan for the full diagnosis).

    The catch-up itself still occurs.  H15 doesn't
    eliminate the "tail wagging the dog" architectural
    smell (that would require a proper deviation-halt
    + plan-repair pass); it just makes the catch-up
    HONEST about which stages it abandoned without
    meeting their goals.

    A deviation event is also recorded in session_data
    (if provided) under "plan_deviations" — a structured
    log that future architectural work can drive
    decisions from.

    Args:
      program_name: str or None.
      structure_model: optional, used to evaluate
        success_criteria on catch-up intermediate stages.
        When None, intermediate stages fall back to
        legacy mark-COMPLETE behavior (graceful
        degradation when ai_agent.py is pre-H15).
      session_data: optional dict to record deviation
        events into.  Modified in place.
      cycle_number: optional, tagged onto deviation
        events for forensics.
    """
    curr = self.current_stage()
    if curr:
      if program_name is not None:
        # STOP never counts
        if program_name.upper() == "STOP":
          return
        if curr.programs:
          # Does it match THIS stage?
          if _program_matches_phase(
            program_name, curr.programs
          ):
            curr.cycles_used += 1
            return
          # Does it match a LATER stage?  PLAN DEVIATION.
          curr_idx = self.current_stage_index
          for i in range(
            curr_idx + 1, len(self.stages)
          ):
            s = self.stages[i]
            if (s.programs
                and _program_matches_phase(
                  program_name, s.programs)):
              logger.info(
                "Plan deviation: program %s matches "
                "stage '%s' (index %d), advancing "
                "from '%s' (index %d) with strict "
                "criteria check"
                % (program_name, s.id, i,
                   curr.id, curr_idx))
              # v119.H15 Item 1: structured deviation
              # event for session JSON
              self.record_plan_deviation(
                session_data,
                from_stage_id=curr.id,
                to_stage_id=s.id,
                program_name=program_name,
                cycle_number=cycle_number,
              )
              # Walk forward with strict semantics.  Each
              # intermediate stage's criteria are checked;
              # stages that don't meet criteria → FAILED,
              # not COMPLETE.  The first stage we advance
              # INTO (the new active one) gets entered_via
              # tagged so the gate evaluator and resume
              # logic can see this was a deviation.
              first_advance = True
              while (self.current_stage_index < i
                     and self.advance(
                       force=False,
                       structure_model=structure_model,
                       via=("deviation"
                            if first_advance else None),
                     )):
                first_advance = False
              # Count for the new current stage if it
              # matches the program (defensive against
              # overshoot e.g. if target was SKIPPED).
              new_curr = self.current_stage()
              if (new_curr
                  and new_curr.programs
                  and _program_matches_phase(
                    program_name, new_curr.programs)):
                new_curr.cycles_used += 1
              return
          # Does it match an EARLIER stage? Skip.
          for i in range(0, curr_idx):
            s = self.stages[i]
            if (s.programs
                and _program_matches_phase(
                  program_name, s.programs)):
              return
          # Matches no stage → reactive, count
          curr.cycles_used += 1
        else:
          curr.cycles_used += 1
      else:
        curr.cycles_used += 1

  def mark_phase_complete(self, cycle_number,
                          metrics=None):
    """Mark the current stage as complete.

    Args:
      cycle_number: int.
      metrics: dict of final metrics for this stage.
    """
    curr = self.current_stage()
    if curr:
      curr.status = STAGE_COMPLETE
      curr.end_cycle = int(cycle_number)
      if metrics:
        curr.result_metrics = dict(metrics)

  def is_complete(self):
    """Check if all stages are done (complete/skipped).

    Returns:
      True if no pending or active stages remain.

    Note: returns True even when one or more stages are in
    STAGE_FAILED state.  Callers that need to distinguish
    "all clean" from "exhausted with failures" should use
    get_failed_stage_ids() in conjunction with this check.
    """
    for stage in self.stages:
      if stage.status in (
        STAGE_PENDING, STAGE_ACTIVE
      ):
        return False
    return len(self.stages) > 0

  def get_failed_stage_ids(self):
    """List the IDs of stages currently in STAGE_FAILED.

    v119.H15 H1 fix: enables the gate evaluator to emit a
    stop reason that distinguishes "plan ran to completion"
    from "plan ran out of stages with one or more failed."
    Without this, the gate's stop event always said "all
    stages complete" even when a stage had status FAILED
    (e.g., abandoned via H15 Item 1's strict catch-up).

    Pure query — no mutation.

    Returns:
      list of str.  Empty list when no stages are FAILED.
    """
    return [
      s.id for s in self.stages
      if s.status == STAGE_FAILED
    ]

  def get_stage_by_id(self, stage_id):
    """Look up a stage by ID.

    Returns:
      StageDef or None.
    """
    for stage in self.stages:
      if stage.id == stage_id:
        return stage
    return None

  def get_previous_phase(self):
    """Get the stage before the current one.

    Used by retreat logic to determine the one-step-
    back target.

    Returns:
      StageDef or None.
    """
    idx = self.current_stage_index - 1
    if 0 <= idx < len(self.stages):
      return self.stages[idx]
    return None

  def can_retreat(self, cycle_number,
                  max_retreats=2, cooldown=2):
    """Check if retreat is allowed right now.

    Enforces anti-oscillation safeguards from the
    gate evaluator spec:
    1. Retreat counter: max N retreats total
    2. Retreat cooldown: at least N cycles since last

    Args:
      cycle_number: int, current cycle.
      max_retreats: int, maximum total retreats.
      cooldown: int, minimum cycles between retreats.

    Returns:
      (allowed: bool, reason: str).
      reason explains why retreat is blocked.
    """
    if self.retreat_count >= max_retreats:
      return (False,
        "max retreats reached (%d/%d)"
        % (self.retreat_count, max_retreats))
    if self.last_retreat_cycle is not None:
      elapsed = cycle_number - self.last_retreat_cycle
      if elapsed < cooldown:
        return (False,
          "retreat cooldown (%d/%d cycles)"
          % (elapsed, cooldown))
    return (True, "")

  # ── Directive translation ───────────────────────

  def to_directives(self):
    """Convert current stage to reactive agent directives.

    Produces a dict compatible with the existing
    directives system consumed by ai_agent.py and
    the PLAN/BUILD nodes.

    Returns:
      dict with keys:
        workflow_preferences: {prefer_programs, ...}
        stop_conditions: {after_program, ...}
        program_settings: {program: {key: val}, ...}

    Returns empty dict if no current stage.
    """
    curr = self.current_stage()
    if curr is None:
      return {}

    result = {}

    # --- workflow_preferences ---
    wf = {}
    if curr.programs:
      wf["prefer_programs"] = list(curr.programs)
    # Merge any explicit directives from the stage
    phase_wf = curr.directives.get(
      "workflow_preferences", {}
    )
    if phase_wf:
      wf.update(phase_wf)
    if wf:
      result["workflow_preferences"] = wf

    # --- stop_conditions ---
    sc = {}
    # If the stage has a single program, set
    # after_program so the reactive agent knows
    # what this stage is working toward.
    if len(curr.programs) == 1:
      sc["after_program"] = curr.programs[0]
    phase_sc = curr.directives.get(
      "stop_conditions", {}
    )
    if phase_sc:
      sc.update(phase_sc)
    if sc:
      result["stop_conditions"] = sc

    # --- program_settings ---
    ps = {}
    if curr.strategy:
      # Apply strategy settings to each program
      for prog in curr.programs:
        ps[prog] = dict(curr.strategy)
    phase_ps = curr.directives.get(
      "program_settings", {}
    )
    if phase_ps:
      for prog, settings in phase_ps.items():
        if prog not in ps:
          ps[prog] = {}
        ps[prog].update(settings)
    if ps:
      result["program_settings"] = ps

    return result

  # ── Hash computation ────────────────────────────

  def compute_hash(self):
    """Fingerprint of current stage + directives.

    Used to detect plan changes and trigger
    advice_changed. Changes when:
    - Current stage index changes
    - Stage directives change
    - Stages are added/removed

    Returns:
      str (hex digest, 12 chars).
    """
    try:
      data = {
        "current_stage_index": (
          self.current_stage_index
        ),
        "phase_ids": [
          p.id for p in self.stages
        ],
        "current_directives": self.to_directives(),
      }
      raw = json.dumps(
        data, sort_keys=True, default=str
      )
      return hashlib.md5(
        raw.encode("utf-8")
      ).hexdigest()[:12]
    except Exception:
      return ""

  # ── Display ─────────────────────────────────────

  def get_display_phases(self):
    """Stage list for GUI display.

    Returns:
      list of {id, description, programs, status,
        cycles_used, max_cycles, success_criteria}.
      Status is one of: pending, active, complete,
        skipped, failed.
    """
    result = []
    for stage in self.stages:
      result.append({
        "id": stage.id,
        "description": stage.description,
        "programs": list(stage.programs),
        "status": stage.status,
        "cycles_used": stage.cycles_used,
        "max_cycles": stage.max_cycles,
        "success_criteria": dict(
          stage.success_criteria
        ),
      })
    return result

  def format_plan_header(self, stop_after=None):
    """Format the plan as a compact text header.

    Used for the progress panel fixed header and
    for logging at session start.

    Args:
      stop_after: str or None.  If set (a program name
        like "phenix.resolve_cryo_em"), only display
        stages up to and including the one that contains
        this program.  Adds a "will stop here" indicator.

    Display rules:
      - Stages with STAGE_SKIPPED status are HIDDEN
        (e.g. prerequisites skipped via skip_to_program
        because the user asked to start later in the
        plan).  Showing "—" for skipped stages adds
        noise when the user didn't ask for them.
      - The first displayed stage gets a
        "▸ Starting here" annotation when any stages
        were skipped (explains non-1 stage numbering).
      - The last displayed stage gets "▸ Will stop here"
        when stop_after truncation is active.
      - Both annotations can appear on the same stage
        (e.g. "predict and stop" → Stage 2 is both
        start and stop).

    Returns:
      str, multi-line text.
    """
    lines = []
    bar = "=" * 50

    # If stop_after is set, find the last stage to show
    _stop_stage_idx = None
    if stop_after:
      for i, stage in enumerate(self.stages):
        if stop_after in stage.programs:
          _stop_stage_idx = i
          break

    # Where does "Starting here" point?  At the stage the
    # agent is currently about to work on — current_stage_index.
    # NOT at the first displayed stage, which on a resumed
    # session can be a COMPLETE stage from a previous run.
    # We only show the annotation when SKIPPED stages exist
    # before the current stage (this is what motivates the
    # annotation — it explains why the visible numbering
    # doesn't start at Stage 1).
    _start_here_idx = None
    if (self.current_stage_index is not None
            and 0 <= self.current_stage_index
                  < len(self.stages)):
      _curr = self.current_stage_index
      # Only annotate if the current stage will actually be
      # displayed (not past the stop point, not SKIPPED).
      _curr_displayable = (
        (_stop_stage_idx is None
         or _curr <= _stop_stage_idx)
        and self.stages[_curr].status != STAGE_SKIPPED)
      # Only annotate if there's at least one SKIPPED stage
      # before the current one — that's what the annotation
      # is for.
      _any_skipped_before = any(
        self.stages[i].status == STAGE_SKIPPED
        for i in range(_curr))
      if _curr_displayable and _any_skipped_before:
        _start_here_idx = _curr

    lines.append(bar)
    lines.append(
      " STRATEGY PLAN: %s" % self.goal
    )
    lines.append(bar)
    for i, stage in enumerate(self.stages):
      # Display truncation: skip stages after stop
      if _stop_stage_idx is not None and i > _stop_stage_idx:
        break
      # Hide skipped stages — user didn't ask for them,
      # and showing "skipped" adds noise to the plan
      # display.
      if stage.status == STAGE_SKIPPED:
        continue
      # Status indicator
      if stage.status == STAGE_COMPLETE:
        indicator = "✓"
      elif stage.status == STAGE_ACTIVE:
        indicator = "●"
      elif stage.status == STAGE_FAILED:
        indicator = "✗"
      else:
        indicator = "○"
      # Stage description (human-readable) or
      # formatted id as fallback
      phase_label = (
        stage.description
        or stage.id.replace("_", " ").title()
      )
      line = " %s Stage %d: %s" % (
        indicator, i + 1, phase_label
      )
      lines.append(line)
      # Success criteria (human-readable)
      if stage.success_criteria:
        goal_text = self._format_goal_readable(
          stage.success_criteria
        )
        if goal_text:
          lines.append(
            "          Goal: %s" % goal_text
          )
      # Combined start/stop annotation — both can fire
      # on the same stage (e.g. "predict and stop" where
      # the predict stage is both where the agent starts
      # and where it stops).
      _is_start_here = (
        _start_here_idx is not None
        and i == _start_here_idx)
      _is_stop_here = (
        _stop_stage_idx is not None
        and i == _stop_stage_idx)
      if _is_start_here and _is_stop_here:
        lines.append(
          "          ▸ Starting here — Will stop here"
        )
      elif _is_start_here:
        lines.append(
          "          ▸ Starting here"
        )
      elif _is_stop_here:
        lines.append(
          "          ▸ Will stop here"
        )
    lines.append(bar)
    return "\n".join(lines)

  @staticmethod
  def _format_goal_readable(criteria):
    """Format success criteria for non-expert display.

    Translates internal metric names and boolean
    flags into plain-language goal descriptions.

    Args:
      criteria: dict, e.g. {"r_free": "<0.30",
        "xtriage_completed": "true"}.

    Returns:
      str, human-readable goal text.
    """
    # Readable names for common metrics
    _METRIC_NAMES = {
      "r_free": "R-free",
      "r_work": "R-work",
      "model_map_cc": "map-model CC",
      "tfz": "MR score (TFZ)",
      "llg": "log-likelihood gain",
      "ligand_cc": "ligand fit score",
    }
    # Boolean completion flags
    _COMPLETION_NAMES = {
      "xtriage_completed":
        "Data quality analysis",
      "mtriage_completed":
        "Map quality analysis",
      "autosol_completed":
        "Experimental phasing",
      "autobuild_completed":
        "Model building",
      "predict_completed":
        "Model prediction",
      "phaser_completed":
        "Molecular replacement",
    }

    parts = []
    for key, value in criteria.items():
      if (str(value).lower() in ("true", "1")
          and key in _COMPLETION_NAMES):
        parts.append(
          "%s complete"
          % _COMPLETION_NAMES[key]
        )
      elif key in _COMPLETION_NAMES:
        parts.append(
          "%s complete"
          % _COMPLETION_NAMES[key]
        )
      elif key in _METRIC_NAMES:
        parts.append(
          "%s %s"
          % (_METRIC_NAMES[key], value)
        )
      else:
        parts.append(
          "%s %s"
          % (key.replace("_", " "), value)
        )

    return (
      ", ".join(parts) if parts else ""
    )

  # ── Serialization ───────────────────────────────

  def to_dict(self):
    """Serialize to JSON-safe dict."""
    try:
      return {
        "goal": self.goal,
        "stages": [
          p.to_dict() for p in self.stages
        ],
        "current_stage_index": (
          self.current_stage_index
        ),
        "strategy_hash": self.strategy_hash,
        "created_at_cycle": self.created_at_cycle,
        "revised_at_cycle": self.revised_at_cycle,
        "revision_reason": self.revision_reason,
        "template_id": self.template_id,
        "retreat_count": self.retreat_count,
        "last_retreat_cycle": (
          self.last_retreat_cycle
        ),
        "_version": 1,
      }
    except Exception:
      logger.debug(
        "StructurePlan.to_dict failed",
        exc_info=True,
      )
      return {}

  @classmethod
  def from_dict(cls, d):
    """Deserialize from session data.

    Tolerant of missing/extra keys.

    Args:
      d: dict from to_dict().

    Returns:
      StructurePlan instance.
    """
    if not isinstance(d, dict):
      return cls()
    try:
      # Accept both "stages" (v114.1+) and "phases"
      # (pre-v114.1 sessions) for backward compat.
      stages_data = (
        d.get("stages") or d.get("phases") or [])
      stages = [
        StageDef.from_dict(p)
        for p in stages_data
        if isinstance(p, dict)
      ]
      plan = cls(
        goal=d.get("goal", ""),
        stages=stages,
        current_stage_index=d.get(
          "current_stage_index",
          d.get("current_phase_index", 0)
        ),
        strategy_hash=d.get("strategy_hash", ""),
        created_at_cycle=d.get(
          "created_at_cycle", 0
        ),
        revised_at_cycle=d.get("revised_at_cycle"),
        revision_reason=d.get(
          "revision_reason", ""
        ),
        template_id=d.get("template_id", ""),
      )
      plan.retreat_count = int(
        d.get("retreat_count", 0)
      )
      plan.last_retreat_cycle = d.get(
        "last_retreat_cycle"
      )
      return plan
    except Exception:
      logger.debug(
        "StructurePlan.from_dict failed",
        exc_info=True,
      )
      return cls()


# ── Merge directives helper ─────────────────────────

def merge_directives(plan_directives,
                     user_directives):
  """Merge plan-generated and user-extracted directives.

  User directives take priority. Returns a merged dict
  and a list of conflict warnings.

  Args:
    plan_directives: dict from StructurePlan.to_directives()
    user_directives: dict from directive extraction

  Returns:
    tuple of (merged_dict, conflict_warnings)
    where conflict_warnings is a list of strings.
  """
  if not plan_directives:
    return (dict(user_directives or {}), [])
  if not user_directives:
    return (dict(plan_directives), [])

  merged = copy.deepcopy(plan_directives)
  warnings = []

  ud = user_directives

  # --- workflow_preferences ---
  plan_wf = merged.get(
    "workflow_preferences", {}
  )
  user_wf = ud.get("workflow_preferences", {})
  if user_wf:
    # User skip_programs may conflict with plan
    user_skip = user_wf.get("skip_programs", [])
    plan_prefer = plan_wf.get(
      "prefer_programs", []
    )
    for prog in user_skip:
      if prog in plan_prefer:
        warnings.append(
          "[Plan Conflict] User skips %s, but "
          "the plan requires it for this stage."
          % prog
        )
    # Merge: user overrides plan
    if "workflow_preferences" not in merged:
      merged["workflow_preferences"] = {}
    merged["workflow_preferences"].update(user_wf)
    # But keep plan's prefer_programs if user
    # didn't specify their own
    if "prefer_programs" not in user_wf and \
        plan_prefer:
      merged["workflow_preferences"][
        "prefer_programs"
      ] = plan_prefer

  # --- stop_conditions ---
  user_sc = ud.get("stop_conditions", {})
  if user_sc:
    # User after_program overrides plan
    if "stop_conditions" not in merged:
      merged["stop_conditions"] = {}
    plan_ap = merged.get(
      "stop_conditions", {}
    ).get("after_program")
    user_ap = user_sc.get("after_program")
    if user_ap and plan_ap and user_ap != plan_ap:
      warnings.append(
        "[Plan Conflict] User after_program=%s "
        "overrides plan after_program=%s."
        % (user_ap, plan_ap)
      )
    merged["stop_conditions"].update(user_sc)

  # --- program_settings ---
  user_ps = ud.get("program_settings", {})
  if user_ps:
    if "program_settings" not in merged:
      merged["program_settings"] = {}
    for prog, settings in user_ps.items():
      if prog not in merged["program_settings"]:
        merged["program_settings"][prog] = {}
      # User settings override plan settings
      merged["program_settings"][prog].update(
        settings
      )

  # Pass through any other user directive keys
  for key in ud:
    if key not in (
      "workflow_preferences",
      "stop_conditions",
      "program_settings",
    ):
      merged[key] = ud[key]

  return (merged, warnings)
