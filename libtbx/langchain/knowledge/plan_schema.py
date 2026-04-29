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

  def advance(self):
    """Move to next stage.

    Marks the current stage as complete and activates
    the next pending stage. Skips stages whose skip_if
    condition is noted (actual evaluation is done by
    the gate evaluator, which calls skip_stage()).

    Returns:
      True if advanced to a new stage.
      False if already at the last stage (plan done).
    """
    if not self.stages:
      return False
    # Mark current as complete
    curr = self.current_stage()
    if curr and curr.status == STAGE_ACTIVE:
      curr.status = STAGE_COMPLETE
    # Find next pending stage
    for i in range(
      self.current_stage_index + 1, len(self.stages)
    ):
      if self.stages[i].status == STAGE_PENDING:
        self.current_stage_index = i
        self.stages[i].status = STAGE_ACTIVE
        return True
    # No more pending stages — plan is done
    self.current_stage_index = len(self.stages)
    return False

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

  def record_stage_cycle(self, program_name=None):
    """Increment the cycle counter for current stage.

    Counting rules:
    - program_name is None → always count (compat)
    - program is STOP → never count
    - program matches current stage → count
    - program matches NO stage → count (reactive
      program like pdbtools running during stage)
    - program matches a LATER stage → advance the
      plan to that stage and count there (the agent
      ran ahead of the plan tracker)
    - program matches an EARLIER stage → skip
      (e.g. xtriage re-run during refinement)

    Matching includes variant programs: e.g.
    phenix.autobuild_denmod counts as autobuild.

    Args:
      program_name: str or None.
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
          # Does it match a LATER stage? If so, the
          # agent ran ahead — advance the plan to
          # catch up.  Intermediate stages are marked
          # complete (the agent effectively handled
          # them).
          curr_idx = self.current_stage_index
          for i in range(
            curr_idx + 1, len(self.stages)
          ):
            s = self.stages[i]
            if (s.programs
                and _program_matches_phase(
                  program_name, s.programs)):
              # Advance through all intermediate
              # stages to reach this one
              logger.info(
                "Plan catch-up: program %s matches "
                "stage '%s' (index %d), advancing "
                "from '%s' (index %d)"
                % (program_name, s.id, i,
                   curr.id, curr_idx))
              while (self.current_stage_index < i
                     and self.advance()):
                pass
              # Now count for the new current stage,
              # but only if it actually matches the
              # program.  If advance() overshot (e.g.
              # the target stage was SKIPPED), don't
              # mis-count on an unrelated stage.
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
    """
    for stage in self.stages:
      if stage.status in (
        STAGE_PENDING, STAGE_ACTIVE
      ):
        return False
    return len(self.stages) > 0

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

  def format_plan_header(self):
    """Format the plan as a compact text header.

    Used for the progress panel fixed header and
    for logging at session start.

    Returns:
      str, multi-line text.
    """
    lines = []
    bar = "=" * 50
    lines.append(bar)
    lines.append(
      " STRATEGY PLAN: %s" % self.goal
    )
    lines.append(bar)
    for i, stage in enumerate(self.stages):
      # Status indicator
      if stage.status == STAGE_COMPLETE:
        indicator = "✓"
      elif stage.status == STAGE_ACTIVE:
        indicator = "●"
      elif stage.status == STAGE_SKIPPED:
        indicator = "—"
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
