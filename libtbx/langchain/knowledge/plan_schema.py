"""
Plan Schema for the Goal-Directed Agent.

Defines the data structures for multi-phase strategy
plans: PhaseDef (one phase in a plan) and StructurePlan
(the full plan with navigation and serialization).

The plan is generated once at session start (Phase 2.3),
consumed by the gate evaluator each cycle (Phase 3),
and displayed by the explanation engine (Phase 5).

Entry points: PhaseDef, StructurePlan classes.

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import copy
import hashlib
import json
import logging

logger = logging.getLogger(__name__)


# ── Phase status constants ──────────────────────────

PHASE_PENDING = "pending"
PHASE_ACTIVE = "active"
PHASE_COMPLETE = "complete"
PHASE_SKIPPED = "skipped"
PHASE_FAILED = "failed"

_VALID_PHASE_STATUSES = frozenset([
  PHASE_PENDING, PHASE_ACTIVE, PHASE_COMPLETE,
  PHASE_SKIPPED, PHASE_FAILED,
])


# ── PhaseDef ────────────────────────────────────────

class PhaseDef(object):
  """One phase in a structure determination plan.

  Defines what programs to run, when to advance, when
  to retreat, and what success looks like. Each phase
  maps to a segment of the crystallographic workflow
  (e.g. "molecular replacement", "initial refinement").

  The reactive agent sees only the directives derived
  from the current phase — it does not know about the
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
    # What this phase provides to downstream phases
    # e.g. ["resolution", "twinning_status"]
    self.provides = (
      list(provides) if provides else []
    )
    # Repair rules if this phase is skipped
    # e.g. {"resolution": "extract from refine log"}
    self.if_skipped = (
      dict(if_skipped) if if_skipped else {}
    )

    # --- Runtime state (updated by gate evaluator) ---
    self.status = PHASE_PENDING
    self.cycles_used = 0
    self.start_cycle = None    # cycle number when entered
    self.end_cycle = None      # cycle number when exited
    self.result_metrics = {}   # metrics at phase completion

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
    phase = cls(
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
    status = d.get("status", PHASE_PENDING)
    if status in _VALID_PHASE_STATUSES:
      phase.status = status
    phase.cycles_used = int(
      d.get("cycles_used", 0)
    )
    phase.start_cycle = d.get("start_cycle")
    phase.end_cycle = d.get("end_cycle")
    phase.result_metrics = dict(
      d.get("result_metrics", {})
    )
    return phase

  def is_exhausted(self):
    """Check if this phase has used all its cycles.

    Returns:
      True if cycles_used >= max_cycles.
    """
    return self.cycles_used >= self.max_cycles

  def __repr__(self):
    return (
      "PhaseDef(id=%r, status=%r, programs=%r)"
      % (self.id, self.status, self.programs)
    )


# ── Program matching helper ────────────────────────

def _program_matches_phase(program_name, phase_programs):
  """Check if a program name matches any phase program.

  Exact match first, then variant match: a program
  is a variant if a phase program is a prefix of the
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
  # Variant match: phase program + "_" is a prefix
  for pp in phase_programs:
    if program_name.startswith(pp + "_"):
      return True
  return False


# ── StructurePlan ───────────────────────────────────

class StructurePlan(object):
  """Multi-phase strategy plan for structure determination.

  Generated once at session start from templates
  (Step 2.3). Consumed by the gate evaluator after
  each cycle to decide advance/continue/retreat.
  Persisted to session JSON for resume support.

  Navigation:
    current_phase() — get the active PhaseDef
    advance() — move to next phase (returns False at end)
    retreat_to(id) — go back to a named phase
    mark_phase_complete() — mark current as done

  Directive translation:
    to_directives() — convert current phase to the
      reactive agent's directive format

  Display:
    get_display_phases() — list of {id, status, ...}
      for the GUI progress indicator
  """

  def __init__(self, goal="", phases=None,
               current_phase_index=0,
               strategy_hash="",
               created_at_cycle=0,
               revised_at_cycle=None,
               revision_reason="",
               template_id=""):
    self.goal = str(goal)
    self.phases = (
      list(phases) if phases else []
    )
    self.current_phase_index = max(
      0, int(current_phase_index)
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

  def current_phase(self):
    """Get the active PhaseDef.

    Returns:
      PhaseDef or None if index out of range or
      plan has no phases.
    """
    if not self.phases:
      return None
    idx = self.current_phase_index
    if 0 <= idx < len(self.phases):
      return self.phases[idx]
    return None

  def advance(self):
    """Move to next phase.

    Marks the current phase as complete and activates
    the next pending phase. Skips phases whose skip_if
    condition is noted (actual evaluation is done by
    the gate evaluator, which calls skip_phase()).

    Returns:
      True if advanced to a new phase.
      False if already at the last phase (plan done).
    """
    if not self.phases:
      return False
    # Mark current as complete
    curr = self.current_phase()
    if curr and curr.status == PHASE_ACTIVE:
      curr.status = PHASE_COMPLETE
    # Find next pending phase
    for i in range(
      self.current_phase_index + 1, len(self.phases)
    ):
      if self.phases[i].status == PHASE_PENDING:
        self.current_phase_index = i
        self.phases[i].status = PHASE_ACTIVE
        return True
    # No more pending phases — plan is done
    self.current_phase_index = len(self.phases)
    return False

  def retreat_to(self, phase_id, cycle_number=None):
    """Go back to a named phase.

    Resets the target phase to active, and resets all
    phases AFTER the target to pending (they need to
    re-run with the new strategy). Skipped phases are
    preserved.

    The gate evaluator handles blacklisting the failed
    strategy before calling retreat_to.

    Args:
      phase_id: str, the id of the target phase.
      cycle_number: int or None. Current cycle, used
        for retreat cooldown tracking.

    Returns:
      True if retreat succeeded.
      False if phase_id not found.
    """
    target_idx = None
    for i, phase in enumerate(self.phases):
      if phase.id == phase_id:
        target_idx = i
        break
    if target_idx is None:
      return False

    # Reset target phase
    target = self.phases[target_idx]
    target.status = PHASE_ACTIVE
    target.cycles_used = 0
    target.start_cycle = None
    target.end_cycle = None
    target.result_metrics = {}

    # Reset all phases AFTER the target to pending.
    # They may need to re-run with the new strategy.
    # Skipped phases are preserved — they were skipped
    # for structural reasons independent of the failed
    # strategy.
    for i in range(target_idx + 1, len(self.phases)):
      p = self.phases[i]
      if p.status in (
        PHASE_COMPLETE, PHASE_FAILED, PHASE_ACTIVE
      ):
        p.status = PHASE_PENDING
        p.cycles_used = 0
        p.start_cycle = None
        p.end_cycle = None
        p.result_metrics = {}

    self.current_phase_index = target_idx
    self.retreat_count += 1
    if cycle_number is not None:
      self.last_retreat_cycle = int(cycle_number)
    return True

  def skip_phase(self, phase_id, reason=""):
    """Skip a phase (e.g. skip_if condition met).

    Args:
      phase_id: str, phase to skip.
      reason: str, why skipped.

    Returns:
      True if skipped, False if not found.
    """
    for phase in self.phases:
      if phase.id == phase_id:
        phase.status = PHASE_SKIPPED
        return True
    return False

  def mark_phase_started(self, cycle_number):
    """Mark the current phase as started at a cycle.

    Called by the gate evaluator when a phase begins
    execution.

    Args:
      cycle_number: int.
    """
    curr = self.current_phase()
    if curr:
      curr.status = PHASE_ACTIVE
      if curr.start_cycle is None:
        curr.start_cycle = int(cycle_number)

  def record_phase_cycle(self, program_name=None):
    """Increment the cycle counter for current phase.

    Counting rules:
    - program_name is None → always count (compat)
    - program is STOP → never count
    - program matches current phase → count
    - program matches NO phase → count (reactive
      program like pdbtools running during phase)
    - program matches a DIFFERENT phase → skip
      (e.g. ligandfit during refinement phase)

    Matching includes variant programs: e.g.
    phenix.autobuild_denmod counts as autobuild.

    Args:
      program_name: str or None.
    """
    curr = self.current_phase()
    if curr:
      if program_name is not None:
        # STOP never counts
        if program_name.upper() == "STOP":
          return
        if curr.programs:
          # Does it match THIS phase?
          if _program_matches_phase(
            program_name, curr.programs
          ):
            curr.cycles_used += 1
            return
          # Does it match ANY other phase?
          for p in self.phases:
            if p is curr:
              continue
            if (p.programs
                and _program_matches_phase(
                  program_name, p.programs)):
              return
          # Matches no phase → reactive, count
          curr.cycles_used += 1
        else:
          curr.cycles_used += 1
      else:
        curr.cycles_used += 1

  def mark_phase_complete(self, cycle_number,
                          metrics=None):
    """Mark the current phase as complete.

    Args:
      cycle_number: int.
      metrics: dict of final metrics for this phase.
    """
    curr = self.current_phase()
    if curr:
      curr.status = PHASE_COMPLETE
      curr.end_cycle = int(cycle_number)
      if metrics:
        curr.result_metrics = dict(metrics)

  def is_complete(self):
    """Check if all phases are done (complete/skipped).

    Returns:
      True if no pending or active phases remain.
    """
    for phase in self.phases:
      if phase.status in (
        PHASE_PENDING, PHASE_ACTIVE
      ):
        return False
    return len(self.phases) > 0

  def get_phase_by_id(self, phase_id):
    """Look up a phase by ID.

    Returns:
      PhaseDef or None.
    """
    for phase in self.phases:
      if phase.id == phase_id:
        return phase
    return None

  def get_previous_phase(self):
    """Get the phase before the current one.

    Used by retreat logic to determine the one-step-
    back target.

    Returns:
      PhaseDef or None.
    """
    idx = self.current_phase_index - 1
    if 0 <= idx < len(self.phases):
      return self.phases[idx]
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
    """Convert current phase to reactive agent directives.

    Produces a dict compatible with the existing
    directives system consumed by ai_agent.py and
    the PLAN/BUILD nodes.

    Returns:
      dict with keys:
        workflow_preferences: {prefer_programs, ...}
        stop_conditions: {after_program, ...}
        program_settings: {program: {key: val}, ...}

    Returns empty dict if no current phase.
    """
    curr = self.current_phase()
    if curr is None:
      return {}

    result = {}

    # --- workflow_preferences ---
    wf = {}
    if curr.programs:
      wf["prefer_programs"] = list(curr.programs)
    # Merge any explicit directives from the phase
    phase_wf = curr.directives.get(
      "workflow_preferences", {}
    )
    if phase_wf:
      wf.update(phase_wf)
    if wf:
      result["workflow_preferences"] = wf

    # --- stop_conditions ---
    sc = {}
    # If the phase has a single program, set
    # after_program so the reactive agent knows
    # what this phase is working toward.
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
    """Fingerprint of current phase + directives.

    Used to detect plan changes and trigger
    advice_changed. Changes when:
    - Current phase index changes
    - Phase directives change
    - Phases are added/removed

    Returns:
      str (hex digest, 12 chars).
    """
    try:
      data = {
        "current_phase_index": (
          self.current_phase_index
        ),
        "phase_ids": [
          p.id for p in self.phases
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
    """Phase list for GUI display.

    Returns:
      list of {id, description, programs, status,
        cycles_used, max_cycles, success_criteria}.
      Status is one of: pending, active, complete,
        skipped, failed.
    """
    result = []
    for phase in self.phases:
      result.append({
        "id": phase.id,
        "description": phase.description,
        "programs": list(phase.programs),
        "status": phase.status,
        "cycles_used": phase.cycles_used,
        "max_cycles": phase.max_cycles,
        "success_criteria": dict(
          phase.success_criteria
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
    for i, phase in enumerate(self.phases):
      # Status indicator
      if phase.status == PHASE_COMPLETE:
        indicator = "✓"
      elif phase.status == PHASE_ACTIVE:
        indicator = "●"
      elif phase.status == PHASE_SKIPPED:
        indicator = "—"
      elif phase.status == PHASE_FAILED:
        indicator = "✗"
      else:
        indicator = "○"
      # Phase description (human-readable) or
      # formatted id as fallback
      phase_label = (
        phase.description
        or phase.id.replace("_", " ").title()
      )
      line = " %s Phase %d: %s" % (
        indicator, i + 1, phase_label
      )
      lines.append(line)
      # Success criteria (human-readable)
      if phase.success_criteria:
        goal_text = self._format_goal_readable(
          phase.success_criteria
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
        "phases": [
          p.to_dict() for p in self.phases
        ],
        "current_phase_index": (
          self.current_phase_index
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
      phases_data = d.get("phases", [])
      phases = [
        PhaseDef.from_dict(p)
        for p in phases_data
        if isinstance(p, dict)
      ]
      plan = cls(
        goal=d.get("goal", ""),
        phases=phases,
        current_phase_index=d.get(
          "current_phase_index", 0
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
          "the plan requires it for this phase."
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
