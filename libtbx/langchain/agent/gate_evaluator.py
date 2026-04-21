"""
Gate Evaluator for the Goal-Directed Agent.

Evaluates stage progress against plan criteria after
each cycle. Decides advance/continue/retreat/stop using
purely deterministic logic — no LLM involvement.

Uses the Structure Model (ground truth) and
ValidationHistory for metric extraction.

Entry point: GateEvaluator.evaluate(plan, structure_model,
  validation_history, cycle_number)

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import logging
import re

logger = logging.getLogger(__name__)


# ── GateResult ──────────────────────────────────────

class GateResult(object):
  """Result of a gate evaluation.

  Attributes:
    action: str. One of:
      "continue" — keep going in current stage
      "advance" — move to next stage (success)
      "skip" — skip current stage, move to next
      "retreat" — go back to an earlier stage
      "fallback" — try within-stage fallback
      "stop" — end the session
    reason: str. Human-readable explanation.
    new_phase_id: str or None. Target stage for
      advance/retreat.
    blacklist_entry: dict or None. Strategy to
      blacklist on retreat.
    details: list of dicts. Per-criterion results.
  """

  def __init__(self, action="continue", reason="",
               new_phase_id=None,
               blacklist_entry=None,
               details=None):
    self.action = str(action)
    self.reason = str(reason)
    self.new_phase_id = new_phase_id
    self.blacklist_entry = blacklist_entry
    self.details = details if details else []

  def __repr__(self):
    return (
      "GateResult(action=%r, reason=%r)"
      % (self.action, self.reason)
    )


# ── Criterion parsing ───────────────────────────────

_CRITERION_RE = re.compile(
  r'^([<>=!]+)\s*(-?[0-9.]+)$'
)

_BOOLEAN_VALUES = frozenset([
  "true", "false", "yes", "no",
])


def parse_criterion(criterion_str):
  """Parse a criterion string like "<0.35" or ">8".

  Args:
    criterion_str: str, e.g. "<0.35", ">8",
      "<=0.30", "true".

  Returns:
    (operator, value) where operator is str and
    value is float. For booleans, returns
    ("==", True/False).
    Returns (None, None) if unparseable.
  """
  s = str(criterion_str).strip()
  if not s:
    return (None, None)
  # Boolean check
  if s.lower() in _BOOLEAN_VALUES:
    return ("==", s.lower() in ("true", "yes"))
  # Numeric with operator
  m = _CRITERION_RE.match(s)
  if m:
    op = m.group(1)
    try:
      val = float(m.group(2))
      return (op, val)
    except ValueError:
      return (None, None)
  return (None, None)


def evaluate_criterion(actual, operator, target):
  """Evaluate a single criterion.

  Args:
    actual: float, int, or bool. Current value.
    operator: str ("<", ">", "<=", ">=", "==", "!=").
    target: float or bool.

  Returns:
    True if criterion is met.
  """
  if actual is None:
    return False
  if isinstance(target, bool):
    return bool(actual) == target
  try:
    actual_f = float(actual)
    target_f = float(target)
  except (ValueError, TypeError):
    return False
  if operator == "<":
    return actual_f < target_f
  elif operator == "<=":
    return actual_f <= target_f
  elif operator == ">":
    return actual_f > target_f
  elif operator == ">=":
    return actual_f >= target_f
  elif operator == "==":
    return abs(actual_f - target_f) < 1e-9
  elif operator == "!=":
    return abs(actual_f - target_f) >= 1e-9
  return False


def apply_hysteresis(target_value, operator,
                     buffer_fraction=0.015):
  """Apply hysteresis buffer to a threshold.

  Tightens the threshold slightly to prevent
  oscillation on noisy metrics.

  For "<" targets: tighten (0.35 → ~0.345)
  For ">" targets: tighten (8.0 → ~8.12)

  Args:
    target_value: float.
    operator: str.
    buffer_fraction: float (default 1.5%).

  Returns:
    float, adjusted threshold.
  """
  if isinstance(target_value, bool):
    return target_value
  try:
    tv = float(target_value)
  except (ValueError, TypeError):
    return target_value
  buffer = abs(tv) * buffer_fraction
  if operator in ("<", "<="):
    return tv - buffer
  elif operator in (">", ">="):
    return tv + buffer
  return tv


# ── Gate condition parsing ──────────────────────────

_GATE_AFTER_RE = re.compile(
  r'after\s+(\d+)\s+cycle'
)


def parse_gate_condition(condition_str):
  """Parse a gate condition string.

  Format: "metric op value [after N cycles]"
  Examples:
    "r_free > 0.45 after 2 cycles"
    "tfz < 5"
    "r_free 0.35-0.45 after 3 cycles" (range)

  Returns:
    dict with:
      metric: str
      operator: str
      value: float
      after_cycles: int or None
    Or None if unparseable.
  """
  s = str(condition_str).strip()
  if not s:
    return None

  # Extract "after N cycles"
  after_cycles = None
  m_after = _GATE_AFTER_RE.search(s)
  if m_after:
    after_cycles = int(m_after.group(1))
    s = s[:m_after.start()].strip()

  # Try range format "metric low-high"
  range_m = re.match(
    r'^(\w+)\s+([0-9.]+)-([0-9.]+)$', s
  )
  if range_m:
    metric = range_m.group(1)
    low = float(range_m.group(2))
    high = float(range_m.group(3))
    return {
      "metric": metric,
      "operator": "range",
      "value": (low, high),
      "after_cycles": after_cycles,
    }

  # Try "metric op value"
  op_m = re.match(
    r'^(\w+)\s*([<>=!]+)\s*(-?[0-9.]+)$', s
  )
  if op_m:
    return {
      "metric": op_m.group(1),
      "operator": op_m.group(2),
      "value": float(op_m.group(3)),
      "after_cycles": after_cycles,
    }

  return None


def evaluate_gate_condition(parsed, actual_value,
                            cycles_in_phase):
  """Evaluate a parsed gate condition.

  Args:
    parsed: dict from parse_gate_condition.
    actual_value: float or None.
    cycles_in_phase: int.

  Returns:
    True if the condition fires (gate triggered).
  """
  if parsed is None or actual_value is None:
    return False

  # Check after_cycles requirement
  after = parsed.get("after_cycles")
  if after is not None and cycles_in_phase < after:
    return False

  op = parsed["operator"]
  target = parsed["value"]

  if op == "range":
    low, high = target
    try:
      v = float(actual_value)
      return low <= v <= high
    except (ValueError, TypeError):
      return False

  return evaluate_criterion(actual_value, op, target)


# ── GateEvaluator ───────────────────────────────────

class GateEvaluator(object):
  """Evaluates stage progress against plan criteria.

  Called after each cycle. Uses the Structure Model
  (ground truth) to check success criteria without
  any LLM involvement. Purely deterministic.

  Usage:
    evaluator = GateEvaluator()
    result = evaluator.evaluate(
      plan, structure_model,
      validation_history, cycle_number)
  """

  def evaluate(self, plan, structure_model,
               validation_history, cycle_number):
    """Evaluate current stage progress.

    Order of checks:
    1. Is plan complete? → stop
    2. Is current stage None? → stop
    3. Does skip_if apply? → skip
    4. Are success criteria met? → advance
    5. Do gate conditions fire? → retreat/stop/fallback
    6. Do fallback conditions fire? → fallback
    7. Is stage exhausted (max_cycles)? → advance
    8. Otherwise → continue

    Gates are checked BEFORE exhaustion because gates
    may include critical stop/retreat conditions that
    must not be skipped at max_cycles.

    Args:
      plan: StructurePlan.
      structure_model: StructureModel (or dict with
        get_metric method, or None).
      validation_history: ValidationHistory (or None).
      cycle_number: int.

    Returns:
      GateResult.

    Never raises.
    """
    try:
      return self._evaluate_inner(
        plan, structure_model,
        validation_history, cycle_number,
      )
    except Exception:
      logger.debug(
        "GateEvaluator.evaluate failed",
        exc_info=True,
      )
      return GateResult(
        action="continue",
        reason="gate evaluation failed — continuing",
      )

  def _evaluate_inner(self, plan, structure_model,
                      validation_history,
                      cycle_number):
    """Inner evaluation. May raise."""
    if plan is None:
      return GateResult(
        action="continue",
        reason="no plan",
      )

    # Is plan already complete?
    if plan.is_complete():
      return GateResult(
        action="stop",
        reason="all stages complete",
      )

    stage = plan.current_stage()
    if stage is None:
      return GateResult(
        action="stop",
        reason="no active stage",
      )

    # --- Check skip_if (stage entry only) ---
    # skip_if is a pre-entry check: "don't enter this
    # stage if the condition is already met." Once the
    # stage has started (cycles_used > 0), skip_if no
    # longer applies — success criteria take over.
    if (stage.skip_if and structure_model
        and stage.cycles_used == 0):
      skip_result = self._check_skip(
        stage, structure_model
      )
      if skip_result is not None:
        return skip_result

    # --- Check success criteria ---
    if stage.success_criteria and structure_model:
      all_met, details = self._check_success(
        stage.success_criteria, structure_model
      )
      if all_met:
        return GateResult(
          action="advance",
          reason="success criteria met",
          details=details,
        )

    # --- Anomalous signal guard (v115.05) ──────────
    # After the first autosol attempt in
    # experimental_phasing, verify that anomalous
    # signal was actually usable.  Without this, the
    # agent wastes cycles running autosol on noise
    # when anomalous_measurability is negligible.
    #
    # Placed AFTER success criteria: if autosol
    # somehow succeeded, let the advance happen.
    # Uses cycles_used <= 1 (not == 0) because the
    # gate evaluates AFTER each cycle, so the first
    # evaluation of this stage has cycles_used == 1.
    if (stage.id == "experimental_phasing"
        and stage.cycles_used <= 1
        and structure_model):
      _anom_m = _get_metric_value(
        structure_model,
        "anomalous_measurability",
      )
      _has_anom = _get_metric_value(
        structure_model, "has_anomalous",
      )
      # Block only when has_anomalous is not explicitly
      # True.  When the user genuinely has anomalous data,
      # allow autosol even with weak measurability.
      if (_anom_m is not None
          and _has_anom is not True):
        try:
          _anom_f = float(_anom_m)
        except (ValueError, TypeError):
          _anom_f = None
        if _anom_f is not None and _anom_f < 0.05:
          return GateResult(
            action="stop",
            reason=(
              "Anomalous measurability %.3f "
              "< 0.05 — signal too weak for "
              "SAD phasing. Consider "
              "predict_and_build or MR "
              "instead." % _anom_f
            ),
          )

    # --- Check gate conditions ---
    # Gates must fire BEFORE exhaustion check,
    # because gates can include critical stop/retreat
    # conditions that shouldn't be skipped when
    # max_cycles is reached.
    if stage.gate_conditions and structure_model:
      gate_result = self._check_gates(
        plan, stage, structure_model,
        cycle_number,
      )
      if gate_result is not None:
        return gate_result

    # --- Check fallback conditions ---
    if stage.fallbacks and structure_model:
      fb_result = self._check_fallbacks(
        stage, structure_model,
      )
      if fb_result is not None:
        return fb_result

    # --- Check stage exhaustion ---
    if stage.is_exhausted():
      return self._handle_exhaustion(
        plan, stage, structure_model,
        cycle_number,
      )

    return GateResult(
      action="continue",
      reason="stage in progress",
    )

  # ── Hypothesis evaluation ───────────────────────

  def evaluate_hypotheses(self, structure_model,
                          cycle_number):
    """Evaluate active hypotheses and revalidate
    confirmed ones.

    Called after the main stage evaluation.
    Separate from evaluate() because hypothesis
    lifecycle is orthogonal to stage transitions.

    Args:
      structure_model: StructureModel.
      cycle_number: int.

    Returns:
      list of HypothesisResult. Empty if no
      hypotheses or imports unavailable.

    Never raises.
    """
    try:
      try:
        from libtbx.langchain.agent \
          .hypothesis_evaluator import (
            evaluate_hypotheses,
            revalidate_confirmed,
          )
      except ImportError:
        from agent.hypothesis_evaluator import (
          evaluate_hypotheses,
          revalidate_confirmed,
        )
      results = evaluate_hypotheses(
        structure_model, cycle_number,
      )
      reval = revalidate_confirmed(
        structure_model,
      )
      results.extend(reval)
      return results
    except ImportError:
      return []
    except Exception:
      logger.debug(
        "GateEvaluator.evaluate_hypotheses failed",
        exc_info=True,
      )
      return []

  # ── Success criteria ────────────────────────────

  def _check_success(self, criteria, structure_model):
    """Check if all success criteria are met.

    Args:
      criteria: dict, e.g. {"r_free": "<0.35"}.
      structure_model: StructureModel.

    Returns:
      (all_met: bool, details: list of dicts).
    """
    details = []
    all_met = True
    for metric_name, criterion_str in criteria.items():
      op, target = parse_criterion(criterion_str)
      if op is None:
        # Unparseable criterion — treat as not met
        details.append({
          "criterion": metric_name,
          "target": criterion_str,
          "actual": None,
          "met": False,
        })
        all_met = False
        continue

      actual = _get_metric_value(
        structure_model, metric_name
      )

      # Apply hysteresis for advancement
      if not isinstance(target, bool):
        adj_target = apply_hysteresis(target, op)
        met = evaluate_criterion(actual, op, adj_target)
      else:
        met = evaluate_criterion(actual, op, target)

      if not met:
        all_met = False
      details.append({
        "criterion": metric_name,
        "target": criterion_str,
        "actual": actual,
        "met": met,
      })
    return (all_met, details)

  # ── Skip conditions ─────────────────────────────

  def _check_skip(self, stage, structure_model):
    """Check if stage should be skipped.

    Args:
      stage: StageDef with skip_if.
      structure_model: StructureModel.

    Returns:
      GateResult with action="skip" if skip, None
      otherwise. The caller must call plan.skip_stage()
      then plan.advance() to handle this correctly.
    """
    parsed = parse_gate_condition(stage.skip_if)
    if parsed is None:
      return None
    actual = _get_metric_value(
      structure_model, parsed["metric"]
    )
    if actual is None:
      return None
    if evaluate_gate_condition(parsed, actual, 0):
      return GateResult(
        action="skip",
        reason="skip_if met: %s (actual=%s)"
          % (stage.skip_if, actual),
      )
    return None

  # ── Gate conditions (retreat triggers) ──────────

  def _check_gates(self, plan, stage,
                   structure_model, cycle_number):
    """Check gate conditions for retreat/stop.

    Args:
      plan: StructurePlan.
      stage: StageDef.
      structure_model: StructureModel.
      cycle_number: int.

    Returns:
      GateResult or None.
    """
    for gate in stage.gate_conditions:
      condition_str = gate.get("if", "")
      action_str = gate.get("action", "")
      parsed = parse_gate_condition(condition_str)
      if parsed is None:
        continue
      actual = _get_metric_value(
        structure_model, parsed["metric"]
      )
      if not evaluate_gate_condition(
        parsed, actual, stage.cycles_used
      ):
        continue

      # Gate condition fired
      if action_str == "stop_report":
        return GateResult(
          action="stop",
          reason="gate: %s (actual=%s)"
            % (condition_str, actual),
        )

      if action_str.startswith("retreat_to "):
        target_id = action_str.split(
          "retreat_to ", 1
        )[1].strip()
        return self._evaluate_retreat(
          plan, target_id,
          condition_str, actual,
          cycle_number, structure_model,
          explicit_target=True,
        )

      if action_str == "try_rebuilding":
        # Advance to the model_rebuilding stage if it
        # exists in the plan.  This is stronger than a
        # fallback hint — it skips remaining refine
        # cycles and goes straight to autobuild.
        rebuild_stage = None
        for s in plan.stages:
          if s.id == "model_rebuilding":
            rebuild_stage = s
            break
        if rebuild_stage is not None:
          return GateResult(
            action="advance",
            reason="gate: %s → advancing to "
              "model_rebuilding (autobuild)"
              % condition_str,
          )
        # No model_rebuilding stage — fall back to hint
        return GateResult(
          action="fallback",
          reason="gate suggests rebuilding: %s"
            % condition_str,
        )

      # Generic action — return as fallback
      return GateResult(
        action="fallback",
        reason="gate: %s → %s"
          % (condition_str, action_str),
      )
    return None

  # ── Fallback conditions ─────────────────────────

  def _check_fallbacks(self, stage, structure_model):
    """Check within-stage fallback conditions.

    Returns GateResult with action="fallback" or None.
    """
    for fb in stage.fallbacks:
      condition_str = fb.get("if", "")
      action_str = fb.get("action", "")
      parsed = parse_gate_condition(condition_str)
      if parsed is None:
        continue
      actual = _get_metric_value(
        structure_model, parsed["metric"]
      )
      if evaluate_gate_condition(
        parsed, actual, stage.cycles_used
      ):
        return GateResult(
          action="fallback",
          reason="fallback: %s → %s (actual=%s)"
            % (condition_str, action_str, actual),
        )
    return None

  # ── Stage exhaustion ────────────────────────────

  def _handle_exhaustion(self, plan, stage,
                         structure_model,
                         cycle_number):
    """Handle stage that used all its cycles.

    If some success criteria are met, advance anyway
    (partial success is enough to move forward).
    Otherwise, advance if we're making progress, or
    retreat if stalled.

    Returns GateResult.
    """
    # Check if any success criteria are met
    if stage.success_criteria and structure_model:
      _, details = self._check_success(
        stage.success_criteria, structure_model,
      )
      n_met = sum(1 for d in details if d["met"])
      n_total = len(details)
      if n_met > 0:
        return GateResult(
          action="advance",
          reason=(
            "All steps completed (%d/%d cycles), "
            "%d/%d criteria met — advancing"
            % (stage.cycles_used, stage.max_cycles,
               n_met, n_total)
          ),
          details=details,
        )

    # No criteria met — check if we should retreat
    # Default: advance anyway (the next stage may
    # improve things). Only retreat if a gate
    # condition explicitly says to.
    return GateResult(
      action="advance",
      reason=(
        "All steps completed (%d/%d cycles) — "
        "advancing to next stage"
        % (stage.cycles_used, stage.max_cycles)
      ),
    )

  # ── Retreat evaluation ──────────────────────────

  def _evaluate_retreat(self, plan, target_id,
                        condition_str, actual_value,
                        cycle_number,
                        structure_model,
                        explicit_target=False):
    """Evaluate whether a retreat should proceed.

    Anti-oscillation safeguards:
    1. Strategy Blacklist — is target blacklisted?
    2. Retreat counter — max 2 total
    3. Monotonic progress gate — only retreat if the
       key metric is WORSE than at stage start
    4. Retreat cooldown — at least 2 cycles since last
    5. Retreat depth limit — max 1 stage back (unless
       explicit_target from template gate rule)

    Args:
      plan: StructurePlan.
      target_id: str, stage to retreat to.
      condition_str: str, the triggering condition.
      actual_value: the metric value that triggered.
      cycle_number: int.
      structure_model: StructureModel.
      explicit_target: bool. If True, depth limit
        doesn't apply (template author override).

    Returns:
      GateResult (retreat or continue).
    """
    # --- Safeguard 1: Blacklist ---
    if structure_model and hasattr(
      structure_model, "is_blacklisted"
    ):
      if structure_model.is_blacklisted(target_id):
        return GateResult(
          action="continue",
          reason=(
            "retreat to %s blocked: strategy "
            "blacklisted" % target_id
          ),
        )

    # --- Safeguard 2+4: Counter + cooldown ---
    allowed, block_reason = plan.can_retreat(
      cycle_number
    )
    if not allowed:
      return GateResult(
        action="continue",
        reason=(
          "retreat blocked: %s" % block_reason
        ),
      )

    # --- Safeguard 3: Monotonic progress gate ---
    # Only retreat if the triggering metric is worse
    # than (or equal to) its value at stage start.
    # If the stage is improving, let it continue.
    stage = plan.current_stage()
    if (stage and stage.start_cycle is not None
        and actual_value is not None
        and structure_model):
      progress_block = self._check_monotonic(
        structure_model, stage, actual_value,
        condition_str,
      )
      if progress_block is not None:
        return progress_block

    # --- Safeguard 5: Depth limit ---
    if not explicit_target:
      prev = plan.get_previous_phase()
      if prev and prev.id != target_id:
        # Retreat target is deeper than one stage
        # back. Use one-stage-back instead.
        target_id = prev.id

    # --- Build blacklist entry ---
    bl_entry = None
    if stage:
      metrics = {}
      if structure_model:
        for m in ("r_free", "r_work", "tfz", "llg",
                   "clashscore", "model_map_cc"):
          v = _get_metric_value(structure_model, m)
          if v is not None:
            metrics[m] = v
      bl_entry = {
        "strategy_id": "%s_cycle%d" % (
          stage.id, cycle_number
        ),
        "reason": "gate: %s (actual=%s)"
          % (condition_str, actual_value),
        "metrics_at_retreat": metrics,
      }

    return GateResult(
      action="retreat",
      reason="gate: %s (actual=%s)"
        % (condition_str, actual_value),
      new_phase_id=target_id,
      blacklist_entry=bl_entry,
    )

  def _check_monotonic(self, structure_model,
                       stage, actual_value,
                       condition_str):
    """Check if stage is making progress.

    If the triggering metric has improved since
    stage start, block the retreat.

    Returns:
      GateResult to block retreat, or None to allow.
    """
    # Extract metric name from condition
    parsed = parse_gate_condition(condition_str)
    if parsed is None:
      return None
    metric = parsed.get("metric")
    if metric is None:
      return None

    # Get value at stage start from progress history
    start_val = None
    if hasattr(structure_model, "progress"):
      for entry in structure_model.progress:
        if entry.get("cycle") == stage.start_cycle:
          start_val = entry.get(metric)
          break
    if start_val is None:
      return None  # can't compare, allow retreat

    try:
      start_f = float(start_val)
      actual_f = float(actual_value)
    except (ValueError, TypeError):
      return None

    # Determine if metric improved. For most metrics
    # lower is better (r_free, clashscore), but for
    # some higher is better (tfz, llg, model_map_cc).
    _HIGHER_BETTER = frozenset([
      "tfz", "llg", "model_map_cc",
      "rama_favored", "ligand_cc",
    ])
    if metric in _HIGHER_BETTER:
      improved = actual_f > start_f
    else:
      improved = actual_f < start_f

    if improved:
      return GateResult(
        action="continue",
        reason=(
          "retreat blocked: %s improving "
          "(%.3f → %.3f since stage start)"
          % (metric, start_f, actual_f)
        ),
      )
    return None


# ── Metric extraction helper ───────────────────────

def _get_metric_value(structure_model, metric_name):
  """Get a metric from the structure model.

  Handles StructureModel objects, dicts with
  get_metric method, and plain dicts.

  Returns float, int, bool, or None.
  """
  if structure_model is None:
    return None
  # StructureModel or object with get_metric
  if hasattr(structure_model, "get_metric"):
    return structure_model.get_metric(metric_name)
  # Plain dict (for testing)
  if isinstance(structure_model, dict):
    return structure_model.get(metric_name)
  return None
