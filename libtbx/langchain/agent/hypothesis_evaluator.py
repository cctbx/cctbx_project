"""
Hypothesis Evaluator for the Goal-Directed Agent.

Evaluates active hypotheses after each cycle using
the Structure Model as ground truth. Manages the
hypothesis lifecycle: proposed → testing → pending →
confirmed/refuted/abandoned.

Also provides the hypothesis section for the THINK
prompt (Step 4.2) and revalidation of confirmed
hypotheses (Step 4.1 spec).

Entry points:
  evaluate_hypotheses(structure_model, cycle) -> list
  revalidate_confirmed(structure_model) -> list
  build_hypothesis_prompt(structure_model) -> str

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import logging
import re

logger = logging.getLogger(__name__)

try:
  from libtbx.langchain.agent.structure_model \
    import Hypothesis
except ImportError:
  try:
    from agent.structure_model import Hypothesis
  except ImportError:
    Hypothesis = None

try:
  from libtbx.langchain.agent.gate_evaluator \
    import (
      parse_criterion, evaluate_criterion,
      _get_metric_value,
    )
except ImportError:
  try:
    from agent.gate_evaluator import (
      parse_criterion, evaluate_criterion,
      _get_metric_value,
    )
  except ImportError:
    parse_criterion = None
    evaluate_criterion = None
    _get_metric_value = None


# ── Hypothesis evaluation result ────────────────────

class HypothesisResult(object):
  """Result of evaluating one hypothesis.

  Attributes:
    hypothesis_id: str.
    action: str. One of:
      "confirmed" — criteria met
      "refuted" — refute criteria met
      "abandoned" — inconclusive after timeout
      "countdown" — test_cycles_remaining decremented
      "waiting" — not yet time to evaluate
      "revalidate" — confirmed hypothesis lost validity
      "still_valid" — confirmed hypothesis still good
    reason: str. Human-readable explanation.
    details: list of per-criterion results.
  """

  def __init__(self, hypothesis_id, action,
               reason="", details=None):
    self.hypothesis_id = str(hypothesis_id)
    self.action = str(action)
    self.reason = str(reason)
    self.details = details if details else []

  def __repr__(self):
    return (
      "HypothesisResult(id=%r, action=%r)"
      % (self.hypothesis_id, self.action)
    )


# ── Main evaluation entry point ─────────────────────

# Maximum cycles after test_cycles_remaining reaches 0
# before abandoning as inconclusive.
_ABANDON_TIMEOUT = 2


def evaluate_hypotheses(structure_model,
                        cycle_number):
  """Evaluate all active hypotheses.

  Called after each cycle by the gate evaluator.
  Processes hypotheses in this order:
  1. "testing" → decrement or move to "pending"
  2. "pending" with remaining cycles → decrement
  3. "pending" with 0 remaining → evaluate

  Args:
    structure_model: StructureModel.
    cycle_number: int.

  Returns:
    list of HypothesisResult. Empty if no active
    hypotheses or structure_model is None.

  Never raises.
  """
  if structure_model is None:
    return []
  try:
    return _evaluate_hypotheses_inner(
      structure_model, cycle_number
    )
  except Exception:
    logger.debug(
      "evaluate_hypotheses failed",
      exc_info=True,
    )
    return []


def _evaluate_hypotheses_inner(structure_model,
                               cycle_number):
  """Inner evaluation. May raise."""
  results = []
  hypotheses = getattr(
    structure_model, "hypotheses", []
  )

  for h in hypotheses:
    if h.status == "testing":
      # Test program just ran — transition to pending
      h.status = "pending"
      if h.test_cycles_remaining > 0:
        h.test_cycles_remaining -= 1
        results.append(HypothesisResult(
          h.id, "countdown",
          reason=(
            "test running, %d cycle(s) remaining"
            % h.test_cycles_remaining
          ),
        ))
      else:
        # Ready for evaluation next cycle
        results.append(HypothesisResult(
          h.id, "waiting",
          reason="test complete, evaluating "
                 "next cycle",
        ))

    elif h.status == "pending":
      if h.test_cycles_remaining > 0:
        h.test_cycles_remaining -= 1
        results.append(HypothesisResult(
          h.id, "countdown",
          reason=(
            "stabilizing, %d cycle(s) remaining"
            % h.test_cycles_remaining
          ),
        ))
      else:
        # Evaluate confirm/refute criteria
        result = _evaluate_single_hypothesis(
          h, structure_model, cycle_number
        )
        results.append(result)

  return results


def _evaluate_single_hypothesis(hypothesis,
                                structure_model,
                                cycle_number):
  """Evaluate confirm/refute criteria for one hypothesis.

  Order:
  1. Check confirm_if → confirmed
  2. Check refute_if → refuted
  3. If inconclusive and past timeout → abandoned
  4. Otherwise → waiting (try again next cycle)

  Args:
    hypothesis: Hypothesis object.
    structure_model: StructureModel.
    cycle_number: int.

  Returns:
    HypothesisResult.
  """
  h = hypothesis

  # --- Check confirmation ---
  if h.confirm_if:
    confirmed, conf_details = _check_criteria(
      h.confirm_if, structure_model
    )
    if confirmed:
      h.status = "confirmed"
      h.resolved_at_cycle = cycle_number
      return HypothesisResult(
        h.id, "confirmed",
        reason="confirm criteria met: %s"
          % h.confirm_if,
        details=conf_details,
      )

  # --- Check refutation ---
  if h.refute_if:
    refuted, ref_details = _check_criteria(
      h.refute_if, structure_model
    )
    if refuted:
      h.status = "refuted"
      h.resolved_at_cycle = cycle_number
      return HypothesisResult(
        h.id, "refuted",
        reason="refute criteria met: %s"
          % h.refute_if,
        details=ref_details,
      )

  # --- Check abandon timeout ---
  cycles_since = cycle_number - h.proposed_at_cycle
  if cycles_since > h.test_cycles_remaining + \
      _ABANDON_TIMEOUT + 2:
    h.status = "abandoned"
    h.resolved_at_cycle = cycle_number
    return HypothesisResult(
      h.id, "abandoned",
      reason="inconclusive after %d cycles"
        % cycles_since,
    )

  return HypothesisResult(
    h.id, "waiting",
    reason="criteria not yet met, continuing",
  )


# ── Criteria checking ───────────────────────────────

# Criteria string format:
#   "metric op value [AND metric op value ...]"
# Examples:
#   "r_free < 0.30"
#   "ligand_cc > 0.7 AND clashscore < 5"
#   "anomalous_peak > 5"

_AND_SPLIT = re.compile(
  r'\s+AND\s+', re.IGNORECASE
)

_CRITERION_PART = re.compile(
  r'^(\w+)\s*([<>=!]+)\s*(-?[0-9.]+)$'
)


def _check_criteria(criteria_str, structure_model):
  """Check a compound criteria string.

  All clauses must be met (AND logic).

  Args:
    criteria_str: str, e.g.
      "ligand_cc > 0.7 AND clashscore < 5"
    structure_model: StructureModel.

  Returns:
    (all_met: bool, details: list of dicts).
  """
  if not criteria_str or parse_criterion is None:
    return (False, [])

  clauses = _AND_SPLIT.split(criteria_str.strip())
  if not clauses:
    return (False, [])

  details = []
  all_met = True

  for clause in clauses:
    clause = clause.strip()
    m = _CRITERION_PART.match(clause)
    if not m:
      details.append({
        "clause": clause,
        "met": False,
        "reason": "unparseable",
      })
      all_met = False
      continue

    metric = m.group(1)
    op = m.group(2)
    try:
      target = float(m.group(3))
    except ValueError:
      all_met = False
      continue

    actual = _get_metric_value(
      structure_model, metric
    )
    met = evaluate_criterion(actual, op, target)
    if not met:
      all_met = False
    details.append({
      "clause": clause,
      "metric": metric,
      "operator": op,
      "target": target,
      "actual": actual,
      "met": met,
    })

  return (all_met, details)


# ── Revalidation of confirmed hypotheses ────────────

_REVALIDATION_THRESHOLDS = {
  # Placed feature B-factor too high
  "high_b_factor": 80.0,
  # Ligand RSCC dropped
  "low_ligand_rscc": 0.5,
  # R-free spike after feature added
  "rfree_spike": 0.02,
}


def revalidate_confirmed(structure_model):
  """Revalidate all confirmed hypotheses.

  Called once per cycle. Checks if any confirmed
  hypothesis has lost validity based on model
  quality indicators.

  Args:
    structure_model: StructureModel.

  Returns:
    list of HypothesisResult with action
    "revalidate" or "still_valid".

  Never raises.
  """
  if structure_model is None:
    return []
  try:
    return _revalidate_inner(structure_model)
  except Exception:
    logger.debug(
      "revalidate_confirmed failed",
      exc_info=True,
    )
    return []


def _revalidate_inner(structure_model):
  """Inner revalidation. May raise."""
  results = []
  hypotheses = getattr(
    structure_model, "hypotheses", []
  )

  for h in hypotheses:
    if h.status != "confirmed":
      continue

    # --- Check for degradation signals ---
    reasons = []

    # Check ligand RSCC if hypothesis involved
    # ligand placement
    if "ligand" in h.statement.lower() or \
        "rscc" in h.confirm_if.lower():
      rscc = _get_metric_value(
        structure_model, "ligand_cc"
      )
      if rscc is not None and rscc < \
          _REVALIDATION_THRESHOLDS["low_ligand_rscc"]:
        reasons.append(
          "ligand RSCC dropped to %.2f" % rscc
        )

    # Check R-free for general degradation
    # (simplified — full implementation would
    # track R-free at confirmation time and compare)
    r_free = _get_metric_value(
      structure_model, "r_free"
    )
    if r_free is not None and r_free > 0.45:
      reasons.append(
        "R-free degraded to %.3f" % r_free
      )

    if reasons:
      h.status = "proposed"
      h.revalidation_reason = "; ".join(reasons)
      h.resolved_at_cycle = None
      results.append(HypothesisResult(
        h.id, "revalidate",
        reason="re-evaluation needed: %s"
          % "; ".join(reasons),
      ))
    else:
      results.append(HypothesisResult(
        h.id, "still_valid",
        reason="confirmed hypothesis unchanged",
      ))

  return results


# ── THINK prompt builder (Step 4.2) ─────────────────

def build_hypothesis_prompt(structure_model):
  """Build the hypothesis section for the THINK prompt.

  Returns a text block that is inserted into the
  advanced-mode THINK prompt after the Structure Model
  summary. The block either:
  - Shows the active hypothesis and instructs the LLM
    not to propose a new one
  - Invites the LLM to propose a hypothesis based on
    current problems/unmodeled density

  Args:
    structure_model: StructureModel or None.

  Returns:
    str. Empty string if no hypothesis context or
    structure_model is None.

  Never raises.
  """
  if structure_model is None:
    return ""
  try:
    return _build_prompt_inner(structure_model)
  except Exception:
    return ""


def _build_prompt_inner(structure_model):
  """Inner prompt builder. May raise."""
  active = structure_model.get_active_hypothesis()

  if active is not None:
    # Active hypothesis exists — show it and block
    # new proposals
    lines = [
      "",
      "=== ACTIVE HYPOTHESIS ===",
      "Hypothesis: %s" % active.statement,
      "Status: %s" % active.status,
      "Test: %s %s" % (
        active.test_program,
        _format_params(active.test_parameters),
      ),
    ]
    if active.test_cycles_remaining > 0:
      lines.append(
        "Waiting: %d cycle(s) before evaluation"
        % active.test_cycles_remaining
      )
    lines.append(
      "Confirm if: %s" % active.confirm_if
    )
    lines.append(
      "Refute if: %s" % active.refute_if
    )
    if active.revalidation_reason:
      lines.append(
        "RE-EVALUATING because: %s"
        % active.revalidation_reason
      )
    lines.append("")
    lines.append(
      "DO NOT propose a new hypothesis. The active"
      " hypothesis must be resolved first."
    )
    lines.append(
      "========================="
    )
    return "\n".join(lines)

  # No active hypothesis — check if conditions
  # warrant proposing one
  problems = structure_model.get_current_problems()
  has_peaks = False
  dp = structure_model.model_state.get(
    "diff_peaks", {}
  )
  pos_peaks = dp.get("positive", [])
  if len(pos_peaks) >= 2:
    has_peaks = True

  if not problems and not has_peaks:
    return ""  # Nothing to hypothesize about

  lines = [
    "",
    "=== HYPOTHESIS OPPORTUNITY ===",
    "No active hypothesis. Based on the current "
    "model state:",
  ]

  # Summarize what's available for hypothesis
  if has_peaks:
    lines.append(
      "- %d positive difference density peak(s) "
      "> 4 sigma" % len(pos_peaks)
    )
    for peak in pos_peaks[:3]:
      lines.append(
        "  %.1f sigma near %s"
        % (peak.get("height", 0),
           peak.get("near_residue", "?"))
      )

  for prob in problems[:3]:
    lines.append(
      "- %s" % prob.get("problem", "")
    )

  lines.append("")
  lines.append(
    "If you identify a testable hypothesis, "
    "respond with:"
  )
  lines.append(
    '{"hypothesis": "...", "test_program": "...",'
    ' "test_parameters": {...},'
  )
  lines.append(
    ' "confirm_if": "metric op value",'
    ' "refute_if": "metric op value"}'
  )
  lines.append(
    "Otherwise, focus on the current plan stage."
  )
  lines.append(
    "==============================="
  )
  return "\n".join(lines)


def _format_params(params):
  """Format test parameters as a compact string."""
  if not params:
    return ""
  parts = [
    "%s=%s" % (k, v)
    for k, v in params.items()
  ]
  return "(%s)" % ", ".join(parts)
