"""
Explanation Engine for the Goal-Directed Agent.

Produces crystallographer-level commentary on what's
happening and why. Three levels:

1. Per-cycle commentary (deterministic template)
2. Stage transition summaries
3. Final/stopped reports (template + data)

No LLM calls for per-cycle commentary. Phase and final
reports use structured templates filled from the
Structure Model (ground truth) to prevent hallucination.

Entry points:
  generate_cycle_commentary(...) -> str
  generate_stage_summary(...) -> str
  generate_final_report(...) -> str
  generate_stopped_report(...) -> str

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import logging

logger = logging.getLogger(__name__)


# ── Per-cycle commentary (Step 5.1) ─────────────────

def generate_cycle_commentary(structure_model,
                              cycle_number,
                              program_name,
                              metrics_before,
                              metrics_after):
  """One-paragraph summary of what happened this cycle.

  Deterministic template — no LLM call needed.
  Reads from the Structure Model and metric deltas.

  Args:
    structure_model: StructureModel or None.
    cycle_number: int.
    program_name: str.
    metrics_before: dict of {metric: value} or None.
    metrics_after: dict of {metric: value} or None.

  Returns:
    str. One paragraph, or empty string if no
    meaningful data.

  Never raises.
  """
  try:
    return _generate_cycle_commentary_inner(
      structure_model, cycle_number,
      program_name, metrics_before, metrics_after,
    )
  except Exception:
    return ""


def _generate_cycle_commentary_inner(
  structure_model, cycle_number, program_name,
  metrics_before, metrics_after,
):
  """Inner commentary generator. May raise."""
  parts = []
  before = metrics_before or {}
  after = metrics_after or {}
  prog = program_name or "unknown"

  # Header
  parts.append(
    "Cycle %d (%s):" % (cycle_number, prog)
  )

  # R-factor changes
  rf_text = _format_metric_change(
    before, after, "r_free", "R-free", lower=True
  )
  if rf_text:
    parts.append(rf_text)

  rw_text = _format_metric_change(
    before, after, "r_work", "R-work", lower=True
  )

  # Model-map CC (cryo-EM)
  cc_text = _format_metric_change(
    before, after, "model_map_cc", "model-map CC",
    lower=False,
  )
  if cc_text:
    parts.append(cc_text)

  # Geometry summary
  geom_parts = []
  cs_text = _format_metric_change(
    before, after, "clashscore", "clashscore",
    lower=True,
  )
  if cs_text:
    geom_parts.append(cs_text)
  rama_text = _format_metric_change(
    before, after, "rama_favored",
    "Ramachandran favored", lower=False,
  )
  if rama_text:
    geom_parts.append(rama_text)
  if geom_parts:
    parts.append(
      "Geometry: %s." % "; ".join(geom_parts)
    )

  # Model contents from Structure Model
  if structure_model is not None:
    contents = _format_model_contents(
      structure_model
    )
    if contents:
      parts.append(contents)

  # Problems
  if structure_model is not None:
    prob_text = _format_problems_brief(
      structure_model
    )
    if prob_text:
      parts.append(prob_text)

  if len(parts) <= 1:
    return ""  # Only header, nothing meaningful

  return " ".join(parts)


def _format_metric_change(before, after, key, label,
                          lower=True):
  """Format a metric change as text.

  Args:
    before, after: dicts with metric values.
    key: str, metric key.
    label: str, human-readable name.
    lower: bool, True if lower is better.

  Returns:
    str or empty.
  """
  val_b = before.get(key)
  val_a = after.get(key)
  if val_a is None:
    return ""
  if val_b is not None:
    try:
      fb = float(val_b)
      fa = float(val_a)
      diff = fa - fb
      if abs(diff) < 1e-6:
        return ""
      if lower:
        direction = "improved" if diff < 0 \
          else "worsened"
      else:
        direction = "improved" if diff > 0 \
          else "worsened"
      return (
        "%s %s from %.3f to %.3f"
        % (label, direction, fb, fa)
      )
    except (ValueError, TypeError):
      pass
  # No before value — just report current
  try:
    return "%s = %.3f" % (label, float(val_a))
  except (ValueError, TypeError):
    return ""


def _format_model_contents(structure_model):
  """Format model contents as a brief clause."""
  ms = structure_model.model_state
  parts = []
  chains = ms.get("chains", [])
  if chains:
    n = len(chains)
    parts.append("%d chain(s)" % n)
  waters = ms.get("waters", 0)
  if waters:
    parts.append("%d waters" % waters)
  ligands = ms.get("ligands", [])
  if ligands:
    lig_strs = []
    for lig in ligands[:3]:
      name = lig.get("name", "?")
      rscc = lig.get("rscc")
      if rscc is not None:
        lig_strs.append(
          "%s (RSCC=%.2f)" % (name, rscc)
        )
      else:
        lig_strs.append(name)
    parts.append(
      "Ligands: %s" % ", ".join(lig_strs)
    )
  ions = ms.get("ions", [])
  if ions:
    parts.append(
      "%d ion(s)" % len(ions)
    )
  if not parts:
    return ""
  return "Model: %s." % ", ".join(parts)


def _format_problems_brief(structure_model):
  """Format current problems as a brief clause."""
  problems = structure_model.get_current_problems()
  if not problems:
    return ""
  if len(problems) == 1:
    return "Issue: %s." % problems[0].get(
      "problem", ""
    )
  return "Issues: %s." % "; ".join(
    p.get("problem", "") for p in problems[:3]
  )


# ── Stage transition summary (Step 5.2) ─────────────

def generate_stage_summary(structure_model, plan,
                           completed_stage,
                           next_stage):
  """Summary of a completed stage.

  Deterministic template summarizing what was
  achieved, what problems remain, and why the agent
  is moving to the next step.

  Args:
    structure_model: StructureModel or None.
    plan: StructurePlan or None.
    completed_stage: StageDef that just finished.
    next_stage: StageDef coming next, or None.

  Returns:
    str. 2-4 sentences.

  Never raises.
  """
  try:
    return _generate_stage_summary_inner(
      structure_model, plan,
      completed_stage, next_stage,
    )
  except Exception:
    return ""


def _generate_stage_summary_inner(
  structure_model, plan,
  completed_stage, next_stage,
):
  """Inner stage summary. May raise."""
  parts = []
  cp = completed_stage
  if cp is None:
    return ""

  # What stage completed
  status = cp.status
  if status == "complete":
    parts.append(
      "Stage '%s' completed in %d cycle(s)."
      % (cp.id, cp.cycles_used)
    )
  elif status == "skipped":
    parts.append(
      "Stage '%s' was skipped." % cp.id
    )
  elif status == "failed":
    parts.append(
      "Stage '%s' failed after %d cycle(s)."
      % (cp.id, cp.cycles_used)
    )
  else:
    parts.append(
      "Stage '%s' ended (%s)."
      % (cp.id, status)
    )

  # Metrics at completion
  rm = cp.result_metrics
  if rm:
    metric_strs = []
    for k, v in rm.items():
      try:
        metric_strs.append(
          "%s=%.3f" % (k, float(v))
        )
      except (ValueError, TypeError):
        metric_strs.append(
          "%s=%s" % (k, v)
        )
    if metric_strs:
      parts.append(
        "Final metrics: %s."
        % ", ".join(metric_strs)
      )

  # Success criteria comparison
  sc = cp.success_criteria
  if sc and rm:
    met = []
    not_met = []
    for k, target in sc.items():
      actual = rm.get(k)
      if actual is not None:
        met.append("%s %s (got %.3f)" % (
          k, target, float(actual)
        ))
      else:
        not_met.append(
          "%s %s (no data)" % (k, target)
        )
    # We don't need to re-evaluate here;
    # just report the comparison.

  # What comes next — omitted from summary because
  # the GUI gate transition block already shows
  # "ADVANCING TO: <next stage>" separately.

  # Current problems
  if structure_model is not None:
    prob_text = _format_problems_brief(
      structure_model
    )
    if prob_text:
      parts.append(prob_text)

  return " ".join(parts)


# ── Final report (Step 5.3) ─────────────────────────

def generate_final_report(structure_model, plan):
  """Comprehensive report at session completion.

  Template-based: skeleton filled from Structure Model
  and plan data. No LLM needed for basic report.

  Args:
    structure_model: StructureModel or None.
    plan: StructurePlan or None.

  Returns:
    str. Multi-section report text.

  Never raises.
  """
  try:
    return _generate_final_report_inner(
      structure_model, plan
    )
  except Exception:
    logger.debug(
      "generate_final_report failed",
      exc_info=True,
    )
    return "Report generation failed."


def _generate_final_report_inner(structure_model,
                                 plan):
  """Inner report generator. May raise."""
  bar = "=" * 50
  lines = [
    bar,
    " STRUCTURE DETERMINATION REPORT",
    bar,
    "",
  ]

  # --- Data characteristics ---
  if structure_model is not None:
    dc = structure_model.data_characteristics
    lines.append("DATA")
    lines.append("----")
    res = dc.get("resolution")
    if res is not None:
      lines.append(
        "  Resolution: %.2f A" % float(res)
      )
    sg = dc.get("space_group")
    if sg:
      lines.append("  Space group: %s" % sg)
    uc = dc.get("unit_cell")
    if uc:
      lines.append("  Unit cell: %s" % str(uc))
    exp = dc.get("experiment_type")
    if exp:
      lines.append("  Experiment: %s" % exp)
    tw = dc.get("twinning", {})
    if tw.get("is_twinned"):
      tf = tw.get("twin_fraction")
      lines.append(
        "  Twinning: yes (fraction=%.2f)"
        % (tf or 0)
      )
    # MR metrics
    tfz = dc.get("mr_tfz")
    if tfz is not None:
      llg = dc.get("mr_llg")
      mr_line = "  MR solution: TFZ=%.1f" % tfz
      if llg is not None:
        mr_line += ", LLG=%.0f" % llg
      lines.append(mr_line)
    # Phasing metrics (SAD/MAD)
    fom = dc.get("phasing_fom")
    if fom is not None:
      ph_line = "  Phasing: FOM=%.3f" % fom
      bcc = dc.get("phasing_bayes_cc")
      if bcc is not None:
        ph_line += ", BAYES-CC=%.1f" % bcc
      sites = dc.get("sites_found")
      if sites is not None:
        ph_line += ", %d sites" % sites
      lines.append(ph_line)
    lines.append("")

  # --- Final model ---
  if structure_model is not None:
    ms = structure_model.model_state
    lines.append("FINAL MODEL")
    lines.append("-----------")
    rf = ms.get("r_free")
    rw = ms.get("r_work")
    if rf is not None:
      lines.append(
        "  R-free: %.4f" % float(rf)
      )
    if rw is not None:
      lines.append(
        "  R-work: %.4f" % float(rw)
      )
    cc = ms.get("model_map_cc")
    if cc is not None:
      lines.append(
        "  Model-map CC: %.3f" % float(cc)
      )
    geom = ms.get("geometry", {})
    cs = geom.get("clashscore")
    if cs is not None:
      lines.append(
        "  Clashscore: %.1f" % float(cs)
      )
    rf_pct = geom.get("rama_favored")
    if rf_pct is not None:
      lines.append(
        "  Ramachandran favored: %.1f%%"
        % (float(rf_pct) * 100)
      )
    # Contents
    chains = ms.get("chains", [])
    if chains:
      lines.append(
        "  Chains: %d" % len(chains)
      )
      for ch in chains[:6]:
        cid = ch.get("chain_id", "?")
        comp = ch.get("completeness")
        if comp is not None:
          lines.append(
            "    %s: %.0f%% complete"
            % (cid, comp * 100)
          )
    waters = ms.get("waters", 0)
    if waters:
      lines.append(
        "  Waters: %d" % waters
      )
    ligands = ms.get("ligands", [])
    if ligands:
      for lig in ligands[:5]:
        name = lig.get("name", "?")
        rscc = lig.get("rscc")
        chain = lig.get("chain", "?")
        resid = lig.get("resid", "?")
        if rscc is not None:
          lines.append(
            "  Ligand: %s (%s/%s, RSCC=%.2f)"
            % (name, chain, resid, rscc)
          )
        else:
          lines.append(
            "  Ligand: %s (%s/%s)"
            % (name, chain, resid)
          )
    ions = ms.get("ions", [])
    if ions:
      lines.append(
        "  Ions: %d" % len(ions)
      )
    lines.append("")

  # --- Phase timeline ---
  if plan is not None and plan.stages:
    lines.append("STAGE TIMELINE")
    lines.append("--------------")
    for i, stage in enumerate(plan.stages):
      if stage.status == "complete":
        marker = "[done]"
      elif stage.status == "skipped":
        marker = "[skip]"
      elif stage.status == "failed":
        marker = "[FAIL]"
      elif stage.status == "active":
        marker = "[....]"
      else:
        marker = "[    ]"
      line = "  %d. %-6s %s" % (
        i + 1, marker, stage.id
      )
      if stage.cycles_used > 0:
        line += " (%d cycle%s)" % (
          stage.cycles_used,
          "s" if stage.cycles_used != 1 else "",
        )
      lines.append(line)
      # Phase metrics
      rm = stage.result_metrics
      if rm:
        m_strs = []
        for k, v in rm.items():
          try:
            m_strs.append(
              "%s=%.3f" % (k, float(v))
            )
          except (ValueError, TypeError):
            pass
        if m_strs:
          lines.append(
            "         %s" % ", ".join(m_strs)
          )
    lines.append("")

  # --- Hypotheses ---
  if structure_model is not None:
    hyps = structure_model.hypotheses
    if hyps:
      lines.append("HYPOTHESES TESTED")
      lines.append("-----------------")
      for h in hyps:
        lines.append(
          "  [%s] %s"
          % (h.status, h.statement[:60])
        )
        if h.resolved_at_cycle is not None:
          lines.append(
            "    Resolved at cycle %d"
            % h.resolved_at_cycle
          )
      lines.append("")

  # --- Outstanding issues ---
  if structure_model is not None:
    problems = (
      structure_model.get_current_problems()
    )
    if problems:
      lines.append("OUTSTANDING ISSUES")
      lines.append("------------------")
      for p in problems[:10]:
        lines.append(
          "  - %s" % p.get("problem", "")
        )
        action = p.get("suggested_action", "")
        if action:
          lines.append(
            "    Suggested: %s" % action
          )
      lines.append("")

  # --- Blacklisted strategies ---
  if structure_model is not None:
    bl = structure_model.strategy_blacklist
    if bl:
      lines.append("STRATEGIES TRIED AND FAILED")
      lines.append("---------------------------")
      for entry in bl:
        lines.append(
          "  - %s" % entry.get("reason", "")
        )
      lines.append("")

  lines.append(bar)
  return "\n".join(lines)


def generate_stopped_report(structure_model, plan):
  """Report when agent stops before completion.

  Emphasizes what went wrong and what to try next.
  More important than the completion report because
  the user needs actionable guidance.

  Args:
    structure_model: StructureModel or None.
    plan: StructurePlan or None.

  Returns:
    str. Multi-section report text.

  Never raises.
  """
  try:
    return _generate_stopped_report_inner(
      structure_model, plan
    )
  except Exception:
    logger.debug(
      "generate_stopped_report failed",
      exc_info=True,
    )
    return "Stopped report generation failed."


def _generate_stopped_report_inner(structure_model,
                                   plan):
  """Inner stopped report. May raise."""
  bar = "=" * 50
  lines = [
    bar,
    " SESSION STOPPED - INCOMPLETE",
    bar,
    "",
  ]

  # --- Best achieved ---
  lines.append("BEST METRICS ACHIEVED")
  lines.append("---------------------")
  if structure_model is not None:
    ms = structure_model.model_state
    rf = ms.get("r_free")
    rw = ms.get("r_work")
    cc_val = ms.get("model_map_cc")
    if rf is not None:
      lines.append(
        "  R-free: %.4f" % float(rf)
      )
    if rw is not None:
      lines.append(
        "  R-work: %.4f" % float(rw)
      )
    if cc_val is not None:
      lines.append(
        "  Model-map CC: %.3f" % float(cc_val)
      )
    if rf is None and rw is None and cc_val is None:
      lines.append(
        "  No R-factors or CC computed."
      )
  else:
    lines.append("  No structure model available.")
  lines.append("")

  # --- What was tried ---
  if plan is not None and plan.stages:
    lines.append("STAGES ATTEMPTED")
    lines.append("----------------")
    for stage in plan.stages:
      if stage.status == "pending":
        continue  # never reached
      lines.append(
        "  %s: %s (%d cycle%s)"
        % (stage.id, stage.status,
           stage.cycles_used,
           "s" if stage.cycles_used != 1
           else "")
      )
      rm = stage.result_metrics
      if rm:
        m_strs = [
          "%s=%.3f" % (k, float(v))
          for k, v in rm.items()
          if _is_numeric(v)
        ]
        if m_strs:
          lines.append(
            "    Metrics: %s"
            % ", ".join(m_strs)
          )
    lines.append("")

  # --- What failed and why ---
  if structure_model is not None:
    bl = structure_model.strategy_blacklist
    if bl:
      lines.append("STRATEGIES THAT FAILED")
      lines.append("----------------------")
      for entry in bl:
        lines.append(
          "  - %s"
          % entry.get("reason", "unknown")
        )
        metrics = entry.get(
          "metrics_at_retreat", {}
        )
        if metrics:
          m_strs = [
            "%s=%.3f" % (k, float(v))
            for k, v in metrics.items()
            if _is_numeric(v)
          ]
          if m_strs:
            lines.append(
              "    At retreat: %s"
              % ", ".join(m_strs)
            )
      lines.append("")

  # --- Recommendations ---
  lines.append("RECOMMENDATIONS")
  lines.append("---------------")
  recommendations = _generate_recommendations(
    structure_model, plan
  )
  if recommendations:
    for rec in recommendations:
      lines.append("  - %s" % rec)
  else:
    lines.append(
      "  No specific recommendations available."
    )
  lines.append("")

  # --- Outstanding issues ---
  if structure_model is not None:
    problems = (
      structure_model.get_current_problems()
    )
    if problems:
      lines.append("OUTSTANDING ISSUES")
      lines.append("------------------")
      for p in problems[:10]:
        lines.append(
          "  - %s" % p.get("problem", "")
        )
      lines.append("")

  lines.append(bar)
  return "\n".join(lines)


def _generate_recommendations(structure_model,
                              plan):
  """Generate actionable recommendations.

  Based on the current state and what was tried,
  suggest concrete next steps for the user.

  Returns list of str.
  """
  recs = []

  if structure_model is None:
    return ["Examine program logs for errors."]

  ms = structure_model.model_state
  dc = structure_model.data_characteristics

  # R-free based recommendations
  rf = ms.get("r_free")
  if rf is not None:
    rf_val = float(rf)
    if rf_val > 0.45:
      recs.append(
        "R-free is very high (%.3f). Consider "
        "trying a different search model or "
        "checking data quality." % rf_val
      )
    elif rf_val > 0.35:
      recs.append(
        "R-free is %.3f. Try manual model "
        "rebuilding in Coot, focusing on "
        "regions with poor density." % rf_val
      )

  # Twinning
  tw = dc.get("twinning", {})
  if tw.get("is_twinned"):
    bl = structure_model.strategy_blacklist
    twin_tried = any(
      "twin" in entry.get("reason", "").lower()
      for entry in bl
    )
    if not twin_tried:
      recs.append(
        "Data appears twinned. Ensure "
        "refinement uses the correct twin law."
      )

  # Ligand issues
  ligands = ms.get("ligands", [])
  for lig in ligands:
    rscc = lig.get("rscc")
    if rscc is not None and float(rscc) < 0.6:
      recs.append(
        "Ligand %s has poor density "
        "(RSCC=%.2f). Verify ligand identity "
        "and placement in Coot."
        % (lig.get("name", "?"), rscc)
      )

  # Clashscore
  geom = ms.get("geometry", {})
  cs = geom.get("clashscore")
  if cs is not None and float(cs) > 15:
    recs.append(
      "Clashscore is high (%.1f). Review "
      "side chain conformations and crystal "
      "contacts." % float(cs)
    )

  if not recs:
    recs.append(
      "Review the model in Coot with the "
      "difference density map to identify "
      "remaining issues."
    )

  return recs


def _is_numeric(v):
  """Check if a value can be formatted as float."""
  try:
    float(v)
    return True
  except (ValueError, TypeError):
    return False
