"""
Derive knowledge base lookup tags from agent state.

Maps the current validation results, workflow stage,
and metrics into category + tags that the KB loader
uses for tag-based matching.

Entry point: derive_tags(...) -> (category, tags)

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function


def derive_tags(
  validation_result=None,
  log_metrics=None,
  workflow_stage=None,
  program_name=None,
  xtriage_results=None,
  resolution=None,
  experiment_type=None,
  r_free_trend=None,
):
  """Derive KB category and tags from current state.

  Maps the agent's current state into the tag
  vocabulary used by the expert knowledge base.
  The returned (category, tags) are passed directly
  to ExpertKnowledgeBase.query().

  Args:
    validation_result: dict from run_validation()
      or None.
    log_metrics: dict with r_work, r_free, etc.
    workflow_stage: str (e.g. "refinement").
    program_name: str (e.g. "phenix.refine").
    xtriage_results: dict from xtriage analysis
      or None.
    resolution: float (data resolution in A).
    experiment_type: "xray" or "cryoem".
    r_free_trend: list of float (R-free per cycle).

  Returns:
    (category, tags) where category is a str and
    tags is a list of str.
  """
  tags = []
  category = _stage_to_category(workflow_stage)

  # Workflow stage as a tag — matches v2 entries
  # that use "refinement", "phasing", etc.
  if workflow_stage and workflow_stage != "unknown":
    tags.append(workflow_stage)

  # Coerce resolution to float; reject invalid
  if resolution is not None:
    try:
      resolution = float(resolution)
    except (ValueError, TypeError):
      resolution = None
  if resolution is not None and resolution > 0:
    tags.append("resolution")
    tags.extend(_resolution_tags(resolution))

  # Experiment type
  if experiment_type == "cryoem":
    tags.append("cryoem")
    category = "cryoem"

  # Program-specific
  if program_name:
    prog = program_name.replace("phenix.", "")
    tags.append(prog)
    # Higher-level tags for v2 KB matching
    tags.extend(
      _program_context_tags(prog)
    )

  # R-free diagnostics (coerce values inside)
  if log_metrics:
    tags.extend(_r_free_tags(log_metrics))

  # R-free trend: stuck?
  if r_free_trend and len(r_free_trend) >= 3:
    trend_tags = _trend_tags(r_free_trend)
    tags.extend(trend_tags)
    if "plateau" in trend_tags:
      tags.append("convergence")

  # Geometry diagnostics
  if validation_result:
    geom_tags = _geometry_tags(validation_result)
    tags.extend(geom_tags)
    if geom_tags:
      tags.append("geometry")

  # Model contents
  if validation_result:
    tags.extend(
      _contents_tags(validation_result)
    )

  # Xtriage
  if xtriage_results:
    tags.extend(
      _xtriage_tags(xtriage_results)
    )

  # Deduplicate while preserving order
  seen = set()
  unique_tags = []
  for t in tags:
    if t not in seen:
      seen.add(t)
      unique_tags.append(t)

  return category, unique_tags


# =========================================================
# Internal tag derivation helpers
# =========================================================

def _stage_to_category(stage):
  """Map workflow stage to KB category."""
  if not stage:
    return "refinement"
  mapping = {
    "data_analysis": "data_assessment",
    "initial": "data_assessment",
    "molecular_replacement": "phasing",
    "phasing": "phasing",
    "refinement": "refinement",
    "building": "refinement",
    "autobuild": "refinement",
    "ligand_fitting": "ligand_fitting",
    "validation": "validation",
    "cryoem_refinement": "cryoem",
    "cryoem_building": "cryoem",
    "cryoem": "cryoem",
  }
  return mapping.get(stage, "refinement")


def _program_context_tags(prog):
  """Higher-level context tags from program name.

  Maps the short program name to broader conceptual
  tags that the v2 KB uses for matching.
  """
  _PROGRAM_TAGS = {
    "refine": [
      "parameters",
    ],
    "real_space_refine": [
      "parameters",
    ],
    "phaser": [
      "molecular_replacement", "phasing",
    ],
    "autosol": [
      "phasing",
    ],
    "autobuild": [
      "model_building",
    ],
    "predict_and_build": [
      "model_building",
    ],
    "map_to_model": [
      "model_building",
    ],
    "ligandfit": [
      "ligand_fitting",
    ],
    "xtriage": [
      "data_assessment",
    ],
    "polder": [
      "ligand_fitting", "validation",
    ],
    "model_vs_data": [
      "validation",
    ],
    "map_sharpening": [
      "map", "cryoem",
    ],
    "resolve_cryo_em": [
      "map", "cryoem",
    ],
  }
  return _PROGRAM_TAGS.get(prog, [])


def _resolution_tags(resolution):
  """Resolution-band tags."""
  if resolution < 1.5:
    return ["atomic_resolution"]
  elif resolution < 2.0:
    return ["high_resolution"]
  elif resolution < 3.0:
    return ["medium_resolution"]
  elif resolution < 4.0:
    return ["low_resolution"]
  else:
    return ["very_low_resolution"]


def _r_free_tags(metrics):
  """R-free diagnostic tags."""
  tags = []
  rf = metrics.get("r_free")
  if rf is None:
    return tags
  try:
    rf = float(rf)
  except (ValueError, TypeError):
    return tags
  if rf < 0 or rf > 1:
    return tags
  tags.append("r_free")
  if rf > 0.50:
    tags.append("r_free_very_high")
  elif rf >= 0.40:
    tags.append("r_free_high")
  elif rf < 0.25:
    tags.append("r_free_good")

  # R-factor gap
  rw = metrics.get("r_work")
  if rw is not None:
    try:
      rw = float(rw)
    except (ValueError, TypeError):
      rw = None
  if rw is not None and 0 <= rw <= 1:
    gap = rf - rw
    if gap > 0.08:
      tags.append("r_factor_gap_large")
      tags.append("overfitting")
    elif gap > 0.05:
      tags.append("r_factor_gap")
  return tags


def _trend_tags(r_free_trend):
  """Tags from R-free trend analysis."""
  tags = []
  last3 = r_free_trend[-3:]
  # Guard against None and string values that could slip in
  # from JSON round-tripping
  coerced = []
  for v in last3:
    if v is None:
      return tags
    try:
      coerced.append(float(v))
    except (ValueError, TypeError):
      return tags
  last3 = coerced
  diffs = [
    round(abs(last3[i] - last3[i + 1]), 4)
    for i in range(len(last3) - 1)
  ]
  if all(d <= 0.005 for d in diffs):
    tags.append("r_free_stuck")
    tags.append("plateau")

  # Improving quickly?
  if len(r_free_trend) >= 2:
    first = r_free_trend[0]
    last = r_free_trend[-1]
    try:
      first = float(first) if first is not None else None
      last = float(last) if last is not None else None
    except (ValueError, TypeError):
      first = last = None
    if first is not None and last is not None:
      total_drop = first - last
      if total_drop > 0.10:
        tags.append("improving")

  return tags


def _geometry_tags(validation_result):
  """Tags from geometry validation."""
  tags = []
  geom = validation_result.get("geometry", {})

  cs = geom.get("clashscore")
  if cs is not None:
    if cs > 40:
      tags.append("geometry_terrible")
      tags.append("clashscore")
    elif cs > 20:
      tags.append("geometry_poor")
      tags.append("clashscore")
    elif cs < 5:
      tags.append("geometry_good")

  ro = geom.get("rama_outliers")
  if ro is not None and ro > 0.01:
    tags.append("ramachandran")
    tags.append("outliers")

  rot = geom.get("rotamer_outliers")
  if rot is not None and rot > 0.05:
    tags.append("rotamer_outliers")

  return tags


def _contents_tags(validation_result):
  """Tags from model contents."""
  tags = []
  mc = validation_result.get("model_contents", {})
  if mc.get("ligands"):
    tags.append("ligand_present")
  if mc.get("waters", 0) > 0:
    tags.append("waters")
  if mc.get("ions"):
    tags.append("ions")
  return tags


def _xtriage_tags(xtriage_results):
  """Tags from xtriage analysis."""
  tags = []
  tf = xtriage_results.get("twin_fraction", 0)
  if tf and tf > 0.05:
    tags.append("twinning")
  if xtriage_results.get("tncs_detected"):
    tags.append("tncs")
  if xtriage_results.get("anisotropy_detected"):
    tags.append("anisotropy")
  if xtriage_results.get("ice_rings"):
    tags.append("ice_rings")
  return tags
