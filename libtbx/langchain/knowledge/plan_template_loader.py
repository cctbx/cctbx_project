"""
Plan Template Loader for the Goal-Directed Agent.

Loads plan_templates.yaml, resolves template inheritance
(extends/overrides/insert_after), and produces StructurePlan
instances ready for the plan generator.

Entry points:
  load_templates(yaml_path) -> dict of {id: template_dict}
  select_template(templates, context) -> template_id
  build_plan_from_template(template, context) -> StructurePlan

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import copy
import logging
import os

logger = logging.getLogger(__name__)

try:
  from libtbx.langchain.knowledge.plan_schema import (
    StageDef, StructurePlan,
  )
except ImportError:
  try:
    from knowledge.plan_schema import (
      StageDef, StructurePlan,
    )
  except ImportError:
    StageDef = None
    StructurePlan = None


# ── Template loading ────────────────────────────────

_templates_cache = None
_templates_path = None


def load_templates(yaml_path=None):
  """Load and resolve all plan templates from YAML.

  Caches the result — subsequent calls return the
  same dict unless the path changes.

  Args:
    yaml_path: str or None. Path to plan_templates.yaml.
      If None, uses the default location alongside
      this module.

  Returns:
    dict of {template_id: resolved_template_dict}.
    Each template_dict has a flat "stages" list with
    all inheritance resolved. Returns {} on failure.

  Never raises.
  """
  global _templates_cache, _templates_path
  if yaml_path is None:
    yaml_path = os.path.join(
      os.path.dirname(os.path.abspath(__file__)),
      "plan_templates.yaml",
    )
  if (
    _templates_cache is not None
    and _templates_path == yaml_path
  ):
    return _templates_cache

  try:
    result = _load_templates_inner(yaml_path)
    _templates_cache = result
    _templates_path = yaml_path
    return result
  except Exception:
    logger.debug(
      "load_templates failed", exc_info=True,
    )
    return {}


def _load_templates_inner(yaml_path):
  """Inner loader. May raise."""
  import yaml
  with open(yaml_path, "r") as f:
    raw = yaml.safe_load(f)
  if not raw or not isinstance(raw, dict):
    return {}
  templates_raw = raw.get("templates", {})
  if not isinstance(templates_raw, dict):
    return {}

  # First pass: store raw templates
  raw_by_id = {}
  for tid, tdata in templates_raw.items():
    if isinstance(tdata, dict):
      raw_by_id[tid] = tdata

  # Second pass: resolve inheritance
  resolved = {}
  for tid in raw_by_id:
    resolved[tid] = _resolve_template(
      tid, raw_by_id
    )
  return resolved


def _resolve_template(template_id, raw_by_id,
                      depth=0):
  """Resolve extends/overrides/insert_after.

  Returns a template dict with a flat stages list.
  Handles up to 5 levels of inheritance to prevent
  infinite recursion.
  """
  if depth > 5:
    logger.warning(
      "Template inheritance too deep: %s",
      template_id,
    )
    return raw_by_id.get(template_id, {})

  raw = raw_by_id.get(template_id, {})
  if not isinstance(raw, dict):
    return {}

  extends = raw.get("extends")
  if not extends or extends not in raw_by_id:
    # Base template — return as-is
    return copy.deepcopy(raw)

  # Resolve parent first
  parent = _resolve_template(
    extends, raw_by_id, depth + 1
  )
  if not parent:
    return copy.deepcopy(raw)

  # Start with parent's stages
  result = copy.deepcopy(parent)

  # Copy over child metadata
  for key in (
    "description", "applicable_when", "priority",
  ):
    if key in raw:
      result[key] = copy.deepcopy(raw[key])

  # Apply insert_after: insert child stages after
  # the named parent stage
  insert_after = raw.get("insert_after")
  child_phases = raw.get("stages") or raw.get("phases") or []
  if insert_after and child_phases:
    parent_phases = result.get("stages") or result.get("phases") or []
    insert_idx = None
    for i, p in enumerate(parent_phases):
      if p.get("id") == insert_after:
        insert_idx = i + 1
        break
    if insert_idx is not None:
      for j, cp in enumerate(child_phases):
        parent_phases.insert(
          insert_idx + j,
          copy.deepcopy(cp),
        )
    else:
      # Insert point not found — append
      parent_phases.extend(
        copy.deepcopy(cp) for cp in child_phases
      )
    result["stages"] = parent_phases

  # Apply overrides: merge into existing stages
  overrides = raw.get("overrides", {})
  if overrides and isinstance(overrides, dict):
    for stage_id, override_data in overrides.items():
      if not isinstance(override_data, dict):
        continue
      for stage in result.get("stages") or result.get("phases") or []:
        if stage.get("id") == stage_id:
          _merge_phase_override(stage, override_data)
          break

  # Store the template_id and parent reference
  result["_resolved_from"] = template_id
  result["_extends"] = extends

  return result


def _merge_phase_override(stage, override):
  """Merge override data into a stage dict.

  For dicts (success, strategy): merge keys.
  For lists (gate, fallbacks): replace entirely.
  For scalars (max_cycles): replace.
  """
  for key, val in override.items():
    if key in ("success", "strategy"):
      # Merge dict keys
      if key not in stage:
        stage[key] = {}
      if isinstance(val, dict):
        stage[key].update(val)
      else:
        stage[key] = val
    elif key in ("gate", "fallbacks"):
      # Replace list entirely
      stage[key] = (
        list(val) if isinstance(val, list) else val
      )
    else:
      # Replace scalar
      stage[key] = val


# ── Template selection ──────────────────────────────

def select_template(templates, context):
  """Select the best-matching template for a session.

  Matching is deterministic: score each template by
  how well its applicable_when conditions match the
  context, then pick the highest-scoring one. Ties
  broken by priority field.

  Args:
    templates: dict from load_templates().
    context: dict with:
      experiment_type: "xray" or "cryoem"
      has_search_model: bool
      has_sequence: bool
      has_ligand_code: bool
      resolution: float or None
      is_twinned: bool or None

  Returns:
    str (template_id) or None if nothing matches.
  """
  if not templates:
    return None

  best_id = None
  best_score = -1
  best_priority = -1

  for tid, tdata in templates.items():
    if not isinstance(tdata, dict):
      continue
    conditions = tdata.get("applicable_when", {})
    if not isinstance(conditions, dict):
      continue
    score = _score_match(conditions, context)
    if score < 0:
      continue  # hard mismatch
    priority = int(tdata.get("priority", 0))
    # Prefer higher score, then higher priority
    if (score > best_score or
        (score == best_score and
         priority > best_priority)):
      best_id = tid
      best_score = score
      best_priority = priority

  return best_id


def _score_match(conditions, context):
  """Score how well conditions match context.

  Returns:
    int >= 0 for match (higher = more specific).
    -1 for hard mismatch (required condition fails).
  """
  score = 0
  for key, expected in conditions.items():
    actual = context.get(key)

    # Resolution comparison
    if key == "resolution" and isinstance(
      expected, str
    ):
      if actual is None:
        # Unknown resolution = cannot confirm the
        # template's resolution requirement.
        # Hard mismatch — prevents lowres/highres
        # templates from matching when resolution
        # is unknown.
        return -1
      try:
        actual_f = float(actual)
      except (ValueError, TypeError):
        return -1
      if expected.startswith(">"):
        threshold = float(expected[1:])
        if actual_f > threshold:
          score += 1
        else:
          return -1  # hard mismatch
      elif expected.startswith("<"):
        threshold = float(expected[1:])
        if actual_f < threshold:
          score += 1
        else:
          return -1
      continue

    # Boolean conditions
    if isinstance(expected, bool):
      if actual is None:
        if expected is True:
          # Template REQUIRES this flag to be true.
          # Unknown = can't confirm = no match.
          return -1
        else:
          # Template requires flag to be false.
          # Unknown = don't penalize (absence is
          # closer to "not present" than "present").
          continue
      if bool(actual) != expected:
        return -1
      score += 1
      continue

    # String equality
    if actual is None:
      if key == "experiment_type":
        return -1  # required
      continue
    if str(actual).lower() != str(expected).lower():
      return -1
    score += 1

  return score


# ── Plan building ───────────────────────────────────

def build_plan_from_template(template_id, templates,
                             context,
                             cycle_number=0):
  """Build a StructurePlan from a resolved template.

  Converts template stages to StageDef objects and
  creates a StructurePlan ready for the gate evaluator.

  Args:
    template_id: str, key into templates dict.
    templates: dict from load_templates().
    context: dict with session context (for goal text).
    cycle_number: int, current cycle for created_at.

  Returns:
    StructurePlan or None on failure.

  Never raises.
  """
  if StageDef is None or StructurePlan is None:
    return None

  try:
    return _build_plan_inner(
      template_id, templates, context,
      cycle_number,
    )
  except Exception:
    logger.debug(
      "build_plan_from_template failed",
      exc_info=True,
    )
    return None


def _build_plan_inner(template_id, templates,
                      context, cycle_number):
  """Inner plan builder. May raise."""
  tdata = templates.get(template_id)
  if not tdata or not isinstance(tdata, dict):
    return None

  raw_phases = tdata.get("stages") or tdata.get("phases") or []
  if not raw_phases:
    return None

  stages = []
  for rp in raw_phases:
    if not isinstance(rp, dict):
      continue
    stage = StageDef(
      id=rp.get("id", "unknown"),
      programs=rp.get("programs", []),
      max_cycles=rp.get("max_cycles", 5),
      success_criteria=rp.get("success", {}),
      gate_conditions=rp.get("gate", []),
      fallbacks=rp.get("fallbacks", []),
      skip_if=rp.get("skip_if", ""),
      strategy=rp.get("strategy", {}),
      description=rp.get("description", ""),
      provides=rp.get("provides", []),
      if_skipped=rp.get("if_skipped", {}),
    )
    # Prerequisites stored in directives for now
    prereqs = rp.get("prerequisites", {})
    if prereqs:
      stage.directives["prerequisites"] = prereqs
    stages.append(stage)

  # Build goal text
  exp_type = context.get(
    "experiment_type", "unknown"
  )
  resolution = context.get("resolution")
  desc = tdata.get("description", template_id)
  goal_parts = [desc]
  if resolution:
    goal_parts.append("at %.1fÅ" % float(resolution))
  goal = " ".join(goal_parts)

  plan = StructurePlan(
    goal=goal,
    stages=stages,
    created_at_cycle=cycle_number,
    template_id=template_id,
  )
  plan.strategy_hash = plan.compute_hash()
  return plan
