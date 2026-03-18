"""
Plan Generator for the Goal-Directed Agent.

Generates a StructurePlan at session start by selecting
the best-matching template from plan_templates.yaml and
customizing it based on available files, user advice,
and data characteristics.

Template selection is deterministic (rule-based).
The LLM is NOT used — the template constrains the
plan structure, ensuring correctness.

Entry points:
  generate_plan(...) -> StructurePlan
  plan_to_directives(plan) -> dict
  check_plan_revision(plan, session_data) -> bool
  repair_plan(plan, user_directives) -> list

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import logging
import os
import re

logger = logging.getLogger(__name__)

try:
  from libtbx.langchain.knowledge.plan_schema import (
    StructurePlan, merge_directives,
  )
  from libtbx.langchain.knowledge.plan_template_loader \
    import (
      load_templates, select_template,
      build_plan_from_template,
    )
except ImportError:
  try:
    from knowledge.plan_schema import (
      StructurePlan, merge_directives,
    )
    from knowledge.plan_template_loader import (
      load_templates, select_template,
      build_plan_from_template,
    )
  except ImportError:
    StructurePlan = None
    merge_directives = None
    load_templates = None
    select_template = None
    build_plan_from_template = None


# ── Main entry point ────────────────────────────────

def generate_plan(data_characteristics=None,
                  user_advice="",
                  available_files=None,
                  directives=None,
                  structure_model=None,
                  cycle_number=0):
  """Generate a StructurePlan for the current session.

  Step 1: Build context from available information.
  Step 2: Select the best-matching template.
  Step 3: Build and customize the plan.

  The plan generator can NEVER produce an invalid plan —
  the template guarantees structural correctness.

  Args:
    data_characteristics: dict from StructureModel or
      session. May be None at session start.
    user_advice: str, user's project advice.
    available_files: list of str, file paths.
    directives: dict, extracted user directives.
    structure_model: StructureModel or None.
    cycle_number: int, current cycle.

  Returns:
    StructurePlan or None if plan generation fails
    (e.g. templates not loadable, no matching template).

  Never raises.
  """
  if StructurePlan is None or load_templates is None:
    return None

  try:
    return _generate_plan_inner(
      data_characteristics, user_advice,
      available_files, directives,
      structure_model, cycle_number,
    )
  except Exception:
    logger.debug(
      "generate_plan failed", exc_info=True,
    )
    return None


def _generate_plan_inner(
  data_characteristics, user_advice,
  available_files, directives,
  structure_model, cycle_number,
):
  """Inner plan generation. May raise."""
  # --- Step 1: Build context ---
  context = _build_context(
    data_characteristics=data_characteristics,
    user_advice=user_advice,
    available_files=available_files,
    directives=directives,
    structure_model=structure_model,
  )

  # --- Step 2: Select template ---
  templates = load_templates()
  if not templates:
    logger.warning(
      "No plan templates loaded — cannot generate"
      " plan"
    )
    return None

  template_id = select_template(templates, context)
  if template_id is None:
    logger.warning(
      "No matching template for context: %s",
      {k: v for k, v in context.items()
       if v is not None},
    )
    return None

  # --- Step 3: Build plan ---
  plan = build_plan_from_template(
    template_id, templates, context,
    cycle_number=cycle_number,
  )
  if plan is None:
    return None

  # --- Step 4: Customize ---
  _customize_plan(plan, context, user_advice)

  # Compute initial strategy hash
  plan.strategy_hash = plan.compute_hash()

  logger.debug(
    "Generated plan: template=%s goal=%s "
    "stages=%d",
    template_id, plan.goal, len(plan.stages),
  )
  return plan


# ── Context building ────────────────────────────────

def _build_context(data_characteristics=None,
                   user_advice="",
                   available_files=None,
                   directives=None,
                   structure_model=None):
  """Build template-selection context from all sources.

  Merges information from files, advice, directives,
  and the structure model (if available from a previous
  run or resumed session).

  Returns:
    dict with keys matching template applicable_when:
      experiment_type, has_search_model, has_sequence,
      has_ligand_code, has_anomalous_atoms, resolution,
      is_twinned.
  """
  ctx = {
    "experiment_type": None,
    "has_search_model": False,
    "has_sequence": False,
    "has_ligand_code": False,
    "has_anomalous_atoms": False,
    "wants_mr_sad": False,
    "wants_polder": False,
    "wants_validation_only": False,
    "model_is_placed": False,
    "resolution": None,
    "is_twinned": None,
  }

  # --- From available files ---
  has_pdb = False
  pdb_files = []
  cif_files = []
  files = available_files or []

  # Ligand PDB hints, split into two tiers:
  #
  # Tier 1 (unambiguous): These words in a PDB filename
  # always mean a ligand coordinate file, never a protein
  # model.  Safe to skip has_search_model.
  _ligand_pdb_certain = ("ligand", "lig_", "lig.")
  #
  # Tier 2 (ambiguous): These MIGHT indicate a ligand
  # file but also appear in protein model names (e.g.
  # drug_resistant_mutant.pdb, kinase_inhibitor_complex.pdb).
  # These set _has_ligand_pdb for the advice-based check
  # below but do NOT block has_search_model on their own.
  _ligand_pdb_broad = (
    "random", "compound", "drug", "inhibitor",
  )

  for f in files:
    ext = os.path.splitext(f)[1].lower()
    bn = os.path.basename(f).lower()
    if ext in (".pdb", ".ent"):
      has_pdb = True
      pdb_files.append(bn)
      # Check if this PDB is unambiguously a ligand file.
      # Only the "certain" hints (ligand, lig_, lig.) are
      # safe to use here — broader hints like "drug" or
      # "inhibitor" can appear in protein model filenames
      # (drug_resistant_mutant.pdb, inhibitor_complex.pdb).
      if any(h in bn for h in _ligand_pdb_certain):
        ctx["has_ligand_code"] = True
        # Don't set has_search_model for ligand PDBs
      else:
        ctx["has_search_model"] = True
    elif ext == ".cif":
      cif_files.append(bn)
    elif ext in (".mtz", ".sca", ".hkl"):
      if ctx["experiment_type"] is None:
        ctx["experiment_type"] = "xray"
    elif ext in (".mrc", ".ccp4", ".map"):
      ctx["experiment_type"] = "cryoem"
    elif ext in (
      ".seq", ".fa", ".fasta", ".pir", ".dat",
    ):
      ctx["has_sequence"] = True

  # Classify CIF files after scanning all files
  # so we know whether a PDB model is present.
  for bn in cif_files:
    if _is_ligand_cif(bn, has_pdb):
      ctx["has_ligand_code"] = True
    else:
      ctx["has_search_model"] = True

  # Detect placed-model signals from filenames.
  # PDB files named *_mr_solution* are already-placed
  # MR models (not search models).  Boundary-aware:
  # require _ before mr_solution (or at start of name)
  # so nmr_solution.pdb doesn't match.
  for bn in pdb_files:
    if ('_mr_solution' in bn or
        bn.startswith('mr_solution')):
      ctx["model_is_placed"] = True
      break

  # Track whether any PDB matched ligand hints (used
  # below by the advice-based ligand detection).
  # Certain hints match regardless; broad hints only
  # match when there are 2+ PDB files (the original
  # guard — prevents "drug_resistant_mutant.pdb" alone
  # from triggering ligand detection, while still
  # catching model.pdb + drug.pdb as model + ligand).
  _has_ligand_pdb = any(
    any(h in pbn for h in _ligand_pdb_certain)
    for pbn in pdb_files
  )
  if not _has_ligand_pdb and len(pdb_files) >= 2:
    _has_ligand_pdb = any(
      any(h in pbn for h in _ligand_pdb_broad)
      for pbn in pdb_files
    )

  # --- From directives ---
  d = directives or {}
  wf = d.get("workflow_preferences", {})
  ps = d.get("program_settings", {})
  # Ligand code in directives
  lf_settings = ps.get("phenix.ligandfit", {})
  if lf_settings.get("ligand_code"):
    ctx["has_ligand_code"] = True
  # Resolution from directives
  if ctx["resolution"] is None:
    default_settings = ps.get("default", {})
    if isinstance(default_settings, dict):
      res_str = default_settings.get("resolution")
      if res_str is not None:
        try:
          ctx["resolution"] = float(res_str)
        except (ValueError, TypeError):
          pass
  # Anomalous atom type → SAD experiment
  autosol_settings = ps.get(
    "phenix.autosol", {})
  if isinstance(autosol_settings, dict):
    if autosol_settings.get("atom_type"):
      ctx["has_anomalous_atoms"] = True

  # --- From user advice (text heuristics) ---
  advice = (user_advice or "").lower()
  # Ligand detection from advice: only set
  # has_ligand_code when there is FILE EVIDENCE
  # of a ligand (CIF already detected, or a PDB
  # with ligand-like name). This prevents
  # impossible ligandfit phases when user says
  # "fit ATP" but provides no ligand file (S11A),
  # while still enabling ligand fitting when files
  # are present but not named with .cif extension
  # (1J4R-ligand tutorial).
  if not ctx["has_ligand_code"]:
    if _has_ligand_pdb or _mentions_ligand(advice):
      # Ligand PDB detected or advice mentions
      # ligand — but only set if there's some
      # file evidence (ligand PDB or multiple PDBs)
      if _has_ligand_pdb:
        ctx["has_ligand_code"] = True
      elif cif_files:
        # CIF files exist but weren't classified
        # as ligand — trust the advice
        ctx["has_ligand_code"] = True
  # Experiment type from advice
  if ctx["experiment_type"] is None:
    if any(w in advice for w in (
      "cryo-em", "cryoem", "cryo em",
      "electron microscopy", "map_to_model",
    )):
      ctx["experiment_type"] = "cryoem"
    elif any(w in advice for w in (
      "x-ray", "xray", "diffraction",
      "molecular replacement", "phaser",
      "refinement", "refine",
    )):
      ctx["experiment_type"] = "xray"

  # Anomalous from advice
  if not ctx["has_anomalous_atoms"]:
    if any(w in advice for w in (
      "sad", "mad", "anomalous", "selenium",
      "se-met", "s-sad", "sulfur",
      "atom_type", "heavy atom",
    )):
      ctx["has_anomalous_atoms"] = True

  # MR-SAD intent: user explicitly wants molecular
  # replacement combined with anomalous phasing.
  # Without this, pure "SAD phasing" advice will
  # select sad_phasing (not mr_sad) even when a
  # search model is present.
  _wants_mr_sad = False
  if any(w in advice for w in (
    "mr-sad", "mr sad", "mrsad", "mr_sad",
  )):
    _wants_mr_sad = True
  elif (
    any(w in advice for w in (
      "molecular replacement", "phaser",
      "mr "))
    and any(w in advice for w in (
      "sad", "anomalous", "mad"))
  ):
    _wants_mr_sad = True
  # If the user explicitly says to run autosol or
  # do SAD phasing AND a search model is present,
  # that's MR-SAD (autosol does MR-SAD when given
  # a partial model).  This is more targeted than
  # checking has_anomalous_atoms (which fires on
  # any data with anomalous signal).
  elif (
    ctx["has_search_model"]
    and any(w in advice for w in (
      "autosol", "sad phasing", "run sad",
      "use sad", "do sad"))
  ):
    _wants_mr_sad = True
  # Also from directives
  if wf.get("use_mr_sad"):
    _wants_mr_sad = True
  ctx["wants_mr_sad"] = _wants_mr_sad

  # Model placement from advice: if user says
  # "refine", "fit ligand", etc. WITHOUT mentioning
  # "molecular replacement", "phaser", "solve", etc.,
  # the model is already placed.
  _mr_keywords = (
    "molecular replacement", "phaser", "solve",
    "mr ", "autosol", "predict",
    "sad ", "sad,", "mad ", "mad,",
    "anomalous", "phasing",
    "heavy atom", "selenium",
    # Docking keywords — if the user mentions
    # docking, the model is NOT yet placed.
    "dock", "docking", "dock_in_map",
    "place the model", "place model",
    "fit into map", "fit model into",
  )
  _placed_keywords = (
    "refine", "refinement", "fit ligand",
    "fit the ligand", "ligandfit",
    "fit atp", "fit nad", "fit fad", "fit heme",
    "polder", "validate", "molprobity",
    "just refine", "only refine",
    "further refine", "continue refin",
  )
  _mentions_mr = any(w in advice for w in _mr_keywords)
  _mentions_placed = any(
    w in advice for w in _placed_keywords
  )
  if _mentions_placed and not _mentions_mr:
    ctx["model_is_placed"] = True

  # Polder override (v115.05): polder ALWAYS implies
  # the model is placed (you need a placed model to
  # compute omit maps).  Override even when "solve"
  # is present, because "Solve omit maps with polder"
  # is a task request, not an MR request.
  if "polder" in advice:
    ctx["model_is_placed"] = True
    ctx["wants_polder"] = True

  # Also check directives for placement signals
  d = directives or {}
  _d_constraints = d.get("constraints", [])
  for c in _d_constraints:
    c_lower = (c.lower()
               if isinstance(c, str) else "")
    if any(kw in c_lower for kw in (
      "refine", "ligandfit", "fit ligand",
      "polder", "validate",
    )):
      ctx["model_is_placed"] = True
      break
  # user_wants_ligandfit → model is placed
  wf = d.get("workflow_preferences", {})
  sc = d.get("stop_conditions", {})
  _after = (sc.get("after_program") or "").lower()
  _prefer = wf.get("prefer_programs", [])
  if ("ligandfit" in _after or
      any("ligandfit" in p.lower()
          for p in _prefer)):
    ctx["model_is_placed"] = True
  # after_program is a refinement program
  if _after in (
    "phenix.refine", "phenix.real_space_refine",
    "phenix.ligandfit", "phenix.polder",
    "phenix.molprobity",
  ):
    ctx["model_is_placed"] = True

  # Polder intent from after_program directive.
  # (The advice-based detection already happened above
  # in the polder override block.)
  if "polder" in _after:
    ctx["wants_polder"] = True

  # v115.09: Validation-only intent from directives
  if wf.get("wants_validation_only"):
    ctx["wants_validation_only"] = True

  # --- From data characteristics ---
  dc = data_characteristics or {}
  if isinstance(dc, dict):
    if dc.get("experiment_type"):
      ctx["experiment_type"] = dc["experiment_type"]
    res = dc.get("resolution")
    if res is not None:
      try:
        ctx["resolution"] = float(res)
      except (ValueError, TypeError):
        pass
    tw = dc.get("twinning", {})
    if isinstance(tw, dict):
      if tw.get("is_twinned") is not None:
        ctx["is_twinned"] = bool(tw["is_twinned"])
    elif dc.get("is_twinned") is not None:
      ctx["is_twinned"] = bool(dc["is_twinned"])
    # Placement detection from session (v114.1)
    if dc.get("model_is_placed"):
      ctx["model_is_placed"] = True

  # --- From structure model (resumed session) ---
  if structure_model is not None:
    sm_dc = None
    if hasattr(structure_model, "data_characteristics"):
      sm_dc = structure_model.data_characteristics
    elif isinstance(structure_model, dict):
      sm_dc = structure_model.get(
        "data_characteristics"
      )
    if sm_dc and isinstance(sm_dc, dict):
      if sm_dc.get("experiment_type"):
        ctx["experiment_type"] = sm_dc[
          "experiment_type"
        ]
      if sm_dc.get("resolution") is not None:
        try:
          ctx["resolution"] = float(
            sm_dc["resolution"]
          )
        except (ValueError, TypeError):
          pass
      tw = sm_dc.get("twinning", {})
      if isinstance(tw, dict) and \
          tw.get("is_twinned") is not None:
        ctx["is_twinned"] = bool(
          tw["is_twinned"]
        )

  # Default experiment type
  if ctx["experiment_type"] is None:
    ctx["experiment_type"] = "xray"

  return ctx


_LIGAND_PATTERN = re.compile(
  r'\b(?:'
  r'ligand|ligandfit|atp|adp|nad|fad|heme|'
  r'cofactor|inhibitor|substrate|drug|'
  r'small.?molecule|binding.?site|'
  r'fit.?ligand|place.?ligand|'
  r'ligand.?fitting'
  r')\b',
  re.IGNORECASE,
)

# Patterns in CIF filenames that indicate ligand
# restraints rather than a coordinate model.
_LIGAND_CIF_MARKERS = re.compile(
  r'elbow|restraint|ligand|_lig[_.]|'
  r'^[a-z0-9]{1,3}\.cif$',
  re.IGNORECASE,
)

# Patterns that indicate a model CIF (not ligand).
_MODEL_CIF_MARKERS = re.compile(
  r'model|refine|autobuild|predicted|'
  r'^\d[a-z0-9]{3}\.cif$',
  re.IGNORECASE,
)


def _is_ligand_cif(basename, has_pdb_model):
  """Classify a CIF file as ligand or model.

  Heuristics (in priority order):
  1. Filename contains ligand markers (elbow,
     restraint, etc.) → ligand.
  2. Filename matches model markers (model, refine,
     PDB code pattern) → model.
  3. Short name (≤3 chars before .cif) → ligand.
  4. If a PDB model is already present, CIF is
     more likely ligand restraints.
  5. Default: model.

  Args:
    basename: str, lowercase filename.
    has_pdb_model: bool, whether a .pdb file is
      present in the file list.

  Returns:
    True if likely a ligand CIF.
  """
  if _LIGAND_CIF_MARKERS.search(basename):
    return True
  if _MODEL_CIF_MARKERS.search(basename):
    return False
  # Short name → ligand (ATP.cif, LIG.cif)
  name_no_ext = os.path.splitext(basename)[0]
  if len(name_no_ext) <= 3:
    return True
  # If PDB already present, CIF is likely restraints
  if has_pdb_model:
    return True
  return False


def _mentions_ligand(text):
  """Check if user advice mentions ligand fitting.

  Returns bool.
  """
  return bool(_LIGAND_PATTERN.search(text or ""))


# ── Plan customization ──────────────────────────────

def _customize_plan(plan, context, user_advice):
  """Customize the plan based on context.

  Adjusts resolution-dependent thresholds for
  intermediate resolutions not covered by the
  lowres/highres templates.

  Args:
    plan: StructurePlan (modified in place).
    context: dict from _build_context.
    user_advice: str.
  """
  resolution = context.get("resolution")
  if resolution is None:
    return

  # Intermediate resolution adjustments (1.5-3.0 Å
  # range, where standard templates apply but
  # thresholds can be tuned).
  if 2.5 <= resolution <= 3.0:
    # Slightly relaxed targets for 2.5-3.0 Å
    for stage in plan.stages:
      sc = stage.success_criteria
      if stage.id == "initial_refinement":
        if sc.get("r_free") == "<0.35":
          sc["r_free"] = "<0.37"
      elif stage.id == "final_refinement":
        if sc.get("r_free") == "<0.25":
          sc["r_free"] = "<0.27"


# ── Plan-to-directives translator (Step 2.4) ───────

def plan_to_directives(plan):
  """Convert current step to reactive agent directives.

  Thin wrapper around StructurePlan.to_directives().
  Exists as a module-level function so ai_agent.py
  can call it without importing StructurePlan.

  Args:
    plan: StructurePlan.

  Returns:
    dict compatible with the directives system.
    Empty dict if plan is None.
  """
  if plan is None:
    return {}
  try:
    return plan.to_directives()
  except Exception:
    return {}


# ── Strategy hash and advice_changed (Step 2.6) ────

def check_plan_revision(plan, session_data):
  """Check if the plan has changed since last stored.

  Compares the plan's current hash against the one
  stored in session_data. If different, updates the
  stored hash and sets advice_changed=True.

  Args:
    plan: StructurePlan.
    session_data: dict (session.data), modified in
      place.

  Returns:
    True if the plan changed, False otherwise.
  """
  if plan is None or not isinstance(
    session_data, dict
  ):
    return False
  try:
    new_hash = plan.compute_hash()
    old_hash = plan.strategy_hash
    if new_hash != old_hash:
      plan.strategy_hash = new_hash
      session_data["advice_changed"] = True
      logger.debug(
        "Plan hash changed: %s -> %s",
        old_hash, new_hash,
      )
      return True
    return False
  except Exception:
    return False


# ── Plan repair (Step 2.5) ──────────────────────────

def repair_plan(plan, user_directives):
  """Repair the plan when user directives conflict.

  Checks if any user skip_programs conflict with
  stages that provide essential data. Applies repair
  rules from the template's if_skipped definitions.

  Args:
    plan: StructurePlan.
    user_directives: dict, user's extracted directives.

  Returns:
    list of repair/warning messages (strings).
    Empty list if no conflicts.

  Never raises.
  """
  if plan is None or not user_directives:
    return []
  try:
    return _repair_plan_inner(
      plan, user_directives
    )
  except Exception:
    logger.debug(
      "repair_plan failed", exc_info=True,
    )
    return []


def _repair_plan_inner(plan, user_directives):
  """Inner plan repair. May raise."""
  messages = []
  ud_wf = user_directives.get(
    "workflow_preferences", {}
  )
  skip_progs = ud_wf.get("skip_programs", [])
  if not skip_progs:
    return messages

  for stage in plan.stages:
    if not stage.provides:
      continue
    # Check if any skip_programs match this phase
    phase_progs = set(stage.programs)
    skipped = phase_progs & set(skip_progs)
    if not skipped:
      continue
    # This phase's programs are being skipped
    prog_str = ", ".join(sorted(skipped))
    if stage.if_skipped:
      # Apply repair rules
      for provided, repair in stage.if_skipped.items():
        messages.append(
          "[Plan Repair] User skips %s. "
          "Substituting: %s"
          % (prog_str, repair)
        )
    else:
      # No repair available — warning only
      provides_str = ", ".join(stage.provides)
      messages.append(
        "[Plan Conflict] User skips %s, but "
        "the plan requires it for %s. "
        "Proceeding without — downstream "
        "decisions may be affected."
        % (prog_str, provides_str)
      )
  return messages
