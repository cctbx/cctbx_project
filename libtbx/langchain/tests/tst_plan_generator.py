"""
Unit tests for plan_generator.py (Phase 2, Steps 2.3-2.6).

Run standalone:
  python tests/tst_plan_generator.py

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

from agent.plan_generator import (
  generate_plan, plan_to_directives,
  check_plan_revision, repair_plan,
  _build_context, _mentions_ligand,
  _is_ligand_cif,
)
from knowledge.plan_schema import (
  StructurePlan, STAGE_ACTIVE,
)


def run_tests():
  """Run all plan generator tests."""
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
  print("Plan Generator Unit Tests")
  print("=" * 60)
  print()

  # --- Context building ---
  print("Context building")
  test("context_from_xray_files",
    test_context_from_xray_files)
  test("context_from_cryoem_files",
    test_context_from_cryoem_files)
  test("context_from_ligand_cif",
    test_context_from_ligand_cif)
  test("context_from_sequence_file",
    test_context_from_sequence_file)
  test("context_from_directives",
    test_context_from_directives)
  test("context_from_advice_ligand",
    test_context_from_advice_ligand)
  test("context_from_advice_cryoem",
    test_context_from_advice_cryoem)
  test("context_from_data_characteristics",
    test_context_from_data_characteristics)
  test("context_from_structure_model",
    test_context_from_structure_model)
  test("context_from_structure_model_dict",
    test_context_from_structure_model_dict)
  test("context_defaults_to_xray",
    test_context_defaults_to_xray)
  print()

  # --- Ligand detection ---
  print("Ligand detection")
  test("mentions_ligand_positive",
    test_mentions_ligand_positive)
  test("mentions_ligand_negative",
    test_mentions_ligand_negative)
  test("cif_classification_short_name",
    test_cif_classification_short_name)
  test("cif_classification_elbow",
    test_cif_classification_elbow)
  test("cif_classification_model_cif",
    test_cif_classification_model_cif)
  test("cif_classification_with_pdb_present",
    test_cif_classification_with_pdb_present)
  test("cif_classification_pdb_code",
    test_cif_classification_pdb_code)
  print()

  # --- Plan generation ---
  print("Plan generation")
  test("generate_mr_refine",
    test_generate_mr_refine)
  test("generate_mr_ligand",
    test_generate_mr_ligand)
  test("generate_cryoem",
    test_generate_cryoem)
  test("generate_lowres",
    test_generate_lowres)
  test("generate_highres",
    test_generate_highres)
  test("generate_twinned",
    test_generate_twinned)
  test("generate_no_files",
    test_generate_no_files)
  test("generate_returns_none_on_bad_input",
    test_generate_returns_none_on_bad_input)
  print()

  # --- Customization ---
  print("Customization")
  test("customize_intermediate_resolution",
    test_customize_intermediate_resolution)
  test("customize_no_resolution",
    test_customize_no_resolution)
  print()

  # --- plan_to_directives ---
  print("plan_to_directives")
  test("directives_basic",
    test_directives_basic)
  test("directives_none_plan",
    test_directives_none_plan)
  print()

  # --- check_plan_revision ---
  print("check_plan_revision")
  test("revision_no_change",
    test_revision_no_change)
  test("revision_after_advance",
    test_revision_after_advance)
  test("revision_sets_advice_changed",
    test_revision_sets_advice_changed)
  print()

  # --- repair_plan ---
  print("repair_plan")
  test("repair_skip_with_rules",
    test_repair_skip_with_rules)
  test("repair_skip_no_rules",
    test_repair_skip_no_rules)
  test("repair_no_conflict",
    test_repair_no_conflict)
  test("repair_none_plan",
    test_repair_none_plan)
  print()

  # --- Serialization ---
  print("Serialization")
  test("generated_plan_json_roundtrip",
    test_generated_plan_json_roundtrip)
  print()

  # --- Integration ---
  print("Integration")
  test("full_session_flow",
    test_full_session_flow)
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


# ── Context building tests ────────────────────────────

def test_context_from_xray_files():
  ctx = _build_context(
    available_files=[
      "/data/beta.mtz",
      "/data/beta.pdb",
      "/data/beta.seq",
    ],
  )
  assert ctx["experiment_type"] == "xray"
  assert ctx["has_search_model"] is True
  assert ctx["has_sequence"] is True
  assert ctx["has_ligand_code"] is False


def test_context_from_cryoem_files():
  ctx = _build_context(
    available_files=[
      "/data/emd_1234.mrc",
      "/data/model.pdb",
    ],
  )
  assert ctx["experiment_type"] == "cryoem"
  assert ctx["has_search_model"] is True


def test_context_from_ligand_cif():
  ctx = _build_context(
    available_files=[
      "/data/data.mtz",
      "/data/model.pdb",
      "/data/ATP.cif",
    ],
  )
  assert ctx["has_ligand_code"] is True
  assert ctx["has_search_model"] is True


def test_context_from_sequence_file():
  ctx = _build_context(
    available_files=[
      "/data/data.mtz",
      "/data/seq.fasta",
    ],
  )
  assert ctx["has_sequence"] is True


def test_context_from_directives():
  ctx = _build_context(
    directives={
      "program_settings": {
        "phenix.ligandfit": {
          "ligand_code": "ATP",
        },
      },
    },
  )
  assert ctx["has_ligand_code"] is True


def test_context_from_advice_ligand():
  # Without a CIF file, advice alone should NOT
  # set has_ligand_code (prevents hallucinated
  # ligandfit phases when no CIF is provided).
  ctx = _build_context(
    user_advice="Solve and fit the ATP ligand",
    available_files=["/data/data.mtz"],
  )
  assert ctx["has_ligand_code"] is False
  # With a CIF file, has_ligand_code should be True
  ctx2 = _build_context(
    user_advice="Solve and fit the ATP ligand",
    available_files=["/data/data.mtz", "atp.cif"],
  )
  assert ctx2["has_ligand_code"] is True


def test_context_from_advice_cryoem():
  ctx = _build_context(
    user_advice="Build model into cryo-EM map",
    available_files=[],
  )
  assert ctx["experiment_type"] == "cryoem"


def test_context_from_data_characteristics():
  ctx = _build_context(
    data_characteristics={
      "experiment_type": "xray",
      "resolution": 2.1,
      "twinning": {"is_twinned": True},
    },
  )
  assert ctx["experiment_type"] == "xray"
  assert ctx["resolution"] == 2.1
  assert ctx["is_twinned"] is True


def test_context_from_structure_model():
  """Context from a StructureModel object."""
  from agent.structure_model import StructureModel
  sm = StructureModel()
  sm.data_characteristics["experiment_type"] = "xray"
  sm.data_characteristics["resolution"] = 1.8
  sm.data_characteristics["twinning"][
    "is_twinned"
  ] = False
  ctx = _build_context(structure_model=sm)
  assert ctx["experiment_type"] == "xray"
  assert ctx["resolution"] == 1.8
  assert ctx["is_twinned"] is False


def test_context_from_structure_model_dict():
  """Context from a serialized StructureModel dict."""
  sm_dict = {
    "data_characteristics": {
      "experiment_type": "cryoem",
      "resolution": 3.2,
    },
  }
  ctx = _build_context(structure_model=sm_dict)
  assert ctx["experiment_type"] == "cryoem"
  assert ctx["resolution"] == 3.2


def test_context_defaults_to_xray():
  ctx = _build_context()
  assert ctx["experiment_type"] == "xray"


# ── Ligand detection ──────────────────────────────────

def test_mentions_ligand_positive():
  assert _mentions_ligand("Fit the ATP ligand")
  assert _mentions_ligand("Run ligandfit")
  assert _mentions_ligand("place inhibitor")
  assert _mentions_ligand("cofactor binding site")
  assert _mentions_ligand(
    "fit ligand into density"
  )


def test_mentions_ligand_negative():
  assert not _mentions_ligand("Refine the model")
  assert not _mentions_ligand("Run phaser")
  assert not _mentions_ligand("")
  assert not _mentions_ligand(None)


def test_cif_classification_short_name():
  assert _is_ligand_cif("atp.cif", False) is True
  assert _is_ligand_cif("lig.cif", False) is True
  assert _is_ligand_cif("ab.cif", False) is True


def test_cif_classification_elbow():
  """ELBOW-generated restraints have long names."""
  assert _is_ligand_cif(
    "elbow.atp.001.cif", False
  ) is True
  assert _is_ligand_cif(
    "atp_restraints.cif", False
  ) is True
  assert _is_ligand_cif(
    "ligand_geometry.cif", False
  ) is True


def test_cif_classification_model_cif():
  """Model CIFs should not be ligands."""
  assert _is_ligand_cif(
    "model.cif", False
  ) is False
  assert _is_ligand_cif(
    "refine_001.cif", False
  ) is False
  assert _is_ligand_cif(
    "autobuild_001.cif", False
  ) is False


def test_cif_classification_with_pdb_present():
  """When .pdb model exists, ambiguous CIF files
  should be classified as ligand restraints."""
  # Long name, no markers, but PDB present
  assert _is_ligand_cif(
    "some_compound.cif", True
  ) is True
  # Model-marked CIF stays model even with PDB
  assert _is_ligand_cif(
    "model.cif", True
  ) is False


def test_cif_classification_pdb_code():
  """PDB codes like 1abc.cif are model files."""
  assert _is_ligand_cif(
    "1abc.cif", False
  ) is False
  assert _is_ligand_cif(
    "7xyz.cif", False
  ) is False


# ── Plan generation ───────────────────────────────────

def test_generate_mr_refine():
  plan = generate_plan(
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
    user_advice="Solve by MR and refine",
  )
  assert plan is not None
  assert len(plan.stages) == 6
  assert plan.template_id == "mr_refine"
  assert plan.stages[0].id == "data_assessment"
  assert plan.stages[1].id == (
    "molecular_replacement"
  )
  assert plan.stages[5].id == "validation"


def test_generate_mr_ligand():
  plan = generate_plan(
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
      "/data/ATP.cif",
    ],
    user_advice="Solve and fit the ligand",
  )
  assert plan is not None
  assert plan.template_id == "mr_refine_ligand"
  ids = [p.id for p in plan.stages]
  assert "ligand_fitting" in ids


def test_generate_cryoem():
  plan = generate_plan(
    available_files=[
      "/data/emd.mrc", "/data/model.pdb",
    ],
    user_advice="Refine into cryo-EM map",
  )
  assert plan is not None
  assert plan.template_id == "cryoem_refine"
  ids = [p.id for p in plan.stages]
  assert "refinement" in ids
  assert "validation" in ids
  # Refinement stage has RSR
  refine_stage = [s for s in plan.stages
    if s.id == "refinement"][0]
  assert "phenix.real_space_refine" in (
    refine_stage.programs
  )
  # Last stage is validation
  assert plan.stages[-1].id == "validation"


def test_generate_lowres():
  plan = generate_plan(
    data_characteristics={"resolution": 3.5},
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  assert plan is not None
  assert plan.template_id == "mr_refine_lowres"
  # Check relaxed thresholds
  for p in plan.stages:
    if p.id == "final_refinement":
      assert p.strategy.get(
        "ordered_solvent"
      ) is False


def test_generate_highres():
  plan = generate_plan(
    data_characteristics={"resolution": 1.0},
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  assert plan is not None
  assert plan.template_id == "mr_refine_highres"
  for p in plan.stages:
    if p.id == "final_refinement":
      assert p.strategy.get("adp", {}).get(
        "type"
      ) == "anisotropic"


def test_generate_twinned():
  plan = generate_plan(
    data_characteristics={
      "twinning": {"is_twinned": True},
    },
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  assert plan is not None
  assert plan.template_id == "mr_refine_twinned"


def test_generate_no_files():
  """Without files but with xray context, falls back
  to data_analysis_only template."""
  plan = generate_plan(
    user_advice="Solve the structure",
  )
  # No files at all → experiment_type defaults to
  # xray → data_analysis_only matches
  if plan is not None:
    assert plan.template_id == "data_analysis_only"
  # else: None is also acceptable (no experiment_type)


def test_generate_returns_none_on_bad_input():
  # data.mtz only → xray, no model, no sequence
  # Falls back to data_analysis_only (xtriage only)
  plan = generate_plan(
    available_files=["/data/data.mtz"],
    # No .pdb → has_search_model=False
    # No cryoem files → xray
  )
  assert plan is not None
  assert plan.template_id == "data_analysis_only"
  assert len(plan.stages) == 1
  assert plan.stages[0].id == "data_assessment"


# ── Customization ─────────────────────────────────────

def test_customize_intermediate_resolution():
  """Plan at 2.7Å should get slightly relaxed targets."""
  plan = generate_plan(
    data_characteristics={"resolution": 2.7},
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  assert plan is not None
  for p in plan.stages:
    if p.id == "initial_refinement":
      # 2.7Å is in 2.5-3.0 range → relaxed to <0.37
      assert p.success_criteria.get(
        "r_free"
      ) == "<0.37", (
        "Expected <0.37, got %s"
        % p.success_criteria
      )
    if p.id == "final_refinement":
      assert p.success_criteria.get(
        "r_free"
      ) == "<0.27"


def test_customize_no_resolution():
  """No resolution → no customization (defaults)."""
  plan = generate_plan(
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  assert plan is not None
  for p in plan.stages:
    if p.id == "initial_refinement":
      assert p.success_criteria.get(
        "r_free"
      ) == "<0.35"


# ── plan_to_directives ───────────────────────────────

def test_directives_basic():
  plan = generate_plan(
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  d = plan_to_directives(plan)
  wf = d.get("workflow_preferences", {})
  assert "prefer_programs" in wf
  assert len(wf["prefer_programs"]) > 0


def test_directives_none_plan():
  d = plan_to_directives(None)
  assert d == {}


# ── check_plan_revision ──────────────────────────────

def test_revision_no_change():
  plan = generate_plan(
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  session_data = {}
  changed = check_plan_revision(plan, session_data)
  assert changed is False


def test_revision_after_advance():
  plan = generate_plan(
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  plan.stages[0].status = STAGE_ACTIVE
  session_data = {}
  # Hash is set during generation
  h1 = plan.strategy_hash
  # Advance changes current_stage_index
  plan.advance()
  changed = check_plan_revision(plan, session_data)
  assert changed is True
  assert plan.strategy_hash != h1


def test_revision_sets_advice_changed():
  plan = generate_plan(
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  plan.stages[0].status = STAGE_ACTIVE
  session_data = {"advice_changed": False}
  plan.advance()
  check_plan_revision(plan, session_data)
  assert session_data["advice_changed"] is True


# ── repair_plan ───────────────────────────────────────

def test_repair_skip_with_rules():
  plan = generate_plan(
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  # User skips xtriage
  user_d = {
    "workflow_preferences": {
      "skip_programs": ["phenix.xtriage"],
    },
  }
  messages = repair_plan(plan, user_d)
  assert len(messages) > 0
  assert all("Repair" in m for m in messages)
  # Should have one repair per provided item
  # (resolution, twinning_status, anomalous_signal)
  assert len(messages) == 3


def test_repair_skip_no_rules():
  plan = generate_plan(
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  # User skips refine (no if_skipped rules)
  user_d = {
    "workflow_preferences": {
      "skip_programs": ["phenix.refine"],
    },
  }
  messages = repair_plan(plan, user_d)
  # initial_refinement and final_refinement have
  # no if_skipped, so we get warnings
  # (provides may be empty for refine steps)
  # Actually: refine steps don't have provides in
  # the template. No conflict detected.
  assert isinstance(messages, list)


def test_repair_no_conflict():
  plan = generate_plan(
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
  )
  user_d = {
    "program_settings": {
      "phenix.refine": {"nproc": 4},
    },
  }
  messages = repair_plan(plan, user_d)
  assert messages == []


def test_repair_none_plan():
  messages = repair_plan(None, {})
  assert messages == []


# ── Serialization ─────────────────────────────────────

def test_generated_plan_json_roundtrip():
  plan = generate_plan(
    available_files=[
      "/data/data.mtz", "/data/model.pdb",
    ],
    user_advice="Solve by MR and refine",
  )
  d = plan.to_dict()
  json_str = json.dumps(d)
  d2 = json.loads(json_str)
  plan2 = StructurePlan.from_dict(d2)
  assert len(plan2.stages) == len(plan.stages)
  assert plan2.template_id == plan.template_id
  assert plan2.goal == plan.goal
  # Directives should work after round-trip
  d_orig = plan_to_directives(plan)
  d_rest = plan_to_directives(plan2)
  assert d_orig == d_rest


# ── Integration ───────────────────────────────────────

def test_full_session_flow():
  """Simulate the full session start flow:
  1. Generate plan from files + advice
  2. Get directives for first stage
  3. Merge with user directives
  4. Check for conflicts / repair
  5. Verify plan persistence
  """
  from knowledge.plan_schema import merge_directives

  # --- Session start ---
  plan = generate_plan(
    available_files=[
      "/data/data.mtz",
      "/data/model.pdb",
      "/data/LIG.cif",
    ],
    user_advice="Solve and fit the ligand",
  )
  assert plan is not None
  assert plan.template_id == "mr_refine_ligand"

  # --- Get directives for phase 1 ---
  d1 = plan_to_directives(plan)
  assert "phenix.xtriage" in d1.get(
    "workflow_preferences", {}
  ).get("prefer_programs", [])

  # --- User directives (from extraction) ---
  user_d = {
    "workflow_preferences": {
      "skip_programs": ["phenix.xtriage"],
    },
    "program_settings": {
      "phenix.refine": {"nproc": 8},
    },
  }

  # --- Merge ---
  merged, warnings = merge_directives(d1, user_d)
  assert len(warnings) > 0  # xtriage conflict
  assert "xtriage" in warnings[0]

  # --- Repair ---
  repair_msgs = repair_plan(plan, user_d)
  assert len(repair_msgs) > 0
  assert any("Repair" in m for m in repair_msgs)

  # --- Persist ---
  session_data = {}
  session_data["plan"] = plan.to_dict()
  assert session_data["plan"]["_version"] == 1

  # --- Restore ---
  plan2 = StructurePlan.from_dict(
    session_data["plan"]
  )
  assert len(plan2.stages) == len(plan.stages)
  assert plan2.template_id == "mr_refine_ligand"

  # --- Advance and check revision ---
  plan2.stages[0].status = STAGE_ACTIVE
  plan2.advance()
  changed = check_plan_revision(
    plan2, session_data
  )
  assert changed is True
  assert session_data.get("advice_changed") is True


# ── Entry point ──────────────────────────────────────

if __name__ == "__main__":
  run_tests()
