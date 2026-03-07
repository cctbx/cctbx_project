"""
Tests for v112.13 fixes:
- Companion file discovery (refine, autobuild, pdbtools)
- Intermediate file filtering (TEMP0/, EDITED_*)
- Output validation (output.file_name excluded from validation)
- File categorization (with_ligand, EDITED exclusion, refine_map_coeffs)
- has_model_for_mr phaser condition
- Program input priorities (pdbtools, refine)
- Workflow YAML conditions (phaser model_for_mr)
"""

from __future__ import absolute_import, division, print_function
import os
import sys
import shutil
import tempfile
import re
import yaml

# Ensure test utilities are available
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tests.tst_utils import (
    assert_true, assert_equal, assert_in,
    run_tests_with_fail_fast
)

# Paths to knowledge files
KNOWLEDGE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "knowledge")


def _load_categories():
    with open(os.path.join(KNOWLEDGE_DIR, "file_categories.yaml")) as f:
        return yaml.safe_load(f)


def _load_programs():
    with open(os.path.join(KNOWLEDGE_DIR, "programs.yaml")) as f:
        return yaml.safe_load(f)


def _load_workflows():
    with open(os.path.join(KNOWLEDGE_DIR, "workflows.yaml")) as f:
        return yaml.safe_load(f)


# ---- Inline reimplementations of tested functions ----
# (graph_nodes.py can't be imported without libtbx)

def _filter_intermediate_files(available_files):
    """Inline copy of graph_nodes._filter_intermediate_files for testing."""
    temp_dir_markers = ("/TEMP", "/temp", "/TEMP0/", "/scratch/")
    temp_file_prefixes = ("EDITED_", "TEMP_")

    filtered = []
    for f in available_files:
        if not f:
            continue
        if any(marker in f for marker in temp_dir_markers):
            continue
        basename = os.path.basename(f)
        if any(basename.startswith(prefix) for prefix in temp_file_prefixes):
            continue
        filtered.append(f)

    return filtered


# =============================================================================
# 1. PERCEIVE DOES NOT SCAN INPUT DIRECTORIES (v112.79)
# =============================================================================

def test_perceive_no_input_dir_scanning():
  """Verify that the perceive node's file pipeline does NOT scan
  directories of user-supplied files.

  The old _discover_companion_files() would look in the directory of
  every available file.  If a user supplied /data/project/data.mtz,
  and /data/project/ also contained mask.ccp4, the agent could pick
  up mask.ccp4 even though the user never supplied it.

  After v112.79, the perceive node only:
    - Injects output files from history entries
    - Filters intermediate/temp files
  No directory scanning occurs in the graph pipeline.  Companion
  file discovery is handled on the client side by
  session.get_available_files().
  """
  print("Test: perceive_no_input_dir_scanning")

  tmpdir = tempfile.mkdtemp()
  try:
    # Create a user input directory with extra files the agent
    # should NOT discover
    user_dir = os.path.join(tmpdir, "user_data")
    os.makedirs(user_dir)

    # User-supplied file
    data_mtz = os.path.join(user_dir, "data.mtz")
    with open(data_mtz, "w") as f:
      f.write("test")

    # Extra files in the same directory — NOT supplied by user
    for extra in ["mask.ccp4", "old_model.pdb",
                  "refine_001.mtz", "refine_001_data.mtz"]:
      with open(os.path.join(user_dir, extra), "w") as f:
        f.write("test")

    # Simulate what perceive does: start with available_files,
    # inject from history, then filter intermediates.
    # NO directory scanning should occur.
    available_files = [data_mtz]

    # Simulate history injection (as perceive does)
    history = [
      {"output_files": [
        os.path.join(tmpdir, "sub_01_xtriage", "xtriage.log")
      ]},
    ]
    seen_basenames = {
      os.path.basename(f) for f in available_files if f
    }
    for hist_entry in history:
      if isinstance(hist_entry, dict):
        for f in hist_entry.get("output_files", []):
          if f:
            bn = os.path.basename(f)
            if bn not in seen_basenames:
              available_files = list(available_files) + [f]
              seen_basenames.add(bn)

    # Apply intermediate file filter (as perceive does)
    available_files = _filter_intermediate_files(available_files)

    # Verify: only user-supplied + history files are present
    basenames = {os.path.basename(f) for f in available_files}
    assert_true(
      "data.mtz" in basenames,
      "User-supplied file should be present")
    assert_true(
      "xtriage.log" in basenames,
      "History output file should be present")

    # The key assertion: files from the input directory that
    # were NOT supplied or in history must NOT appear
    assert_true(
      "mask.ccp4" not in basenames,
      "mask.ccp4 from input dir should NOT be discovered")
    assert_true(
      "old_model.pdb" not in basenames,
      "old_model.pdb from input dir should NOT be discovered")
    assert_true(
      "refine_001.mtz" not in basenames,
      "refine_001.mtz from input dir should NOT be discovered")
    assert_true(
      "refine_001_data.mtz" not in basenames,
      "refine_001_data.mtz from input dir should NOT "
      "be discovered")
    assert_equal(
      len(available_files), 2,
      "Should have exactly 2 files (1 user + 1 history)")
  finally:
    shutil.rmtree(tmpdir)

  print("  PASSED")


def test_history_injection_still_works():
  """Verify that output files from history entries are still
  injected into available_files by the perceive pipeline."""
  print("Test: history_injection_still_works")

  # Simulate perceive's history injection logic
  available_files = ["/p/data.mtz", "/p/sequence.fa"]
  history = [
    {"output_files": []},  # xtriage: no output files
    {"output_files": [     # refine: has outputs
      "/p/sub_02_refine/refine_001.pdb",
      "/p/sub_02_refine/refine_001_data.mtz",
    ]},
  ]

  seen_basenames = {
    os.path.basename(f) for f in available_files if f
  }
  history_additions = []
  for hist_entry in history:
    if isinstance(hist_entry, dict):
      for f in hist_entry.get("output_files", []):
        if f:
          bn = os.path.basename(f)
          if bn not in seen_basenames:
            history_additions.append(f)
            seen_basenames.add(bn)
  if history_additions:
    available_files = list(available_files) + history_additions

  basenames = {os.path.basename(f) for f in available_files}
  assert_equal(len(available_files), 4,
    "Should have 4 files (2 original + 2 from history)")
  assert_in("refine_001.pdb", basenames,
    "Refine model should be injected from history")
  assert_in("refine_001_data.mtz", basenames,
    "Refine data MTZ should be injected from history")

  print("  PASSED")


# =============================================================================
# 2. INTERMEDIATE FILE FILTERING
# =============================================================================

def test_filter_temp_directory_files():
    """Files in TEMP0/ directories should be filtered out."""
    print("Test: filter_temp_directory_files")

    files = [
        "/p/7qz0.pdb",
        "/p/ligandfit/LigandFit_run_1_/ligand_fit_1.pdb",
        "/p/ligandfit/LigandFit_run_1_/TEMP0/EDITED_7qz0.pdb",
        "/p/ligandfit/LigandFit_run_1_/TEMP0/other.pdb",
    ]

    result = _filter_intermediate_files(files)
    basenames = [os.path.basename(f) for f in result]

    assert_in("7qz0.pdb", basenames, "Original file should be kept")
    assert_in("ligand_fit_1.pdb", basenames, "Ligandfit output should be kept")
    assert_true("EDITED_7qz0.pdb" not in basenames,
               "TEMP0 files should be filtered")
    assert_true("other.pdb" not in basenames,
               "Files in TEMP0 dir should be filtered")
    assert_equal(len(result), 2, "Only 2 files should remain")

    print("  PASSED")


def test_filter_edited_prefix_files():
    """Files with EDITED_ prefix should be filtered even outside TEMP dirs."""
    print("Test: filter_edited_prefix_files")

    files = [
        "/p/7qz0.pdb",
        "/p/EDITED_7qz0.pdb",
        "/p/TEMP_scratch.pdb",
    ]

    result = _filter_intermediate_files(files)
    assert_equal(len(result), 1, "Only original PDB should remain")
    assert_equal(os.path.basename(result[0]), "7qz0.pdb")

    print("  PASSED")


def test_filter_preserves_legitimate_files():
    """Legitimate files should not be filtered."""
    print("Test: filter_preserves_legitimate_files")

    files = [
        "/p/refine_001.pdb",
        "/p/ligand_fit_1.pdb",
        "/p/overall_best.pdb",
        "/p/7qz0_with_ligand.pdb",
        "/p/sequence.fa",
        "/p/data.mtz",
    ]

    result = _filter_intermediate_files(files)
    assert_equal(len(result), len(files), "No legitimate files should be filtered")

    print("  PASSED")


# =============================================================================
# 3. FILE CATEGORIZATION
# =============================================================================

def test_with_ligand_categorization():
    """*_with_ligand.pdb should be in with_ligand -> model, also_in refined."""
    print("Test: with_ligand_categorization")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    files = ["/p/refine_001_001_with_ligand.pdb"]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    path = "/p/refine_001_001_with_ligand.pdb"
    assert_in(path, cats.get("with_ligand", []),
             "Should be in with_ligand subcategory")
    assert_in(path, cats.get("model", []),
             "Should bubble to model parent")
    assert_in(path, cats.get("refined", []),
             "Should be also_in refined (via with_ligand config)")

    print("  PASSED")


def test_ligand_fit_output_parent_is_ligand():
    """ligand_fit_1.pdb should be in ligand_fit_output -> ligand, NOT model."""
    print("Test: ligand_fit_output_parent_is_ligand")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    files = ["/p/ligand_fit_1.pdb"]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    path = "/p/ligand_fit_1.pdb"
    assert_in(path, cats.get("ligand_fit_output", []),
             "Should be in ligand_fit_output")
    assert_in(path, cats.get("ligand", []),
             "Should bubble to ligand parent")
    assert_true(path not in cats.get("model", []),
               "ligand_fit_1.pdb should NOT be in model category")

    print("  PASSED")


def test_edited_pdb_excluded_from_model():
    """EDITED_*.pdb files should be excluded from unclassified_pdb/model."""
    print("Test: edited_pdb_excluded_from_model")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    files = ["/p/EDITED_7qz0.pdb", "/p/7qz0.pdb"]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    model_basenames = [os.path.basename(f) for f in cats.get("model", [])]
    assert_true("EDITED_7qz0.pdb" not in model_basenames,
               "EDITED_7qz0.pdb should NOT be in model category")
    assert_in("7qz0.pdb", model_basenames,
             "7qz0.pdb should be in model category")

    print("  PASSED")


def test_refine_map_coeffs_bare_pattern():
    """refine_001.mtz (bare prefix) should be in refine_map_coeffs."""
    print("Test: refine_map_coeffs_bare_pattern")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    files = ["/p/refine_001.mtz", "/p/refine_001_data.mtz"]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    refine_mc = [os.path.basename(f) for f in cats.get("refine_map_coeffs", [])]
    map_mc = [os.path.basename(f) for f in cats.get("map_coeffs_mtz", [])]

    assert_in("refine_001.mtz", refine_mc,
              "refine_001.mtz should be in refine_map_coeffs")
    assert_in("refine_001.mtz", map_mc,
              "refine_001.mtz should bubble to map_coeffs_mtz")
    assert_true("refine_001_data.mtz" not in refine_mc,
               "refine_001_data.mtz should NOT be in refine_map_coeffs")

    print("  PASSED")


def test_refine_data_mtz_is_data_not_map_coeffs():
    """refine_001_data.mtz -> data_mtz, not map_coeffs_mtz."""
    print("Test: refine_data_mtz_is_data_not_map_coeffs")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    files = ["/p/refine_001_data.mtz"]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    assert_in("/p/refine_001_data.mtz", cats.get("data_mtz", []),
             "refine_001_data.mtz should be in data_mtz")
    assert_true("/p/refine_001_data.mtz" not in cats.get("map_coeffs_mtz", []),
               "refine_001_data.mtz should NOT be in map_coeffs_mtz")

    print("  PASSED")


# =============================================================================
# 4. PHASER model_for_mr CONDITION
# =============================================================================

def test_has_model_for_mr_with_search_model():
    """has_model_for_mr should be True when only a search_model is available."""
    print("Test: has_model_for_mr_with_search_model")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    files = ["/p/search_model.pdb", "/p/data.mtz"]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    has_search = bool(cats.get("search_model"))
    has_model_for_mr = bool(cats.get("model") or cats.get("search_model"))

    assert_true(has_search,
               "search_model.pdb should be in search_model category")
    assert_true(has_model_for_mr,
               "has_model_for_mr should be True (search_model OR model)")

    print("  PASSED")


def test_has_model_for_mr_with_generic_pdb():
    """has_model_for_mr should be True when a generic PDB is available."""
    print("Test: has_model_for_mr_with_generic_pdb")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    files = ["/p/7qz0.pdb", "/p/data.mtz"]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    has_model = bool(cats.get("model"))
    has_model_for_mr = bool(cats.get("model") or cats.get("search_model"))

    assert_true(has_model, "7qz0.pdb should be in model category")
    assert_true(has_model_for_mr, "has_model_for_mr should be True")

    print("  PASSED")


def test_phaser_condition_in_workflows_yaml():
    """workflows.yaml should use model_for_mr for phaser condition."""
    print("Test: phaser_condition_in_workflows_yaml")

    workflows = _load_workflows()
    xray = workflows.get("xray", {})
    steps = xray.get("steps") or xray.get("phases") or {}

    # Search all phases for phaser with model_for_mr condition
    found_model_for_mr = False
    for step_name, step_def in steps.items():
        programs = step_def.get("programs", [])
        for prog in programs:
            if isinstance(prog, dict):
                prog_name = prog.get("program", "")
                if "phaser" in str(prog_name).lower():
                    conditions = prog.get("conditions", [])
                    for cond in conditions:
                        if isinstance(cond, dict) and cond.get("has") == "model_for_mr":
                            found_model_for_mr = True

    assert_true(found_model_for_mr,
               "phaser condition should use 'has: model_for_mr' in workflows.yaml")

    print("  PASSED")


# =============================================================================
# 5. OUTPUT VALIDATION
# =============================================================================

def test_output_filename_excluded_from_validation():
    """output.file_name=X.pdb should not be validated as an input file."""
    print("Test: output_filename_excluded_from_validation")

    command = ("phenix.pdbtools /p/overall_best.pdb /p/ligand_fit_1.pdb "
               "output.file_name=overall_best_with_ligand.pdb")

    output_pattern = r'output\.\w+=\S+'
    command_for_validation = re.sub(output_pattern, '', command)

    file_pattern = r'[\w\-\.\/]+\.(?:pdb|mtz|fa|fasta|seq|cif|mrc|ccp4|map)\b'
    referenced_files = re.findall(file_pattern, command_for_validation, re.IGNORECASE)

    ref_basenames = [os.path.basename(f) for f in referenced_files]
    assert_true("overall_best_with_ligand.pdb" not in ref_basenames,
               "Output filename should NOT be extracted for validation")
    assert_in("overall_best.pdb", ref_basenames,
             "Input files should still be extracted")
    assert_in("ligand_fit_1.pdb", ref_basenames,
             "Input files should still be extracted")

    print("  PASSED")


def test_output_prefix_excluded_from_validation():
    """output.prefix=refine_002 should not cause validation issues."""
    print("Test: output_prefix_excluded_from_validation")

    command = ("phenix.refine /p/model.pdb /p/data.mtz "
               "output.prefix=refine_002 output.serial=1")

    output_pattern = r'output\.\w+=\S+'
    command_for_validation = re.sub(output_pattern, '', command)

    file_pattern = r'[\w\-\.\/]+\.(?:pdb|mtz|fa|fasta|seq|cif|mrc|ccp4|map)\b'
    referenced_files = re.findall(file_pattern, command_for_validation, re.IGNORECASE)

    assert_equal(len(referenced_files), 2, "Only 2 input files should be extracted")

    print("  PASSED")


def test_multiple_output_args_stripped():
    """Multiple output.X=Y arguments should all be stripped."""
    print("Test: multiple_output_args_stripped")

    command = ("phenix.refine /p/model.pdb /p/data.mtz "
               "output.prefix=refine_002 output.serial=1 "
               "output.file_name=result.pdb output.suffix=_final")

    output_pattern = r'output\.\w+=\S+'
    command_for_validation = re.sub(output_pattern, '', command)

    file_pattern = r'[\w\-\.\/]+\.(?:pdb|mtz|fa|fasta|seq|cif|mrc|ccp4|map)\b'
    referenced_files = re.findall(file_pattern, command_for_validation, re.IGNORECASE)

    ref_basenames = [os.path.basename(f) for f in referenced_files]
    assert_true("result.pdb" not in ref_basenames,
               "output.file_name value should be stripped")
    assert_equal(len(ref_basenames), 2, "Only 2 input files should remain")

    print("  PASSED")


# =============================================================================
# 6. PROGRAM INPUT PRIORITIES (YAML)
# =============================================================================

def test_pdbtools_protein_prefers_autobuild():
    """pdbtools protein slot should prefer autobuild_output over generic PDB."""
    print("Test: pdbtools_protein_prefers_autobuild")

    programs = _load_programs()
    pdbtools = programs["phenix.pdbtools"]
    protein_priorities = pdbtools["input_priorities"]["protein"]

    assert_in("autobuild_output", protein_priorities["prefer_subcategories"],
             "autobuild_output should be in pdbtools protein prefer_subcategories")
    assert_in("ligand_fit_output", protein_priorities["exclude_categories"],
             "ligand_fit_output should be excluded from pdbtools protein")

    print("  PASSED")


def test_pdbtools_excludes_edited_pattern():
    """pdbtools should exclude EDITED files from protein input."""
    print("Test: pdbtools_excludes_edited_pattern")

    programs = _load_programs()
    protein_def = programs["phenix.pdbtools"]["inputs"]["required"]["protein"]

    assert_in("EDITED", protein_def.get("exclude_patterns", []),
             "EDITED should be in pdbtools protein exclude_patterns")

    print("  PASSED")


def test_refine_excludes_ligand_categories():
    """phenix.refine should exclude ligand categories from model slot."""
    print("Test: refine_excludes_ligand_categories")

    programs = _load_programs()
    refine = programs["phenix.refine"]
    model_priorities = refine["input_priorities"]["model"]

    assert_in("ligand_fit_output", model_priorities["exclude_categories"],
             "ligand_fit_output should be excluded from refine model")
    assert_in("ligand", model_priorities["exclude_categories"],
             "ligand should be excluded from refine model")

    print("  PASSED")


def test_refine_prefers_with_ligand_first():
    """phenix.refine should prefer with_ligand models first."""
    print("Test: refine_prefers_with_ligand_first")

    programs = _load_programs()
    refine = programs["phenix.refine"]
    model_priorities = refine["input_priorities"]["model"]

    assert_equal(model_priorities["prefer_subcategories"][0], "with_ligand",
                "with_ligand should be first in refine model prefer_subcategories")

    print("  PASSED")


# =============================================================================
# 7. END-TO-END: POST-PDBTOOLS REFINEMENT FILE SELECTION
# =============================================================================

def test_post_pdbtools_refine_selects_with_ligand():
    """After pdbtools, refinement should pick *_with_ligand.pdb, not ligand_fit_1.pdb."""
    print("Test: post_pdbtools_refine_selects_with_ligand")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    programs = _load_programs()

    files = [
        "/p/7qz0.pdb", "/p/7qz0.mtz", "/p/7qz0.fa", "/p/7qz0_ligand.pdb",
        "/p/refine_001_data.mtz", "/p/refine_001.mtz", "/p/refine_001.pdb",
        "/p/refine_001_001.cif",
        "/p/ligand_fit_1.pdb",
        "/p/refine_001_001_with_ligand.pdb",
    ]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    refine_model_priorities = programs["phenix.refine"]["input_priorities"]["model"]
    exclude_cats = refine_model_priorities["exclude_categories"]
    prefer_subcats = refine_model_priorities["prefer_subcategories"]

    excluded_files = set()
    for cat in exclude_cats:
        for f in cats.get(cat, []):
            excluded_files.add(f)

    assert_true("/p/ligand_fit_1.pdb" in excluded_files,
               "ligand_fit_1.pdb should be excluded from refine model")

    selected = None
    for subcat in prefer_subcats:
        for f in cats.get(subcat, []):
            if f not in excluded_files:
                selected = f
                break
        if selected:
            break

    assert_equal(os.path.basename(selected), "refine_001_001_with_ligand.pdb",
                "Refine should select the with_ligand model, got: %s" %
                (os.path.basename(selected) if selected else "None"))

    print("  PASSED")


def test_best_model_ligand_fit_excluded_from_refine():
    """If best_files points to ligand_fit_1.pdb, refine should NOT use it."""
    print("Test: best_model_ligand_fit_excluded_from_refine")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    programs = _load_programs()

    files = ["/p/ligand_fit_1.pdb", "/p/refine_001_001_with_ligand.pdb"]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    best_model = "/p/ligand_fit_1.pdb"
    best_bn = os.path.basename(best_model)

    refine_model_priorities = programs["phenix.refine"]["input_priorities"]["model"]
    exclude_cats = refine_model_priorities["exclude_categories"]

    is_excluded = False
    for cat in exclude_cats:
        cat_files = cats.get(cat, [])
        if any(os.path.basename(f) == best_bn for f in cat_files):
            is_excluded = True
            break

    assert_true(is_excluded,
               "best_model ligand_fit_1.pdb should be detected as excluded for refine")

    print("  PASSED")


# =============================================================================
# 8. COMBINE LIGAND PHASE
# =============================================================================

def test_combine_ligand_phase_exists():
    """The combine_ligand step should be defined in xray workflow."""
    print("Test: combine_ligand_phase_exists")

    workflows = _load_workflows()
    xray = workflows.get("xray", {})
    steps = xray.get("steps") or xray.get("phases") or {}

    assert_in("combine_ligand", steps,
             "combine_ligand step should exist in xray workflow")

    cl_step = steps["combine_ligand"]
    programs = cl_step.get("programs", [])

    prog_names = []
    for p in programs:
        if isinstance(p, str):
            prog_names.append(p)
        elif isinstance(p, dict):
            prog_names.extend(p.keys())

    assert_in("phenix.pdbtools", prog_names,
             "combine_ligand step should include phenix.pdbtools")

    transitions = cl_step.get("transitions", {})
    assert_equal(transitions.get("on_complete"), "refine",
                "combine_ligand should transition to refine on complete")

    print("  PASSED")


# =============================================================================
# 9. SHARPENED MAP CATEGORIZATION
# =============================================================================

def test_sharpened_half_map_not_in_half_map():
    """Sharpened maps with 'half' in name should NOT be in half_map category."""
    print("Test: sharpened_half_map_not_in_half_map")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    files = [
        "/p/emd_20026_half_map_2_box_sharpened.ccp4",  # Sharpened output
        "/p/emd_20026_half_map_1.ccp4",                 # Original half map 1
        "/p/emd_20026_half_map_2.ccp4",                 # Original half map 2
    ]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    half_basenames = [os.path.basename(f) for f in cats.get("half_map", [])]
    map_basenames = [os.path.basename(f) for f in cats.get("map", [])]

    # Sharpened file should NOT be in half_map
    assert_true("emd_20026_half_map_2_box_sharpened.ccp4" not in half_basenames,
               "Sharpened map should NOT be in half_map")

    # But it should be in map
    assert_in("emd_20026_half_map_2_box_sharpened.ccp4", map_basenames,
             "Sharpened map should be in map category")

    # Original half maps should still be in half_map
    assert_in("emd_20026_half_map_1.ccp4", half_basenames,
             "Original half map 1 should be in half_map")
    assert_in("emd_20026_half_map_2.ccp4", half_basenames,
             "Original half map 2 should be in half_map")

    print("  PASSED")


def test_mtriage_finds_sharpened_map_as_full_map():
    """mtriage input priorities should pick sharpened map for full_map slot."""
    print("Test: mtriage_finds_sharpened_map_as_full_map")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    programs = _load_programs()

    files = [
        "/p/emd_20026_half_map_2_box_sharpened.ccp4",
        "/p/emd_20026_half_map_1.ccp4",
        "/p/emd_20026_half_map_2.ccp4",
    ]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    # Simulate mtriage full_map input selection
    mtriage = programs["phenix.mtriage"]
    full_map_prio = mtriage["input_priorities"]["full_map"]
    search_cats = full_map_prio["categories"]      # [full_map, map]
    exclude_cats = full_map_prio["exclude_categories"]  # [half_map]

    # Build excluded file set
    excluded = set()
    for ecat in exclude_cats:
        for f in cats.get(ecat, []):
            excluded.add(f)

    # Search through categories in priority order
    selected = None
    for cat in search_cats:
        for f in cats.get(cat, []):
            if f not in excluded:
                selected = f
                break
        if selected:
            break

    assert_true(selected is not None, "Should find a map for full_map slot")
    assert_equal(os.path.basename(selected), "emd_20026_half_map_2_box_sharpened.ccp4",
                "mtriage should select sharpened map for full_map, got: %s" %
                (os.path.basename(selected) if selected else "None"))

    print("  PASSED")


# =============================================================================
# 10. DIRECTIVE DONE-FLAG CHECKS (YAML CONFIG)
# =============================================================================

def test_run_once_programs_have_done_flags():
    """All run_once programs should have a flag defined for done checking."""
    print("Test: run_once_programs_have_done_flags")

    programs = _load_programs()

    for prog_name, prog_def in programs.items():
        tracking = prog_def.get("done_tracking", {})
        if tracking.get("run_once") or tracking.get("strategy") == "run_once":
            flag = tracking.get("flag", "")
            assert_true(bool(flag),
                       "%s has run_once but no flag defined" % prog_name)

    print("  PASSED")


def test_start_with_program_respects_done_flags_config():
    """The start_with_program directive path should be able to detect done programs.
    Verify that map_symmetry has the necessary done_tracking config."""
    print("Test: start_with_program_respects_done_flags_config")

    programs = _load_programs()

    # These programs can be requested via start_with_program directive
    # and should have run_once tracking to prevent re-running
    run_once_programs = ["phenix.map_symmetry", "phenix.mtriage",
                        "phenix.xtriage"]

    for prog_name in run_once_programs:
        prog_def = programs.get(prog_name, {})
        tracking = prog_def.get("done_tracking", {})
        has_run_once = (tracking.get("run_once") or
                       tracking.get("strategy") == "run_once")
        flag = tracking.get("flag", "")

        assert_true(has_run_once,
                   "%s should have run_once tracking" % prog_name)
        assert_true(bool(flag),
                   "%s should have a done flag" % prog_name)
        assert_true(flag.endswith("_done"),
                   "%s done flag should end with '_done', got '%s'" % (prog_name, flag))

    print("  PASSED")


# =============================================================================
# 11. MAP_SYMMETRY REQUIRES NON-HALF MAP
# =============================================================================

def test_map_symmetry_requires_non_half_map():
    """map_symmetry in workflows.yaml should require has: non_half_map."""
    print("Test: map_symmetry_requires_non_half_map")

    workflows = _load_workflows()
    cryoem = workflows.get("cryoem", {})
    steps = cryoem.get("steps") or cryoem.get("phases") or {}

    # Find map_symmetry in any step
    found_condition = False
    for step_name, step_def in steps.items():
        programs = step_def.get("programs", [])
        for prog in programs:
            if isinstance(prog, dict) and prog.get("program") == "phenix.map_symmetry":
                conditions = prog.get("conditions", [])
                for cond in conditions:
                    if isinstance(cond, dict) and cond.get("has") == "non_half_map":
                        found_condition = True

    assert_true(found_condition,
               "map_symmetry should have 'has: non_half_map' condition in workflows.yaml")

    print("  PASSED")


def test_non_half_map_with_only_half_maps():
    """has_non_half_map should be False when only half maps are available."""
    print("Test: non_half_map_with_only_half_maps")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    files = [
        "/p/emd_20026_half_map_1_box.ccp4",
        "/p/emd_20026_half_map_2_box.ccp4",
    ]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    map_set = set(cats.get("map", []))
    half_set = set(cats.get("half_map", []))
    has_non_half = bool(map_set - half_set)

    assert_true(not has_non_half,
               "has_non_half_map should be False with only half maps")

    print("  PASSED")


def test_non_half_map_with_sharpened_map():
    """has_non_half_map should be True when a sharpened map is available."""
    print("Test: non_half_map_with_sharpened_map")

    from agent.workflow_state import _categorize_files_yaml, _bubble_up_to_parents

    rules = _load_categories()
    files = [
        "/p/emd_20026_half_map_1_box.ccp4",
        "/p/emd_20026_half_map_2_box.ccp4",
        "/p/emd_20026_half_map_2_box_sharpened.ccp4",  # Sharpened output
    ]

    cats = _categorize_files_yaml(files, rules)
    cats = _bubble_up_to_parents(cats, rules)

    map_set = set(cats.get("map", []))
    half_set = set(cats.get("half_map", []))
    has_non_half = bool(map_set - half_set)

    assert_true(has_non_half,
               "has_non_half_map should be True when sharpened map is present")

    # Verify the sharpened file is what makes the difference
    non_half = map_set - half_set
    non_half_basenames = [os.path.basename(f) for f in non_half]
    assert_in("emd_20026_half_map_2_box_sharpened.ccp4", non_half_basenames,
             "Sharpened map should be the non-half map")

    print("  PASSED")


def test_pdbtools_modified_recognized_as_with_ligand():
    """Pdbtools *_modified.pdb must be treated as with_ligand, not refined.

    Regression: After ligandfit + pdbtools combination, the output file
    (e.g. refine_001_001_modified.pdb) was classified as 'refined' because
    it contained 'refine' in the name.  This caused best_model to stay on
    the old ligand-free refined model, and the command builder would
    override the LLM's correct file choice.
    """
    print("Test: pdbtools_modified_recognized_as_with_ligand")

    from agent.best_files_tracker import BestFilesTracker

    tracker = BestFilesTracker()

    # First, register the refined model (no ligand) with good metrics
    tracker.evaluate_file(
        "/p/refine_001_001.pdb", cycle=3,
        metrics={"r_free": 0.28}, stage="refined")

    # Then register the pdbtools combined file (model+ligand)
    # This file has no metrics yet — just came from pdbtools
    tracker.evaluate_file(
        "/p/refine_001_001_modified.pdb", cycle=5,
        metrics={}, stage="with_ligand")

    best = tracker.get_best_dict()
    best_model = best.get("model", "")

    assert_in("modified", os.path.basename(best_model),
             "best model should be the _modified (with_ligand) file, "
             "got: %s" % best_model)

    # Also verify _classify_stage recognizes _modified as with_ligand
    stage = tracker._classify_stage(
        "/p/refine_001_001_modified.pdb", "model")
    assert_equal("with_ligand", stage,
                "classify_stage should return with_ligand for _modified.pdb")

    print("  PASSED")


def test_pdbtools_modified_not_overridden_by_best_model():
    """Command builder should NOT override LLM's _modified.pdb choice.

    When the LLM correctly picks the pdbtools combined model+ligand file,
    the build step must not replace it with the old best_model.
    """
    print("Test: pdbtools_modified_not_overridden_by_best_model")

    from agent.best_files_tracker import BestFilesTracker

    tracker = BestFilesTracker()

    # Register refined model
    tracker.evaluate_file(
        "/p/refine_001_001.pdb", cycle=3,
        metrics={"r_free": 0.28}, stage="refined")

    # Register _modified with explicit with_ligand stage (as session.py
    # would set it after pdbtools via _infer_stage_from_program)
    tracker.evaluate_file(
        "/p/refine_001_001_modified.pdb", cycle=5,
        metrics={}, stage="with_ligand")

    best = tracker.get_best_dict()
    best_model = os.path.basename(best.get("model", ""))

    # The with_ligand file should have higher stage score (110 > 100)
    # even without metrics, so it should be the best model
    assert_equal("refine_001_001_modified.pdb", best_model,
                "with_ligand (110) should beat refined (100+metric) — "
                "got: %s" % best_model)

    print("  PASSED")


# =============================================================================
# RUN ALL TESTS
# =============================================================================

def run_all_tests():
    """Run all v112.13 tests."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
