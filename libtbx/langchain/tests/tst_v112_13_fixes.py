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
import glob
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

def _discover_companion_files(available_files):
    """Inline copy of graph_nodes._discover_companion_files for testing."""
    seen = {os.path.basename(f) for f in available_files if f}
    new_files = []

    for f in available_files:
        if not f:
            continue
        basename = os.path.basename(f)

        if basename.endswith("_data.mtz") and "refine" in basename.lower():
            prefix = basename[:-9]
            output_dir = os.path.dirname(os.path.abspath(f))

            for companion in [prefix + ".mtz", prefix + "_001.mtz"]:
                full_path = os.path.join(output_dir, companion)
                if companion not in seen and os.path.exists(full_path):
                    new_files.append(full_path)
                    seen.add(companion)
                    break

            for companion in [prefix + ".pdb", prefix + "_001.pdb"]:
                full_path = os.path.join(output_dir, companion)
                if companion not in seen and os.path.exists(full_path):
                    new_files.append(full_path)
                    seen.add(companion)
                    break

            found_basenames = {os.path.basename(x) for x in new_files}
            for mtz in glob.glob(os.path.join(output_dir, prefix + "*.mtz")):
                bn = os.path.basename(mtz)
                if (bn not in seen and
                    bn not in found_basenames and
                    not bn.endswith("_data.mtz")):
                    new_files.append(mtz)
                    seen.add(bn)

        elif "overall_best" in basename.lower() and basename.endswith(".mtz"):
            output_dir = os.path.dirname(os.path.abspath(f))
            for companion in ["overall_best.pdb", "overall_best_refine_001.pdb"]:
                full_path = os.path.join(output_dir, companion)
                if companion not in seen and os.path.exists(full_path):
                    new_files.append(full_path)
                    seen.add(companion)
                    break

    # pdbtools output discovery
    agent_dirs = set()
    for f in available_files:
        if not f:
            continue
        abs_f = os.path.abspath(f)
        parts = abs_f.replace("\\", "/").split("/")
        for i, part in enumerate(parts):
            if part.startswith("sub_") and "_" in part[4:]:
                agent_dir = os.sep.join(parts[:i])
                if os.path.isdir(agent_dir):
                    agent_dirs.add(agent_dir)
                break

    for agent_dir in agent_dirs:
        for entry in glob.glob(os.path.join(agent_dir, "sub_*_pdbtools")):
            if os.path.isdir(entry):
                for pdb in glob.glob(os.path.join(entry, "*_with_ligand.pdb")):
                    bn = os.path.basename(pdb)
                    if bn not in seen:
                        new_files.append(pdb)
                        seen.add(bn)

    if new_files:
        return list(available_files) + new_files
    return available_files


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
# 1. COMPANION FILE DISCOVERY
# =============================================================================

def test_discover_refine_companions():
    """After refinement, map coefficients MTZ and refined PDB should be discovered
    even if the client only tracked _data.mtz and .cif."""
    print("Test: discover_refine_companions")

    tmpdir = tempfile.mkdtemp()
    try:
        refine_dir = os.path.join(tmpdir, "sub_02_refine")
        os.makedirs(refine_dir)

        for name in ["refine_001.pdb", "refine_001.mtz",
                     "refine_001_data.mtz", "refine_001_001.cif"]:
            with open(os.path.join(refine_dir, name), "w") as f:
                f.write("test")

        available = [
            os.path.join(refine_dir, "refine_001_data.mtz"),
            os.path.join(refine_dir, "refine_001_001.cif"),
        ]

        result = _discover_companion_files(available)
        basenames = [os.path.basename(f) for f in result]

        assert_in("refine_001.mtz", basenames,
                  "Map coefficients MTZ should be discovered")
        assert_in("refine_001.pdb", basenames,
                  "Refined model PDB should be discovered")
        assert_equal(len(result), 4,
                    "Should have 4 files total (2 original + 2 discovered)")
    finally:
        shutil.rmtree(tmpdir)

    print("  PASSED")


def test_discover_refine_001_format():
    """Discover companion files with _001 naming (refine_001_001.mtz format)."""
    print("Test: discover_refine_001_format")

    tmpdir = tempfile.mkdtemp()
    try:
        refine_dir = os.path.join(tmpdir, "sub_02_refine")
        os.makedirs(refine_dir)

        for name in ["refine_001_001.pdb", "refine_001_001.mtz",
                     "refine_001_data.mtz"]:
            with open(os.path.join(refine_dir, name), "w") as f:
                f.write("test")

        available = [os.path.join(refine_dir, "refine_001_data.mtz")]
        result = _discover_companion_files(available)
        basenames = [os.path.basename(f) for f in result]

        assert_in("refine_001_001.mtz", basenames,
                  "_001 format map coefficients should be discovered")
        assert_in("refine_001_001.pdb", basenames,
                  "_001 format model should be discovered")
    finally:
        shutil.rmtree(tmpdir)

    print("  PASSED")


def test_discover_autobuild_model():
    """Discover overall_best.pdb when only autobuild MTZ files are tracked."""
    print("Test: discover_autobuild_model")

    tmpdir = tempfile.mkdtemp()
    try:
        ab_dir = os.path.join(tmpdir, "AutoBuild_run_1_")
        os.makedirs(ab_dir)

        for name in ["overall_best.pdb", "overall_best_refine_data.mtz",
                     "overall_best_denmod_map_coeffs.mtz"]:
            with open(os.path.join(ab_dir, name), "w") as f:
                f.write("test")

        available = [
            os.path.join(ab_dir, "overall_best_refine_data.mtz"),
            os.path.join(ab_dir, "overall_best_denmod_map_coeffs.mtz"),
        ]

        result = _discover_companion_files(available)
        basenames = [os.path.basename(f) for f in result]
        assert_in("overall_best.pdb", basenames,
                  "Autobuild model PDB should be discovered")
    finally:
        shutil.rmtree(tmpdir)

    print("  PASSED")


def test_discover_pdbtools_output():
    """Discover *_with_ligand.pdb from pdbtools output directory."""
    print("Test: discover_pdbtools_output")

    tmpdir = tempfile.mkdtemp()
    try:
        agent_dir = os.path.join(tmpdir, "ai_agent_directory")
        refine_dir = os.path.join(agent_dir, "sub_02_refine")
        pdbtools_dir = os.path.join(agent_dir, "sub_04_pdbtools")
        os.makedirs(refine_dir)
        os.makedirs(pdbtools_dir)

        with open(os.path.join(refine_dir, "refine_001_data.mtz"), "w") as f:
            f.write("test")
        with open(os.path.join(pdbtools_dir, "refine_001_001_with_ligand.pdb"), "w") as f:
            f.write("test")

        available = [os.path.join(refine_dir, "refine_001_data.mtz")]
        result = _discover_companion_files(available)

        basenames = [os.path.basename(f) for f in result]
        assert_in("refine_001_001_with_ligand.pdb", basenames,
                  "pdbtools with_ligand output should be discovered")
    finally:
        shutil.rmtree(tmpdir)

    print("  PASSED")


def test_discover_no_duplicates():
    """Discovery should not add files that are already in the list."""
    print("Test: discover_no_duplicates")

    tmpdir = tempfile.mkdtemp()
    try:
        refine_dir = os.path.join(tmpdir, "sub_02_refine")
        os.makedirs(refine_dir)

        for name in ["refine_001.pdb", "refine_001.mtz", "refine_001_data.mtz"]:
            with open(os.path.join(refine_dir, name), "w") as f:
                f.write("test")

        available = [
            os.path.join(refine_dir, "refine_001_data.mtz"),
            os.path.join(refine_dir, "refine_001.mtz"),
            os.path.join(refine_dir, "refine_001.pdb"),
        ]

        result = _discover_companion_files(available)
        assert_equal(len(result), 3, "No new files should be added")
    finally:
        shutil.rmtree(tmpdir)

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
    phases = xray.get("phases", {})

    # Search all phases for phaser with model_for_mr condition
    found_model_for_mr = False
    for phase_name, phase_def in phases.items():
        programs = phase_def.get("programs", [])
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
    """The combine_ligand phase should be defined in xray workflow."""
    print("Test: combine_ligand_phase_exists")

    workflows = _load_workflows()
    xray = workflows.get("xray", {})
    phases = xray.get("phases", {})

    assert_in("combine_ligand", phases,
             "combine_ligand phase should exist in xray workflow")

    cl_phase = phases["combine_ligand"]
    programs = cl_phase.get("programs", [])

    prog_names = []
    for p in programs:
        if isinstance(p, str):
            prog_names.append(p)
        elif isinstance(p, dict):
            prog_names.extend(p.keys())

    assert_in("phenix.pdbtools", prog_names,
             "combine_ligand phase should include phenix.pdbtools")

    transitions = cl_phase.get("transitions", {})
    assert_equal(transitions.get("on_complete"), "refine",
                "combine_ligand should transition to refine on complete")

    print("  PASSED")


# =============================================================================
# RUN ALL TESTS
# =============================================================================

def run_all_tests():
    """Run all v112.13 tests."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
