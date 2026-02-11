"""
Tests for File Categorization Logic.

This tests the critical file categorization that determines:
- Which files are PDBs, MTZs, maps, sequences, etc.
- Which files are refined, predicted, phaser output, etc.
- Priority ordering for file selection

Run with: python tests/test_file_categorization.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Mock libtbx imports
import types
libtbx = types.ModuleType('libtbx')
libtbx.langchain = types.ModuleType('libtbx.langchain')
libtbx.langchain.agent = types.ModuleType('libtbx.langchain.agent')
libtbx.langchain.knowledge = types.ModuleType('libtbx.langchain.knowledge')
sys.modules['libtbx'] = libtbx
sys.modules['libtbx.langchain'] = libtbx.langchain
sys.modules['libtbx.langchain.agent'] = libtbx.langchain.agent
sys.modules['libtbx.langchain.knowledge'] = libtbx.langchain.knowledge

import knowledge.yaml_loader
libtbx.langchain.knowledge.yaml_loader = knowledge.yaml_loader
sys.modules['libtbx.langchain.knowledge.yaml_loader'] = knowledge.yaml_loader

# Import file categorization
from agent.best_files_tracker import BestFilesTracker


# =============================================================================
# TEST UTILITIES
# =============================================================================

def assert_equal(a, b, msg=""):
    if a != b:
        raise AssertionError("%s: %r != %r" % (msg, a, b) if msg else "%r != %r" % (a, b))


def assert_true(condition, msg="Assertion failed"):
    if not condition:
        raise AssertionError(msg)


def assert_in(item, container, msg=""):
    if item not in container:
        raise AssertionError("%s: %r not in %r" % (msg, item, container) if msg else "%r not in %r" % (item, container))


def assert_not_in(item, container, msg=""):
    if item in container:
        raise AssertionError("%s: %r should not be in %r" % (msg, item, container) if msg else "%r should not be in %r" % (item, container))


def assert_not_none(value, msg=""):
    if value is None:
        raise AssertionError(msg if msg else "Expected non-None value")


# =============================================================================
# FILE EXTENSION TESTS
# =============================================================================

def test_pdb_file_detection():
    """Test detection of PDB/model files."""
    print("Test: pdb_file_detection")

    # Various PDB file patterns that should be recognized
    pdb_files = [
        "/path/to/model.pdb",
        "/path/to/refined.pdb",
        "/path/to/PHASER.1.pdb",
        "/path/to/structure.cif",
        "/path/to/refine_001_001.pdb",
        "/path/to/rsr_001_real_space_refined_000.pdb",
    ]

    for f in pdb_files:
        ext = os.path.splitext(f)[1].lower()
        assert_in(ext, ['.pdb', '.cif'], "Should recognize as model: %s" % f)

    print("  PASSED")


def test_mtz_file_detection():
    """Test detection of MTZ/reflection files."""
    print("Test: mtz_file_detection")

    mtz_files = [
        "/path/to/data.mtz",
        "/path/to/refine_001_001.mtz",
        "/path/to/scaled_data.mtz",
    ]

    for f in mtz_files:
        ext = os.path.splitext(f)[1].lower()
        assert_equal(ext, '.mtz', "Should recognize as MTZ: %s" % f)

    print("  PASSED")


def test_map_file_detection():
    """Test detection of map files."""
    print("Test: map_file_detection")

    map_files = [
        "/path/to/map.mrc",
        "/path/to/map.ccp4",
        "/path/to/density.map",
        "/path/to/emd_12345.map",
    ]

    map_extensions = ['.mrc', '.ccp4', '.map']

    for f in map_files:
        ext = os.path.splitext(f)[1].lower()
        assert_in(ext, map_extensions, "Should recognize as map: %s" % f)

    print("  PASSED")


def test_sequence_file_detection():
    """Test detection of sequence files."""
    print("Test: sequence_file_detection")

    seq_files = [
        "/path/to/sequence.fa",
        "/path/to/protein.fasta",
        "/path/to/seq.seq",
    ]

    seq_extensions = ['.fa', '.fasta', '.seq']

    for f in seq_files:
        ext = os.path.splitext(f)[1].lower()
        assert_in(ext, seq_extensions, "Should recognize as sequence: %s" % f)

    print("  PASSED")


# =============================================================================
# FILE STAGE CLASSIFICATION TESTS
# =============================================================================

def test_refined_model_detection():
    """Test detection of refined model files."""
    print("Test: refined_model_detection")

    tracker = BestFilesTracker()

    # Files that should be classified as "refined"
    refined_patterns = [
        "refine_001_001.pdb",
        "refine_002_001.pdb",
        "refined_model.pdb",
        "model_refined.pdb",
    ]

    for filename in refined_patterns:
        stage = tracker._classify_stage("/path/to/" + filename, "model")
        assert_equal(stage, "refined", "Should classify as refined: %s (got %s)" % (filename, stage))

    print("  PASSED")


def test_rsr_output_detection():
    """Test detection of real_space_refine output files."""
    print("Test: rsr_output_detection")

    tracker = BestFilesTracker()

    # Files that should be classified as "rsr_output"
    rsr_patterns = [
        "rsr_001_real_space_refined_000.pdb",
        "rsr_001_real_space_refined_001.pdb",
        "rsr_002_real_space_refined_000.pdb",
        "model_real_space_refined.pdb",
    ]

    for filename in rsr_patterns:
        stage = tracker._classify_stage("/path/to/" + filename, "model")
        assert_equal(stage, "rsr_output", "Should classify as rsr_output: %s (got %s)" % (filename, stage))

    print("  PASSED")


def test_phaser_output_detection():
    """Test detection of phaser output files."""
    print("Test: phaser_output_detection")

    tracker = BestFilesTracker()

    # Files that should be classified as "phaser_output"
    phaser_patterns = [
        "PHASER.1.pdb",
        "PHASER.2.pdb",
        "phaser_output.pdb",
    ]

    for filename in phaser_patterns:
        stage = tracker._classify_stage("/path/to/" + filename, "model")
        assert_equal(stage, "phaser_output", "Should classify as phaser_output: %s (got %s)" % (filename, stage))

    print("  PASSED")


def test_predicted_model_detection():
    """Test detection of predicted model files."""
    print("Test: predicted_model_detection")

    tracker = BestFilesTracker()

    # Files that should be classified as "predicted"
    # Note: With semantic categories, predicted is a subcategory of search_model
    predicted_patterns = [
        "predicted_model.pdb",
        "alphafold_prediction.pdb",
    ]

    for filename in predicted_patterns:
        stage = tracker._classify_stage("/path/to/" + filename, "search_model")
        assert_equal(stage, "predicted",
                  "Should classify as predicted: %s (got %s)" % (filename, stage))

    print("  PASSED")


def test_autobuild_output_detection():
    """Test detection of autobuild output files."""
    print("Test: autobuild_output_detection")

    tracker = BestFilesTracker()

    # Files that should be classified as "autobuild_output"
    autobuild_patterns = [
        "overall_best.pdb",
        "autobuild_output.pdb",
    ]

    for filename in autobuild_patterns:
        stage = tracker._classify_stage("/path/to/" + filename, "model")
        assert_equal(stage, "autobuild_output",
                  "Should classify as autobuild_output: %s (got %s)" % (filename, stage))

    print("  PASSED")


# =============================================================================
# BEST FILE SELECTION TESTS
# =============================================================================

def test_best_model_selection_prefers_refined():
    """Test that best model selection prefers refined over input."""
    print("Test: best_model_selection_prefers_refined")

    tracker = BestFilesTracker()

    # Add files in order: input, then refined
    tracker.evaluate_file("/path/to/input_model.pdb", cycle=1, category="model")
    tracker.evaluate_file("/path/to/refine_001_001.pdb", cycle=2, category="model")

    best = tracker.get_best_path("model")
    assert_true(best is not None, "Should have a best model")
    assert_in("refine", best.lower(), "Should prefer refined model: %s" % best)

    print("  PASSED")


def test_best_model_selection_prefers_recent():
    """Test that best model selection prefers more recent refined models."""
    print("Test: best_model_selection_prefers_recent")

    tracker = BestFilesTracker()

    # Add multiple refined models
    tracker.evaluate_file("/path/to/refine_001_001.pdb", cycle=1, category="model")
    tracker.evaluate_file("/path/to/refine_002_001.pdb", cycle=2, category="model")
    tracker.evaluate_file("/path/to/refine_003_001.pdb", cycle=3, category="model")

    best = tracker.get_best_path("model")
    assert_true(best is not None, "Should have a best model")
    # Should prefer the most recent (refine_003)
    assert_in("003", best, "Should prefer most recent refined model: %s" % best)

    print("  PASSED")


def test_best_mtz_locked():
    """Test that best data_mtz respects R-free locking."""
    print("Test: best_mtz_locked")

    tracker = BestFilesTracker()

    # Add data_mtz files - data.mtz should lock on R-free
    tracker.evaluate_file("/path/to/data.mtz", cycle=1, category="data_mtz",
                          metrics={"has_rfree_flags": True})
    # This should NOT replace because data_mtz locks on R-free
    tracker.evaluate_file("/path/to/refine_001_data.mtz", cycle=2, category="data_mtz",
                          metrics={"has_rfree_flags": True})

    # Best data_mtz should be the locked one (first with R-free)
    best = tracker.get_best_path("data_mtz")
    assert_equal(best, "/path/to/data.mtz", "Should use locked R-free data_mtz")

    print("  PASSED")


# =============================================================================
# HALF-MAP DETECTION TESTS
# =============================================================================

def test_half_map_detection():
    """Test detection of half-map files."""
    print("Test: half_map_detection")

    # Half-map patterns that should be excluded from full map selection
    half_map_patterns = [
        "half_map_1.mrc",
        "half_map_2.mrc",
        "map_half1.ccp4",
        "map_half2.ccp4",
        "emd_12345_half_map_1.map",
        "emd_12345_half_map_2.map",
    ]

    for filename in half_map_patterns:
        basename = filename.lower()
        is_half = ("half" in basename or "_1" in basename or "_2" in basename)
        assert_true(is_half, "Should detect as half-map: %s" % filename)

    # Full map patterns that should NOT be excluded
    full_map_patterns = [
        "full_map.mrc",
        "density.ccp4",
        "denmod_map.ccp4",
    ]

    for filename in full_map_patterns:
        basename = filename.lower()
        is_half = ("half" in basename)
        assert_true(not is_half, "Should NOT detect as half-map: %s" % filename)

    print("  PASSED")


def test_single_half_map_promoted_to_full_map():
    """Test that a single half-map is promoted to full_map for usability."""
    print("Test: single_half_map_promoted_to_full_map")

    # Use _categorize_files which applies post-processing for half-map promotion
    from agent.workflow_state import _categorize_files

    # Case 1: Single file with half-map naming should be treated as full_map
    files_single_half = ["/path/to/emd-20026_auto_sharpen_A.ccp4"]
    result = _categorize_files(files_single_half)

    assert_in("/path/to/emd-20026_auto_sharpen_A.ccp4", result["full_map"],
              "Single map with _A suffix should be promoted to full_map")
    assert_equal(len(result["half_map"]), 0,
                 "half_map should be empty when single map is promoted")

    # Case 2: Two half-maps should remain as half-maps
    files_two_halves = [
        "/path/to/map_half_1.mrc",
        "/path/to/map_half_2.mrc",
    ]
    result = _categorize_files(files_two_halves)

    assert_equal(len(result["half_map"]), 2,
                 "Two half-maps should stay as half_maps")
    assert_equal(len(result["full_map"]), 0,
                 "full_map should be empty when we have two half-maps")

    # Case 3: One full map and one half-map - half-map stays as half-map
    files_mixed = [
        "/path/to/full_density.mrc",
        "/path/to/half_map_1.mrc",
    ]
    result = _categorize_files(files_mixed)

    assert_equal(len(result["full_map"]), 1,
                 "Full map should be detected")
    assert_equal(len(result["half_map"]), 1,
                 "Single half-map with full map present should stay as half_map")

    print("  PASSED")


# =============================================================================
# FILE NUMBER EXTRACTION TESTS
# =============================================================================

def test_file_number_extraction():
    """Test extraction of numbers from filenames for sorting."""
    print("Test: file_number_extraction")

    import re

    def get_file_number(path):
        basename = os.path.basename(path)
        numbers = re.findall(r'(\d+)', basename)
        if numbers:
            return int(numbers[-1])
        return 0

    # Test various filename patterns
    test_cases = [
        ("rsr_001_real_space_refined_000.pdb", 0),
        ("rsr_001_real_space_refined_001.pdb", 1),
        ("rsr_001_real_space_refined_002.pdb", 2),
        ("refine_001_001.pdb", 1),
        ("refine_002_001.pdb", 1),
        ("refine_002_003.pdb", 3),
        ("model.pdb", 0),
    ]

    for filename, expected in test_cases:
        result = get_file_number(filename)
        assert_equal(result, expected, "Number extraction for %s" % filename)

    print("  PASSED")


def test_file_sorting_by_number():
    """Test that files are sorted correctly by number."""
    print("Test: file_sorting_by_number")

    import re

    def get_file_number(path):
        basename = os.path.basename(path)
        numbers = re.findall(r'(\d+)', basename)
        if numbers:
            return int(numbers[-1])
        return 0

    files = [
        "/path/to/rsr_001_real_space_refined_000.pdb",
        "/path/to/rsr_001_real_space_refined_002.pdb",
        "/path/to/rsr_001_real_space_refined_001.pdb",
    ]

    # Sort by number descending (most recent first)
    sorted_files = sorted(files, key=get_file_number, reverse=True)

    assert_in("002", sorted_files[0], "First should be 002")
    assert_in("001", sorted_files[1], "Second should be 001")
    assert_in("000", sorted_files[2], "Third should be 000")

    print("  PASSED")


# =============================================================================
# CATEGORY EXCLUSION TESTS
# =============================================================================

def test_predicted_excluded_from_refinement():
    """Test that predicted models are excluded from refinement input."""
    print("Test: predicted_excluded_from_refinement")

    # The input_priorities for real_space_refine should exclude predicted
    from agent.program_registry import ProgramRegistry

    registry = ProgramRegistry()
    priorities = registry.get_input_priorities("phenix.real_space_refine", "model")

    exclude_cats = priorities.get("exclude_categories", [])
    assert_in("predicted", exclude_cats,
              "predicted should be excluded from RSR model input")
    assert_in("processed_predicted", exclude_cats,
              "processed_predicted should be excluded from RSR model input")

    print("  PASSED")


def test_intermediate_mr_excluded():
    """Test that intermediate MR files are excluded."""
    print("Test: intermediate_mr_excluded")

    from agent.program_registry import ProgramRegistry

    registry = ProgramRegistry()
    priorities = registry.get_input_priorities("phenix.real_space_refine", "model")

    exclude_cats = priorities.get("exclude_categories", [])
    assert_in("intermediate_mr", exclude_cats,
              "intermediate_mr should be excluded from RSR model input")

    print("  PASSED")


# =============================================================================
# INTEGRATION TESTS
# =============================================================================

def test_cryo_em_file_selection_workflow():
    """Test file selection for a typical cryo-EM workflow."""
    print("Test: cryo_em_file_selection_workflow")

    tracker = BestFilesTracker()

    # Test that stage classification is correct for cryo-EM files
    # Note: predicted_model.pdb should use search_model category
    model_stages = {
        "/data/PHASER.1.pdb": "phaser_output",
        "/data/rsr_001_real_space_refined_000.pdb": "rsr_output",
        "/data/rsr_001_real_space_refined_001.pdb": "rsr_output",
    }

    for path, expected_stage in model_stages.items():
        stage = tracker._classify_stage(path, "model")
        assert_equal(stage, expected_stage, "Stage for %s" % path)

    # Predicted model uses search_model category
    stage = tracker._classify_stage("/data/predicted_model.pdb", "search_model")
    assert_equal(stage, "predicted", "Stage for predicted_model.pdb")

    # Verify RSR outputs are tracked correctly when evaluated
    tracker.evaluate_file("/data/rsr_001_real_space_refined_000.pdb", cycle=4, category="model")
    tracker.evaluate_file("/data/rsr_001_real_space_refined_001.pdb", cycle=5, category="model")

    best = tracker.get_best_path("model")
    assert_true(best is not None, "Should have a best model")
    # Most recent RSR output should win
    assert_in("001", best, "Should prefer most recent RSR output: %s" % best)

    print("  PASSED")


def test_xray_file_selection_workflow():
    """Test file selection for a typical X-ray workflow."""
    print("Test: xray_file_selection_workflow")

    tracker = BestFilesTracker()

    # Test that stage classification is correct for X-ray files
    # Note: With semantic categories, generic model.pdb gets 'model' as default stage
    model_stages = {
        "/data/refine_001_001.pdb": "refined",
        "/data/refine_002_001.pdb": "refined",
    }

    for path, expected_stage in model_stages.items():
        stage = tracker._classify_stage(path, "model")
        assert_equal(stage, expected_stage, "Stage for %s" % path)

    # Generic model.pdb gets default 'model' stage
    stage = tracker._classify_stage("/data/model.pdb", "model")
    assert_equal(stage, "model", "Stage for generic model.pdb")

    # Test data_mtz stages
    data_mtz_stages = {
        "/data/data.mtz": "data_mtz",           # Generic goes to data_mtz
        "/data/refine_001_data.mtz": "original_data_mtz",  # Data output
    }

    for path, expected_stage in data_mtz_stages.items():
        stage = tracker._classify_stage(path, "data_mtz")
        assert_equal(stage, expected_stage, "Stage for %s" % path)

    # Test map_coeffs_mtz stages
    coeffs_stages = {
        "/data/refine_001_001.mtz": "refine_map_coeffs",
        "/data/denmod_map_coeffs.mtz": "denmod_map_coeffs",
    }

    for path, expected_stage in coeffs_stages.items():
        stage = tracker._classify_stage(path, "map_coeffs_mtz")
        assert_equal(stage, expected_stage, "Stage for %s" % path)

    # Verify refined models are tracked correctly
    tracker.evaluate_file("/data/refine_001_001.pdb", cycle=2, category="model")
    tracker.evaluate_file("/data/refine_002_001.pdb", cycle=3, category="model")

    best_model = tracker.get_best_path("model")
    assert_in("refine_002", best_model, "Should use most recent refined model: %s" % best_model)

    # Test data_mtz R-free locking
    tracker.evaluate_file("/data/data.mtz", cycle=1, category="data_mtz",
                          metrics={"has_rfree_flags": True})
    tracker.evaluate_file("/data/refine_001_data.mtz", cycle=2, category="data_mtz",
                          metrics={"has_rfree_flags": True})

    best_mtz = tracker.get_best_path("data_mtz")
    assert_equal(best_mtz, "/data/data.mtz", "Should use locked R-free data_mtz")

    print("  PASSED")


def test_cryo_em_map_categorization():
    """Test categorization of cryo-EM maps including resolve_cryo_em outputs (v74)."""
    print("Test: cryo_em_map_categorization")

    from agent.best_files_tracker import BestFilesTracker

    tracker = BestFilesTracker()

    # Test map classification via _classify_stage
    map_stages = {
        # Optimized maps from resolve_cryo_em
        "/path/rcm_0/denmod_map.ccp4": "optimized_full_map",
        "/path/density_modified_map.mrc": "optimized_full_map",
        "/path/sharpened_map.ccp4": "sharpened",

        # Intermediate maps (should be deprioritized)
        "/path/rcm_0/initial_map.ccp4": "intermediate_map",
        "/path/initial.mrc": "intermediate_map",

        # Regular full maps
        "/path/reconstruction.mrc": "full_map",
        "/path/full_map.ccp4": "full_map",

        # Half maps (use realistic naming conventions)
        "/path/map_half1.mrc": "half_map",
        "/path/emd_1234_half_map_1.mrc": "half_map",
    }

    for path, expected_stage in map_stages.items():
        stage = tracker._classify_stage(path, "map")
        assert_equal(stage, expected_stage, "Stage for %s" % os.path.basename(path))

    print("  PASSED")


def test_optimized_full_map_preferred():
    """Test that optimized_full_map is preferred over initial_map."""
    print("Test: optimized_full_map_preferred")

    from agent.best_files_tracker import BestFilesTracker

    tracker = BestFilesTracker()

    # Add intermediate map first
    tracker.evaluate_file("/path/rcm_0/initial_map.ccp4", cycle=1, stage="intermediate_map")

    # Add optimized map
    tracker.evaluate_file("/path/rcm_0/denmod_map.ccp4", cycle=1, stage="optimized_full_map")

    best = tracker.get_best("map")
    assert_not_none(best, "Should have best map")
    assert_in("denmod", best.path, "Should prefer denmod_map over initial_map")

    print("  PASSED")


def test_docked_model_bubbles_to_model():
    """Test that docked models bubble up to the model category."""
    print("Test: docked_model_bubbles_to_model")

    from agent.workflow_state import _categorize_files

    # Test file list with docked model
    files = [
        "/path/sequence.fa",
        "/path/map.mrc",
        "/path/placed_model.pdb",
        "/path/predicted.pdb",
    ]

    categorized = _categorize_files(files)

    # placed_model should be in docked category
    assert_in("docked", categorized, "Should have docked category")
    assert_true(any("placed_model" in f for f in categorized.get("docked", [])),
               "placed_model.pdb should be in docked category")

    # placed_model should also bubble up to model category
    assert_in("model", categorized, "Should have model category")
    assert_true(any("placed_model" in f for f in categorized.get("model", [])),
               "placed_model.pdb should bubble up to model category")

    print("  PASSED")


def test_predict_and_build_output_is_model():
    """
    Test that PredictAndBuild output files are categorized as model, not search_model.

    Bug fix: PredictAndBuild_0_overall_best.pdb was incorrectly being categorized
    as search_model (because it contains "predict") instead of model.

    These are fully built/refined models, not search models for MR.
    """
    print("test_predict_and_build_output_is_model:", end=" ")

    from agent.workflow_state import _categorize_files

    files = [
        "/path/PredictAndBuild_0_overall_best.pdb",
        "/path/PredictAndBuild_0_overall_best_superposed_predicted_models.pdb",
        "/path/PredictAndBuild_0_superposed_predicted_untrimmed_cycle_1.pdb",
        "/path/7qz0.fa",
        "/path/7qz0.mtz",
    ]

    categorized = _categorize_files(files)

    # PredictAndBuild_0_overall_best.pdb should be in model category (not just search_model)
    assert_in("model", categorized, "Should have model category")
    model_files = categorized.get("model", [])
    assert_true(any("overall_best.pdb" in f for f in model_files),
               "PredictAndBuild_0_overall_best.pdb should be in model category, got: %s" % model_files)

    # It should also be in predict_and_build_output subcategory
    assert_in("predict_and_build_output", categorized, "Should have predict_and_build_output category")
    pab_files = categorized.get("predict_and_build_output", [])
    assert_true(any("overall_best.pdb" in f for f in pab_files),
               "PredictAndBuild_0_overall_best.pdb should be in predict_and_build_output category")

    # The intermediate predicted models should still be in search_model
    search_files = categorized.get("search_model", [])
    assert_true(any("superposed_predicted_untrimmed" in f for f in search_files),
               "Intermediate predicted models should still be in search_model")

    # Critically: overall_best should NOT be in search_model
    assert_true(not any("overall_best.pdb" in f for f in search_files),
               "overall_best.pdb should NOT be in search_model, but found in: %s" % search_files)

    print("  PASSED")


def test_ligand_file_patterns():
    """
    Test that various ligand file naming patterns are correctly detected.

    Bug fix: 7qz0_ligand.pdb wasn't being detected as a ligand file because
    the pattern was ligand_*.pdb (ligand at start) but not *_ligand.pdb.
    """
    print("test_ligand_file_patterns:", end=" ")

    from agent.workflow_state import _categorize_files

    files = [
        "/path/lig.pdb",           # Standard short name
        "/path/ligand.pdb",        # Standard name
        "/path/ligand_001.pdb",    # Numbered ligand
        "/path/7qz0_ligand.pdb",   # Protein code + ligand suffix
        "/path/abc_lig.pdb",       # Protein code + lig suffix
        "/path/model.pdb",         # NOT a ligand
    ]

    categorized = _categorize_files(files)

    # Check ligand category exists
    assert_in("ligand", categorized, "Should have ligand category")
    ligand_files = categorized.get("ligand", [])

    # All ligand patterns should be detected
    assert_true(any("lig.pdb" in f for f in ligand_files),
               "lig.pdb should be in ligand category")
    assert_true(any("ligand.pdb" in f for f in ligand_files),
               "ligand.pdb should be in ligand category")
    assert_true(any("7qz0_ligand.pdb" in f for f in ligand_files),
               "7qz0_ligand.pdb should be in ligand category, got: %s" % ligand_files)
    assert_true(any("abc_lig.pdb" in f for f in ligand_files),
               "abc_lig.pdb should be in ligand category")

    # model.pdb should NOT be in ligand
    assert_true(not any("model.pdb" in f for f in ligand_files),
               "model.pdb should NOT be in ligand category")

    print("  PASSED")


def test_predict_and_build_mtz_detection():
    """
    Test that predict_and_build map coefficient files are correctly categorized.

    Bug fix: PredictAndBuild_0_overall_best_map_coeffs.mtz was not being
    recognized for ligandfit. The _refinement.mtz doesn't have map coefficients.
    """
    print("test_predict_and_build_mtz_detection:", end=" ")

    from agent.workflow_state import _categorize_files

    files = [
        "/path/PredictAndBuild_0_overall_best_map_coeffs.mtz",  # Has FP/PHIFP
        "/path/PredictAndBuild_0_overall_best_refinement.mtz", # Has data, NOT map coeffs
        "/path/refine_001_001.mtz",  # Standard refined MTZ (map coefficients)
    ]

    categorized = _categorize_files(files)

    # Check predict_build_map_coeffs category exists (renamed from predict_and_build_mtz)
    assert_in("predict_build_map_coeffs", categorized,
             "Should have predict_build_map_coeffs category, got: %s" % list(categorized.keys()))

    pab_mtz = categorized.get("predict_build_map_coeffs", [])
    assert_true(any("map_coeffs" in f for f in pab_mtz),
               "map_coeffs.mtz should be in predict_build_map_coeffs category")

    # The _refinement.mtz should NOT be in predict_build_map_coeffs (it's data, not map coeffs)
    assert_true(not any("refinement.mtz" in f for f in pab_mtz),
               "_refinement.mtz should NOT be in predict_build_map_coeffs")

    # The _refinement.mtz should be in data_mtz (via original_data_mtz)
    assert_in("data_mtz", categorized)
    data_mtz = categorized.get("data_mtz", [])
    assert_true(any("refinement.mtz" in f for f in data_mtz),
               "_refinement.mtz should be in data_mtz")

    print("  PASSED")


# =============================================================================
# RUN ALL TESTS
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.tst_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
