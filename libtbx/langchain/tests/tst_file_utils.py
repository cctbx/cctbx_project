"""
Unit tests for the file_utils module.

Tests the shared file classification utilities used across the codebase.

Run with: python tests/test_file_utils.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_false,
    run_tests_with_fail_fast
)

from agent.file_utils import (
    classify_mtz_type,
    get_mtz_stage,
    is_mtz_file,
    is_model_file,
    is_map_file,
    is_sequence_file,
)


# =============================================================================
# MTZ CLASSIFICATION TESTS
# =============================================================================

def test_classify_mtz_refine_map_coeffs():
    """Refine output with maps is classified as map_coeffs_mtz."""
    print("Test: classify_mtz_refine_map_coeffs")

    # Standard refine map output
    assert_equal(classify_mtz_type("/path/to/refine_001_001.mtz"), "map_coeffs_mtz")
    assert_equal(classify_mtz_type("/path/to/refine_002_001.mtz"), "map_coeffs_mtz")
    assert_equal(classify_mtz_type("/path/to/refine_123_001.mtz"), "map_coeffs_mtz")

    print("  PASSED")


def test_classify_mtz_numbered_output():
    """Other _001.mtz outputs are classified as map_coeffs_mtz."""
    print("Test: classify_mtz_numbered_output")

    assert_equal(classify_mtz_type("/path/to/model_refine_001.mtz"), "map_coeffs_mtz")
    assert_equal(classify_mtz_type("/path/to/output_001.mtz"), "map_coeffs_mtz")

    print("  PASSED")


def test_classify_mtz_explicit_map_coeffs():
    """Files with 'map_coeffs' in name are classified as map_coeffs_mtz."""
    print("Test: classify_mtz_explicit_map_coeffs")

    assert_equal(classify_mtz_type("/path/to/map_coeffs.mtz"), "map_coeffs_mtz")
    assert_equal(classify_mtz_type("/path/to/overall_best_map_coeffs.mtz"), "map_coeffs_mtz")
    assert_equal(classify_mtz_type("/path/to/predict_map_coeffs.mtz"), "map_coeffs_mtz")

    print("  PASSED")


def test_classify_mtz_denmod():
    """Density-modified files are classified as map_coeffs_mtz."""
    print("Test: classify_mtz_denmod")

    assert_equal(classify_mtz_type("/path/to/denmod.mtz"), "map_coeffs_mtz")
    assert_equal(classify_mtz_type("/path/to/model_denmod.mtz"), "map_coeffs_mtz")
    assert_equal(classify_mtz_type("/path/to/density_mod.mtz"), "map_coeffs_mtz")

    print("  PASSED")


def test_classify_mtz_data_files():
    """Data MTZ files are classified correctly."""
    print("Test: classify_mtz_data_files")

    # Explicit data files
    assert_equal(classify_mtz_type("/path/to/data.mtz"), "data_mtz")
    assert_equal(classify_mtz_type("/path/to/refine_001_data.mtz"), "data_mtz")
    assert_equal(classify_mtz_type("/path/to/refinement_data.mtz"), "data_mtz")

    # Generic MTZ files default to data
    assert_equal(classify_mtz_type("/path/to/input.mtz"), "data_mtz")
    assert_equal(classify_mtz_type("/path/to/p9.mtz"), "data_mtz")

    print("  PASSED")


def test_classify_mtz_case_insensitive():
    """MTZ classification is case-insensitive."""
    print("Test: classify_mtz_case_insensitive")

    assert_equal(classify_mtz_type("/path/to/REFINE_001_001.MTZ"), "map_coeffs_mtz")
    assert_equal(classify_mtz_type("/path/to/Map_Coeffs.mtz"), "map_coeffs_mtz")
    assert_equal(classify_mtz_type("/path/to/DENMOD.MTZ"), "map_coeffs_mtz")

    print("  PASSED")


# =============================================================================
# MTZ STAGE TESTS
# =============================================================================

def test_get_mtz_stage_data():
    """Data MTZ stages are identified correctly."""
    print("Test: get_mtz_stage_data")

    assert_equal(get_mtz_stage("/path/to/refine_001_data.mtz", "data_mtz"), "original_data_mtz")
    assert_equal(get_mtz_stage("/path/to/phased.mtz", "data_mtz"), "phased_data_mtz")
    assert_equal(get_mtz_stage("/path/to/input.mtz", "data_mtz"), "data_mtz")

    print("  PASSED")


def test_get_mtz_stage_map_coeffs():
    """Map coeffs MTZ stages are identified correctly."""
    print("Test: get_mtz_stage_map_coeffs")

    assert_equal(get_mtz_stage("/path/to/denmod.mtz", "map_coeffs_mtz"), "denmod_map_coeffs")
    assert_equal(get_mtz_stage("/path/to/predict_map_coeffs.mtz", "map_coeffs_mtz"), "predict_build_map_coeffs")
    assert_equal(get_mtz_stage("/path/to/overall_best_map_coeffs.mtz", "map_coeffs_mtz"), "predict_build_map_coeffs")
    assert_equal(get_mtz_stage("/path/to/refine_001_001.mtz", "map_coeffs_mtz"), "refine_map_coeffs")

    print("  PASSED")


# =============================================================================
# FILE TYPE DETECTION TESTS
# =============================================================================

def test_is_mtz_file():
    """MTZ file detection works."""
    print("Test: is_mtz_file")

    assert_true(is_mtz_file("/path/to/data.mtz"))
    assert_true(is_mtz_file("/path/to/DATA.MTZ"))
    assert_false(is_mtz_file("/path/to/model.pdb"))
    assert_false(is_mtz_file("/path/to/map.mrc"))

    print("  PASSED")


def test_is_model_file():
    """Model file detection works."""
    print("Test: is_model_file")

    assert_true(is_model_file("/path/to/model.pdb"))
    assert_true(is_model_file("/path/to/structure.ent"))
    assert_true(is_model_file("/path/to/model.cif"))
    assert_false(is_model_file("/path/to/ligand.cif"))  # Has 'ligand' in name
    assert_false(is_model_file("/path/to/data.mtz"))

    print("  PASSED")


def test_is_map_file():
    """Map file detection works."""
    print("Test: is_map_file")

    assert_true(is_map_file("/path/to/density.map"))
    assert_true(is_map_file("/path/to/emd_1234.mrc"))
    assert_true(is_map_file("/path/to/map.ccp4"))
    assert_false(is_map_file("/path/to/model.pdb"))

    print("  PASSED")


def test_is_sequence_file():
    """Sequence file detection works."""
    print("Test: is_sequence_file")

    assert_true(is_sequence_file("/path/to/sequence.fa"))
    assert_true(is_sequence_file("/path/to/seq.fasta"))
    assert_true(is_sequence_file("/path/to/protein.seq"))
    assert_true(is_sequence_file("/path/to/seq.dat"))
    assert_false(is_sequence_file("/path/to/model.pdb"))

    print("  PASSED")


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
