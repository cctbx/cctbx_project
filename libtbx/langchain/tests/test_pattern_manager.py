#!/usr/bin/env python
"""
Unit tests for pattern_manager.py

Tests:
- Pattern loading and interpolation
- Primitive replacement
- Match/search/findall operations
- Metric extraction
- Pattern validation
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agent.pattern_manager import get_patterns


# =============================================================================
# TEST FIXTURES
# =============================================================================

def get_test_manager():
    """Get a PatternManager instance for testing."""
    return get_patterns()


# =============================================================================
# PRIMITIVE TESTS
# =============================================================================

def test_primitives_loaded():
    """Test that primitives are loaded from YAML."""
    print("Test: primitives_loaded")

    pm = get_test_manager()
    primitives = pm.list_primitives()

    assert "FLOAT" in primitives, "Should have FLOAT primitive"
    assert "INT" in primitives, "Should have INT primitive"
    assert "OPT_WS" in primitives, "Should have OPT_WS primitive"

    print("  PASSED")


def test_primitive_interpolation():
    """Test that primitives are correctly interpolated."""
    print("Test: primitive_interpolation")

    pm = get_test_manager()

    # Test single primitive
    result = pm._interpolate("{INT}")
    assert "\\d" in result, "INT should contain digit pattern"

    # Test multiple primitives
    result = pm._interpolate("value{OPT_WS}={OPT_WS}{FLOAT}")
    assert "\\s*" in result, "Should have optional whitespace"
    assert "{FLOAT}" not in result, "FLOAT should be replaced"

    print("  PASSED")


# =============================================================================
# PATTERN LOADING TESTS
# =============================================================================

def test_patterns_loaded():
    """Test that patterns are loaded from YAML."""
    print("Test: patterns_loaded")

    pm = get_test_manager()

    # Check sections exist
    sections = pm.list_sections()
    assert "metrics" in sections, "Should have metrics section"
    assert "system" in sections, "Should have system section"

    # Check some patterns exist
    assert pm.has_pattern("metrics.r_free"), "Should have r_free pattern"
    assert pm.has_pattern("system.refine_output"), "Should have refine_output pattern"

    print("  PASSED")


def test_pattern_description():
    """Test getting pattern descriptions."""
    print("Test: pattern_description")

    pm = get_test_manager()

    desc = pm.get_description("metrics.r_free")
    assert "R-free" in desc or "refinement" in desc.lower(), \
        "r_free description should mention R-free or refinement"

    print("  PASSED")


# =============================================================================
# MATCHING TESTS
# =============================================================================

def test_match_refine_output():
    """Test matching refinement output filenames."""
    print("Test: match_refine_output")

    pm = get_test_manager()

    # Should match
    match = pm.match("system.refine_output", "refine_001_002.pdb")
    assert match is not None, "Should match refine_001_002.pdb"
    assert match.group(1) == "001", "First group should be 001"
    assert match.group(2) == "002", "Second group should be 002"

    # Should not match
    match = pm.match("system.refine_output", "model.pdb")
    assert match is None, "Should not match model.pdb"

    print("  PASSED")


def test_search_r_free():
    """Test searching for R-free in log text."""
    print("Test: search_r_free")

    pm = get_test_manager()

    # Various formats
    tests = [
        ("R-free = 0.2534", "0.2534"),
        ("R-free: 0.25", "0.25"),
        ("Rfree=0.30", "0.30"),
        ("Final R-free = 0.28", "0.28"),
    ]

    for text, expected in tests:
        match = pm.search("metrics.r_free", text)
        assert match is not None, "Should match: %s" % text
        assert match.group(1) == expected, \
            "Expected %s, got %s for: %s" % (expected, match.group(1), text)

    print("  PASSED")


def test_search_case_insensitive():
    """Test case-insensitive searching."""
    print("Test: search_case_insensitive")

    pm = get_test_manager()
    import re

    # Should match with IGNORECASE
    match = pm.search("metrics.r_free", "r-free = 0.25", flags=re.IGNORECASE)
    assert match is not None, "Should match lowercase r-free with IGNORECASE"

    print("  PASSED")


def test_findall():
    """Test finding all matches."""
    print("Test: findall")

    pm = get_test_manager()

    text = "R-free = 0.30, then R-free = 0.28, finally R-free = 0.25"
    matches = pm.findall("metrics.r_free", text)

    assert len(matches) == 3, "Should find 3 R-free values"
    assert "0.30" in matches, "Should find 0.30"
    assert "0.28" in matches, "Should find 0.28"
    assert "0.25" in matches, "Should find 0.25"

    print("  PASSED")


# =============================================================================
# EXTRACTION TESTS
# =============================================================================

def test_extract_string():
    """Test extracting string values."""
    print("Test: extract_string")

    pm = get_test_manager()

    value = pm.extract("metrics.r_free", "R-free = 0.2534")
    assert value == "0.2534", "Should extract 0.2534"

    # With default
    value = pm.extract("metrics.r_free", "no r-free here", default="N/A")
    assert value == "N/A", "Should return default when no match"

    print("  PASSED")


def test_extract_float():
    """Test extracting float values."""
    print("Test: extract_float")

    pm = get_test_manager()

    value = pm.extract_float("metrics.r_free", "R-free = 0.2534")
    assert abs(value - 0.2534) < 0.0001, "Should extract 0.2534 as float"

    # With default
    value = pm.extract_float("metrics.r_free", "no match", default=1.0)
    assert value == 1.0, "Should return default"

    print("  PASSED")


def test_extract_int():
    """Test extracting integer values."""
    print("Test: extract_int")

    pm = get_test_manager()

    value = pm.extract_int("system.refine_output", "refine_003_001.pdb", group=1)
    assert value == 3, "Should extract 3 as int"

    print("  PASSED")


# =============================================================================
# METRIC EXTRACTION TESTS
# =============================================================================

def test_extract_metrics_from_log():
    """Test extracting multiple metrics from a log."""
    print("Test: extract_metrics_from_log")

    pm = get_test_manager()

    log = """
    Refinement cycle 5 completed.
    R-work = 0.2100
    R-free = 0.2534
    Resolution: 2.5 A
    Clashscore = 5.2
    Ramachandran outliers = 0.5 %
    """

    r_free = pm.extract_float("metrics.r_free", log)
    r_work = pm.extract_float("metrics.r_work", log)
    resolution = pm.extract_float("metrics.resolution", log)
    clashscore = pm.extract_float("metrics.clashscore", log)

    assert abs(r_free - 0.2534) < 0.0001, "R-free should be 0.2534"
    assert abs(r_work - 0.2100) < 0.0001, "R-work should be 0.2100"
    assert abs(resolution - 2.5) < 0.0001, "Resolution should be 2.5"
    assert abs(clashscore - 5.2) < 0.0001, "Clashscore should be 5.2"

    print("  PASSED")


def test_extract_phaser_metrics():
    """Test extracting Phaser metrics."""
    print("Test: extract_phaser_metrics")

    pm = get_test_manager()

    log = """
    SOLU SET  RFZ=8.2 TFZ==24.5 PAK=0 LLG=1500
    """

    tfz = pm.extract_float("metrics.tfz_score", log)
    llg = pm.extract_float("metrics.llg", log)

    assert tfz is not None, "Should find TFZ"
    assert abs(tfz - 24.5) < 0.1, "TFZ should be 24.5"
    assert llg is not None, "Should find LLG"
    assert abs(llg - 1500) < 1, "LLG should be 1500"

    print("  PASSED")


# =============================================================================
# VALIDATION TESTS
# =============================================================================

def test_validate_pattern():
    """Test pattern self-validation."""
    print("Test: validate_pattern")

    pm = get_test_manager()

    passed, failed = pm.validate_pattern("metrics.r_free")
    assert len(failed) == 0, "r_free validation should pass: %s" % failed
    assert len(passed) > 0, "r_free should have some tests"

    print("  PASSED")


def test_validate_all():
    """Test validating all patterns."""
    print("Test: validate_all")

    pm = get_test_manager()

    results = pm.validate_all()

    all_passed = True
    for name, (passed, failed) in results.items():
        if failed:
            print("  FAILED: %s - %s" % (name, failed))
            all_passed = False

    assert all_passed, "All pattern validations should pass"

    print("  PASSED")


# =============================================================================
# EDGE CASE TESTS
# =============================================================================

def test_unknown_pattern():
    """Test handling of unknown pattern names."""
    print("Test: unknown_pattern")

    pm = get_test_manager()

    try:
        pm.get_compiled("nonexistent.pattern")
        assert False, "Should raise KeyError"
    except KeyError:
        pass  # Expected

    print("  PASSED")


def test_empty_text():
    """Test handling of empty text."""
    print("Test: empty_text")

    pm = get_test_manager()

    match = pm.search("metrics.r_free", "")
    assert match is None, "Should return None for empty text"

    value = pm.extract_float("metrics.r_free", "", default=0.0)
    assert value == 0.0, "Should return default for empty text"

    print("  PASSED")


# =============================================================================
# RUN ALL TESTS
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.test_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
