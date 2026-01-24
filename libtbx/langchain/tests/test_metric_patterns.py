#!/usr/bin/env python
"""
Tests for knowledge/metric_patterns.py
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_pattern_loading():
    """Test that patterns load from YAML."""
    from knowledge.metric_patterns import get_all_metric_patterns

    patterns = get_all_metric_patterns()

    # Should have patterns for multiple programs
    assert len(patterns) > 0, "Should load patterns for at least one program"

    # Check specific programs exist
    assert "phenix.map_symmetry" in patterns, "Should have map_symmetry patterns"
    assert "phenix.refine" in patterns, "Should have refine patterns"
    assert "phenix.mtriage" in patterns, "Should have mtriage patterns"

    print("  test_pattern_loading PASSED")


def test_map_symmetry_extraction():
    """Test map_symmetry metric extraction."""
    from knowledge.metric_patterns import extract_metrics_for_program

    log_text = """
    Best NCS type is:
      SCORE    CC   OPERATORS     SYMMETRY
       3.40   0.93    14          D7 (a)  Best NCS type

    Final symmetry obtained:
    NCS type: D7 (a)
    Correlation of symmetry-related regions: 0.93   Copies: 14
    """

    metrics = extract_metrics_for_program(log_text, "phenix.map_symmetry")

    assert "symmetry_type" in metrics, "Should extract symmetry_type"
    assert metrics["symmetry_type"] == "D7 (a)", f"Expected 'D7 (a)', got {metrics['symmetry_type']}"
    assert "ncs_copies" in metrics, "Should extract ncs_copies"
    assert metrics["ncs_copies"] == 14, f"Expected 14, got {metrics['ncs_copies']}"
    assert "ncs_cc" in metrics, "Should extract ncs_cc"
    assert abs(metrics["ncs_cc"] - 0.93) < 0.01, f"Expected 0.93, got {metrics['ncs_cc']}"

    print("  test_map_symmetry_extraction PASSED")


def test_no_symmetry_extraction():
    """Test that 'no symmetry' case is handled."""
    from knowledge.metric_patterns import extract_metrics_for_program

    log_text = "No suitable symmetry found"

    metrics = extract_metrics_for_program(log_text, "phenix.map_symmetry")

    assert "symmetry_type" in metrics, "Should extract symmetry_type for no-match case"
    assert metrics["symmetry_type"] == "None", f"Expected 'None', got {metrics['symmetry_type']}"

    print("  test_no_symmetry_extraction PASSED")


def test_refine_extraction():
    """Test phenix.refine metric extraction."""
    from knowledge.metric_patterns import extract_metrics_for_program

    log_text = """
    phenix.refine

    Final R-work = 0.2150
    Final R-free = 0.2580
    Bond RMSD: 0.0123
    Angle RMSD: 1.456
    """

    metrics = extract_metrics_for_program(log_text, "phenix.refine")

    assert "r_free" in metrics, "Should extract r_free"
    assert abs(metrics["r_free"] - 0.258) < 0.001, f"Expected 0.258, got {metrics['r_free']}"
    assert "r_work" in metrics, "Should extract r_work"
    assert abs(metrics["r_work"] - 0.215) < 0.001, f"Expected 0.215, got {metrics['r_work']}"

    print("  test_refine_extraction PASSED")


def test_mtriage_extraction():
    """Test phenix.mtriage metric extraction."""
    from knowledge.metric_patterns import extract_metrics_for_program

    log_text = """
    phenix.mtriage

    d_fsc_model_05 = 3.21
    CC_mask = 0.856
    """

    metrics = extract_metrics_for_program(log_text, "phenix.mtriage")

    assert "resolution" in metrics, "Should extract resolution"
    assert abs(metrics["resolution"] - 3.21) < 0.01, f"Expected 3.21, got {metrics['resolution']}"
    assert "map_cc" in metrics, "Should extract map_cc"
    assert abs(metrics["map_cc"] - 0.856) < 0.001, f"Expected 0.856, got {metrics['map_cc']}"

    print("  test_mtriage_extraction PASSED")


def test_molprobity_extraction():
    """Test phenix.molprobity metric extraction."""
    from knowledge.metric_patterns import extract_metrics_for_program

    log_text = """
    phenix.molprobity

    Clashscore = 4.56
    Ramachandran favored = 96.5
    Ramachandran outliers = 0.3
    Rotamer outliers = 1.2
    MolProbity score = 1.85
    """

    metrics = extract_metrics_for_program(log_text, "phenix.molprobity")

    assert "clashscore" in metrics, "Should extract clashscore"
    assert abs(metrics["clashscore"] - 4.56) < 0.01, f"Expected 4.56, got {metrics['clashscore']}"
    assert "ramachandran_favored" in metrics, "Should extract ramachandran_favored"
    assert "molprobity_score" in metrics, "Should extract molprobity_score"

    print("  test_molprobity_extraction PASSED")


def test_phaser_extraction():
    """Test phenix.phaser metric extraction."""
    from knowledge.metric_patterns import extract_metrics_for_program

    log_text = """
    phenix.phaser

    TFZ=12.5
    LLG=1234
    RFZ=8.2
    SOLU SPAC P 21 21 21
    """

    metrics = extract_metrics_for_program(log_text, "phenix.phaser")

    assert "tfz" in metrics, "Should extract tfz"
    assert abs(metrics["tfz"] - 12.5) < 0.1, f"Expected 12.5, got {metrics['tfz']}"
    assert "llg" in metrics, "Should extract llg"
    assert "space_group" in metrics, "Should extract space_group"

    print("  test_phaser_extraction PASSED")


def test_format_metric_value():
    """Test metric value formatting."""
    from knowledge.metric_patterns import format_metric_value

    # Test R-free formatting
    formatted = format_metric_value("phenix.refine", "r_free", 0.258)
    assert "0.2580" in formatted, f"Expected '0.2580' in '{formatted}'"

    # Test symmetry formatting
    formatted = format_metric_value("phenix.map_symmetry", "symmetry_type", "D7 (a)")
    assert "D7 (a)" in formatted, f"Expected 'D7 (a)' in '{formatted}'"

    print("  test_format_metric_value PASSED")


def test_display_config():
    """Test getting display configuration."""
    from knowledge.metric_patterns import get_metric_display_config

    config = get_metric_display_config("phenix.refine", "r_free")

    assert "display_name" in config, "Should have display_name"
    assert config["display_name"] == "R-free", f"Expected 'R-free', got {config['display_name']}"

    print("  test_display_config PASSED")


def test_empty_log():
    """Test handling of empty log."""
    from knowledge.metric_patterns import extract_metrics_for_program

    metrics = extract_metrics_for_program("", "phenix.refine")
    assert metrics == {}, "Empty log should return empty dict"

    metrics = extract_metrics_for_program(None, "phenix.refine")
    assert metrics == {}, "None log should return empty dict"

    print("  test_empty_log PASSED")


def test_unknown_program():
    """Test handling of unknown program."""
    from knowledge.metric_patterns import extract_metrics_for_program

    metrics = extract_metrics_for_program("some log text", "phenix.unknown_program")
    assert metrics == {}, "Unknown program should return empty dict"

    print("  test_unknown_program PASSED")


def run_all_tests():
    """Run all tests."""
    print("\nRunning metric_patterns tests...")

    test_pattern_loading()
    test_map_symmetry_extraction()
    test_no_symmetry_extraction()
    test_refine_extraction()
    test_mtriage_extraction()
    test_molprobity_extraction()
    test_phaser_extraction()
    test_format_metric_value()
    test_display_config()
    test_empty_log()
    test_unknown_program()

    print("\nAll metric_patterns tests PASSED!")


if __name__ == "__main__":
    run_all_tests()
