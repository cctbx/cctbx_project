#!/usr/bin/env python
"""
Tests for knowledge/summary_display.py
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_quality_table_config_loading():
    """Test that quality table config loads from YAML."""
    from knowledge.summary_display import get_quality_table_config

    config = get_quality_table_config()

    # Should have multiple rows configured
    assert len(config) > 0, "Should load quality table config"

    # Check structure
    for row in config:
        assert 'metrics' in row, "Each row should have metrics"
        assert 'label' in row, "Each row should have label"
        assert 'format' in row, "Each row should have format"

    print("  test_quality_table_config_loading PASSED")


def test_step_metrics_config_loading():
    """Test that step metrics config loads from YAML."""
    from knowledge.summary_display import get_step_metrics_config

    config = get_step_metrics_config()

    # Should have configs for multiple programs
    assert len(config) > 0, "Should load step metrics config"

    # Check specific programs exist
    assert 'phenix.refine' in config, "Should have phenix.refine config"
    assert 'phenix.real_space_refine' in config, "Should have phenix.real_space_refine config"
    assert 'phenix.map_symmetry' in config, "Should have phenix.map_symmetry config"

    print("  test_step_metrics_config_loading PASSED")


def test_format_quality_table_xray():
    """Test formatting quality table for X-ray data."""
    from knowledge.summary_display import format_quality_table_rows

    metrics = {
        'r_free': 0.2345,
        'r_work': 0.2012,
        'clashscore': 3.5,
        'resolution': 2.1,
    }

    rows = format_quality_table_rows(metrics, 'xray')

    # Should have rows
    assert len(rows) > 0, "Should format some rows"

    # Check R-free row
    r_free_rows = [r for r in rows if r['label'] == 'R-free']
    assert len(r_free_rows) == 1, "Should have R-free row"
    assert '0.2345' in r_free_rows[0]['value'], f"R-free value should contain 0.2345, got {r_free_rows[0]['value']}"

    # Check clashscore row
    clash_rows = [r for r in rows if r['label'] == 'Clashscore']
    assert len(clash_rows) == 1, "Should have Clashscore row"

    print("  test_format_quality_table_xray PASSED")


def test_format_quality_table_cryoem():
    """Test formatting quality table for cryo-EM data."""
    from knowledge.summary_display import format_quality_table_rows

    metrics = {
        'map_cc': 0.856,
        'symmetry_type': 'D7 (a)',
        'ncs_copies': 14,
        'ncs_cc': 0.93,
        'clashscore': 2.1,
        'resolution': 3.2,
    }

    rows = format_quality_table_rows(metrics, 'cryoem')

    # Should have rows
    assert len(rows) > 0, "Should format some rows"

    # Check Map CC row
    cc_rows = [r for r in rows if r['label'] == 'Map CC']
    assert len(cc_rows) == 1, "Should have Map CC row"

    # Check Symmetry row
    sym_rows = [r for r in rows if r['label'] == 'Symmetry']
    assert len(sym_rows) == 1, "Should have Symmetry row"
    assert 'D7 (a)' in sym_rows[0]['value'], "Symmetry value should contain D7 (a)"
    assert '14 copies' in sym_rows[0]['value'], "Symmetry should include copies"
    assert 'CC' in sym_rows[0]['detail'], "Symmetry detail should include CC"

    print("  test_format_quality_table_cryoem PASSED")


def test_format_quality_table_missing_metrics():
    """Test that missing metrics don't create rows."""
    from knowledge.summary_display import format_quality_table_rows

    # Only resolution, no R-free or map_cc
    metrics = {
        'resolution': 2.5,
    }

    rows = format_quality_table_rows(metrics, 'xray')

    # Should have resolution but not R-free
    labels = [r['label'] for r in rows]
    assert 'R-free' not in labels, "Should not have R-free row without r_free metric"
    assert 'Resolution' in labels, "Should have Resolution row"

    print("  test_format_quality_table_missing_metrics PASSED")


def test_format_step_metric_refine():
    """Test step metric formatting for phenix.refine."""
    from knowledge.summary_display import format_step_metric

    metrics = {'r_free': 0.258}
    result = format_step_metric('phenix.refine', metrics)

    assert result is not None, "Should return formatted string"
    assert 'R-free' in result or 'r_free' in result.lower(), f"Should contain R-free, got {result}"
    assert '0.258' in result, f"Should contain 0.258, got {result}"

    print("  test_format_step_metric_refine PASSED")


def test_format_step_metric_rsr():
    """Test step metric formatting for phenix.real_space_refine."""
    from knowledge.summary_display import format_step_metric

    metrics = {'map_cc': 0.823}
    result = format_step_metric('phenix.real_space_refine', metrics)

    assert result is not None, "Should return formatted string"
    assert 'CC' in result, f"Should contain CC, got {result}"
    assert '0.823' in result, f"Should contain 0.823, got {result}"

    print("  test_format_step_metric_rsr PASSED")


def test_format_step_metric_map_symmetry():
    """Test step metric formatting for phenix.map_symmetry."""
    from knowledge.summary_display import format_step_metric

    metrics = {'symmetry_type': 'D7 (a)', 'ncs_copies': 14}
    result = format_step_metric('phenix.map_symmetry', metrics)

    assert result is not None, "Should return formatted string"
    assert 'D7 (a)' in result, f"Should contain D7 (a), got {result}"
    assert '14' in result, f"Should contain 14, got {result}"

    print("  test_format_step_metric_map_symmetry PASSED")


def test_format_step_metric_fallback():
    """Test step metric fallback when metrics missing."""
    from knowledge.summary_display import format_step_metric

    # No metrics provided
    result = format_step_metric('phenix.refine', {})

    # Should return fallback
    assert result is not None, "Should return fallback string"
    assert result == 'Refined' or 'Refined' in result, f"Should return 'Refined' fallback, got {result}"

    print("  test_format_step_metric_fallback PASSED")


def test_format_step_metric_unknown_program():
    """Test step metric for unknown program."""
    from knowledge.summary_display import format_step_metric

    result = format_step_metric('phenix.unknown_program', {})

    # Should use default
    assert result is not None, "Should return default string"

    print("  test_format_step_metric_unknown_program PASSED")


def test_experiment_type_filtering():
    """Test that experiment_type filters rows correctly."""
    from knowledge.summary_display import format_quality_table_rows

    metrics = {
        'r_free': 0.25,
        'map_cc': 0.85,
        'clashscore': 3.0,
    }

    # X-ray should have R-free but not Map CC
    xray_rows = format_quality_table_rows(metrics, 'xray')
    xray_labels = [r['label'] for r in xray_rows]
    assert 'R-free' in xray_labels, "X-ray should have R-free"
    assert 'Map CC' not in xray_labels, "X-ray should not have Map CC"

    # Cryo-EM should have Map CC but not R-free
    cryoem_rows = format_quality_table_rows(metrics, 'cryoem')
    cryoem_labels = [r['label'] for r in cryoem_rows]
    assert 'Map CC' in cryoem_labels, "Cryo-EM should have Map CC"
    assert 'R-free' not in cryoem_labels, "Cryo-EM should not have R-free"

    # Both should have clashscore (no experiment_type filter)
    assert 'Clashscore' in xray_labels, "X-ray should have Clashscore"
    assert 'Clashscore' in cryoem_labels, "Cryo-EM should have Clashscore"

    print("  test_experiment_type_filtering PASSED")


def test_get_metric_assessment():
    """Test metric quality assessment."""
    from knowledge.summary_display import get_metric_assessment

    # R-free tests
    assert get_metric_assessment('r_free', 0.20) == 'Good', "R-free 0.20 should be Good"
    assert get_metric_assessment('r_free', 0.25) == 'Acceptable', "R-free 0.25 should be Acceptable"
    assert get_metric_assessment('r_free', 0.35) == 'Needs Improvement', "R-free 0.35 should be Needs Improvement"

    # Clashscore tests (minimize)
    assert get_metric_assessment('clashscore', 1.5) == 'Good', "Clashscore 1.5 should be Good"

    # Map CC tests (maximize)
    assert get_metric_assessment('map_cc', 0.85) == 'Good', "Map CC 0.85 should be Good"

    print("  test_get_metric_assessment PASSED")


def run_all_tests():
    """Run all tests."""
    print("\nRunning summary_display tests...")

    test_quality_table_config_loading()
    test_step_metrics_config_loading()
    test_format_quality_table_xray()
    test_format_quality_table_cryoem()
    test_format_quality_table_missing_metrics()
    test_format_step_metric_refine()
    test_format_step_metric_rsr()
    test_format_step_metric_map_symmetry()
    test_format_step_metric_fallback()
    test_format_step_metric_unknown_program()
    test_experiment_type_filtering()
    test_get_metric_assessment()

    print("\nAll summary_display tests PASSED!")


if __name__ == "__main__":
    run_all_tests()
