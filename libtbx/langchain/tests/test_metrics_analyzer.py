#!/usr/bin/env python
"""
Unit tests for metrics_analyzer.py

Tests:
- derive_metrics_from_history with various history formats
- analyze_metrics_trend for plateau detection
- Edge cases: empty history, missing metrics, non-monotonic R-free
"""

from __future__ import absolute_import, division, print_function
import sys
import os

# Add agent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'agent'))

from metrics_analyzer import (
    derive_metrics_from_history,
    analyze_metrics_trend,
    get_latest_resolution,
    get_best_r_free,
    get_latest_r_free,
    format_trend_for_prompt
)


def test_empty_history():
    """Test with empty history."""
    print("Test: empty_history")
    
    metrics = derive_metrics_from_history([])
    assert metrics == [], "Expected empty list for empty history"
    
    trend = analyze_metrics_trend([])
    assert trend["should_stop"] == False
    assert trend["trend_summary"] == "No refinement history"
    
    print("  PASSED")


def test_derive_metrics_dict_format():
    """Test metrics extraction from dict history entries."""
    print("Test: derive_metrics_dict_format")
    
    history = [
        {
            "cycle_number": 1,
            "program": "phenix.xtriage",
            "analysis": {"resolution": 2.5},
            "result": "SUCCESS"
        },
        {
            "cycle_number": 2,
            "program": "phenix.phaser",
            "analysis": {"tfz": 24.5, "llg": 544},
            "result": "SUCCESS TFZ=24.5"
        },
        {
            "cycle_number": 3,
            "program": "phenix.refine",
            "analysis": {"r_free": 0.35, "r_work": 0.30},
            "result": "R-free=0.35"
        }
    ]
    
    metrics = derive_metrics_from_history(history)
    
    assert len(metrics) == 3, "Expected 3 metrics entries"
    assert metrics[0]["resolution"] == 2.5
    assert metrics[1]["tfz"] == 24.5
    assert metrics[1]["llg"] == 544
    assert metrics[2]["r_free"] == 0.35
    assert metrics[2]["r_work"] == 0.30
    
    print("  PASSED")


def test_derive_metrics_string_format():
    """Test metrics extraction from legacy string history entries."""
    print("Test: derive_metrics_string_format")
    
    history = [
        "phenix.xtriage resolution=2.5",
        "phenix.phaser TFZ=24.5 LLG=544",
        "phenix.refine R-free: 0.35 R-work: 0.30"
    ]
    
    metrics = derive_metrics_from_history(history)
    
    assert len(metrics) == 3
    assert metrics[0]["resolution"] == 2.5
    assert metrics[1]["tfz"] == 24.5
    assert metrics[2]["r_free"] == 0.35
    
    print("  PASSED")


def test_derive_metrics_fallback_to_result():
    """Test that metrics are extracted from result field when analysis is empty."""
    print("Test: derive_metrics_fallback_to_result")
    
    history = [
        {
            "cycle_number": 1,
            "program": "phenix.refine",
            "analysis": {},  # Empty analysis
            "result": "Final R-free: 0.28, R-work: 0.24"
        }
    ]
    
    metrics = derive_metrics_from_history(history)
    
    assert len(metrics) == 1
    assert metrics[0]["r_free"] == 0.28
    assert metrics[0]["r_work"] == 0.24
    
    print("  PASSED")


def test_trend_improving():
    """Test trend analysis with improving R-free."""
    print("Test: trend_improving")
    
    metrics = [
        {"cycle": 1, "program": "phenix.refine", "r_free": 0.40},
        {"cycle": 2, "program": "phenix.refine", "r_free": 0.35},
        {"cycle": 3, "program": "phenix.refine", "r_free": 0.30},
    ]
    
    trend = analyze_metrics_trend(metrics)
    
    assert trend["should_stop"] == False
    assert trend["improvement_rate"] > 10  # >10% improvement
    assert "improving" in trend["trend_summary"].lower()
    assert trend["consecutive_refines"] == 3
    
    print("  PASSED")


def test_trend_plateau():
    """Test plateau detection (0.3% improvement for 3+ cycles)."""
    print("Test: trend_plateau")
    
    # Need 4 values to detect plateau over last 3 improvements
    # Each improvement < 0.3%
    metrics = [
        {"cycle": 1, "program": "phenix.refine", "r_free": 0.300},
        {"cycle": 2, "program": "phenix.refine", "r_free": 0.2995},  # 0.17% improvement
        {"cycle": 3, "program": "phenix.refine", "r_free": 0.2991},  # 0.13% improvement
        {"cycle": 4, "program": "phenix.refine", "r_free": 0.2988},  # 0.10% improvement
        {"cycle": 5, "program": "phenix.refine", "r_free": 0.2986},  # 0.07% improvement
    ]
    
    trend = analyze_metrics_trend(metrics)
    
    assert trend["should_stop"] == True
    assert "PLATEAU" in trend["reason"]
    assert trend["recommendation"] == "stop"
    
    print("  PASSED")


def test_trend_success():
    """Test success detection (R-free well below target)."""
    print("Test: trend_success")
    
    # Test 1: Success without validation - should recommend validation, not stop
    metrics_no_validation = [
        {"cycle": 1, "program": "phenix.refine", "r_free": 0.30},
        {"cycle": 2, "program": "phenix.refine", "r_free": 0.22},
    ]
    
    trend = analyze_metrics_trend(metrics_no_validation, resolution=2.5)
    
    assert trend["should_stop"] == False, "Should not stop without validation"
    assert trend["recommendation"] == "validate", "Should recommend validation"
    assert trend["suggest_validation"] == True
    
    # Test 2: Success with validation - should stop
    metrics_with_validation = [
        {"cycle": 1, "program": "phenix.refine", "r_free": 0.30},
        {"cycle": 2, "program": "phenix.refine", "r_free": 0.22},
        {"cycle": 3, "program": "phenix.molprobity"},  # Validation done
    ]
    
    trend2 = analyze_metrics_trend(metrics_with_validation, resolution=2.5)
    
    assert trend2["should_stop"] == True
    assert "SUCCESS" in trend2["reason"]
    
    print("  PASSED")


def test_trend_borderline():
    """Test borderline R-free (close to target but not auto-stop)."""
    print("Test: trend_borderline")
    
    # Resolution 2.5Å -> target 0.25, success threshold 0.23
    # R-free 0.24 is borderline - should NOT auto-stop, let LLM decide
    metrics = [
        {"cycle": 1, "program": "phenix.refine", "r_free": 0.30},
        {"cycle": 2, "program": "phenix.refine", "r_free": 0.24},
    ]
    
    trend = analyze_metrics_trend(metrics, resolution=2.5)
    
    # Should NOT auto-stop for borderline cases
    assert trend["should_stop"] == False
    assert "could try autobuild" in trend["trend_summary"]
    
    print("  PASSED")


def test_trend_excessive_refines():
    """Test excessive refinement detection (8+ cycles)."""
    print("Test: trend_excessive_refines")
    
    # Need 8+ consecutive refines to trigger excessive
    metrics = [
        {"cycle": i, "program": "phenix.refine", "r_free": 0.30 - i*0.005}
        for i in range(9)
    ]
    
    trend = analyze_metrics_trend(metrics)
    
    assert trend["should_stop"] == True
    assert "EXCESSIVE" in trend["reason"]
    assert trend["consecutive_refines"] == 9
    
    print("  PASSED")


def test_trend_worsening():
    """Test detection of worsening R-free."""
    print("Test: trend_worsening")
    
    metrics = [
        {"cycle": 1, "program": "phenix.refine", "r_free": 0.28},
        {"cycle": 2, "program": "phenix.refine", "r_free": 0.30},  # Got worse
    ]
    
    trend = analyze_metrics_trend(metrics)
    
    assert trend["should_stop"] == False  # Don't auto-stop, but warn
    assert trend["improvement_rate"] < 0
    assert trend["recommendation"] == "consider_stopping"
    
    print("  PASSED")


def test_trend_non_monotonic():
    """Test with non-monotonic R-free (temporary increase)."""
    print("Test: trend_non_monotonic")
    
    metrics = [
        {"cycle": 1, "program": "phenix.refine", "r_free": 0.35},
        {"cycle": 2, "program": "phenix.refine", "r_free": 0.32},
        {"cycle": 3, "program": "phenix.refine", "r_free": 0.34},  # Temporary increase
        {"cycle": 4, "program": "phenix.refine", "r_free": 0.29},  # Recovered
    ]
    
    trend = analyze_metrics_trend(metrics)
    
    # Should see improvement from 0.34 to 0.29
    assert trend["improvement_rate"] > 10
    assert trend["should_stop"] == False
    
    print("  PASSED")


def test_trend_mixed_programs():
    """Test with mixed programs (non-refine cycles interleaved)."""
    print("Test: trend_mixed_programs")
    
    metrics = [
        {"cycle": 1, "program": "phenix.xtriage", "resolution": 2.5},
        {"cycle": 2, "program": "phenix.phaser", "tfz": 24},
        {"cycle": 3, "program": "phenix.refine", "r_free": 0.35},
        {"cycle": 4, "program": "phenix.ligandfit", "r_free": None},  # No R-free
        {"cycle": 5, "program": "phenix.refine", "r_free": 0.30},
    ]
    
    trend = analyze_metrics_trend(metrics)
    
    # Should only count 2 refine cycles for R-free trend
    assert len(trend["r_free_trend"]) == 2
    # Consecutive refines should be 1 (only the last one)
    assert trend["consecutive_refines"] == 1
    
    print("  PASSED")


def test_trend_cryoem():
    """Test cryo-EM trend analysis (CC-based)."""
    print("Test: trend_cryoem")
    
    metrics = [
        {"cycle": 1, "program": "phenix.mtriage", "resolution": 3.0},
        {"cycle": 2, "program": "phenix.real_space_refine", "map_cc": 0.65},
        {"cycle": 3, "program": "phenix.real_space_refine", "map_cc": 0.72},
        {"cycle": 4, "program": "phenix.real_space_refine", "map_cc": 0.78},
    ]
    
    trend = analyze_metrics_trend(metrics, experiment_type="cryoem")
    
    assert trend["should_stop"] == True
    assert "SUCCESS" in trend["reason"]
    assert "0.75" in trend["reason"] or "CC" in trend["reason"]
    
    print("  PASSED")


def test_dynamic_resolution_target():
    """Test dynamic R-free target based on resolution."""
    print("Test: dynamic_resolution_target")
    
    # Low resolution (3.5Å) -> target 0.30, success threshold 0.28
    # R-free 0.27 should trigger validation recommendation (not stop without validation)
    metrics1 = [
        {"cycle": 1, "program": "phenix.refine", "r_free": 0.27},
    ]
    trend1 = analyze_metrics_trend(metrics1, resolution=3.5)
    # Without validation, should recommend validate, not stop
    assert trend1["should_stop"] == False, "Should not stop without validation"
    assert trend1["recommendation"] == "validate"
    
    # With validation, should stop
    metrics1_val = [
        {"cycle": 1, "program": "phenix.refine", "r_free": 0.27},
        {"cycle": 2, "program": "phenix.molprobity"},
    ]
    trend1_val = analyze_metrics_trend(metrics1_val, resolution=3.5)
    assert trend1_val["should_stop"] == True  # 0.27 < 0.28 (success threshold) + validated
    
    # Low resolution (3.5Å) -> R-free 0.29 is borderline (below target but above success threshold)
    metrics1b = [
        {"cycle": 1, "program": "phenix.refine", "r_free": 0.29},
    ]
    trend1b = analyze_metrics_trend(metrics1b, resolution=3.5)
    assert trend1b["should_stop"] == False  # 0.29 > 0.28 (borderline, not auto-stop)
    assert "could try autobuild" in trend1b["trend_summary"]
    
    # High resolution (1.5Å) -> target 0.20, success threshold 0.18
    metrics2 = [{"cycle": 1, "program": "phenix.refine", "r_free": 0.22}]
    trend2 = analyze_metrics_trend(metrics2, resolution=1.5)
    assert trend2["should_stop"] == False  # 0.22 > 0.20 (above target)
    
    print("  PASSED")


def test_get_latest_resolution():
    """Test resolution extraction from metrics."""
    print("Test: get_latest_resolution")
    
    metrics = [
        {"cycle": 1, "program": "phenix.xtriage", "resolution": 2.5},
        {"cycle": 2, "program": "phenix.refine", "resolution": None},
    ]
    
    res = get_latest_resolution(metrics)
    assert res == 2.5
    
    print("  PASSED")


def test_get_best_r_free():
    """Test best R-free extraction."""
    print("Test: get_best_r_free")
    
    metrics = [
        {"cycle": 1, "program": "phenix.refine", "r_free": 0.35},
        {"cycle": 2, "program": "phenix.refine", "r_free": 0.28},
        {"cycle": 3, "program": "phenix.refine", "r_free": 0.30},  # Got worse
    ]
    
    best = get_best_r_free(metrics)
    assert best == 0.28
    
    print("  PASSED")


def test_format_trend_for_prompt():
    """Test prompt formatting."""
    print("Test: format_trend_for_prompt")
    
    trend = {
        "should_stop": True,
        "reason": "PLATEAU: R-free improvement stalled",
        "trend_summary": "R-free: 0.30 -> 0.29 -> 0.289",
        "consecutive_refines": 3,
        "consecutive_rsr": 0,
        "recommendation": "stop"
    }
    
    formatted = format_trend_for_prompt(trend)
    
    assert "REFINEMENT PROGRESS" in formatted
    assert "0.30" in formatted
    assert "STOP RECOMMENDED" in formatted
    assert "PLATEAU" in formatted
    
    print("  PASSED")


def test_first_refinement():
    """Test trend with only one refinement cycle."""
    print("Test: first_refinement")
    
    metrics = [
        {"cycle": 1, "program": "phenix.xtriage", "resolution": 2.5},
        {"cycle": 2, "program": "phenix.refine", "r_free": 0.35},
    ]
    
    trend = analyze_metrics_trend(metrics)
    
    assert trend["should_stop"] == False
    assert "first refinement" in trend["trend_summary"].lower()
    assert trend["improvement_rate"] == 0.0  # No previous to compare
    
    print("  PASSED")


def test_autobuild_rfree_extraction():
    """Test extraction of R-free from autobuild table format."""
    print("Test: autobuild_rfree_extraction")
    
    # Simulate autobuild output in result text
    autobuild_result = """
    Running autobuild...
    
    SOLUTION  CYCLE     R        RFREE     BUILT   PLACED
    1         1      0.21        0.25      0       0
    
    Best model: overall_best.pdb
    """
    
    history = [
        {
            "cycle_number": 1,
            "program": "phenix.autobuild",
            "result": autobuild_result,
        }
    ]
    
    metrics = derive_metrics_from_history(history)
    assert len(metrics) == 1
    assert metrics[0]["program"] == "phenix.autobuild"
    assert metrics[0]["r_free"] == 0.25, f"Expected R-free 0.25, got {metrics[0]['r_free']}"
    
    print("  PASSED")


def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("METRICS ANALYZER TESTS")
    print("=" * 60)
    
    test_empty_history()
    test_derive_metrics_dict_format()
    test_derive_metrics_string_format()
    test_derive_metrics_fallback_to_result()
    test_trend_improving()
    test_trend_plateau()
    test_trend_success()
    test_trend_borderline()
    test_trend_excessive_refines()
    test_trend_worsening()
    test_trend_non_monotonic()
    test_trend_mixed_programs()
    test_trend_cryoem()
    test_dynamic_resolution_target()
    test_get_latest_resolution()
    test_get_best_r_free()
    test_format_trend_for_prompt()
    test_first_refinement()
    test_autobuild_rfree_extraction()
    
    print("=" * 60)
    print("ALL METRICS ANALYZER TESTS PASSED!")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
