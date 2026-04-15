#!/usr/bin/env python
"""
Tests for knowledge/program_registration.py
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_trackable_programs():
    """Test that trackable programs are loaded from YAML."""
    from knowledge.program_registration import get_trackable_programs

    trackable = get_trackable_programs()

    # Should have at least the known run_once programs
    assert len(trackable) >= 3, f"Should have at least 3 trackable programs, got {len(trackable)}"

    # Check specific programs
    assert "phenix.xtriage" in trackable, "xtriage should be trackable"
    assert "phenix.mtriage" in trackable, "mtriage should be trackable"
    assert "phenix.map_symmetry" in trackable, "map_symmetry should be trackable"

    print("  test_trackable_programs PASSED")


def test_done_flag_format():
    """Test that done flags have correct format."""
    from knowledge.program_registration import get_trackable_programs

    trackable = get_trackable_programs()

    for prog_name, info in trackable.items():
        # Check short_name
        assert 'short_name' in info, f"{prog_name} should have short_name"
        assert 'phenix.' not in info['short_name'], f"short_name should not contain 'phenix.'"

        # Check done_flag
        assert 'done_flag' in info, f"{prog_name} should have done_flag"
        assert info['done_flag'].endswith('_done'), f"done_flag should end with '_done'"
        assert info['done_flag'] == f"{info['short_name']}_done", \
            f"done_flag format mismatch for {prog_name}: {info['done_flag']} != {info['short_name']}_done"
        # Note: This holds for current run_once programs (xtriage, mtriage,
        # map_symmetry) but the flag name actually comes from YAML
        # done_tracking.flag, so irregular names are supported.

    print("  test_done_flag_format PASSED")


def test_initial_flags():
    """Test that initial flags are all False."""
    from knowledge.program_registration import get_initial_history_flags

    flags = get_initial_history_flags()

    # Should have flags
    assert len(flags) >= 3, f"Should have at least 3 flags, got {len(flags)}"

    # All should be False
    for flag_name, value in flags.items():
        assert value == False, f"{flag_name} should be False initially"

    print("  test_initial_flags PASSED")


def test_detect_xtriage():
    """Test detection of xtriage in history."""
    from knowledge.program_registration import detect_programs_in_history, get_initial_history_flags

    history = [
        {"program": "phenix.xtriage", "command": "phenix.xtriage data.mtz", "result": "SUCCESS"},
    ]

    flags = get_initial_history_flags()
    detected = detect_programs_in_history(history, flags)

    assert detected.get("xtriage_done") == True, "Should detect xtriage"
    assert detected.get("mtriage_done") == False, "Should not detect mtriage"

    print("  test_detect_xtriage PASSED")


def test_detect_mtriage():
    """Test detection of mtriage in history."""
    from knowledge.program_registration import detect_programs_in_history, get_initial_history_flags

    history = [
        {"program": "phenix.mtriage", "command": "phenix.mtriage map.mrc", "result": "SUCCESS"},
    ]

    flags = get_initial_history_flags()
    detected = detect_programs_in_history(history, flags)

    assert detected.get("mtriage_done") == True, "Should detect mtriage"
    assert detected.get("xtriage_done") == False, "Should not detect xtriage"

    print("  test_detect_mtriage PASSED")


def test_detect_map_symmetry():
    """Test detection of map_symmetry in history."""
    from knowledge.program_registration import detect_programs_in_history, get_initial_history_flags

    history = [
        {"program": "phenix.map_symmetry", "command": "phenix.map_symmetry map.mrc", "result": "SUCCESS"},
    ]

    flags = get_initial_history_flags()
    detected = detect_programs_in_history(history, flags)

    assert detected.get("map_symmetry_done") == True, "Should detect map_symmetry"

    print("  test_detect_map_symmetry PASSED")


def test_detect_multiple():
    """Test detection of multiple programs in history."""
    from knowledge.program_registration import detect_programs_in_history, get_initial_history_flags

    history = [
        {"program": "phenix.mtriage", "command": "phenix.mtriage map.mrc"},
        {"program": "phenix.map_symmetry", "command": "phenix.map_symmetry map.mrc"},
    ]

    flags = get_initial_history_flags()
    detected = detect_programs_in_history(history, flags)

    assert detected.get("mtriage_done") == True, "Should detect mtriage"
    assert detected.get("map_symmetry_done") == True, "Should detect map_symmetry"
    assert detected.get("xtriage_done") == False, "Should not detect xtriage"

    print("  test_detect_multiple PASSED")


def test_detection_patterns():
    """Test that various pattern formats are detected."""
    from knowledge.program_registration import detect_programs_in_history, get_initial_history_flags

    # Test different ways xtriage might appear
    test_cases = [
        [{"program": "phenix.xtriage"}],
        [{"program": "xtriage"}],
        [{"command": "phenix.xtriage data.mtz"}],
    ]

    for history in test_cases:
        flags = get_initial_history_flags()
        detected = detect_programs_in_history(history, flags)
        assert detected.get("xtriage_done") == True, f"Should detect xtriage in {history}"

    print("  test_detection_patterns PASSED")


def test_empty_history():
    """Test with empty history."""
    from knowledge.program_registration import detect_programs_in_history, get_initial_history_flags

    flags = get_initial_history_flags()
    detected = detect_programs_in_history([], flags)

    # All should remain False
    for flag_name, value in detected.items():
        assert value == False, f"{flag_name} should be False with empty history"

    print("  test_empty_history PASSED")


def test_get_all_done_flags():
    """Test getting list of all done flags."""
    from knowledge.program_registration import get_all_done_flags

    flags = get_all_done_flags()

    assert isinstance(flags, list), "Should return a list"
    assert len(flags) >= 3, "Should have at least 3 flags"
    assert "xtriage_done" in flags, "Should include xtriage_done"
    assert "mtriage_done" in flags, "Should include mtriage_done"

    print("  test_get_all_done_flags PASSED")


def test_is_program_trackable():
    """Test checking if program is trackable."""
    from knowledge.program_registration import is_program_trackable

    assert is_program_trackable("phenix.xtriage") == True
    assert is_program_trackable("phenix.mtriage") == True
    assert is_program_trackable("phenix.refine") == False  # No run_once
    assert is_program_trackable("phenix.unknown") == False

    print("  test_is_program_trackable PASSED")


def test_get_done_flag_for_program():
    """Test getting done flag name for a program."""
    from knowledge.program_registration import get_done_flag_for_program

    assert get_done_flag_for_program("phenix.xtriage") == "xtriage_done"
    assert get_done_flag_for_program("phenix.mtriage") == "mtriage_done"
    assert get_done_flag_for_program("phenix.map_symmetry") == "map_symmetry_done"
    assert get_done_flag_for_program("phenix.refine") is None

    print("  test_get_done_flag_for_program PASSED")


def test_workflow_state_integration():
    """Test that workflow_state uses auto-registration."""
    from agent.workflow_state import _analyze_history

    history = [
        {"program": "phenix.xtriage", "command": "phenix.xtriage data.mtz", "result": "SUCCESS"},
        {"program": "phenix.mtriage", "command": "phenix.mtriage map.mrc", "result": "SUCCESS"},
    ]

    info = _analyze_history(history)

    # Check auto-generated flags
    assert info.get("xtriage_done") == True, "Should detect xtriage via auto-registration"
    assert info.get("mtriage_done") == True, "Should detect mtriage via auto-registration"

    # Check that complex flags still work
    assert "phaser_done" in info, "Should still have phaser_done"
    assert "refine_count" in info, "Should still have refine_count"

    print("  test_workflow_state_integration PASSED")


def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.tst_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
