#!/usr/bin/env python
"""
Tests for session_tools.py - Session management tool.

Run with: python tests/tst_session_tools.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import json
import tempfile
import shutil

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def create_test_session(session_dir):
    """Create a test session with some cycles."""
    session_file = os.path.join(session_dir, "agent_session.json")

    session_data = {
        "session_id": "2026-01-21_10-00-00",
        "project_advice": "Test advice",
        "original_files": ["/path/to/data.mtz", "/path/to/seq.fa"],
        "cycles": [
            {
                "cycle_number": 1,
                "program": "phenix.xtriage",
                "command": "phenix.xtriage data.mtz",
                "result": "SUCCESS: Analysis complete",
                "output_files": ["/path/to/xtriage.log"],
                "metrics": {"resolution": 2.5}
            },
            {
                "cycle_number": 2,
                "program": "phenix.phaser",
                "command": "phenix.phaser data.mtz model.pdb",
                "result": "SUCCESS: MR solution found",
                "output_files": ["/path/to/phaser.pdb"],
                "metrics": {"llg": 500}
            },
            {
                "cycle_number": 3,
                "program": "phenix.refine",
                "command": "phenix.refine phaser.pdb data.mtz",
                "result": "SUCCESS: Refinement complete",
                "output_files": ["/path/to/refine.pdb"],
                "metrics": {"r_free": 0.28, "r_work": 0.24}
            }
        ],
        "resolution": 2.5,
        "experiment_type": "xray"
    }

    with open(session_file, 'w') as f:
        json.dump(session_data, f, indent=2)

    return session_file, session_data


def test_load_session():
    """Test loading a session file."""
    print("Test: load_session")

    from agent.session_tools import load_session

    temp_dir = tempfile.mkdtemp()
    try:
        session_file, original_data = create_test_session(temp_dir)

        data, loaded_file = load_session(temp_dir)

        assert data is not None, "Should load session data"
        assert loaded_file == session_file, "Should return correct file path"
        assert data["session_id"] == original_data["session_id"], "Should have correct session_id"
        assert len(data["cycles"]) == 3, "Should have 3 cycles"
    finally:
        shutil.rmtree(temp_dir)

    print("  PASSED")


def test_load_session_missing():
    """Test loading a non-existent session."""
    print("Test: load_session_missing")

    from agent.session_tools import load_session

    temp_dir = tempfile.mkdtemp()
    try:
        data, session_file = load_session(temp_dir)

        assert data is None, "Should return None for missing session"
    finally:
        shutil.rmtree(temp_dir)

    print("  PASSED")


def test_save_session():
    """Test saving a session file."""
    print("Test: save_session")

    from agent.session_tools import load_session, save_session

    temp_dir = tempfile.mkdtemp()
    try:
        session_file, original_data = create_test_session(temp_dir)

        # Modify and save
        data, _ = load_session(temp_dir)
        data["project_advice"] = "Modified advice"
        save_session(data, session_file)

        # Reload and verify
        data2, _ = load_session(temp_dir)
        assert data2["project_advice"] == "Modified advice", "Should save modifications"
    finally:
        shutil.rmtree(temp_dir)

    print("  PASSED")


def test_show_session():
    """Test displaying session summary."""
    print("Test: show_session")

    from agent.session_tools import load_session, show_session

    temp_dir = tempfile.mkdtemp()
    try:
        create_test_session(temp_dir)
        data, _ = load_session(temp_dir)

        # Capture output (show_session prints to stdout)
        # Just verify it doesn't crash
        show_session(data, detailed=False)
        show_session(data, detailed=True)
    finally:
        shutil.rmtree(temp_dir)

    print("  PASSED")


def test_remove_cycles():
    """Test removing cycles from session."""
    print("Test: remove_cycles")

    from agent.session_tools import load_session, save_session, remove_last_cycles

    temp_dir = tempfile.mkdtemp()
    try:
        session_file, _ = create_test_session(temp_dir)

        data, _ = load_session(temp_dir)
        assert len(data["cycles"]) == 3, "Should start with 3 cycles"

        # Remove last 2 cycles
        data, removed = remove_last_cycles(data, 2)
        assert removed == 2, "Should remove 2 cycles"
        save_session(data, session_file)

        # Reload and verify
        data2, _ = load_session(temp_dir)
        assert len(data2["cycles"]) == 1, "Should have 1 cycle after removing 2"
        assert data2["cycles"][0]["cycle_number"] == 1, "Should keep first cycle"
    finally:
        shutil.rmtree(temp_dir)

    print("  PASSED")


def test_remove_cycles_all():
    """Test removing all cycles (more than exist)."""
    print("Test: remove_cycles_all")

    from agent.session_tools import load_session, remove_last_cycles

    temp_dir = tempfile.mkdtemp()
    try:
        create_test_session(temp_dir)
        data, _ = load_session(temp_dir)

        # Remove more cycles than exist
        data, removed = remove_last_cycles(data, 10)

        assert len(data["cycles"]) == 0, "Should have 0 cycles"
    finally:
        shutil.rmtree(temp_dir)

    print("  PASSED")


def test_reset_session():
    """Test resetting a session."""
    print("Test: reset_session")

    from agent.session_tools import load_session, save_session, reset_session

    temp_dir = tempfile.mkdtemp()
    try:
        session_file, _ = create_test_session(temp_dir)

        data, _ = load_session(temp_dir)
        assert len(data["cycles"]) == 3, "Should start with 3 cycles"

        # Reset session
        data = reset_session(data)
        save_session(data, session_file)

        # Reload and verify
        data2, _ = load_session(temp_dir)
        assert len(data2["cycles"]) == 0, "Should have 0 cycles after reset"
        assert data2["project_advice"] == "Test advice", "Should preserve project_advice"
    finally:
        shutil.rmtree(temp_dir)

    print("  PASSED")


def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.tst_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
