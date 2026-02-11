#!/usr/bin/env python
"""
Template for PHENIX AI Agent tests.

This file shows the standard test pattern used in this project.
Copy this file and modify for new test modules.

Key principles:
1. Each test is a simple function named test_*()
2. Tests use assert statements - crash on failure, pass silently
3. Each test prints its name and "PASSED" on success
4. run_all_tests() calls all test functions
5. Can be run standalone or via run_all_tests.py

Run with: python tests/tst_template.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# =============================================================================
# Test Functions
# =============================================================================

def test_basic_functionality():
    """Test description goes here."""
    print("Test: basic_functionality")

    # Setup
    value = 1 + 1

    # Assertions - use assert, crash on failure
    assert value == 2, "Basic math should work"
    assert value != 3, "Should not equal 3"

    print("  PASSED")


def test_with_expected_exception():
    """Test that an exception is raised."""
    print("Test: with_expected_exception")

    exception_raised = False
    try:
        int("not a number")
    except ValueError:
        exception_raised = True

    assert exception_raised, "Should raise ValueError for invalid int"

    print("  PASSED")


def test_with_setup_and_teardown():
    """Test with temporary resources."""
    print("Test: with_setup_and_teardown")

    import tempfile
    import shutil

    # Setup
    temp_dir = tempfile.mkdtemp()
    try:
        # Test logic
        test_file = os.path.join(temp_dir, "test.txt")
        with open(test_file, 'w') as f:
            f.write("test content")

        assert os.path.exists(test_file), "File should exist"

        with open(test_file, 'r') as f:
            content = f.read()
        assert content == "test content", "Content should match"

    finally:
        # Teardown - always clean up
        shutil.rmtree(temp_dir)

    print("  PASSED")


def test_checking_list_contents():
    """Test list membership and properties."""
    print("Test: checking_list_contents")

    items = ["apple", "banana", "cherry"]

    assert len(items) == 3, "Should have 3 items"
    assert "banana" in items, "Should contain banana"
    assert "grape" not in items, "Should not contain grape"
    assert items[0] == "apple", "First item should be apple"

    print("  PASSED")


def test_checking_dict_contents():
    """Test dictionary contents."""
    print("Test: checking_dict_contents")

    data = {"name": "test", "value": 42}

    assert "name" in data, "Should have 'name' key"
    assert data["name"] == "test", "Name should be 'test'"
    assert data.get("missing") is None, "Missing key should return None"
    assert data.get("missing", "default") == "default", "Should return default"

    print("  PASSED")


# =============================================================================
# Test Runner
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.tst_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
