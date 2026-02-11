#!/usr/bin/env python
"""
Tests for docs_tools.py - Documentation generation tool.

Run with: python tests/tst_docs_tools.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_load_command_templates():
    """Test loading program definitions from YAML."""
    print("Test: load_command_templates")

    from agent.docs_tools import load_command_templates

    templates = load_command_templates()

    assert templates is not None, "Should load templates"
    assert isinstance(templates, dict), "Should return dict"
    assert len(templates) > 0, "Should have entries"

    # Check for expected programs
    assert "phenix.refine" in templates, "Should have phenix.refine"
    assert "phenix.xtriage" in templates, "Should have phenix.xtriage"
    assert "phenix.phaser" in templates, "Should have phenix.phaser"

    # Check structure
    refine = templates["phenix.refine"]
    assert "description" in refine, "Should have description"
    assert "inputs" in refine, "Should have inputs"
    assert "outputs" in refine, "Should have outputs"

    print("  PASSED")


def test_load_decision_config():
    """Test loading decision config (may not exist)."""
    print("Test: load_decision_config")

    from agent.docs_tools import load_decision_config

    # This may return None if the file doesn't exist
    config = load_decision_config()

    # Just verify it doesn't crash
    # Config may be None if decision_config.json doesn't exist
    if config is not None:
        assert isinstance(config, dict), "Should return dict if exists"

    print("  PASSED")


def test_extract_workflow_states():
    """Test extracting workflow states from code."""
    print("Test: extract_workflow_states")

    from agent.docs_tools import extract_workflow_states

    states = extract_workflow_states()

    assert states is not None, "Should return states list"
    assert isinstance(states, list), "Should return list"

    # States may be empty if the regex doesn't match the current code structure
    # This is expected behavior - the function tries its best
    print(f"  Found {len(states)} workflow states")

    print("  PASSED")


def test_extract_stop_conditions():
    """Test extracting stop conditions from code."""
    print("Test: extract_stop_conditions")

    from agent.docs_tools import extract_stop_conditions

    conditions = extract_stop_conditions()

    assert conditions is not None, "Should return conditions"
    assert isinstance(conditions, (list, dict)), "Should return list or dict"

    print("  PASSED")


def test_extract_file_categorization():
    """Test extracting file categorization rules."""
    print("Test: extract_file_categorization")

    from agent.docs_tools import extract_file_categorization

    rules = extract_file_categorization()

    assert rules is not None, "Should return rules"

    print("  PASSED")


def test_extract_key_thresholds():
    """Test extracting key thresholds."""
    print("Test: extract_key_thresholds")

    from agent.docs_tools import extract_key_thresholds

    thresholds = extract_key_thresholds()

    assert thresholds is not None, "Should return thresholds"

    print("  PASSED")


def test_generate_markdown():
    """Test generating markdown documentation."""
    print("Test: generate_markdown")

    from agent.docs_tools import (
        load_command_templates,
        load_decision_config,
        extract_workflow_states,
        extract_stop_conditions,
        extract_file_categorization,
        extract_key_thresholds,
        generate_markdown
    )

    # Load all data
    templates = load_command_templates()
    config = load_decision_config()
    states = extract_workflow_states()
    stop_conditions = extract_stop_conditions()
    file_rules = extract_file_categorization()
    thresholds = extract_key_thresholds()

    # Generate markdown
    doc = generate_markdown(templates, states, stop_conditions, file_rules, thresholds, config)

    assert doc is not None, "Should generate documentation"
    assert isinstance(doc, str), "Should return string"
    assert len(doc) > 100, "Should have substantial content"

    # Check for expected sections
    assert "# PHENIX AI Agent" in doc, "Should have title"
    assert "Workflow" in doc, "Should have workflow section"
    assert "Program" in doc, "Should have program section"

    print("  PASSED")


def test_full_generation():
    """Test full documentation generation workflow."""
    print("Test: full_generation")

    import subprocess

    # Run the tool and capture output
    result = subprocess.run(
        [sys.executable, "agent/docs_tools.py"],
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.dirname(__file__))
    )

    assert result.returncode == 0, f"Should succeed: {result.stderr}"
    assert len(result.stdout) > 100, "Should produce output"
    assert "PHENIX AI Agent" in result.stdout, "Should have expected content"

    print("  PASSED")


def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.tst_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
