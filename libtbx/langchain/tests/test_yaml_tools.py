#!/usr/bin/env python
"""
Tests for yaml_tools.py - YAML configuration management tool.

Run with: python tests/test_yaml_tools.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import tempfile
import shutil

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_find_yaml_files():
    """Test finding YAML files in knowledge directory."""
    print("Test: find_yaml_files")

    from agent.yaml_tools import find_yaml_files

    # Find from agent directory (should find ../knowledge/)
    agent_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "agent")
    yaml_files = find_yaml_files(agent_dir)

    assert len(yaml_files) > 0, "Should find at least one YAML file"

    # Check that expected files are found
    filenames = [os.path.basename(f) for f in yaml_files]
    assert "programs.yaml" in filenames, "Should find programs.yaml"
    assert "workflows.yaml" in filenames, "Should find workflows.yaml"
    assert "metrics.yaml" in filenames, "Should find metrics.yaml"
    assert "file_categories.yaml" in filenames, "Should find file_categories.yaml"

    print("  PASSED")


def test_load_yaml_file():
    """Test loading YAML files."""
    print("Test: load_yaml_file")

    from agent.yaml_tools import load_yaml_file, find_yaml_files

    agent_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "agent")
    yaml_files = find_yaml_files(agent_dir)

    # Load programs.yaml
    programs_file = [f for f in yaml_files if "programs.yaml" in f][0]
    data, error = load_yaml_file(programs_file)

    assert error is None, f"Should load without error: {error}"
    assert data is not None, "Should return data"
    assert isinstance(data, dict), "Should return dict"
    assert len(data) > 0, "Should have entries"

    # Check for expected programs
    assert "phenix.refine" in data, "Should have phenix.refine"
    assert "phenix.xtriage" in data, "Should have phenix.xtriage"

    print("  PASSED")


def test_load_yaml_file_invalid():
    """Test loading invalid YAML file."""
    print("Test: load_yaml_file_invalid")

    from agent.yaml_tools import load_yaml_file

    # Test with non-existent file
    data, error = load_yaml_file("/nonexistent/file.yaml")
    assert error is not None, "Should return error for missing file"
    assert data is None, "Should return None for missing file"

    # Test with invalid YAML
    temp_dir = tempfile.mkdtemp()
    try:
        invalid_yaml = os.path.join(temp_dir, "invalid.yaml")
        with open(invalid_yaml, 'w') as f:
            f.write("invalid: yaml: content: [")  # Invalid YAML

        data, error = load_yaml_file(invalid_yaml)
        assert error is not None, "Should return error for invalid YAML"
    finally:
        shutil.rmtree(temp_dir)

    print("  PASSED")


def test_validate_programs_yaml():
    """Test validation of programs.yaml structure."""
    print("Test: validate_programs_yaml")

    from agent.yaml_tools import validate_yaml_structure, load_yaml_file, find_yaml_files

    agent_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "agent")
    yaml_files = find_yaml_files(agent_dir)

    programs_file = [f for f in yaml_files if "programs.yaml" in f][0]
    data, _ = load_yaml_file(programs_file)

    issues = validate_yaml_structure(data, programs_file)

    # Should have no errors (warnings are OK)
    errors = [i for i in issues if i[0] == 'error']
    assert len(errors) == 0, f"Should have no errors: {errors}"

    print("  PASSED")


def test_validate_workflows_yaml():
    """Test validation of workflows.yaml structure."""
    print("Test: validate_workflows_yaml")

    from agent.yaml_tools import validate_yaml_structure, load_yaml_file, find_yaml_files

    agent_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "agent")
    yaml_files = find_yaml_files(agent_dir)

    workflows_file = [f for f in yaml_files if "workflows.yaml" in f][0]
    data, _ = load_yaml_file(workflows_file)

    issues = validate_yaml_structure(data, workflows_file)

    # Should have no errors
    errors = [i for i in issues if i[0] == 'error']
    assert len(errors) == 0, f"Should have no errors: {errors}"

    print("  PASSED")


def test_validate_metrics_yaml():
    """Test validation of metrics.yaml structure."""
    print("Test: validate_metrics_yaml")

    from agent.yaml_tools import validate_yaml_structure, load_yaml_file, find_yaml_files

    agent_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "agent")
    yaml_files = find_yaml_files(agent_dir)

    metrics_file = [f for f in yaml_files if "metrics.yaml" in f][0]
    data, _ = load_yaml_file(metrics_file)

    issues = validate_yaml_structure(data, metrics_file)

    # Should have no errors
    errors = [i for i in issues if i[0] == 'error']
    assert len(errors) == 0, f"Should have no errors: {errors}"

    print("  PASSED")


def test_validate_file_categories_yaml():
    """Test validation of file_categories.yaml structure."""
    print("Test: validate_file_categories_yaml")

    from agent.yaml_tools import validate_yaml_structure, load_yaml_file, find_yaml_files

    agent_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "agent")
    yaml_files = find_yaml_files(agent_dir)

    categories_file = [f for f in yaml_files if "file_categories.yaml" in f][0]
    data, _ = load_yaml_file(categories_file)

    issues = validate_yaml_structure(data, categories_file)

    # Should have no errors
    errors = [i for i in issues if i[0] == 'error']
    assert len(errors) == 0, f"Should have no errors: {errors}"

    print("  PASSED")


def test_validate_all_files():
    """Test that validate_files succeeds on all YAML files."""
    print("Test: validate_all_files")

    from agent.yaml_tools import validate_files

    # validate_files returns True if all pass
    success = validate_files([])  # Empty list uses default directory

    assert success, "All YAML files should validate successfully"

    print("  PASSED")


def test_extract_terms():
    """Test extraction of terms from YAML files."""
    print("Test: extract_terms")

    from agent.yaml_tools import collect_all_terms, find_yaml_files

    agent_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "agent")
    yaml_files = find_yaml_files(agent_dir)

    terms = collect_all_terms(yaml_files)

    # Should have various term categories
    assert "file_categories" in terms, "Should have file_categories"
    assert "strategy_keys" in terms, "Should have strategy_keys"
    assert "metrics" in terms, "Should have metrics"
    assert "workflow_phases" in terms, "Should have workflow_phases"

    # Check some expected values
    assert len(terms["file_categories"]) > 0, "Should have file categories"
    assert len(terms["strategy_keys"]) > 0, "Should have strategy keys"

    print("  PASSED")


def run_all_tests():
    """Run all yaml_tools tests."""
    print("=" * 60)
    print("YAML TOOLS TESTS")
    print("=" * 60)

    test_find_yaml_files()
    test_load_yaml_file()
    test_load_yaml_file_invalid()
    test_validate_programs_yaml()
    test_validate_workflows_yaml()
    test_validate_metrics_yaml()
    test_validate_file_categories_yaml()
    test_validate_all_files()
    test_extract_terms()

    print()
    print("=" * 60)
    print("ALL YAML TOOLS TESTS PASSED!")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
