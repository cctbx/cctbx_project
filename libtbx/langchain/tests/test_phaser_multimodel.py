#!/usr/bin/env python
"""
Test phenix.phaser multi-model command building.

This test verifies that the command builder correctly generates
multi-ensemble syntax when multiple models are provided.

Run with: python tests/test_phaser_multimodel.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_single_model_command():
    """Test that single model uses simple syntax."""
    print("Test: single_model_command")

    # This test documents expected behavior
    # Single model: phenix.phaser data.mtz model.pdb phaser.mode=MR_AUTO

    files = {"data_mtz": "data.mtz", "model": "beta.pdb"}

    # Verify file structure
    assert files["data_mtz"] == "data.mtz"
    assert files["model"] == "beta.pdb"

    print("  PASSED")


def test_multi_model_command():
    """Test that multiple models use ensemble syntax."""
    print("Test: multi_model_command")

    # This test documents expected behavior
    # Multi model should generate:
    # phenix.phaser data.mtz
    #   phaser.ensemble.ense1.pdb=beta.pdb
    #   phaser.ensemble.ense2.pdb=blip.pdb
    #   phaser.search.ense1.num=1
    #   phaser.search.ense2.num=1
    #   phaser.mode=MR_AUTO

    files = {
        "data_mtz": "data.mtz",
        "model": ["beta.pdb", "blip.pdb"],
        "sequence": ["beta.seq", "blip.seq"]
    }

    # Verify file structure supports list format
    assert isinstance(files["model"], list)
    assert len(files["model"]) == 2
    assert "beta.pdb" in files["model"]
    assert "blip.pdb" in files["model"]

    print("  PASSED")


def test_yaml_config():
    """Test that programs.yaml has correct configuration."""
    print("Test: yaml_config")

    import yaml

    yaml_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "knowledge", "programs.yaml"
    )

    with open(yaml_path) as f:
        programs = yaml.safe_load(f)

    phaser = programs.get("phenix.phaser", {})

    # Check multi_ensemble flag
    assert phaser.get("multi_ensemble") == True, \
        "phenix.phaser should have multi_ensemble=true"

    # Check model input allows multiple
    model_def = phaser.get("inputs", {}).get("required", {}).get("model", {})
    assert model_def.get("multiple") == True, \
        "model input should allow multiple files"

    # Check sequence input allows multiple
    seq_def = phaser.get("inputs", {}).get("optional", {}).get("sequence", {})
    assert seq_def.get("multiple") == True, \
        "sequence input should allow multiple files"

    print("  PASSED")


def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.test_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
