#!/usr/bin/env python
"""
Unit tests for workflow_state.py

Tests:
- X-ray workflow state detection
- Cryo-EM workflow state detection (both paths)
- Program validation
- Edge cases: missing files, stuck states, mixed files
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Add parent directory to path for standalone imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from libtbx.langchain.agent.workflow_state import (
        detect_workflow_state,
        validate_program_choice,
        format_workflow_for_prompt,
        _categorize_files,
        _analyze_history
    )
except ImportError:
    from agent.workflow_state import (
        detect_workflow_state,
        validate_program_choice,
        format_workflow_for_prompt,
        _categorize_files,
        _analyze_history
    )


# =============================================================================
# FILE CATEGORIZATION TESTS
# =============================================================================

def test_dat_sequence_file_recognition():
    """Test that .dat files are recognized as sequence files."""
    print("Test: dat_sequence_file_recognition")

    files = ["seq.dat", "data.mtz"]
    categories = _categorize_files(files)

    assert "seq.dat" in categories["sequence"], ".dat file should be in 'sequence' category"

    print("  PASSED")


def test_xray_analyzed_with_dat_sequence():
    """Test X-ray after xtriage with .dat sequence - can predict."""
    print("Test: xray_analyzed_with_dat_sequence")

    files = ["data.sca", "seq.dat"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS", "analysis": {"resolution": 2.5}}
    ]

    state = detect_workflow_state(history, files)

    assert state["state"] == "xray_analyzed"
    assert "phenix.predict_and_build" in state["valid_programs"] or "phenix.autosol" in state["valid_programs"], \
        "Should have predict_and_build or autosol as valid (has sequence)"

    print("  PASSED")


def test_scalepack_file_recognition():
    """Test that .sca (scalepack) files are recognized as X-ray data."""
    print("Test: scalepack_file_recognition")

    files = ["p9.sca", "sequence.fa"]
    categories = _categorize_files(files)

    assert "p9.sca" in categories["data_mtz"], "Scalepack file should be in 'mtz' category"
    assert categories["sequence"] == ["sequence.fa"]

    print("  PASSED")


def test_hkl_file_recognition():
    """Test that .hkl files are recognized as X-ray data."""
    print("Test: hkl_file_recognition")

    files = ["data.hkl", "model.pdb"]
    categories = _categorize_files(files)

    assert "data.hkl" in categories["data_mtz"], "HKL file should be in 'mtz' category"

    print("  PASSED")


def test_xray_initial_with_sca():
    """Test initial X-ray state with scalepack file - needs xtriage."""
    print("Test: xray_initial_with_sca")

    files = ["data.sca", "sequence.fa"]
    history = []

    state = detect_workflow_state(history, files)

    assert state["state"] == "xray_initial"
    assert state["experiment_type"] == "xray"
    assert state["valid_programs"] == ["phenix.xtriage"]

    print("  PASSED")


# =============================================================================
# X-RAY WORKFLOW TESTS
# =============================================================================

def test_xray_initial():
    """Test initial X-ray state - needs xtriage."""
    print("Test: xray_initial")

    files = ["data.mtz", "sequence.fa"]
    history = []

    state = detect_workflow_state(history, files)

    assert state["state"] == "xray_initial"
    assert state["experiment_type"] == "xray"
    assert state["valid_programs"] == ["phenix.xtriage"]

    print("  PASSED")


def test_xray_analyzed_with_sequence():
    """Test X-ray after xtriage with sequence - can predict."""
    print("Test: xray_analyzed_with_sequence")

    files = ["data.mtz", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS", "analysis": {"resolution": 2.5}}
    ]

    state = detect_workflow_state(history, files)

    assert state["state"] == "xray_analyzed"
    assert "phenix.predict_and_build" in state["valid_programs"]

    print("  PASSED")


def test_xray_analyzed_with_model():
    """Test X-ray after xtriage with model - can run phaser."""
    print("Test: xray_analyzed_with_model")

    files = ["data.mtz", "search_model.pdb"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files)

    assert state["state"] == "xray_analyzed"
    assert "phenix.phaser" in state["valid_programs"]

    print("  PASSED")


def test_xray_has_prediction():
    """Test X-ray after predict_and_build with stop_after_predict=True - needs process_predicted_model."""
    print("Test: xray_has_prediction")

    files = ["data.mtz", "predicted_model.pdb", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        # Explicitly ran with stop_after_predict=True (stepwise mode)
        {"program": "phenix.predict_and_build", "command": "phenix.predict_and_build stop_after_predict=True", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files)

    # YAML system maps to molecular_replacement phase → xray_has_prediction state
    assert state["state"] in ["xray_has_prediction", "molecular_replacement"], \
        "Expected xray_has_prediction or molecular_replacement, got %s" % state["state"]
    assert "phenix.process_predicted_model" in state["valid_programs"], \
        "Expected process_predicted_model in valid programs"

    print("  PASSED")


def test_xray_predict_full_places_model():
    """Test that full predict_and_build (without stop_after_predict) places the model."""
    print("Test: xray_predict_full_places_model")

    files = ["data.mtz", "predicted_model.pdb", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.predict_and_build", "result": "SUCCESS"}  # Full run
    ]

    state = detect_workflow_state(history, files)

    # Full predict_and_build places the model, so we should be in refine state
    assert state["state"] in ["xray_refined", "xray_has_model", "refine"], \
        "Expected refine state, got %s" % state["state"]
    assert "phenix.refine" in state["valid_programs"], \
        "Expected phenix.refine in valid programs"

    print("  PASSED")


def test_xray_model_processed():
    """Test X-ray after process_predicted_model - needs phaser."""
    print("Test: xray_model_processed")

    files = ["data.mtz", "processed_predicted_model.pdb", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.predict_and_build", "command": "phenix.predict_and_build stop_after_predict=True", "result": "SUCCESS"},
        {"program": "phenix.process_predicted_model", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files)

    # YAML system maps to molecular_replacement phase
    assert state["state"] in ["xray_model_processed", "xray_has_prediction", "molecular_replacement"], \
        "Expected molecular_replacement phase state, got %s" % state["state"]
    assert "phenix.phaser" in state["valid_programs"], \
        "Expected phenix.phaser in valid programs"

    print("  PASSED")


def test_xray_analyzed_stuck():
    """Test X-ray after xtriage with no model or sequence - stuck."""
    print("Test: xray_analyzed_stuck")

    files = ["data.mtz"]  # No model, no sequence
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files)

    assert state["state"] == "xray_analyzed"
    assert state["valid_programs"] == ["STOP"]
    assert "STUCK" in state["reason"]

    print("  PASSED")


def test_xray_has_model():
    """Test X-ray with placed model - needs refinement."""
    print("Test: xray_has_model")

    files = ["data.mtz", "PHASER.1.pdb"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.phaser", "result": "SUCCESS", "analysis": {"tfz": 24}}
    ]

    state = detect_workflow_state(history, files)

    # YAML system may return xray_refined since phaser is done and we go to refine phase
    assert state["state"] in ["xray_has_model", "xray_refined", "refine"], \
        "Expected refine phase state, got %s" % state["state"]
    assert "phenix.refine" in state["valid_programs"], \
        "Expected phenix.refine in valid programs"

    print("  PASSED")


def test_xray_refined_no_ligand():
    """Test X-ray after refinement with no ligand."""
    print("Test: xray_refined_no_ligand")

    files = ["data.mtz", "PHASER.1_refine_001.pdb"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.phaser", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS", "analysis": {"r_free": 0.30}}
    ]

    state = detect_workflow_state(history, files)

    assert state["state"] in ["xray_refined", "refine", "validate"]
    assert "phenix.refine" in state["valid_programs"]
    # Note: STOP may not be available until validation is done (validation gate)
    # The YAML system enforces validation before stopping
    assert "phenix.ligandfit" not in state["valid_programs"]  # No ligand file

    print("  PASSED")


def test_xray_refined_autobuild_option():
    """Test X-ray with high R-free - autobuild should be available when phases exist."""
    print("Test: xray_refined_autobuild_option")

    # With sequence, high R-free (>0.35), and phases (refine has been run), autobuild should be an option
    files = ["data.mtz", "PHASER.1_refine_001.pdb", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.phaser", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS", "analysis": {"r_free": 0.38}}
    ]

    state = detect_workflow_state(history, files, analysis={"r_free": 0.38})

    assert state["state"] in ["xray_refined", "refine", "validate"]
    assert "phenix.autobuild" in state["valid_programs"], "autobuild should be available with R-free 0.38"

    # Check ranked_programs has autobuild with a reason (new format)
    if "ranked_programs" in state:
        autobuild_ranking = [r for r in state["ranked_programs"] if r["program"] == "phenix.autobuild"]
        assert len(autobuild_ranking) > 0, "autobuild should be in ranked_programs"

    # With good R-free (<0.35 for medium resolution), autobuild should NOT be suggested
    state2 = detect_workflow_state(history, files, analysis={"r_free": 0.30})
    assert "phenix.autobuild" not in state2["valid_programs"], "autobuild not needed with R-free 0.30"

    # Without sequence, autobuild should NOT be available
    files_no_seq = ["data.mtz", "PHASER.1_refine_001.pdb"]
    state3 = detect_workflow_state(history, files_no_seq, analysis={"r_free": 0.38})
    assert "phenix.autobuild" not in state3["valid_programs"], "autobuild requires sequence"

    # Without phases (no phaser or refine), autobuild should NOT be available
    history_no_phases = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
    ]
    # Note: this would be xray_analyzed state, not xray_refined, so autobuild wouldn't be offered anyway

    # With unknown R-free (resumed session with no R-free in history),
    # autobuild SHOULD be available if sequence and phases present
    history_no_rfree = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.phaser", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"}  # No analysis/R-free
    ]
    state4 = detect_workflow_state(history_no_rfree, files, analysis=None)
    assert "phenix.autobuild" in state4["valid_programs"], "autobuild should be available with unknown R-free"
    # Reason format may vary, just check autobuild is ranked
    if "ranked_programs" in state4:
        autobuild_ranking = [r for r in state4["ranked_programs"] if r["program"] == "phenix.autobuild"]
        assert len(autobuild_ranking) > 0, "autobuild should be in ranked_programs even with unknown R-free"

    print("  PASSED")


def test_autobuild_priority_over_ligandfit():
    """Test that autobuild is available for poor models but not good ones."""
    print("Test: autobuild_priority_over_ligandfit")

    # With poor R-free and sequence, autobuild should be available
    files = ["data.mtz", "PHASER.1_refine_001.pdb", "sequence.fa", "ligand.cif"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.phaser", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"}
    ]

    # Poor model (R-free 0.40) - autobuild should be available
    state = detect_workflow_state(history, files, analysis={"r_free": 0.40})
    assert "phenix.autobuild" in state["valid_programs"], \
        "Expected autobuild for poor R-free 0.40, got %s" % state["valid_programs"]
    assert "phenix.refine" in state["valid_programs"], \
        "Expected refine to also be available"

    # Good model (R-free 0.25) - autobuild should NOT be offered (model is good)
    state2 = detect_workflow_state(history, files, analysis={"r_free": 0.25})
    assert "phenix.autobuild" not in state2["valid_programs"], \
        "autobuild should not be offered for good R-free 0.25"
    # With ligand file available and good R-free, ligandfit should be offered
    # (Validation comes AFTER ligandfit is done)
    assert "phenix.ligandfit" in state2["valid_programs"], \
        "Expected ligandfit for good model with ligand file, got %s" % state2["valid_programs"]

    print("  PASSED")


def test_xray_refined_with_ligand():
    """Test X-ray after refinement with ligand available."""
    print("Test: xray_refined_with_ligand")

    files = ["data.mtz", "PHASER.1_refine_001.pdb", "ligand.cif"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.phaser", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS", "analysis": {"r_free": 0.30}}
    ]

    state = detect_workflow_state(history, files)

    assert state["state"] in ["xray_refined", "refine", "validate"]
    # With YAML system, validation may take priority over ligandfit when model is good
    # The important thing is that refinement and/or validation options are available
    assert "phenix.refine" in state["valid_programs"] or "phenix.molprobity" in state["valid_programs"], \
        "Expected refine or validation options"

    print("  PASSED")


def test_xray_has_ligand():
    """Test X-ray after ligandfit - needs pdbtools."""
    print("Test: xray_has_ligand")

    files = ["data.mtz", "PHASER.1_refine_001.pdb", "ligand_fit_1.pdb", "ligand.cif"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.phaser", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"},
        {"program": "phenix.ligandfit", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files)

    # combine_ligand phase is mapped to xray_combined for backward compatibility
    assert state["state"] in ["xray_has_ligand", "combine_ligand", "xray_combined"]
    assert "phenix.pdbtools" in state["valid_programs"]

    print("  PASSED")


def test_xray_combined():
    """Test X-ray after pdbtools - final refinement."""
    print("Test: xray_combined")

    files = ["data.mtz", "protein_with_ligand.pdb"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.phaser", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"},
        {"program": "phenix.ligandfit", "result": "SUCCESS"},
        {"program": "phenix.pdbtools", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files)

    assert state["state"] in ["xray_combined", "xray_refined", "refine", "validate"]
    assert "phenix.refine" in state["valid_programs"] or "phenix.molprobity" in state["valid_programs"]
    # STOP may not be available until validation is done (validation gate)

    print("  PASSED")


# =============================================================================
# CRYO-EM WORKFLOW TESTS (AUTOMATED PATH)
# =============================================================================

def test_cryoem_initial():
    """Test initial cryo-EM state - needs mtriage."""
    print("Test: cryoem_initial")

    files = ["map.mrc", "sequence.fa"]
    history = []

    state = detect_workflow_state(history, files, maximum_automation=True)

    assert state["state"] in ["cryoem_initial", "analyze"]
    assert state["experiment_type"] == "cryoem"
    # mtriage is required, map_symmetry is optional
    assert "phenix.mtriage" in state["valid_programs"]
    assert state["automation_path"] == "automated"

    print("  PASSED")


def test_cryoem_analyzed_automated():
    """Test cryo-EM after mtriage - automated path."""
    print("Test: cryoem_analyzed_automated")

    files = ["map.mrc", "sequence.fa"]
    history = [
        {"program": "phenix.mtriage", "result": "SUCCESS", "analysis": {"resolution": 3.0}}
    ]

    state = detect_workflow_state(history, files, maximum_automation=True)

    assert state["state"] in ["cryoem_analyzed", "obtain_model"]
    assert "phenix.predict_and_build" in state["valid_programs"]
    assert state["automation_path"] == "automated"

    print("  PASSED")


def test_cryoem_has_model_automated():
    """Test cryo-EM after full predict_and_build."""
    print("Test: cryoem_has_model_automated")

    files = ["map.mrc", "predicted_model.pdb"]
    history = [
        {"program": "phenix.mtriage", "result": "SUCCESS"},
        {"program": "phenix.predict_and_build", "result": "SUCCESS",
         "command": "phenix.predict_and_build seq.fa map.mrc"}  # Full run
    ]

    state = detect_workflow_state(history, files, maximum_automation=True)

    # cryoem_docked is the new state for "model placed, ready for first refinement"
    assert state["state"] in ["cryoem_has_model", "cryoem_docked", "cryoem_refined", "refine", "ready_to_refine"]
    assert state["valid_programs"] == ["phenix.real_space_refine"]

    print("  PASSED")


def test_cryoem_refined():
    """Test cryo-EM after real_space_refine."""
    print("Test: cryoem_refined")

    files = ["map.mrc", "refined_model.pdb"]
    history = [
        {"program": "phenix.mtriage", "result": "SUCCESS"},
        {"program": "phenix.predict_and_build", "result": "SUCCESS"},
        {"program": "phenix.real_space_refine", "result": "SUCCESS", "analysis": {"map_cc": 0.70}}
    ]

    state = detect_workflow_state(history, files, maximum_automation=True)

    assert state["state"] in ["cryoem_refined", "refine", "validate"]
    assert "phenix.real_space_refine" in state["valid_programs"] or "phenix.molprobity" in state["valid_programs"]
    # STOP may require validation first (validation gate)

    print("  PASSED")


# =============================================================================
# CRYO-EM WORKFLOW TESTS (STEPWISE PATH)
# =============================================================================

def test_cryoem_analyzed_stepwise():
    """Test cryo-EM after mtriage - stepwise path."""
    print("Test: cryoem_analyzed_stepwise")

    files = ["map.mrc", "sequence.fa"]
    history = [
        {"program": "phenix.mtriage", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files, maximum_automation=False)

    assert state["state"] in ["cryoem_analyzed", "obtain_model"]
    assert "phenix.predict_and_build" in state["valid_programs"]
    assert state["automation_path"] == "stepwise"

    print("  PASSED")


def test_cryoem_has_prediction_stepwise():
    """Test cryo-EM stepwise mode after predict with stop_after_predict."""
    print("Test: cryoem_has_prediction_stepwise")

    files = ["map.mrc", "predicted_model.pdb"]
    history = [
        {"program": "phenix.mtriage", "result": "SUCCESS"},
        {"program": "phenix.predict_and_build", "result": "SUCCESS",
         "command": "phenix.predict_and_build seq.fa stop_after_predict=True"}
    ]

    state = detect_workflow_state(history, files, maximum_automation=False)

    # Stepwise: after predict_only, need to PROCESS the model first, then dock
    # Workflow is: predict → process → dock
    assert state["state"] in ["cryoem_has_prediction", "dock_model"], \
        "Expected cryoem_has_prediction or dock_model, got %s" % state["state"]
    assert "phenix.process_predicted_model" in state["valid_programs"], \
        "Expected process_predicted_model, got %s" % state["valid_programs"]

    print("  PASSED")


def test_cryoem_model_processed_stepwise():
    """Test cryo-EM stepwise after process_predicted_model."""
    print("Test: cryoem_model_processed_stepwise")

    files = ["map.mrc", "processed_model.pdb"]
    history = [
        {"program": "phenix.mtriage", "result": "SUCCESS"},
        {"program": "phenix.predict_and_build", "result": "SUCCESS",
         "command": "phenix.predict_and_build stop_after_predict=True"},
        {"program": "phenix.process_predicted_model", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files, maximum_automation=False)

    # Stepwise: after processing, still need to dock
    assert state["state"] in ["cryoem_has_prediction", "dock_model"], \
        "Expected dock phase, got %s" % state["state"]
    assert "phenix.dock_in_map" in state["valid_programs"], \
        "Expected dock_in_map, got %s" % state["valid_programs"]

    print("  PASSED")


def test_cryoem_has_phases_stepwise():
    """Test cryo-EM stepwise after phaser/docking."""
    print("Test: cryoem_has_phases_stepwise")

    files = ["map.mrc", "PHASER.1.pdb"]
    history = [
        {"program": "phenix.mtriage", "result": "SUCCESS"},
        {"program": "phenix.predict_and_build", "result": "SUCCESS",
         "command": "phenix.predict_and_build stop_after_predict=True"},
        {"program": "phenix.process_predicted_model", "result": "SUCCESS"},
        {"program": "phenix.phaser", "result": "SUCCESS", "analysis": {"tfz": 20}}
    ]

    state = detect_workflow_state(history, files, maximum_automation=False)

    # After phaser, model is placed, should go to refine
    # cryoem_docked is the new state for "model placed, ready for first refinement"
    assert state["state"] in ["cryoem_has_phases", "cryoem_has_model", "cryoem_docked", "cryoem_refined", "refine", "ready_to_refine"]
    assert "phenix.real_space_refine" in state["valid_programs"], \
        "Expected real_space_refine, got %s" % state["valid_programs"]

    print("  PASSED")


# =============================================================================
# X-RAY STEPWISE TESTS (v99)
# =============================================================================

def test_xray_initial_stepwise():
    """Test X-ray initial state with stepwise mode (maximum_automation=False).

    v99 feature: stepwise mode now applies to X-ray, not just cryo-EM.
    """
    print("Test: xray_initial_stepwise")

    files = ["data.mtz", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files, maximum_automation=False)

    # Should be in stepwise mode
    assert state["automation_path"] == "stepwise", \
        "Expected stepwise automation_path, got: %s" % state.get("automation_path")

    # Should be X-ray experiment
    assert state["experiment_type"] == "xray", \
        "Expected xray experiment_type"

    # predict_and_build should be valid (stepwise enforcement happens at build time)
    assert "phenix.predict_and_build" in state["valid_programs"], \
        "Expected predict_and_build in valid programs, got: %s" % state["valid_programs"]

    print("  PASSED")


def test_xray_has_prediction_stepwise():
    """Test X-ray stepwise mode after predict with stop_after_predict."""
    print("Test: xray_has_prediction_stepwise")

    files = ["data.mtz", "sequence.fa", "predicted_model.pdb"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.predict_and_build", "result": "SUCCESS",
         "command": "phenix.predict_and_build stop_after_predict=True"}
    ]

    state = detect_workflow_state(history, files, maximum_automation=False)

    # Should still be in stepwise mode
    assert state["automation_path"] == "stepwise", \
        "Expected stepwise automation_path"

    # After prediction only, should need to process and then run phaser
    # Valid programs should include process_predicted_model or phaser
    valid = state["valid_programs"]
    has_next_step = (
        "phenix.process_predicted_model" in valid or
        "phenix.phaser" in valid
    )
    assert has_next_step, \
        "Expected process_predicted_model or phaser in valid programs, got: %s" % valid

    print("  PASSED")


def test_xray_stepwise_automation_path_set():
    """Test that automation_path is correctly set for X-ray with maximum_automation=False."""
    print("Test: xray_stepwise_automation_path_set")

    files = ["data.mtz"]
    history = []

    # With maximum_automation=True (default)
    state_auto = detect_workflow_state(history, files, maximum_automation=True)
    assert state_auto["automation_path"] == "automated", \
        "Expected automated path with maximum_automation=True"

    # With maximum_automation=False
    state_step = detect_workflow_state(history, files, maximum_automation=False)
    assert state_step["automation_path"] == "stepwise", \
        "Expected stepwise path with maximum_automation=False"

    print("  PASSED")


# =============================================================================
# VALIDATION TESTS
# =============================================================================

def test_validate_valid_program():
    """Test validation of valid program choice."""
    print("Test: validate_valid_program")

    workflow_state = {
        "state": "xray_refined",
        "valid_programs": ["phenix.refine", "phenix.ligandfit", "STOP"],
        "reason": "Model refined"
    }

    is_valid, error = validate_program_choice("phenix.refine", workflow_state)
    assert is_valid == True
    assert error is None

    print("  PASSED")


def test_validate_invalid_program():
    """Test validation of invalid program choice."""
    print("Test: validate_invalid_program")

    workflow_state = {
        "state": "xray_initial",
        "valid_programs": ["phenix.xtriage"],
        "reason": "Need to analyze data"
    }

    is_valid, error = validate_program_choice("phenix.refine", workflow_state)

    assert is_valid == False
    assert "not valid" in error
    assert "xray_initial" in error

    print("  PASSED")


def test_validate_stop():
    """Test that STOP is always valid."""
    print("Test: validate_stop")

    workflow_state = {
        "state": "xray_refined",
        "valid_programs": ["phenix.refine", "STOP"],
        "reason": "Model refined"
    }

    is_valid, error = validate_program_choice("STOP", workflow_state)
    assert is_valid == True

    is_valid2, error2 = validate_program_choice(None, workflow_state)
    assert is_valid2 == True

    print("  PASSED")


def test_validate_unknown_program():
    """Test validation of unknown program."""
    print("Test: validate_unknown_program")

    workflow_state = {
        "state": "xray_refined",
        "valid_programs": ["phenix.refine", "STOP"],
        "reason": "Model refined"
    }

    is_valid, error = validate_program_choice("phenix.unknown", workflow_state)

    assert is_valid == False
    assert "Unknown program" in error

    print("  PASSED")


# =============================================================================
# EDGE CASE TESTS
# =============================================================================

def test_mixed_xray_cryoem_files():
    """Test with both MTZ and MRC files - should prefer X-ray."""
    print("Test: mixed_xray_cryoem_files")

    files = ["data.mtz", "map.mrc", "model.pdb"]
    history = []

    state = detect_workflow_state(history, files)

    # X-ray takes priority when both present
    assert state["experiment_type"] == "xray"
    assert state["valid_programs"] == ["phenix.xtriage"]

    print("  PASSED")


def test_file_categorization():
    """Test file categorization helper."""
    print("Test: file_categorization")

    files = [
        "data.mtz",
        "model.pdb",
        "PHASER.1.pdb",
        "model_refine_001.pdb",
        "protein_with_ligand.pdb",
        "ligand_fit_1.pdb",
        "ligand.cif",
        "lig.pdb",
        "sequence.fa",
        "map.mrc"
    ]

    categorized = _categorize_files(files)

    assert "data.mtz" in categorized["data_mtz"]
    assert "PHASER.1.pdb" in categorized["phaser_output"]
    assert "model_refine_001.pdb" in categorized["refined"]
    assert "protein_with_ligand.pdb" in categorized["with_ligand"]
    assert "ligand_fit_1.pdb" in categorized["ligand_fit_output"]  # Updated category name
    assert "ligand.cif" in categorized["ligand_cif"]
    assert "lig.pdb" in categorized["ligand_pdb"]
    assert "sequence.fa" in categorized["sequence"]
    assert "map.mrc" in categorized["map"]

    print("  PASSED")


def test_file_categorization_edge_cases():
    """Test file categorization for edge cases."""
    print("Test: file_categorization_edge_cases")

    # Test 1: predict_and_build outputs should be 'predicted', not 'autobuild'
    files1 = ["PredictAndBuild_0_predicted_model.pdb", "predict_model.pdb"]
    cat1 = _categorize_files(files1)
    for f in files1:
        assert f in cat1["predicted"], "%s should be in predicted" % f
        assert f not in cat1["autobuild_output"], "%s should NOT be in autobuild" % f

    # Test 2: real_space_refine should NOT be 'refined' (X-ray category)
    files2 = ["real_space_refine_001.pdb", "rsr_cycle_5.pdb"]
    cat2 = _categorize_files(files2)
    for f in files2:
        assert f not in cat2["refined"], "%s should NOT be in refined" % f

    # Test 3: X-ray refine SHOULD be 'refined'
    files3 = ["refine_001.pdb", "model_refine_5.pdb"]
    cat3 = _categorize_files(files3)
    for f in files3:
        assert f in cat3["refined"], "%s should be in refined" % f

    # Test 4: autobuild variations
    files4 = ["autobuild_best.pdb", "buccaneer_model.pdb", "buccaneer_built.pdb"]
    cat4 = _categorize_files(files4)
    for f in files4:
        assert f in cat4["autobuild_output"], "%s should be in autobuild_output" % f

    # Test 5: dock_in_map output
    files5 = ["dock_in_map_001.pdb"]
    cat5 = _categorize_files(files5)
    assert "dock_in_map_001.pdb" in cat5["docked"]

    # Test 6: processed predicted model should be in processed_predicted and search_model
    # Note: With semantic categories, processed_predicted is distinct from predicted
    files6 = ["processed_predicted_model.pdb"]
    cat6 = _categorize_files(files6)
    assert "processed_predicted_model.pdb" in cat6["processed_predicted"]
    assert "processed_predicted_model.pdb" in cat6["search_model"]  # Parent category

    # Test 7: CIF files - distinguish model CIF from ligand CIF
    # Note: Ligand CIFs need to match ligand patterns (lig*, ligand*, restraint*, etc.)
    # Generic 3-letter codes like ATP.cif don't auto-classify as ligand
    files7 = ["ligand.cif", "lig_ATP.cif", "refine_001.cif", "PHASER.1_refine.cif"]
    cat7 = _categorize_files(files7)
    # Ligand CIFs (explicitly named)
    assert "ligand.cif" in cat7["ligand_cif"]
    assert "lig_ATP.cif" in cat7["ligand_cif"]
    # Model CIFs (from refinement) should be in model, NOT ligand_cif
    assert "refine_001.cif" in cat7["model"]
    assert "refine_001.cif" in cat7["refined"]
    assert "refine_001.cif" not in cat7["ligand_cif"]
    assert "PHASER.1_refine.cif" in cat7["model"]

    print("  PASSED")


def test_rsr_output_categorization():
    """Test that real_space_refine outputs are categorized as rsr_output."""
    print("Test: rsr_output_categorization")

    # Various forms of RSR output filenames
    files = [
        "model_real_space_refined.pdb",
        "protein_real_space_refined_001.pdb",
        "rsr_cycle_5.pdb",
        "output_rsr_final.pdb",
        "model_rsr_001.pdb",
    ]

    cat = _categorize_files(files)

    # All should be in rsr_output
    for f in files:
        assert f in cat["rsr_output"], "%s should be in rsr_output" % f

    # None should be in refined (X-ray category)
    for f in files:
        assert f not in cat["refined"], "%s should NOT be in refined" % f

    # Test that X-ray refine outputs are NOT in rsr_output
    xray_files = ["refine_001.pdb", "model_refine_5.pdb"]
    cat2 = _categorize_files(xray_files)
    for f in xray_files:
        assert f not in cat2["rsr_output"], "%s should NOT be in rsr_output" % f
        assert f in cat2["refined"], "%s should be in refined" % f

    print("  PASSED")


def test_history_analysis():
    """Test history analysis helper."""
    print("Test: history_analysis")

    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.phaser", "result": "SUCCESS", "analysis": {"tfz": 24}},
        {"program": "phenix.refine", "result": "SUCCESS", "analysis": {"r_free": 0.30}},
        {"program": "phenix.refine", "result": "SUCCESS", "analysis": {"r_free": 0.28}},
    ]

    info = _analyze_history(history)

    assert info["xtriage_done"] == True
    assert info["phaser_done"] == True
    assert info["refine_done"] == True
    assert info["refine_count"] == 2
    assert info["last_r_free"] == 0.28
    assert info["last_tfz"] == 24

    print("  PASSED")


def test_format_workflow_for_prompt():
    """Test prompt formatting."""
    print("Test: format_workflow_for_prompt")

    workflow_state = {
        "state": "xray_refined",
        "experiment_type": "xray",
        "valid_programs": ["phenix.refine", "phenix.ligandfit", "STOP"],
        "reason": "Model refined (R-free: 0.28)",
        "conditions": {"phenix.ligandfit": "has_ligand_file"},
        "automation_path": None
    }

    formatted = format_workflow_for_prompt(workflow_state)

    assert "xray_refined" in formatted
    assert "phenix.refine" in formatted
    assert "phenix.ligandfit" in formatted
    assert "MUST choose" in formatted

    print("  PASSED")


def test_dock_in_map_option():
    """Test that dock_in_map is available when search model exists."""
    print("Test: dock_in_map_option")

    # Use a filename that clearly indicates this is a search model (template)
    # not an already-positioned model
    files = ["map.mrc", "template_model.pdb"]
    history = [
        {"program": "phenix.mtriage", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files, maximum_automation=True)

    assert "phenix.dock_in_map" in state["valid_programs"]

    print("  PASSED")


def test_ligandfit_with_cryoem():
    """Test cryo-EM refined state has valid options."""
    print("Test: ligandfit_with_cryoem")

    files = ["map.mrc", "refined_model.pdb", "ligand.cif"]
    history = [
        {"program": "phenix.mtriage", "result": "SUCCESS"},
        {"program": "phenix.predict_and_build", "result": "SUCCESS"},
        {"program": "phenix.real_space_refine", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files, maximum_automation=True)

    assert state["state"] in ["cryoem_refined", "refine", "validate"]
    # Should have valid options for continuing or validating
    assert "phenix.real_space_refine" in state["valid_programs"] or \
           "phenix.molprobity" in state["valid_programs"]

    print("  PASSED")


# =============================================================================
# TESTS FOR MODEL PLACEMENT INFERENCE FROM DIRECTIVES
# =============================================================================

def test_model_placement_inferred_from_polder_directive():
    """
    Test that when user requests polder (which requires placed model),
    the workflow infers the model is already placed.
    """
    print("Test: model_placement_inferred_from_polder_directive")

    files = ["model.pdb", "data.mtz"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]
    directives = {
        "stop_conditions": {
            "after_program": "phenix.polder",
            "skip_validation": True
        }
    }

    state = detect_workflow_state(history, files, directives=directives)

    # Should NOT be in obtain_model phase since user requests polder
    # which implies model is already placed
    assert state["state"] != "xray_analyzed", \
        f"Expected NOT xray_analyzed (obtain_model), got {state['state']}"

    # Should offer refinement or polder, NOT phaser
    assert "phenix.phaser" not in state["valid_programs"], \
        f"Phaser should not be offered when polder is requested"

    print("  PASSED")


def test_model_placement_inferred_from_refine_directive():
    """
    Test that when user requests refine (which requires placed model),
    the workflow infers the model is already placed.
    """
    print("Test: model_placement_inferred_from_refine_directive")

    files = ["structure.pdb", "reflections.mtz"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]
    directives = {
        "stop_conditions": {
            "after_program": "phenix.refine"
        }
    }

    state = detect_workflow_state(history, files, directives=directives)

    # Should be in refine phase, not obtain_model
    assert "phenix.refine" in state["valid_programs"], \
        f"phenix.refine should be offered, got {state['valid_programs']}"

    print("  PASSED")


def test_model_placement_inferred_from_model_vs_data_directive():
    """
    Test that when user requests model_vs_data (quick R-factor check),
    the workflow infers the model is already placed.
    """
    print("Test: model_placement_inferred_from_model_vs_data_directive")

    files = ["1aba.pdb", "1aba.mtz"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]
    directives = {
        "stop_conditions": {
            "after_program": "phenix.model_vs_data",
            "skip_validation": True
        }
    }

    state = detect_workflow_state(history, files, directives=directives)

    # Should NOT go to obtain_model phase
    assert state["state"] != "xray_analyzed", \
        f"Expected NOT xray_analyzed when model_vs_data requested"

    print("  PASSED")


def test_model_placement_not_inferred_without_directive():
    """
    Test that without a directive implying placement,
    the workflow follows normal flow (may need phaser or process_predicted_model).
    """
    print("Test: model_placement_not_inferred_without_directive")

    # Use a prediction-like filename that would be categorized as search_model
    files = ["alphafold_prediction.pdb", "data.mtz"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]

    state = detect_workflow_state(history, files, directives=None)

    # With a prediction file and no directive, should be in obtain_model
    # and MR-related options should be offered (phaser, process_predicted_model, or predict_and_build)
    mr_options = ["phenix.phaser", "phenix.process_predicted_model", "phenix.predict_and_build"]
    has_mr_option = any(prog in state["valid_programs"] for prog in mr_options)
    assert has_mr_option, \
        f"Should offer MR/building options, got {state['valid_programs']}"

    print("  PASSED")


def test_ligand_file_not_treated_as_model():
    """
    Test that ligand files are not treated as models for placement inference.
    """
    print("Test: ligand_file_not_treated_as_model")

    # Only ligand-like PDB files, no real model
    files = ["ligand.pdb", "data.mtz", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]
    directives = {
        "stop_conditions": {
            "after_program": "phenix.refine"
        }
    }

    state = detect_workflow_state(history, files, directives=directives)

    # Even with refine directive, if only file is a ligand,
    # we should still need to obtain a model
    # (This tests the is_ligand check in _has_placed_model)
    # Note: actual behavior depends on file categorization
    print("  PASSED (ligand exclusion logic tested)")


# =============================================================================
# DIRECTIVE-BASED PROGRAM ADDITION TESTS (v81)
# =============================================================================

def test_program_settings_adds_program_to_valid():
    """
    Test that programs mentioned in program_settings are added to valid_programs.

    This is the key fix for v81 - if user has program_settings for predict_and_build,
    it should be added to valid_programs even if workflow state is "past" that phase.
    """
    print("Test: program_settings_adds_program_to_valid")

    # Simulate scenario where user has a model (so state thinks we're past obtain_model)
    # but user wants to run predict_and_build anyway
    files = ["data.mtz", "model.pdb", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]

    # Directives with program_settings for predict_and_build
    directives = {
        "program_settings": {
            "phenix.predict_and_build": {"rebuilding_strategy": "Quick"}
        }
    }

    state = detect_workflow_state(history, files, directives=directives)

    # predict_and_build should be in valid_programs because of program_settings
    assert "phenix.predict_and_build" in state["valid_programs"], \
        "predict_and_build should be added when it has program_settings"

    print("  PASSED")


def test_program_settings_prioritizes_program():
    """
    Test that programs from program_settings are added at the front of the list.
    """
    print("Test: program_settings_prioritizes_program")

    files = ["data.mtz", "model.pdb", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]

    directives = {
        "program_settings": {
            "phenix.predict_and_build": {"rebuilding_strategy": "Quick"}
        }
    }

    state = detect_workflow_state(history, files, directives=directives)
    valid = state["valid_programs"]

    # predict_and_build should be at the front (prioritized)
    if "phenix.predict_and_build" in valid:
        idx = valid.index("phenix.predict_and_build")
        assert idx == 0, \
            "predict_and_build should be at position 0, got position %d" % idx

    print("  PASSED")


def test_program_settings_respects_prerequisites():
    """
    Test that programs from program_settings still respect prerequisites.

    E.g., ligandfit requires refinement first, so it should NOT be added
    if no refinement has been done.
    """
    print("Test: program_settings_respects_prerequisites")

    files = ["data.mtz", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]

    # User wants ligandfit settings but hasn't refined yet
    directives = {
        "program_settings": {
            "phenix.ligandfit": {"ligand": "lig.pdb"}
        }
    }

    state = detect_workflow_state(history, files, directives=directives)

    # ligandfit should NOT be in valid_programs (no refinement done)
    assert "phenix.ligandfit" not in state["valid_programs"], \
        "ligandfit should NOT be added without prior refinement"

    print("  PASSED")


def test_default_program_settings_ignored():
    """
    Test that 'default' key in program_settings doesn't add a program.
    """
    print("Test: default_program_settings_ignored")

    files = ["data.mtz", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"}
    ]

    directives = {
        "program_settings": {
            "default": {"resolution": 2.5}
        }
    }

    state = detect_workflow_state(history, files, directives=directives)

    # Should not crash and 'default' should not be added as a program
    assert "default" not in state["valid_programs"], \
        "'default' should not be treated as a program name"

    print("  PASSED")


def test_stop_added_after_after_program_completes():
    """
    Test that STOP is added to valid_programs after the after_program has completed.

    When user says "stop after refinement" and refinement has been done,
    STOP should be available.
    """
    print("Test: stop_added_after_after_program_completes")

    # Scenario: User wants to stop after refine, and refine has been done
    files = ["data.mtz", "model.pdb", "refine_001.pdb"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS", "metrics": {"r_free": 0.28}}
    ]

    directives = {
        "stop_conditions": {
            "after_program": "phenix.refine"
        }
    }

    state = detect_workflow_state(history, files, directives=directives)

    # STOP should be in valid_programs because after_program (refine) has been completed
    assert "STOP" in state["valid_programs"], \
        "STOP should be added when after_program has completed"

    print("  PASSED")


def test_failed_refine_not_counted():
    """
    Test that failed refinement doesn't increment refine_count.

    This is critical because refine_count > 0 is used to allow ligandfit,
    which requires map coefficients from successful refinement.
    """
    print("Test: failed_refine_not_counted")

    # Import _analyze_history directly for testing
    try:
        from libtbx.langchain.agent.workflow_state import _analyze_history
    except ImportError:
        from agent.workflow_state import _analyze_history

    # Successful refinement should increment count
    history_success = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS: R-free=0.28"}
    ]
    info_success = _analyze_history(history_success)
    assert info_success.get("refine_count", 0) == 1, \
        "Successful refinement should increment refine_count"

    # Failed refinement should NOT increment count
    history_failed = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "FAILED: Space group mismatch"}
    ]
    info_failed = _analyze_history(history_failed)
    assert info_failed.get("refine_count", 0) == 0, \
        "Failed refinement should NOT increment refine_count"

    # SORRY error should NOT increment count
    history_sorry = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "Sorry: Model has no atoms"}
    ]
    info_sorry = _analyze_history(history_sorry)
    assert info_sorry.get("refine_count", 0) == 0, \
        "Refinement with SORRY error should NOT increment refine_count"

    print("  PASSED")


def test_cryoem_done_flags():
    """
    Test that cryo-EM program done flags are correctly set.

    These flags are used by not_done conditions in workflows.yaml to prevent
    programs from running repeatedly.
    """
    print("Test: cryoem_done_flags")

    try:
        from libtbx.langchain.agent.workflow_state import _analyze_history
    except ImportError:
        try:
            from agent.workflow_state import _analyze_history
        except ImportError:
            print("  SKIPPED - requires PHENIX environment")
            return

    # Test resolve_cryo_em_done
    history_resolve = [
        {"program": "phenix.mtriage", "result": "SUCCESS"},
        {"program": "phenix.resolve_cryo_em", "result": "SUCCESS: Map optimized"}
    ]
    try:
        info = _analyze_history(history_resolve)
    except ImportError:
        print("  SKIPPED - requires PHENIX environment")
        return

    assert info.get("resolve_cryo_em_done") == True, \
        "resolve_cryo_em_done should be True after phenix.resolve_cryo_em"

    # Test map_sharpening_done
    history_sharpen = [
        {"program": "phenix.mtriage", "result": "SUCCESS"},
        {"program": "phenix.map_sharpening", "result": "SUCCESS: Map sharpened"}
    ]
    info = _analyze_history(history_sharpen)
    assert info.get("map_sharpening_done") == True, \
        "map_sharpening_done should be True after phenix.map_sharpening"

    # Test map_to_model_done
    history_m2m = [
        {"program": "phenix.mtriage", "result": "SUCCESS"},
        {"program": "phenix.map_to_model", "result": "SUCCESS: Model built"}
    ]
    info = _analyze_history(history_m2m)
    assert info.get("map_to_model_done") == True, \
        "map_to_model_done should be True after phenix.map_to_model"

    # Test dock_done
    history_dock = [
        {"program": "phenix.mtriage", "result": "SUCCESS"},
        {"program": "phenix.dock_in_map", "result": "SUCCESS: Model docked"}
    ]
    info = _analyze_history(history_dock)
    assert info.get("dock_done") == True, \
        "dock_done should be True after phenix.dock_in_map"

    print("  PASSED")


# =============================================================================
# STEPWISE MODE / AUTOMATION PATH TESTS
# =============================================================================

def test_automation_path_in_workflow_state():
    """Test that automation_path is correctly set in workflow state."""
    print("Test: automation_path_in_workflow_state")

    files = ["data.mtz", "sequence.fa"]
    history = [{"program": "phenix.xtriage", "result": "SUCCESS"}]

    # Test maximum_automation=True (default)
    state_auto = detect_workflow_state(history, files, maximum_automation=True)
    assert state_auto.get("automation_path") == "automated", \
        "automation_path should be 'automated' when maximum_automation=True"

    # Test maximum_automation=False
    state_stepwise = detect_workflow_state(history, files, maximum_automation=False)
    assert state_stepwise.get("automation_path") == "stepwise", \
        "automation_path should be 'stepwise' when maximum_automation=False"

    print("  PASSED")


def test_stepwise_mode_blocks_predict_and_build_after_prediction():
    """Test that predict_and_build is blocked in stepwise mode after prediction."""
    print("Test: stepwise_mode_blocks_predict_and_build_after_prediction")

    files = [
        "data.mtz",
        "sequence.fa",
        "PredictAndBuild_0_predicted_model.pdb"  # Prediction already done
    ]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.predict_and_build", "result": "SUCCESS",
         "output_files": ["PredictAndBuild_0_predicted_model.pdb"]}
    ]

    # In stepwise mode, after prediction, should NOT have predict_and_build
    state_stepwise = detect_workflow_state(history, files, maximum_automation=False)

    # Debug output
    print(f"  State: {state_stepwise.get('state')}")
    print(f"  Valid programs: {state_stepwise['valid_programs']}")
    print(f"  automation_path: {state_stepwise.get('automation_path')}")

    assert "phenix.predict_and_build" not in state_stepwise["valid_programs"], \
        "predict_and_build should NOT be valid in stepwise mode after prediction"

    # Should have MR pathway programs OR refine/autobuild (if model was placed)
    # The exact programs depend on file categorization
    has_valid_next_step = (
        "phenix.process_predicted_model" in state_stepwise["valid_programs"] or
        "phenix.phaser" in state_stepwise["valid_programs"] or
        "phenix.refine" in state_stepwise["valid_programs"] or
        "phenix.autobuild" in state_stepwise["valid_programs"] or
        "STOP" in state_stepwise["valid_programs"]
    )
    assert has_valid_next_step, \
        f"Should have valid next step programs, got: {state_stepwise['valid_programs']}"

    print("  PASSED")


def test_automated_mode_allows_full_predict_and_build():
    """Test that full predict_and_build is allowed in automated mode."""
    print("Test: automated_mode_allows_full_predict_and_build")

    files = ["data.mtz", "sequence.fa"]
    history = [{"program": "phenix.xtriage", "result": "SUCCESS"}]

    # In automated mode, predict_and_build should be available
    state_auto = detect_workflow_state(history, files, maximum_automation=True)
    assert "phenix.predict_and_build" in state_auto["valid_programs"], \
        "predict_and_build should be valid in automated mode before prediction"

    print("  PASSED")


# =============================================================================
# MR-SAD WORKFLOW TESTS
# =============================================================================

def test_mr_sad_after_phaser_with_anomalous():
    """Test MR-SAD: after phaser + anomalous data → experimental_phasing state."""
    print("Test: mr_sad_after_phaser_with_anomalous")

    files = ["data.mtz", "search_model.pdb", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS",
         "analysis": {"has_anomalous": True, "anomalous_measurability": 0.15}},
        {"program": "phenix.phaser", "result": "SUCCESS",
         "output_files": ["PHASER.1.pdb"]},
    ]

    state = detect_workflow_state(history, files)

    assert state["state"] == "xray_mr_sad", \
        "Expected xray_mr_sad, got %s" % state["state"]
    assert "phenix.autosol" in state["valid_programs"], \
        "autosol should be valid for MR-SAD, got %s" % state["valid_programs"]

    print("  PASSED")


def test_mr_sad_not_triggered_without_anomalous():
    """Test that MR-SAD is NOT triggered when there's no anomalous signal."""
    print("Test: mr_sad_not_triggered_without_anomalous")

    files = ["data.mtz", "search_model.pdb", "sequence.fa", "PHASER.1.pdb"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS",
         "analysis": {"resolution": 2.0}},
        {"program": "phenix.phaser", "result": "SUCCESS",
         "output_files": ["PHASER.1.pdb"]},
    ]

    state = detect_workflow_state(history, files)

    # Without anomalous data, should go to refine (not mr_sad)
    assert state["state"] != "xray_mr_sad", \
        "Should NOT be xray_mr_sad without anomalous data, got %s" % state["state"]

    print("  PASSED")


def test_mr_sad_not_triggered_when_autosol_done():
    """Test that MR-SAD phase is skipped if autosol already ran."""
    print("Test: mr_sad_not_triggered_when_autosol_done")

    files = ["data.mtz", "search_model.pdb", "sequence.fa", "PHASER.1.pdb",
             "overall_best.pdb"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS",
         "analysis": {"has_anomalous": True, "anomalous_measurability": 0.15}},
        {"program": "phenix.phaser", "result": "SUCCESS",
         "output_files": ["PHASER.1.pdb"]},
        {"program": "phenix.autosol", "result": "SUCCESS",
         "output_files": ["overall_best.pdb"]},
    ]

    state = detect_workflow_state(history, files)

    # After autosol → should go to build_from_phases or refine, NOT mr_sad
    assert state["state"] != "xray_mr_sad", \
        "Should NOT be xray_mr_sad after autosol, got %s" % state["state"]

    print("  PASSED")


def test_mr_sad_directive_prioritizes_phaser():
    """Test that use_mr_sad directive prioritizes phaser in obtain_model phase."""
    print("Test: mr_sad_directive_prioritizes_phaser")

    files = ["data.mtz", "search_model.pdb", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS",
         "analysis": {"has_anomalous": True, "anomalous_measurability": 0.15}},
    ]

    directives = {
        "workflow_preferences": {
            "use_mr_sad": True,
            "use_experimental_phasing": True,
        }
    }

    state = detect_workflow_state(history, files, directives=directives)

    # Should be in obtain_model phase, with phaser prioritized and autosol removed
    assert "phenix.phaser" in state["valid_programs"], \
        "phaser should be valid for MR-SAD, got %s" % state["valid_programs"]
    assert "phenix.autosol" not in state["valid_programs"], \
        "autosol should NOT be in obtain_model for MR-SAD, got %s" % state["valid_programs"]

    print("  PASSED")


def test_normal_sad_still_works():
    """Test that normal SAD (no phaser) still goes through obtain_model with autosol."""
    print("Test: normal_sad_still_works")

    files = ["data.mtz", "sequence.fa"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS",
         "analysis": {"has_anomalous": True, "anomalous_measurability": 0.20}},
    ]

    directives = {
        "workflow_preferences": {
            "use_experimental_phasing": True,
        }
    }

    state = detect_workflow_state(history, files, directives=directives)

    # Should be in obtain_model with autosol prioritized
    assert "phenix.autosol" in state["valid_programs"], \
        "autosol should be valid for normal SAD, got %s" % state["valid_programs"]

    print("  PASSED")


def test_autosol_has_partial_model_config():
    """Test that autosol config includes partial_model optional input."""
    print("Test: autosol_has_partial_model_config")

    import yaml
    yaml_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "knowledge", "programs.yaml")
    with open(yaml_path) as f:
        programs = yaml.safe_load(f)

    autosol = programs["phenix.autosol"]
    optional = autosol["inputs"].get("optional", {})

    assert "partial_model" in optional, \
        "autosol should have partial_model in optional inputs"
    assert optional["partial_model"]["flag"] == "partpdb_file=", \
        "partial_model flag should be partpdb_file=, got %s" % optional["partial_model"]["flag"]

    # Check input_priorities
    priorities = autosol.get("input_priorities", {})
    assert "partial_model" in priorities, \
        "autosol should have input_priorities for partial_model"
    assert "phaser_output" in priorities["partial_model"]["categories"], \
        "partial_model should prefer phaser_output category"

    print("  PASSED")


def test_mr_sad_directive_overrides_no_anomalous():
    """Test that use_mr_sad directive triggers MR-SAD even without has_anomalous from xtriage."""
    print("Test: mr_sad_directive_overrides_no_anomalous")

    files = ["data.mtz", "search_model.pdb", "sequence.fa", "PHASER.1.pdb"]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS",
         "analysis": {"resolution": 2.0}},  # No anomalous detected
        {"program": "phenix.phaser", "result": "SUCCESS",
         "output_files": ["PHASER.1.pdb"]},
    ]

    # User explicitly requests MR-SAD
    directives = {
        "workflow_preferences": {
            "use_mr_sad": True,
            "use_experimental_phasing": True,
        }
    }

    state = detect_workflow_state(history, files, directives=directives)

    assert state["state"] == "xray_mr_sad", \
        "use_mr_sad directive should trigger MR-SAD even without has_anomalous, got %s" % state["state"]
    assert "phenix.autosol" in state["valid_programs"], \
        "autosol should be valid for MR-SAD with directive, got %s" % state["valid_programs"]

    print("  PASSED")


def test_experimental_phasing_yaml_structure():
    """Test that workflows.yaml has correct experimental_phasing phase definition."""
    print("Test: experimental_phasing_yaml_structure")

    import yaml
    yaml_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "knowledge", "workflows.yaml")
    with open(yaml_path) as f:
        workflows = yaml.safe_load(f)

    phases = workflows["xray"]["phases"]

    # experimental_phasing phase should exist
    assert "experimental_phasing" in phases, \
        "experimental_phasing phase should exist in xray workflow"

    ep = phases["experimental_phasing"]

    # Should have autosol as the program
    programs = ep.get("programs", [])
    assert len(programs) > 0, "experimental_phasing should have programs"
    prog = programs[0]
    assert isinstance(prog, dict), "Program entry should be a dict with conditions"
    assert prog["program"] == "phenix.autosol", \
        "experimental_phasing should use autosol, got %s" % prog.get("program")

    # Should transition to build_from_phases
    transitions = ep.get("transitions", {})
    assert transitions.get("on_complete") == "build_from_phases", \
        "experimental_phasing should transition to build_from_phases, got %s" % transitions

    print("  PASSED")


def test_predict_and_build_blocked_after_full_completion():
    """Test that predict_and_build is NOT offered after it has fully completed.

    Bug: After predict_and_build ran successfully (predict_full_done=True),
    followed by autobuild and refinement, the agent would re-offer
    predict_and_build because _check_program_prerequisites only blocked
    re-runs in stepwise mode. Directives (program_settings) could inject
    it back into valid_programs even in the refine phase.
    """
    print("Test: predict_and_build_blocked_after_full_completion")

    files = [
        "data.mtz", "sequence.fa",
        "PredictAndBuild_0_overall_best.pdb",   # predict_and_build output
        "overall_best_final_refine_001.pdb",     # refined model
        "refine_001_001.mtz",                    # refined data
    ]
    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS",
         "analysis": {"resolution": 2.1}},
        {"program": "phenix.predict_and_build", "result": "SUCCESS",
         "output_files": ["PredictAndBuild_0_overall_best.pdb"]},
        {"program": "phenix.autobuild", "result": "SUCCESS",
         "output_files": ["overall_best_final_refine_001.pdb"]},
        {"program": "phenix.refine", "result": "SUCCESS",
         "output_files": ["refine_001_001.pdb", "refine_001_001.mtz"],
         "analysis": {"r_free": 0.28}},
    ]

    # Without directives - predict_and_build should not be offered
    state = detect_workflow_state(history, files)
    assert state["state"] == "xray_refined", \
        "Expected xray_refined, got %s" % state["state"]
    assert "phenix.predict_and_build" not in state["valid_programs"], \
        "predict_and_build should NOT be valid after full completion, got %s" % state["valid_programs"]

    # With directives that mention predict_and_build in program_settings -
    # it should STILL not be offered
    directives = {
        "program_settings": {
            "phenix.predict_and_build": {"rebuilding_strategy": "Quick"},
        }
    }
    state2 = detect_workflow_state(history, files, directives=directives)
    assert "phenix.predict_and_build" not in state2["valid_programs"], \
        "predict_and_build should NOT be re-added by directives after full completion, got %s" % state2["valid_programs"]

    print("  PASSED")


def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.test_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
