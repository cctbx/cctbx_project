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

from libtbx.langchain.agent.workflow_state import (
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

    assert "p9.sca" in categories["mtz"], "Scalepack file should be in 'mtz' category"
    assert categories["sequence"] == ["sequence.fa"]

    print("  PASSED")


def test_hkl_file_recognition():
    """Test that .hkl files are recognized as X-ray data."""
    print("Test: hkl_file_recognition")

    files = ["data.hkl", "model.pdb"]
    categories = _categorize_files(files)

    assert "data.hkl" in categories["mtz"], "HKL file should be in 'mtz' category"

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
    # Validation should be suggested for good model
    assert "phenix.molprobity" in state2["valid_programs"], \
        "Expected validation for good model, got %s" % state2["valid_programs"]

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
    assert state["valid_programs"] == ["phenix.mtriage"]
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

    assert state["state"] in ["cryoem_has_model", "cryoem_refined", "refine"]
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
    assert state["state"] in ["cryoem_has_phases", "cryoem_has_model", "cryoem_refined", "refine"]
    assert "phenix.real_space_refine" in state["valid_programs"], \
        "Expected real_space_refine, got %s" % state["valid_programs"]

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

    assert "data.mtz" in categorized["mtz"]
    assert "PHASER.1.pdb" in categorized["phaser_output"]
    assert "model_refine_001.pdb" in categorized["refined"]
    assert "protein_with_ligand.pdb" in categorized["with_ligand"]
    assert "ligand_fit_1.pdb" in categorized["ligand_fit"]
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

    # Test 6: processed predicted model should be in BOTH predicted and processed
    files6 = ["processed_predicted_model.pdb"]
    cat6 = _categorize_files(files6)
    assert "processed_predicted_model.pdb" in cat6["predicted"]
    assert "processed_predicted_model.pdb" in cat6["processed_predicted"]

    # Test 7: CIF files - distinguish model CIF from ligand CIF
    files7 = ["ligand.cif", "ATP.cif", "refine_001.cif", "PHASER.1_refine.cif"]
    cat7 = _categorize_files(files7)
    # Ligand CIFs
    assert "ligand.cif" in cat7["ligand_cif"]
    assert "ATP.cif" in cat7["ligand_cif"]
    # Model CIFs (from refinement) should be in pdb/refined, NOT ligand_cif
    assert "refine_001.cif" in cat7["pdb"]
    assert "refine_001.cif" in cat7["refined"]
    assert "refine_001.cif" not in cat7["ligand_cif"]
    assert "PHASER.1_refine.cif" in cat7["pdb"]

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
    """Test that dock_in_map is available when model exists."""
    print("Test: dock_in_map_option")

    files = ["map.mrc", "existing_model.pdb"]
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


def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("WORKFLOW STATE TESTS")
    print("=" * 60)

    # File categorization tests (new)
    test_dat_sequence_file_recognition()
    test_xray_analyzed_with_dat_sequence()
    test_scalepack_file_recognition()
    test_hkl_file_recognition()
    test_xray_initial_with_sca()

    # X-ray tests
    test_xray_initial()
    test_xray_analyzed_with_sequence()
    test_xray_analyzed_with_model()
    test_xray_has_prediction()
    test_xray_predict_full_places_model()
    test_xray_model_processed()
    test_xray_analyzed_stuck()
    test_xray_has_model()
    test_xray_refined_no_ligand()
    test_xray_refined_with_ligand()
    test_xray_has_ligand()
    test_xray_combined()
    test_xray_refined_autobuild_option()
    test_autobuild_priority_over_ligandfit()

    # Cryo-EM automated tests
    test_cryoem_initial()
    test_cryoem_analyzed_automated()
    test_cryoem_has_model_automated()
    test_cryoem_refined()

    # Cryo-EM stepwise tests
    test_cryoem_analyzed_stepwise()
    test_cryoem_has_prediction_stepwise()
    test_cryoem_model_processed_stepwise()
    test_cryoem_has_phases_stepwise()

    # Validation tests
    test_validate_valid_program()
    test_validate_invalid_program()
    test_validate_stop()
    test_validate_unknown_program()

    # Edge case tests
    test_mixed_xray_cryoem_files()
    test_file_categorization()
    test_file_categorization_edge_cases()
    test_rsr_output_categorization()
    test_history_analysis()
    test_format_workflow_for_prompt()
    test_dock_in_map_option()
    test_ligandfit_with_cryoem()

    print("=" * 60)
    print("ALL WORKFLOW STATE TESTS PASSED!")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
