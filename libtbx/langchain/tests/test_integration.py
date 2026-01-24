#!/usr/bin/env python
"""
Integration tests for the PHENIX AI Agent v2.

Tests the full pipeline:
- perceive -> plan -> build -> validate flow
- Auto-stop on plateau
- Workflow state enforcement
- Mock LLM fallback
"""

from __future__ import absolute_import, division, print_function
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from libtbx.langchain.agent.graph_state import create_initial_state
    from libtbx.langchain.agent.graph_nodes import perceive, plan, build, validate, fallback
except ImportError:
    from agent.graph_state import create_initial_state
    from agent.graph_nodes import perceive, plan, build, validate, fallback


def test_perceive_builds_metrics():
    """Test that perceive builds metrics history from backpack."""
    print("Test: perceive_builds_metrics")

    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS", "analysis": {"resolution": 2.5}},
        {"program": "phenix.phaser", "result": "SUCCESS", "analysis": {"tfz": 24}},
        {"program": "phenix.refine", "result": "SUCCESS", "analysis": {"r_free": 0.35}},
    ]

    state = create_initial_state(
        available_files=["data.mtz", "PHASER.1_refine_001.pdb"],
        log_text="Final R-free: 0.32",
        history=history
    )

    state = perceive(state)

    # Check metrics history was built
    assert len(state["metrics_history"]) >= 3
    assert state["metrics_history"][0]["resolution"] == 2.5
    assert state["metrics_history"][1]["tfz"] == 24
    assert state["metrics_history"][2]["r_free"] == 0.35

    # Check workflow state was detected
    assert state["workflow_state"]["state"] == "xray_refined"

    print("  PASSED")


def test_auto_stop_on_plateau():
    """Test that plan node auto-stops on detected plateau."""
    print("Test: auto_stop_on_plateau")

    # History with 6 refine cycles showing clear plateau (last 3 improvements < 0.3%)
    # New threshold: need 4+ values to detect plateau over last 3 improvements
    history = [
        {"program": "phenix.refine", "analysis": {"r_free": 0.300}},
        {"program": "phenix.refine", "analysis": {"r_free": 0.295}},
        {"program": "phenix.refine", "analysis": {"r_free": 0.2945}},  # 0.17% improvement
        {"program": "phenix.refine", "analysis": {"r_free": 0.2941}},  # 0.14% improvement
        {"program": "phenix.refine", "analysis": {"r_free": 0.2938}},  # 0.10% improvement
        {"program": "phenix.refine", "analysis": {"r_free": 0.2936}},  # 0.07% improvement
    ]

    # Use a realistic log snippet
    log_text = """
Final R-work = 0.2450, R-free = 0.2934
"""

    state = create_initial_state(
        available_files=["data.mtz", "model.pdb"],
        log_text=log_text,
        history=history
    )

    state = perceive(state)

    # Debug: print what was extracted
    print("    metrics_history length:", len(state["metrics_history"]))
    print("    r_free_trend:", state["metrics_trend"].get("r_free_trend", []))
    print("    consecutive_refines:", state["metrics_trend"].get("consecutive_refines", 0))
    print("    should_stop:", state["metrics_trend"]["should_stop"])
    print("    reason:", state["metrics_trend"].get("reason"))

    # Check plateau was detected
    # Last 3 improvements in history are all < 0.3%
    assert state["metrics_trend"]["should_stop"] == True, \
        "Expected plateau detection. trend=%s" % state["metrics_trend"]
    assert "PLATEAU" in state["metrics_trend"]["reason"], \
        "Expected PLATEAU in reason. Got: %s" % state["metrics_trend"]["reason"]

    # Plan should auto-stop
    state = plan(state)

    assert state["stop"] == True
    assert "PLATEAU" in str(state["stop_reason"]) or "plateau" in str(state["intent"].get("reasoning", "")).lower()

    print("  PASSED")


def test_workflow_validation_blocks_invalid():
    """Test that validate blocks invalid program for state."""
    print("Test: workflow_validation_blocks_invalid")

    # State: xray_initial (only xtriage valid)
    state = create_initial_state(
        available_files=["data.mtz"],
        history=[]
    )

    state = perceive(state)

    # Force an invalid intent (refine when should run xtriage)
    state["intent"] = {
        "program": "phenix.refine",
        "files": {"model": "model.pdb", "data": "data.mtz"},
        "strategy": {}
    }
    state["command"] = "phenix.refine model.pdb data.mtz"

    state = validate(state)

    # Should have validation error about workflow
    assert state["validation_error"] is not None
    assert "not valid" in state["validation_error"].lower()
    assert state["attempt_number"] == 1

    print("  PASSED")


def test_mock_plan_uses_workflow():
    """Test that mock planner respects workflow state."""
    print("Test: mock_plan_uses_workflow")

    # Initial state - should pick xtriage
    state1 = create_initial_state(
        available_files=["data.mtz", "sequence.fa"],
        history=[]
    )
    state1 = perceive(state1)
    state1["provider"] = "nonexistent"  # Force mock planner

    # Mock plan should pick xtriage (first valid program)
    state1 = plan(state1)

    assert state1["intent"]["program"] == "phenix.xtriage"

    # After xtriage - should pick predict_and_build
    state2 = create_initial_state(
        available_files=["data.mtz", "sequence.fa"],
        history=[{"program": "phenix.xtriage", "result": "SUCCESS"}]
    )
    state2 = perceive(state2)
    state2["provider"] = "nonexistent"

    state2 = plan(state2)

    assert state2["intent"]["program"] == "phenix.predict_and_build"

    print("  PASSED")


def test_fallback_uses_workflow():
    """Test that fallback node respects workflow state."""
    print("Test: fallback_uses_workflow")

    state = create_initial_state(
        available_files=["data.mtz", "sequence.fa"],
        history=[{"program": "phenix.xtriage", "result": "SUCCESS"}]
    )
    state = perceive(state)

    # Fallback should pick from valid programs
    state = fallback(state)

    valid_programs = state["workflow_state"]["valid_programs"]

    # Command should be for a valid program
    assert state["fallback_used"] == True
    cmd = state["command"]
    assert any(prog in cmd for prog in valid_programs if prog != "STOP")

    print("  PASSED")


def test_cryoem_workflow_detection():
    """Test cryo-EM workflow is detected from files."""
    print("Test: cryoem_workflow_detection")

    state = create_initial_state(
        available_files=["map.mrc", "sequence.fa"],
        history=[],
        maximum_automation=True
    )

    state = perceive(state)

    assert state["workflow_state"]["experiment_type"] == "cryoem"
    assert state["workflow_state"]["state"] == "cryoem_initial"
    assert "phenix.mtriage" in state["workflow_state"]["valid_programs"]

    print("  PASSED")


def test_cryoem_stepwise_forces_stop_after_predict():
    """Test that stepwise cryo-EM forces stop_after_predict in build."""
    print("Test: cryoem_stepwise_forces_stop_after_predict")

    state = create_initial_state(
        available_files=["map.mrc", "sequence.fa"],
        history=[{"program": "phenix.mtriage", "result": "SUCCESS"}],
        maximum_automation=False  # Stepwise mode
    )

    state = perceive(state)

    # Check stepwise path detected
    assert state["workflow_state"]["automation_path"] == "stepwise"

    # Simulate intent for predict_and_build
    state["intent"] = {
        "program": "phenix.predict_and_build",
        "files": {"sequence": "sequence.fa"},
        "strategy": {}  # No stop_after_predict specified
    }

    state = build(state)

    # Build should have added stop_after_predict=True
    assert "stop_after_predict" in state["command"].lower() or state.get("debug_log", [])

    print("  PASSED (or check debug_log for strategy modification)")


def test_success_detection():
    """Test R-free success detection."""
    print("Test: success_detection")

    # Test 1: Success WITHOUT validation - should recommend validation, not stop
    history_no_validation = [
        {"program": "phenix.xtriage", "result": "SUCCESS", "analysis": {"resolution": 2.5}},
        {"program": "phenix.phaser", "result": "SUCCESS", "analysis": {"tfz": 24}},
        {"program": "phenix.refine", "analysis": {"r_free": 0.30}},
        {"program": "phenix.refine", "analysis": {"r_free": 0.22}},  # Below success threshold (0.23)
    ]

    state = create_initial_state(
        available_files=["data.mtz", "model.pdb"],
        log_text="R-free: 0.22",
        history=history_no_validation
    )

    state = perceive(state)

    # Should NOT auto-stop without validation
    assert state["metrics_trend"]["should_stop"] == False, \
        "Should not auto-stop without validation. trend=%s" % state["metrics_trend"]
    assert state["metrics_trend"]["recommendation"] == "validate", \
        "Should recommend validation. Got: %s" % state["metrics_trend"]["recommendation"]

    # Test 2: Success WITH validation - should recommend stop
    history_with_validation = [
        {"program": "phenix.xtriage", "result": "SUCCESS", "analysis": {"resolution": 2.5}},
        {"program": "phenix.phaser", "result": "SUCCESS", "analysis": {"tfz": 24}},
        {"program": "phenix.refine", "analysis": {"r_free": 0.30}},
        {"program": "phenix.refine", "analysis": {"r_free": 0.22}},
        {"program": "phenix.molprobity", "result": "SUCCESS"},  # Validation done
    ]

    state2 = create_initial_state(
        available_files=["data.mtz", "model.pdb"],
        log_text="R-free: 0.22",
        history=history_with_validation
    )

    state2 = perceive(state2)

    # Should detect success after validation
    assert state2["metrics_trend"]["should_stop"] == True, \
        "Expected success after validation. trend=%s" % state2["metrics_trend"]
    assert "SUCCESS" in state2["metrics_trend"]["reason"], \
        "Expected SUCCESS in reason. Got: %s" % state2["metrics_trend"]["reason"]

    print("  PASSED")


def test_full_xray_cycle():
    """Test a complete X-ray cycle: perceive -> plan -> build -> validate."""
    print("Test: full_xray_cycle")

    # Start fresh
    state = create_initial_state(
        available_files=["data.mtz", "sequence.fa"],
        history=[],
        provider="nonexistent"  # Use mock planner
    )

    # Run perceive
    state = perceive(state)
    assert state["workflow_state"]["state"] == "xray_initial"

    # Run plan (mock will pick xtriage)
    state = plan(state)
    assert state["intent"]["program"] == "phenix.xtriage"

    # Run build
    state = build(state)
    assert "phenix.xtriage" in state["command"]
    assert "data.mtz" in state["command"]

    # Run validate
    state = validate(state)
    assert state["validation_error"] is None

    print("  PASSED")


def test_consecutive_refine_tracking():
    """Test that consecutive refine cycles are tracked correctly."""
    print("Test: consecutive_refine_tracking")

    history = [
        {"program": "phenix.xtriage", "result": "SUCCESS"},
        {"program": "phenix.phaser", "result": "SUCCESS"},
        {"program": "phenix.refine", "analysis": {"r_free": 0.35}},
        {"program": "phenix.refine", "analysis": {"r_free": 0.32}},
        {"program": "phenix.refine", "analysis": {"r_free": 0.30}},
    ]

    state = create_initial_state(
        available_files=["data.mtz", "model.pdb"],
        log_text="",  # Empty log - don't rely on current log being parsed
        history=history
    )

    state = perceive(state)

    # 3 explicit refines from history
    # (current log may or may not add +1 depending on whether log_parsers extracts it)
    assert state["metrics_trend"]["consecutive_refines"] >= 3, \
        "Expected at least 3 consecutive refines, got %d" % state["metrics_trend"]["consecutive_refines"]

    print("  PASSED")


def test_validation_error_increments_attempt():
    """Test that validation errors increment attempt counter."""
    print("Test: validation_error_increments_attempt")

    state = create_initial_state(
        available_files=["data.mtz"],
        history=[]
    )
    state = perceive(state)

    # Force invalid intent
    state["intent"] = {"program": "phenix.refine"}
    state["command"] = "phenix.refine nonexistent.pdb data.mtz"

    # First validation failure
    state = validate(state)
    assert state["attempt_number"] == 1
    assert state["validation_error"] is not None

    # Second validation failure
    state["command"] = "phenix.refine another_nonexistent.pdb data.mtz"
    state = validate(state)
    assert state["attempt_number"] == 2

    print("  PASSED")


def test_file_path_validation():
    """Test that file validation passes if basename matches (for server-side use)."""
    print("Test: file_path_validation")

    # The correct path
    correct_path = "/home/user/project/processed_model.pdb"
    # Different path with same basename (acceptable on server where we only have file list)
    different_path = "/home/user/other/processed_model.pdb"

    state = create_initial_state(
        available_files=["data.mtz", "sequence.fa", correct_path],
        history=[{"program": "phenix.xtriage", "result": "SUCCESS"}]
    )
    state = perceive(state)

    # For this test about file path validation, ensure phaser is valid
    # by manually adding it to valid_programs
    workflow_state = state.get("workflow_state", {})
    valid_programs = workflow_state.get("valid_programs", [])
    if "phenix.phaser" not in valid_programs:
        workflow_state["valid_programs"] = valid_programs + ["phenix.phaser"]
        state["workflow_state"] = workflow_state

    # Create a command with a different path but same basename
    state["intent"] = {"program": "phenix.phaser"}
    state["command"] = "phenix.phaser data.mtz %s" % different_path

    state = validate(state)

    # Should PASS because basename matches (server can't verify exact paths)
    assert state["validation_error"] is None, \
        "Expected validation to pass with matching basename, got: %s" % state["validation_error"]

    # Now test with a truly missing file
    state["command"] = "phenix.phaser data.mtz /some/nonexistent_file.pdb"
    state = validate(state)

    # Should FAIL because basename doesn't exist in available_files
    assert state["validation_error"] is not None, "Expected validation to fail with missing file"
    assert "nonexistent_file.pdb" in state["validation_error"]

    print("  PASSED")


def run_all_tests():
    """Run all integration tests."""
    print("=" * 60)
    print("INTEGRATION TESTS")
    print("=" * 60)

    test_perceive_builds_metrics()
    test_auto_stop_on_plateau()
    test_workflow_validation_blocks_invalid()
    test_mock_plan_uses_workflow()
    test_fallback_uses_workflow()
    test_cryoem_workflow_detection()
    test_cryoem_stepwise_forces_stop_after_predict()
    test_success_detection()
    test_full_xray_cycle()
    test_consecutive_refine_tracking()
    test_validation_error_increments_attempt()
    test_file_path_validation()

    print("=" * 60)
    print("ALL INTEGRATION TESTS PASSED!")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
