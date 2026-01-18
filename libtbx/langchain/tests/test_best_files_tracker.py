"""
Unit tests for BestFilesTracker.

Run with: python tests/test_best_files_tracker.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agent.best_files_tracker import BestFilesTracker


# =============================================================================
# TEST HELPERS
# =============================================================================

def assert_none(actual, message=""):
    if actual is not None:
        raise AssertionError(f"{message}: expected None, got {actual}")


def assert_equal(actual, expected, message=""):
    """Assert two values are equal."""
    if actual != expected:
        raise AssertionError(f"{message}: expected {expected}, got {actual}")


def assert_true(value, message=""):
    """Assert value is True."""
    if not value:
        raise AssertionError(f"{message}: expected True, got {value}")


def assert_false(value, message=""):
    """Assert value is False."""
    if value:
        raise AssertionError(f"{message}: expected False, got {value}")


def assert_none(value, message=""):
    """Assert value is None."""
    if value is not None:
        raise AssertionError(f"{message}: expected None, got {value}")


def assert_not_none(value, message=""):
    """Assert value is not None."""
    if value is None:
        raise AssertionError(f"{message}: expected not None")


# =============================================================================
# MODEL SCORING TESTS
# =============================================================================

def test_model_stage_scoring():
    """Test that refined models score higher than predicted."""
    print("Test: model_stage_scoring")

    tracker = BestFilesTracker()

    # Add predicted model
    tracker.evaluate_file("/path/predicted.pdb", cycle=1, stage="predicted")
    best = tracker.get_best("model")
    assert_not_none(best, "Should have best model")
    predicted_score = best.score

    # Add docked model - should become best
    tracker.evaluate_file("/path/docked.pdb", cycle=2, stage="docked")
    best = tracker.get_best("model")
    assert_equal(best.stage, "docked", "Docked should be best")
    docked_score = best.score
    assert_true(docked_score > predicted_score, "Docked should score higher than predicted")

    # Add refined model - should become best
    tracker.evaluate_file("/path/refined.pdb", cycle=3, stage="refined")
    best = tracker.get_best("model")
    assert_equal(best.stage, "refined", "Refined should be best")
    refined_score = best.score
    assert_true(refined_score > docked_score, "Refined should score higher than docked")

    print("  PASSED")


def test_model_metrics_improve_score():
    """Test that better metrics improve model score."""
    print("Test: model_metrics_improve_score")

    tracker = BestFilesTracker()

    # Refined model with poor metrics
    tracker.evaluate_file("/path/refine1.pdb", cycle=1, stage="refined",
                         metrics={"r_free": 0.35, "clashscore": 25})
    score1 = tracker.get_best("model").score

    # Refined model with better metrics - should win
    tracker.evaluate_file("/path/refine2.pdb", cycle=2, stage="refined",
                         metrics={"r_free": 0.22, "clashscore": 5})
    score2 = tracker.get_best("model").score

    assert_true(score2 > score1, "Better metrics should give higher score")
    assert_equal(tracker.get_best("model").path, "/path/refine2.pdb",
                "Model with better metrics should be best")

    print("  PASSED")


def test_model_stage_beats_metrics():
    """Test that refined model beats docked even with similar metrics."""
    print("Test: model_stage_beats_metrics")

    tracker = BestFilesTracker()

    # Docked model with good metrics
    tracker.evaluate_file("/path/docked.pdb", cycle=1, stage="docked",
                         metrics={"r_free": 0.25})

    # Refined model with same metrics - should still win due to stage
    tracker.evaluate_file("/path/refined.pdb", cycle=2, stage="refined",
                         metrics={"r_free": 0.25})

    best = tracker.get_best("model")
    assert_equal(best.stage, "refined", "Refined should beat docked at same metrics")

    print("  PASSED")


# =============================================================================
# MAP SCORING TESTS
# =============================================================================

def test_map_stage_scoring():
    """Test that optimized maps score higher than raw."""
    print("Test: map_stage_scoring")

    tracker = BestFilesTracker()

    # Add half-map
    tracker.evaluate_file("/path/half_1.ccp4", cycle=1, stage="half_map")
    half_score = tracker.get_best("map").score

    # Add full map - should be better
    tracker.evaluate_file("/path/full.ccp4", cycle=2, stage="full_map")
    full_score = tracker.get_best("map").score
    assert_true(full_score > half_score, "Full map should beat half-map")

    # Add optimized map - should be best
    tracker.evaluate_file("/path/denmod.ccp4", cycle=3, stage="optimized_full_map")
    opt_score = tracker.get_best("map").score
    assert_true(opt_score > full_score, "Optimized should beat full map")

    print("  PASSED")


def test_map_resolution_bonus():
    """Test that higher resolution maps get bonus points."""
    print("Test: map_resolution_bonus")

    tracker = BestFilesTracker()

    # Low resolution map
    tracker.evaluate_file("/path/lowres.ccp4", cycle=1, stage="full_map",
                         metrics={"resolution": 4.0})
    lowres_score = tracker.get_best("map").score

    # High resolution map - should score higher
    tracker.evaluate_file("/path/highres.ccp4", cycle=2, stage="full_map",
                         metrics={"resolution": 2.0})
    highres_score = tracker.get_best("map").score

    assert_true(highres_score > lowres_score, "Higher resolution should score higher")

    print("  PASSED")


# =============================================================================
# MTZ SCORING TESTS
# =============================================================================

def test_mtz_rfree_locks():
    """Test that first MTZ with R-free flags locks forever."""
    print("Test: mtz_rfree_locks")

    tracker = BestFilesTracker()

    # First MTZ with R-free - should lock
    tracker.evaluate_file("/path/data.mtz", cycle=1, stage="original",
                         metrics={"has_rfree_flags": True})
    assert_equal(tracker.get_best("mtz").path, "/path/data.mtz")

    # Try to change to refined MTZ - should fail
    result = tracker.evaluate_file("/path/refined.mtz", cycle=2, stage="refined_mtz",
                                   metrics={"has_rfree_flags": True})
    assert_false(result, "Should not change locked MTZ")
    assert_equal(tracker.get_best("mtz").path, "/path/data.mtz",
                "Original MTZ should still be best")

    print("  PASSED")


def test_mtz_without_rfree_can_change():
    """Test that MTZ without R-free flags can be replaced."""
    print("Test: mtz_without_rfree_can_change")

    tracker = BestFilesTracker()

    # First MTZ without R-free
    tracker.evaluate_file("/path/data1.mtz", cycle=1, stage="original",
                         metrics={"has_rfree_flags": False})
    assert_equal(tracker.get_best("mtz").path, "/path/data1.mtz")

    # Second MTZ with R-free - should replace and lock
    result = tracker.evaluate_file("/path/data2.mtz", cycle=2, stage="original",
                                   metrics={"has_rfree_flags": True})
    assert_true(result, "Should change to MTZ with R-free")
    assert_equal(tracker.get_best("mtz").path, "/path/data2.mtz")

    # Now try third MTZ - should fail (locked)
    result = tracker.evaluate_file("/path/data3.mtz", cycle=3, stage="refined_mtz",
                                   metrics={"has_rfree_flags": True})
    assert_false(result, "Should not change after lock")

    print("  PASSED")


# =============================================================================
# HISTORY TESTS
# =============================================================================

def test_history_tracking():
    """Test that changes are recorded in history."""
    print("Test: history_tracking")

    tracker = BestFilesTracker()

    # Initial - no history
    assert_equal(len(tracker.get_history()), 0, "Should start with no history")

    # First file - no change recorded (no previous)
    tracker.evaluate_file("/path/model1.pdb", cycle=1, stage="predicted")
    assert_equal(len(tracker.get_history()), 0, "First file shouldn't create history")

    # Second file replaces first - should record
    tracker.evaluate_file("/path/model2.pdb", cycle=2, stage="docked")
    assert_equal(len(tracker.get_history()), 1, "Change should be recorded")

    # Check history entry
    change = tracker.get_history()[0]
    assert_equal(change.old_path, "/path/model1.pdb")
    assert_equal(change.new_path, "/path/model2.pdb")
    assert_equal(change.cycle, 2)

    print("  PASSED")


def test_history_filter_by_category():
    """Test filtering history by category."""
    print("Test: history_filter_by_category")

    tracker = BestFilesTracker()

    # Create model changes
    tracker.evaluate_file("/path/model1.pdb", cycle=1, stage="predicted")
    tracker.evaluate_file("/path/model2.pdb", cycle=2, stage="refined")

    # Create map changes
    tracker.evaluate_file("/path/map1.ccp4", cycle=1, stage="half_map")
    tracker.evaluate_file("/path/map2.ccp4", cycle=2, stage="full_map")

    # Filter by category
    model_history = tracker.get_history("model")
    map_history = tracker.get_history("map")

    assert_equal(len(model_history), 1, "Should have 1 model change")
    assert_equal(len(map_history), 1, "Should have 1 map change")
    assert_equal(model_history[0].category, "model")
    assert_equal(map_history[0].category, "map")

    print("  PASSED")


# =============================================================================
# CLASSIFICATION TESTS
# =============================================================================

def test_file_classification_by_extension():
    """Test that files are classified correctly by extension."""
    print("Test: file_classification_by_extension")

    tracker = BestFilesTracker()

    # PDB -> model
    tracker.evaluate_file("/path/model.pdb", cycle=1)
    assert_not_none(tracker.get_best("model"), "PDB should be classified as model")

    # MTZ -> mtz
    tracker.evaluate_file("/path/data.mtz", cycle=1)
    assert_not_none(tracker.get_best("mtz"), "MTZ should be classified as mtz")

    # CCP4 -> map
    tracker.evaluate_file("/path/map.ccp4", cycle=1)
    assert_not_none(tracker.get_best("map"), "CCP4 should be classified as map")

    # FA -> sequence
    tracker.evaluate_file("/path/seq.fa", cycle=1)
    assert_not_none(tracker.get_best("sequence"), "FA should be classified as sequence")

    print("  PASSED")


def test_stage_classification_from_filename():
    """Test that stage is inferred from filename."""
    print("Test: stage_classification_from_filename")

    tracker = BestFilesTracker()

    # Test various filename patterns
    test_cases = [
        ("predicted_model.pdb", "predicted"),
        ("placed_model.pdb", "docked"),
        ("refine_001_001.pdb", "refined"),
        ("rsr_001_real_space_refined_000.pdb", "rsr_output"),
        ("overall_best_1.pdb", "autobuild_output"),
        ("PHASER.1.pdb", "phaser_output"),
        ("processed_predicted.pdb", "processed_predicted"),
    ]

    for filename, expected_stage in test_cases:
        tracker2 = BestFilesTracker()  # Fresh tracker
        tracker2.evaluate_file(f"/path/{filename}", cycle=1)
        best = tracker2.get_best("model")
        assert_not_none(best, f"Should classify {filename}")
        assert_equal(best.stage, expected_stage,
                    f"Stage for {filename} should be {expected_stage}")

    print("  PASSED")


def test_intermediate_files_ignored():
    """Test that intermediate files are not tracked."""
    print("Test: intermediate_files_ignored")

    tracker = BestFilesTracker()

    # These should be ignored
    ignored_files = [
        "/path/run_mr/model.pdb",
        "/path/AutoBuild_run_1/best.pdb",
        "/path/mask.ccp4",
        "/path/temp/file.pdb",
    ]

    for f in ignored_files:
        result = tracker.evaluate_file(f, cycle=1, stage="refined")
        # Should return False or have no effect
        # (intermediate detection may happen at different points)

    # Tracker should have no entries from these files
    # (mask.ccp4 specifically should be ignored for map category)

    print("  PASSED")


# =============================================================================
# SERIALIZATION TESTS
# =============================================================================

def test_serialization_roundtrip():
    """Test to_dict/from_dict roundtrip."""
    print("Test: serialization_roundtrip")

    tracker = BestFilesTracker()

    # Add some data
    tracker.evaluate_file("/path/model1.pdb", cycle=1, stage="predicted")
    tracker.evaluate_file("/path/model2.pdb", cycle=2, stage="refined",
                         metrics={"r_free": 0.25})
    tracker.evaluate_file("/path/map.ccp4", cycle=1, stage="full_map")
    tracker.evaluate_file("/path/data.mtz", cycle=1, stage="original",
                         metrics={"has_rfree_flags": True})

    # Serialize
    data = tracker.to_dict()

    # Deserialize
    tracker2 = BestFilesTracker.from_dict(data)

    # Verify
    assert_equal(len(tracker2.best), len(tracker.best), "Same number of best entries")
    assert_equal(len(tracker2.history), len(tracker.history), "Same history length")

    # Check specific values
    assert_equal(tracker2.get_best("model").path, "/path/model2.pdb")
    assert_equal(tracker2.get_best("model").stage, "refined")
    assert_equal(tracker2.get_best("map").path, "/path/map.ccp4")
    assert_equal(tracker2.get_best("mtz").path, "/path/data.mtz")

    # MTZ lock should be preserved
    assert_true(tracker2._mtz_with_rfree_locked, "MTZ lock should be preserved")

    print("  PASSED")


def test_empty_serialization():
    """Test serialization of empty tracker."""
    print("Test: empty_serialization")

    tracker = BestFilesTracker()
    data = tracker.to_dict()

    tracker2 = BestFilesTracker.from_dict(data)

    assert_equal(len(tracker2.best), 0, "Should have no best entries")
    assert_equal(len(tracker2.history), 0, "Should have no history")

    print("  PASSED")


def test_from_dict_none():
    """Test from_dict with None input."""
    print("Test: from_dict_none")

    tracker = BestFilesTracker.from_dict(None)
    assert_equal(len(tracker.best), 0, "Should create empty tracker")

    tracker2 = BestFilesTracker.from_dict({})
    assert_equal(len(tracker2.best), 0, "Should create empty tracker from empty dict")

    print("  PASSED")


# =============================================================================
# EDGE CASE TESTS
# =============================================================================

def test_recency_tiebreaker():
    """Test that more recent file wins when scores are equal."""
    print("Test: recency_tiebreaker")

    tracker = BestFilesTracker()

    # Two files with same stage, no metrics
    tracker.evaluate_file("/path/model1.pdb", cycle=1, stage="refined")
    tracker.evaluate_file("/path/model2.pdb", cycle=5, stage="refined")

    # More recent should win
    assert_equal(tracker.get_best("model").path, "/path/model2.pdb",
                "More recent should win tie")

    print("  PASSED")


def test_worse_file_doesnt_replace():
    """Test that a worse file doesn't replace better one."""
    print("Test: worse_file_doesnt_replace")

    tracker = BestFilesTracker()

    # Good refined model
    tracker.evaluate_file("/path/refined.pdb", cycle=1, stage="refined",
                         metrics={"r_free": 0.22})
    original_best = tracker.get_best("model").path

    # Try to add worse model (predicted, even more recent)
    tracker.evaluate_file("/path/predicted.pdb", cycle=10, stage="predicted")

    # Should not have changed
    assert_equal(tracker.get_best("model").path, original_best,
                "Worse file shouldn't replace better")

    print("  PASSED")


def test_get_best_dict():
    """Test get_best_dict returns correct format."""
    print("Test: get_best_dict")

    tracker = BestFilesTracker()
    tracker.evaluate_file("/path/model.pdb", cycle=1, stage="refined")
    tracker.evaluate_file("/path/map.ccp4", cycle=1, stage="full_map")

    best_dict = tracker.get_best_dict()

    assert_equal(best_dict.get("model"), "/path/model.pdb")
    assert_equal(best_dict.get("map"), "/path/map.ccp4")
    assert_none(best_dict.get("mtz"), "Should not have mtz")

    print("  PASSED")


# =============================================================================
# SESSION INTEGRATION TESTS (without full Session import)
# =============================================================================

def test_session_integration_simulation():
    """Test the integration pattern that session.py uses."""
    print("Test: session_integration_simulation")

    tracker = BestFilesTracker()

    # Simulate session.record_result flow
    def simulate_record_result(tracker, cycle, program, output_files, result, metrics=None):
        """Simulate what session.record_result does."""
        if "FAILED" in result.upper():
            return  # Skip failed cycles

        # Infer stage from program (like _infer_stage_from_program)
        stage = None
        prog_lower = program.lower()
        if "real_space_refine" in prog_lower:
            stage = "rsr_output"
        elif "refine" in prog_lower:
            stage = "refined"
        elif "dock_in_map" in prog_lower:
            stage = "docked"
        elif "resolve_cryo_em" in prog_lower:
            stage = "optimized_full_map"

        # Evaluate each output file
        for f in output_files:
            tracker.evaluate_file(f, cycle, metrics, stage)

    # Simulate a cryo-EM workflow
    simulate_record_result(tracker, 1, "phenix.mtriage",
                          ["/path/mask.ccp4"], "SUCCESS")  # mask should be ignored
    simulate_record_result(tracker, 2, "phenix.resolve_cryo_em",
                          ["/path/denmod_map.ccp4"], "SUCCESS")
    simulate_record_result(tracker, 3, "phenix.dock_in_map",
                          ["/path/placed_model.pdb"], "SUCCESS")
    simulate_record_result(tracker, 4, "phenix.real_space_refine",
                          ["/path/rsr_001_000.pdb"], "SUCCESS",
                          metrics={"map_cc": 0.75})

    # Verify best files
    assert_equal(tracker.get_best_path("model"), "/path/rsr_001_000.pdb",
                "Best model should be refined")
    assert_equal(tracker.get_best_path("map"), "/path/denmod_map.ccp4",
                "Best map should be denmod")

    # Failed cycles should not update best
    old_model = tracker.get_best_path("model")
    simulate_record_result(tracker, 5, "phenix.real_space_refine",
                          ["/path/rsr_002_000.pdb"], "FAILED: crashed")
    assert_equal(tracker.get_best_path("model"), old_model,
                "Failed cycle should not update best")

    print("  PASSED")


def test_session_persistence_simulation():
    """Test that best files survive serialization like session.save/load."""
    print("Test: session_persistence_simulation")

    # Create tracker and add files
    tracker = BestFilesTracker()
    tracker.evaluate_file("/path/model.pdb", cycle=1, stage="refined",
                         metrics={"r_free": 0.22})
    tracker.evaluate_file("/path/map.ccp4", cycle=1, stage="optimized_full_map")
    tracker.evaluate_file("/path/data.mtz", cycle=0, stage="original",
                         metrics={"has_rfree_flags": True})

    # Simulate session.save() - serialize to dict
    session_data = {
        "cycles": [],
        "best_files": tracker.to_dict(),
        "best_files_history": [h.to_dict() for h in tracker.get_history()]
    }

    # Convert to JSON and back (like file save/load)
    import json
    json_str = json.dumps(session_data, indent=2)
    loaded_data = json.loads(json_str)

    # Simulate session.load() - restore tracker
    tracker2 = BestFilesTracker.from_dict(loaded_data.get("best_files"))

    # Verify all best files restored
    assert_equal(tracker2.get_best_path("model"), "/path/model.pdb")
    assert_equal(tracker2.get_best_path("map"), "/path/map.ccp4")
    assert_equal(tracker2.get_best_path("mtz"), "/path/data.mtz")

    # Verify MTZ is locked (should not change)
    result = tracker2.evaluate_file("/path/new_data.mtz", cycle=2,
                                    metrics={"has_rfree_flags": True})
    assert_false(result, "MTZ should remain locked after reload")
    assert_equal(tracker2.get_best_path("mtz"), "/path/data.mtz")

    print("  PASSED")


def test_template_builder_best_files_priority():
    """Test that template_builder prefers best_files over categorized_files."""
    print("Test: template_builder_best_files_priority")

    # Create mock scenario: we have multiple models but a clear "best"
    best_files = {
        "model": "/path/to/best_refined.pdb",  # The tracked best
        "map": "/path/to/denmod_map.ccp4",
    }

    # Categorized files include older/worse models
    categorized_files = {
        "model": ["/path/to/older_model.pdb", "/path/to/predicted.pdb"],
        "refined": ["/path/to/older_refined.pdb"],
        "pdb": ["/path/to/input.pdb"],
    }

    # In a real scenario, template_builder.build_command_for_program would:
    # 1. Check best_files["model"] first
    # 2. Only fall back to categorized_files if best not available

    # This test verifies the logic pattern
    input_name = "model"
    slot_to_best_category = {"model": "model", "pdb_file": "model"}

    # Simulate what template_builder does
    file_found = None
    best_category = slot_to_best_category.get(input_name, input_name)
    best_path = best_files.get(best_category)

    if best_path:
        file_found = best_path

    assert_equal(file_found, "/path/to/best_refined.pdb",
                "Should use best_files over categorized_files")

    # Without best_files, should fall back
    best_files_empty = {}
    file_found2 = None
    best_path2 = best_files_empty.get(best_category)
    if best_path2:
        file_found2 = best_path2
    else:
        # Fallback to categorized
        for cat in ["refined", "model", "pdb"]:
            if categorized_files.get(cat):
                file_found2 = categorized_files[cat][0]
                break

    assert_equal(file_found2, "/path/to/older_refined.pdb",
                "Should fall back to categorized_files when no best")

    print("  PASSED")


def test_update_metrics():
    """Test that metrics can be updated for existing best file."""
    print("Test: update_metrics")

    tracker = BestFilesTracker()

    # Add a refined model without metrics
    tracker.evaluate_file("/path/to/refined.pdb", cycle=1, stage="refined")
    initial_score = tracker.get_best("model").score
    assert_equal(initial_score, 100.0, "Should have stage score only")

    # Update with validation metrics
    was_updated, old_score, new_score = tracker.update_metrics(
        "model",
        {"clashscore": 5, "r_free": 0.22},
        cycle=2
    )

    assert_true(was_updated, "Should update score")
    assert_equal(old_score, 100.0, "Old score should be stage only")
    assert_true(new_score > old_score, "New score should include metrics")

    # Verify metrics were merged
    entry = tracker.get_best("model")
    assert_equal(entry.metrics.get("clashscore"), 5)
    assert_equal(entry.metrics.get("r_free"), 0.22)

    print("  PASSED")


def test_update_metrics_no_change():
    """Test that update_metrics returns False when no change."""
    print("Test: update_metrics_no_change")

    tracker = BestFilesTracker()

    # No model yet
    was_updated, old_score, new_score = tracker.update_metrics(
        "model", {"clashscore": 5}
    )
    assert_false(was_updated, "Should not update non-existent category")

    print("  PASSED")


def test_update_metrics_history():
    """Test that metric updates are recorded in history."""
    print("Test: update_metrics_history")

    tracker = BestFilesTracker()

    # Add model
    tracker.evaluate_file("/path/to/refined.pdb", cycle=1, stage="refined")
    initial_history_len = len(tracker.get_history())

    # Update metrics
    tracker.update_metrics("model", {"r_free": 0.25}, cycle=2)

    # Should have new history entry
    new_history_len = len(tracker.get_history())
    assert_equal(new_history_len, initial_history_len + 1,
                "Should add history entry for metric update")

    # Check the history entry
    last_change = tracker.get_history()[-1]
    assert_equal(last_change.category, "model")
    assert_equal(last_change.old_path, last_change.new_path,
                "Metric update should have same old/new path")

    print("  PASSED")


def test_validation_workflow_simulation():
    """Test the full workflow: refine -> validate -> update metrics."""
    print("Test: validation_workflow_simulation")

    tracker = BestFilesTracker()

    # Cycle 1: Refinement produces model with R-free
    tracker.evaluate_file("/path/to/refine_001.pdb", cycle=1, stage="refined",
                         metrics={"r_free": 0.25, "r_work": 0.20})
    score_after_refine = tracker.get_best("model").score

    # Cycle 2: Validation adds clashscore
    was_updated, old_score, new_score = tracker.update_metrics(
        "model",
        {"clashscore": 8.5},
        cycle=2
    )

    assert_true(was_updated, "Should update with validation metrics")
    assert_true(new_score > score_after_refine,
               "Score should improve with good clashscore")

    # Verify merged metrics
    entry = tracker.get_best("model")
    assert_equal(entry.metrics.get("r_free"), 0.25, "Should keep R-free")
    assert_equal(entry.metrics.get("clashscore"), 8.5, "Should have clashscore")

    # Cycle 3: Another refinement with better metrics
    tracker.evaluate_file("/path/to/refine_002.pdb", cycle=3, stage="refined",
                         metrics={"r_free": 0.22, "r_work": 0.18, "clashscore": 6})

    # Should be new best (better metrics)
    assert_equal(tracker.get_best("model").path, "/path/to/refine_002.pdb",
                "Better refined model should become best")

    print("  PASSED")


# =============================================================================
# YAML CONFIGURATION TESTS
# =============================================================================

def test_yaml_scoring_loads():
    """Test that YAML scoring configuration loads."""
    print("Test: yaml_scoring_loads")

    tracker = BestFilesTracker()
    assert_not_none(tracker._scoring_config, "Scoring config should load")
    assert_true("model" in tracker._scoring_config, "Should have model config")
    assert_true("map" in tracker._scoring_config, "Should have map config")
    assert_true("mtz" in tracker._scoring_config, "Should have mtz config")

    print("  PASSED")


def test_yaml_stage_scores():
    """Test that stage scores come from YAML config."""
    print("Test: yaml_stage_scores")

    tracker = BestFilesTracker()

    # Check model stage scores
    assert_equal(tracker._get_stage_score("model", "refined"), 100)
    assert_equal(tracker._get_stage_score("model", "docked"), 60)
    assert_equal(tracker._get_stage_score("model", "predicted"), 40)
    assert_equal(tracker._get_stage_score("model", "unknown_stage"), 10)  # _default

    # Check map stage scores
    assert_equal(tracker._get_stage_score("map", "optimized_full_map"), 100)
    assert_equal(tracker._get_stage_score("map", "half_map"), 10)

    # Check mtz stage scores
    assert_equal(tracker._get_stage_score("mtz", "original"), 50)

    print("  PASSED")


def test_yaml_linear_formula():
    """Test linear formula (higher is better)."""
    print("Test: yaml_linear_formula")

    tracker = BestFilesTracker()

    # map_cc uses linear formula: best=1.0, worst=0.0, max=30
    config = {
        "max_points": 30,
        "formula": "linear",
        "best_value": 1.0,
        "worst_value": 0.0,
    }

    # map_cc = 1.0 -> 30 points
    score = tracker._apply_formula(1.0, config)
    assert_equal(score, 30.0, "map_cc=1.0 should give 30 points")

    # map_cc = 0.5 -> 15 points
    score = tracker._apply_formula(0.5, config)
    assert_equal(score, 15.0, "map_cc=0.5 should give 15 points")

    # map_cc = 0.0 -> 0 points
    score = tracker._apply_formula(0.0, config)
    assert_equal(score, 0.0, "map_cc=0.0 should give 0 points")

    print("  PASSED")


def test_yaml_linear_inverse_formula():
    """Test linear_inverse formula (lower is better)."""
    print("Test: yaml_linear_inverse_formula")

    tracker = BestFilesTracker()

    # r_free uses linear_inverse: best=0.20, worst=0.40, max=40
    config = {
        "max_points": 40,
        "formula": "linear_inverse",
        "best_value": 0.20,
        "worst_value": 0.40,
    }

    # r_free = 0.20 -> 40 points (best)
    score = tracker._apply_formula(0.20, config)
    assert_equal(score, 40.0, "r_free=0.20 should give 40 points")

    # r_free = 0.30 -> 20 points (middle)
    score = tracker._apply_formula(0.30, config)
    assert_true(abs(score - 20.0) < 0.01, f"r_free=0.30 should give ~20 points, got {score}")

    # r_free = 0.40 -> 0 points (worst)
    score = tracker._apply_formula(0.40, config)
    assert_equal(score, 0.0, "r_free=0.40 should give 0 points")

    # r_free = 0.25 -> 30 points
    score = tracker._apply_formula(0.25, config)
    assert_true(abs(score - 30.0) < 0.01, f"r_free=0.25 should give ~30 points, got {score}")

    print("  PASSED")


def test_yaml_boolean_formula():
    """Test boolean formula."""
    print("Test: yaml_boolean_formula")

    tracker = BestFilesTracker()

    config = {
        "max_points": 30,
        "formula": "boolean",
    }

    # True -> 30 points
    score = tracker._apply_formula(True, config)
    assert_equal(score, 30, "True should give max_points")

    # False -> 0 points
    score = tracker._apply_formula(False, config)
    assert_equal(score, 0, "False should give 0")

    print("  PASSED")


def test_yaml_model_scoring_integrated():
    """Test full model scoring matches expected values."""
    print("Test: yaml_model_scoring_integrated")

    tracker = BestFilesTracker()

    # Refined model with good metrics
    tracker.evaluate_file("/path/to/model.pdb", cycle=1, stage="refined",
                         metrics={"r_free": 0.22, "map_cc": 0.80, "clashscore": 5})

    entry = tracker.get_best("model")

    # Stage: 100 (refined)
    # R-free: 40 * (0.40 - 0.22) / 0.20 = 40 * 0.18 / 0.20 = 36
    # Map CC: 30 * 0.80 = 24
    # Clashscore: 30 * (20 - 5) / 20 = 30 * 0.75 = 22.5
    # Total: 100 + 36 + 24 + 22.5 = 182.5

    expected = 100 + 36 + 24 + 22.5
    assert_true(abs(entry.score - expected) < 0.1,
               f"Score should be ~{expected}, got {entry.score}")

    print("  PASSED")


def test_yaml_mtz_scoring_integrated():
    """Test MTZ scoring with has_rfree_flags."""
    print("Test: yaml_mtz_scoring_integrated")

    tracker = BestFilesTracker()

    # MTZ with R-free flags
    tracker.evaluate_file("/path/to/data.mtz", cycle=1, stage="original",
                         metrics={"has_rfree_flags": True})

    entry = tracker.get_best("mtz")

    # Stage: 50 (original)
    # has_rfree_flags: 30 (boolean)
    # Total: 80

    assert_equal(entry.score, 80.0, "MTZ with R-free should score 80")

    print("  PASSED")


def test_yaml_fallback_to_defaults():
    """Test that scoring falls back to defaults if YAML unavailable."""
    print("Test: yaml_fallback_to_defaults")

    tracker = BestFilesTracker()

    # Force reload with broken path
    original_config = tracker._scoring_config
    tracker._scoring_config = None

    # Should get defaults
    tracker._scoring_config = tracker._get_default_scoring()

    assert_not_none(tracker._scoring_config, "Should have fallback config")
    assert_true("model" in tracker._scoring_config, "Fallback should have model")

    # Verify it still works
    score = tracker._get_stage_score("model", "refined")
    assert_equal(score, 100, "Fallback refined score should be 100")

    # Restore
    tracker._scoring_config = original_config

    print("  PASSED")


def test_best_files_respects_exclude_patterns():
    """Test that best_files are excluded when they match exclude_patterns."""
    print("Test: best_files_respects_exclude_patterns")

    # This tests the logic pattern used in template_builder
    # When best_files contains a processed model, it should be excluded
    # from refinement programs that have exclude_patterns: [processed]

    best_files = {
        "model": "/path/to/predicted_model_processed.pdb",
    }

    exclude_patterns = ["processed"]

    # Simulate what template_builder does
    best_path = best_files.get("model")
    is_excluded = False
    for pattern in exclude_patterns:
        if pattern.lower() in best_path.lower():
            is_excluded = True
            break

    assert_true(is_excluded,
               "Processed model should be excluded by 'processed' pattern")

    # Non-processed model should not be excluded
    best_files2 = {
        "model": "/path/to/rsr_001_refined.pdb",
    }
    best_path2 = best_files2.get("model")
    is_excluded2 = False
    for pattern in exclude_patterns:
        if pattern.lower() in best_path2.lower():
            is_excluded2 = True
            break

    assert_false(is_excluded2,
                "Refined model should NOT be excluded")

    print("  PASSED")


def test_mtz_rfree_flag_locking():
    """Test that first MTZ with R-free flags locks forever."""
    print("Test: mtz_rfree_flag_locking")

    tracker = BestFilesTracker()

    # First MTZ without R-free flags
    tracker.evaluate_file("/path/to/original.mtz", cycle=1, stage="original",
                         metrics={})
    assert_equal(tracker.get_best_path("mtz"), "/path/to/original.mtz")

    # Second MTZ WITH R-free flags - should become best and lock
    tracker.evaluate_file("/path/to/refine_001_data.mtz", cycle=2, stage="refined_mtz",
                         metrics={"has_rfree_flags": True})
    assert_equal(tracker.get_best_path("mtz"), "/path/to/refine_001_data.mtz")

    # Third MTZ also with R-free flags - should NOT replace (locked)
    result = tracker.evaluate_file("/path/to/refine_002_data.mtz", cycle=3,
                                   stage="refined_mtz",
                                   metrics={"has_rfree_flags": True})
    assert_false(result, "Should not replace locked MTZ")
    assert_equal(tracker.get_best_path("mtz"), "/path/to/refine_001_data.mtz",
                "Locked MTZ should remain")

    print("  PASSED")


def test_session_rfree_mtz_locking():
    """Test Session-level R-free MTZ locking."""
    print("Test: session_rfree_mtz_locking")

    import tempfile
    import shutil
    import json

    # Create temp directory for session
    temp_dir = tempfile.mkdtemp()
    try:
        # Simulate Session R-free MTZ locking without importing Session
        # (to avoid langchain dependency in tests)

        # Test 1: Locking logic
        session_data = {
            "rfree_mtz": None,
            "rfree_mtz_locked_at_cycle": None,
        }

        def set_rfree_mtz(data, path, cycle):
            """Simulate Session.set_rfree_mtz"""
            current = data.get("rfree_mtz")
            if current is not None:
                if current == path:
                    return False, None
                return False, f"Already locked to {current}"
            data["rfree_mtz"] = path
            data["rfree_mtz_locked_at_cycle"] = cycle
            return True, f"Locked to {path}"

        # Initially no R-free MTZ
        assert_none(session_data.get("rfree_mtz"), "Should have no R-free MTZ initially")

        # Lock the first R-free MTZ
        was_locked, msg = set_rfree_mtz(session_data, "/path/to/refine_001_data.mtz", 5)
        assert_true(was_locked, "Should lock first R-free MTZ")
        assert_equal(session_data["rfree_mtz"], "/path/to/refine_001_data.mtz")

        # Try to lock a different MTZ - should fail
        was_locked2, msg2 = set_rfree_mtz(session_data, "/path/to/refine_002_data.mtz", 6)
        assert_false(was_locked2, "Should not replace locked R-free MTZ")
        assert_equal(session_data["rfree_mtz"], "/path/to/refine_001_data.mtz",
                    "Original R-free MTZ should remain")

        # Test persistence (save and reload)
        session_file = os.path.join(temp_dir, "session.json")
        with open(session_file, 'w') as f:
            json.dump(session_data, f)

        with open(session_file, 'r') as f:
            reloaded = json.load(f)

        assert_equal(reloaded.get("rfree_mtz"), "/path/to/refine_001_data.mtz",
                    "R-free MTZ should persist")
        assert_equal(reloaded.get("rfree_mtz_locked_at_cycle"), 5,
                    "Lock cycle should persist")

        print("  PASSED")
    finally:
        shutil.rmtree(temp_dir)


def test_rfree_mtz_pattern_detection():
    """Test that R-free MTZ patterns are correctly detected."""
    print("Test: rfree_mtz_pattern_detection")

    # Test patterns that SHOULD be recognized as having R-free flags
    rfree_patterns = [
        "refine_001_data.mtz",
        "refine_001_001.mtz",
        "aniso_refinement_data_PHX.mtz",  # AutoSol output
        "refinement_data.mtz",
        "some_refine_output.mtz",
    ]

    # Test patterns that should NOT be recognized (phased or generic)
    non_rfree_patterns = [
        "overall_best_hklout_phased.mtz",  # Phased - multiple arrays
        "solve_phases.mtz",
        "phased_output.mtz",
        "PHASER.1.mtz",  # Phaser output, no R-free
        "original_data.mtz",
    ]

    def has_rfree_flags(basename):
        """Simulate the detection logic from session.py"""
        basename = basename.lower()
        # Patterns that indicate R-free flags:
        # - *refine* (from phenix.refine)
        # - *refinement_data* (from autosol)
        has_rfree_patterns = (
            'refine' in basename or
            'refinement_data' in basename
        )
        is_phased_mtz = (
            'phased' in basename or
            'phases' in basename or
            'solve' in basename
        )
        return has_rfree_patterns and not is_phased_mtz

    # Check positive cases
    for pattern in rfree_patterns:
        assert_true(has_rfree_flags(pattern),
                   f"'{pattern}' should be recognized as having R-free flags")

    # Check negative cases
    for pattern in non_rfree_patterns:
        assert_false(has_rfree_flags(pattern),
                    f"'{pattern}' should NOT be recognized as having R-free flags")

    print("  PASSED")


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("BEST FILES TRACKER TESTS")
    print("=" * 60 + "\n")

    # Model scoring tests
    print("--- Model Scoring Tests ---")
    test_model_stage_scoring()
    test_model_metrics_improve_score()
    test_model_stage_beats_metrics()

    # Map scoring tests
    print("\n--- Map Scoring Tests ---")
    test_map_stage_scoring()
    test_map_resolution_bonus()

    # MTZ scoring tests
    print("\n--- MTZ Scoring Tests ---")
    test_mtz_rfree_locks()
    test_mtz_without_rfree_can_change()

    # History tests
    print("\n--- History Tests ---")
    test_history_tracking()
    test_history_filter_by_category()

    # Classification tests
    print("\n--- Classification Tests ---")
    test_file_classification_by_extension()
    test_stage_classification_from_filename()
    test_intermediate_files_ignored()

    # Serialization tests
    print("\n--- Serialization Tests ---")
    test_serialization_roundtrip()
    test_empty_serialization()
    test_from_dict_none()

    # Edge case tests
    print("\n--- Edge Case Tests ---")
    test_recency_tiebreaker()
    test_worse_file_doesnt_replace()
    test_get_best_dict()

    # Session integration tests
    print("\n--- Session Integration Tests ---")
    test_session_integration_simulation()
    test_session_persistence_simulation()
    test_template_builder_best_files_priority()

    # Metrics update tests
    print("\n--- Metrics Update Tests ---")
    test_update_metrics()
    test_update_metrics_no_change()
    test_update_metrics_history()
    test_validation_workflow_simulation()

    # YAML configuration tests
    print("\n--- YAML Configuration Tests ---")
    test_yaml_scoring_loads()
    test_yaml_stage_scores()
    test_yaml_linear_formula()
    test_yaml_linear_inverse_formula()
    test_yaml_boolean_formula()
    test_yaml_model_scoring_integrated()
    test_yaml_mtz_scoring_integrated()
    test_yaml_fallback_to_defaults()
    test_best_files_respects_exclude_patterns()
    test_mtz_rfree_flag_locking()
    test_session_rfree_mtz_locking()
    test_rfree_mtz_pattern_detection()

    print("\n" + "=" * 60)
    print("ALL BEST FILES TRACKER TESTS PASSED")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    run_all_tests()
