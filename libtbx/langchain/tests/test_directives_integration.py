"""
Integration tests for the complete directives system.

These tests verify that directives flow correctly through the entire system,
from extraction to validation to workflow state modification.

Run with: python tests/test_directives_integration.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os
import unittest

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agent.directive_extractor import (
    extract_directives_simple,
    validate_directives,
    get_program_settings,
    check_stop_conditions,
)
from agent.directive_validator import (
    validate_intent,
    augment_intent_with_directives,
    _check_stop_conditions as validator_check_stop,
)


class TestDirectivesEndToEnd(unittest.TestCase):
    """End-to-end tests for directive processing."""

    def test_resolution_flow_simple_extraction(self):
        """Test resolution directive from extraction to validation."""
        # Step 1: Extract directives from user advice
        user_advice = "use resolution 3.0 in autosol but 2.5 in other programs"

        # Using simple extraction (no LLM)
        directives = extract_directives_simple(user_advice)

        # Simple extraction only gets global resolution
        # For program-specific, we'd need LLM extraction
        # But we can test with manually constructed directives
        directives = {
            "program_settings": {
                "phenix.autosol": {"resolution": 3.0},
                "default": {"resolution": 2.5}
            }
        }

        # Step 2: Validate directives
        validated = validate_directives(directives)
        self.assertEqual(validated["program_settings"]["phenix.autosol"]["resolution"], 3.0)
        self.assertEqual(validated["program_settings"]["default"]["resolution"], 2.5)

        # Step 3: Get program settings
        autosol_settings = get_program_settings(validated, "phenix.autosol")
        refine_settings = get_program_settings(validated, "phenix.refine")

        self.assertEqual(autosol_settings["resolution"], 3.0)
        self.assertEqual(refine_settings["resolution"], 2.5)

        # Step 4: Validate LLM intent
        # Case A: LLM correctly uses autosol resolution
        intent_a = {
            "program": "phenix.autosol",
            "strategy": {"resolution": 3.0}
        }
        result_a = validate_intent(intent_a, validated, cycle_number=1)
        self.assertEqual(result_a["validated_intent"]["strategy"]["resolution"], 3.0)
        self.assertEqual(len(result_a["warnings"]), 0)

        # Case B: LLM forgets resolution for refine
        intent_b = {
            "program": "phenix.refine",
            "strategy": {}
        }
        result_b = validate_intent(intent_b, validated, cycle_number=2)
        self.assertEqual(result_b["validated_intent"]["strategy"]["resolution"], 2.5)
        self.assertTrue(len(result_b["modifications"]) > 0)

        # Case C: LLM uses wrong resolution for autosol
        intent_c = {
            "program": "phenix.autosol",
            "strategy": {"resolution": 2.5}  # Wrong!
        }
        result_c = validate_intent(intent_c, validated, cycle_number=1)
        self.assertEqual(result_c["validated_intent"]["strategy"]["resolution"], 3.0)  # Corrected
        self.assertTrue(len(result_c["warnings"]) > 0)

    def test_stop_after_first_refinement(self):
        """Test 'stop after first refinement' directive."""
        # Extract from simple patterns
        directives = extract_directives_simple("stop after the first refinement")

        self.assertIn("stop_conditions", directives)
        self.assertEqual(directives["stop_conditions"]["after_program"], "phenix.refine")
        self.assertEqual(directives["stop_conditions"]["max_refine_cycles"], 1)

        # Validate directive
        validated = validate_directives(directives)

        # Check stop conditions using validator function (which takes history)
        # Before refine runs
        history_before = [
            {"program": "phenix.phaser"}
        ]
        should_stop, reason = validator_check_stop(
            validated, cycle_number=2, last_program="phenix.phaser", history=history_before, log=lambda x: None
        )
        self.assertFalse(should_stop)

        # After refine completes
        history_after = [
            {"program": "phenix.phaser"},
            {"program": "phenix.refine"}
        ]
        should_stop, reason = validator_check_stop(
            validated, cycle_number=3, last_program="phenix.refine", history=history_after, log=lambda x: None
        )
        self.assertTrue(should_stop)
        self.assertIn("refine", reason.lower())

    def test_stop_at_cycle_n(self):
        """Test 'stop at cycle N' directive."""
        directives = extract_directives_simple("stop at cycle 4")

        self.assertIn("stop_conditions", directives)
        self.assertEqual(directives["stop_conditions"]["after_cycle"], 4)

        validated = validate_directives(directives)

        # Cycle 3 - should not stop
        should_stop, _ = check_stop_conditions(validated, cycle_number=3, last_program="phenix.refine")
        self.assertFalse(should_stop)

        # Cycle 4 - should stop
        should_stop, reason = check_stop_conditions(validated, cycle_number=4, last_program="phenix.refine")
        self.assertTrue(should_stop)
        self.assertIn("4", reason)

    def test_anisotropic_refinement(self):
        """Test 'run anisotropic refinement' directive."""
        directives = extract_directives_simple("run anisotropic refinement")

        self.assertIn("program_settings", directives)
        self.assertTrue(directives["program_settings"]["phenix.refine"]["anisotropic_adp"])

        validated = validate_directives(directives)

        # LLM intent without anisotropic
        intent = {
            "program": "phenix.refine",
            "strategy": {"resolution": 2.0}
        }

        result = validate_intent(intent, validated, cycle_number=1)

        # Should have anisotropic_adp added
        self.assertTrue(result["validated_intent"]["strategy"]["anisotropic_adp"])
        self.assertTrue(any("anisotropic" in m.lower() for m in result["modifications"]))

    def test_skip_validation_directive(self):
        """Test skip_validation directive allows stopping."""
        directives = extract_directives_simple("stop after the first refinement, skip validation")

        # Add skip_validation manually since simple extraction may not get it
        directives["stop_conditions"]["skip_validation"] = True

        validated = validate_directives(directives)

        # LLM wants to stop
        intent = {
            "program": None,
            "stop": True,
            "stop_reason": "User requested"
        }

        result = validate_intent(intent, validated, cycle_number=2, history=[])

        # validate_intent no longer returns should_stop (checked post-execution)
        # Just verify it returns the intent unchanged
        self.assertIn("validated_intent", result)

    def test_multiple_directives_combined(self):
        """Test multiple directives working together."""
        # Construct complex directives
        directives = {
            "program_settings": {
                "phenix.autosol": {"resolution": 3.0, "atom_type": "Se"},
                "phenix.refine": {"anisotropic_adp": True},
                "default": {"resolution": 2.5}
            },
            "stop_conditions": {
                "max_refine_cycles": 2,
                "r_free_target": 0.25
            },
            "workflow_preferences": {
                "skip_programs": ["phenix.autobuild"]
            },
            "constraints": [
                "User will do ligand fitting manually"
            ]
        }

        validated = validate_directives(directives)

        # Check autosol settings
        autosol_settings = get_program_settings(validated, "phenix.autosol")
        self.assertEqual(autosol_settings["resolution"], 3.0)
        self.assertEqual(autosol_settings["atom_type"], "Se")

        # Check refine settings (should merge default + specific)
        refine_settings = get_program_settings(validated, "phenix.refine")
        self.assertEqual(refine_settings["resolution"], 2.5)  # From default
        self.assertTrue(refine_settings["anisotropic_adp"])  # From specific

        # Check r_free target stop using validator function
        history = [
            {"program": "phenix.refine", "metrics": {"r_free": 0.24}}
        ]
        should_stop, reason = validator_check_stop(
            validated, cycle_number=3, last_program="phenix.refine", history=history, log=lambda x: None
        )
        self.assertTrue(should_stop)
        self.assertIn("R-free", reason)

    def test_empty_directives_passthrough(self):
        """Test that empty directives don't affect normal operation."""
        intent = {
            "program": "phenix.refine",
            "strategy": {"resolution": 2.5},
            "files": {"mtz": "data.mtz"}
        }

        # Empty directives
        result = validate_intent(intent, {}, cycle_number=1)

        # Intent should pass through unchanged
        self.assertEqual(result["validated_intent"]["program"], "phenix.refine")
        self.assertEqual(result["validated_intent"]["strategy"]["resolution"], 2.5)
        self.assertEqual(result["modifications"], [])
        self.assertEqual(result["warnings"], [])
        # NOTE: should_stop no longer returned from validate_intent
        # Stop conditions are checked post-execution in ai_agent.py

    def test_directive_validation_with_history(self):
        """Test directive validation applies program settings correctly."""
        # NOTE: Stop condition checking has been moved to post-execution.
        # This test now verifies that validate_intent still applies settings.
        directives = {
            "program_settings": {
                "phenix.refine": {"resolution": 2.0}
            }
        }

        validated = validate_directives(directives)

        # History doesn't affect program settings
        history = [
            {"program": "phenix.phaser"},
            {"program": "phenix.refine"},
        ]

        intent = {"program": "phenix.refine", "strategy": {}}
        result = validate_intent(intent, validated, cycle_number=3, history=history)

        # Should apply the resolution setting
        self.assertEqual(result["validated_intent"]["strategy"]["resolution"], 2.0)


class TestDirectiveEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def test_invalid_resolution_ignored(self):
        """Invalid resolution values should be ignored."""
        directives = {
            "program_settings": {
                "default": {"resolution": "invalid"}
            }
        }

        # Should not crash, just ignore invalid value
        validated = validate_directives(directives)
        # The invalid value should be converted or ignored
        # Based on implementation, it tries float() which will fail
        self.assertNotIn("program_settings", validated)

    def test_unknown_program_fixed(self):
        """Unknown program names should be fixed if possible."""
        directives = {
            "program_settings": {
                "refine": {"resolution": 2.5}  # Missing "phenix." prefix
            }
        }

        validated = validate_directives(directives)

        # Should be fixed to phenix.refine
        self.assertIn("phenix.refine", validated.get("program_settings", {}))

    def test_contradictory_stop_conditions(self):
        """Test handling of potentially contradictory stop conditions."""
        directives = {
            "stop_conditions": {
                "after_cycle": 3,
                "after_program": "phenix.refine",
                "r_free_target": 0.20
            }
        }

        validated = validate_directives(directives)

        # At cycle 3 with refine done and r_free not at target
        history = [
            {"program": "phenix.phaser"},
            {"program": "phenix.refine", "metrics": {"r_free": 0.30}},
            {"program": "phenix.refine", "metrics": {"r_free": 0.28}},
        ]

        # after_cycle should trigger first (using simple check_stop_conditions)
        should_stop, reason = check_stop_conditions(
            validated, cycle_number=3, last_program="phenix.refine"
        )
        self.assertTrue(should_stop)

    def test_none_values_handled(self):
        """Test that None values don't cause crashes."""
        intent = {
            "program": "phenix.refine",
            "strategy": None,  # Could be None
            "files": None
        }

        directives = {
            "program_settings": {
                "default": {"resolution": 2.5}
            }
        }

        # Should not crash
        result = validate_intent(intent, directives, cycle_number=1)
        self.assertIsNotNone(result["validated_intent"])
        # Resolution should still be added
        self.assertEqual(result["validated_intent"]["strategy"]["resolution"], 2.5)

    def test_augment_preserves_existing(self):
        """augment_intent_with_directives should not override existing values."""
        intent = {
            "program": "phenix.refine",
            "strategy": {"resolution": 3.0}  # Already set
        }

        directives = {
            "program_settings": {
                "default": {"resolution": 2.5, "cycles": 5}
            }
        }

        result = augment_intent_with_directives(intent, directives)

        # Resolution should NOT be changed (already set)
        self.assertEqual(result["strategy"]["resolution"], 3.0)
        # Cycles should be added (was missing)
        self.assertEqual(result["strategy"]["cycles"], 5)


class TestWorkflowPreferences(unittest.TestCase):
    """Test workflow preference directives."""

    def test_skip_programs_in_directives(self):
        """Test that skip_programs is properly validated."""
        directives = {
            "workflow_preferences": {
                "skip_programs": ["phenix.autobuild", "phenix.ligandfit"]
            }
        }

        validated = validate_directives(directives)

        self.assertIn("workflow_preferences", validated)
        self.assertEqual(len(validated["workflow_preferences"]["skip_programs"]), 2)
        self.assertIn("phenix.autobuild", validated["workflow_preferences"]["skip_programs"])

    def test_prefer_programs_in_directives(self):
        """Test that prefer_programs is properly validated."""
        directives = {
            "workflow_preferences": {
                "prefer_programs": ["phenix.phaser", "phenix.refine"]
            }
        }

        validated = validate_directives(directives)

        self.assertIn("workflow_preferences", validated)
        self.assertEqual(len(validated["workflow_preferences"]["prefer_programs"]), 2)

    def test_experimental_phasing_preference(self):
        """Test experimental phasing preference."""
        directives = {
            "workflow_preferences": {
                "use_experimental_phasing": True
            }
        }

        validated = validate_directives(directives)

        self.assertTrue(validated["workflow_preferences"]["use_experimental_phasing"])


def run_all_tests():
    """Run all tests and raise exception if any fail."""
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromModule(sys.modules[__name__])
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    if not result.wasSuccessful():
        raise Exception("Some tests failed")


if __name__ == "__main__":
    unittest.main()
