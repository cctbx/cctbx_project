"""
Tests for the clean decision flow architecture.

These tests verify that:
1. Directives properly modify valid_programs (workflow_engine)
2. Validation gate is a simple "in list?" check (graph_nodes)
3. Stop conditions are checked only post-execution (ai_agent)

Run with: python tests/test_decision_flow.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os
import unittest

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Try to import workflow_engine - may fail without libtbx
try:
    from agent.workflow_engine import WorkflowEngine
    HAS_WORKFLOW_ENGINE = True
except ImportError:
    HAS_WORKFLOW_ENGINE = False
    WorkflowEngine = None

from agent.directive_validator import validate_intent
from agent.directive_extractor import check_stop_conditions


@unittest.skipUnless(HAS_WORKFLOW_ENGINE, "Requires libtbx for WorkflowEngine")
class TestWorkflowEngineDirectives(unittest.TestCase):
    """Test that workflow_engine._apply_directives correctly modifies valid_programs."""

    def setUp(self):
        self.engine = WorkflowEngine()

    def test_skip_validation_adds_stop(self):
        """skip_validation=true should add STOP to valid_programs."""
        valid_programs = ["phenix.refine", "phenix.molprobity"]
        directives = {
            "stop_conditions": {"skip_validation": True}
        }

        result = self.engine._apply_directives(valid_programs, directives, "refine")

        self.assertIn("STOP", result)

    def test_skip_validation_no_duplicate_stop(self):
        """Should not add duplicate STOP if already present."""
        valid_programs = ["phenix.refine", "STOP"]
        directives = {
            "stop_conditions": {"skip_validation": True}
        }

        result = self.engine._apply_directives(valid_programs, directives, "refine")

        self.assertEqual(result.count("STOP"), 1)

    def test_after_program_adds_target(self):
        """after_program should add the target program to valid_programs."""
        valid_programs = ["phenix.predict_and_build"]
        directives = {
            "stop_conditions": {
                "after_program": "phenix.resolve_cryo_em",
                "skip_validation": True
            }
        }

        result = self.engine._apply_directives(valid_programs, directives, "obtain_model")

        self.assertIn("phenix.resolve_cryo_em", result)
        # Should be at front (high priority)
        self.assertEqual(result[0], "phenix.resolve_cryo_em")

    def test_after_program_no_duplicate(self):
        """Should not add duplicate if target already in list."""
        valid_programs = ["phenix.resolve_cryo_em", "phenix.predict_and_build"]
        directives = {
            "stop_conditions": {
                "after_program": "phenix.resolve_cryo_em"
            }
        }

        result = self.engine._apply_directives(valid_programs, directives, "obtain_model")

        self.assertEqual(result.count("phenix.resolve_cryo_em"), 1)

    def test_skip_programs(self):
        """skip_programs should remove programs from list."""
        valid_programs = ["phenix.refine", "phenix.molprobity", "phenix.autobuild"]
        directives = {
            "workflow_preferences": {
                "skip_programs": ["phenix.autobuild"]
            }
        }

        result = self.engine._apply_directives(valid_programs, directives, "refine")

        self.assertNotIn("phenix.autobuild", result)
        self.assertIn("phenix.refine", result)
        self.assertIn("phenix.molprobity", result)

    def test_prefer_programs(self):
        """prefer_programs should move programs to front."""
        valid_programs = ["phenix.refine", "phenix.molprobity", "phenix.autobuild"]
        directives = {
            "workflow_preferences": {
                "prefer_programs": ["phenix.autobuild"]
            }
        }

        result = self.engine._apply_directives(valid_programs, directives, "refine")

        self.assertEqual(result[0], "phenix.autobuild")

    def test_empty_directives(self):
        """Empty directives should return list unchanged."""
        valid_programs = ["phenix.refine", "phenix.molprobity"]

        result = self.engine._apply_directives(valid_programs, {}, "refine")

        self.assertEqual(result, valid_programs)

    def test_none_directives(self):
        """None directives should return list unchanged."""
        valid_programs = ["phenix.refine", "phenix.molprobity"]

        result = self.engine._apply_directives(valid_programs, None, "refine")

        self.assertEqual(result, valid_programs)


class TestValidateIntentSimplified(unittest.TestCase):
    """Test that validate_intent no longer handles stop conditions."""

    def test_no_should_stop_in_result(self):
        """validate_intent should not return should_stop anymore."""
        intent = {"program": "phenix.refine", "strategy": {}}
        directives = {
            "stop_conditions": {"after_cycle": 1}  # Would have triggered stop before
        }

        result = validate_intent(intent, directives, cycle_number=1)

        # should_stop should not be in the result
        self.assertNotIn("should_stop", result)
        self.assertNotIn("stop_reason", result)

    def test_still_applies_program_settings(self):
        """validate_intent should still apply program settings."""
        intent = {"program": "phenix.refine", "strategy": {}}
        directives = {
            "program_settings": {
                "default": {"resolution": 2.5}
            }
        }

        result = validate_intent(intent, directives, cycle_number=1)

        self.assertEqual(result["validated_intent"]["strategy"]["resolution"], 2.5)

    def test_returns_modifications(self):
        """validate_intent should still return modifications list."""
        intent = {"program": "phenix.refine", "strategy": {}}
        directives = {
            "program_settings": {
                "default": {"resolution": 2.5}
            }
        }

        result = validate_intent(intent, directives, cycle_number=1)

        self.assertIn("modifications", result)
        self.assertTrue(len(result["modifications"]) > 0)


class TestPostExecutionStopCheck(unittest.TestCase):
    """Test that stop conditions are checked correctly post-execution."""

    def test_after_program_triggers_stop(self):
        """after_program should trigger stop after program completes."""
        directives = {
            "stop_conditions": {"after_program": "phenix.resolve_cryo_em"}
        }

        # Program just completed
        should_stop, reason = check_stop_conditions(
            directives,
            cycle_number=2,
            last_program="phenix.resolve_cryo_em"
        )

        self.assertTrue(should_stop)
        self.assertIn("resolve_cryo_em", reason)

    def test_after_program_normalized_match(self):
        """Should match with or without phenix. prefix."""
        directives = {
            "stop_conditions": {"after_program": "phenix.xtriage"}
        }

        # Program without prefix
        should_stop, reason = check_stop_conditions(
            directives,
            cycle_number=1,
            last_program="xtriage"
        )

        self.assertTrue(should_stop)

    def test_after_program_no_match(self):
        """Should not stop if different program ran."""
        directives = {
            "stop_conditions": {"after_program": "phenix.resolve_cryo_em"}
        }

        should_stop, reason = check_stop_conditions(
            directives,
            cycle_number=2,
            last_program="phenix.mtriage"
        )

        self.assertFalse(should_stop)

    def test_after_cycle_triggers_stop(self):
        """after_cycle should trigger stop at specified cycle."""
        directives = {
            "stop_conditions": {"after_cycle": 3}
        }

        # At cycle 3
        should_stop, reason = check_stop_conditions(
            directives,
            cycle_number=3,
            last_program="phenix.refine"
        )

        self.assertTrue(should_stop)
        self.assertIn("cycle", reason.lower())

    def test_after_cycle_before_target(self):
        """Should not stop before target cycle."""
        directives = {
            "stop_conditions": {"after_cycle": 3}
        }

        should_stop, reason = check_stop_conditions(
            directives,
            cycle_number=2,
            last_program="phenix.refine"
        )

        self.assertFalse(should_stop)

    def test_empty_directives_no_stop(self):
        """Empty directives should not trigger stop."""
        should_stop, reason = check_stop_conditions(
            {},
            cycle_number=5,
            last_program="phenix.refine"
        )

        self.assertFalse(should_stop)

    def test_no_stop_conditions_no_stop(self):
        """Directives without stop_conditions should not trigger stop."""
        directives = {
            "program_settings": {"default": {"resolution": 2.5}}
        }

        should_stop, reason = check_stop_conditions(
            directives,
            cycle_number=5,
            last_program="phenix.refine"
        )

        self.assertFalse(should_stop)


@unittest.skipUnless(HAS_WORKFLOW_ENGINE, "Requires libtbx for WorkflowEngine")
class TestDecisionFlowIntegration(unittest.TestCase):
    """Integration tests for the complete decision flow."""

    def setUp(self):
        self.engine = WorkflowEngine()

    def test_tutorial_flow_resolve_cryo_em(self):
        """
        Tutorial: "Run density modification and stop"
        Should add resolve_cryo_em to valid_programs and allow stop after.
        """
        # Initial valid programs (from workflow phase - might not include resolve_cryo_em)
        valid_programs = ["phenix.predict_and_build"]

        # User directive for tutorial
        directives = {
            "stop_conditions": {
                "after_program": "phenix.resolve_cryo_em",
                "skip_validation": True
            }
        }

        # Step 1: Workflow engine adds required program
        result = self.engine._apply_directives(valid_programs, directives, "obtain_model")
        self.assertIn("phenix.resolve_cryo_em", result)
        self.assertIn("STOP", result)

        # Step 2: LLM can now choose resolve_cryo_em (it's in valid_programs)
        # (graph_nodes validation would pass)

        # Step 3: After resolve_cryo_em completes, check stop condition
        should_stop, reason = check_stop_conditions(
            directives,
            cycle_number=2,
            last_program="phenix.resolve_cryo_em"
        )
        self.assertTrue(should_stop)

    def test_normal_workflow_no_directives(self):
        """Normal workflow without directives should work unchanged."""
        valid_programs = ["phenix.refine", "phenix.molprobity"]

        # No directives
        result = self.engine._apply_directives(valid_programs, {}, "refine")

        # Should be unchanged
        self.assertEqual(result, valid_programs)
        self.assertNotIn("STOP", result)

        # No stop condition triggered
        should_stop, _ = check_stop_conditions({}, 1, "phenix.refine")
        self.assertFalse(should_stop)


if __name__ == "__main__":
    unittest.main()
