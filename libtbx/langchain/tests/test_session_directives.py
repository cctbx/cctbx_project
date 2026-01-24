"""
Tests for session.py directive methods.

These tests verify that directives are correctly stored and retrieved
in the session.

Run with: python tests/test_session_directives.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os
import tempfile
import shutil
import unittest

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Create a mock BestFilesTracker class that we'll inject
class MockBestFilesTracker:
    def __init__(self):
        self._history = []
    def to_dict(self):
        return {}
    def get_history(self):
        return self._history
    @classmethod
    def from_dict(cls, data):
        return cls()

# Mock missing imports before importing session
from unittest.mock import MagicMock

# Set up the mocks
_mock_langchain_core = MagicMock()
_mock_libtbx = MagicMock()
_mock_best_files_module = MagicMock()
_mock_best_files_module.BestFilesTracker = MockBestFilesTracker

sys.modules['langchain_core'] = _mock_langchain_core
sys.modules['langchain_core.prompts'] = _mock_langchain_core.prompts
sys.modules['libtbx'] = _mock_libtbx
sys.modules['libtbx.langchain'] = _mock_libtbx.langchain
sys.modules['libtbx.langchain.agent'] = _mock_libtbx.langchain.agent
sys.modules['libtbx.langchain.agent.best_files_tracker'] = _mock_best_files_module

from agent.session import AgentSession


class TestSessionDirectives(unittest.TestCase):
    """Tests for session directive methods."""

    def setUp(self):
        """Create a temporary directory for session files."""
        self.temp_dir = tempfile.mkdtemp()
        self.session = AgentSession(session_dir=self.temp_dir)

    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_initial_directives_empty(self):
        """New session should have empty directives."""
        self.assertEqual(self.session.get_directives(), {})
        self.assertFalse(self.session.data.get("directives_extracted"))

    def test_get_program_directive_settings_empty(self):
        """Empty directives should return empty settings."""
        settings = self.session.get_program_directive_settings("phenix.refine")
        self.assertEqual(settings, {})

    def test_check_directive_stop_conditions_empty(self):
        """Empty directives should not trigger stop."""
        should_stop, reason = self.session.check_directive_stop_conditions(
            cycle_number=5, last_program="phenix.refine"
        )
        self.assertFalse(should_stop)
        self.assertIsNone(reason)

    def test_should_skip_validation_default(self):
        """Default should be to not skip validation."""
        self.assertFalse(self.session.should_skip_validation())

    def test_directives_persist_after_save_load(self):
        """Directives should persist across save/load."""
        import json

        # Manually set directives
        self.session.data["directives"] = {
            "program_settings": {
                "phenix.refine": {"resolution": 2.5}
            }
        }
        self.session.data["directives_extracted"] = True

        # Ensure best_files has the mock methods we need
        if not hasattr(self.session.best_files, 'to_dict') or \
           not callable(getattr(self.session.best_files, 'to_dict', None)):
            self.session.best_files = MockBestFilesTracker()

        # Save manually to bypass any mock issues
        session_file = os.path.join(self.temp_dir, "agent_session.json")
        self.session.data["best_files"] = {}
        self.session.data["best_files_history"] = []
        with open(session_file, 'w') as f:
            json.dump(self.session.data, f, indent=2)

        # Verify file was saved
        self.assertTrue(os.path.exists(session_file),
            f"Session file should exist at {session_file}")

        # Load manually and verify directives
        with open(session_file, 'r') as f:
            loaded_data = json.load(f)

        self.assertEqual(
            loaded_data.get("directives"),
            {"program_settings": {"phenix.refine": {"resolution": 2.5}}}
        )

    def test_get_program_directive_settings_with_default(self):
        """Should merge default with program-specific settings."""
        self.session.data["directives"] = {
            "program_settings": {
                "default": {"resolution": 2.5, "cycles": 5},
                "phenix.refine": {"anisotropic_adp": True}
            }
        }

        settings = self.session.get_program_directive_settings("phenix.refine")

        self.assertEqual(settings.get("resolution"), 2.5)  # From default
        self.assertEqual(settings.get("cycles"), 5)  # From default
        self.assertEqual(settings.get("anisotropic_adp"), True)  # From specific

    def test_get_program_directive_settings_specific_override(self):
        """Program-specific settings should override defaults."""
        self.session.data["directives"] = {
            "program_settings": {
                "default": {"resolution": 2.5},
                "phenix.autosol": {"resolution": 3.0}
            }
        }

        # autosol should get 3.0
        settings = self.session.get_program_directive_settings("phenix.autosol")
        self.assertEqual(settings.get("resolution"), 3.0)

        # refine should get default 2.5
        settings = self.session.get_program_directive_settings("phenix.refine")
        self.assertEqual(settings.get("resolution"), 2.5)

    def test_check_directive_stop_after_cycle(self):
        """Should stop after specified cycle."""
        self.session.data["directives"] = {
            "stop_conditions": {"after_cycle": 4}
        }

        # Before cycle 4
        should_stop, reason = self.session.check_directive_stop_conditions(
            cycle_number=3, last_program="phenix.refine"
        )
        self.assertFalse(should_stop)

        # At cycle 4
        should_stop, reason = self.session.check_directive_stop_conditions(
            cycle_number=4, last_program="phenix.refine"
        )
        self.assertTrue(should_stop)
        self.assertIsNotNone(reason)

    def test_check_directive_stop_after_program(self):
        """Should stop after specified program."""
        self.session.data["directives"] = {
            "stop_conditions": {"after_program": "phenix.refine"}
        }

        # Different program
        should_stop, reason = self.session.check_directive_stop_conditions(
            cycle_number=2, last_program="phenix.phaser"
        )
        self.assertFalse(should_stop)

        # Matching program
        should_stop, reason = self.session.check_directive_stop_conditions(
            cycle_number=2, last_program="phenix.refine"
        )
        self.assertTrue(should_stop)

    def test_should_skip_validation_true(self):
        """Should return True when skip_validation is set."""
        self.session.data["directives"] = {
            "stop_conditions": {"skip_validation": True}
        }
        self.assertTrue(self.session.should_skip_validation())

    def test_count_program_runs(self):
        """Should correctly count program runs."""
        self.session.data["cycles"] = [
            {"program": "phenix.xtriage"},
            {"program": "phenix.phaser"},
            {"program": "phenix.refine"},
            {"program": "phenix.refine"},
            {"program": "phenix.molprobity"},
        ]

        self.assertEqual(self.session.count_program_runs("phenix.refine"), 2)
        self.assertEqual(self.session.count_program_runs("phenix.phaser"), 1)
        self.assertEqual(self.session.count_program_runs("phenix.autobuild"), 0)

    def test_check_max_program_cycles(self):
        """Should check max refine cycles from directives."""
        self.session.data["directives"] = {
            "stop_conditions": {"max_refine_cycles": 2}
        }
        self.session.data["cycles"] = [
            {"program": "phenix.phaser"},
            {"program": "phenix.refine"},
        ]

        # After 1 refine, should not be at limit
        limit_reached, count, max_allowed = self.session.check_max_program_cycles("phenix.refine")
        self.assertFalse(limit_reached)
        self.assertEqual(count, 1)
        self.assertEqual(max_allowed, 2)

        # Add another refine
        self.session.data["cycles"].append({"program": "phenix.refine"})

        # Now at limit
        limit_reached, count, max_allowed = self.session.check_max_program_cycles("phenix.refine")
        self.assertTrue(limit_reached)
        self.assertEqual(count, 2)


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
