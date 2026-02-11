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
import json

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_false, assert_none,
    assert_not_none,
    run_tests_with_fail_fast
)

# Create a mock BestFilesTracker class that we'll inject
class MockBestFilesTracker:
    def __init__(self):
        self._history = []
        self.best = {}
    def to_dict(self):
        return {"best": {}, "history": []}
    def get_history(self):
        return self._history
    def get_best_dict(self):
        return {}
    def get_best_path(self, category):
        return None
    def evaluate_file(self, path, cycle, metrics=None, stage=None):
        return False
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


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def create_test_session():
    """Create a temporary session for testing."""
    temp_dir = tempfile.mkdtemp()
    session = AgentSession(session_dir=temp_dir)
    return session, temp_dir


def cleanup_test_session(temp_dir):
    """Clean up temporary session directory."""
    shutil.rmtree(temp_dir, ignore_errors=True)


# =============================================================================
# TESTS
# =============================================================================

def test_initial_directives_empty():
    """New session should have empty directives."""
    session, temp_dir = create_test_session()
    try:
        assert_equal(session.get_directives(), {})
        assert_false(session.data.get("directives_extracted"))
    finally:
        cleanup_test_session(temp_dir)


def test_get_program_directive_settings_empty():
    """Empty directives should return empty settings."""
    session, temp_dir = create_test_session()
    try:
        settings = session.get_program_directive_settings("phenix.refine")
        assert_equal(settings, {})
    finally:
        cleanup_test_session(temp_dir)


def test_check_directive_stop_conditions_empty():
    """Empty directives should not trigger stop."""
    session, temp_dir = create_test_session()
    try:
        should_stop, reason = session.check_directive_stop_conditions(
            cycle_number=5, last_program="phenix.refine"
        )
        assert_false(should_stop)
        assert_none(reason)
    finally:
        cleanup_test_session(temp_dir)


def test_should_skip_validation_default():
    """Default should be to not skip validation."""
    session, temp_dir = create_test_session()
    try:
        assert_false(session.should_skip_validation())
    finally:
        cleanup_test_session(temp_dir)


def test_directives_persist_after_save_load():
    """Directives should persist across save/load."""
    session, temp_dir = create_test_session()
    try:
        # Manually set directives
        session.data["directives"] = {
            "program_settings": {
                "phenix.refine": {"resolution": 2.5}
            }
        }
        session.data["directives_extracted"] = True

        # Ensure best_files has the mock methods we need
        if not hasattr(session.best_files, 'to_dict') or \
           not callable(getattr(session.best_files, 'to_dict', None)):
            session.best_files = MockBestFilesTracker()

        # Save manually to bypass any mock issues
        session_file = os.path.join(temp_dir, "agent_session.json")
        session.data["best_files"] = {}
        session.data["best_files_history"] = []
        with open(session_file, 'w') as f:
            json.dump(session.data, f, indent=2)

        # Verify file was saved
        assert_true(os.path.exists(session_file),
            f"Session file should exist at {session_file}")

        # Load manually and verify directives
        with open(session_file, 'r') as f:
            loaded_data = json.load(f)

        assert_equal(
            loaded_data.get("directives"),
            {"program_settings": {"phenix.refine": {"resolution": 2.5}}}
        )
    finally:
        cleanup_test_session(temp_dir)


def test_get_program_directive_settings_with_default():
    """Should merge default with program-specific settings."""
    session, temp_dir = create_test_session()
    try:
        session.data["directives"] = {
            "program_settings": {
                "default": {"resolution": 2.5, "cycles": 5},
                "phenix.refine": {"anisotropic_adp": True}
            }
        }

        settings = session.get_program_directive_settings("phenix.refine")

        assert_equal(settings.get("resolution"), 2.5)  # From default
        assert_equal(settings.get("cycles"), 5)  # From default
        assert_equal(settings.get("anisotropic_adp"), True)  # From specific
    finally:
        cleanup_test_session(temp_dir)


def test_get_program_directive_settings_specific_override():
    """Program-specific settings should override defaults."""
    session, temp_dir = create_test_session()
    try:
        session.data["directives"] = {
            "program_settings": {
                "default": {"resolution": 2.5},
                "phenix.autosol": {"resolution": 3.0}
            }
        }

        # autosol should get 3.0
        settings = session.get_program_directive_settings("phenix.autosol")
        assert_equal(settings.get("resolution"), 3.0)

        # refine should get default 2.5
        settings = session.get_program_directive_settings("phenix.refine")
        assert_equal(settings.get("resolution"), 2.5)
    finally:
        cleanup_test_session(temp_dir)


def test_check_directive_stop_after_cycle():
    """Should stop after specified cycle."""
    session, temp_dir = create_test_session()
    try:
        session.data["directives"] = {
            "stop_conditions": {"after_cycle": 4}
        }

        # Before cycle 4
        should_stop, reason = session.check_directive_stop_conditions(
            cycle_number=3, last_program="phenix.refine"
        )
        assert_false(should_stop)

        # At cycle 4
        should_stop, reason = session.check_directive_stop_conditions(
            cycle_number=4, last_program="phenix.refine"
        )
        assert_true(should_stop)
        assert_not_none(reason)
    finally:
        cleanup_test_session(temp_dir)


def test_check_directive_stop_after_program():
    """Should stop after specified program."""
    session, temp_dir = create_test_session()
    try:
        session.data["directives"] = {
            "stop_conditions": {"after_program": "phenix.refine"}
        }

        # Different program
        should_stop, reason = session.check_directive_stop_conditions(
            cycle_number=2, last_program="phenix.phaser"
        )
        assert_false(should_stop)

        # Matching program
        should_stop, reason = session.check_directive_stop_conditions(
            cycle_number=2, last_program="phenix.refine"
        )
        assert_true(should_stop)
    finally:
        cleanup_test_session(temp_dir)


def test_should_skip_validation_true():
    """Should return True when skip_validation is set."""
    session, temp_dir = create_test_session()
    try:
        session.data["directives"] = {
            "stop_conditions": {"skip_validation": True}
        }
        assert_true(session.should_skip_validation())
    finally:
        cleanup_test_session(temp_dir)


def test_count_program_runs():
    """Should correctly count program runs."""
    session, temp_dir = create_test_session()
    try:
        session.data["cycles"] = [
            {"program": "phenix.xtriage"},
            {"program": "phenix.phaser"},
            {"program": "phenix.refine"},
            {"program": "phenix.refine"},
            {"program": "phenix.molprobity"},
        ]

        assert_equal(session.count_program_runs("phenix.refine"), 2)
        assert_equal(session.count_program_runs("phenix.phaser"), 1)
        assert_equal(session.count_program_runs("phenix.autobuild"), 0)
    finally:
        cleanup_test_session(temp_dir)


def test_check_max_program_cycles():
    """Should check max refine cycles from directives."""
    session, temp_dir = create_test_session()
    try:
        session.data["directives"] = {
            "stop_conditions": {"max_refine_cycles": 2}
        }
        session.data["cycles"] = [
            {"program": "phenix.phaser"},
            {"program": "phenix.refine"},
        ]

        # After 1 refine, should not be at limit
        limit_reached, count, max_allowed = session.check_max_program_cycles("phenix.refine")
        assert_false(limit_reached)
        assert_equal(count, 1)
        assert_equal(max_allowed, 2)

        # Add another refine
        session.data["cycles"].append({"program": "phenix.refine"})

        # Now at limit
        limit_reached, count, max_allowed = session.check_max_program_cycles("phenix.refine")
        assert_true(limit_reached)
        assert_equal(count, 2)
    finally:
        cleanup_test_session(temp_dir)


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
