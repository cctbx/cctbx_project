"""
Tests for event logging system.
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agent.event_log import (
    EventType, Verbosity, EventLog,
    get_event_log, emit_event, migrate_debug_log,
    verbosity_index,
)
from agent.event_formatter import EventFormatter


def test_event_type_constants():
    """Test that EventType constants are defined."""
    print("Test: EventType constants")

    required_types = [
        "CYCLE_START", "CYCLE_COMPLETE", "STATE_DETECTED",
        "METRICS_EXTRACTED", "METRICS_TREND", "SANITY_CHECK",
        "PROGRAM_SELECTED", "PROGRAM_MODIFIED", "STOP_DECISION",
        "FILES_SELECTED", "COMMAND_BUILT", "ERROR", "WARNING",
        "DEBUG", "THOUGHT",
    ]

    for type_name in required_types:
        assert hasattr(EventType, type_name), f"Missing EventType.{type_name}"
        value = getattr(EventType, type_name)
        assert isinstance(value, str), f"EventType.{type_name} should be string"

    print("  PASSED")


def test_verbosity_constants():
    """Test that Verbosity constants are defined."""
    print("Test: Verbosity constants")

    # There are 3 verbosity levels (no DEBUG level - DEBUG events use VERBOSE)
    assert Verbosity.QUIET == "quiet"
    assert Verbosity.NORMAL == "normal"
    assert Verbosity.VERBOSE == "verbose"

    # Test ordering
    assert verbosity_index(Verbosity.QUIET) < verbosity_index(Verbosity.NORMAL)
    assert verbosity_index(Verbosity.NORMAL) < verbosity_index(Verbosity.VERBOSE)

    print("  PASSED")


def test_event_log_emit():
    """Test EventLog.emit()."""
    print("Test: EventLog.emit()")

    events = []
    log = EventLog(events)

    event = log.emit(EventType.STATE_DETECTED, workflow_state="test", foo="bar")

    assert len(events) == 1
    assert events[0]["type"] == EventType.STATE_DETECTED
    assert events[0]["workflow_state"] == "test"
    assert events[0]["foo"] == "bar"
    assert event is events[0]

    print("  PASSED")


def test_event_log_convenience_methods():
    """Test EventLog convenience methods."""
    print("Test: EventLog convenience methods")

    events = []
    log = EventLog(events)

    log.emit_debug("Debug message")
    log.emit_warning("Warning message")
    log.emit_error("Error message", details="Details here")
    log.emit_thought("Chain of thought")

    assert len(events) == 4
    assert events[0]["type"] == EventType.DEBUG
    assert events[1]["type"] == EventType.WARNING
    assert events[2]["type"] == EventType.ERROR
    assert events[3]["type"] == EventType.THOUGHT

    assert log.has_errors()
    assert log.has_warnings()

    print("  PASSED")


def test_event_log_filtering():
    """Test EventLog.get_events() filtering."""
    print("Test: EventLog filtering")

    events = []
    log = EventLog(events)

    log.emit(EventType.CYCLE_START, cycle=1)  # QUIET level
    log.emit(EventType.STATE_DETECTED, state="test")  # NORMAL level
    log.emit(EventType.FILES_SELECTED, files=[])  # VERBOSE level
    log.emit(EventType.DEBUG, message="debug")  # VERBOSE level (DEBUG events are at VERBOSE)

    # Filter by type
    state_events = log.get_events(event_type=EventType.STATE_DETECTED)
    assert len(state_events) == 1

    # Filter by verbosity
    quiet_events = log.get_events(max_verbosity=Verbosity.QUIET)
    assert len(quiet_events) == 1  # Only CYCLE_START

    normal_events = log.get_events(max_verbosity=Verbosity.NORMAL)
    assert len(normal_events) == 2  # CYCLE_START + STATE_DETECTED

    verbose_events = log.get_events(max_verbosity=Verbosity.VERBOSE)
    assert len(verbose_events) == 4  # All events (FILES_SELECTED and DEBUG are both VERBOSE)

    print("  PASSED")


def test_get_event_log_helper():
    """Test get_event_log() helper function."""
    print("Test: get_event_log() helper")

    # With no events in state
    state = {"foo": "bar"}
    log, events = get_event_log(state)

    assert events == []
    assert isinstance(log, EventLog)

    # Events list should be usable
    log.emit(EventType.CYCLE_START, cycle=1)
    assert len(events) == 1

    # With existing events
    state2 = {"events": [{"type": "existing"}]}
    log2, events2 = get_event_log(state2)
    assert len(events2) == 1

    print("  PASSED")


def test_emit_event_helper():
    """Test emit_event() helper function."""
    print("Test: emit_event() helper")

    state = {"events": [{"type": "existing"}]}

    new_state = emit_event(state, EventType.STATE_DETECTED, workflow_state="test")

    # Original state should be unchanged
    assert len(state["events"]) == 1

    # New state should have new event
    assert len(new_state["events"]) == 2
    assert new_state["events"][1]["type"] == EventType.STATE_DETECTED

    print("  PASSED")


def test_migrate_debug_log():
    """Test migrate_debug_log() function."""
    print("Test: migrate_debug_log()")

    state = {
        "debug_log": ["Message 1", "Message 2"],
        "events": [{"type": "existing"}],
    }

    new_state = migrate_debug_log(state)

    # Should have 3 events now
    assert len(new_state["events"]) == 3
    assert new_state["events"][1]["type"] == EventType.DEBUG
    assert new_state["events"][1]["message"] == "Message 1"

    print("  PASSED")


def test_formatter_normal_verbosity():
    """Test EventFormatter at normal verbosity."""
    print("Test: EventFormatter normal verbosity")

    events = [
        {"type": EventType.CYCLE_START, "cycle_number": 1},
        {"type": EventType.STATE_DETECTED, "workflow_state": "xray_refining"},
        {"type": EventType.PROGRAM_SELECTED, "program": "phenix.refine", "reasoning": "Test"},
        {"type": EventType.COMMAND_BUILT, "command": "phenix.refine model.pdb"},
    ]

    formatter = EventFormatter(verbosity=Verbosity.NORMAL)
    output = formatter.format_cycle(events, cycle_number=1)

    assert "CYCLE 1" in output
    assert "xray_refining" in output
    assert "phenix.refine" in output
    assert "model.pdb" in output

    print("  PASSED")


def test_formatter_quiet_verbosity():
    """Test EventFormatter at quiet verbosity."""
    print("Test: EventFormatter quiet verbosity")

    events = [
        {"type": EventType.PROGRAM_SELECTED, "program": "phenix.refine"},
    ]

    formatter = EventFormatter(verbosity=Verbosity.QUIET)
    output = formatter.format_cycle(events, cycle_number=3)

    assert output == "CYCLE 3: phenix.refine"

    # Test with error
    events_with_error = [
        {"type": EventType.PROGRAM_SELECTED, "program": "phenix.refine"},
        {"type": EventType.ERROR, "message": "Failed"},
    ]

    output_error = formatter.format_cycle(events_with_error, cycle_number=3)
    assert "ERROR" in output_error
    assert "Failed" in output_error

    print("  PASSED")


def test_formatter_verbose_includes_files():
    """Test that verbose mode includes file selection."""
    print("Test: EventFormatter verbose includes files")

    events = [
        {"type": EventType.FILES_SELECTED, "selections": {
            "model": {"selected": "model.pdb", "score": 85},
        }},
    ]

    # Normal should not show
    formatter_normal = EventFormatter(verbosity=Verbosity.NORMAL)
    output_normal = formatter_normal.format_cycle(events, cycle_number=1)
    assert "model.pdb" not in output_normal

    # Verbose should show
    formatter_verbose = EventFormatter(verbosity=Verbosity.VERBOSE)
    output_verbose = formatter_verbose.format_cycle(events, cycle_number=1)
    assert "model.pdb" in output_verbose

    print("  PASSED")


def test_formatter_full_reasoning():
    """Test that reasoning is not truncated."""
    print("Test: EventFormatter full reasoning (no truncation)")

    long_reasoning = "Word" * 200  # 800 chars of continuous text

    events = [
        {"type": EventType.PROGRAM_SELECTED, "program": "test", "reasoning": long_reasoning},
    ]

    formatter = EventFormatter(verbosity=Verbosity.NORMAL)
    output = formatter.format_cycle(events, cycle_number=1)

    # The full reasoning should appear (may be word-wrapped but not truncated)
    # Count total "Word" occurrences - should be 200, not truncated to ~50
    count = output.count("Word")
    assert count == 200, f"Expected 200 occurrences, got {count} - reasoning was truncated"

    # Also verify no "..." truncation marker
    assert "..." not in output or output.count("...") == 0, "Found truncation marker"

    print("  PASSED")


def test_user_request_invalid_event():
    """Test USER_REQUEST_INVALID event formatting."""
    print("Test: USER_REQUEST_INVALID event formatting")

    events = [
        {
            "type": EventType.USER_REQUEST_INVALID,
            "requested_program": "phenix.xxx",
            "reason": "Unknown program (not supported by PHENIX AI Agent)",
            "selected_instead": "phenix.refine",
            "valid_programs": ["phenix.refine", "phenix.molprobity"],
            "suggestion": "Check the program name",
        },
        {"type": EventType.PROGRAM_SELECTED, "program": "phenix.refine"},
    ]

    # Should be shown at all verbosity levels (QUIET level)
    formatter_quiet = EventFormatter(verbosity=Verbosity.QUIET)
    output_quiet = formatter_quiet.format_cycle(events, cycle_number=1)
    assert "WARNING" in output_quiet
    assert "phenix.xxx" in output_quiet
    assert "not available" in output_quiet.lower() or "requested program" in output_quiet.lower()

    # Should also be prominent at normal level
    formatter_normal = EventFormatter(verbosity=Verbosity.NORMAL)
    output_normal = formatter_normal.format_cycle(events, cycle_number=1)
    assert "WARNING" in output_normal
    assert "phenix.xxx" in output_normal
    assert "phenix.refine" in output_normal

    print("  PASSED")


def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.tst_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
