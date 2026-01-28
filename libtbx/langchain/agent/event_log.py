"""
Event logging for PHENIX AI Agent transparency.

This module provides structured event logging that captures
the agent's decision-making process for display to users.

Design Notes:
- Events are stored as a list of dicts in state["events"] for JSON serialization
- EventLog is a helper class that wraps the list with convenience methods
- Each node retrieves events from state, adds new events, returns updated list
- This maintains LangGraph's state immutability pattern

Usage in LangGraph nodes:

    from libtbx.langchain.agent.event_log import EventType, get_event_log

    def perceive(state):
        # Get or create events list and wrapper
        log, events = get_event_log(state)

        # Add events
        log.emit(EventType.STATE_DETECTED, state="xray_refining", ...)

        # Return updated state (events list is modified in place)
        return {**state, "events": events}
"""

from __future__ import absolute_import, division, print_function


# =============================================================================
# EVENT TYPES - Strict constants for filtering and categorization
# =============================================================================

class EventType:
    """
    Event type constants.

    Using a class with string constants allows:
    - IDE autocomplete
    - Easy filtering by category
    - Clear documentation of valid types
    """
    # High-level cycle events
    CYCLE_START = "cycle_start"
    CYCLE_COMPLETE = "cycle_complete"

    # State and analysis
    STATE_DETECTED = "state_detected"
    METRICS_EXTRACTED = "metrics_extracted"
    METRICS_TREND = "metrics_trend"
    SANITY_CHECK = "sanity_check"

    # Decision events
    PROGRAM_SELECTED = "program_selected"
    PROGRAM_MODIFIED = "program_modified"
    STOP_DECISION = "stop_decision"
    DIRECTIVE_APPLIED = "directive_applied"
    USER_REQUEST_INVALID = "user_request_invalid"  # User requested unavailable program

    # File selection events
    FILES_SELECTED = "files_selected"
    FILE_SCORED = "file_scored"

    # Command events
    COMMAND_BUILT = "command_built"

    # Error/diagnostic events
    ERROR = "error"
    WARNING = "warning"

    # Debug/trace events (maps to old debug_log)
    DEBUG = "debug"
    THOUGHT = "thought"  # LLM chain-of-thought/reasoning traces


# =============================================================================
# VERBOSITY LEVELS
# =============================================================================

class Verbosity:
    """
    Verbosity level constants.

    Three levels:
      QUIET:   Errors, warnings, and final result only
      NORMAL:  Key decisions, metrics, reasoning (default)
      VERBOSE: Full detail including file selection and internal state
    """
    QUIET = "quiet"      # Errors and final result only
    NORMAL = "normal"    # Key decisions and metrics
    VERBOSE = "verbose"  # Full detail including file selection, internal state


# Event type to minimum verbosity level mapping
# Events at a given level are shown at that level and above
EVENT_VERBOSITY = {
    EventType.CYCLE_START: Verbosity.QUIET,
    EventType.CYCLE_COMPLETE: Verbosity.QUIET,
    EventType.STATE_DETECTED: Verbosity.NORMAL,
    EventType.METRICS_EXTRACTED: Verbosity.NORMAL,
    EventType.METRICS_TREND: Verbosity.NORMAL,
    EventType.SANITY_CHECK: Verbosity.NORMAL,
    EventType.PROGRAM_SELECTED: Verbosity.NORMAL,
    EventType.PROGRAM_MODIFIED: Verbosity.NORMAL,
    EventType.STOP_DECISION: Verbosity.NORMAL,
    EventType.DIRECTIVE_APPLIED: Verbosity.NORMAL,
    EventType.USER_REQUEST_INVALID: Verbosity.QUIET,  # Always show - user needs to know
    EventType.FILES_SELECTED: Verbosity.VERBOSE,
    EventType.FILE_SCORED: Verbosity.VERBOSE,
    EventType.COMMAND_BUILT: Verbosity.NORMAL,
    EventType.ERROR: Verbosity.QUIET,
    EventType.WARNING: Verbosity.NORMAL,
    EventType.DEBUG: Verbosity.VERBOSE,    # DEBUG events shown at VERBOSE level
    EventType.THOUGHT: Verbosity.VERBOSE,  # THOUGHT events shown at VERBOSE level
}

# Verbosity level ordering (for filtering)
VERBOSITY_ORDER = [Verbosity.QUIET, Verbosity.NORMAL, Verbosity.VERBOSE]


def verbosity_index(level):
    """Get numeric index for a verbosity level."""
    try:
        return VERBOSITY_ORDER.index(level)
    except ValueError:
        return 1  # Default to NORMAL


# =============================================================================
# EVENT LOG CLASS
# =============================================================================

class EventLog:
    """
    Helper class for working with the events list in state.

    This class wraps a list of event dicts and provides convenience methods.
    It does NOT own the list - it operates on the list stored in state["events"].

    Usage in LangGraph nodes:

        def perceive(state):
            # Get or create events list
            events = state.get("events", [])
            log = EventLog(events)

            # Add events
            log.emit(EventType.STATE_DETECTED, state="xray_refining", ...)

            # Return updated state (events list is modified in place)
            return {**state, "events": events}

    For immutable patterns, use a copy:

        def perceive(state):
            events = list(state.get("events", []))  # Copy
            log = EventLog(events)
            log.emit(EventType.STATE_DETECTED, ...)
            return {**state, "events": events}
    """

    def __init__(self, events_list=None):
        """
        Initialize with an events list.

        Args:
            events_list: List to store events in. If None, creates new list.
        """
        self._events = events_list if events_list is not None else []

    @property
    def events(self):
        """Get the underlying events list."""
        return self._events

    def emit(self, event_type, **kwargs):
        """
        Record an event.

        Args:
            event_type: One of EventType constants
            **kwargs: Event-specific data

        Returns:
            The created event dict
        """
        event = {
            "type": event_type,
            **kwargs
        }
        self._events.append(event)
        return event

    def emit_debug(self, message):
        """
        Emit a debug event (replacement for old _log/debug_log).

        Args:
            message: Debug message string

        Returns:
            The created event dict
        """
        return self.emit(EventType.DEBUG, message=message)

    def emit_thought(self, thought):
        """
        Emit a thought/reasoning trace event.

        Args:
            thought: LLM reasoning or chain-of-thought text

        Returns:
            The created event dict
        """
        return self.emit(EventType.THOUGHT, content=thought)

    def emit_error(self, message, details=None):
        """
        Emit an error event.

        Args:
            message: Error message
            details: Optional additional details

        Returns:
            The created event dict
        """
        return self.emit(EventType.ERROR, message=message, details=details)

    def emit_warning(self, message):
        """
        Emit a warning event.

        Args:
            message: Warning message

        Returns:
            The created event dict
        """
        return self.emit(EventType.WARNING, message=message)

    def get_events(self, event_type=None, max_verbosity=None):
        """
        Get events, optionally filtered.

        Args:
            event_type: Filter to specific type (or list of types)
            max_verbosity: Filter to events at or below this verbosity level

        Returns:
            List of matching events
        """
        result = self._events

        # Filter by type
        if event_type:
            if isinstance(event_type, (list, tuple)):
                result = [e for e in result if e.get("type") in event_type]
            else:
                result = [e for e in result if e.get("type") == event_type]

        # Filter by verbosity
        if max_verbosity:
            max_idx = verbosity_index(max_verbosity)
            result = [
                e for e in result
                if verbosity_index(EVENT_VERBOSITY.get(e.get("type"), Verbosity.NORMAL)) <= max_idx
            ]

        return result

    def has_errors(self):
        """Check if any error events were recorded."""
        return any(e.get("type") == EventType.ERROR for e in self._events)

    def has_warnings(self):
        """Check if any warning events were recorded."""
        return any(e.get("type") == EventType.WARNING for e in self._events)

    def get_last(self, event_type):
        """Get the most recent event of a given type."""
        for event in reversed(self._events):
            if event.get("type") == event_type:
                return event
        return None

    def get_first(self, event_type):
        """Get the first event of a given type."""
        for event in self._events:
            if event.get("type") == event_type:
                return event
        return None

    def clear(self):
        """Clear all events."""
        self._events.clear()

    def __len__(self):
        return len(self._events)

    def __iter__(self):
        return iter(self._events)

    def __bool__(self):
        return len(self._events) > 0


# =============================================================================
# HELPER FUNCTIONS FOR USE IN GRAPH NODES
# =============================================================================

def get_event_log(state):
    """
    Get an EventLog wrapper for the state's events list.

    Creates the events list if it doesn't exist.

    Args:
        state: LangGraph state dict

    Returns:
        tuple: (EventLog instance, events list)

    Example:
        def my_node(state):
            log, events = get_event_log(state)
            log.emit(EventType.STATE_DETECTED, ...)
            return {**state, "events": events}
    """
    events = state.get("events")
    if events is None:
        events = []
    log = EventLog(events)
    return log, events


def emit_event(state, event_type, **kwargs):
    """
    Convenience function to emit a single event and return updated state.

    This creates a copy of the events list for immutability.

    Args:
        state: LangGraph state dict
        event_type: Event type constant
        **kwargs: Event data

    Returns:
        Updated state dict with new event

    Example:
        state = emit_event(state, EventType.ERROR, message="Something failed")
    """
    events = list(state.get("events", []))  # Copy for immutability
    events.append({"type": event_type, **kwargs})
    return {**state, "events": events}


def migrate_debug_log(state):
    """
    Migrate old debug_log entries to new event system.

    Call this at the start of processing to preserve backward compatibility.
    Converts each debug_log string to a DEBUG event.

    Args:
        state: LangGraph state dict

    Returns:
        Updated state with debug_log entries converted to DEBUG events
    """
    debug_log = state.get("debug_log", [])
    if not debug_log:
        return state

    events = list(state.get("events", []))

    for msg in debug_log:
        if isinstance(msg, str):
            events.append({"type": EventType.DEBUG, "message": msg})

    return {**state, "events": events}


# =============================================================================
# VALIDATION
# =============================================================================

if __name__ == "__main__":
    # Basic validation
    print("EventLog module validation")
    print("=" * 40)

    # Test EventType constants
    print("\nEvent Types:")
    for attr in dir(EventType):
        if not attr.startswith("_"):
            print(f"  {attr}: {getattr(EventType, attr)}")

    # Test Verbosity constants
    print("\nVerbosity Levels:")
    for level in VERBOSITY_ORDER:
        print(f"  {level}")

    # Test EventLog
    print("\nEventLog test:")
    events = []
    log = EventLog(events)

    log.emit(EventType.CYCLE_START, cycle_number=1)
    log.emit(EventType.STATE_DETECTED, state="xray_refining", experiment_type="X-ray")
    log.emit(EventType.METRICS_EXTRACTED, r_free=0.24, r_work=0.21)
    log.emit(EventType.PROGRAM_SELECTED, program="phenix.refine", reasoning="Continue refinement")
    log.emit_debug("This is a debug message")
    log.emit_error("Test error", details="Some details")

    print(f"  Total events: {len(log)}")
    print(f"  Has errors: {log.has_errors()}")
    print(f"  Events at NORMAL verbosity: {len(log.get_events(max_verbosity=Verbosity.NORMAL))}")
    print(f"  Events at VERBOSE verbosity: {len(log.get_events(max_verbosity=Verbosity.VERBOSE))}")

    # Test helper functions
    print("\nHelper function test:")
    state = {"cycle_number": 1}
    log, events = get_event_log(state)
    log.emit(EventType.CYCLE_START, cycle_number=1)
    state = {**state, "events": events}
    print(f"  State has {len(state['events'])} events")

    # Test emit_event
    state = emit_event(state, EventType.STATE_DETECTED, workflow_state="test")
    print(f"  After emit_event: {len(state['events'])} events")

    print("\n" + "=" * 40)
    print("Validation complete!")
