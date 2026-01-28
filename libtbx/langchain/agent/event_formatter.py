"""
Event formatter for PHENIX AI Agent output.

Converts structured events to human-readable text for display.
Produces consistent output for both local and remote execution.
"""

from __future__ import absolute_import, division, print_function

from libtbx.langchain.agent.event_log import (
    EventType,
    Verbosity,
    EVENT_VERBOSITY,
    verbosity_index,
)


class EventFormatter:
    """
    Formats events for human-readable display.

    Supports multiple verbosity levels and produces consistent
    output for both local and remote execution.

    Usage:
        formatter = EventFormatter(verbosity=Verbosity.NORMAL)
        output = formatter.format_cycle(events, cycle_number=3)
        print(output, file=logger)
    """

    def __init__(self, verbosity=Verbosity.NORMAL, width=80):
        """
        Initialize formatter.

        Args:
            verbosity: Verbosity level (quiet, normal, verbose, debug)
            width: Line width for formatting (default 80)
        """
        self.verbosity = verbosity
        self._width = width

    def format_cycle(self, events, cycle_number=None):
        """
        Format all events for a cycle.

        Args:
            events: List of event dicts
            cycle_number: Cycle number for header (optional)

        Returns:
            str: Formatted output
        """
        if not events:
            return ""

        # For quiet mode, generate one-liner from ALL events (unfiltered)
        if self.verbosity == Verbosity.QUIET:
            return self._format_quiet(events, cycle_number)

        # Filter events by verbosity for other modes
        filtered = self._filter_by_verbosity(events)

        if not filtered:
            return ""

        lines = []

        # Header
        lines.append(self._header(cycle_number))

        # Format each event
        for event in filtered:
            formatted = self._format_event(event)
            if formatted:
                lines.append(formatted)

        # Footer
        lines.append(self._separator("-"))

        return "\n".join(lines)

    def format_single(self, event):
        """
        Format a single event.

        Args:
            event: Event dict

        Returns:
            str: Formatted output, or empty string if filtered out
        """
        # Check verbosity
        event_type = event.get("type", "")
        event_verbosity = EVENT_VERBOSITY.get(event_type, Verbosity.NORMAL)

        if verbosity_index(event_verbosity) > verbosity_index(self.verbosity):
            return ""

        return self._format_event(event) or ""

    def _filter_by_verbosity(self, events):
        """Filter events based on current verbosity level."""
        current_idx = verbosity_index(self.verbosity)

        return [
            e for e in events
            if verbosity_index(EVENT_VERBOSITY.get(e.get("type"), Verbosity.NORMAL)) <= current_idx
        ]

    def _header(self, cycle_number):
        """Format cycle header."""
        if cycle_number:
            title = " CYCLE %d" % cycle_number
        else:
            title = " AGENT DECISION"
        return self._separator("=") + "\n" + title + "\n" + self._separator("=")

    def _separator(self, char):
        """Create a separator line."""
        return char * self._width

    def _format_quiet(self, events, cycle_number):
        """Format for quiet verbosity (single line, but warnings shown)."""
        program = None
        error = None
        user_request_invalid = None

        for e in events:
            etype = e.get("type")
            if etype == EventType.PROGRAM_SELECTED:
                program = e.get("program", "unknown")
            elif etype == EventType.ERROR:
                error = e.get("message", "error")
            elif etype == EventType.USER_REQUEST_INVALID:
                user_request_invalid = e

        lines = []

        # Always show user request invalid warning prominently
        if user_request_invalid:
            lines.append(self._format_user_request_invalid(user_request_invalid))

        if error:
            lines.append("CYCLE %s: ERROR - %s" % (cycle_number or "?", error))
        elif program:
            lines.append("CYCLE %s: %s" % (cycle_number or "?", program))

        return "\n".join(lines) if lines else ""

    def _format_event(self, event):
        """Format a single event based on its type."""
        event_type = event.get("type", "")

        # Map event types to formatter methods
        formatters = {
            EventType.STATE_DETECTED: self._format_state_detected,
            EventType.METRICS_EXTRACTED: self._format_metrics_extracted,
            EventType.METRICS_TREND: self._format_metrics_trend,
            EventType.SANITY_CHECK: self._format_sanity_check,
            EventType.PROGRAM_SELECTED: self._format_program_selected,
            EventType.PROGRAM_MODIFIED: self._format_program_modified,
            EventType.USER_REQUEST_INVALID: self._format_user_request_invalid,
            EventType.FILES_SELECTED: self._format_files_selected,
            EventType.FILE_SCORED: self._format_file_scored,
            EventType.COMMAND_BUILT: self._format_command_built,
            EventType.STOP_DECISION: self._format_stop_decision,
            EventType.DIRECTIVE_APPLIED: self._format_directive_applied,
            EventType.ERROR: self._format_error,
            EventType.WARNING: self._format_warning,
            EventType.DEBUG: self._format_debug,
            EventType.THOUGHT: self._format_thought,
        }

        formatter = formatters.get(event_type)
        if formatter:
            return formatter(event)
        return None

    # -------------------------------------------------------------------------
    # Individual event formatters
    # -------------------------------------------------------------------------

    def _format_state_detected(self, event):
        """Format state detection event."""
        lines = [
            "",
            "State: %s" % event.get("workflow_state", event.get("state", "unknown")),
        ]
        if event.get("experiment_type"):
            lines.append("  Experiment type: %s" % event["experiment_type"])
        if event.get("reason"):
            lines.append("  Reason: %s" % event["reason"])
        return "\n".join(lines)

    def _format_metrics_extracted(self, event):
        """Format metrics extraction event."""
        lines = ["", "Metrics:"]

        # Define metrics with their display names and formats
        metrics = [
            ("r_free", "R-free", "%.4f"),
            ("r_work", "R-work", "%.4f"),
            ("resolution", "Resolution", "%.2f A"),
            ("map_cc", "Map CC", "%.4f"),
            ("clashscore", "Clashscore", "%.1f"),
            ("tfz", "TFZ", "%.1f"),
            ("llg", "LLG", "%.0f"),
            ("ramachandran_outliers", "Ramachandran outliers", "%.1f%%"),
            ("rotamer_outliers", "Rotamer outliers", "%.1f%%"),
        ]

        has_metrics = False
        for key, label, fmt in metrics:
            if event.get(key) is not None:
                has_metrics = True
                value = event[key]
                prev_key = key + "_prev"
                prev = event.get(prev_key)

                if prev is not None:
                    diff = value - prev
                    # For R-factors, clashscore, outliers: lower is better
                    if key in ("r_free", "r_work", "clashscore", "ramachandran_outliers", "rotamer_outliers"):
                        direction = "improved" if diff < 0 else "worsened"
                    else:
                        direction = "improved" if diff > 0 else "worsened"

                    # Format with change
                    val_str = fmt % value
                    prev_str = (fmt % prev).split()[0]  # Remove units for prev
                    lines.append("  %s: %s -> %s (%s by %.4f)" % (label, prev_str, val_str, direction, abs(diff)))
                else:
                    lines.append("  %s: %s" % (label, fmt % value))

        return "\n".join(lines) if has_metrics else None

    def _format_metrics_trend(self, event):
        """Format metrics trend event."""
        if event.get("plateau"):
            cycles = event.get("cycles_analyzed", "?")
            return "  Trend: PLATEAU detected after %s cycles" % cycles
        elif event.get("improving"):
            return "  Trend: Improving (no plateau detected)"
        else:
            return "  Trend: Not improving"

    def _format_sanity_check(self, event):
        """Format sanity check event."""
        if event.get("passed"):
            return None  # Don't show if passed

        lines = ["", "Sanity Check Issues:"]
        for flag in event.get("red_flags", []):
            lines.append("  [RED FLAG] %s" % flag)
        for warning in event.get("warnings", []):
            lines.append("  [WARNING] %s" % warning)

        return "\n".join(lines) if len(lines) > 1 else None

    def _format_program_selected(self, event):
        """Format program selection event."""
        lines = [
            "",
            "Decision: %s" % event.get("program", "unknown"),
        ]

        # Full reasoning (no truncation!)
        reasoning = event.get("reasoning", "")
        if reasoning:
            wrapped = self._wrap_text(str(reasoning), width=76, indent="  ")
            lines.append("  Reasoning: %s" % wrapped)

        source = event.get("source", "")
        if source:
            provider = event.get("provider", "")
            if provider and source == "llm":
                lines.append("  Source: %s (%s)" % (source, provider))
            else:
                lines.append("  Source: %s" % source)

        return "\n".join(lines)

    def _format_program_modified(self, event):
        """Format program modification event."""
        original = event.get("original", "?")
        selected = event.get("selected", "?")
        reason = event.get("reason", "")
        return "  Note: Changed from '%s' to '%s' (%s)" % (original, selected, reason)

    def _format_user_request_invalid(self, event):
        """Format user request invalid event (always shown prominently)."""
        requested = event.get("requested_program", "unknown")
        reason = event.get("reason", "not available")
        selected = event.get("selected_instead", "")
        valid_programs = event.get("valid_programs", [])

        lines = [
            "",
            "=" * 60,
            "  WARNING: Requested program not available",
            "=" * 60,
            "  You requested: %s" % requested,
            "  Reason: %s" % reason,
        ]

        if selected:
            lines.append("  Running instead: %s" % selected)

        if valid_programs:
            lines.append("  Available programs: %s" % ", ".join(valid_programs))

        suggestion = event.get("suggestion", "")
        if suggestion:
            lines.append("  Suggestion: %s" % suggestion)

        lines.append("=" * 60)
        return "\n".join(lines)

    def _format_files_selected(self, event):
        """Format file selection event (verbose only)."""
        lines = ["", "File Selection:"]

        selections = event.get("selections", {})
        for input_type, info in selections.items():
            if isinstance(info, dict):
                selected = info.get("selected", "none")
                lines.append("  %s: %s" % (input_type.title(), selected))

                if info.get("reason"):
                    lines.append("    Reason: %s" % info["reason"])

                if info.get("score") is not None:
                    lines.append("    Score: %s" % info["score"])

                candidates = info.get("candidates", [])
                if candidates and len(candidates) > 1:
                    others = []
                    for c in candidates[1:4]:  # Show top 3 alternatives
                        if isinstance(c, dict):
                            others.append("%s (score=%s)" % (c.get("file", "?"), c.get("score", "?")))
                        else:
                            others.append(str(c))
                    if others:
                        lines.append("    Other candidates: %s" % ", ".join(others))
            else:
                # Simple string value
                lines.append("  %s: %s" % (input_type.title(), info))

        return "\n".join(lines) if len(lines) > 1 else None

    def _format_file_scored(self, event):
        """Format file scoring event (verbose only)."""
        filename = event.get("file", "?")
        score = event.get("score", "?")
        components = event.get("components", {})

        lines = ["  Scored: %s = %s" % (filename, score)]
        if components:
            comp_strs = ["%s=%s" % (k, v) for k, v in components.items()]
            lines.append("    Components: %s" % ", ".join(comp_strs))

        return "\n".join(lines)

    def _format_command_built(self, event):
        """Format command built event."""
        command = event.get("command", "")
        if not command:
            return None
        return "\nCommand:\n  %s" % command

    def _format_stop_decision(self, event):
        """Format stop decision event."""
        if event.get("stop"):
            reason = event.get("reason", "unknown reason")
            return "\nStopping: %s" % reason
        return None

    def _format_directive_applied(self, event):
        """Format directive application event."""
        directive = event.get("directive", "")
        action = event.get("action", "")
        return "  Directive applied: %s -> %s" % (directive, action)

    def _format_error(self, event):
        """Format error event."""
        message = event.get("message", "Unknown error")
        lines = ["", "ERROR: %s" % message]
        if event.get("details"):
            lines.append("  Details: %s" % event["details"])
        return "\n".join(lines)

    def _format_warning(self, event):
        """Format warning event."""
        message = event.get("message", "")
        return "WARNING: %s" % message

    def _format_debug(self, event):
        """Format debug event (only shown at debug verbosity)."""
        message = event.get("message", "")
        return "[DEBUG] %s" % message

    def _format_thought(self, event):
        """Format thought/reasoning trace event."""
        content = event.get("content", "")
        if not content:
            return None
        wrapped = self._wrap_text(str(content), width=78, indent="  ")
        return "\n[THOUGHT]\n  %s" % wrapped

    # -------------------------------------------------------------------------
    # Utility methods
    # -------------------------------------------------------------------------

    def _wrap_text(self, text, width=76, indent=""):
        """Word-wrap text at specified width."""
        if not text:
            return ""

        words = text.split()
        if not words:
            return ""

        lines = []
        current_line = []
        current_length = 0

        for word in words:
            word_len = len(word)
            if current_length + word_len + 1 <= width:
                current_line.append(word)
                current_length += word_len + 1
            else:
                if current_line:
                    lines.append(" ".join(current_line))
                current_line = [word]
                current_length = word_len

        if current_line:
            lines.append(" ".join(current_line))

        if len(lines) > 1:
            return lines[0] + "\n" + "\n".join(indent + line for line in lines[1:])
        return lines[0]


# =============================================================================
# VALIDATION
# =============================================================================

if __name__ == "__main__":
    print("EventFormatter validation")
    print("=" * 60)

    # Create sample events
    sample_events = [
        {"type": EventType.CYCLE_START, "cycle_number": 3},
        {"type": EventType.STATE_DETECTED, "workflow_state": "xray_refining",
         "experiment_type": "X-ray crystallography", "reason": "Has model and MTZ data"},
        {"type": EventType.METRICS_EXTRACTED, "r_free": 0.24, "r_free_prev": 0.26,
         "r_work": 0.21, "r_work_prev": 0.22, "resolution": 2.5},
        {"type": EventType.METRICS_TREND, "improving": True, "plateau": False, "cycles_analyzed": 3},
        {"type": EventType.PROGRAM_SELECTED, "program": "phenix.refine",
         "reasoning": "R-free continues to improve (0.26 -> 0.24). No plateau detected after 3 cycles. Continue refinement to optimize geometry and reduce R-factors. The model shows good density fit with reasonable geometry.",
         "source": "llm", "provider": "google"},
        {"type": EventType.FILES_SELECTED, "selections": {
            "model": {"selected": "refine_002_001.pdb", "score": 85, "reason": "Highest cycle, best R-free"},
            "data": {"selected": "data.mtz", "reason": "Original data file"},
        }},
        {"type": EventType.COMMAND_BUILT, "command": "phenix.refine refine_002_001.pdb data.mtz output.prefix=refine_003"},
        {"type": EventType.DEBUG, "message": "This is a debug message"},
    ]

    # Test different verbosity levels
    for level in [Verbosity.QUIET, Verbosity.NORMAL, Verbosity.VERBOSE]:
        print("\n" + "-" * 60)
        print("VERBOSITY: %s" % level)
        print("-" * 60)

        formatter = EventFormatter(verbosity=level)
        output = formatter.format_cycle(sample_events, cycle_number=3)
        print(output)

    print("\n" + "=" * 60)
    print("Validation complete!")
