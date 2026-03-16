"""
Event formatter for PHENIX AI Agent output.

Converts structured events to human-readable text for display.
Produces consistent output for both local and remote execution.
"""

from __future__ import absolute_import, division, print_function

try:
    from libtbx.langchain.agent.event_log import (
        EventType,
        Verbosity,
        EVENT_VERBOSITY,
        verbosity_index,
    )
except ImportError:
    from agent.event_log import (
        EventType,
        Verbosity,
        EVENT_VERBOSITY,
        verbosity_index,
    )


class EventFormatter:
    """
    Formats events for human-readable display.

    Supports multiple verbosity levels and produces consistent
    output for both local and remote execution. Decision blocks
    are wrapped in box frames for visual clarity and log
    searchability.

    Usage:
        formatter = EventFormatter(verbosity=Verbosity.NORMAL)
        output = formatter.format_cycle(events, cycle_number=3)
        print(output, file=logger)

        # After program execution:
        result_output = formatter.format_result_block(
            status="SUCCESS", metrics={"r_free": 0.230}, ...)
        print(result_output, file=logger)
    """

    # Tag markers for log searchability (grep-friendly)
    TAG_DECISION = "[DECISION]"
    TAG_REASONING = "[REASONING]"
    TAG_COMMAND = "[COMMAND]"
    TAG_RESULT = "[RESULT]"
    TAG_GATE = "[GATE]"
    TAG_EXPERT = "[EXPERT]"
    TAG_STOP = "[STOP]"
    TAG_STATE = "[STATE]"

    # Box-drawing characters
    _BOX_TL = "\u250c"  # ┌
    _BOX_TR = "\u2510"  # ┐
    _BOX_BL = "\u2514"  # └
    _BOX_BR = "\u2518"  # ┘
    _BOX_H = "\u2500"   # ─
    _BOX_V = "\u2502"   # │
    _BOX_ML = "\u251c"  # ├
    _BOX_MR = "\u2524"  # ┤

    # Default inner width for box content (72 fits in 80-char
    # terminal with box chars + padding)
    _BOX_INNER = 70

    def __init__(self, verbosity=Verbosity.NORMAL, width=80):
        """
        Initialize formatter.

        Args:
            verbosity: Verbosity level (quiet, normal, verbose, debug)
            width: Line width for formatting (default 80)
        """
        self.verbosity = verbosity
        self._width = width
        self._BOX_INNER = width - 4  # 2 for "│ " + 2 for " │"

    def format_cycle(self, events, cycle_number=None):
        """
        Format all events for a cycle.

        At NORMAL and VERBOSE verbosity, produces a boxed
        decision block with tagged lines.  At QUIET, produces
        a one-liner summary.

        Args:
            events: List of event dicts
            cycle_number: Cycle number for header (optional)

        Returns:
            str: Formatted output
        """
        if not events:
            return ""

        # Use boxed format for normal/verbose
        return self.format_decision_block(
            events, cycle_number=cycle_number)

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
        expert_stop = None
        stop_reason = None

        for e in events:
            etype = e.get("type")
            if etype == EventType.PROGRAM_SELECTED:
                program = e.get("program", "unknown")
            elif etype == EventType.ERROR:
                error = e.get("message", "error")
            elif etype == EventType.USER_REQUEST_INVALID:
                user_request_invalid = e
            elif etype == EventType.EXPERT_ASSESSMENT:
                if e.get("action") == "stop":
                    expert_stop = e
            elif etype == EventType.STOP_DECISION:
                if e.get("stop"):
                    stop_reason = e.get("reason", "")

        lines = []

        # Always show user request invalid warning prominently
        if user_request_invalid:
            lines.append(self._format_user_request_invalid(user_request_invalid))

        if error:
            lines.append("CYCLE %s: ERROR - %s" % (cycle_number or "?", error))
        elif expert_stop:
            # Expert recommended stopping — always show.
            # Use the analysis text directly (stop_reason
            # has an "expert:" prefix which would double up).
            reason = (
                expert_stop.get("analysis", "")[:80]
                or stop_reason
                or "expert recommendation"
            )
            lines.append(
                "CYCLE %s: STOP (expert: %s)"
                % (cycle_number or "?", reason)
            )
        elif stop_reason:
            lines.append(
                "CYCLE %s: STOP (%s)"
                % (cycle_number or "?", stop_reason)
            )
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
            EventType.EXPERT_ASSESSMENT: self._format_expert_assessment,
            EventType.FILES_SELECTED: self._format_files_selected,
            EventType.FILE_SCORED: self._format_file_scored,
            EventType.COMMAND_BUILT: self._format_command_built,
            EventType.STOP_DECISION: self._format_stop_decision,
            EventType.DIRECTIVE_APPLIED: self._format_directive_applied,
            EventType.ERROR: self._format_error,
            EventType.WARNING: self._format_warning,
            EventType.NOTICE: self._format_notice,
            EventType.DEBUG: self._format_debug,
            EventType.THOUGHT: self._format_thought,
        }

        formatter = formatters.get(event_type)
        if formatter:
            return formatter(event)

        # String-based dispatch for new event types
        # not yet in the EventType enum.
        if event_type == "phase_transition":
            return self._format_phase_transition(event)

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
                    try:
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
                    except (TypeError, ValueError):
                        # Guard: metric value may be a string (e.g. twin_law, space_group)
                        lines.append("  %s: %s" % (label, str(value)))
                else:
                    try:
                        lines.append("  %s: %s" % (label, fmt % value))
                    except (TypeError, ValueError):
                        # Guard: metric value may be a string (e.g. twin_law, space_group)
                        lines.append("  %s: %s" % (label, str(value)))

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

    def _format_expert_assessment(self, event):
        """Format expert assessment event (thinking agent output).

        Shown at NORMAL verbosity. Presents the expert's
        analysis alongside any structural validation data
        so the user understands what the agent saw and why
        it made its decision.
        """
        action = event.get("action", "let_run")
        confidence = event.get("confidence", "medium")
        thinking_level = event.get("thinking_level", "")
        analysis = event.get("analysis", "")
        guidance = event.get("guidance", "")
        concerns = event.get("concerns", [])
        validation = event.get("validation_summary", "")
        validation_skip = event.get(
            "validation_skip_reason", "")
        validated_model = event.get("validated_model", "")
        rfree_trend = event.get("rfree_trend", "")
        kb_matched = event.get("kb_rules_matched", False)

        # --- Header ---
        header = "Expert Assessment"
        if confidence:
            header += " (%s confidence)" % confidence
        if action == "stop":
            header += " -- STOP RECOMMENDED"
        lines = ["", header]

        # --- Structure Model summary (v114) ---
        sm_summary = event.get(
            "structure_model_summary", "")
        if sm_summary:
            lines.append("  Structure:")
            for sline in sm_summary.split("\n"):
                sline = sline.strip()
                if sline:
                    lines.append("    %s" % sline)
            # Current problems (nested under Structure)
            problems = event.get(
                "current_problems", [])
            if problems and isinstance(problems, list):
                p_strs = [
                    str(p.get("problem", ""))
                    for p in problems[:3] if p
                ]
                p_strs = [s for s in p_strs if s]
                if p_strs:
                    lines.append(
                        "    Problems: %s"
                        % "; ".join(p_strs))

        # --- Structural validation (advanced mode) ---
        if validation:
            label = "Structural validation"
            if validated_model:
                label += " (%s)" % validated_model
            lines.append("  %s:" % label)
            for vline in validation.split("\n"):
                vline = vline.strip()
                if vline:
                    lines.append("    %s" % vline)
        elif validation_skip:
            lines.append(
                "  Validation: not available (%s)"
                % validation_skip)

        # --- R-free trend ---
        if rfree_trend:
            lines.append("  R-free trend: %s" % rfree_trend)

        # --- Analysis ---
        if analysis:
            wrapped = self._wrap_text(
                str(analysis), width=72, indent="    ")
            lines.append("  Analysis: %s" % wrapped)

        # --- Guidance ---
        if guidance:
            label = "Recommendation" if action == "stop" \
                else "Guidance"
            wrapped = self._wrap_text(
                str(guidance), width=72, indent="    ")
            lines.append("  %s: %s" % (label, wrapped))

        # --- Concerns ---
        if concerns and isinstance(concerns, list):
            shown = [str(c) for c in concerns[:3] if c]
            if shown:
                lines.append(
                    "  Concerns: %s" % "; ".join(shown))

        # --- Mode indicator (subtle) ---
        extras = []
        if thinking_level:
            extras.append("mode=%s" % thinking_level)
        if kb_matched:
            extras.append("KB rules consulted")
        if extras:
            lines.append("  [%s]" % ", ".join(extras))

        # If we only have the header (LLM parse failure
        # or missing fields), add a minimal note so the
        # user knows the expert was consulted.
        if len(lines) <= 2 and not analysis:
            lines.append("  (no detailed assessment available)")

        return "\n".join(lines)

    def _format_phase_transition(self, event):
        """Format phase transition event (v114).

        Shown for advance, retreat, and skip transitions.
        Uses distinctive visual separator for visibility.
        """
        transition = event.get("transition", "advance")
        from_phase = event.get("from_phase", "?")
        to_phase = event.get("to_phase", "?")
        reason = event.get("reason", "")
        blacklisted = event.get("blacklisted", "")

        bar = "=" * 50
        lines = ["", bar]

        if transition == "retreat":
            lines.append(
                " RETREAT: %s -> %s"
                % (from_phase, to_phase))
        elif transition == "skip":
            lines.append(
                " SKIP: %s (-> %s)"
                % (from_phase, to_phase))
        else:
            lines.append(
                " PHASE TRANSITION: %s -> %s"
                % (from_phase, to_phase))

        if reason:
            wrapped = self._wrap_text(
                str(reason), width=46,
                indent="   ")
            lines.append(" Reason: %s" % wrapped)

        if blacklisted:
            lines.append(
                " Strategy blacklisted: %s"
                % blacklisted)

        lines.append(bar)
        return "\n".join(lines)

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

    def _format_notice(self, event):
        """Format notice event — important info the user must see."""
        message = event.get("message", "")
        details = event.get("details", "")
        lines = [
            "",
            "=" * 60,
            "  NOTICE",
            "=" * 60,
            "  %s" % message,
        ]
        if details:
            lines.append("  %s" % details)
        lines.append("=" * 60)
        return "\n".join(lines)

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
    # Box formatting
    # -------------------------------------------------------------------------

    def _box_top(self, title=None):
        """Top border of a box, optionally with a title."""
        if title:
            # ┌─ TITLE ──────────────────────────────┐
            inner = " %s " % title
            fill = self._BOX_INNER - len(inner) + 1
            return (self._BOX_TL + self._BOX_H
                    + inner
                    + self._BOX_H * max(fill, 0)
                    + self._BOX_TR)
        return (self._BOX_TL
                + self._BOX_H * (self._BOX_INNER + 2)
                + self._BOX_TR)

    def _box_mid(self):
        """Middle separator line."""
        return (self._BOX_ML
                + self._BOX_H * (self._BOX_INNER + 2)
                + self._BOX_MR)

    def _box_bottom(self):
        """Bottom border of a box."""
        return (self._BOX_BL
                + self._BOX_H * (self._BOX_INNER + 2)
                + self._BOX_BR)

    def _box_line(self, text=""):
        """A single content line inside a box."""
        # Truncate if too long (shouldn't happen after
        # wrapping, but safety net)
        if len(text) > self._BOX_INNER:
            text = text[:self._BOX_INNER - 1] + "\u2026"
        padding = self._BOX_INNER - len(text)
        return "%s %s%s %s" % (
            self._BOX_V, text, " " * padding, self._BOX_V)

    def _box_wrap_lines(self, text, indent=0):
        """Wrap text to fit inside a box, return list of
        box_line() results."""
        if not text:
            return []
        content_width = self._BOX_INNER - indent
        prefix = " " * indent
        words = str(text).split()
        lines = []
        current = []
        length = 0
        for word in words:
            if length + len(word) + 1 > content_width and current:
                lines.append(
                    self._box_line(prefix + " ".join(current)))
                current = [word]
                length = len(word)
            else:
                current.append(word)
                length += len(word) + 1
        if current:
            lines.append(
                self._box_line(prefix + " ".join(current)))
        return lines

    def _box_command(self, command):
        """Format a command inside a box with backslash
        continuation for long lines, suitable for copy-paste."""
        if not command:
            return []
        lines = []
        max_w = self._BOX_INNER - 2  # indent for "  "
        cmd = str(command).strip()
        if len(cmd) <= max_w:
            lines.append(self._box_line("  " + cmd))
        else:
            # Split at spaces, add backslash continuation
            parts = cmd.split()
            current = parts[0]
            for part in parts[1:]:
                if len(current) + 1 + len(part) > max_w - 2:
                    lines.append(
                        self._box_line("  " + current + " \\"))
                    current = "    " + part
                else:
                    current += " " + part
            lines.append(self._box_line("  " + current))
        return lines

    def format_decision_block(self, events, cycle_number=None):
        """Format all decision events for a cycle as a boxed block.

        Produces a visually distinct decision block with tagged
        lines for log searchability.

        Args:
            events: List of event dicts from the graph pipeline
            cycle_number: Cycle number for the header

        Returns:
            str: Boxed decision block, or empty string in quiet mode
        """
        if self.verbosity == Verbosity.QUIET:
            return self._format_quiet(events, cycle_number)

        filtered = self._filter_by_verbosity(events)
        if not filtered:
            return ""

        lines = []
        title = "CYCLE %d" % cycle_number if cycle_number else "AGENT DECISION"
        lines.append(self._box_top(title))

        for event in filtered:
            etype = event.get("type", "")

            if etype == EventType.STATE_DETECTED:
                state = event.get("workflow_state", "?")
                reason = event.get("reason", "")
                lines.append(self._box_line(
                    "%s  %s" % (self.TAG_STATE, state)))
                if reason:
                    for rl in self._box_wrap_lines(
                            reason, indent=2):
                        lines.append(rl)
                lines.append(self._box_line(""))

            elif etype == EventType.METRICS_EXTRACTED:
                mline = self._format_metrics_line(event)
                if mline:
                    lines.append(self._box_line(mline))

            elif etype == EventType.EXPERT_ASSESSMENT:
                analysis = event.get("analysis", "")
                confidence = event.get(
                    "confidence", "")
                if analysis:
                    label = "%s" % self.TAG_EXPERT
                    if confidence:
                        label += "  (%s confidence)" % confidence
                    lines.append(self._box_line(label))
                    for al in self._box_wrap_lines(
                            analysis, indent=2):
                        lines.append(al)
                    lines.append(self._box_line(""))

            elif etype == EventType.PROGRAM_SELECTED:
                prog = event.get("program", "unknown")
                reasoning = event.get("reasoning", "")
                source = event.get("source", "")
                provider = event.get("provider", "")

                lines.append(self._box_line(
                    "%s  %s" % (self.TAG_DECISION, prog)))
                if reasoning:
                    lines.append(self._box_line(
                        "%s" % self.TAG_REASONING))
                    for rl in self._box_wrap_lines(
                            reasoning, indent=2):
                        lines.append(rl)
                if source:
                    src = "%s (%s)" % (source, provider) \
                        if provider and source == "llm" \
                        else source
                    lines.append(self._box_line(
                        "  Source: %s" % src))
                lines.append(self._box_line(""))

            elif etype == EventType.STOP_DECISION:
                if event.get("stop"):
                    reason = event.get("reason", "")
                    lines.append(self._box_line(
                        "%s" % self.TAG_STOP))
                    for rl in self._box_wrap_lines(
                            reason, indent=2):
                        lines.append(rl)
                    lines.append(self._box_line(""))

            elif etype == EventType.FILES_SELECTED:
                selections = event.get("selections", {})
                if selections:
                    lines.append(self._box_line("Files:"))
                    for slot, info in selections.items():
                        if isinstance(info, dict):
                            fname = info.get("selected", "?")
                            reason = info.get("reason", "")
                            entry = "  %s: %s" % (
                                slot.title(), fname)
                            if reason:
                                entry += " (%s)" % reason
                            lines.append(self._box_line(entry))
                        else:
                            lines.append(self._box_line(
                                "  %s: %s" % (slot.title(), info)))
                    lines.append(self._box_line(""))

            elif etype == EventType.COMMAND_BUILT:
                command = event.get("command", "")
                if command:
                    lines.append(self._box_line(
                        "%s" % self.TAG_COMMAND))
                    for cl in self._box_command(command):
                        lines.append(cl)
                    lines.append(self._box_line(""))

            elif etype == EventType.DIRECTIVE_APPLIED:
                directive = event.get("directive", "")
                action = event.get("action", "")
                lines.append(self._box_line(
                    "  Directive: %s -> %s"
                    % (directive, action)))

            elif etype == EventType.PROGRAM_MODIFIED:
                original = event.get("original", "?")
                selected = event.get("selected", "?")
                reason = event.get("reason", "")
                lines.append(self._box_line(
                    "  Changed: %s -> %s (%s)"
                    % (original, selected, reason)))

            elif etype == EventType.USER_REQUEST_INVALID:
                formatted = self._format_user_request_invalid(
                    event)
                if formatted:
                    for fl in formatted.split("\n"):
                        lines.append(self._box_line(fl))

            elif etype == EventType.WARNING:
                msg = event.get("message", "")
                lines.append(self._box_line(
                    "WARNING: %s" % msg))

            elif etype == EventType.ERROR:
                msg = event.get("message", "")
                lines.append(self._box_line(
                    "ERROR: %s" % msg))

            elif etype == EventType.SANITY_CHECK:
                if not event.get("passed"):
                    for flag in event.get("red_flags", []):
                        lines.append(self._box_line(
                            "  [RED FLAG] %s" % flag))
                    for warn in event.get("warnings", []):
                        lines.append(self._box_line(
                            "  [WARNING] %s" % warn))

            elif etype == EventType.METRICS_TREND:
                if event.get("plateau"):
                    cycles = event.get(
                        "cycles_analyzed", "?")
                    lines.append(self._box_line(
                        "  Trend: PLATEAU after %s cycles"
                        % cycles))
                elif not event.get("improving"):
                    lines.append(self._box_line(
                        "  Trend: Not improving"))
                # "Improving" is the happy path — skip it
                # to reduce noise.

            elif etype == EventType.NOTICE:
                msg = event.get("message", "")
                details = event.get("details", "")
                lines.append(self._box_line(
                    "NOTICE: %s" % msg))
                if details:
                    for dl in self._box_wrap_lines(
                            details, indent=2):
                        lines.append(dl)

            # Skip CYCLE_START, CYCLE_COMPLETE, DEBUG,
            # THOUGHT, FILE_SCORED in the box — they're
            # redundant with the header or too detailed.

        lines.append(self._box_bottom())
        return "\n".join(lines)

    def _format_metrics_line(self, event):
        """Compact one-line metrics for inside a box."""
        parts = []
        r_free = event.get("r_free")
        r_free_prev = event.get("r_free_prev")
        if r_free is not None:
            try:
                r_free = float(r_free)
                s = "R-free: %.4f" % r_free
                if r_free_prev is not None:
                    r_free_prev = float(r_free_prev)
                    diff = r_free - r_free_prev
                    arrow = "\u2193" if diff < 0 else (
                        "\u2191" if diff > 0 else "\u2192")
                    s += " %s (was %.4f)" % (arrow, r_free_prev)
                parts.append(s)
            except (ValueError, TypeError):
                parts.append("R-free: %s" % r_free)
        resolution = event.get("resolution")
        if resolution is not None:
            try:
                parts.append("Res: %.2f \u00c5" % float(resolution))
            except (ValueError, TypeError):
                pass
        cc = event.get("map_cc")
        if cc is not None:
            try:
                parts.append("CC: %.3f" % float(cc))
            except (ValueError, TypeError):
                pass
        if parts:
            return "  " + "  |  ".join(parts)
        return None

    def format_result_block(self, status, metrics=None,
                            error_text=None, program=None):
        """Format a post-execution result as a small boxed block.

        Called by ai_agent.py after a program finishes.

        Args:
            status: "SUCCESS" or "FAILED" or "ERROR"
            metrics: dict of metric name → value (optional)
            error_text: error message for failures (optional)
            program: program name for context (optional)

        Returns:
            str: Boxed result block

        Never raises — returns a minimal block on any error.
        """
        try:
            return self._format_result_block_inner(
                status, metrics, error_text, program)
        except Exception:
            return "%s  %s" % (self.TAG_RESULT, status or "?")

    def _format_result_block_inner(self, status, metrics,
                                   error_text, program):
        """Inner implementation of format_result_block."""
        if self.verbosity == Verbosity.QUIET:
            tag = self.TAG_RESULT
            s = "%s  %s" % (tag, status or "?")
            if metrics and metrics.get("r_free") is not None:
                s += "  R-free=%.4f" % metrics["r_free"]
            return s

        lines = []
        # Result header
        result_label = status or "?"
        if program:
            result_label = "%s: %s" % (program, status)
        lines.append(self._box_top("RESULT"))
        lines.append(self._box_line(
            "%s  %s" % (self.TAG_RESULT, result_label)))

        if status == "FAILED" or status == "ERROR":
            if error_text:
                # Show first 3 lines of error
                err_lines = str(error_text).strip().splitlines()
                for el in err_lines[:3]:
                    for wl in self._box_wrap_lines(
                            el.strip(), indent=2):
                        lines.append(wl)
                if len(err_lines) > 3:
                    lines.append(self._box_line(
                        "  ... (%d more lines)"
                        % (len(err_lines) - 3)))
            else:
                lines.append(self._box_line(
                    "  (no error details available)"))
        elif metrics:
            for key, val in metrics.items():
                if val is not None:
                    if isinstance(val, float):
                        lines.append(self._box_line(
                            "  %s: %.4f" % (key, val)))
                    else:
                        lines.append(self._box_line(
                            "  %s: %s" % (key, val)))

        lines.append(self._box_bottom())
        return "\n".join(lines)

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
