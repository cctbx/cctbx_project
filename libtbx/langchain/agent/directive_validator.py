"""
Directive Validator for PHENIX AI Agent.

This module validates and augments LLM planning decisions against stored
directives. It ensures that user instructions are consistently followed
across all cycles.

The validator:
1. Fills in missing strategy parameters from directives
2. Overrides LLM choices that contradict directives
3. Checks stop conditions
4. Logs modifications and warnings

Usage:
    from libtbx.langchain.agent.directive_validator import validate_intent

    result = validate_intent(
        intent=llm_intent,
        directives=session.get_directives(),
        cycle_number=3,
        history=session.get_cycles()
    )

    validated_intent = result["validated_intent"]
    if result["should_stop"]:
        # Handle stop condition
"""

from __future__ import absolute_import, division, print_function


def validate_intent(intent, directives, cycle_number, history=None, log_func=None):
    """
    Validate and modify LLM intent based on user directives.

    This function applies program-specific settings and constraints from
    user directives. It does NOT check stop conditions - that happens
    post-execution in ai_agent._run_single_cycle().

    Args:
        intent: LLM's planned action dict with keys:
            - program: Program name (e.g., "phenix.refine")
            - files: Dict of input files
            - strategy: Dict of strategy parameters
            - reasoning: LLM's reasoning string
            - stop: Boolean indicating if LLM wants to stop
        directives: Stored directive dict from session
        cycle_number: Current cycle number
        history: List of previous cycle dicts (optional)
        log_func: Optional logging function

    Returns:
        dict: {
            "validated_intent": Modified intent dict,
            "modifications": List of changes made,
            "warnings": List of potential issues
        }

    Note:
        Stop condition checking has been moved to post-execution.
        This function no longer returns should_stop or stop_reason.
    """
    def log(msg):
        if log_func:
            log_func(msg)

    if not intent:
        intent = {}

    if not directives:
        # No directives - return intent unchanged
        return {
            "validated_intent": intent,
            "modifications": [],
            "warnings": []
        }

    modifications = []
    warnings = []

    # Make a copy to avoid modifying original
    validated = dict(intent)
    validated["strategy"] = dict(intent.get("strategy") or {})
    validated["files"] = dict(intent.get("files") or {})

    program = validated.get("program")

    # 1. Apply program-specific settings
    if program and program != "STOP":
        setting_mods, setting_warns = _apply_program_settings(
            validated, directives, program, log
        )
        modifications.extend(setting_mods)
        warnings.extend(setting_warns)

    # 2. Check file preferences
    file_warns = _check_file_preferences(validated, directives, log)
    warnings.extend(file_warns)

    # 3. Check workflow preferences
    workflow_warns = _check_workflow_preferences(validated, directives, log)
    warnings.extend(workflow_warns)

    # NOTE: Stop condition checking has been removed from here.
    # Stop conditions are now checked ONLY in ai_agent._run_single_cycle()
    # after the program executes. This ensures clean separation of concerns:
    # - workflow_engine: determines what programs are valid
    # - directive_validator: applies program settings
    # - ai_agent: checks stop conditions after execution

    return {
        "validated_intent": validated,
        "modifications": modifications,
        "warnings": warnings
    }


def _apply_program_settings(intent, directives, program, log):
    """
    Apply directive settings to intent strategy.

    Fills in missing parameters from directives and warns/overrides
    when LLM contradicts directives.

    Args:
        intent: Intent dict (modified in place)
        directives: Directives dict
        program: Program name
        log: Logging function

    Returns:
        tuple: (modifications list, warnings list)
    """
    modifications = []
    warnings = []

    prog_settings = directives.get("program_settings", {})
    if not prog_settings:
        return modifications, warnings

    # Get default settings
    default_settings = prog_settings.get("default", {})

    # Get program-specific settings (override defaults)
    specific_settings = prog_settings.get(program, {})

    # Merge: specific overrides default
    directive_settings = {**default_settings, **specific_settings}

    if not directive_settings:
        return modifications, warnings

    strategy = intent.get("strategy", {})

    # Key aliases: map common LLM variations to canonical names
    # This handles cases where LLM uses different names than directives
    key_aliases = {
        # cycles variations
        "number_of_cycles": "cycles",
        "macro_cycles": "cycles",
        "num_cycles": "cycles",
        # resolution variations
        "high_resolution": "resolution",
        "d_min": "resolution",
        "resolution_limit": "resolution",
    }

    # Normalize strategy keys to canonical names
    normalized_strategy = {}
    for key, value in strategy.items():
        canonical_key = key_aliases.get(key, key)
        normalized_strategy[canonical_key] = value

    # Also normalize directive keys (in case directive uses non-canonical names)
    normalized_directives = {}
    for key, value in directive_settings.items():
        canonical_key = key_aliases.get(key, key)
        normalized_directives[canonical_key] = value

    for key, directive_value in normalized_directives.items():
        if key not in normalized_strategy:
            # LLM didn't specify - add from directive
            strategy[key] = directive_value
            modifications.append(
                "Added %s=%s from directive (program: %s)" % (key, directive_value, program)
            )
            log("VALIDATOR: Added %s=%s for %s" % (key, directive_value, program))

        elif normalized_strategy[key] != directive_value:
            # LLM specified different value - warn and override
            llm_value = normalized_strategy[key]
            warnings.append(
                "LLM used %s=%s but directive says %s=%s for %s" % (
                    key, llm_value, key, directive_value, program
                )
            )
            log("VALIDATOR WARNING: Overriding %s from %s to %s (directive)" % (
                key, llm_value, directive_value
            ))
            strategy[key] = directive_value
            modifications.append(
                "Override %s: %s -> %s (directive for %s)" % (
                    key, llm_value, directive_value, program
                )
            )

    # Remove aliased keys from strategy (keep only canonical names)
    for alias, canonical in key_aliases.items():
        if alias in strategy and canonical in strategy:
            del strategy[alias]  # Remove the alias, keep canonical

    intent["strategy"] = strategy
    return modifications, warnings


def _check_file_preferences(intent, directives, log):
    """
    Check if intent respects file preferences from directives.

    Currently only warns - doesn't modify file choices.

    Args:
        intent: Intent dict
        directives: Directives dict
        log: Logging function

    Returns:
        list: Warnings
    """
    warnings = []

    file_prefs = directives.get("file_preferences", {})
    if not file_prefs:
        return warnings

    intent_files = intent.get("files", {})

    # Check for excluded files
    excluded = file_prefs.get("exclude", [])
    for file_type, filepath in intent_files.items():
        if filepath:
            # Check if any excluded pattern matches
            for excl in excluded:
                if excl in str(filepath):
                    warnings.append(
                        "Using file %s which matches excluded pattern '%s'" % (filepath, excl)
                    )
                    log("VALIDATOR WARNING: File %s matches exclude pattern %s" % (filepath, excl))

    # Check preferred files
    for file_type in ("model", "sequence", "mtz"):
        preferred = file_prefs.get(file_type)
        if preferred and file_type in intent_files:
            actual = intent_files[file_type]
            if actual and preferred not in str(actual):
                warnings.append(
                    "Using %s=%s but directive prefers %s" % (file_type, actual, preferred)
                )

    return warnings


def _check_workflow_preferences(intent, directives, log):
    """
    Check if intent respects workflow preferences.

    Args:
        intent: Intent dict
        directives: Directives dict
        log: Logging function

    Returns:
        list: Warnings
    """
    warnings = []

    workflow_prefs = directives.get("workflow_preferences", {})
    if not workflow_prefs:
        return warnings

    program = intent.get("program")
    if not program:
        return warnings

    # Check skip_programs
    skip_programs = workflow_prefs.get("skip_programs", [])
    if program in skip_programs:
        warnings.append(
            "Running %s which is in skip_programs directive" % program
        )
        log("VALIDATOR WARNING: %s is in skip_programs but was selected" % program)

    return warnings


def _check_stop_conditions(directives, cycle_number, last_program, history, log):
    """
    Check if any stop conditions from directives are met.

    Args:
        directives: Directives dict
        cycle_number: Current cycle number
        last_program: Program that would run (or just ran)
        history: List of previous cycle dicts
        log: Logging function

    Returns:
        tuple: (should_stop: bool, reason: str or None)
    """
    stop_cond = directives.get("stop_conditions", {})
    if not stop_cond:
        return False, None

    # Check after_cycle
    if "after_cycle" in stop_cond:
        target_cycle = stop_cond["after_cycle"]
        if cycle_number >= target_cycle:
            reason = "Reached cycle %d (directive: stop after cycle %d)" % (
                cycle_number, target_cycle
            )
            log("VALIDATOR: Stop condition met - %s" % reason)
            return True, reason

    # Check after_program - this triggers AFTER the program completes
    # So we check if the program just ran in history
    if "after_program" in stop_cond and history:
        target_program = stop_cond["after_program"]
        # Check if last completed cycle was the target program
        if history:
            last_cycle = history[-1] if history else {}
            if last_cycle.get("program") == target_program:
                reason = "Completed %s (directive: stop after %s)" % (
                    target_program, target_program
                )
                log("VALIDATOR: Stop condition met - %s" % reason)
                return True, reason

    # Check r_free_target
    if "r_free_target" in stop_cond and history:
        target_r_free = stop_cond["r_free_target"]
        # Look for most recent r_free in history
        for cycle in reversed(history):
            metrics = cycle.get("metrics", {})
            r_free = metrics.get("r_free")
            if r_free is not None:
                if r_free <= target_r_free:
                    reason = "R-free %.3f reached target %.3f (directive)" % (
                        r_free, target_r_free
                    )
                    log("VALIDATOR: Stop condition met - %s" % reason)
                    return True, reason
                break  # Only check most recent r_free

    # Check map_cc_target
    if "map_cc_target" in stop_cond and history:
        target_cc = stop_cond["map_cc_target"]
        for cycle in reversed(history):
            metrics = cycle.get("metrics", {})
            map_cc = metrics.get("map_cc")
            if map_cc is not None:
                if map_cc >= target_cc:
                    reason = "Map CC %.3f reached target %.3f (directive)" % (
                        map_cc, target_cc
                    )
                    log("VALIDATOR: Stop condition met - %s" % reason)
                    return True, reason
                break

    return False, None


def _check_max_program_cycles(directives, program, history, log):
    """
    Check if max cycles for a specific program have been reached.

    Args:
        directives: Directives dict
        program: Program about to run
        history: List of previous cycle dicts
        log: Logging function

    Returns:
        tuple: (should_stop: bool, reason: str or None)
    """
    stop_cond = directives.get("stop_conditions", {})

    # Check max_refine_cycles for refinement programs
    if program in ("phenix.refine", "phenix.real_space_refine"):
        max_cycles = stop_cond.get("max_refine_cycles")
        if max_cycles is not None:
            # Count how many times refine has run
            refine_count = 0
            for cycle in (history or []):
                if cycle.get("program") in ("phenix.refine", "phenix.real_space_refine"):
                    refine_count += 1

            if refine_count >= max_cycles:
                reason = "Reached max refinement cycles (%d/%d) per directive" % (
                    refine_count, max_cycles
                )
                log("VALIDATOR: Stop condition met - %s" % reason)
                return True, reason

    return False, None


def augment_intent_with_directives(intent, directives, program=None):
    """
    Simple function to add directive settings to an intent.

    Use this when you just want to add missing settings without
    full validation.

    Args:
        intent: Intent dict (modified in place)
        directives: Directives dict
        program: Program name (uses intent["program"] if not provided)

    Returns:
        dict: Modified intent
    """
    if not directives:
        return intent

    if program is None:
        program = intent.get("program")

    if not program or program == "STOP":
        return intent

    prog_settings = directives.get("program_settings", {})
    default_settings = prog_settings.get("default", {})
    specific_settings = prog_settings.get(program, {})

    directive_settings = {**default_settings, **specific_settings}

    strategy = intent.get("strategy", {})
    for key, value in directive_settings.items():
        if key not in strategy:
            strategy[key] = value

    intent["strategy"] = strategy
    return intent


def get_stop_reason_from_directives(directives, cycle_number, last_program, history=None):
    """
    Get the stop reason if directives say to stop.

    Convenience function for checking stop conditions without full validation.

    Args:
        directives: Directives dict
        cycle_number: Current cycle number
        last_program: Last program that ran
        history: Previous cycle history

    Returns:
        str or None: Stop reason if should stop, None otherwise
    """
    should_stop, reason = _check_stop_conditions(
        directives, cycle_number, last_program, history, lambda x: None
    )
    return reason if should_stop else None


def format_validation_result(result):
    """
    Format validation result for display.

    Args:
        result: Result dict from validate_intent

    Returns:
        str: Formatted string
    """
    lines = []

    if result.get("modifications"):
        lines.append("Modifications applied:")
        for mod in result["modifications"]:
            lines.append("  - %s" % mod)

    if result.get("warnings"):
        lines.append("Warnings:")
        for warn in result["warnings"]:
            lines.append("  ! %s" % warn)

    if result.get("should_stop"):
        lines.append("Stop triggered: %s" % result.get("stop_reason", "Unknown"))

    return "\n".join(lines) if lines else "No changes from directives"
