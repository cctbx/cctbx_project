"""
LangGraph Node Functions for PHENIX AI Agent.

Nodes:
- perceive: Extract metrics from log, build metrics history
- plan: Decide next action (LLM), with auto-stop on plateau
- build: Construct command from intent (templates)
- validate: Validate command (files, workflow state, duplicates)
- fallback: Mechanical auto-selection when LLM fails
- output_node: Prepare response for client

NEW: Set USE_NEW_COMMAND_BUILDER = True to use unified CommandBuilder.
"""

from __future__ import absolute_import, division, print_function
import json
import re
import os

# Centralized pattern utilities
from libtbx.langchain.agent.pattern_manager import (
    patterns,
    extract_cycle_number,
    extract_all_numbers,
)

# Event logging for transparency
from libtbx.langchain.agent.event_log import EventType

# Import from libtbx.langchain
from libtbx.langchain.agent.template_builder import TemplateBuilder

from libtbx.langchain.agent.metrics_analyzer import (
    derive_metrics_from_history,
    analyze_metrics_trend,
    get_latest_resolution,
)

from libtbx.langchain.agent.workflow_state import (
    detect_workflow_state,
    validate_program_choice
)

from libtbx.langchain.knowledge.prompts_hybrid import get_planning_prompt

from libtbx.langchain.agent.rules_selector import select_action_by_rules

# Feature flag: Set to True to use new CommandBuilder
USE_NEW_COMMAND_BUILDER = True

def _get_command_builder():
    """Lazy-load CommandBuilder and CommandContext classes to avoid circular imports."""
    from libtbx.langchain.agent.command_builder import CommandBuilder, CommandContext
    return CommandBuilder, CommandContext


def _get_most_recent_file(file_list):
    """
    Get the most recent file from a list based on cycle numbering.

    For files with numbered suffixes (e.g., rsr_011_real_space_refined_000.pdb),
    prefers higher cycle numbers (the first number, e.g., 011).
    Falls back to modification time.

    Uses centralized cycle extraction from pattern_manager for consistency
    across all parts of the agent.

    Args:
        file_list: List of file paths

    Returns:
        str: Most recent file path, or None if list is empty
    """
    if not file_list:
        return None
    if len(file_list) == 1:
        return file_list[0]

    # Sort by cycle number (descending), then all numbers, then mtime
    def sort_key(path):
        # Primary: cycle number (for cross-cycle comparison)
        cycle = extract_cycle_number(path, default=0)
        # Secondary: all numbers as tuple (for within-cycle comparison)
        all_nums = extract_all_numbers(path)
        try:
            mtime = os.path.getmtime(path) if os.path.exists(path) else 0
        except OSError:
            mtime = 0
        return (cycle, all_nums, mtime)

    sorted_files = sorted(file_list, key=sort_key, reverse=True)
    return sorted_files[0]


# =============================================================================
# YAML MODE CONFIGURATION
# =============================================================================
# These are now True by default to use YAML-driven logic.
# Set to False to use legacy hardcoded logic.

USE_YAML_WORKFLOW = True     # Use workflows.yaml for state detection
USE_YAML_METRICS = True      # Use metrics.yaml for trend analysis
USE_RULES_SELECTOR = False   # Use rules-based selection (no LLM) - off by default

def set_yaml_mode(enabled=True):
    """
    Enable or disable YAML-driven mode globally.

    When enabled:
    - Workflow state detection uses workflows.yaml
    - Metrics trend analysis uses metrics.yaml thresholds

    Note: This does NOT affect LLM vs rules selection.
    Use set_rules_only_mode() for that.

    Args:
        enabled: True to enable YAML mode, False for legacy hardcoded mode
    """
    global USE_YAML_WORKFLOW, USE_YAML_METRICS
    USE_YAML_WORKFLOW = enabled
    USE_YAML_METRICS = enabled


def set_rules_only_mode(enabled=True):
    """
    Enable or disable rules-only mode (no LLM).

    When enabled, the agent uses deterministic rules from YAML
    to select programs instead of calling an LLM.

    Args:
        enabled: True for rules-only, False for LLM-based selection
    """
    global USE_RULES_SELECTOR
    USE_RULES_SELECTOR = enabled


builder = TemplateBuilder()

# Global LLM instance - initialized on first use
_llm = None
_llm_provider = None

# Supported providers
SUPPORTED_PROVIDERS = ["google", "openai", "ollama"]


# =============================================================================
# LLM MANAGEMENT
# =============================================================================

def validate_provider(provider):
    """
    Validate that the provider is supported and properly configured.

    Returns:
        tuple: (is_valid, error_message)
    """
    if provider not in SUPPORTED_PROVIDERS:
        return False, "Unknown provider '%s'. Supported: %s" % (
            provider, ", ".join(SUPPORTED_PROVIDERS)
        )

    # Check API keys
    if provider == "google":
        if not os.getenv("GOOGLE_API_KEY"):
            return False, (
                "Provider 'google' requires GOOGLE_API_KEY environment variable.\n"
                "Set it with: export GOOGLE_API_KEY='your-key-here'"
            )
    elif provider == "openai":
        if not os.getenv("OPENAI_API_KEY"):
            return False, (
                "Provider 'openai' requires OPENAI_API_KEY environment variable.\n"
                "Set it with: export OPENAI_API_KEY='your-key-here'"
            )
    elif provider == "ollama":
        # Check if Ollama is reachable
        ollama_url = os.getenv("OLLAMA_BASE_URL", "http://localhost:11434")
        try:
            import urllib.request
            req = urllib.request.Request(ollama_url, method='HEAD')
            req.timeout = 2
            urllib.request.urlopen(req, timeout=2)
        except Exception:
            return False, (
                "Provider 'ollama' requires Ollama server running at %s.\n"
                "Either:\n"
                "  1. Start Ollama: ollama serve\n"
                "  2. Set OLLAMA_BASE_URL to point to your server\n"
                "  3. Use a different provider: provider='google' or provider='openai'"
            ) % ollama_url

    return True, None


def get_planning_llm(provider=None):
    """
    Get or create the LLM for planning.

    Uses the existing get_llm_and_embeddings or get_expensive_llm from langchain_tools.
    The LLM is cached for reuse.

    Args:
        provider: LLM provider ("google", "openai", "ollama")

    Returns:
        tuple: (llm, error_message) - llm is None if there's an error
    """
    global _llm, _llm_provider

    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "google")

    # Validate provider first
    is_valid, error = validate_provider(provider)
    if not is_valid:
        return None, error

    # Return cached LLM if same provider
    if _llm is not None and _llm_provider == provider:
        return _llm, None

    try:
        from libtbx.langchain.core.llm import get_expensive_llm

        # For Ollama, we need json_mode for structured output
        if provider == 'ollama':
            _llm, _ = get_expensive_llm(
                provider=provider,
                timeout=120,
                json_mode=True
            )
        else:
            # Google and OpenAI handle JSON well without special mode
            _llm, _ = get_expensive_llm(
                provider=provider,
                timeout=120,
                json_mode=False
            )

        _llm_provider = provider
        return _llm, None

    except ImportError as e:
        return None, "Missing LLM package for provider '%s': %s" % (provider, str(e))
    except Exception as e:
        return None, "Failed to initialize LLM for provider '%s': %s" % (provider, str(e))


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _log(state, msg):
    """
    Helper to add debug message to state.

    Writes to both:
    - debug_log (list of strings) for backward compatibility
    - events (list of dicts) for new event system
    """
    debug_log = list(state.get("debug_log", []))
    debug_log.append(msg)

    # Also emit as DEBUG event for new system
    events = list(state.get("events", []))
    events.append({"type": EventType.DEBUG, "message": msg})

    return {**state, "debug_log": debug_log, "events": events}


def _emit(state, event_type, **kwargs):
    """
    Helper to emit a structured event to state.

    Args:
        state: Current state dict
        event_type: EventType constant
        **kwargs: Event-specific data

    Returns:
        Updated state with new event
    """
    events = list(state.get("events", []))
    events.append({"type": event_type, **kwargs})
    return {**state, "events": events}


def _increment_attempt(state):
    """Record failed attempt and increment counter."""
    attempt = {
        "intent": state.get("intent", {}),
        "command": state.get("command", ""),
        "error": state.get("validation_error", "")
    }
    previous = list(state.get("previous_attempts", []))

    return {
        **state,
        "attempt_number": state.get("attempt_number", 0) + 1,
        "previous_attempts": previous + [attempt]
    }


# =============================================================================
# NODE: PERCEIVE
# =============================================================================

def perceive(state):
    """
    Node A: Extract metrics from log text and build metrics history.

    This node:
    1. Derives metrics_history from the history backpack
    2. Extracts metrics from the current log
    3. Appends current metrics to history
    4. Analyzes trends for plateau detection
    5. Detects workflow state
    6. Runs sanity checks for red flag detection
    """
    # Initialize events list if not present
    if "events" not in state:
        state = {**state, "events": []}

    # Emit CYCLE_START event
    cycle_number = state.get("cycle_number", 1)
    state = _emit(state, EventType.CYCLE_START, cycle_number=cycle_number)

    log_text = state.get("log_text", "")
    history = state.get("history", [])

    # 1. Derive metrics from history
    metrics_history = derive_metrics_from_history(history)

    # 2. Extract metrics from current log
    try:
        from phenix.phenix_ai.log_parsers import extract_all_metrics
        analysis = extract_all_metrics(log_text)
    except ImportError:
        # Fallback if log_parsers not available
        analysis = _fallback_extract_metrics(log_text)

    # 2b. CRITICAL: Ensure we have the LAST r_free/r_work values
    # Refinement logs contain multiple R-free values (one per macro cycle).
    # Some extractors return the FIRST match, but we need the FINAL value.
    # Always re-extract using our robust last-match method to be safe.
    if log_text:
        last_r_free = patterns.extract_last_float("metrics.r_free", log_text, flags=re.IGNORECASE)
        last_r_work = patterns.extract_last_float("metrics.r_work", log_text, flags=re.IGNORECASE)

        if last_r_free is not None:
            # Only override if we found a value and it differs
            # (the last value should be <= first value after refinement converges)
            if analysis is None:
                analysis = {}
            analysis["r_free"] = last_r_free

        if last_r_work is not None:
            if analysis is None:
                analysis = {}
            analysis["r_work"] = last_r_work

    # 3. Append current cycle metrics if we have useful data
    if analysis and any(v is not None for k, v in analysis.items() if k != "program"):
        current_metrics = {
            "cycle": len(metrics_history) + 1,
            "program": analysis.get("program", "unknown"),
            "r_free": analysis.get("r_free"),
            "r_work": analysis.get("r_work"),
            "tfz": analysis.get("tfz"),
            "llg": analysis.get("llg"),
            "resolution": analysis.get("resolution"),
            "map_cc": analysis.get("map_cc"),
        }
        metrics_history.append(current_metrics)

    # 4. Analyze metrics trend
    resolution = get_latest_resolution(metrics_history)

    # Get available_files first (needed for both experiment type detection and workflow state)
    available_files = state.get("available_files", [])

    # Detect experiment type for proper trend analysis
    # PRIORITY: Use session's locked experiment_type first, then detect from files
    session_info = state.get("session_info", {})
    session_experiment_type = session_info.get("experiment_type")

    if session_experiment_type:
        # Use session's locked experiment type (already validated)
        experiment_type = session_experiment_type
        if experiment_type == "cryoEM":  # Normalize case
            experiment_type = "cryoem"
    else:
        # Fall back to file-based detection
        # DEFENSIVE: Filter out None values from available_files
        available_files_clean = [f for f in available_files if f is not None]
        if len(available_files_clean) != len(available_files):
            state = _log(state, "PERCEIVE: WARNING - Filtered %d None values from available_files" %
                        (len(available_files) - len(available_files_clean)))
            available_files = available_files_clean
        has_map = any(f.lower().endswith(('.mrc', '.ccp4', '.map')) for f in available_files)
        has_mtz = any(f.lower().endswith(('.mtz', '.sca', '.hkl')) for f in available_files)
        experiment_type = "cryoem" if has_map and not has_mtz else "xray"

    # Check for YAML mode (from state or global setting)
    use_yaml = state.get("use_yaml_mode", USE_YAML_METRICS)

    metrics_trend = analyze_metrics_trend(
        metrics_history, resolution, experiment_type,
        use_yaml_evaluator=use_yaml
    )

    # 5. Detect workflow state
    use_yaml_workflow = state.get("use_yaml_mode", USE_YAML_WORKFLOW)

    # Debug: check history for mtriage
    mtriage_in_history = any(
        "mtriage" in (h.get("program", "") + h.get("command", "")).lower()
        for h in history if isinstance(h, dict)
    )
    state = _log(state, "PERCEIVE: History has %d entries, mtriage_in_history=%s" % (
        len(history), mtriage_in_history))

    workflow_state = detect_workflow_state(
        history=history,
        available_files=available_files,
        analysis=analysis,
        maximum_automation=state.get("maximum_automation", True),
        use_yaml_engine=use_yaml_workflow,
        directives=state.get("directives", {})
    )

    # Debug categorization (helps diagnose model vs search_model issues)
    categorized = workflow_state.get("categorized_files", {})
    model_files = categorized.get("model", [])
    search_model_files = categorized.get("search_model", [])
    if search_model_files and not model_files:
        state = _log(state, "PERCEIVE: Has search_model but no model - sanity check may trigger")
        for f in search_model_files[:3]:
            state = _log(state, "PERCEIVE:   search_model: %s" % os.path.basename(f))

    # Override detected experiment_type with locked session experiment_type if available
    session_info = state.get("session_info", {})
    session_experiment_type = session_info.get("experiment_type")
    if session_experiment_type:
        state = _log(state, "PERCEIVE: Using locked experiment_type from session: %s (detected was: %s)" % (
            session_experiment_type, workflow_state.get("experiment_type")))
        workflow_state["experiment_type"] = session_experiment_type

    # === EMIT STATE_DETECTED EVENT ===
    state = _emit(state, EventType.STATE_DETECTED,
        workflow_state=workflow_state.get("state", "unknown"),
        experiment_type=workflow_state.get("experiment_type", experiment_type),
        reason=workflow_state.get("reason", ""),
        valid_programs=workflow_state.get("valid_programs", []))

    # === EMIT METRICS_EXTRACTED EVENT ===
    # Get previous metrics for comparison
    prev_metrics = metrics_history[-2] if len(metrics_history) >= 2 else None
    metrics_event_data = {}
    if analysis:
        for key in ["r_free", "r_work", "resolution", "map_cc", "tfz", "llg", "clashscore"]:
            if analysis.get(key) is not None:
                metrics_event_data[key] = analysis[key]
                # Add previous value if available
                if prev_metrics and prev_metrics.get(key) is not None:
                    metrics_event_data[key + "_prev"] = prev_metrics[key]
    if metrics_event_data:
        state = _emit(state, EventType.METRICS_EXTRACTED, **metrics_event_data)

    # Log the results
    if analysis:
        metrics_str = ", ".join([
            "%s=%s" % (k, v) for k, v in sorted(analysis.items())
            if k != "program" and v is not None
        ])
        state = _log(state, "PERCEIVE: Extracted metrics: %s" % metrics_str)
    else:
        state = _log(state, "PERCEIVE: No metrics extracted")

    state = _log(state, "PERCEIVE: Workflow state: %s" % workflow_state["state"])

    # Log valid_programs for debugging
    valid_progs = workflow_state.get("valid_programs", [])
    state = _log(state, "PERCEIVE: Valid programs: %s" % ", ".join(valid_progs))

    # Log forced_program if set
    forced_prog = workflow_state.get("forced_program")
    if forced_prog:
        state = _log(state, "PERCEIVE: Forced program (after_program directive): %s" % forced_prog)

    # Log file categorization for debugging
    categorized = workflow_state.get("categorized_files", {})
    if categorized:
        cat_summary = ", ".join([
            "%s=%d" % (k, len(v) if isinstance(v, list) else 1)
            for k, v in categorized.items() if v
        ])
        state = _log(state, "PERCEIVE: Categorized files: %s" % cat_summary)
    state = _log(state, "PERCEIVE: Available files count: %d" % len(available_files))

    # Log metrics trend details for debugging
    consecutive_rsr = metrics_trend.get("consecutive_rsr", 0)
    consecutive_refines = metrics_trend.get("consecutive_refines", 0)
    if consecutive_rsr > 0:
        state = _log(state, "PERCEIVE: Consecutive RSR cycles: %d (stop at 8)" % consecutive_rsr)
    if consecutive_refines > 0:
        state = _log(state, "PERCEIVE: Consecutive refine cycles: %d (stop at 8)" % consecutive_refines)

    trend_summary = metrics_trend.get("trend_summary", "")
    if trend_summary:
        state = _log(state, "PERCEIVE: Metrics trend: %s" % trend_summary)

    if metrics_trend.get("should_stop"):
        state = _log(state, "PERCEIVE: Auto-stop recommended: %s" % metrics_trend["reason"])

    # === EMIT METRICS_TREND EVENT ===
    state = _emit(state, EventType.METRICS_TREND,
        improving=metrics_trend.get("improving", False),
        plateau=metrics_trend.get("plateau_detected", metrics_trend.get("should_stop", False)),
        cycles_analyzed=len(metrics_history),
        trend_summary=trend_summary,
        should_stop=metrics_trend.get("should_stop", False),
        stop_reason=metrics_trend.get("reason", ""))

    # 6. Run sanity checks for red flag detection
    session_info = state.get("session_info", {})
    abort_on_red_flags = state.get("abort_on_red_flags", True)
    abort_on_warnings = state.get("abort_on_warnings", False)

    # Use session_info experiment_type (locked) first, then workflow_state, then detected
    session_experiment_type = session_info.get("experiment_type")
    effective_experiment_type = (
        session_experiment_type or
        workflow_state.get("experiment_type") or
        experiment_type
    )

    # Log best_files status and validate against available files
    best_files = session_info.get("best_files", {})
    if best_files:
        best_files_summary = ", ".join([
            "%s=%s" % (k, os.path.basename(v)) for k, v in best_files.items() if v
        ])
        state = _log(state, "PERCEIVE: Best files: %s" % best_files_summary)

    # VALIDATION: Check if best_files seem stale compared to available files
    # If we have a refined model in available_files but best_files shows phaser output,
    # we need to recalculate best_files on-the-fly
    if best_files:
        current_best_model = best_files.get("model")
        if current_best_model:
            current_best_basename = os.path.basename(current_best_model).lower()
            # Check if best model is a phaser/early-stage model
            if 'phaser' in current_best_basename or 'mr_' in current_best_basename:
                # Check if we have refined models in available_files
                for f in available_files:
                    basename = os.path.basename(f).lower()
                    # If we find a refined model, best_files is stale
                    if 'refine' in basename and basename.endswith('.pdb'):
                        state = _log(state, "PERCEIVE: WARNING - best_files may be stale (has phaser but refined exists: %s)" % os.path.basename(f))
                        # Recalculate best_files on the fly
                        try:
                            from libtbx.langchain.agent.best_files_tracker import BestFilesTracker
                            temp_tracker = BestFilesTracker()
                            for avail_f in available_files:
                                temp_tracker.evaluate_file(avail_f, cycle=1, metrics=None, stage=None)
                            new_best = temp_tracker.get_best_dict()
                            if new_best.get("model") != current_best_model:
                                state = _log(state, "PERCEIVE: Recalculated best_files: model=%s" % os.path.basename(new_best.get("model", "none")))
                                # Update session_info with recalculated best_files
                                session_info["best_files"] = new_best
                                state = {**state, "session_info": session_info}
                        except Exception as e:
                            state = _log(state, "PERCEIVE: Could not recalculate best_files: %s" % str(e))
                        break

    # Build sanity check context
    categorized_files = workflow_state.get("categorized_files", {})
    session_resolution = state.get("session_resolution") or resolution
    sanity_context = {
        "experiment_type": effective_experiment_type,
        # SEMANTIC: 'model' = positioned models (phaser_output, refined, docked)
        # SEMANTIC: 'search_model' = templates NOT yet positioned (predicted, pdb_template)
        "has_model": bool(categorized_files.get("model")),
        "has_search_model": bool(categorized_files.get("search_model")),
        "has_mtz": bool(categorized_files.get("mtz")),
        # Check all map categories for cryo-EM
        "has_map": bool(
            categorized_files.get("map") or
            categorized_files.get("full_map") or
            categorized_files.get("half_map")
        ),
        # Resolution is known if we have a value from mtriage/xtriage
        "has_resolution": session_resolution is not None,
        "state": workflow_state.get("state", ""),
        "history": history,
        "metrics_history": metrics_history,
        "resolution": session_resolution,
        "categorized_files": categorized_files,
    }

    # Debug: log key sanity context values
    state = _log(state, "PERCEIVE: Sanity context: exp_type=%s, session_exp_type=%s, has_mtz=%s, has_model=%s, has_search_model=%s, has_map=%s" % (
        sanity_context.get("experiment_type"),
        session_info.get("experiment_type"),
        sanity_context.get("has_mtz"),
        sanity_context.get("has_model"),
        sanity_context.get("has_search_model"),
        sanity_context.get("has_map"),
    ))

    # Run sanity checks
    try:
        from libtbx.langchain.agent.sanity_checker import SanityChecker
        checker = SanityChecker()
    except ImportError:
        checker = None

    if checker:
        sanity_result = checker.check(
            sanity_context,
            session_info,
            abort_on_red_flags=abort_on_red_flags,
            abort_on_warnings=abort_on_warnings
        )

        # === EMIT SANITY_CHECK EVENT ===
        red_flags = [i.message for i in sanity_result.issues if i.severity == "red_flag"]
        warnings = [i.message for i in sanity_result.issues if i.severity == "warning"]
        state = _emit(state, EventType.SANITY_CHECK,
            passed=not sanity_result.should_abort,
            red_flags=red_flags,
            warnings=warnings)

        # Log any warnings
        for issue in sanity_result.issues:
            if issue.severity == "warning":
                state = _log(state, "PERCEIVE: SANITY WARNING [%s]: %s" % (issue.code, issue.message))

        # Handle abort
        if sanity_result.should_abort:
            # Log the specific issues
            for issue in sanity_result.issues:
                state = _log(state, "PERCEIVE: RED FLAG [%s]: %s" % (issue.code, issue.message))
            state = _log(state, "PERCEIVE: SANITY ABORT - Stopping due to red flag(s)")

            # Build detailed abort message with diagnostic info
            abort_msg = sanity_result.abort_message
            diag_info = []
            diag_info.append("Diagnostic info:")
            diag_info.append("  - Available files: %d" % len(state.get("available_files", [])))
            diag_info.append("  - has_mtz: %s" % sanity_context.get("has_mtz"))
            diag_info.append("  - has_model: %s" % sanity_context.get("has_model"))
            diag_info.append("  - has_map: %s" % sanity_context.get("has_map"))
            diag_info.append("  - experiment_type: %s" % sanity_context.get("experiment_type"))
            diag_info.append("  - workflow_state: %s" % sanity_context.get("state"))
            cat_files = sanity_context.get("categorized_files", {})
            if cat_files:
                diag_info.append("  - categorized: %s" % ", ".join([
                    "%s=%d" % (k, len(v) if isinstance(v, list) else 1)
                    for k, v in cat_files.items() if v
                ]))
            else:
                diag_info.append("  - categorized: (empty)")
            abort_msg = abort_msg + "\n\n" + "\n".join(diag_info)

            return {
                **state,
                "log_analysis": analysis or {},
                "metrics_history": metrics_history,
                "metrics_trend": metrics_trend,
                "workflow_state": workflow_state,
                "command": "STOP",
                "stop": True,
                "stop_reason": "red_flag",
                "red_flag_issues": [i.to_dict() for i in sanity_result.issues],
                "abort_message": abort_msg,
            }

    return {
        **state,
        "log_analysis": analysis or {},
        "metrics_history": metrics_history,
        "metrics_trend": metrics_trend,
        "workflow_state": workflow_state,
    }


def _fallback_extract_metrics(log_text):
    """
    Fallback metric extraction using centralized patterns.

    This is used when the full log_parsers module is not available.
    It uses the PatternManager for consistent regex patterns.

    NOTE: This is a FALLBACK. The primary extraction should use
    phenix.phenix_ai.log_parsers.extract_all_metrics() which has
    program-specific parsing logic.
    """
    analysis = {}

    if not log_text:
        return analysis

    # Use centralized patterns for metric extraction
    # patterns is imported at module level

    # R-free and R-work: Use extract_last_float because refinement logs
    # contain multiple values (one per macro cycle) and we want the FINAL one
    r_free = patterns.extract_last_float("metrics.r_free", log_text, flags=re.IGNORECASE)
    if r_free is not None:
        analysis["r_free"] = r_free

    r_work = patterns.extract_last_float("metrics.r_work", log_text, flags=re.IGNORECASE)
    if r_work is not None:
        analysis["r_work"] = r_work

    tfz = patterns.extract_float("metrics.tfz_score", log_text)
    if tfz is not None:
        analysis["tfz"] = tfz

    llg = patterns.extract_float("metrics.llg", log_text)
    if llg is not None:
        analysis["llg"] = llg

    map_cc = patterns.extract_float("metrics.map_cc", log_text, flags=re.IGNORECASE)
    if map_cc is not None:
        analysis["map_cc"] = map_cc

    resolution = patterns.extract_float("metrics.resolution_from_range", log_text, flags=re.IGNORECASE)
    if resolution is None:
        resolution = patterns.extract_float("metrics.resolution", log_text, flags=re.IGNORECASE)
    if resolution is not None:
        analysis["resolution"] = resolution

    # Error detection using patterns
    if patterns.has_pattern("system.phenix_failed"):
        error_match = patterns.search("system.phenix_failed", log_text)
        if error_match:
            analysis["error"] = error_match.group(1).strip()

    if "error" not in analysis and patterns.has_pattern("system.phenix_sorry"):
        sorry_match = patterns.search("system.phenix_sorry", log_text)
        if sorry_match:
            analysis["error"] = sorry_match.group(1).strip()

    # Program detection (simple string matching, not regex-dependent)
    log_lower = log_text.lower()
    if "phenix.real_space_refine" in log_lower:
        analysis["program"] = "phenix.real_space_refine"
    elif "phenix.refine" in log_lower:
        analysis["program"] = "phenix.refine"
    elif "phenix.phaser" in log_lower or "SOLU SET" in log_text:
        analysis["program"] = "phenix.phaser"
    elif "phenix.xtriage" in log_lower:
        analysis["program"] = "phenix.xtriage"
    elif "phenix.mtriage" in log_lower:
        analysis["program"] = "phenix.mtriage"
    elif "predict_and_build" in log_lower:
        analysis["program"] = "phenix.predict_and_build"
    elif "phenix.ligandfit" in log_lower:
        analysis["program"] = "phenix.ligandfit"

    return analysis


# =============================================================================
# NODE: PLAN
# =============================================================================

def parse_intent_json(llm_output):
    """Clean and parse LLM JSON output."""
    text = llm_output.strip()

    # Remove markdown code blocks if present
    if "```json" in text:
        text = text.split("```json")[1].split("```")[0]
    elif "```" in text:
        parts = text.split("```")
        if len(parts) >= 2:
            text = parts[1]

    # Remove any leading/trailing whitespace
    text = text.strip()

    try:
        intent = json.loads(text)
    except json.JSONDecodeError:
        # Fallback: try to find JSON object boundaries
        start = text.find("{")
        end = text.rfind("}") + 1
        if start >= 0 and end > start:
            intent = json.loads(text[start:end])
        else:
            raise ValueError("Could not parse JSON from LLM response: %s..." % text[:100])

    # Sanitize intent fields to ensure correct types
    # This prevents crashes when LLM returns unexpected types
    if not isinstance(intent.get("files"), dict):
        intent["files"] = {}
    if not isinstance(intent.get("strategy"), dict):
        intent["strategy"] = {}
    if not isinstance(intent.get("program"), (str, type(None))):
        intent["program"] = str(intent.get("program")) if intent.get("program") else None
    if not isinstance(intent.get("reasoning"), str):
        intent["reasoning"] = str(intent.get("reasoning", ""))

    return intent


def plan(state):
    """
    Node B: Decide Intent (LLM), with auto-stop on plateau.

    This node:
    0. Checks for forced retry (from error recovery) - bypasses all other logic
    1. Checks for auto-stop conditions (plateau, success)
       BUT: Suppresses auto-stop if there's an after_program directive that hasn't been fulfilled
    2. Optionally uses rules-only mode (no LLM)
    3. Constructs prompt with workflow state and metrics trend
    4. Calls LLM to decide next action
    5. Falls back to rules selector if LLM unavailable
    """
    # 0. Check for forced retry (from error recovery)
    # This takes priority over everything else - when recovery is triggered,
    # we want to immediately retry the failed program with recovery flags
    session_info = state.get("session_info", {})
    force_retry = session_info.get("force_retry_program")

    # Debug: log what we received
    if session_info.get("recovery_strategies"):
        state = _log(state, f"PLAN: Recovery strategies present for files: {list(session_info['recovery_strategies'].keys())}")

    if force_retry:
        state = _log(state, f"PLAN: Forced retry of {force_retry} (error recovery)")

        # Emit event for transparency
        state = _emit(state, EventType.PROGRAM_SELECTED,
                      program=force_retry,
                      reason="error_recovery_retry",
                      forced=True)

        return {
            **state,
            "intent": {
                "program": force_retry,
                "reasoning": "Retrying with recovery flags after recoverable error",
                "strategy": {},  # Strategy will be merged from recovery_strategies in BUILD
                "files": {},
                "forced_retry": True
            },
            # Signal to clear the flag after this cycle
            "clear_force_retry": True
        }

    # 1. Check for auto-stop from metrics trend
    metrics_trend = state.get("metrics_trend", {})

    if metrics_trend.get("should_stop"):
        reason = metrics_trend.get("reason", "Metrics indicate completion")

        # BUT: Check if there's an after_program directive that hasn't been fulfilled
        # In that case, we should NOT auto-stop until that program runs
        directives = state.get("directives", {})
        stop_cond = directives.get("stop_conditions", {})
        after_program = stop_cond.get("after_program")

        if after_program:
            # Check if the after_program has been run in history
            history = state.get("history", [])
            after_program_done = any(
                after_program.replace("phenix.", "") in
                (h.get("program", "") + h.get("command", "")).lower().replace("phenix.", "")
                for h in history if isinstance(h, dict)
            )

            if not after_program_done:
                # Suppress auto-stop - user wants a specific program to run first
                state = _log(state, "PLAN: Suppressing AUTO-STOP because after_program=%s hasn't run yet" % after_program)
                # Continue to LLM planning instead of stopping
            else:
                # after_program has run, allow auto-stop
                state = _log(state, "PLAN: AUTO-STOP triggered: %s" % reason)
                # Emit STOP_DECISION event
                state = _emit(state, EventType.STOP_DECISION,
                    stop=True,
                    reason=reason)
                return {
                    **state,
                    "intent": {
                        "program": None,
                        "stop": True,
                        "stop_reason": reason,
                        "reasoning": "Automatically stopping: %s. %s" % (
                            reason, metrics_trend.get("trend_summary", "")
                        ),
                        "files": {},
                        "strategy": {},
                    },
                    "stop": True,
                    "stop_reason": reason,
                }
        else:
            # No after_program directive, allow auto-stop
            state = _log(state, "PLAN: AUTO-STOP triggered: %s" % reason)
            # Emit STOP_DECISION event
            state = _emit(state, EventType.STOP_DECISION,
                stop=True,
                reason=reason)
            return {
                **state,
                "intent": {
                    "program": None,
                    "stop": True,
                    "stop_reason": reason,
                    "reasoning": "Automatically stopping: %s. %s" % (
                        reason, metrics_trend.get("trend_summary", "")
                    ),
                    "files": {},
                    "strategy": {},
                },
                "stop": True,
                "stop_reason": reason,
            }

    # 2. Check for rules-only mode (no LLM)
    use_rules_only = state.get("use_rules_only", False) or state.get("use_yaml_mode", USE_RULES_SELECTOR)
    if use_rules_only:
        state = _log(state, "PLAN: Rules-only mode enabled, skipping LLM")
        return _mock_plan(state)

    # 3. Get provider from state
    provider = state.get("provider")
    if not provider:
        provider = os.getenv("LLM_PROVIDER", "google")

    # 3. Construct Prompt with workflow and metrics context
    workflow_state = state.get("workflow_state", {})

    # Pre-filter valid programs based on user advice
    # This ensures the LLM only sees programs the user wants
    user_advice = state.get("user_advice", "")
    if user_advice and workflow_state.get("valid_programs"):
        from libtbx.langchain.agent.rules_selector import RulesSelector
        selector = RulesSelector()
        filtered_programs = selector._apply_user_advice(
            workflow_state["valid_programs"], user_advice
        )
        if filtered_programs != workflow_state["valid_programs"]:
            state = _log(state, "PLAN: User advice filtered programs to: %s" % filtered_programs)
            # Create a copy to avoid modifying original
            workflow_state = dict(workflow_state)
            workflow_state["valid_programs"] = filtered_programs

    # Get directives from state
    directives = state.get("directives", {})

    # Get best_files from session_info
    session_info = state.get("session_info", {})
    best_files = session_info.get("best_files", {})

    system_msg, user_msg = get_planning_prompt(
        history=state.get("history", []),
        analysis=state.get("log_analysis", {}),
        available_files=state.get("available_files", []),
        previous_attempts=state.get("previous_attempts", []),
        user_advice=user_advice,
        metrics_trend=metrics_trend,
        workflow_state=workflow_state,
        directives=directives,
        best_files=best_files
    )

    # 4. Get LLM (with validation)
    llm, error = get_planning_llm(provider)

    if llm is None:
        state = _log(state, "PLAN: LLM error - %s" % error)
        state = _log(state, "PLAN: Falling back to mock planner")
        return _mock_plan(state)

    # 5. Call LLM with rate limit handling
    try:
        state = _log(state, "PLAN: Calling LLM (provider=%s)..." % provider)

        from langchain_core.messages import SystemMessage, HumanMessage

        messages = [
            SystemMessage(content=system_msg),
            HumanMessage(content=user_msg)
        ]

        # Try to use rate limit handler with pre-configured settings
        handler = None
        try:
            from libtbx.langchain.agent.rate_limit_handler import (
                get_google_handler, get_openai_handler, get_anthropic_handler
            )
            if provider == "google":
                handler = get_google_handler()
            elif provider == "openai":
                handler = get_openai_handler()
            elif provider == "anthropic":
                handler = get_anthropic_handler()
        except ImportError:
            try:
                # Try relative import for standalone testing
                from agent.rate_limit_handler import (
                    get_google_handler, get_openai_handler, get_anthropic_handler
                )
                if provider == "google":
                    handler = get_google_handler()
                elif provider == "openai":
                    handler = get_openai_handler()
                elif provider == "anthropic":
                    handler = get_anthropic_handler()
            except ImportError:
                pass  # No rate limit handler available

        if handler:
            # Create a list to capture state updates from log wrapper
            log_messages = []
            def log_wrapper(msg):
                log_messages.append(msg)

            def make_call():
                return llm.invoke(messages)

            response = handler.call_with_retry(make_call, log_wrapper)

            # Add any rate limit log messages to state
            for msg in log_messages:
                state = _log(state, "PLAN: %s" % msg)
        else:
            # No rate limit handler available, call directly
            response = llm.invoke(messages)

        content = response.content

        state = _log(state, "PLAN: LLM response received (%d chars)" % len(content))

    except Exception as e:
        state = _log(state, "PLAN: LLM call failed - %s" % str(e))
        return _mock_plan(state)

    # 6. Parse Response
    try:
        intent = parse_intent_json(content)

        chosen_program = intent.get("program")
        state = _log(state, "PLAN: LLM selected program=%s" % chosen_program)

        # === VALIDATE INTENT AGAINST USER DIRECTIVES ===
        directives = state.get("directives", {})
        if directives:
            try:
                from libtbx.langchain.agent.directive_validator import validate_intent

                cycle_number = state.get("cycle_number", 1)
                history = state.get("history", [])
                attempt_number = state.get("attempt_number", 0)

                validation_result = validate_intent(
                    intent=intent,
                    directives=directives,
                    cycle_number=cycle_number,
                    history=history,
                    attempt_number=attempt_number,
                    log_func=lambda msg: None  # Suppress internal logs, we'll log summary
                )

                # Apply validated intent
                intent = validation_result["validated_intent"]

                # Log modifications
                for mod in validation_result.get("modifications", []):
                    state = _log(state, "PLAN: DIRECTIVE - %s" % mod)

                # Log warnings
                for warn in validation_result.get("warnings", []):
                    state = _log(state, "PLAN: DIRECTIVE WARNING - %s" % warn)

                # NOTE: Stop condition checking removed from here.
                # Stop conditions are checked post-execution in ai_agent._run_single_cycle()

                # Update chosen_program in case it was modified
                chosen_program = intent.get("program")

            except ImportError:
                state = _log(state, "PLAN: directive_validator not available, skipping validation")
            except Exception as e:
                state = _log(state, "PLAN: Directive validation error: %s" % str(e))

        # === VALIDATE LLM CHOICE AGAINST VALID_PROGRAMS ===
        # This is a simple check - all directive logic is handled upstream
        # by workflow_engine._apply_directives() which modifies valid_programs
        valid_programs = workflow_state.get("valid_programs", [])

        # Normalize: treat STOP request same as choosing "STOP" program
        # BUT: Only if the LLM didn't select a valid non-STOP program
        # This prevents cases where LLM says {program: "phenix.autosol", stop: true}
        # from incorrectly stopping instead of running the program
        if chosen_program == "STOP" or (intent.get("stop") and chosen_program not in valid_programs):
            chosen_program = "STOP"
            intent["program"] = "STOP"
        elif intent.get("stop") and chosen_program in valid_programs and chosen_program != "STOP":
            # LLM chose a valid program but also set stop=true
            # Trust the program choice - the stop might have been set for "after this program"
            state = _log(state, "PLAN: LLM set stop=true but also chose valid program %s, running program" % chosen_program)
            intent["stop"] = False  # Clear the stop flag since we have a valid program to run

        # Simple validation: is the chosen program in valid_programs?
        if chosen_program is None or chosen_program not in valid_programs:
            if valid_programs:
                # Override to first valid program (or first non-STOP if available)
                override_program = next((p for p in valid_programs if p != "STOP"), valid_programs[0])
                original_reasoning = intent.get("reasoning", "")

                state = _log(state, "PLAN: LLM choice '%s' not in valid_programs: %s" % (
                    chosen_program or "none", valid_programs))
                state = _log(state, "PLAN: Selecting %s instead" % override_program)

                # Check if this was a user-requested program (from user_advice)
                user_advice = state.get("user_advice", "")
                user_requested = False
                if chosen_program and user_advice:
                    # Check if user explicitly mentioned this program
                    if chosen_program.lower() in user_advice.lower():
                        user_requested = True
                    # Also check for partial matches like "run refine" -> "phenix.refine"
                    program_short = chosen_program.replace("phenix.", "")
                    if program_short.lower() in user_advice.lower():
                        user_requested = True

                # Emit USER_REQUEST_INVALID event if user asked for this
                if user_requested:
                    # Determine why it's invalid
                    all_known_programs = [
                        "phenix.xtriage", "phenix.mtriage", "phenix.predict_and_build",
                        "phenix.phaser", "phenix.refine", "phenix.real_space_refine",
                        "phenix.ligandfit", "phenix.pdbtools", "phenix.dock_in_map",
                        "phenix.autobuild", "phenix.autosol", "phenix.process_predicted_model",
                        "phenix.molprobity", "phenix.resolve_cryo_em", "phenix.map_sharpening",
                        "phenix.map_to_model", "phenix.map_symmetry", "phenix.validation_cryoem",
                    ]

                    if chosen_program in all_known_programs:
                        reason = "Not valid in current workflow state '%s'" % workflow_state.get("state", "unknown")
                        suggestion = "This program requires different conditions. %s" % workflow_state.get("reason", "")
                    else:
                        reason = "Unknown program (not supported by PHENIX AI Agent)"
                        suggestion = "Check the program name or use one of the available programs"

                    state = _emit(state, EventType.USER_REQUEST_INVALID,
                        requested_program=chosen_program,
                        reason=reason,
                        selected_instead=override_program,
                        valid_programs=valid_programs,
                        suggestion=suggestion)

                intent["program"] = override_program
                intent["stop"] = False
                intent["stop_reason"] = None
                intent["reasoning"] = "Selected %s ('%s' not available). [Original: %s]" % (
                    override_program,
                    chosen_program or "no choice",
                    original_reasoning[:200] if original_reasoning else "none"
                )
            else:
                # No valid programs - shouldn't happen but handle gracefully
                state = _log(state, "PLAN: ERROR - no valid programs available")
                intent["stop"] = True
                intent["stop_reason"] = "No valid programs available"

        # If STOP is the chosen program and it's valid, set stop flags
        elif chosen_program == "STOP":
            state = _log(state, "PLAN: STOP is valid, allowing stop")
            intent["stop"] = True
            if not intent.get("stop_reason"):
                intent["stop_reason"] = "Agent decided to stop"

        # === EMIT PROGRAM_SELECTED EVENT ===
        final_program = intent.get("program")
        state = _emit(state, EventType.PROGRAM_SELECTED,
            program=final_program,
            reasoning=intent.get("reasoning", ""),
            source="llm",
            provider=provider)

        # Check if program was modified from LLM choice
        if chosen_program != final_program and chosen_program is not None:
            state = _emit(state, EventType.PROGRAM_MODIFIED,
                original=chosen_program,
                selected=final_program,
                reason="Not in valid_programs")

        # Emit STOP_DECISION if stopping
        if intent.get("stop"):
            state = _emit(state, EventType.STOP_DECISION,
                stop=True,
                reason=intent.get("stop_reason", "Agent decided to stop"))

        return {
            **state,
            "intent": intent,
            "stop": intent.get("stop", False),
            "stop_reason": intent.get("stop_reason")
        }

    except Exception as e:
        state = _log(state, "PLAN: Error parsing LLM output - %s" % str(e))
        state = _log(state, "PLAN: Raw response: %s" % content[:200])
        return {
            **state,
            "validation_error": "LLM JSON Parse Error: %s" % str(e)
        }


def _mock_plan(state):
    """
    Fallback planner when LLM is unavailable.

    Uses RulesSelector if available for intelligent selection,
    otherwise falls back to simple first-valid-program logic.
    """
    state = _log(state, "PLAN: Using fallback planner (no LLM)")

    workflow_state = state.get("workflow_state", {})
    valid_programs = workflow_state.get("valid_programs", [])
    metrics_trend = state.get("metrics_trend", {})

    # Use rules-based selection
    try:
        state = _log(state, "PLAN: Using RulesSelector for intelligent selection")

        # Get categorized files from workflow state
        files = workflow_state.get("categorized_files", {})

        # Build history info from state
        history_info = _extract_history_info(state)

        # Get log analysis for error handling
        log_analysis = state.get("log_analysis", {})

        # Get user advice
        user_advice = state.get("user_advice", "")

        intent = select_action_by_rules(
            workflow_state=workflow_state,
            files=files,
            metrics_trend=metrics_trend,
            history_info=history_info,
            user_advice=user_advice,
            log_analysis=log_analysis
        )

        state = _log(state, "PLAN: RulesSelector chose %s" % intent.get("program"))

        # === EMIT PROGRAM_SELECTED EVENT ===
        state = _emit(state, EventType.PROGRAM_SELECTED,
            program=intent.get("program"),
            reasoning=intent.get("reasoning", ""),
            source="rules")

        if intent.get("stop"):
            state = _emit(state, EventType.STOP_DECISION,
                stop=True,
                reason=intent.get("stop_reason", "Rules-based decision"))

        return {
            **state,
            "intent": intent,
            "stop": intent.get("stop", False),
            "stop_reason": intent.get("stop_reason")
        }
    except Exception as e:
        state = _log(state, "PLAN: RulesSelector failed: %s, using simple fallback" % e)

    # Simple fallback if rules selector fails
    if not valid_programs:
        state = _log(state, "PLAN: No workflow state, cannot determine valid programs")
        return {
            **state,
            "intent": {
                "program": None,
                "reasoning": "Fallback: No workflow state available, stopping.",
                "files": {},
                "strategy": {},
                "stop": True,
                "stop_reason": "no_workflow_state"
            },
            "command": "STOP",
            "stop_reason": "no_workflow_state"
        }

    # Pick first valid program (excluding STOP unless it's the only option)
    if len(valid_programs) == 1:
        program = valid_programs[0]
    else:
        program = next((p for p in valid_programs if p != "STOP"), valid_programs[0])

    if program == "STOP":
        mock_intent = {
            "program": None,
            "reasoning": "Fallback: No valid programs available, stopping.",
            "files": {},
            "strategy": {},
            "stop": True,
            "stop_reason": "no_valid_programs"
        }
    else:
        mock_intent = {
            "program": program,
            "reasoning": "Fallback: LLM unavailable, using first valid program from workflow state.",
            "files": {},  # Empty = let builder auto-select
            "strategy": {},
            "confidence": "low",
            "stop": False,
            "stop_reason": None
        }

    # === EMIT PROGRAM_SELECTED EVENT ===
    final_program = mock_intent.get("program")
    state = _emit(state, EventType.PROGRAM_SELECTED,
        program=final_program,
        reasoning=mock_intent.get("reasoning", ""),
        source="rules" if "RulesSelector" in str(state.get("debug_log", [])) else "fallback")

    if mock_intent.get("stop"):
        state = _emit(state, EventType.STOP_DECISION,
            stop=True,
            reason=mock_intent.get("stop_reason", "Fallback decision"))

    return {
        **state,
        "intent": mock_intent,
        "stop": mock_intent.get("stop", False),
        "stop_reason": mock_intent.get("stop_reason")
    }


def _extract_history_info(state):
    """Extract history info from state for rules selector."""
    history = state.get("history", [])
    history_info = {
        "refine_count": 0,
        "rsr_count": 0,
        "twin_law": None,
        "has_ncs": False,
        "cycle_number": state.get("cycle_number", 1),
    }

    for entry in history:
        if isinstance(entry, dict):
            prog = (entry.get("program") or "").lower()
            if "refine" in prog and "real_space" not in prog:
                history_info["refine_count"] += 1
            elif "real_space" in prog:
                history_info["rsr_count"] += 1

            # Check for twin law in analysis
            analysis = entry.get("analysis", {})
            if isinstance(analysis, dict):
                if analysis.get("twin_law"):
                    history_info["twin_law"] = analysis["twin_law"]
                if analysis.get("ncs_found"):
                    history_info["has_ncs"] = True

    return history_info


# =============================================================================
# NODE: BUILD
# =============================================================================

def _build_with_new_builder(state):
    """
    Build command using the new unified CommandBuilder.

    This is a simpler implementation that delegates all complexity
    to CommandBuilder, which handles:
    - File selection (with priorities)
    - Strategy building (with output_prefix)
    - Invariant application (resolution, R-free flags, etc.)
    - Command assembly
    """
    intent = state.get("intent", {})

    if intent.get("stop"):
        state = _log(state, "BUILD: Stop requested, no command needed")
        return {**state, "command": "STOP"}

    program = intent.get("program", "phenix.refine")
    available_files = state.get("available_files", [])

    # Get the CommandBuilder and CommandContext classes
    CommandBuilder, CommandContext = _get_command_builder()

    # Sanitize LLM strategy
    raw_strategy = intent.get("strategy", {})
    if not isinstance(raw_strategy, dict):
        state = _log(state, "BUILD: LLM returned invalid strategy type, using empty dict")
        raw_strategy = {}

    # Log strategy for debugging
    if raw_strategy:
        strategy_summary = ", ".join(["%s=%s" % (k, v) for k, v in raw_strategy.items()])
        state = _log(state, "BUILD: Strategy from intent: {%s}" % strategy_summary)
    else:
        state = _log(state, "BUILD: No strategy in intent")

    # Correct LLM file paths (map basenames to actual paths)
    llm_files = intent.get("files", {})
    corrected_files = {}
    basename_to_path = {os.path.basename(f): f for f in available_files}

    if llm_files:
        for slot, filepath in llm_files.items():
            if filepath:
                if isinstance(filepath, list):
                    corrected = []
                    for fp in filepath:
                        basename = os.path.basename(fp) if fp else ""
                        if basename in basename_to_path:
                            corrected.append(basename_to_path[basename])
                    if corrected:
                        corrected_files[slot] = corrected
                else:
                    basename = os.path.basename(filepath)
                    if basename in basename_to_path:
                        corrected_files[slot] = basename_to_path[basename]

    # NOTE: Half-map vs full-map categorization is handled by _categorize_files()
    # in workflow_state.py. If a single map was promoted to full_map, trust that decision.

    # Build context from state
    session_info = state.get("session_info", {})
    workflow_state = state.get("workflow_state", {})

    # Handle stepwise cryo-EM
    strategy = dict(raw_strategy)
    if (workflow_state.get("automation_path") == "stepwise" and
        program == "phenix.predict_and_build" and
        workflow_state.get("state") == "cryoem_analyzed"):
        strategy["stop_after_predict"] = True
        state = _log(state, "BUILD: Forcing stop_after_predict=True for stepwise cryo-EM")

    # Debug: log recovery strategies if present
    recovery_strategies = session_info.get("recovery_strategies", {})
    if recovery_strategies:
        state = _log(state, f"BUILD: Recovery strategies available for: {list(recovery_strategies.keys())}")

    # Create log wrapper that updates state
    build_logs = []
    def log_wrapper(msg):
        build_logs.append(msg)

    # Create CommandContext
    context = CommandContext(
        cycle_number=state.get("cycle_number", 1),
        experiment_type=session_info.get("experiment_type", ""),
        resolution=state.get("session_resolution") or workflow_state.get("resolution"),
        best_files=session_info.get("best_files", {}),
        rfree_mtz=session_info.get("rfree_mtz"),
        categorized_files=workflow_state.get("categorized_files", {}),
        workflow_state=workflow_state.get("state", ""),
        history=state.get("history", []),
        llm_files=corrected_files if corrected_files else None,
        llm_strategy=strategy if strategy else None,
        recovery_strategies=recovery_strategies,
        log=log_wrapper,
    )

    # Build command
    builder = CommandBuilder()
    command = builder.build(program, available_files, context)

    # Add builder logs to state
    for log_msg in build_logs:
        state = _log(state, log_msg)

    if not command:
        state = _log(state, "BUILD: Failed to build command for %s" % program)
        state = _emit(state, EventType.ERROR, message="Failed to build command", details=program)
        return {**state, "command": "", "validation_error": "Failed to build command"}

    state = _log(state, "BUILD: Generated command: %s" % command[:150])

    # === EMIT FILES_SELECTED EVENT (verbose level) ===
    selection_details = builder.get_selection_details()
    if selection_details:
        state = _emit(state, EventType.FILES_SELECTED, selections=selection_details)

    # === EMIT COMMAND_BUILT EVENT ===
    state = _emit(state, EventType.COMMAND_BUILT, command=command, program=program)

    return {
        **state,
        "command": command,
        "program": program,
        "corrected_files": corrected_files,
        "strategy": strategy,
    }


def build(state):
    """
    Node C: Build Command from Intent using Templates.

    This is 100% deterministic - no LLM involved.
    The template builder guarantees correct syntax.

    Applies Tier 2 defaults from workflow_state and logs any overrides.

    NEW: If USE_NEW_COMMAND_BUILDER is True, delegates to _build_with_new_builder.
    """
    # Optionally use new unified builder
    if USE_NEW_COMMAND_BUILDER:
        return _build_with_new_builder(state)

    intent = state.get("intent", {})

    if intent.get("stop"):
        state = _log(state, "BUILD: Stop requested, no command needed")
        return {**state, "command": "STOP"}

    program = intent.get("program", "phenix.refine")

    # Handle special strategy for stepwise cryo-EM
    # IMPORTANT: Safely handle case where LLM returns non-dict for strategy
    raw_strategy = intent.get("strategy", {})
    if not isinstance(raw_strategy, dict):
        state = _log(state, "BUILD: LLM returned invalid strategy type (%s), using empty dict" % type(raw_strategy).__name__)
        raw_strategy = {}
    strategy = dict(raw_strategy)

    workflow_state = state.get("workflow_state", {})

    if (workflow_state.get("automation_path") == "stepwise" and
        program == "phenix.predict_and_build" and
        workflow_state.get("state") == "cryoem_analyzed"):
        # Force stop_after_predict for stepwise mode
        strategy["stop_after_predict"] = True
        state = _log(state, "BUILD: Forcing stop_after_predict=True for stepwise cryo-EM")

    # === APPLY TIER 2 DEFAULTS AND LOG OVERRIDES ===
    recommended_strategy = workflow_state.get("recommended_strategy", {})

    if recommended_strategy and program == "phenix.refine":
        for key, info in recommended_strategy.items():
            recommended_value = info.get("value")
            reason = info.get("reason", "")
            override_warning = info.get("override_warning", "")

            if key in strategy:
                # LLM specified a value - check if it differs from recommendation
                llm_value = strategy[key]
                if llm_value != recommended_value:
                    # Log the override with warning
                    state = _log(state, "BUILD: OVERRIDE - %s=%s (recommended: %s)" % (
                        key, llm_value, recommended_value))
                    if override_warning:
                        state = _log(state, "BUILD: Warning - %s" % override_warning)
                else:
                    state = _log(state, "BUILD: %s=%s (matches recommendation)" % (key, llm_value))
            else:
                # LLM didn't specify - apply the default
                if recommended_value is not None:
                    strategy[key] = recommended_value
                    state = _log(state, "BUILD: Applying default %s=%s (%s)" % (
                        key, recommended_value, reason[:60] if reason else "recommended"))

    # === LEGACY: APPLY REFINE_HINTS (backward compatibility) ===
    # This handles cases where recommended_strategy might not be populated
    if program == "phenix.refine" and not recommended_strategy:
        refine_hints = workflow_state.get("refine_hints", [])

        # R-free flags: only generate on first refinement
        if "generate_rfree_flags" in refine_hints:
            if "generate_rfree_flags" not in strategy:
                strategy["generate_rfree_flags"] = True
                state = _log(state, "BUILD: First refinement - will generate R-free flags")
        elif "use_existing_rfree_flags" in refine_hints:
            if strategy.get("generate_rfree_flags"):
                state = _log(state, "BUILD: Warning - LLM requested new R-free flags but this is not first refinement")

        # Water building
        if "add_waters" in refine_hints:
            if "add_waters" not in strategy:
                strategy["add_waters"] = True
                state = _log(state, "BUILD: Adding waters (model ready, resolution ok)")

    # === TWINNING (always check, even with new structure) ===
    if program == "phenix.refine":
        if workflow_state.get("has_twinning") and workflow_state.get("twin_law"):
            if "twin_law" not in strategy:
                strategy["twin_law"] = workflow_state["twin_law"]
                state = _log(state, "BUILD: Adding twin law %s (twinning detected)" % workflow_state["twin_law"])

    # === OUTPUT PREFIX (prevent long filename accumulation) ===
    # Use descriptive prefixes: refine_001, refine_002 for X-ray; rsr_001, rsr_002 for cryo-EM
    # This makes output progression clear: refine_001_001.pdb → refine_002_001.pdb
    if program in ("phenix.refine", "phenix.real_space_refine"):
        if "output_prefix" not in strategy:
            history = state.get("history", [])

            # Count successful runs of THIS specific program
            if program == "phenix.refine":
                prefix_base = "refine"
                matching_runs = [h for h in history
                               if isinstance(h, dict) and
                               h.get("program") == "phenix.refine" and
                               str(h.get("result", "")).startswith("SUCCESS")]
            else:  # phenix.real_space_refine
                prefix_base = "rsr"
                matching_runs = [h for h in history
                               if isinstance(h, dict) and
                               h.get("program") == "phenix.real_space_refine" and
                               str(h.get("result", "")).startswith("SUCCESS")]

            run_count = len(matching_runs)

            # Debug: show which cycles were counted
            if matching_runs:
                cycle_nums = [str(h.get("cycle_number", "?")) for h in matching_runs]
                state = _log(state, "BUILD: Found %d previous %s runs in cycles: %s" % (
                    run_count, program, ", ".join(cycle_nums)))

            next_number = run_count + 1
            strategy["output_prefix"] = "%s_%03d" % (prefix_base, next_number)
            state = _log(state, "BUILD: output_prefix=%s_%03d" % (prefix_base, next_number))

    # === HIGH RESOLUTION SUGGESTIONS (just log, don't auto-apply) ===
    if program == "phenix.refine" and workflow_state.get("high_res_suggestions"):
        for suggestion in workflow_state["high_res_suggestions"]:
            if suggestion == "consider_anisotropic_adp":
                state = _log(state, "BUILD: Note - high resolution data, consider anisotropic_adp=true")
            elif suggestion == "consider_riding_hydrogens":
                state = _log(state, "BUILD: Note - very high resolution, consider adding hydrogens")

    # === APPLY AUTOBUILD RESOLUTION LIMIT ===
    if program == "phenix.autobuild":
        resolution = workflow_state.get("resolution")
        if resolution and resolution < 2.0:
            if "resolution" not in strategy:
                strategy["resolution"] = 2.0
                state = _log(state, "BUILD: Limiting autobuild resolution to 2.0Å (data is %.1fÅ)" % resolution)

    # === APPLY RESOLUTION FOR PROGRAMS THAT NEED IT ===
    # Helper to find resolution from various sources (priority order)
    def find_resolution():
        # 1. From session_resolution (single source of truth - set by xtriage/mtriage)
        session_res = state.get("session_resolution")
        if session_res:
            return session_res, "session"

        # 2. From workflow_state (extracted from history analysis)
        res = workflow_state.get("resolution")
        if res:
            return res, "workflow_state"

        # 3. From current log_analysis
        analysis = state.get("log_analysis", {})
        if analysis.get("resolution"):
            return analysis["resolution"], "log_analysis"

        # 4. From previous commands in history (last resort)
        history = state.get("history", [])
        import re
        for entry in reversed(history):
            cmd = entry.get("command", "") if isinstance(entry, dict) else ""
            res_match = re.search(r'resolution[=\s]+([\d\.]+)', cmd, re.IGNORECASE)
            if res_match:
                try:
                    res = float(res_match.group(1))
                    if 0.5 < res < 20:  # Sanity check
                        return res, "previous_command"
                except ValueError:
                    pass

        return None, None

    # NOTE: Resolution requirements for predict_and_build, real_space_refine,
    # and dock_in_map are now handled by invariants in programs.yaml.
    # The validate_and_fix() call below will auto-fill resolution from context.

    # Correct file paths from LLM - map basenames to actual paths
    available_files = state.get("available_files", [])
    basename_to_path = {os.path.basename(f): f for f in available_files}
    available_set = set(available_files) | set(basename_to_path.keys())

    llm_files = intent.get("files", {})
    corrected_files = {}
    rejected_files = []

    # DEBUG: Log what the LLM actually returned in its files dict
    if llm_files:
        llm_files_summary = ", ".join("%s=%s" % (k, os.path.basename(str(v)) if v else "None")
                                       for k, v in llm_files.items())
        state = _log(state, "BUILD: LLM requested files: {%s}" % llm_files_summary)
    else:
        state = _log(state, "BUILD: LLM returned empty files dict - will auto-select")

    if llm_files:
        for slot, filepath in llm_files.items():
            if filepath:
                # Handle case where LLM returns a list (e.g., for half_map)
                if isinstance(filepath, list):
                    # Preserve list - correct each path in the list
                    corrected_list = []
                    for fp in filepath:
                        if fp:
                            basename = os.path.basename(fp)
                            if basename in basename_to_path:
                                corrected_list.append(basename_to_path[basename])
                            elif fp in available_set:
                                corrected_list.append(fp)
                            else:
                                # File not in available_files - REJECT
                                rejected_files.append("%s=%s" % (slot, basename))
                    if corrected_list:
                        corrected_files[slot] = corrected_list
                else:
                    # Single file
                    basename = os.path.basename(filepath)
                    if basename in basename_to_path:
                        # Use the correct path from available_files
                        corrected_files[slot] = basename_to_path[basename]
                    elif filepath in available_set:
                        corrected_files[slot] = filepath
                    else:
                        # File not in available_files - REJECT
                        rejected_files.append("%s=%s" % (slot, basename))

    # If LLM tried to use files not in available_files, log and reject
    if rejected_files:
        error_msg = "LLM selected files not in available_files: %s" % ", ".join(rejected_files)
        state = _log(state, "BUILD: FILE REJECTION - %s" % error_msg)
        # Don't fail completely - let auto-fill try to find alternatives
        state = _log(state, "BUILD: Will attempt to auto-fill valid files instead")

    # NOTE: Half-map vs full-map categorization is handled by _categorize_files()
    # in workflow_state.py. If a single map was promoted to full_map, trust that decision.

    # === AUTO-FILL MODEL FOR AUTOBUILD ===
    # If LLM selected autobuild but didn't include model, find the best refined model
    if program == "phenix.autobuild" and "model" not in corrected_files:
        # Find refined PDB files (prefer most recent by cycle number)
        import re
        refined_pdbs = []
        for f in available_files:
            basename = os.path.basename(f).lower()
            if f.lower().endswith('.pdb') and 'refine' in basename:
                # Extract cycle number if present
                cycle_match = re.search(r'refine[_.]?(\d+)', basename)
                cycle_num = int(cycle_match.group(1)) if cycle_match else 0
                refined_pdbs.append((f, cycle_num))

        if refined_pdbs:
            # Sort by cycle number descending, take the most recent
            refined_pdbs.sort(key=lambda x: x[1], reverse=True)
            best_model = refined_pdbs[0][0]
            corrected_files["model"] = best_model
            state = _log(state, "BUILD: Auto-added model=%s for autobuild" % os.path.basename(best_model))
        else:
            # Try phaser output
            for f in available_files:
                basename = os.path.basename(f).lower()
                if f.lower().endswith('.pdb') and 'phaser' in basename:
                    corrected_files["model"] = f
                    state = _log(state, "BUILD: Auto-added model=%s for autobuild" % os.path.basename(f))
                    break

    # === OVERRIDE MODEL FOR REFINEMENT PROGRAMS ===
    # For refine/real_space_refine, ensure we use the MOST RECENT refined model
    # The LLM often picks older models incorrectly
    # Priority order is now defined in programs.yaml input_priorities
    if program in ("phenix.refine", "phenix.real_space_refine"):
        workflow_state = state.get("workflow_state", {})
        categorized_files = workflow_state.get("categorized_files", {})

        # Get input priorities from YAML
        priorities = builder._registry.get_input_priorities(program, "model")
        priority_categories = priorities.get("categories", [])
        exclude_categories = priorities.get("exclude_categories", [])

        # Find best model from categories (priority order from YAML)
        best_model = None
        best_source = None

        if priority_categories:
            for cat in priority_categories:
                cat_files = categorized_files.get(cat, [])
                if cat_files:
                    # Check exclusions
                    for f in cat_files:
                        excluded = False
                        for excl_cat in exclude_categories:
                            if f in categorized_files.get(excl_cat, []):
                                excluded = True
                                break
                        if not excluded:
                            best_model = f
                            best_source = cat.replace("_", " ")
                            break
                if best_model:
                    break
        else:
            # Fallback: use hardcoded priorities if YAML not defined
            if program == "phenix.real_space_refine":
                # Cryo-EM: prefer rsr_output > docked > with_ligand
                for cat, label in [("rsr_output", "RSR output"),
                                   ("docked", "docked model"),
                                   ("with_ligand", "model with ligand")]:
                    if categorized_files.get(cat):
                        best_model = categorized_files[cat][0]
                        best_source = label
                        break
            else:
                # X-ray: prefer refined > phaser_output > with_ligand
                for cat, label in [("refined", "refined model"),
                                   ("phaser_output", "Phaser output"),
                                   ("with_ligand", "model with ligand")]:
                    if categorized_files.get(cat):
                        best_model = categorized_files[cat][0]
                        best_source = label
                        break

        # Override LLM's choice if we found a better model
        if best_model:
            llm_model = corrected_files.get("model")
            if llm_model != best_model:
                corrected_files["model"] = best_model
                if llm_model:
                    state = _log(state, "BUILD: Overriding LLM model choice with %s (%s)" % (
                        os.path.basename(best_model), best_source))
                else:
                    state = _log(state, "BUILD: Auto-filled model=%s (%s)" % (
                        os.path.basename(best_model), best_source))

    try:
        if not corrected_files:
            # LLM didn't pick files, use auto-select with categorized files
            workflow_state = state.get("workflow_state", {})
            categorized_files = workflow_state.get("categorized_files", {})
            # Get best_files and rfree_mtz from session_info
            session_info = state.get("session_info", {})
            best_files = session_info.get("best_files", {})
            rfree_mtz = session_info.get("rfree_mtz")
            command = builder.build_command_for_program(
                program,
                available_files,
                categorized_files=categorized_files,
                best_files=best_files,
                rfree_mtz=rfree_mtz
            )
            state = _log(state, "BUILD: Auto-selected files for %s" % program)
        else:
            # LLM picked some files - auto-fill any missing required files
            workflow_state = state.get("workflow_state", {})
            categorized_files = workflow_state.get("categorized_files", {})
            # Get best_files and rfree_mtz from session_info
            session_info = state.get("session_info", {})
            best_files = session_info.get("best_files", {})
            rfree_mtz = session_info.get("rfree_mtz")

            # Get required inputs for this program
            try:
                required_inputs = builder._registry.get_required_inputs(program)
            except Exception:
                required_inputs = []

            # Map input names to best_files categories
            input_to_best_category = {
                "model": "model",
                "pdb_file": "model",
                "map": "map",
                "full_map": "map",
                "mtz": "mtz",
                "hkl_file": "mtz",
                "sequence": "sequence",
                "seq_file": "sequence",
            }

            # Auto-fill missing required files
            # PRIORITY 0: Locked R-free MTZ (for mtz input, X-ray refinement)
            # PRIORITY 1: Check best_files first
            # PRIORITY 2: Use YAML input_priorities
            # PRIORITY 3: Simple fallback
            for input_name in required_inputs:
                if input_name not in corrected_files:
                    file_found = None

                    # PRIORITY 0: Locked R-free MTZ
                    if input_name in ("mtz", "hkl_file") and rfree_mtz:
                        if os.path.exists(rfree_mtz):
                            file_found = rfree_mtz
                            state = _log(state, "BUILD: Using LOCKED R-free MTZ for %s: %s" % (
                                input_name, os.path.basename(rfree_mtz)))

                    # PRIORITY 1: Check best_files
                    if not file_found:
                        best_category = input_to_best_category.get(input_name, input_name)
                        best_path = best_files.get(best_category)
                        if best_path and os.path.exists(best_path):
                            file_found = best_path
                            state = _log(state, "BUILD: Using best_%s for %s: %s" % (
                                best_category, input_name, os.path.basename(best_path)))

                    # PRIORITY 2: Try to get priorities from YAML
                    if not file_found:
                        priorities = builder._registry.get_input_priorities(program, input_name)
                        priority_categories = priorities.get("categories", [])
                        exclude_categories = priorities.get("exclude_categories", [])

                        if priority_categories:
                            # Use YAML-defined priorities
                            for cat in priority_categories:
                                cat_files = categorized_files.get(cat, [])
                                for f in cat_files:
                                    # Check exclusions
                                    excluded = False
                                    for excl_cat in exclude_categories:
                                        if f in categorized_files.get(excl_cat, []):
                                            excluded = True
                                            break
                                    if not excluded:
                                        file_found = f
                                        break
                                if file_found:
                                    break

                    # PRIORITY 3: Simple fallback for common input types
                    # Use _get_most_recent_file to prefer the latest file when multiple match
                    if not file_found:
                        if input_name == "mtz":
                            if categorized_files.get("refined_mtz"):
                                file_found = _get_most_recent_file(categorized_files["refined_mtz"])
                            elif categorized_files.get("mtz"):
                                file_found = _get_most_recent_file(categorized_files["mtz"])
                        elif input_name == "model":
                            for cat in ["refined", "rsr_output", "phaser_output",
                                       "processed_predicted", "pdb"]:
                                if categorized_files.get(cat):
                                    file_found = _get_most_recent_file(categorized_files[cat])
                                    break
                        elif input_name == "sequence" and categorized_files.get("sequence"):
                            file_found = _get_most_recent_file(categorized_files["sequence"])
                        elif input_name == "map" and categorized_files.get("full_map"):
                            file_found = _get_most_recent_file(categorized_files["full_map"])

                    if file_found:
                        corrected_files[input_name] = file_found
                        state = _log(state, "BUILD: Auto-filled %s=%s" % (
                            input_name, os.path.basename(file_found)))

            # Validate invariants and apply fixes (single place for all program constraints)
            # Build context for auto-fills (resolution, etc.)
            # Use find_resolution() to get resolution from all possible sources
            found_resolution, res_source = find_resolution()
            invariant_context = {
                "session_resolution": state.get("session_resolution"),
                "resolution": found_resolution,  # From find_resolution() priority order
                "resolution_source": res_source,
                "workflow_state": workflow_state,
            }
            corrected_files, strategy, warnings = builder.validate_and_fix(
                program, corrected_files, strategy,
                log=lambda msg: None,  # Could log to state if needed
                context=invariant_context
            )
            for warning in warnings:
                state = _log(state, "BUILD: %s" % warning)

            # DEBUG: Log final file selection before building command
            if corrected_files:
                final_files_summary = ", ".join("%s=%s" % (k, os.path.basename(str(v)) if v else "None")
                                                 for k, v in corrected_files.items())
                state = _log(state, "BUILD: Final files for command: {%s}" % final_files_summary)

            command = builder.build_command(
                program,
                corrected_files,
                strategy
            )
            state = _log(state, "BUILD: Used LLM-selected files (paths corrected)")

        if not command:
            error_msg = "Builder returned empty command for %s" % program
            state = _log(state, "BUILD: Failed - %s" % error_msg)
            state = _emit(state, EventType.ERROR, message=error_msg)
            return {**state, "command": "", "validation_error": error_msg}

    except Exception as e:
        state = _log(state, "BUILD: Exception - %s" % str(e))
        state = _emit(state, EventType.ERROR, message="Build exception", details=str(e))
        return {**state, "command": "", "validation_error": str(e)}

    state = _log(state, "BUILD: Command = %s" % command[:80])

    # === EMIT COMMAND_BUILT EVENT ===
    state = _emit(state, EventType.COMMAND_BUILT,
        command=command,
        program=program)

    return {**state, "command": command, "validation_error": None}


# =============================================================================
# NODE: VALIDATE
# =============================================================================

def validate(state):
    """
    Node D: Validate command before execution.

    Checks:
    1. Command is not empty
    2. Program is valid for current workflow state (STRICT)
    3. Referenced files exist in available_files
    4. Command is not a duplicate of previous cycle
    """
    command = state.get("command", "")

    if not command:
        # Already has validation_error from Build
        return _increment_attempt(state)

    if command == "STOP":
        return {**state, "validation_error": None}

    intent = state.get("intent", {})
    chosen_program = intent.get("program")

    # === 1. WORKFLOW STATE VALIDATION (STRICT) ===
    workflow_state = state.get("workflow_state", {})

    if workflow_state and chosen_program:
        is_valid, error = validate_program_choice(chosen_program, workflow_state)
        if not is_valid:
            state = _log(state, "VALIDATE: WORKFLOW VIOLATION - %s" % error)
            state = {**state, "validation_error": error}
            return _increment_attempt(state)

    # === 2. FILE VALIDATION ===
    # On the server, we only have the file list passed from client - no filesystem access
    # So we check that referenced files are in available_files (by path or basename)
    available_files = state.get("available_files", [])
    available_set = set(available_files)
    available_basenames = {os.path.basename(f): f for f in available_set}

    # Debug: log file counts for troubleshooting
    state = _log(state, "VALIDATE: Checking against %d available files" % len(available_files))

    # Extract filenames from command
    file_pattern = r'[\w\-\.\/]+\.(?:pdb|mtz|fa|fasta|seq|cif|mrc|ccp4|map)\b'
    referenced_files = re.findall(file_pattern, command, re.IGNORECASE)

    missing = []
    for f in referenced_files:
        basename = os.path.basename(f)

        # Check if file is available (by exact path OR by basename)
        if f in available_set:
            # Exact path match - OK
            continue
        elif basename in available_basenames:
            # Basename match - OK (paths may differ between client/server)
            continue
        else:
            # Neither path nor basename found
            missing.append(f)

    if missing:
        error = "Files not found in available_files: %s" % ", ".join(missing)
        state = _log(state, "VALIDATE: FILE ERROR - %s" % error)
        state = {**state, "validation_error": error}
        return _increment_attempt(state)

    # === 3. DUPLICATE CHECK ===
    for hist in state.get("history", []):
        if isinstance(hist, dict):
            prev_command = hist.get("command", "")
            prev_cycle = hist.get("cycle_number", "?")
        elif isinstance(hist, str):
            prev_command = hist
            prev_cycle = "?"
        else:
            continue

        if prev_command and prev_command == command:
            error = "Duplicate of previous cycle %s" % prev_cycle
            state = _log(state, "VALIDATE: DUPLICATE - %s" % error)
            state = {**state, "validation_error": error}
            return _increment_attempt(state)

    state = _log(state, "VALIDATE: Passed")
    return {**state, "validation_error": None}


# =============================================================================
# NODE: FALLBACK
# =============================================================================

def _fallback_with_new_builder(state):
    """
    Fallback using the new unified CommandBuilder.

    Simpler implementation that delegates to CommandBuilder.
    """
    state = _log(state, "FALLBACK: Using mechanical auto-selection (new builder)")

    # Get CommandBuilder and CommandContext
    CommandBuilder, CommandContext = _get_command_builder()

    # Get previous commands to avoid duplicates
    previous_commands = set()
    for hist in state.get("history", []):
        if isinstance(hist, dict) and hist.get("command"):
            previous_commands.add(hist["command"])

    # Use workflow state to get valid programs
    workflow_state = state.get("workflow_state", {})
    valid_programs = workflow_state.get("valid_programs", [])

    # Filter out STOP unless it's the only option
    runnable = [p for p in valid_programs if p != "STOP"]

    if not runnable:
        state = _log(state, "FALLBACK: No valid programs, stopping")
        return {
            **state,
            "command": "STOP",
            "stop": True,
            "stop_reason": "no_valid_programs",
            "validation_error": None,
            "fallback_used": True
        }

    # Respect start_with_program directive - prioritize it
    directives = state.get("directives", {})
    stop_cond = directives.get("stop_conditions", {})
    start_with = stop_cond.get("start_with_program")
    if start_with and start_with in runnable:
        # Move start_with_program to front of list
        runnable = [start_with] + [p for p in runnable if p != start_with]
        state = _log(state, "FALLBACK: Prioritizing start_with_program: %s" % start_with)

    # Get available files
    available_files = list(state.get("available_files", []))
    session_info = state.get("session_info", {})

    # Create CommandContext (no LLM hints - pure rule-based)
    context = CommandContext(
        cycle_number=state.get("cycle_number", 1),
        experiment_type=session_info.get("experiment_type", ""),
        resolution=state.get("session_resolution") or workflow_state.get("resolution"),
        best_files=session_info.get("best_files", {}),
        rfree_mtz=session_info.get("rfree_mtz"),
        categorized_files=workflow_state.get("categorized_files", {}),
        workflow_state=workflow_state.get("state", ""),
        history=state.get("history", []),
        llm_files=None,  # No LLM hints
        llm_strategy=None,  # No LLM hints
        log=lambda msg: None,
    )

    builder = CommandBuilder()

    # Try each valid program until one works without duplicating
    for program in runnable:
        cmd = builder.build(program, available_files, context)

        if cmd and cmd not in previous_commands:
            state = _log(state, "FALLBACK: Built command for %s" % program)
            return {
                **state,
                "command": cmd,
                "validation_error": None,
                "fallback_used": True
            }
        elif cmd:
            state = _log(state, "FALLBACK: %s would duplicate, trying next" % program)
        else:
            state = _log(state, "FALLBACK: %s failed to build, trying next" % program)

    # Last resort - stop
    state = _log(state, "FALLBACK: Cannot find non-duplicate command, stopping")
    return {
        **state,
        "command": "STOP",
        "stop": True,
        "stop_reason": "all_commands_duplicate",
        "validation_error": None,
        "fallback_used": True
    }


def fallback(state):
    """
    Node E: Mechanical Fallback.

    Used when LLM fails 3 times. Uses workflow state to pick
    a valid program, then auto-selects files.

    Handles duplicate avoidance by:
    1. Checking previous commands in history
    2. Checking recently run programs (not just commands)
    3. Trying alternative programs if first choice would duplicate
    4. For refine, preferring output files from previous cycles

    NEW: If USE_NEW_COMMAND_BUILDER is True, delegates to _fallback_with_new_builder.
    """
    # Optionally use new unified builder
    if USE_NEW_COMMAND_BUILDER:
        return _fallback_with_new_builder(state)

    state = _log(state, "FALLBACK: Using mechanical auto-selection")

    # Get previous commands to avoid duplicates
    previous_commands = set()
    recent_programs = []  # Track recently run programs (last 3)
    for hist in state.get("history", []):
        if isinstance(hist, dict):
            if hist.get("command"):
                previous_commands.add(hist["command"])
            if hist.get("program"):
                recent_programs.append(hist["program"])

    # Get last 3 programs to avoid immediate repeats
    recent_programs = recent_programs[-3:] if recent_programs else []

    # Use workflow state to get valid programs
    workflow_state = state.get("workflow_state", {})
    valid_programs = workflow_state.get("valid_programs", [])

    # Filter out STOP unless it's the only option
    runnable = [p for p in valid_programs if p != "STOP"]

    # Deprioritize recently run programs (move to end of list)
    if recent_programs:
        # Sort: programs not recently run first
        runnable_prioritized = []
        runnable_deprioritized = []
        for p in runnable:
            if p in recent_programs:
                runnable_deprioritized.append(p)
            else:
                runnable_prioritized.append(p)
        runnable = runnable_prioritized + runnable_deprioritized

        if runnable_deprioritized:
            state = _log(state, "FALLBACK: Deprioritized recently run: %s" % runnable_deprioritized)

    # Respect start_with_program directive - prioritize it (overrides recent program logic)
    directives = state.get("directives", {})
    stop_cond = directives.get("stop_conditions", {})
    start_with = stop_cond.get("start_with_program")
    if start_with and start_with in runnable:
        # Move start_with_program to front of list
        runnable = [start_with] + [p for p in runnable if p != start_with]
        state = _log(state, "FALLBACK: Prioritizing start_with_program: %s" % start_with)

    if not runnable:
        state = _log(state, "FALLBACK: No valid programs, stopping")
        return {
            **state,
            "command": "STOP",
            "stop": True,
            "stop_reason": "no_valid_programs",
            "validation_error": None,
            "fallback_used": True
        }

    # Get available files, prioritizing recent outputs
    available_files = list(state.get("available_files", []))

    # Find output files from history (most recent first)
    output_files = []
    for hist in reversed(state.get("history", [])):
        if isinstance(hist, dict):
            hist_outputs = hist.get("output_files", [])
            for f in hist_outputs:
                if f and f not in output_files:
                    output_files.append(f)

    # Prioritize output files for refinement (use most recent refined model)
    prioritized_files = available_files[:]
    for f in output_files:
        if f not in prioritized_files:
            prioritized_files.insert(0, f)  # Add at front for priority

    # Build context for invariant validation (resolution, etc.)
    # This ensures programs like real_space_refine get resolution auto-filled
    fallback_context = {
        "session_resolution": state.get("session_resolution"),
        "resolution": state.get("session_resolution"),  # Use session resolution
        "workflow_state": workflow_state,
    }

    # Get best_files and rfree_mtz from session_info
    session_info = state.get("session_info", {})
    best_files = session_info.get("best_files", {})
    rfree_mtz = session_info.get("rfree_mtz")

    # Try each valid program until one works without duplicating
    for program in runnable:
        # For refine, try with prioritized files (recent outputs first)
        if "refine" in program.lower():
            cmd = builder.build_command_for_program(program, prioritized_files,
                                                    context=fallback_context,
                                                    best_files=best_files,
                                                    rfree_mtz=rfree_mtz)
        else:
            cmd = builder.build_command_for_program(program, available_files,
                                                    context=fallback_context,
                                                    best_files=best_files,
                                                    rfree_mtz=rfree_mtz)

        if cmd and cmd not in previous_commands:
            state = _log(state, "FALLBACK: Built command for %s" % program)
            return {
                **state,
                "command": cmd,
                "validation_error": None,
                "fallback_used": True
            }
        elif cmd:
            state = _log(state, "FALLBACK: %s would duplicate, trying next" % program)
        else:
            state = _log(state, "FALLBACK: %s failed to build, trying next" % program)

    # If all valid programs would duplicate, try autobuild ONLY if it's valid for this state
    sequence_files = [f for f in available_files if f.endswith(('.fa', '.fasta', '.seq'))]
    if sequence_files and "phenix.autobuild" in valid_programs and "phenix.autobuild" not in runnable:
        cmd = builder.build_command_for_program("phenix.autobuild", available_files,
                                                context=fallback_context,
                                                best_files=best_files,
                                                rfree_mtz=rfree_mtz)
        if cmd and cmd not in previous_commands:
            state = _log(state, "FALLBACK: Trying autobuild as alternative")
            return {
                **state,
                "command": cmd,
                "validation_error": None,
                "fallback_used": True
            }

    # If validation is required, run molprobity
    if "phenix.molprobity" in valid_programs:
        # Find a recent PDB file for molprobity
        pdb_files = [f for f in prioritized_files if f.endswith('.pdb')]
        if pdb_files:
            cmd = "phenix.molprobity %s" % pdb_files[0]
            if cmd not in previous_commands:
                state = _log(state, "FALLBACK: Running validation (molprobity)")
                return {
                    **state,
                    "command": cmd,
                    "validation_error": None,
                    "fallback_used": True
                }

    # Last resort: try basic refine with output files
    state = _log(state, "FALLBACK: All options exhausted, trying phenix.refine with outputs")

    # Find refined PDBs from history outputs
    refined_pdbs = [f for f in output_files if f.endswith('.pdb') and 'refine' in f.lower()]
    if refined_pdbs:
        mtz_files = [f for f in available_files if f.endswith(".mtz")]
        if mtz_files:
            cmd = "phenix.refine %s %s" % (refined_pdbs[0], mtz_files[0])
            if cmd not in previous_commands:
                state = _log(state, "FALLBACK: Using refined output from previous cycle")
                return {
                    **state,
                    "command": cmd,
                    "validation_error": None,
                    "fallback_used": True
                }

    # Absolute last resort - stop
    state = _log(state, "FALLBACK: Cannot find non-duplicate command, stopping")
    return {
        **state,
        "command": "STOP",
        "stop": True,
        "stop_reason": "all_commands_duplicate",
        "validation_error": None,
        "fallback_used": True
    }


# =============================================================================
# NODE: OUTPUT
# =============================================================================

def output_node(state):
    """
    Final node: Prepare response for client.

    Adds final debug message and ensures state is clean.
    """
    if state.get("stop"):
        state = _log(state, "OUTPUT: Workflow stopped - %s" % state.get("stop_reason", "unknown"))
    elif state.get("fallback_used"):
        state = _log(state, "OUTPUT: Used fallback command")
    else:
        state = _log(state, "OUTPUT: Normal completion")

    return state
