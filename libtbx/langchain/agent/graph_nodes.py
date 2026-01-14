"""
LangGraph Node Functions for PHENIX AI Agent.

Nodes:
- perceive: Extract metrics from log, build metrics history
- plan: Decide next action (LLM), with auto-stop on plateau
- build: Construct command from intent (templates)
- validate: Validate command (files, workflow state, duplicates)
- fallback: Mechanical auto-selection when LLM fails
- output_node: Prepare response for client
"""

from __future__ import absolute_import, division, print_function
import json
import re
import os

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
    """Helper to add debug message to state."""
    debug_log = list(state.get("debug_log", []))
    debug_log.append(msg)
    return {**state, "debug_log": debug_log}


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
    """
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

    # Detect experiment type for proper trend analysis
    available_files = state.get("available_files", [])
    has_map = any(f.lower().endswith(('.mrc', '.ccp4', '.map')) for f in available_files)
    has_mtz = any(f.lower().endswith('.mtz') for f in available_files)
    experiment_type = "cryoem" if has_map and not has_mtz else "xray"

    # Check for YAML mode (from state or global setting)
    use_yaml = state.get("use_yaml_mode", USE_YAML_METRICS)

    metrics_trend = analyze_metrics_trend(
        metrics_history, resolution, experiment_type,
        use_yaml_evaluator=use_yaml
    )

    # 5. Detect workflow state
    use_yaml_workflow = state.get("use_yaml_mode", USE_YAML_WORKFLOW)

    workflow_state = detect_workflow_state(
        history=history,
        available_files=available_files,
        analysis=analysis,
        maximum_automation=state.get("maximum_automation", True),
        use_yaml_engine=use_yaml_workflow
    )

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

    if metrics_trend.get("should_stop"):
        state = _log(state, "PERCEIVE: Auto-stop recommended: %s" % metrics_trend["reason"])

    return {
        **state,
        "log_analysis": analysis or {},
        "metrics_history": metrics_history,
        "metrics_trend": metrics_trend,
        "workflow_state": workflow_state,
    }


def _fallback_extract_metrics(log_text):
    """Fallback metric extraction if log_parsers module not available."""
    analysis = {}

    if not log_text:
        return analysis

    def safe_float(match):
        if match:
            value = match.group(1).rstrip('.')
            try:
                return float(value)
            except ValueError:
                return None
        return None

    # Basic metrics
    rfree_match = re.search(r'R-free\s*[:=]\s*([\d\.]+)', log_text, re.IGNORECASE)
    if rfree_match:
        val = safe_float(rfree_match)
        if val is not None:
            analysis["r_free"] = val

    rwork_match = re.search(r'R-work\s*[:=]\s*([\d\.]+)', log_text, re.IGNORECASE)
    if rwork_match:
        val = safe_float(rwork_match)
        if val is not None:
            analysis["r_work"] = val

    # TFZ/LLG from SOLU SET (Phaser)
    solu_set_match = re.search(r'SOLU SET\s+(.+)', log_text)
    if solu_set_match:
        solu_line = solu_set_match.group(1)
        tfz_matches = re.findall(r'TFZ==?([\d\.]+)', solu_line)
        if tfz_matches:
            analysis["tfz"] = float(tfz_matches[-1])
        llg_matches = re.findall(r'LLG=([\d\.]+)', solu_line)
        if llg_matches:
            analysis["llg"] = float(llg_matches[-1])

    # Map-model CC (for cryo-EM)
    cc_match = re.search(r'(?:CC_mask|map.model.CC|CC)\s*[:=]\s*([\d\.]+)', log_text, re.IGNORECASE)
    if cc_match:
        val = safe_float(cc_match)
        if val is not None:
            analysis["map_cc"] = val

    # Resolution
    res_match = re.search(r'Resolution\s*:\s*[\d\.]+\s*-\s*([\d\.]+)', log_text, re.IGNORECASE)
    if res_match:
        val = safe_float(res_match)
        if val is not None:
            analysis["resolution"] = val

    # Error detection
    if "FAILED:" in log_text:
        fail_match = re.search(r'FAILED:[:\s]+(.+?)(?:\n|$)', log_text)
        if fail_match:
            analysis["error"] = fail_match.group(1).strip()
    elif "Sorry:" in log_text:
        sorry_match = re.search(r'Sorry[:\s]+(.+?)(?:\n|$)', log_text)
        if sorry_match:
            analysis["error"] = sorry_match.group(1).strip()

    # Program detection
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
    1. Checks for auto-stop conditions (plateau, success)
    2. Optionally uses rules-only mode (no LLM)
    3. Constructs prompt with workflow state and metrics trend
    4. Calls LLM to decide next action
    5. Falls back to rules selector if LLM unavailable
    """
    # 1. Check for auto-stop from metrics trend
    metrics_trend = state.get("metrics_trend", {})

    if metrics_trend.get("should_stop"):
        reason = metrics_trend.get("reason", "Metrics indicate completion")
        state = _log(state, "PLAN: AUTO-STOP triggered: %s" % reason)

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

    system_msg, user_msg = get_planning_prompt(
        history=state.get("history", []),
        analysis=state.get("log_analysis", {}),
        available_files=state.get("available_files", []),
        previous_attempts=state.get("previous_attempts", []),
        user_advice=user_advice,
        metrics_trend=metrics_trend,
        workflow_state=workflow_state
    )

    # 4. Get LLM (with validation)
    llm, error = get_planning_llm(provider)

    if llm is None:
        state = _log(state, "PLAN: LLM error - %s" % error)
        state = _log(state, "PLAN: Falling back to mock planner")
        return _mock_plan(state)

    # 5. Call LLM
    try:
        state = _log(state, "PLAN: Calling LLM (provider=%s)..." % provider)

        from langchain_core.messages import SystemMessage, HumanMessage

        messages = [
            SystemMessage(content=system_msg),
            HumanMessage(content=user_msg)
        ]

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

        # === VALIDATE LLM CHOICE AGAINST VALID_PROGRAMS ===
        valid_programs = workflow_state.get("valid_programs", [])

        # Check if LLM is trying to STOP when STOP is not allowed
        if intent.get("stop") or chosen_program is None or chosen_program == "STOP":
            if "STOP" not in valid_programs:
                # LLM tried to stop but STOP is not allowed (e.g., validation required)
                state = _log(state, "PLAN: LLM tried to STOP but STOP not in valid_programs: %s" % valid_programs)
                state = _log(state, "PLAN: Overriding LLM choice - validation required before stopping")

                # Find the highest priority valid program (should be molprobity)
                override_program = valid_programs[0] if valid_programs else "phenix.refine"
                intent["program"] = override_program
                intent["stop"] = False
                intent["stop_reason"] = None
                intent["reasoning"] = "Override: %s (validation required before stopping)" % override_program
                state = _log(state, "PLAN: Overriding to %s" % override_program)

        # Also check if chosen program is in valid_programs
        elif chosen_program and chosen_program not in valid_programs and chosen_program != "STOP":
            state = _log(state, "PLAN: LLM chose invalid program %s, valid are: %s" % (chosen_program, valid_programs))
            # Use first valid program instead
            if valid_programs:
                override_program = next((p for p in valid_programs if p != "STOP"), valid_programs[0])
                intent["program"] = override_program
                state = _log(state, "PLAN: Overriding to %s" % override_program)

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
            prog = entry.get("program", "").lower()
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

def build(state):
    """
    Node C: Build Command from Intent using Templates.

    This is 100% deterministic - no LLM involved.
    The template builder guarantees correct syntax.

    Applies Tier 2 defaults from workflow_state and logs any overrides.
    """
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
    # Use cycle number as prefix to avoid PHASER.1_refine_001_refine_001_... filenames
    if program in ("phenix.refine", "phenix.real_space_refine"):
        if "output_prefix" not in strategy:
            history_info = state.get("history_info", {})
            cycle_number = history_info.get("cycle_number", 1)
            # Also count from history if cycle_number not set
            if cycle_number <= 1:
                history = state.get("history", [])
                refine_count = sum(1 for h in history
                                   if isinstance(h, dict) and
                                   h.get("program") in ("phenix.refine", "phenix.real_space_refine") and
                                   h.get("result") == "SUCCESS")
                cycle_number = refine_count + 1
            strategy["output_prefix"] = "%03d" % cycle_number
            state = _log(state, "BUILD: Using output_prefix=%03d to avoid long filenames" % cycle_number)

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

    llm_files = intent.get("files", {})
    corrected_files = {}
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
                            else:
                                corrected_list.append(fp)
                    if corrected_list:
                        corrected_files[slot] = corrected_list
                else:
                    # Single file
                    basename = os.path.basename(filepath)
                    if basename in basename_to_path:
                        # Use the correct path from available_files
                        corrected_files[slot] = basename_to_path[basename]
                    else:
                        # Keep original (might still be valid or will fail validation)
                        corrected_files[slot] = filepath

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
            command = builder.build_command_for_program(
                program,
                available_files,
                categorized_files=categorized_files
            )
            state = _log(state, "BUILD: Auto-selected files for %s" % program)
        else:
            # LLM picked some files - auto-fill any missing required files
            workflow_state = state.get("workflow_state", {})
            categorized_files = workflow_state.get("categorized_files", {})

            # Get required inputs for this program
            try:
                required_inputs = builder._registry.get_required_inputs(program)
            except Exception:
                required_inputs = []

            # Auto-fill missing required files from categorized_files
            # Use YAML input_priorities if defined, otherwise simple defaults
            for input_name in required_inputs:
                if input_name not in corrected_files:
                    file_found = None

                    # Try to get priorities from YAML
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
                    else:
                        # Simple fallback for common input types
                        if input_name == "mtz":
                            if categorized_files.get("refined_mtz"):
                                file_found = categorized_files["refined_mtz"][0]
                            elif categorized_files.get("mtz"):
                                file_found = categorized_files["mtz"][0]
                        elif input_name == "model":
                            for cat in ["refined", "phaser_output", "rsr_output",
                                       "processed_predicted", "pdb"]:
                                if categorized_files.get(cat):
                                    file_found = categorized_files[cat][0]
                                    break
                        elif input_name == "sequence" and categorized_files.get("sequence"):
                            file_found = categorized_files["sequence"][0]
                        elif input_name == "map" and categorized_files.get("full_map"):
                            file_found = categorized_files["full_map"][0]

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

            command = builder.build_command(
                program,
                corrected_files,
                strategy
            )
            state = _log(state, "BUILD: Used LLM-selected files (paths corrected)")

        if not command:
            error_msg = "Builder returned empty command for %s" % program
            state = _log(state, "BUILD: Failed - %s" % error_msg)
            return {**state, "command": "", "validation_error": error_msg}

    except Exception as e:
        state = _log(state, "BUILD: Exception - %s" % str(e))
        return {**state, "command": "", "validation_error": str(e)}

    state = _log(state, "BUILD: Command = %s" % command[:80])
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
    """
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

    # Try each valid program until one works without duplicating
    for program in runnable:
        # For refine, try with prioritized files (recent outputs first)
        if "refine" in program.lower():
            cmd = builder.build_command_for_program(program, prioritized_files,
                                                    context=fallback_context)
        else:
            cmd = builder.build_command_for_program(program, available_files,
                                                    context=fallback_context)

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
                                                context=fallback_context)
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
