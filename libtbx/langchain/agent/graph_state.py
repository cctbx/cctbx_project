"""
Agent State Definition for PHENIX AI Agent.

Defines the TypedDict used as state throughout the LangGraph workflow.
"""

from __future__ import absolute_import, division, print_function
from typing import TypedDict, List, Dict, Optional


class AgentState(TypedDict):
    """
    State object passed through the LangGraph workflow.

    Categories:
    - SESSION CONTEXT: Persists across API calls (passed from client)
    - CURRENT CYCLE DATA: Data for the current cycle
    - CURRENT ATTEMPT: Resets on retry within a cycle
    - WORKFLOW TRACKING: Derived state for decision making
    - CONTROL FLAGS: Stop signals and fallback indicators
    - DEBUG: Logging and diagnostics
    """

    # === SESSION CONTEXT (Persists across API calls) ===
    history: List[Dict]           # The "Backpack" of completed cycles
    available_files: List[str]    # The Client's file list
    cycle_number: int             # Current Cycle (1, 2, ...)
    max_cycles: int               # Safety limit (e.g., 20)
    user_advice: str              # User instructions/preferences
    provider: str                 # LLM provider: "google", "openai", or "ollama"
    session_resolution: Optional[float]  # Resolution from session (single source of truth)

    # === CONFIGURATION ===
    maximum_automation: bool      # If True, use fully automated cryo-EM path
    use_rules_only: bool          # If True, use rules-based selection (no LLM)

    # === CURRENT CYCLE DATA ===
    log_text: str                 # Input log for this cycle
    log_analysis: Dict            # Output of perceive node (metrics, errors)

    # === CURRENT ATTEMPT (Resets each cycle) ===
    intent: Dict                  # Output of plan node (IntentJSON)
    command: str                  # Output of build node (CLI string)
    validation_error: Optional[str]  # Output of validate node
    attempt_number: int           # 0, 1, 2 (resets each cycle)
    previous_attempts: List[Dict] # Failed intents in THIS cycle only

    # === WORKFLOW TRACKING (Transient - computed each call) ===
    metrics_history: List[Dict]   # Derived from history - metrics per cycle
    metrics_trend: Dict           # Output of analyze_metrics_trend()
    workflow_state: Dict          # Output of detect_workflow_state()

    # === CONTROL FLAGS ===
    stop: bool                    # True if workflow complete
    stop_reason: Optional[str]    # "success", "stuck", "plateau", etc.
    fallback_used: bool           # True if fallback node was triggered

    # === DEBUG ===
    debug_log: List[str]          # Accumulated debug messages


def create_initial_state(
    available_files,
    log_text="",
    history=None,
    user_advice="",
    max_cycles=20,
    cycle_number=1,
    provider="google",
    maximum_automation=True,
    session_resolution=None,
    use_rules_only=False
):
    """
    Factory function to create a properly initialized AgentState.

    Args:
        available_files: List of files available on client
        log_text: Log text to analyze
        history: List of previous cycle records
        user_advice: User instructions/preferences
        max_cycles: Maximum cycles before stopping
        cycle_number: Current cycle number
        provider: LLM provider - "google", "openai", or "ollama"
        maximum_automation: If True, use fully automated cryo-EM path
            - True (default): predict_and_build runs full pipeline
            - False: Stepwise control with intermediate checkpoints
        session_resolution: Resolution value from session (single source of truth)
        use_rules_only: If True, use rules-based selection instead of LLM.
            This runs the agent without calling any external LLM API.

    Returns:
        AgentState: Properly initialized state dict
    """
    return {
        # Session context
        "history": history if history is not None else [],
        "available_files": list(available_files),
        "cycle_number": cycle_number,
        "max_cycles": max_cycles,
        "user_advice": user_advice,
        "provider": provider,
        "session_resolution": session_resolution,

        # Configuration
        "maximum_automation": maximum_automation,
        "use_rules_only": use_rules_only,

        # Current cycle data
        "log_text": log_text,
        "log_analysis": {},

        # Current attempt
        "intent": {},
        "command": "",
        "validation_error": None,
        "attempt_number": 0,
        "previous_attempts": [],

        # Workflow tracking (computed in nodes)
        "metrics_history": [],
        "metrics_trend": {},
        "workflow_state": {},

        # Control flags
        "stop": False,
        "stop_reason": None,
        "fallback_used": False,

        # Debug
        "debug_log": []
    }
