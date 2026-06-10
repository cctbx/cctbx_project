"""
Agent State Definition for PHENIX AI Agent.

Defines the TypedDict used as state throughout the LangGraph workflow.
"""

from __future__ import absolute_import, division, print_function
from typing import TypedDict, List, Dict, Optional

# ---------------------------------------------------------------------------
# P1B: Typed stop-reason codes for THINK-node classified stops.
# think_stop_override["code"] must be one of these values.
# Unrecognised codes are logged and ignored (graceful degradation).
# ---------------------------------------------------------------------------
STOP_REASON_CODES = frozenset([
    "WRONG_MTZ",              # MTZ file has wrong content for this step
    "WRONG_SPACE_GROUP",      # Space group mismatch prevents solution
    "MISMATCHED_SEQUENCE",    # Sequence file doesn't match electron density
    "NO_SOLUTION_FOUND",      # Search exhausted, no plausible solution
    "REASONLESS_DIVERGENCE",  # R-free diverging with no recoverable cause
])


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
    session_info: Dict            # Session state including:
                                  #   - experiment_type: "xray" or "cryoem"
                                  #   - best_files: {category: path} from BestFilesTracker
    directives: Dict              # User directives extracted from advice
    bad_inject_params: Dict       # {program: [param_names]} — blacklisted params per program

    # === CONFIGURATION ===
    maximum_automation: bool      # If True, use fully automated cryo-EM path
    use_rules_only: bool          # If True, use rules-based selection (no LLM)
    abort_on_red_flags: bool      # If True, abort on critical sanity check failures
    abort_on_warnings: bool       # If True, also abort on warning-level issues
    thinking_level: Optional[str] # None, "basic", or "advanced"
                                  # ("expert" maps to "advanced" here;
                                  # planning gated in ai_agent.py)

    # === THINKING AGENT (v113) ===
    expert_assessment: Dict       # Output of think node (expert analysis)
    strategy_memory: Dict         # Accumulated scientific understanding across cycles

    # === GOAL-DIRECTED AGENT (v114) ===
    structure_model: Dict         # StructureModel state
    validation_history: Dict      # Per-cycle snapshots

    # === SESSION-LEVEL BLOCKING (P4 fix) ===
    session_blocked_programs: List[str]  # Programs blocked for this session (exceeded failure threshold)

    # === P1B FUTURE KEYS (think-node MTZ overrides) ===
    think_file_overrides: Dict    # {category: path} overrides from THINK node (P1B)
    think_label_hints: Dict       # {path: label_str} MTZ label hints from THINK node (P1B)
    think_stop_override: Optional[Dict]  # Typed stop override from THINK node (P1B)

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
    error_classification: Dict    # Error classification from PERCEIVE (v115)
    failure_count: int            # Consecutive failures of same program (v115)

    # === CONTROL FLAGS ===
    stop: bool                    # True if workflow complete
    stop_reason: Optional[str]    # "success", "stuck", "plateau", "red_flag", etc.
    fallback_used: bool           # True if fallback node was triggered

    # === RED FLAG DETECTION ===
    red_flag_issues: Optional[List[Dict]]  # Issues from sanity checker
    abort_message: Optional[str]           # Formatted abort message

    # === EVENTS / WARNINGS ===
    events: List[Dict]            # Structured events emitted by nodes (for client display)
    warnings: List[str]           # Server-side warnings passed back to client

    # === TRANSIENT WITHIN-CYCLE KEYS ===
    # These are written and read within a single cycle; they don't need round-tripping.
    _prev_valid_programs: List[str]  # Loop detection: valid_programs from previous PERCEIVE
    prerequisite_for: Optional[str]  # When set, this cycle is a prerequisite run for the named program

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
    use_rules_only=False,
    session_info=None,
    abort_on_red_flags=True,
    abort_on_warnings=False,
    directives=None,
    bad_inject_params=None,
    thinking_level=None,
    strategy_memory=None,
    structure_model=None,
    validation_history=None,
    session_blocked_programs=None
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
        session_info: Dict with session state:
            - experiment_type: Locked experiment type ("xray" or "cryoem")
            - best_files: Dict of {category: path} from BestFilesTracker
            - rfree_mtz: Path to locked R-free MTZ file (X-ray only)
        abort_on_red_flags: If True, abort on critical sanity check failures
        abort_on_warnings: If True, also abort on warning-level issues
        directives: Dict of user directives extracted from advice
        thinking_level: Controls expert reasoning depth.
            None: No thinking (default, pass-through).
            "basic": LLM reasoning with log analysis and strategy
            memory. No structural validation or knowledge base.
            "advanced": Full pipeline with structural validation,
            file metadata tracking, and expert knowledge base.
        strategy_memory: Dict of accumulated scientific understanding from
            previous cycles. Passed through session_info for persistence.
        structure_model: Dict from StructureModel.to_dict(). Running
            structural knowledge accumulated across cycles. Restored
            from session on resume. None = create fresh.
        validation_history: Dict from ValidationHistory.to_dict().
            Per-cycle validation snapshots for gate evaluator and
            explanation engine. Restored from session on resume.
            None = create fresh.
        session_blocked_programs: List of program names that have exceeded
            the per-program total-failure threshold this session (P4 fix).
            Persisted client-side and re-injected so get_workflow_state
            can filter them from valid_programs each cycle.
            None / [] = no programs currently blocked.

    Returns:
        AgentState: Properly initialized state dict
    """
    # Validate thinking_level
    # "expert" activates the planning layer in ai_agent.py
    # but maps to "advanced" for the graph (THINK node).
    if thinking_level is not None:
      thinking_level = str(thinking_level).lower()
      if thinking_level == "expert":
        thinking_level = "advanced"
      if thinking_level not in ("basic", "advanced"):
        thinking_level = None
    return {
        # Session context
        "history": list(history) if history is not None else [],
        "available_files": list(available_files),
        "cycle_number": cycle_number,
        "max_cycles": max_cycles,
        "user_advice": user_advice,
        "provider": provider,
        "session_resolution": session_resolution,
        "session_info": session_info if session_info is not None else {},
        "directives": directives if directives is not None else {},
        "bad_inject_params": bad_inject_params if bad_inject_params is not None else {},

        # Configuration
        "maximum_automation": maximum_automation,
        "use_rules_only": use_rules_only,
        "abort_on_red_flags": abort_on_red_flags,
        "abort_on_warnings": abort_on_warnings,
        "thinking_level": thinking_level,

        # Thinking agent (v113)
        "expert_assessment": {},
        "strategy_memory": strategy_memory if strategy_memory is not None else {},

        # Goal-directed agent (v114)
        "structure_model": (
          structure_model
          if structure_model is not None else {}
        ),
        "validation_history": (
          validation_history
          if validation_history is not None
          else {}
        ),

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
        "error_classification": {},
        "failure_count": 0,

        # Control flags
        "stop": False,
        "stop_reason": None,
        "fallback_used": False,

        # Red flag detection
        "red_flag_issues": None,
        "abort_message": None,

        # P4: session-level program block list.
        # Programs that have failed N times total in this session are added
        # here and excluded from valid_programs for the rest of the session.
        # Escape hatch: cleared for a program when user advice changes and
        # explicitly names the blocked program.
        # Per-program thresholds: refine/rsr = 6, all others = 4.
        # Restored from session on resume (like strategy_memory etc.).
        "session_blocked_programs": list(
            session_blocked_programs) if session_blocked_programs else [],

        # P1B: typed file overrides from THINK (one-cycle lifetime).
        # Key = input_name (e.g. 'data_mtz'), value = absolute path.
        # Consumed and cleared after BUILD uses it.
        "think_file_overrides": {},

        # P1B: label hints parallel to think_file_overrides.
        "think_label_hints": {},

        # think_stop_override: dict with 'code' and 'analysis' keys,
        # or None.  THINK sets this to force a stop with a classified reason.
        # Code must be in STOP_REASON_CODES (graph_state.py).
        "think_stop_override": None,

        # Debug
        "debug_log": []
    }
