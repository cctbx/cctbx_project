"""
API Client Adapter for PHENIX AI Agent.

This module provides functions to:
1. Build v2 requests from client data
2. Parse v2 responses into the format expected by existing code
3. Serialize/deserialize for transport

Usage:
    from libtbx.langchain.agent.api_client import (
        build_request_v2,
        parse_response_v2,
        serialize_request,
        deserialize_response,
    )
"""

from __future__ import absolute_import, division, print_function

import json
import os

# Import schema definitions
from libtbx.langchain.knowledge.api_schema import (
    create_request,
    apply_request_defaults,
    apply_response_defaults,
)

# Import transport functions for sanitization
from libtbx.langchain.agent.transport import (
    sanitize_string,
    sanitize_for_transport,
    sanitize_dict_recursive,
    sanitize_request,
    sanitize_response,
    get_transport_config,
)
TRANSPORT_AVAILABLE = True

# Silence unused import warnings - these are used conditionally
assert sanitize_request is not None
assert sanitize_response is not None
assert get_transport_config is not None


# =============================================================================
# SESSION STATE BUILDING
# =============================================================================

def build_session_state(session_info, session_resolution=None):
    """
    Build session_state dict from session_info.

    Used by both LocalAgent and RemoteAgent to avoid code duplication.
    The session_state is passed to build_request_v2 and transmitted to the server.

    Args:
        session_info: Dict with session state fields:
            - experiment_type: Locked experiment type
            - rfree_mtz: Locked R-free MTZ path
            - best_files: Current best files dict
            - directives: Extracted user directives
            - force_retry_program: Program to force retry (error recovery)
            - recovery_strategies: Dict of recovery strategies per file
        session_resolution: Optional resolution value

    Returns:
        dict: Session state ready for build_request_v2
    """
    session_state = {}

    if session_resolution is not None:
        session_state["resolution"] = session_resolution

    if session_info:
        # Standard session fields
        for key in ["experiment_type", "rfree_mtz", "best_files", "directives"]:
            if session_info.get(key):
                session_state[key] = session_info[key]

        # Error recovery fields
        for key in ["force_retry_program", "recovery_strategies"]:
            if session_info.get(key):
                session_state[key] = session_info[key]

        # Explicit program request (from user advice)
        if session_info.get("explicit_program"):
            session_state["explicit_program"] = session_info["explicit_program"]

    return session_state


# =============================================================================
# REQUEST BUILDING (Client -> Server)
# =============================================================================

def build_request_v2(
    files,
    cycle_number,
    log_content="",
    history=None,
    session_state=None,
    user_advice="",
    provider="google",
    abort_on_red_flags=True,
    abort_on_warnings=False,
    client_version=None,
):
    """
    Build a v2 API request from client data.

    This is the main entry point for clients to create requests.

    Args:
        files: List of available file paths
        cycle_number: Current cycle number
        log_content: Log text from previous command
        history: List of previous cycle records (dicts)
        session_state: Dict with session state:
            - resolution: float
            - experiment_type: str ("xray" or "cryoem")
            - rfree_mtz: str (path to locked MTZ)
            - best_files: dict (category -> path)
            - directives: dict (extracted user directives)
        user_advice: User instructions/guidelines
        provider: LLM provider name
        abort_on_red_flags: Whether to abort on critical issues
        abort_on_warnings: Whether to abort on warnings
        client_version: Client software version

    Returns:
        dict: A valid v2 request ready for serialization
    """
    # Normalize files to absolute paths
    normalized_files = []
    for f in (files or []):
        if f and isinstance(f, str):
            normalized_files.append(os.path.abspath(f) if not os.path.isabs(f) else f)

    # Normalize history using transport sanitization
    normalized_history = []

    for h in (history or []):
        if isinstance(h, dict):
            # Sanitize metrics recursively
            # NOTE: get_history_for_agent returns 'analysis' key, but we store as 'metrics'
            raw_metrics = h.get("analysis", h.get("metrics", {}))
            sanitized_metrics = sanitize_dict_recursive(raw_metrics, max_len_per_string=500)

            normalized_history.append({
                "cycle": h.get("cycle_number", h.get("cycle", len(normalized_history) + 1)),
                "program": sanitize_string(h.get("program", ""), max_len=100),
                "command": sanitize_string(h.get("command", ""), max_len=1000),
                "result": sanitize_for_transport(
                    h.get("result", ""),
                    max_len=500,
                    truncate_quotes=True,
                    quote_max_len=200
                ),
                "output_files": h.get("output_files", []),
                "metrics": sanitized_metrics,
            })

    # Build session_state
    normalized_session_state = {}
    if session_state:
        if session_state.get("resolution") is not None:
            normalized_session_state["resolution"] = float(session_state["resolution"])
        if session_state.get("experiment_type"):
            normalized_session_state["experiment_type"] = session_state["experiment_type"]
        if session_state.get("rfree_mtz"):
            normalized_session_state["rfree_mtz"] = session_state["rfree_mtz"]
        if session_state.get("best_files"):
            normalized_session_state["best_files"] = dict(session_state["best_files"])
        if session_state.get("directives"):
            normalized_session_state["directives"] = dict(session_state["directives"])
        # Error recovery state
        if session_state.get("force_retry_program"):
            normalized_session_state["force_retry_program"] = session_state["force_retry_program"]
        if session_state.get("recovery_strategies"):
            normalized_session_state["recovery_strategies"] = dict(session_state["recovery_strategies"])
        # Explicit program request
        if session_state.get("explicit_program"):
            normalized_session_state["explicit_program"] = session_state["explicit_program"]

    # Build settings
    settings = {
        "provider": provider or "google",
        "abort_on_red_flags": bool(abort_on_red_flags),
        "abort_on_warnings": bool(abort_on_warnings),
    }

    # Sanitize log_content with quoted string truncation
    # This handles large data dumps like pdb70_text='...'
    sanitized_log_content = sanitize_for_transport(
        log_content or "",
        max_len=50000,
        truncate_quotes=True,
        quote_max_len=500
    )

    # Also sanitize user_advice (user could input tabs)
    sanitized_user_advice = sanitize_string(user_advice or "", max_len=10000)

    # Create the request using schema helper
    return create_request(
        files=normalized_files,
        cycle_number=cycle_number,
        log_content=sanitized_log_content,
        history=normalized_history,
        session_state=normalized_session_state,
        user_advice=sanitized_user_advice,
        settings=settings,
        client_version=client_version,
    )


def build_request_from_params(params, session=None, log_content="", history=None):
    """
    Build a v2 request from PHENIX params object and session.

    This provides compatibility with existing code that uses params.

    Args:
        params: PHENIX params object
        session: AgentSession object (optional)
        log_content: Log text from previous command
        history: List of history records

    Returns:
        dict: A valid v2 request
    """
    # Extract files from params
    files = []
    if hasattr(params, 'ai_analysis') and params.ai_analysis.original_files:
        files = list(params.ai_analysis.original_files)

    # Extract cycle number
    cycle_number = 1
    if history:
        cycle_number = len(history) + 1

    # Extract user advice
    user_advice = ""
    if hasattr(params, 'ai_analysis') and params.ai_analysis.project_advice:
        user_advice = params.ai_analysis.project_advice

    # Extract provider
    provider = "google"
    if hasattr(params, 'communication') and params.communication.provider:
        provider = params.communication.provider

    # Build session state from session object
    session_state = {}
    if session:
        if hasattr(session, 'get_resolution') and session.get_resolution():
            session_state["resolution"] = session.get_resolution()
        if hasattr(session, 'get_experiment_type') and session.get_experiment_type():
            session_state["experiment_type"] = session.get_experiment_type()
        if hasattr(session, 'get_rfree_mtz') and session.get_rfree_mtz():
            session_state["rfree_mtz"] = session.get_rfree_mtz()
        if hasattr(session, 'get_best_files_dict'):
            session_state["best_files"] = session.get_best_files_dict()
        if hasattr(session, 'get_directives') and session.get_directives():
            session_state["directives"] = session.get_directives()

    # Extract abort settings
    abort_on_red_flags = True
    abort_on_warnings = False
    if session and hasattr(session, 'get_abort_settings'):
        settings = session.get_abort_settings()
        abort_on_red_flags = settings.get("abort_on_red_flags", True)
        abort_on_warnings = settings.get("abort_on_warnings", False)

    return build_request_v2(
        files=files,
        cycle_number=cycle_number,
        log_content=log_content,
        history=history,
        session_state=session_state,
        user_advice=user_advice,
        provider=provider,
        abort_on_red_flags=abort_on_red_flags,
        abort_on_warnings=abort_on_warnings,
    )


# =============================================================================
# RESPONSE PARSING (Server -> Client)
# =============================================================================

def parse_response_v2(response):
    """
    Parse a v2 API response into the format expected by existing client code.

    Args:
        response: Dict containing v2 response data

    Returns:
        dict: Response in legacy format compatible with existing code:
            {
                "next_move": {
                    "program": str,
                    "command": str,
                    "explanation": str,
                    "strategy": dict,
                },
                "program": str,
                "command": str,
                "history_record": dict,
                "debug_log": list,
                "error": str or None,
            }
    """
    if not isinstance(response, dict):
        return {
            "next_move": None,
            "program": "",
            "command": "",
            "history_record": None,
            "debug_log": ["ERROR: Invalid response format"],
            "error": "Invalid response format",
        }

    # Apply defaults to ensure all fields exist
    response = apply_response_defaults(response)

    decision = response.get("decision", {})
    metadata = response.get("metadata", {})
    debug = response.get("debug", {})

    program = decision.get("program", "")
    command = decision.get("command", "")
    reasoning = decision.get("reasoning", "")
    strategy = decision.get("strategy", {})

    # Build next_move in legacy format
    next_move = {
        "program": program,
        "command": command,
        "explanation": reasoning,
        "strategy": strategy,
        "confidence": decision.get("confidence", "unknown"),
    }

    # Build history_record in legacy format
    history_record = {
        "program": program,
        "command": command,
        "reasoning": reasoning,
        "summary": f"Agent selected {program}" if program else "",
        "analysis": "",
        "error": response.get("error"),
        "intent": {
            "program": program,
            "reasoning": reasoning,
            "strategy": strategy,
        },
        "next_move": next_move,
        "debug_log": debug.get("log", []),
        "events": response.get("events", []),  # Include events
        "experiment_type": metadata.get("experiment_type"),
        "workflow_state": metadata.get("workflow_state"),
        "stop": response.get("stop", False),
        "stop_reason": response.get("stop_reason"),
        "red_flag_issues": metadata.get("red_flags", []),
        "warnings": metadata.get("warnings", []),
    }

    # Handle stop case
    if response.get("stop"):
        history_record["stop"] = True
        if response.get("stop_reason") == "red_flag":
            history_record["abort_message"] = reasoning

    return {
        "next_move": next_move,
        "program": program,
        "command": command,
        "history_record": history_record,
        "debug_log": debug.get("log", []),
        "events": response.get("events", []),  # Include events
        "error": response.get("error"),
    }


def parse_response_to_group_args(response):
    """
    Parse a v2 response into a group_args object for compatibility.

    This matches the format returned by run_job_on_server.

    Args:
        response: Dict containing v2 response data

    Returns:
        group_args: Object with expected attributes
    """
    from libtbx import group_args

    parsed = parse_response_v2(response)

    return group_args(
        group_args_type='log summary',
        next_move=parsed["next_move"],
        history_record=parsed["history_record"],
        debug_log=parsed["debug_log"],
        events=parsed.get("events", []),  # Include events
        summary=parsed["history_record"].get("summary", ""),
        analysis=parsed["history_record"].get("analysis", ""),
        error=parsed["error"],
        program=parsed["program"],
    )


# =============================================================================
# SERIALIZATION
# =============================================================================

def serialize_request(request):
    """
    Serialize a v2 request to JSON string for transport.

    Args:
        request: Dict containing v2 request data

    Returns:
        str: JSON-encoded request
    """
    json_str = json.dumps(request, indent=None, separators=(',', ':'))
    # Final safety: replace any remaining JSON tab escape sequences with spaces
    # The transport module should have removed tabs, but this catches edge cases
    json_str = json_str.replace('\\t', ' ')
    return json_str


def deserialize_request(json_str):
    """
    Deserialize a JSON string to v2 request.

    Args:
        json_str: JSON-encoded request string

    Returns:
        dict: Parsed request with defaults applied
    """
    try:
        request = json.loads(json_str)
        return apply_request_defaults(request)
    except (json.JSONDecodeError, TypeError) as e:
        return apply_request_defaults({"error": f"Failed to parse request: {e}"})


def serialize_response(response):
    """
    Serialize a v2 response to JSON string for transport.

    Args:
        response: Dict containing v2 response data

    Returns:
        str: JSON-encoded response
    """
    return json.dumps(response, indent=None, separators=(',', ':'))


def deserialize_response(json_str):
    """
    Deserialize a JSON string to v2 response.

    Args:
        json_str: JSON-encoded response string

    Returns:
        dict: Parsed response with defaults applied
    """
    try:
        response = json.loads(json_str)
        return apply_response_defaults(response)
    except (json.JSONDecodeError, TypeError) as e:
        return apply_response_defaults({"error": f"Failed to parse response: {e}"})


# =============================================================================
# VERSION DETECTION
# =============================================================================

def detect_api_version(data):
    """
    Detect the API version from request or response data.

    Args:
        data: Dict that may contain api_version field

    Returns:
        str: Detected version, or "1.0" if not found
    """
    if isinstance(data, dict):
        return data.get("api_version", "1.0")
    return "1.0"


def is_v2_request(data):
    """
    Check if data is a v2 format request.

    Args:
        data: Request data (dict or string)

    Returns:
        bool: True if this is a v2 request
    """
    if isinstance(data, str):
        try:
            data = json.loads(data)
        except (json.JSONDecodeError, TypeError):
            return False

    if not isinstance(data, dict):
        return False

    version = data.get("api_version", "")
    return version.startswith("2.")


def is_v2_response(data):
    """
    Check if data is a v2 format response.

    Args:
        data: Response data (dict or string)

    Returns:
        bool: True if this is a v2 response
    """
    return is_v2_request(data)  # Same logic


# =============================================================================
# SIMPLE LLM CALLS
# =============================================================================

def call_llm_simple(prompt, provider="google", model=None, temperature=0.1, max_tokens=2000):
    """
    Make a simple LLM call without the full agent infrastructure.

    This is used for one-off tasks like directive extraction where we don't
    need the full graph/session machinery.

    Args:
        prompt: The prompt text to send
        provider: LLM provider ("google", "openai", "anthropic", "ollama")
        model: Model name (uses provider default if None)
        temperature: Sampling temperature (default 0.1 for consistency)
        max_tokens: Maximum tokens in response

    Returns:
        str: The LLM response text, or None on failure
    """
    if provider == "google":
        return _call_google_llm(prompt, model, temperature, max_tokens)
    elif provider == "openai":
        return _call_openai_llm(prompt, model, temperature, max_tokens)
    elif provider == "anthropic":
        return _call_anthropic_llm(prompt, model, temperature, max_tokens)
    elif provider == "ollama":
        return _call_ollama_llm(prompt, model, temperature, max_tokens)
    else:
        raise ValueError(f"Unknown provider: {provider}")


def _call_google_llm(prompt, model=None, temperature=0.1, max_tokens=2000):
    """Call Google's Gemini API."""
    api_key = os.environ.get("GOOGLE_API_KEY")
    if not api_key:
        raise ValueError("GOOGLE_API_KEY environment variable not set")

    # Try new google.genai package first
    try:
        from google import genai
        from google.genai import types

        model_name = model or "gemini-2.0-flash"
        client = genai.Client(api_key=api_key)

        response = client.models.generate_content(
            model=model_name,
            contents=prompt,
            config=types.GenerateContentConfig(
                temperature=temperature,
                max_output_tokens=max_tokens
            )
        )
        return response.text

    except ImportError:
        # Fall back to deprecated google.generativeai
        import google.generativeai as genai_old

        genai_old.configure(api_key=api_key)
        model_name = model or "gemini-2.0-flash"
        gen_model = genai_old.GenerativeModel(model_name)

        response = gen_model.generate_content(
            prompt,
            generation_config=genai_old.types.GenerationConfig(
                temperature=temperature,
                max_output_tokens=max_tokens
            )
        )
        return response.text


def _call_openai_llm(prompt, model=None, temperature=0.1, max_tokens=2000):
    """Call OpenAI's API."""
    import openai

    model_name = model or "gpt-4o-mini"
    client = openai.OpenAI()  # Uses OPENAI_API_KEY env var

    response = client.chat.completions.create(
        model=model_name,
        messages=[{"role": "user", "content": prompt}],
        temperature=temperature,
        max_tokens=max_tokens
    )
    return response.choices[0].message.content


def _call_anthropic_llm(prompt, model=None, temperature=0.1, max_tokens=2000):
    """Call Anthropic's API."""
    import anthropic

    model_name = model or "claude-sonnet-4-20250514"
    client = anthropic.Anthropic()  # Uses ANTHROPIC_API_KEY env var

    response = client.messages.create(
        model=model_name,
        max_tokens=max_tokens,
        messages=[{"role": "user", "content": prompt}]
    )
    return response.content[0].text


def _call_ollama_llm(prompt, model=None, temperature=0.1, max_tokens=2000):
    """Call Ollama API (OpenAI-compatible)."""
    import openai

    model_name = model or "llama3.2"
    base_url = os.environ.get("OLLAMA_BASE_URL", "http://localhost:11434/v1")

    client = openai.OpenAI(
        base_url=base_url,
        api_key="ollama"  # Ollama doesn't require a real API key
    )

    response = client.chat.completions.create(
        model=model_name,
        messages=[{"role": "user", "content": prompt}],
        temperature=temperature,
        max_tokens=max_tokens
    )
    return response.choices[0].message.content
