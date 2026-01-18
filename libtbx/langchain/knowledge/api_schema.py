"""
API Schema Definitions for PHENIX AI Agent Client-Server Communication.

This module defines the versioned request/response schemas for communication
between the client (user's PHENIX installation) and the server (decision-making).

Design Principles:
1. Backwards Compatible: Old clients work with new servers
2. Graceful Degradation: Missing fields have sensible defaults
3. Forward Flexible: Unknown fields are ignored
4. Self-Documenting: Schema serves as API documentation

Usage:
    from knowledge.api_schema import (
        RequestV2, ResponseV2,
        validate_request, validate_response,
        apply_request_defaults, apply_response_defaults
    )
"""

from __future__ import absolute_import, division, print_function

# Current API version
API_VERSION = "2.0"

# Minimum supported API version (for backwards compatibility)
MIN_SUPPORTED_VERSION = "1.0"


# =============================================================================
# REQUEST SCHEMA (Client -> Server)
# =============================================================================

# Required fields for v2 request
REQUEST_V2_REQUIRED = [
    "api_version",
    "files",
    "cycle_number",
]

# All fields with their types and defaults
REQUEST_V2_SCHEMA = {
    # Protocol version
    "api_version": {
        "type": str,
        "default": "2.0",
        "description": "API protocol version",
    },
    "client_version": {
        "type": str,
        "default": None,
        "description": "Client software version (e.g., '1.22.0')",
    },

    # Core data
    "log_content": {
        "type": str,
        "default": "",
        "description": "Log text from the previous command",
    },
    "files": {
        "type": list,
        "default": [],
        "description": "List of available file paths (absolute)",
    },
    "history": {
        "type": list,
        "default": [],
        "description": "List of previous cycle records",
    },

    # Session state (new in v2 - enables better decisions)
    "session_state": {
        "type": dict,
        "default": {},
        "description": "Current session state from client",
        "subfields": {
            "resolution": {
                "type": float,
                "default": None,
                "description": "Session resolution in Angstroms",
            },
            "experiment_type": {
                "type": str,
                "default": None,
                "description": "Locked experiment type: 'xray' or 'cryoem'",
            },
            "rfree_mtz": {
                "type": str,
                "default": None,
                "description": "Path to locked R-free MTZ file",
            },
            "best_files": {
                "type": dict,
                "default": {},
                "description": "Current best files by category",
            },
        },
    },

    # User input
    "user_advice": {
        "type": str,
        "default": "",
        "description": "User instructions/guidelines",
    },

    # Settings
    "settings": {
        "type": dict,
        "default": {},
        "description": "Execution settings",
        "subfields": {
            "provider": {
                "type": str,
                "default": "google",
                "description": "LLM provider: 'google', 'openai', 'anthropic', 'ollama'",
            },
            "abort_on_red_flags": {
                "type": bool,
                "default": True,
                "description": "Abort on critical sanity check failures",
            },
            "abort_on_warnings": {
                "type": bool,
                "default": False,
                "description": "Also abort on warning-level issues",
            },
            "max_cycles": {
                "type": int,
                "default": 20,
                "description": "Maximum number of cycles",
            },
            "use_rules_only": {
                "type": bool,
                "default": False,
                "description": "Use rules-based selection without LLM",
            },
            "maximum_automation": {
                "type": bool,
                "default": True,
                "description": "If True, use fully automated path. If False, use stepwise mode with checkpoints.",
            },
        },
    },

    # Metadata
    "cycle_number": {
        "type": int,
        "default": 1,
        "description": "Current cycle number",
    },
}


# =============================================================================
# RESPONSE SCHEMA (Server -> Client)
# =============================================================================

# Required fields for v2 response
RESPONSE_V2_REQUIRED = [
    "api_version",
    "decision",
]

# All fields with their types and defaults
RESPONSE_V2_SCHEMA = {
    # Protocol version
    "api_version": {
        "type": str,
        "default": "2.0",
        "description": "API protocol version",
    },
    "server_version": {
        "type": str,
        "default": None,
        "description": "Server software version",
    },

    # Core decision
    "decision": {
        "type": dict,
        "default": {},
        "description": "The decision made by the server",
        "subfields": {
            "program": {
                "type": str,
                "default": "",
                "description": "Selected program (e.g., 'phenix.refine')",
            },
            "command": {
                "type": str,
                "default": "",
                "description": "Full command to execute",
            },
            "reasoning": {
                "type": str,
                "default": "",
                "description": "Explanation of why this decision was made",
            },
            "strategy": {
                "type": dict,
                "default": {},
                "description": "Strategy options used (e.g., resolution, output_prefix)",
            },
            "confidence": {
                "type": str,
                "default": "unknown",
                "description": "Confidence level: 'high', 'medium', 'low', 'unknown'",
            },
        },
    },

    # Stop signal
    "stop": {
        "type": bool,
        "default": False,
        "description": "Whether workflow should stop",
    },
    "stop_reason": {
        "type": str,
        "default": None,
        "description": "Reason for stopping: 'converged', 'red_flag', 'max_cycles', etc.",
    },

    # Metadata
    "metadata": {
        "type": dict,
        "default": {},
        "description": "Additional information from server",
        "subfields": {
            "experiment_type": {
                "type": str,
                "default": None,
                "description": "Detected experiment type",
            },
            "workflow_state": {
                "type": str,
                "default": None,
                "description": "Current workflow state",
            },
            "warnings": {
                "type": list,
                "default": [],
                "description": "Warning messages",
            },
            "red_flags": {
                "type": list,
                "default": [],
                "description": "Critical issues detected",
            },
            "final_metrics": {
                "type": dict,
                "default": {},
                "description": "Final quality metrics (when stopping)",
            },
        },
    },

    # Debug info
    "debug": {
        "type": dict,
        "default": {},
        "description": "Debug information",
        "subfields": {
            "log": {
                "type": list,
                "default": [],
                "description": "Debug log messages",
            },
            "timing_ms": {
                "type": int,
                "default": None,
                "description": "Processing time in milliseconds",
            },
        },
    },

    # Error handling
    "error": {
        "type": str,
        "default": None,
        "description": "Error message if request failed",
    },
}


# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

def validate_request(request, strict=False):
    """
    Validate a v2 request.

    Args:
        request: Dict containing request data
        strict: If True, fail on unknown fields; if False, ignore them

    Returns:
        tuple: (is_valid, errors_list)
    """
    errors = []

    if not isinstance(request, dict):
        return False, ["Request must be a dictionary"]

    # Check required fields
    for field in REQUEST_V2_REQUIRED:
        if field not in request:
            errors.append(f"Missing required field: {field}")

    # Check field types
    for field, spec in REQUEST_V2_SCHEMA.items():
        if field in request:
            value = request[field]
            expected_type = spec["type"]

            # Allow None for optional fields
            if value is None and field not in REQUEST_V2_REQUIRED:
                continue

            if not isinstance(value, expected_type):
                errors.append(
                    f"Field '{field}' has wrong type: expected {expected_type.__name__}, "
                    f"got {type(value).__name__}"
                )

    # Check for unknown fields (only in strict mode)
    if strict:
        for field in request:
            if field not in REQUEST_V2_SCHEMA:
                errors.append(f"Unknown field: {field}")

    return len(errors) == 0, errors


def validate_response(response, strict=False):
    """
    Validate a v2 response.

    Args:
        response: Dict containing response data
        strict: If True, fail on unknown fields; if False, ignore them

    Returns:
        tuple: (is_valid, errors_list)
    """
    errors = []

    if not isinstance(response, dict):
        return False, ["Response must be a dictionary"]

    # Check required fields
    for field in RESPONSE_V2_REQUIRED:
        if field not in response:
            errors.append(f"Missing required field: {field}")

    # Check field types
    for field, spec in RESPONSE_V2_SCHEMA.items():
        if field in response:
            value = response[field]
            expected_type = spec["type"]

            # Allow None for optional fields
            if value is None and field not in RESPONSE_V2_REQUIRED:
                continue

            if not isinstance(value, expected_type):
                errors.append(
                    f"Field '{field}' has wrong type: expected {expected_type.__name__}, "
                    f"got {type(value).__name__}"
                )

    # Check for unknown fields (only in strict mode)
    if strict:
        for field in response:
            if field not in RESPONSE_V2_SCHEMA:
                errors.append(f"Unknown field: {field}")

    return len(errors) == 0, errors


# =============================================================================
# DEFAULT APPLICATION
# =============================================================================

def apply_request_defaults(request):
    """
    Apply default values to a request, filling in missing optional fields.

    Args:
        request: Dict containing request data (modified in place)

    Returns:
        dict: The request with defaults applied
    """
    if not isinstance(request, dict):
        request = {}

    for field, spec in REQUEST_V2_SCHEMA.items():
        if field not in request:
            request[field] = spec["default"]
        elif spec["type"] == dict and "subfields" in spec:
            # Apply defaults to nested dict
            if request[field] is None:
                request[field] = {}
            # Only apply subfield defaults if they have non-None defaults
            # Don't add None values - that would overwrite intentionally missing fields
            for subfield, subspec in spec["subfields"].items():
                if subfield not in request[field] and subspec["default"] is not None:
                    request[field][subfield] = subspec["default"]

    return request


def apply_response_defaults(response):
    """
    Apply default values to a response, filling in missing optional fields.

    Args:
        response: Dict containing response data (modified in place)

    Returns:
        dict: The response with defaults applied
    """
    if not isinstance(response, dict):
        response = {}

    for field, spec in RESPONSE_V2_SCHEMA.items():
        if field not in response:
            response[field] = spec["default"]
        elif spec["type"] == dict and "subfields" in spec:
            # Apply defaults to nested dict
            if response[field] is None:
                response[field] = {}
            for subfield, subspec in spec["subfields"].items():
                if subfield not in response[field]:
                    response[field][subfield] = subspec["default"]

    return response


# =============================================================================
# VERSION HELPERS
# =============================================================================

def parse_version(version_str):
    """
    Parse a version string into a tuple.

    Args:
        version_str: Version string like "2.0" or "2.1.3"

    Returns:
        tuple: Version as tuple of ints, e.g., (2, 0) or (2, 1, 3)
    """
    if not version_str:
        return (1, 0)  # Default to v1

    try:
        parts = version_str.split(".")
        return tuple(int(p) for p in parts)
    except (ValueError, AttributeError):
        return (1, 0)


def is_version_supported(version_str):
    """
    Check if a version is supported.

    Args:
        version_str: Version string to check

    Returns:
        bool: True if version is supported
    """
    version = parse_version(version_str)
    min_version = parse_version(MIN_SUPPORTED_VERSION)
    current_version = parse_version(API_VERSION)

    return min_version <= version <= current_version


def get_api_version():
    """Get the current API version."""
    return API_VERSION


# =============================================================================
# CONVENIENCE CONSTRUCTORS
# =============================================================================

def create_request(
    files,
    cycle_number,
    log_content="",
    history=None,
    session_state=None,
    user_advice="",
    settings=None,
    client_version=None,
):
    """
    Create a v2 request with proper defaults.

    Args:
        files: List of available file paths
        cycle_number: Current cycle number
        log_content: Log text from previous command
        history: List of previous cycle records
        session_state: Session state dict
        user_advice: User instructions
        settings: Settings dict
        client_version: Client version string

    Returns:
        dict: A valid v2 request
    """
    request = {
        "api_version": API_VERSION,
        "client_version": client_version,
        "log_content": log_content or "",
        "files": list(files) if files else [],
        "history": list(history) if history else [],
        "session_state": dict(session_state) if session_state else {},
        "user_advice": user_advice or "",
        "settings": dict(settings) if settings else {},
        "cycle_number": cycle_number,
    }

    return apply_request_defaults(request)


def create_response(
    program,
    command,
    reasoning="",
    strategy=None,
    stop=False,
    stop_reason=None,
    experiment_type=None,
    workflow_state=None,
    warnings=None,
    red_flags=None,
    debug_log=None,
    error=None,
    server_version=None,
):
    """
    Create a v2 response with proper defaults.

    Args:
        program: Selected program name
        command: Full command to execute
        reasoning: Explanation of decision
        strategy: Strategy options dict
        stop: Whether to stop workflow
        stop_reason: Reason for stopping
        experiment_type: Detected experiment type
        workflow_state: Current workflow state
        warnings: List of warning messages
        red_flags: List of critical issues
        debug_log: List of debug messages
        error: Error message if failed
        server_version: Server version string

    Returns:
        dict: A valid v2 response
    """
    response = {
        "api_version": API_VERSION,
        "server_version": server_version,
        "decision": {
            "program": program or "",
            "command": command or "",
            "reasoning": reasoning or "",
            "strategy": dict(strategy) if strategy else {},
            "confidence": "unknown",
        },
        "stop": bool(stop),
        "stop_reason": stop_reason,
        "metadata": {
            "experiment_type": experiment_type,
            "workflow_state": workflow_state,
            "warnings": list(warnings) if warnings else [],
            "red_flags": list(red_flags) if red_flags else [],
            "final_metrics": {},
        },
        "debug": {
            "log": list(debug_log) if debug_log else [],
            "timing_ms": None,
        },
        "error": error,
    }

    return apply_response_defaults(response)


def create_error_response(error_message, debug_log=None):
    """
    Create a v2 error response.

    Args:
        error_message: Error description
        debug_log: Optional debug messages

    Returns:
        dict: A v2 response indicating error
    """
    return create_response(
        program="",
        command="",
        reasoning="",
        stop=True,
        stop_reason="error",
        error=error_message,
        debug_log=debug_log,
    )


def create_stop_response(
    stop_reason,
    reasoning="",
    final_metrics=None,
    debug_log=None,
    server_version=None,
):
    """
    Create a v2 stop response.

    Args:
        stop_reason: Why workflow is stopping
        reasoning: Explanation
        final_metrics: Final quality metrics
        debug_log: Debug messages
        server_version: Server version

    Returns:
        dict: A v2 response indicating stop
    """
    response = create_response(
        program="STOP",
        command="STOP",
        reasoning=reasoning,
        stop=True,
        stop_reason=stop_reason,
        debug_log=debug_log,
        server_version=server_version,
    )

    if final_metrics:
        response["metadata"]["final_metrics"] = dict(final_metrics)

    return response
