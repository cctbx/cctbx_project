"""
Backward Compatibility Contract
================================

Single source of truth for the client <-> server interface.

Every field the server reads from session_info, and every field the server
returns to the client, is registered here.  Static analysis tests
(tst_contract_compliance.py) verify that the running code conforms.

See docs/guides/BACKWARD_COMPATIBILITY.md for the full plan.


HOW TO ADD A NEW SESSION_INFO FIELD
-------------------------------------
1. Add the entry to SESSION_INFO_FIELDS below (with the next version number).
2. In ai_agent.py (client), add it to the session_info dict in _query_agent().
3. On the server, access it ONLY via session_info.get("field", default).
   The default MUST match what you register here.
4. Run tst_contract_compliance.py to confirm.


HOW TO ADD A NEW RESPONSE FIELD
---------------------------------
1. Add the entry to RESPONSE_FIELDS below.
2. Ensure it is set in the graph state before output_node().
3. In ai_agent.py (client), access it via .get() or getattr() with a fallback.
4. Old clients that don't know about the field will simply ignore it.
5. NEVER remove an existing field -- old clients depend on them.
"""

# =========================================================================
# PROTOCOL VERSION
# =========================================================================

# Bump CURRENT when the client adds new session_info fields.
# Bump MIN_SUPPORTED only after a deprecation period (see BACKWARD_COMPATIBILITY.md).
CURRENT_PROTOCOL_VERSION = 3
MIN_SUPPORTED_PROTOCOL_VERSION = 1


# =========================================================================
# CLIENT -> SERVER: session_info fields
# =========================================================================
#
# Each entry: (field_name, default_value, added_version, description)
#
# Rules:
#   - The server MUST use  session_info.get(field, default)  for every access.
#   - The default here must produce correct (if degraded) behavior when the
#     field is absent (i.e. when an old client doesn't send it).
#   - Fields must never be removed.  Mark as deprecated if no longer used.
#
# NOTE on mutable defaults: The _normalize function below calls copy() on
# dict/list defaults, so registering {} or [] here is safe.

SESSION_INFO_FIELDS = [
    # --- v1: original fields ------------------------------------------------
    ("experiment_type",        "",    1,
     "Experiment type string: 'xray' or 'cryoem'"),

    ("best_files",             {},    1,
     "Dict of category -> best file path, from BestFilesTracker"),

    ("rfree_mtz",              None,  1,
     "Locked R-free MTZ file path (str or None)"),

    ("directives",             {},    1,
     "Parsed user directives dict (stop conditions, preferences, etc.)"),

    # --- v2: error recovery & user interaction ------------------------------
    ("rfree_resolution",       None,  2,
     "Resolution limit of R-free flags in Angstroms (float or None)"),

    ("force_retry_program",    None,  2,
     "One-shot forced program name for error recovery (str or None)"),

    ("recovery_strategies",    {},    2,
     "Dict of filename -> list of recovery strategy dicts"),

    ("explicit_program",       None,  2,
     "User-requested program extracted from project advice (str or None)"),

    ("advice_changed",         False, 2,
     "True when user provided new advice on resume; suppresses AUTO-STOP for one cycle"),

    ("bad_inject_params",      {},    2,
     "Dict of program_name -> set of parameter names to blacklist from injection"),

    # --- v3: client-side pre-extraction for server parity -------------------
    ("unplaced_model_cell",    None,  3,
     "Pre-extracted CRYST1 unit cell as [a, b, c, alpha, beta, gamma] or None"),

    ("model_hetatm_residues",  None,  3,
     "Pre-extracted HETATM residues as [[chain, resseq, resname], ...] or None"),

    ("client_protocol_version", 1,    3,
     "Protocol version of the sending client (int). "
     "Defaults to 1 for clients that predate this field."),

    # --- v4: goal-directed plan status ------------------------------------
    ("plan_has_pending_stages", False, 4,
     "True when the plan has stages still pending "
     "with programs that could run. Suppresses "
     "AUTO-STOP so the plan can drive the workflow."),

    # --- v5: ASU copy count (v115.04) ------------------------------------
    ("asu_copies",              None,  5,
     "Number of copies in the ASU (int or None). "
     "From user directives or xtriage log analysis."),
]

# Convenience lookup: field_name -> (default, version, description)
_SESSION_INFO_LOOKUP = {
    name: (default, ver, desc)
    for name, default, ver, desc in SESSION_INFO_FIELDS
}


# =========================================================================
# CLIENT -> SERVER: top-level decide_next_step() arguments
# =========================================================================
#
# These are the positional/keyword arguments the client passes.
# The graph entry point unpacks them into the initial state dict.

DECIDE_NEXT_STEP_ARGS = [
    ("log_content",        str,   ""),
    ("history",            list,  []),
    ("files",              list,  []),
    ("guidelines",         str,   ""),
    ("session_resolution", float, None),   # float or None
    ("session_info",       dict,  {}),
    ("abort_on_red_flags", bool,  True),
    ("abort_on_warnings",  bool,  False),
]


# =========================================================================
# SERVER -> CLIENT: history_record response fields
# =========================================================================
#
# Rules:
#   - New fields may be added freely (old clients use .get() / getattr()).
#   - Existing fields must NEVER be removed or renamed.
#   - Semantics of existing fields must not change.

RESPONSE_FIELDS = [
    "next_move",           # dict  -- see NEXT_MOVE_FIELDS
    "debug_log",           # list[str]
    "events",              # list[dict]
    "experiment_type",     # str
    "stop_reason",         # str or None
    "abort_message",       # str or None
    "red_flag_issues",     # list
    "warnings",            # list[str] -- deprecation / advisory messages
]

NEXT_MOVE_FIELDS = [
    "command",             # str: the phenix command to run
    "program",             # str: program name (e.g. "phenix.refine")
    "explanation",         # str: LLM reasoning text
    "process_log",         # str: detailed agent thought-process log
]


# =========================================================================
# CLIENT -> SERVER: history entry fields
# =========================================================================
#
# Each element of the `history` list should have these fields.

HISTORY_ENTRY_FIELDS = [
    ("cycle_number", int,  0),
    ("program",      str,  ""),
    ("command",      str,  ""),
    ("result",       str,  ""),
    ("output_files", list, []),
    ("analysis",     dict, {}),
]


# =========================================================================
# GRAPH STATE: fields copied from session_info to top-level state
# =========================================================================
#
# The graph entry point (REST handler or LocalAgent) copies these fields
# from session_info into top-level state keys so that graph nodes can
# access them as  state.get("directives")  rather than going through
# session_info.  This is documented here for clarity.
#
# If you add a new entry here, the graph entry point must be updated to
# copy it, AND the field must also be in SESSION_INFO_FIELDS above.

STATE_PROMOTED_FIELDS = [
    "directives",          # Copied to top-level for convenience
    "bad_inject_params",   # Copied to top-level for convenience
]


# =========================================================================
# RUNTIME HELPERS
# =========================================================================

def normalize_session_info(session_info):
    """
    Ensure all registered fields exist in session_info with safe defaults.

    Call this once at the top of perceive() to make the rest of the server
    code safe against old clients that omit newer fields.

    Mutable defaults (dict, list) are shallow-copied to avoid cross-request
    contamination.
    """
    import copy
    for field_name, default, _version, _desc in SESSION_INFO_FIELDS:
        if field_name not in session_info:
            # Copy mutable defaults so they aren't shared across calls
            if isinstance(default, (dict, list)):
                session_info[field_name] = copy.copy(default)
            else:
                session_info[field_name] = default
    return session_info


def check_client_version(session_info):
    """
    Protocol version gate.

    Returns:
        None                -- client is acceptable (continue normally)
        (str, str, dict)    -- (stop_reason, message, next_move) for rejection

    The caller (perceive or the REST handler) should check the return value
    and set state["stop"] = True with the returned next_move if not None.
    """
    client_version = session_info.get("client_protocol_version", 1)

    if client_version < MIN_SUPPORTED_PROTOCOL_VERSION:
        message = (
            "This version of PHENIX (AI Agent protocol v%d) is no longer "
            "supported by the server (minimum: v%d). "
            "Please update PHENIX to continue using the AI Agent."
            % (client_version, MIN_SUPPORTED_PROTOCOL_VERSION)
        )
        next_move = {
            "command": "STOP",
            "program": "STOP",
            "explanation": message,
            "process_log": "",
        }
        return ("unsupported_client", message, next_move)

    return None


def get_deprecation_warnings(session_info):
    """
    Return a list of warning strings for the client, or [].

    Warnings are advisory -- the server still processes the request.
    The client should print them to the log so the user sees them.
    """
    warnings = []
    client_version = session_info.get("client_protocol_version", 1)

    if client_version < CURRENT_PROTOCOL_VERSION:
        warnings.append(
            "Your PHENIX AI Agent uses protocol v%d (current: v%d). "
            "Consider updating PHENIX for the best experience."
            % (client_version, CURRENT_PROTOCOL_VERSION)
        )

    return warnings
