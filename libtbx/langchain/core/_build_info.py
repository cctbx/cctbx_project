"""Agent build info -- source for the agent_build response field.

v119.H2.  Two public functions:

- get_agent_build_info() -- returns the three-key dict that becomes
  response["agent_build"].
- inject_agent_build(response) -- mutates a response dict in place
  to add the agent_build field.  Called from
  _build_group_args_response in run_ai_agent.py for every response
  that exits the server.

Both functions never raise.

Module-load time is captured here as _AGENT_STARTED_AT.  In
practice this is approximately equal to server-process-start
because this module is imported during run_ai_agent.py module load,
which itself happens at server boot.
"""
from __future__ import absolute_import, division, print_function

import datetime


# Module-load time, captured once.  Strict UTC ISO 8601:
#   YYYY-MM-DDTHH:MM:SSZ -- no microseconds, no offset.
# Pinning the format here keeps it stable for log aggregators,
# dashboards, and uptime calculations downstream.
_AGENT_STARTED_AT = datetime.datetime.now(
    datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# Fallback dict for degraded environments.  All three subfields
# present as strings -- syntactically valid, semantically "unknown".
_FALLBACK_BUILD_INFO = {
  "version": "unknown",
  "defaults_fingerprint": "",
  "started_at": "",
}


def get_agent_build_info():
  """Return the agent_build dict for the response.

  Three keys, all always strings:
    - version: from core/_version.get_version()
               ("unknown" if VERSION absent)
    - defaults_fingerprint: from core.llm.compute_defaults_fingerprint()
                            ("" if compute fails)
    - started_at: module-load time, strict ISO 8601 UTC

  Never raises.  Each underlying call is guarded.  A degraded
  environment still produces a syntactically valid dict.

  Returns:
    dict: with exactly three keys, all str-valued.
  """
  # Defensive imports.  Each failure -> safe fallback function.
  try:
    from libtbx.langchain.core._version import get_version
  except ImportError:
    try:
      from core._version import get_version
    except ImportError:
      def get_version():
        return "unknown"

  try:
    from libtbx.langchain.core.llm import (
      compute_defaults_fingerprint)
  except ImportError:
    try:
      from core.llm import compute_defaults_fingerprint
    except ImportError:
      def compute_defaults_fingerprint():
        return ""

  try:
    version = get_version()
  except Exception:
    version = "unknown"
  try:
    fingerprint = compute_defaults_fingerprint()
  except Exception:
    fingerprint = ""

  return {
    "version": version,
    "defaults_fingerprint": fingerprint,
    "started_at": _AGENT_STARTED_AT,
  }


def inject_agent_build(response):
  """Mutate a response dict in place to add the agent_build field.

  v119.H2 post-processor injection point.  Called from
  _build_group_args_response in run_ai_agent.py for every response
  that exits the server.

  Mutation is intentional and safe in the call context:
    - Every response dict is built fresh per request by
      create_response (or its variants); no upstream caller
      reuses it.
    - Other mutations of the same dict already happen in
      run_ai_agent.py (timing_ms, red_flags, thinking metadata).
    - The mutation is the LAST mutation before the response is
      serialized to JSON (via serialize_response) inside
      _build_group_args_response.

  Overwrites any pre-existing response["agent_build"] (so the
  {} that apply_response_defaults sets during construction is
  replaced with the populated dict here).

  Never raises.  If get_agent_build_info() somehow fails, falls
  back to a structurally-valid "unknown" dict so the response
  build can proceed.

  Args:
    response: dict.  Modified in place.

  Returns:
    None.  The mutation is the side effect.
  """
  try:
    info = get_agent_build_info()
  except Exception:
    info = dict(_FALLBACK_BUILD_INFO)
  if not isinstance(response, dict):
    # Defensive: if somehow response isn't a dict, do nothing rather
    # than raise.  The caller will detect the type mismatch elsewhere.
    return
  response["agent_build"] = info
