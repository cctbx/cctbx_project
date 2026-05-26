"""Shared helpers for the v119.H3 startup canary tools.

Used by:
  - tests/tst_canary.py            (K_H3a metadata canary)
  - tests/llm/canary_check.py      (H3b LLM smoke canary orchestrator)

Both consumers call load_canary_expected() to read the pinned
agent_version + defaults_fingerprint from tests/canary_expected.json.

Stdlib-only.  No langchain deps.  Never raises.
"""
from __future__ import absolute_import, division, print_function

import json
import os


def canary_config_path():
    """Return absolute path to tests/canary_expected.json.

    This module lives at langchain/tests/canary_utils.py, so its
    sibling is the expected location for canary_expected.json.

    Returns the most likely path even if the file doesn't exist
    (callers check existence via load_canary_expected).
    """
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(here, "canary_expected.json"))


def load_canary_expected():
    """Return (config_dict, error_msg).  Never raises.

    Returns (dict, None) on success.  Returns (None, msg) on any
    load failure -- file missing, invalid JSON, missing keys,
    wrong-typed keys.

    Validates that the file is a dict with two string-typed keys:
      - "agent_version"
      - "defaults_fingerprint"

    Both keys are required.  Extra keys are allowed and ignored.
    """
    path = canary_config_path()
    if not os.path.isfile(path):
        return None, "canary_expected.json not found at %s" % path
    try:
        with open(path) as f:
            cfg = json.load(f)
    except Exception as e:
        return None, "Failed to parse %s: %s" % (path, e)

    if not isinstance(cfg, dict):
        return None, ("canary_expected.json must contain a JSON "
                      "object; got %s" % type(cfg).__name__)

    required = ("agent_version", "defaults_fingerprint")
    missing = [k for k in required if k not in cfg]
    if missing:
        return None, ("canary_expected.json missing required keys: %s"
                      % ", ".join(missing))
    for k in required:
        if not isinstance(cfg[k], str):
            return None, ("canary_expected.json key %r must be str, "
                          "got %s" % (k, type(cfg[k]).__name__))
    return cfg, None
