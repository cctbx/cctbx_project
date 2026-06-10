"""Agent build version source-of-truth.

Reads the VERSION file at the langchain/ root.  Returns 'unknown'
if the file is missing.

Stdlib-only; no langchain deps so it imports cleanly anywhere
core/llm.py imports.

v119.H2.
"""
from __future__ import absolute_import, division, print_function

import os


def _version_file_path():
  """Return absolute path to the VERSION file at langchain/ root."""
  here = os.path.dirname(os.path.abspath(__file__))
  # core/ is a direct child of langchain/.
  return os.path.normpath(os.path.join(here, "..", "VERSION"))


def get_version():
  """Return the agent version string, or 'unknown' if absent.

  Reads on every call (no caching).  Cheap (small file, OS-cached
  after first read).  Avoids stale-import problems in development.

  Never raises.  Failure modes (all return "unknown"):
    - File missing
    - File unreadable (permissions, disk error)
    - File empty (or whitespace-only)
    - File contents undecodable as text (e.g. binary garbage)

  Returns:
    str: version string, or "unknown".
  """
  path = _version_file_path()
  try:
    with open(path) as f:
      return f.read().strip() or "unknown"
  except Exception:
    return "unknown"
