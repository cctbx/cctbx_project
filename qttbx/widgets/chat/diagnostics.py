"""Diagnostics dump - copied/saved from Help > Diagnostics.

Pure-Python so it can also be invoked headlessly (debugging, bug reports).
"""

import datetime
import platform
import sys


def build_diagnostics(profile, storage, log_path=None, log_tail_lines=200):
  lines = []
  lines.append("PhenixChat diagnostics")
  lines.append("=" * 40)
  lines.append("timestamp: %s" % datetime.datetime.now().isoformat(
    timespec="seconds"))
  lines.append("platform:  %s" % platform.platform())
  lines.append("python:    %s" % sys.version.replace("\n", " "))

  lines.append("")
  lines.append("Versions:")
  lines.append("  anthropic: %s" % _safe_version("anthropic"))
  lines.append("  fastmcp:   %s" % _safe_version("fastmcp"))
  lines.append("  PySide6:   %s" % _safe_version("PySide6"))
  lines.append("  PySide2:   %s" % _safe_version("PySide2"))

  lines.append("")
  lines.append("Profile:")
  lines.append("  name:  %s" % getattr(profile, "name", "?"))
  lines.append("  model: %s" % getattr(profile, "model", "?"))

  lines.append("")
  lines.append("Storage:")
  lines.append("  project_dir: %s" % getattr(storage, "project_dir", "?"))
  lines.append("  chat_root:   %s" % getattr(storage, "chat_root", "?"))

  if log_path is not None:
    lines.append("")
    lines.append("Recent log (%s):" % log_path)
    lines.append("-" * 40)
    try:
      with open(log_path, "r") as fh:
        log_lines = fh.readlines()
      for line in log_lines[-log_tail_lines:]:
        lines.append(line.rstrip("\n"))
    except OSError:
      lines.append("(log file not readable)")
  return "\n".join(lines) + "\n"


def _safe_version(modname):
  try:
    m = __import__(modname)
  except ImportError:
    return "not installed"
  v = getattr(m, "__version__", None)
  return v or "(unknown)"
