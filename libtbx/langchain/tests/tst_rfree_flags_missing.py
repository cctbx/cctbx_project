"""Regression tests for the rfree_flags_missing auto-recovery.

phenix.refine aborts with "No array of R-free flags found" when it is given
reflection data with no R-free flag array and neither a flags-bearing MTZ nor
xray_data.r_free_flags.generate=True (observed in AF-bromodomain-ligand,
AIAgent_15: refine ran on the raw 7qz0.mtz, failed, and the agent advanced to
polder on the unrefined model).

The fix makes this a *recoverable* error: a new recoverable_errors.yaml entry
(rfree_flags_missing, resolution: force_retry) plus two error_analyzer.py edits
(a force_retry dispatch + resolver, and a force_retry marker in
_extract_error_info so analyze() doesn't bail before dispatch).  On the forced
retry the server-side BUILD node / command_builder re-selects the locked R-free
MTZ (preserving R-free flag continuity) or falls back to generate=True.

These tests drive the REAL ErrorAnalyzer against the REAL recoverable_errors.yaml
(only the session object is a faithful stand-in), so the YAML entry, the
detection patterns, the force_retry dispatch, the _extract_error_info marker, and
the max_retries loop cap are all exercised end-to-end -- not mocked away.

Behavioral, no PHENIX import; 2-space indent.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import types
import importlib.util

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)                       # .../libtbx/langchain
_AGENT = os.path.join(_ROOT, "agent")
_KNOWLEDGE = os.path.join(_ROOT, "knowledge")

# The exact error text phenix.refine emits (from the AIAgent_15 log).
_RFREE_MISSING_LOG = (
  "No array of R-free flags found.\n"
  "For manual selection define:\n"
  "  miller_array.labels.name.test_flag_value\n"
  "  miller_array.labels.name.disable_suitability_test=True\n"
  "\n"
  "  [OldStyle] execution error: If previously used R-free flags are available "
  "run this command again\n"
  "with the name of the file containing the original flags as an\n"
  "additional input. If the structure was never refined before, or if the\n"
  "original R-free flags are unrecoverable, run this command again with\n"
  "the additional definition:\n"
  "\n"
  "xray_data.r_free_flags.generate=True\n"
)


def _find(env, *parts):
  e = os.environ.get(env)
  if e and os.path.isfile(e):
    return e
  p = os.path.join(*parts)
  return p if os.path.isfile(p) else None


class _FakeSession(object):
  """Minimal stand-in for AgentSession: ErrorAnalyzer only touches .data."""

  def __init__(self):
    self.data = {}


def _load_analyzer():
  """Load the REAL ErrorAnalyzer.

  In a PHENIX environment the module is importable directly as
  libtbx.langchain.agent.error_analyzer -- we use that and touch nothing else,
  so running inside the full suite leaves sys.modules exactly as it found it.

  Only when the normal import is unavailable (e.g. a bare sandbox) do we fall
  back to loading error_analyzer.py from disk by path.  That fallback inserts
  temporary parent-package entries into sys.modules and is careful to RESTORE
  sys.modules to its prior state afterward, so it cannot pollute later tests
  (the bug that made tst_dependencies.py fail only inside the suite).

  Its _load_config reads ../knowledge/recoverable_errors.yaml relative to its
  own __file__, so the on-disk path picks up the REAL YAML automatically.
  Returns an ErrorAnalyzer instance, or None if prerequisites are missing."""
  # ErrorAnalyzer needs PyYAML to load its config; skip gracefully if absent.
  try:
    import yaml
    del yaml
  except ImportError:
    return None

  # Preferred path: the real package import (no sys.modules manipulation).
  try:
    from libtbx.langchain.agent.error_analyzer import ErrorAnalyzer
    return ErrorAnalyzer()
  except Exception:
    pass

  # Fallback: load from disk by path, restoring sys.modules afterward.
  ea_path = _find("ERROR_ANALYZER_PY", _AGENT, "error_analyzer.py")
  yaml_path = _find("RECOVERABLE_ERRORS_YAML", _KNOWLEDGE,
                    "recoverable_errors.yaml")
  if not (ea_path and yaml_path):
    return None

  injected = []          # names we add (to delete on restore)
  saved = {}             # names we overwrite (to restore on restore)
  names = ["libtbx", "libtbx.langchain", "libtbx.langchain.agent",
           "libtbx.langchain.knowledge",
           "libtbx.langchain.agent.error_analyzer"]
  for name in names:
    if name in sys.modules:
      saved[name] = sys.modules[name]
    else:
      injected.append(name)

  try:
    for name, path in [("libtbx", None), ("libtbx.langchain", None),
                       ("libtbx.langchain.agent", _AGENT),
                       ("libtbx.langchain.knowledge", _KNOWLEDGE)]:
      if name not in sys.modules:
        m = types.ModuleType(name)
        if path:
          m.__path__ = [path]
        sys.modules[name] = m

    spec = importlib.util.spec_from_file_location(
      "libtbx.langchain.agent.error_analyzer", ea_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["libtbx.langchain.agent.error_analyzer"] = mod
    spec.loader.exec_module(mod)
    return mod.ErrorAnalyzer()
  finally:
    # Restore sys.modules exactly as we found it — delete what we added,
    # put back what we overwrote.  Prevents suite-wide cross-test pollution.
    for name in injected:
      sys.modules.pop(name, None)
    for name, original in saved.items():
      sys.modules[name] = original


# --- Tests ------------------------------------------------------------------

def test_yaml_entry_present_and_well_formed():
  """The rfree_flags_missing entry exists with the expected force_retry shape."""
  analyzer = _load_analyzer()
  if analyzer is None:
    print("  (skip) error_analyzer.py / recoverable_errors.yaml not found")
    return
  errs = analyzer._config.get("errors", {})
  assert "rfree_flags_missing" in errs, \
    "rfree_flags_missing entry missing from recoverable_errors.yaml"
  e = errs["rfree_flags_missing"]
  assert e.get("resolution") == "force_retry", \
    "rfree_flags_missing must use resolution: force_retry"
  assert e.get("retry_program") == "phenix.refine", \
    "rfree_flags_missing must retry phenix.refine"
  assert e.get("max_retries") == 1, "expected max_retries: 1"


def test_full_analyze_returns_force_retry_recovery():
  """End-to-end analyze() returns a force-retry recovery (proves the
  _extract_error_info marker -- without it analyze() bails before dispatch)."""
  analyzer = _load_analyzer()
  if analyzer is None:
    print("  (skip) error_analyzer.py / recoverable_errors.yaml not found")
    return
  rec = analyzer.analyze(_RFREE_MISSING_LOG, "phenix.refine",
                         context={}, session=_FakeSession())
  assert rec is not None, \
    "analyze() returned None -- _extract_error_info marker missing?"
  assert rec.error_type == "rfree_flags_missing", \
    "wrong error_type: %r" % (rec.error_type,)
  assert rec.retry_program == "phenix.refine", \
    "wrong retry_program: %r" % (rec.retry_program,)
  # force_retry carries NO flag additions and NO flag strips -- the builder
  # does the file re-selection on the rebuilt retry.
  assert rec.flags == {}, "force_retry must not add flags: %r" % (rec.flags,)
  assert list(rec.strip_flags) == [], \
    "force_retry must not strip flags: %r" % (rec.strip_flags,)


def test_detection_is_specific():
  """A benign refine log must not trigger the recovery (no false positives)."""
  analyzer = _load_analyzer()
  if analyzer is None:
    print("  (skip)")
    return
  benign = "Refinement converged. R-free = 0.21, R-work = 0.18. Done."
  rec = analyzer.analyze(benign, "phenix.refine",
                         context={}, session=_FakeSession())
  assert rec is None, "false positive on benign log: %r" % (rec,)


def test_does_not_collide_with_rfree_flags_mismatch():
  """The two R-free entries are disjoint: the 'mismatch' text must NOT match
  rfree_flags_missing, and the 'missing' text must NOT match the mismatch entry."""
  analyzer = _load_analyzer()
  if analyzer is None:
    print("  (skip)")
    return
  # mismatch text -> mismatch entry (strip), never our entry
  mismatch_rec = analyzer.analyze(
    "Sorry: Please resolve the R-free flags mismatch.",
    "phenix.refine", context={}, session=_FakeSession())
  assert mismatch_rec is not None, "mismatch entry should still resolve"
  assert mismatch_rec.error_type == "rfree_flags_mismatch", \
    "mismatch text matched wrong entry: %r" % (mismatch_rec.error_type,)
  # and the mismatch path is a strip (regression guard on the twin)
  assert list(mismatch_rec.strip_flags), \
    "rfree_flags_mismatch must still carry strip_flags"
  # our 'missing' text -> our entry, never the mismatch entry
  missing_rec = analyzer.analyze(_RFREE_MISSING_LOG, "phenix.refine",
                                 context={}, session=_FakeSession())
  assert missing_rec.error_type == "rfree_flags_missing", \
    "missing text matched wrong entry: %r" % (missing_rec.error_type,)


def test_loop_cap_blocks_second_identical_failure():
  """max_retries: 1 -- the first failure recovers, the second returns None so
  the agent proceeds to normal failure handling rather than looping forever."""
  analyzer = _load_analyzer()
  if analyzer is None:
    print("  (skip)")
    return
  sess = _FakeSession()
  first = analyzer.analyze(_RFREE_MISSING_LOG, "phenix.refine",
                           context={}, session=sess)
  assert first is not None, "first attempt should recover"
  second = analyzer.analyze(_RFREE_MISSING_LOG, "phenix.refine",
                            context={}, session=sess)
  assert second is None, \
    "second identical failure must be blocked by max_retries: 1"


_TESTS = [
  test_yaml_entry_present_and_well_formed,
  test_full_analyze_returns_force_retry_recovery,
  test_detection_is_specific,
  test_does_not_collide_with_rfree_flags_mismatch,
  test_loop_cap_blocks_second_identical_failure,
]


def run_all_tests():
  for fn in _TESTS:
    fn()
  print("All %d tests passed." % len(_TESTS))
  return True


if __name__ == "__main__":
  p = f = 0
  for fn in _TESTS:
    try:
      fn()
      print("  PASS: %s" % fn.__name__)
      p += 1
    except AssertionError as e:
      print("  FAIL: %s -- %s" % (fn.__name__, e))
      f += 1
  print("\n%d passed, %d failed" % (p, f))
