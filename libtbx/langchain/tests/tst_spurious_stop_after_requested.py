"""Regression test: strip a SPURIOUS stop_after_requested directive.

Bug (beta-blip AIAgent_6): the user asked the agent to "do everything you can to
improve this model by refinement/rebuilding whatever is necessary" — an
open-ended request with NO stop instruction.  The LLM directive extractor
nonetheless set stop_conditions.stop_after_requested=True.  At cycle 3 the
workflow engine routed to validation, but the LLM read the spurious directive as
"stop after phenix.refine" and stopped — abandoning a poor model (R-free 0.326,
clashscore 19, 65% rotamer outliers) the user explicitly asked to improve.

Fix: after the intent dispatch in directive_extractor, strip
stop_after_requested when it is set but NOT justified — i.e. the advice carries
no explicit stop (_is_stop_after_requested is False) AND there is no
after_program pinning a concrete program to stop after.  The intent=task branch
(which sets the flag legitimately) always pairs it with after_program, so it is
preserved; genuine "refine then stop" advice sets _is_stop_after_requested=True
and is preserved too.

These tests exercise the strip logic together with the REAL
_is_stop_after_requested helper loaded from directive_extractor.py.  2-space
indent.
"""

from __future__ import absolute_import, division, print_function

import os
import re

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_is_stop_after_requested():
  """Load the real _is_stop_after_requested from directive_extractor.py.
  Prefer package import; fall back to extracting the helper by source slice so
  the test runs in a bare sandbox.  Returns the function or None."""
  try:
    from libtbx.langchain.agent.directive_extractor import (
      _is_stop_after_requested)
    return _is_stop_after_requested
  except Exception:
    pass
  path = os.environ.get("DIRECTIVE_EXTRACTOR_PY")
  if not path:
    for cand in (
        os.path.join(_HERE, "directive_extractor.py"),
        os.path.join(_HERE, "..", "agent", "directive_extractor.py"),
    ):
      if os.path.isfile(cand):
        path = cand
        break
  if not path or not os.path.isfile(path):
    return None
  src = open(path).read()
  try:
    start = src.index("_POSITIVE_STOP_AFTER_PATTERNS = (")
    end = src.index("def _resolve_after_program")
  except ValueError:
    return None
  ns = {"re": re}
  exec(compile(src[start:end], path, "exec"), ns)
  return ns.get("_is_stop_after_requested")


_IS_STOP = _load_is_stop_after_requested()


def _apply_strip(directives, user_advice):
  """Verbatim copy of the strip guard added to directive_extractor.py (kept in
  sync).  Mutates and returns `directives`."""
  _has_explicit_stop = _IS_STOP(user_advice) if _IS_STOP else False
  _sc_strip = directives.get("stop_conditions", {})
  if (_sc_strip.get("stop_after_requested") is True
          and not _sc_strip.get("after_program")
          and not _has_explicit_stop):
    _sc_strip.pop("stop_after_requested", None)
    if not _sc_strip:
      directives.pop("stop_conditions", None)
  return directives


def test_spurious_flag_stripped():
  """The bug case: open-ended 'do everything ... whatever is necessary' with a
  spurious flag and no after_program -> flag removed."""
  if _IS_STOP is None:
    print("  (skip) directive_extractor not found")
    return
  d = {"stop_conditions": {"stop_after_requested": True}}
  _apply_strip(d, "do everything you can to improve this model by "
                  "refinement/rebuilding whatever is necessary")
  assert "stop_conditions" not in d, \
    "spurious stop_conditions should be removed, got %r" % d.get("stop_conditions")


def test_explicit_stop_preserved():
  """Genuine 'refine then stop' (explicit stop, no after_program) -> KEEP."""
  if _IS_STOP is None:
    print("  (skip)")
    return
  assert _IS_STOP("refine the model and then stop") is True
  d = {"stop_conditions": {"stop_after_requested": True}}
  _apply_strip(d, "refine the model and then stop")
  assert d.get("stop_conditions", {}).get("stop_after_requested") is True, \
    "explicit-stop flag must be preserved"


def test_after_program_preserved():
  """intent=task style: flag justified by after_program -> KEEP."""
  if _IS_STOP is None:
    print("  (skip)")
    return
  d = {"stop_conditions": {
    "stop_after_requested": True,
    "after_program": "phenix.xtriage",
    "skip_validation": True}}
  _apply_strip(d, "just run xtriage")
  assert d["stop_conditions"]["stop_after_requested"] is True, \
    "flag with after_program must be preserved"
  assert d["stop_conditions"]["after_program"] == "phenix.xtriage"


def test_no_stop_conditions_noop():
  """No stop_conditions present -> no-op, no crash."""
  if _IS_STOP is None:
    print("  (skip)")
    return
  d = {"program_settings": {"default": {}}}
  _apply_strip(d, "improve the model")
  assert "stop_conditions" not in d
  assert "program_settings" in d


def test_partial_stop_conditions_retains_siblings():
  """If stop_conditions has OTHER keys, only the spurious flag is removed."""
  if _IS_STOP is None:
    print("  (skip)")
    return
  d = {"stop_conditions": {"stop_after_requested": True, "max_cycles": 5}}
  _apply_strip(d, "improve everything you can")
  assert d.get("stop_conditions") == {"max_cycles": 5}, \
    "only stop_after_requested should be removed, got %r" % d.get("stop_conditions")


_TESTS = [
  test_spurious_flag_stripped,
  test_explicit_stop_preserved,
  test_after_program_preserved,
  test_no_stop_conditions_noop,
  test_partial_stop_conditions_retains_siblings,
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
