"""Regression test: phenix.refine R-work/R-free must come from the FINAL line.

Bug (beta-blip AIAgent_3): phenix.refine emits both
    Start R-work = 0.2672, R-free = 0.2579
    Final R-work = 0.1639, R-free = 0.3258
extract_all_metrics() ran the generic YAML patterns first (STEP 1), which use
first-match semantics and captured the *Start* values, then merged the
Final-anchored _extract_refine_metrics() result with `if k not in metrics` — so
the correct Final values were discarded.  The agent recorded R-free 0.258 (the
pre-refinement value), concluded "target reached" (0.258 < 0.33 at 3.0 A), and
stopped before validation — when the true R-free was 0.326 (overfit).

Fix: in the refine branch of extract_all_metrics(), when a "Final R-work = ...,
R-free = ..." line is present, let _extract_refine_metrics override the YAML
r_work/r_free (otherwise fall back to fill-if-missing).

These tests drive the real extract_all_metrics() with a stubbed YAML extractor
that reproduces STEP 1's first-match behavior, so they exercise the actual merge
logic.  2-space indent.
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys
import types
import importlib.util

_HERE = os.path.dirname(os.path.abspath(__file__))


def _install_yaml_stub():
  """Install a stub libtbx.langchain.knowledge.metric_patterns that mimics
  STEP 1's broad first-match behavior.  Must stay installed while
  extract_all_metrics() runs (it imports metric_patterns lazily, at call time).
  Returns a dict of saved modules for later restoration."""
  saved = {}
  for name in ("libtbx", "libtbx.langchain",
               "libtbx.langchain.knowledge",
               "libtbx.langchain.knowledge.metric_patterns"):
    saved[name] = sys.modules.get(name)

  def _stub_yaml(log_text, program):
    d = {}
    rw = re.search(r'R-work\s*[:=]\s*([0-9.]+)', log_text, re.IGNORECASE)
    rf = re.search(r'R-free\s*[:=]\s*([0-9.]+)', log_text, re.IGNORECASE)
    if rw:
      d["r_work"] = float(rw.group(1))
    if rf:
      d["r_free"] = float(rf.group(1))
    return d

  pkg = types.ModuleType("libtbx")
  pkg.__path__ = []
  lc = types.ModuleType("libtbx.langchain")
  lc.__path__ = []
  kn = types.ModuleType("libtbx.langchain.knowledge")
  kn.__path__ = []
  mp = types.ModuleType("libtbx.langchain.knowledge.metric_patterns")
  mp.extract_metrics_for_program = _stub_yaml
  sys.modules["libtbx"] = pkg
  sys.modules["libtbx.langchain"] = lc
  sys.modules["libtbx.langchain.knowledge"] = kn
  sys.modules["libtbx.langchain.knowledge.metric_patterns"] = mp
  return saved


def _restore_modules(saved):
  for name, val in saved.items():
    if val is None:
      sys.modules.pop(name, None)
    else:
      sys.modules[name] = val


def _load_with_stubbed_yaml():
  """Load log_parsers and return extract_all_metrics, or None if not found.
  The caller is responsible for keeping the YAML stub installed (via
  _install_yaml_stub) around each extract() call."""
  path = os.environ.get("LOG_PARSERS_PY")
  if not path:
    for cand in (
        os.path.join(_HERE, "log_parsers.py"),
        os.path.join(_HERE, "..", "phenix_ai", "log_parsers.py"),
        os.path.join(_HERE, "..", "..", "phenix", "phenix_ai", "log_parsers.py"),
    ):
      if os.path.isfile(cand):
        path = cand
        break
  if not path or not os.path.isfile(path):
    return None
  saved = _install_yaml_stub()
  try:
    spec = importlib.util.spec_from_file_location("_log_parsers_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.extract_all_metrics
  finally:
    _restore_modules(saved)


def _extract(text, program):
  """Run extract_all_metrics with the YAML stub installed for the call."""
  fn = _load_with_stubbed_yaml()
  if fn is None:
    return None
  saved = _install_yaml_stub()
  try:
    return fn(text, program)
  finally:
    _restore_modules(saved)


_START_AND_FINAL = (
  "Start R-work = 0.2672, R-free = 0.2579\n"
  "Final R-work = 0.1639, R-free = 0.3258\n")


def test_final_overrides_start():
  """The bug case: both Start and Final present -> use Final."""
  m = _extract(_START_AND_FINAL, "phenix.refine")
  if m is None:
    print("  (skip) log_parsers.py not found")
    return
  assert m.get("r_free") == 0.3258, "r_free=%r (expected 0.3258 Final)" % m.get("r_free")
  assert m.get("r_work") == 0.1639, "r_work=%r (expected 0.1639 Final)" % m.get("r_work")


def test_final_only():
  """Only a Final line -> those values."""
  m = _extract("Final R-work = 0.1856, R-free = 0.2234", "phenix.refine")
  if m is None:
    print("  (skip)")
    return
  assert m.get("r_free") == 0.2234, "r_free=%r" % m.get("r_free")
  assert m.get("r_work") == 0.1856, "r_work=%r" % m.get("r_work")


def test_no_final_line_fallback():
  """No Final line -> fallback to available values, no crash, no spurious data."""
  m = _extract("R-work = 0.21, R-free = 0.25", "phenix.refine")
  if m is None:
    print("  (skip)")
    return
  assert m.get("r_free") == 0.25, "r_free=%r" % m.get("r_free")
  assert m.get("r_work") == 0.21, "r_work=%r" % m.get("r_work")


def test_real_space_refine_unaffected():
  """real_space_refine routes through a different branch; the refine-branch
  override must not change its behavior."""
  m = _extract(_START_AND_FINAL, "phenix.real_space_refine")
  if m is None:
    print("  (skip)")
    return
  # The real_space branch does not apply the refine override; we only assert it
  # returns a value and does not raise.
  assert m.get("r_free") is not None


_TESTS = [
  test_final_overrides_start,
  test_final_only,
  test_no_final_line_fallback,
  test_real_space_refine_unaffected,
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
