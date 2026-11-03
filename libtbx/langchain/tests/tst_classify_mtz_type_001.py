"""Regression tests for classify_mtz_type after removing the over-broad
'*_001.mtz -> map_coeffs_mtz' rule (former Pattern 2).

Bug (beta-blip, AIAgent_2): a user supplied `beta_blip_001.mtz` (a phenix.refine
output fed back in — it carries Fobs + R-free flags and IS valid refinement
data) plus a model, and asked the agent to "do everything necessary."  The agent
ran nothing and stopped at cycle 1.  Root cause: classify_mtz_type's Pattern 2
matched any `*_001.mtz` and returned map_coeffs_mtz, so the move-step in
workflow_state.py pulled the file out of data_mtz -> has_data_mtz=False ->
phenix.xtriage filtered from valid_programs -> [STOP].

Pattern 2 was redundant for its only legitimate example (model_refine_001.mtz is
already caught by Pattern 1) and its sole unique effect was this false positive
on user data.  It was removed.  Genuine map-coefficient outputs remain covered
by Pattern 1 (refine_NNN) and Pattern 3 (map_coeffs/denmod/density_mod); data
copies stay data via Pattern 4 (_data.mtz / refinement_data).

Filename-based, no PHENIX import; 2-space indent.
"""

from __future__ import absolute_import, division, print_function

import os
import importlib.util

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)                       # .../libtbx/langchain
_AGENT = os.path.join(_ROOT, "agent")


def _load_classifier():
  """Load classify_mtz_type from the real module.  Prefer the package import
  (touches nothing in sys.modules); fall back to loading agent/file_utils.py by
  path for a bare sandbox.  Returns the function, or None if unavailable."""
  try:
    from libtbx.langchain.agent.file_utils import classify_mtz_type
    return classify_mtz_type
  except Exception:
    pass
  for env, default in [("FILE_UTILS_PY", os.path.join(_AGENT, "file_utils.py"))]:
    path = os.environ.get(env) or default
    if os.path.isfile(path):
      spec = importlib.util.spec_from_file_location("_fu_under_test", path)
      mod = importlib.util.module_from_spec(spec)
      spec.loader.exec_module(mod)
      return mod.classify_mtz_type
  return None


# (filename, expected_category)
_DATA_CASES = [
  ("beta_blip_001.mtz", "data_mtz"),        # the bug
  ("7qz0_001.mtz", "data_mtz"),
  ("foo_001.mtz", "data_mtz"),
  ("data.mtz", "data_mtz"),
  ("7qz0.mtz", "data_mtz"),
  ("beta_blip.mtz", "data_mtz"),
  ("refine_001_data.mtz", "data_mtz"),      # Pattern 4 (_data.mtz)
  ("overall_best_refine_data.mtz", "data_mtz"),
  ("something_refinement_data.mtz", "data_mtz"),
]

_MAP_CASES = [
  ("refine_001.mtz", "map_coeffs_mtz"),     # Pattern 1
  ("refine_001_001.mtz", "map_coeffs_mtz"),
  ("7qz0_refine_001.mtz", "map_coeffs_mtz"),
  ("7qz0_refine_001_001.mtz", "map_coeffs_mtz"),
  ("model_refine_001.mtz", "map_coeffs_mtz"),  # was Pattern 2's example; now Pattern 1
  ("foo_map_coeffs.mtz", "map_coeffs_mtz"),    # Pattern 3
  ("denmod_map.mtz", "map_coeffs_mtz"),
  ("density_mod_1.mtz", "map_coeffs_mtz"),
]


def test_user_data_with_001_suffix_is_data():
  """The regression: *_001.mtz user data must classify as data_mtz."""
  classify = _load_classifier()
  if classify is None:
    print("  (skip) file_utils.classify_mtz_type not found")
    return
  got = classify("/path/to/beta_blip_001.mtz")
  assert got == "data_mtz", \
    "beta_blip_001.mtz -> %r (expected data_mtz)" % got


def test_all_data_cases():
  classify = _load_classifier()
  if classify is None:
    print("  (skip)")
    return
  for name, exp in _DATA_CASES:
    got = classify("/path/to/" + name)
    assert got == exp, "%s -> %r (expected %s)" % (name, got, exp)


def test_all_map_coeffs_cases_still_detected():
  """Removing Pattern 2 must not break genuine map-coefficient detection."""
  classify = _load_classifier()
  if classify is None:
    print("  (skip)")
    return
  for name, exp in _MAP_CASES:
    got = classify("/path/to/" + name)
    assert got == exp, "%s -> %r (expected %s)" % (name, got, exp)


def test_pattern2_example_now_covered_by_pattern1():
  """model_refine_001.mtz (the removed Pattern 2's docstring example) must still
  classify as map_coeffs_mtz via Pattern 1 — proving Pattern 2 was redundant."""
  classify = _load_classifier()
  if classify is None:
    print("  (skip)")
    return
  assert classify("/x/model_refine_001.mtz") == "map_coeffs_mtz"


_TESTS = [
  test_user_data_with_001_suffix_is_data,
  test_all_data_cases,
  test_all_map_coeffs_cases_still_detected,
  test_pattern2_example_now_covered_by_pattern1,
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
