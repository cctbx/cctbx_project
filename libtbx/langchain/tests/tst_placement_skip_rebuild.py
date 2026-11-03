"""Regression test: placement-skip must NOT suppress an explicitly-requested
rebuild stage.

Bug (beta-blip AIAgent_17): the user chose the refine_rebuild_placed template
(wants_rebuild), so the plan contained a model_rebuilding (autobuild) stage.  But
a placement heuristic in ai_agent.py — which skips model_rebuilding when the model
is already decent (R-free < 0.35) because "autobuild on a good model is
destructive" — fired and marked the stage skipped, defeating the user's explicit
request.  (beta-blip R-free = 0.326 < 0.35.)

Fix: gate the `_skip_ids.add("model_rebuilding")` on the template NOT being
`refine_rebuild_placed`.  molecular_replacement / experimental_phasing are still
skipped (those are genuinely destructive on a placed model); only the
user-requested rebuild is preserved.

This test reproduces the skip-set construction logic in isolation (it does not
import the GUI program), exercising the exact condition that was changed.
2-space indent.
"""

from __future__ import absolute_import, division, print_function


def _compute_skip_ids(template_id, ev_metric, ev_value):
  """Mirror of the placement-skip skip-set logic in ai_agent.py (kept in sync).
  Returns the set of stage ids that would be skipped."""
  skip_ids = {"molecular_replacement", "experimental_phasing"}
  explicit_rebuild = (template_id == "refine_rebuild_placed")
  if ev_value is not None and not explicit_rebuild:
    try:
      v = float(ev_value)
      if ((ev_metric == "R-free" and v < 0.35)
          or (ev_metric == "R-work" and v < 0.30)
          or (ev_metric in ("CC", "Map CC") and v > 0.5)):
        skip_ids.add("model_rebuilding")
    except (ValueError, TypeError):
      pass
  return skip_ids


def test_explicit_rebuild_preserved_at_low_rfree():
  """refine_rebuild_placed + R-free 0.326 -> model_rebuilding NOT skipped."""
  skip = _compute_skip_ids("refine_rebuild_placed", "R-free", 0.326)
  assert "model_rebuilding" not in skip, \
    "explicit rebuild must be preserved, got skip=%r" % skip
  # MR / phasing still suppressed on a placed model
  assert "molecular_replacement" in skip
  assert "experimental_phasing" in skip


def test_non_rebuild_template_still_skips_at_low_rfree():
  """A placed-model template that is NOT the rebuild template still skips
  model_rebuilding when the model is already decent (unchanged behavior)."""
  skip = _compute_skip_ids("refine_placed", "R-free", 0.326)
  assert "model_rebuilding" in skip, \
    "non-rebuild template should still skip model_rebuilding at R-free<0.35"


def test_rebuild_template_preserved_even_with_good_cc():
  """Explicit rebuild preserved regardless of which placement metric is good."""
  for metric, value in (("R-free", 0.30), ("R-work", 0.25), ("CC", 0.8)):
    skip = _compute_skip_ids("refine_rebuild_placed", metric, value)
    assert "model_rebuilding" not in skip, \
      "rebuild must survive %s=%s, got %r" % (metric, value, skip)


def test_high_rfree_never_skips_rebuild():
  """When the model is poor (R-free >= 0.35), model_rebuilding is not skipped for
  ANY template (the heuristic only suppresses on a good model)."""
  for tpl in ("refine_placed", "refine_rebuild_placed"):
    skip = _compute_skip_ids(tpl, "R-free", 0.50)
    assert "model_rebuilding" not in skip, \
      "%s at R-free 0.50 should not skip rebuild, got %r" % (tpl, skip)


def test_no_evidence_no_rebuild_skip():
  """No placement evidence -> model_rebuilding not skipped (only MR/phasing)."""
  skip = _compute_skip_ids("refine_placed", "", None)
  assert "model_rebuilding" not in skip
  assert skip == {"molecular_replacement", "experimental_phasing"}


_TESTS = [
  test_explicit_rebuild_preserved_at_low_rfree,
  test_non_rebuild_template_still_skips_at_low_rfree,
  test_rebuild_template_preserved_even_with_good_cc,
  test_high_rfree_never_skips_rebuild,
  test_no_evidence_no_rebuild_skip,
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
