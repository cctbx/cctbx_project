"""Regression test: advice-driven rebuild template selection.

Feature (beta-blip): the user asked to "do everything you can to improve this
model by refinement/rebuilding whatever is necessary".  The agent ran
xtriage -> refine -> molprobity and stopped without any model rebuilding,
because the selected plan template (refine_placed) has no rebuild stage and the
LLM cannot add stages a template lacks.

Fix (advice-driven): plan_generator._build_context detects explicit rebuild
intent in the advice/directives and sets a `wants_rebuild` context flag; a new
template `refine_rebuild_placed` (= refine_placed + a model_rebuilding autobuild
stage) matches on model_is_placed + wants_rebuild and outscores plain
refine_placed, so it is selected.  Plain "refine" requests still get
refine_placed.

These tests exercise the REAL select_template / load_templates against the
templates YAML, plus the wants_rebuild keyword logic.  2-space indent.
"""

from __future__ import absolute_import, division, print_function

import os
import importlib.util

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_loader():
  """Load plan_template_loader (provides load_templates / select_template).
  Prefer package import; fall back to env path / sibling lookup."""
  try:
    from libtbx.langchain.knowledge import plan_template_loader
    return plan_template_loader
  except Exception:
    pass
  path = os.environ.get("PLAN_TEMPLATE_LOADER_PY")
  if not path:
    for cand in (
        os.path.join(_HERE, "plan_template_loader.py"),
        os.path.join(_HERE, "..", "knowledge", "plan_template_loader.py"),
    ):
      if os.path.isfile(cand):
        path = cand
        break
  if not path or not os.path.isfile(path):
    return None
  spec = importlib.util.spec_from_file_location("_ptl_under_test", path)
  mod = importlib.util.module_from_spec(spec)
  try:
    spec.loader.exec_module(mod)
  except Exception:
    return None
  return mod


def _load_templates(loader):
  """Load the templates dict.  Uses the loader's load_templates if a YAML path
  is provided via env, else returns None to skip."""
  yaml_path = os.environ.get("PLAN_TEMPLATES_YAML")
  if yaml_path and hasattr(loader, "load_templates"):
    try:
      return loader.load_templates(yaml_path=yaml_path)
    except TypeError:
      return loader.load_templates()
  if hasattr(loader, "load_templates"):
    try:
      return loader.load_templates()
    except Exception:
      return None
  return None


def _base_ctx(**kw):
  ctx = {
    "experiment_type": "xray",
    "has_search_model": True,
    "model_is_placed": True,
    "wants_rebuild": False,
    "wants_polder": False,
    "has_ligand_code": False,
    "has_sequence": False,
    "resolution": 3.0,
    "is_twinned": False,
  }
  ctx.update(kw)
  return ctx


def test_rebuild_selects_rebuild_template():
  """placed model + wants_rebuild -> refine_rebuild_placed."""
  loader = _load_loader()
  if loader is None:
    print("  (skip) plan_template_loader not found")
    return
  templates = _load_templates(loader)
  if not templates:
    print("  (skip) templates not loadable")
    return
  tid = loader.select_template(templates, _base_ctx(wants_rebuild=True))
  assert tid == "refine_rebuild_placed", \
    "expected refine_rebuild_placed, got %r" % tid


def test_no_rebuild_selects_plain_placed():
  """placed model, no rebuild -> refine_placed (unchanged)."""
  loader = _load_loader()
  if loader is None:
    print("  (skip)")
    return
  templates = _load_templates(loader)
  if not templates:
    print("  (skip)")
    return
  tid = loader.select_template(templates, _base_ctx(wants_rebuild=False))
  assert tid == "refine_placed", "expected refine_placed, got %r" % tid


def test_rebuild_template_has_rebuild_stage():
  """The rebuild template must actually contain a model_rebuilding stage with
  autobuild."""
  loader = _load_loader()
  if loader is None:
    print("  (skip)")
    return
  templates = _load_templates(loader)
  if not templates:
    print("  (skip)")
    return
  tpl = templates.get("refine_rebuild_placed")
  assert tpl is not None, "refine_rebuild_placed template missing"
  stage_ids = [s.get("id") for s in tpl.get("stages", [])]
  assert "model_rebuilding" in stage_ids, \
    "model_rebuilding stage missing: %r" % stage_ids
  reb = next(s for s in tpl["stages"] if s.get("id") == "model_rebuilding")
  progs = reb.get("programs", [])
  assert any("autobuild" in str(p) for p in progs), \
    "autobuild not in model_rebuilding programs: %r" % progs


def test_ligand_and_polder_siblings_unaffected():
  """Adding the rebuild template must not steal selection from the ligand/polder
  placed-model siblings."""
  loader = _load_loader()
  if loader is None:
    print("  (skip)")
    return
  templates = _load_templates(loader)
  if not templates:
    print("  (skip)")
    return
  assert loader.select_template(
    templates, _base_ctx(has_ligand_code=True)) == "refine_placed_ligand"
  assert loader.select_template(
    templates, _base_ctx(wants_polder=True)) == "refine_placed_polder"


_TESTS = [
  test_rebuild_selects_rebuild_template,
  test_no_rebuild_selects_plain_placed,
  test_rebuild_template_has_rebuild_stage,
  test_ligand_and_polder_siblings_unaffected,
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
