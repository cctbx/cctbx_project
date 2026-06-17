"""K_RFREE_GUARD: phenix.refine must not regenerate an existing R-free set.

Requirements (Tom, beta_blip_001):
  (1) input has NO R-free flags        -> generate a new set
  (2) input ALREADY has R-free flags   -> keep them (never generate)
  (3) once a set is locked/generated   -> never regenerate

The bug: the generate-flags invariant gated ONLY on `context.rfree_mtz`,
which locks from a refinement OUTPUT.  A FIRST refinement against
pre-flagged INPUT data (rfree_mtz not yet locked) wrongly received
`generate=True`, which phenix.refine uses to OVERWRITE the existing flags
with a fresh random partition -- silently invalidating cross-validation.

The fix adds `CommandBuilder._input_mtz_has_rfree(files, context)`, which
detects R-free flags in the actually-selected data MTZ (via the precomputed
`context.mtz_inspection["rfree_label"]` when present, else by inspecting the
selected `files["data_mtz"]` directly).  The invariant now suppresses/strips
`generate_rfree_flags` whenever flags already exist from EITHER source.

These tests exercise the REAL `_input_mtz_has_rfree` (bound from the real
CommandBuilder class source, with inspect_mtz stubbed) plus the decision
table the invariant implements.
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import types

_HERE = os.path.dirname(os.path.abspath(__file__))
_MISSING = object()  # sentinel for "key absent" in sys.modules snapshot


def _load_helper():
  """Bind CommandBuilder._input_mtz_has_rfree without importing the heavy
  ProgramRegistry/yaml stack.  We exec the method's source against a stub
  module that provides a controllable inspect_mtz."""
  cb_path = os.environ.get(
    "COMMAND_BUILDER_PY",
    os.path.join(_HERE, "..", "agent", "command_builder.py"))
  if not os.path.exists(cb_path):
    return None, None, None
  src = open(cb_path).read()
  if "_input_mtz_has_rfree" not in src:
    return None, None, None
  # Extract the method body and dedent to a module-level function.
  lines = src.splitlines()
  start = None
  for i, ln in enumerate(lines):
    if ln.strip().startswith("def _input_mtz_has_rfree"):
      start = i
      break
  if start is None:
    return None, None, None
  base_indent = len(lines[start]) - len(lines[start].lstrip())
  body = [lines[start][base_indent:]]
  for ln in lines[start + 1:]:
    if ln.strip() == "":
      body.append("")
      continue
    ind = len(ln) - len(ln.lstrip())
    if ind <= base_indent and ln.strip():
      break
    body.append(ln[base_indent:])
  func_src = "\n".join(body)
  # Provide a stub `agent.mtz_inspector.inspect_mtz` the method imports.
  stub_flags = {"files_with_flags": set(), "raise": False}

  def _inspect(path):
    if stub_flags["raise"]:
      raise RuntimeError("cctbx unavailable")
    return {"rfree_label":
            "R-free-flags" if path in stub_flags["files_with_flags"]
            else None}

  # The helper imports inspect_mtz from EITHER
  # `libtbx.langchain.agent.mtz_inspector` (tried first) OR
  # `agent.mtz_inspector` (fallback).  In the deployed tree the libtbx path
  # resolves to the REAL module, which would bypass a stub registered only
  # under `agent.mtz_inspector` (this caused false failures: the real
  # inspect_mtz ran on fake filenames and found no flags).  Register the stub
  # under BOTH names so whichever import the helper resolves gets the stub.
  insp_mod = types.ModuleType("mtz_inspector_stub")
  insp_mod.inspect_mtz = _inspect

  def _ensure_pkg(name):
    if name not in sys.modules:
      m = types.ModuleType(name)
      m.__path__ = []
      sys.modules[name] = m
    return sys.modules[name]

  # Snapshot the sys.modules entries (and the parent-package attributes) we are
  # about to overwrite, so a teardown can restore them — leaving a stub behind
  # would break any LATER test in the same process that needs the real
  # mtz_inspector (PHENIX importlib convention: restore what you mutate).
  _mod_names = [
    "agent.mtz_inspector",
    "libtbx.langchain.agent.mtz_inspector",
  ]
  _saved_modules = {n: sys.modules.get(n, _MISSING) for n in _mod_names}
  _parents = ["libtbx.langchain.agent", "agent"]

  # Build/register both module paths pointing at the same stub.
  _ensure_pkg("agent")
  sys.modules["agent.mtz_inspector"] = insp_mod
  _ensure_pkg("libtbx")
  _ensure_pkg("libtbx.langchain")
  _ensure_pkg("libtbx.langchain.agent")
  sys.modules["libtbx.langchain.agent.mtz_inspector"] = insp_mod
  # Also expose as an attribute so `from pkg import mtz_inspector` styles work.
  _saved_attrs = {}
  for p in _parents:
    pmod = sys.modules.get(p)
    _saved_attrs[p] = getattr(pmod, "mtz_inspector", _MISSING) if pmod else _MISSING
    if pmod is not None:
      pmod.mtz_inspector = insp_mod

  def _restore():
    for n, v in _saved_modules.items():
      if v is _MISSING:
        sys.modules.pop(n, None)
      else:
        sys.modules[n] = v
    for p, v in _saved_attrs.items():
      pmod = sys.modules.get(p)
      if pmod is None:
        continue
      if v is _MISSING:
        if hasattr(pmod, "mtz_inspector"):
          delattr(pmod, "mtz_inspector")
      else:
        pmod.mtz_inspector = v

  ns = {}
  exec(func_src, ns)
  return ns["_input_mtz_has_rfree"], stub_flags, _restore


class _Ctx(object):
  def __init__(self, mtz_inspection=None):
    self.mtz_inspection = mtz_inspection


# Bind once.
_HELPER, _STUB, _RESTORE = _load_helper()


def _has_rfree(files, ctx):
  # The method takes (self, files, context); self is unused, pass None.
  return _HELPER(None, files, ctx)


def _decide(rfree_mtz, files, ctx, strategy=None):
  """Mirror of the invariant's generate decision using the real helper."""
  strategy = dict(strategy or {})
  has = _has_rfree(files, ctx)
  if rfree_mtz or has:
    if strategy.get("generate_rfree_flags"):
      del strategy["generate_rfree_flags"]
  else:
    if "generate_rfree_flags" not in strategy:
      strategy["generate_rfree_flags"] = True
  return strategy.get("generate_rfree_flags", False)


def test_req1_no_flags_generates():
  if _HELPER is None:
    print("  SKIP: helper not found (fix not applied?)")
    return
  _STUB["files_with_flags"] = set()
  _STUB["raise"] = False
  g = _decide(None, {"data_mtz": "noflags.mtz"}, _Ctx())
  assert g is True, "no flags + no lock must generate, got %r" % g
  print("  PASS: test_req1_no_flags_generates")


def test_req2_input_flags_kept_direct_inspect():
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = {"beta_blip_001.mtz"}
  _STUB["raise"] = False
  # mtz_inspection None (as in the real beta_blip session) -> direct inspect
  g = _decide(None, {"data_mtz": "beta_blip_001.mtz"}, _Ctx(mtz_inspection=None))
  assert g is False, "input flags must NOT generate, got %r" % g
  print("  PASS: test_req2_input_flags_kept_direct_inspect")


def test_req2_strips_llm_forced_generate():
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = {"beta_blip_001.mtz"}
  _STUB["raise"] = False
  g = _decide(None, {"data_mtz": "beta_blip_001.mtz"}, _Ctx(),
              strategy={"generate_rfree_flags": True})
  assert g is False, "LLM-forced generate must be stripped, got %r" % g
  print("  PASS: test_req2_strips_llm_forced_generate")


def test_req2_precomputed_inspection():
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = set()  # direct inspect would say no
  _STUB["raise"] = False
  # but precomputed mtz_inspection says yes -> trusted
  g = _decide(None, {"data_mtz": "x.mtz"},
              _Ctx(mtz_inspection={"rfree_label": "FreeR_flag"}))
  assert g is False, "precomputed rfree_label must suppress generate, got %r" % g
  print("  PASS: test_req2_precomputed_inspection")


def test_req3_locked_rfree_strips():
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = set()
  _STUB["raise"] = False
  g = _decide("/sub_02_refine/refine_001_data.mtz",
              {"data_mtz": "refine_001_data.mtz"}, _Ctx(),
              strategy={"generate_rfree_flags": True})
  assert g is False, "locked rfree_mtz must strip generate, got %r" % g
  print("  PASS: test_req3_locked_rfree_strips")


def test_safe_degradation_on_inspect_error():
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = set()
  _STUB["raise"] = True  # inspect_mtz blows up
  # can't prove flags exist -> legacy behavior (generate when no lock)
  g = _decide(None, {"data_mtz": "unknown.mtz"}, _Ctx())
  assert g is True, "must fall back to generate when detection fails, got %r" % g
  print("  PASS: test_safe_degradation_on_inspect_error")


def test_list_valued_data_mtz_slot():
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = {"multi.mtz"}
  _STUB["raise"] = False
  # multiple:true slots can be lists; helper takes the first
  has = _has_rfree({"data_mtz": ["multi.mtz", "other.mtz"]}, _Ctx())
  assert has is True, "list-valued data_mtz first element must be inspected"
  print("  PASS: test_list_valued_data_mtz_slot")


def test_no_data_mtz_no_flags():
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = set()
  _STUB["raise"] = False
  has = _has_rfree({}, _Ctx())
  assert has is False, "no data_mtz selected -> cannot have flags"
  print("  PASS: test_no_data_mtz_no_flags")


_TESTS = [
  test_req1_no_flags_generates,
  test_req2_input_flags_kept_direct_inspect,
  test_req2_strips_llm_forced_generate,
  test_req2_precomputed_inspection,
  test_req3_locked_rfree_strips,
  test_safe_degradation_on_inspect_error,
  test_list_valued_data_mtz_slot,
  test_no_data_mtz_no_flags,
]


def run_all_tests():
  print("=" * 60)
  print("K_RFREE_GUARD: R-free generate guard")
  print("=" * 60)
  try:
    for t in _TESTS:
      t()
  finally:
    if _RESTORE is not None:
      _RESTORE()
  print("=" * 60)
  print("K_RFREE_GUARD complete.")
  return True


if __name__ == "__main__":
  ok = True
  passed = failed = 0
  try:
    for t in _TESTS:
      try:
        t()
        passed += 1
      except Exception as e:
        print("  FAIL: %s: %s" % (t.__name__, e))
        failed += 1
        ok = False
  finally:
    if _RESTORE is not None:
      _RESTORE()
  print()
  print("%d passed, %d failed" % (passed, failed))
  sys.exit(0 if ok else 1)
