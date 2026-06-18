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

v120 (server parity): the original fix inspected the MTZ with inspect_mtz,
which reads the file -- a no-op on the SERVER, where client paths don't
exist (files_local=False).  The fix is now CLIENT-EXTRACTED and tri-state:
the client ships `input_mtz_has_rfree` (True/False/None) in session_info, and
`CommandBuilder._input_mtz_rfree_state(files, context)` returns:
  True  -> flags exist  (strip/suppress generate)
  False -> flags absent (generate -- genuine first refinement)
  None  -> undetermined (strip/suppress -- conservative: a recoverable
           rfree_flags_missing error beats silently overwriting a test set)
Sources, in order: (0) context.input_mtz_has_rfree (client fact, the only one
that works server-side); (1) context.mtz_inspection positive only;
(2) inspect_mtz on the selected data MTZ, LOCAL runs only (gated on
files_local so a server read can't masquerade as a confident False).

These tests exercise the REAL `_input_mtz_rfree_state` (bound from the real
CommandBuilder class source, with inspect_mtz stubbed) plus the three-valued
generate decision the invariant implements.
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import types

_HERE = os.path.dirname(os.path.abspath(__file__))
_MISSING = object()  # sentinel for "key absent" in sys.modules snapshot


def _load_helper():
  """Bind CommandBuilder._input_mtz_rfree_state without importing the heavy
  ProgramRegistry/yaml stack.  We exec the method's source against a stub
  module that provides a controllable inspect_mtz."""
  cb_path = os.environ.get(
    "COMMAND_BUILDER_PY",
    os.path.join(_HERE, "..", "agent", "command_builder.py"))
  if not os.path.exists(cb_path):
    return None, None, None
  src = open(cb_path).read()
  if "_input_mtz_rfree_state" not in src:
    return None, None, None
  # Extract the method body and dedent to a module-level function.
  lines = src.splitlines()
  start = None
  for i, ln in enumerate(lines):
    if ln.strip().startswith("def _input_mtz_rfree_state"):
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
  #   files_with_flags : paths whose inspection reports an R-free label
  #   raise            : inspect_mtz raises (cctbx unavailable / crash)
  #   returns_none     : paths for which inspect_mtz returns None
  #                      (unreadable / not an MTZ) -> tri-state "undetermined"
  stub_flags = {"files_with_flags": set(), "raise": False,
                "returns_none": set()}

  def _inspect(path):
    if stub_flags["raise"]:
      raise RuntimeError("cctbx unavailable")
    if path in stub_flags["returns_none"]:
      return None
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
  return ns["_input_mtz_rfree_state"], stub_flags, _restore


class _Ctx(object):
  def __init__(self, mtz_inspection=None, input_mtz_has_rfree=None,
               files_local=True):
    self.mtz_inspection = mtz_inspection
    # v120: client-extracted tri-state fact (authoritative source 0).
    self.input_mtz_has_rfree = input_mtz_has_rfree
    # v120: when False (server mode), the helper must NOT read the file
    # (it cannot) and must return None rather than a wrong False.
    self.files_local = files_local


# Bind once.
_HELPER, _STUB, _RESTORE = _load_helper()


def _state(files, ctx):
  # The method takes (self, files, context); self is unused, pass None.
  # Returns the TRI-STATE: True | False | None.
  return _HELPER(None, files, ctx)


def _decide(rfree_mtz, files, ctx, strategy=None):
  """Mirror of the invariant's THREE-valued generate decision using the real
  _input_mtz_rfree_state.  Generates ONLY when the state is positively False;
  strips/suppresses when rfree_mtz is locked, state is True, or state is None
  (undetermined -> conservative)."""
  strategy = dict(strategy or {})
  if rfree_mtz:
    rfree_state = True
  else:
    rfree_state = _state(files, ctx)
  if rfree_state is False:
    if "generate_rfree_flags" not in strategy:
      strategy["generate_rfree_flags"] = True
  else:
    if strategy.get("generate_rfree_flags"):
      del strategy["generate_rfree_flags"]
  return strategy.get("generate_rfree_flags", False)


def test_req1_no_flags_generates():
  if _HELPER is None:
    print("  SKIP: helper not found (fix not applied?)")
    return
  _STUB["files_with_flags"] = set()
  _STUB["raise"] = False
  _STUB["returns_none"] = set()
  # Local read confirms no flags -> state False -> generate.
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
  _STUB["returns_none"] = set()
  _STUB["raise"] = True  # inspect_mtz blows up
  # v120 (Gemini review): detection failure -> state None (undetermined).  We do
  # NOT generate when undetermined — a recoverable rfree_flags_missing error is
  # far preferable to silently overwriting an existing test set.  So an
  # LLM-forced generate is STRIPPED, and none is added.
  g = _decide(None, {"data_mtz": "unknown.mtz"}, _Ctx(),
              strategy={"generate_rfree_flags": True})
  assert g is False, "undetermined must NOT generate (strip), got %r" % g
  st = _state({"data_mtz": "unknown.mtz"}, _Ctx())
  assert st is None, "inspect failure must yield None, got %r" % st
  print("  PASS: test_safe_degradation_on_inspect_error")


def test_list_valued_data_mtz_slot():
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = {"multi.mtz"}
  _STUB["raise"] = False
  _STUB["returns_none"] = set()
  # multiple:true slots can be lists; helper takes the first
  st = _state({"data_mtz": ["multi.mtz", "other.mtz"]}, _Ctx())
  assert st is True, "list-valued data_mtz first element must be inspected"
  print("  PASS: test_list_valued_data_mtz_slot")


def test_no_data_mtz_no_flags():
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = set()
  _STUB["raise"] = False
  _STUB["returns_none"] = set()
  # No data MTZ selected locally -> cannot determine -> None (undetermined),
  # NOT a confident False (which would wrongly invite generation).
  st = _state({}, _Ctx())
  assert st is None, "no data_mtz selected -> undetermined (None), got %r" % st
  print("  PASS: test_no_data_mtz_no_flags")


def test_client_fact_true_strips():
  # Source 0: client-shipped fact True is authoritative -> strip.  This is the
  # AIAgent_37 server scenario: server can't read the file but the client
  # shipped input_mtz_has_rfree=True.
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = set()
  _STUB["raise"] = False
  _STUB["returns_none"] = set()
  ctx = _Ctx(input_mtz_has_rfree=True, files_local=False)
  st = _state({"data_mtz": "/client/beta_blip_001.mtz"}, ctx)
  assert st is True, "client fact True must win, got %r" % st
  g = _decide(None, {"data_mtz": "/client/beta_blip_001.mtz"}, ctx,
              strategy={"generate_rfree_flags": True})
  assert g is False, "client fact True must strip generate, got %r" % g
  print("  PASS: test_client_fact_true_strips")


def test_client_fact_false_generates():
  # Source 0: client confirmed NO flags -> generate (genuine first refinement),
  # even server-side where the file can't be read.
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = set()
  _STUB["raise"] = False
  _STUB["returns_none"] = set()
  ctx = _Ctx(input_mtz_has_rfree=False, files_local=False)
  st = _state({"data_mtz": "/client/data.mtz"}, ctx)
  assert st is False, "client fact False must be respected, got %r" % st
  g = _decide(None, {"data_mtz": "/client/data.mtz"}, ctx)
  assert g is True, "client fact False must generate, got %r" % g
  print("  PASS: test_client_fact_false_generates")


def test_server_no_fact_is_undetermined_and_strips():
  # KEY (Gemini pt 2): server mode (files_local=False) with NO client fact and
  # no precomputed inspection -> the helper must NOT read the (absent) file and
  # must return None.  The decision must then STRIP, not generate — a recoverable
  # error beats silently overwriting an existing test set.
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = {"/client/beta_blip_001.mtz"}  # would be True IF read
  _STUB["raise"] = False
  _STUB["returns_none"] = set()
  ctx = _Ctx(input_mtz_has_rfree=None, files_local=False)
  st = _state({"data_mtz": "/client/beta_blip_001.mtz"}, ctx)
  assert st is None, "server + no fact must be undetermined (None), got %r" % st
  g = _decide(None, {"data_mtz": "/client/beta_blip_001.mtz"}, ctx,
              strategy={"generate_rfree_flags": True})
  assert g is False, "undetermined on server must strip generate, got %r" % g
  print("  PASS: test_server_no_fact_is_undetermined_and_strips")


def test_client_fact_overrides_local_read():
  # When both a client fact and a (contradictory) local read are available, the
  # client fact (source 0) wins — it is the authoritative input-MTZ answer.
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = set()  # local read would say False
  _STUB["raise"] = False
  _STUB["returns_none"] = set()
  ctx = _Ctx(input_mtz_has_rfree=True, files_local=True)
  st = _state({"data_mtz": "in.mtz"}, ctx)
  assert st is True, "client fact must override local read, got %r" % st
  print("  PASS: test_client_fact_overrides_local_read")


def test_unreadable_mtz_is_undetermined():
  # Local mode, file present but inspect_mtz returns None (unreadable / not an
  # MTZ) -> undetermined (None), NOT a confident False.
  if _HELPER is None:
    print("  SKIP")
    return
  _STUB["files_with_flags"] = set()
  _STUB["raise"] = False
  _STUB["returns_none"] = {"corrupt.mtz"}
  st = _state({"data_mtz": "corrupt.mtz"}, _Ctx(files_local=True))
  assert st is None, "unreadable MTZ must be undetermined (None), got %r" % st
  # And the decision must not force-generate on an undetermined input.
  g = _decide(None, {"data_mtz": "corrupt.mtz"}, _Ctx(files_local=True),
              strategy={"generate_rfree_flags": True})
  assert g is False, "undetermined must strip, got %r" % g
  print("  PASS: test_unreadable_mtz_is_undetermined")


_TESTS = [
  test_req1_no_flags_generates,
  test_req2_input_flags_kept_direct_inspect,
  test_req2_strips_llm_forced_generate,
  test_req2_precomputed_inspection,
  test_req3_locked_rfree_strips,
  test_safe_degradation_on_inspect_error,
  test_list_valued_data_mtz_slot,
  test_no_data_mtz_no_flags,
  test_client_fact_true_strips,
  test_client_fact_false_generates,
  test_server_no_fact_is_undetermined_and_strips,
  test_client_fact_overrides_local_read,
  test_unreadable_mtz_is_undetermined,
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
