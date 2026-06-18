"""K_RFREE_PLUMBING: R-free session_info fields survive the client->server wire.

Covers both the scalar input_mtz_has_rfree and the v120.2 per-file mtz_rfree_map.

The R-free generate decision (command_builder._input_mtz_rfree_state) relies on
the CLIENT-extracted fact `input_mtz_has_rfree` reaching the server-side builder.
That fact is tri-state (True / False / None) and travels:

    ai_agent.py (client)              session_info["input_mtz_has_rfree"]
      -> api_client.build_session_state            (hop 1)
      -> api_client.build_request_v2 whitelist     (hop 2; the wire chokepoint,
         via knowledge.api_schema.create_request + apply_request_defaults)
      -> JSON transport
      -> run_ai_agent.py map-back                  (hop 3)
      -> command_builder.from_state -> CommandContext.input_mtz_has_rfree

The CRITICAL invariant these tests pin: every hop uses an explicit
`is not None` guard (NOT a truthy guard).  A truthy guard would drop a
confirmed-`False` (treating it like absent/None), which would flip the server's
decision from "generate -- genuine first refinement" to "do not generate".  So
this test asserts:
    True   -> arrives True
    False  -> arrives False   (must NOT be dropped)
    None   -> dropped on the wire -> arrives absent (server reads None ->
              conservative "do not generate")

It exercises the REAL build_session_state / build_request_v2 / create_request /
apply_request_defaults (with only the unrelated transport + core.llm sanitizer
helpers stubbed), then replicates run_ai_agent's map-back exactly.
"""
from __future__ import absolute_import, division, print_function

import json
import os
import sys
import types

_HERE = os.path.dirname(os.path.abspath(__file__))
_MISSING = object()


def _candidate(*parts):
  return os.path.join(_HERE, "..", *parts)


def _load():
  """Import the real api_client + api_schema, stubbing only the unrelated
  transport / core.llm sanitizer helpers they import.  Returns
  (build_session_state, build_request_v2, restore) or (None, None, None)."""
  api_client_py = os.environ.get(
    "API_CLIENT_PY", _candidate("agent", "api_client.py"))
  api_schema_py = os.environ.get(
    "API_SCHEMA_PY", _candidate("knowledge", "api_schema.py"))
  if not (os.path.exists(api_client_py) and os.path.exists(api_schema_py)):
    return None, None, None

  saved = {}
  created = []

  def _ensure_pkg(name):
    if name not in sys.modules:
      saved[name] = sys.modules.get(name, _MISSING)
      m = types.ModuleType(name)
      m.__path__ = []
      sys.modules[name] = m
      created.append(name)
    return sys.modules[name]

  def _install(name, module):
    saved[name] = sys.modules.get(name, _MISSING)
    sys.modules[name] = module
    created.append(name)

  for pkg in ("libtbx", "libtbx.langchain", "libtbx.langchain.agent",
              "libtbx.langchain.knowledge", "libtbx.langchain.core"):
    _ensure_pkg(pkg)

  # Stub the transport sanitizers (identity) + core.llm helpers.  These are
  # unrelated to the plumbing under test; stubbing avoids importing the heavy
  # provider stack.
  transport = types.ModuleType("libtbx.langchain.agent.transport")
  for fn in ("encode_request", "decode_request", "sanitize_string",
             "sanitize_for_transport", "sanitize_dict_recursive",
             "sanitize_request", "sanitize_response", "get_transport_config"):
    setattr(transport, fn, (lambda *a, **k: (a[0] if a else None))
            if fn not in ("get_transport_config",) else (lambda *a, **k: {}))
  _install("libtbx.langchain.agent.transport", transport)

  core_llm = types.ModuleType("libtbx.langchain.core.llm")
  for fn in ("get_llm", "call_llm", "default_model_for_provider",
             "normalize_ollama_openai_base_url", "resolve_model_for_provider"):
    setattr(core_llm, fn, (lambda *a, **k: "model"))
  _install("libtbx.langchain.core.llm", core_llm)

  # Load api_schema then api_client from source under their package names.
  import importlib.util

  def _load_mod(modname, path):
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    saved[modname] = sys.modules.get(modname, _MISSING)
    sys.modules[modname] = mod
    created.append(modname)
    spec.loader.exec_module(mod)
    return mod

  def _restore():
    for name in reversed(created):
      v = saved.get(name, _MISSING)
      if v is _MISSING:
        sys.modules.pop(name, None)
      else:
        sys.modules[name] = v

  try:
    _load_mod("libtbx.langchain.knowledge.api_schema", api_schema_py)
    api_client = _load_mod(
      "libtbx.langchain.agent.api_client", api_client_py)
  except Exception as e:
    print("  (load error: %s: %s)" % (type(e).__name__, e))
    _restore()
    return None, None, None

  return (api_client.build_session_state,
          api_client.build_request_v2,
          _restore)


_BSS, _BRV2, _RESTORE = _load()


def _map_back(session_state):
  """Exact replica of run_ai_agent.py's map-back for the R-free fields."""
  session_info = {}
  if session_state.get("input_mtz_has_rfree") is not None:
    session_info["input_mtz_has_rfree"] = session_state["input_mtz_has_rfree"]
  if session_state.get("mtz_rfree_map") is not None:
    session_info["mtz_rfree_map"] = session_state["mtz_rfree_map"]
  return session_info


def _roundtrip_map(rfree_map, present=True):
  """Run a per-file map through all hops; return it as seen after the map-back,
  or '<ABSENT>' if it never arrived."""
  session_info = {"experiment_type": "xray", "best_files": {},
                  "client_protocol_version": 8}
  if present:
    session_info["mtz_rfree_map"] = rfree_map
  session_state = _BSS(session_info)                       # hop 1
  request = _BRV2(files={"data_mtz": "x.mtz"}, cycle_number=1,
                  session_state=session_state, client_version=8)  # hop 2
  request = json.loads(json.dumps(request))                # JSON transport
  nss = request.get("session_state", {})
  session_info_back = _map_back(nss)                       # hop 3
  return session_info_back.get("mtz_rfree_map", "<ABSENT>")


def _roundtrip(fact, present=True):
  """Run a fact value through all hops; return the value as seen at the end of
  the map-back, or the sentinel string '<ABSENT>' if it never arrived."""
  session_info = {"experiment_type": "xray", "best_files": {},
                  "client_protocol_version": 7}
  if present:
    session_info["input_mtz_has_rfree"] = fact

  session_state = _BSS(session_info)                       # hop 1
  request = _BRV2(files={"data_mtz": "x.mtz"}, cycle_number=1,
                  session_state=session_state, client_version=7)  # hop 2
  request = json.loads(json.dumps(request))                # JSON transport
  nss = request.get("session_state", {})
  session_info_back = _map_back(nss)                       # hop 3
  return session_info_back.get("input_mtz_has_rfree", "<ABSENT>")


def test_true_survives():
  if _BSS is None:
    print("  SKIP: api_client/api_schema not importable")
    return
  got = _roundtrip(True)
  assert got is True, "True must arrive as True, got %r" % got
  print("  PASS: test_true_survives")


def test_false_survives_not_dropped():
  # THE critical case: a truthy guard would drop False.
  if _BSS is None:
    print("  SKIP")
    return
  got = _roundtrip(False)
  assert got is False, "False must arrive as False (not dropped), got %r" % got
  print("  PASS: test_false_survives_not_dropped")


def test_none_is_absent():
  if _BSS is None:
    print("  SKIP")
    return
  got = _roundtrip(None)
  assert got == "<ABSENT>", "None must be dropped -> absent, got %r" % got
  print("  PASS: test_none_is_absent")


def test_missing_field_is_absent():
  if _BSS is None:
    print("  SKIP")
    return
  got = _roundtrip(None, present=False)
  assert got == "<ABSENT>", "absent field must stay absent, got %r" % got
  print("  PASS: test_missing_field_is_absent")


def test_rfree_map_survives():
  # v120.2: the per-file map (dict with mixed True/False) must arrive intact,
  # including the False entries (a truthy-collapsing bug would corrupt them).
  if _BSS is None:
    print("  SKIP")
    return
  m = {"beta_blip_001.mtz": True, "PHASER.1.mtz": False,
       "beta_blip_P3221.mtz": False}
  got = _roundtrip_map(m)
  assert got == m, "map must survive intact, got %r" % got
  print("  PASS: test_rfree_map_survives")


def test_rfree_map_absent_stays_absent():
  if _BSS is None:
    print("  SKIP")
    return
  got = _roundtrip_map(None, present=False)
  assert got == "<ABSENT>", "absent map must stay absent, got %r" % got
  print("  PASS: test_rfree_map_absent_stays_absent")


def test_rfree_map_empty_survives():
  # An empty map is distinguishable from absent (is-not-None guard).
  if _BSS is None:
    print("  SKIP")
    return
  got = _roundtrip_map({})
  assert got == {}, "empty map must survive as {}, got %r" % got
  print("  PASS: test_rfree_map_empty_survives")


_TESTS = [
  test_true_survives,
  test_false_survives_not_dropped,
  test_none_is_absent,
  test_missing_field_is_absent,
  test_rfree_map_survives,
  test_rfree_map_absent_stays_absent,
  test_rfree_map_empty_survives,
]


def run_all_tests():
  print("=" * 60)
  print("K_RFREE_PLUMBING: input_mtz_has_rfree wire round-trip")
  print("=" * 60)
  try:
    for t in _TESTS:
      t()
  finally:
    if _RESTORE is not None:
      _RESTORE()
  print("=" * 60)
  print("K_RFREE_PLUMBING complete.")
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
