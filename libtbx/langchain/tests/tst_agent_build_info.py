"""K_H2: Agent build info regression tests.

v119.H2.  Validates the server build metadata channel:
  - VERSION file at langchain/ root
  - core/_version.get_version()
  - core/llm.compute_defaults_fingerprint()
  - core/_build_info.get_agent_build_info()
  - core/_build_info.inject_agent_build()
  - knowledge/api_schema.py schema entry
  - backwards-compat: no existing field changed
"""
from __future__ import absolute_import, division, print_function

import os
import re


# ---------- Import / file-discovery helpers ------------------------

def _import_core_llm():
  try:
    from libtbx.langchain.core import llm
  except ImportError:
    from core import llm
  return llm


def _import_api_schema():
  try:
    from libtbx.langchain.knowledge import api_schema
  except ImportError:
    from knowledge import api_schema
  return api_schema


def _import_build_info():
  try:
    from libtbx.langchain.core import _build_info
  except ImportError:
    from core import _build_info
  return _build_info


def _import_version_module():
  try:
    from libtbx.langchain.core import _version
  except ImportError:
    from core import _version
  return _version


def _source_path(rel_path):
  """Return absolute path to a file under the langchain/ tree."""
  here = os.path.dirname(os.path.abspath(__file__))
  candidates = [
    os.path.join(here, "..", rel_path),
    os.path.join(here, "..", "..", rel_path),
  ]
  for c in candidates:
    c = os.path.normpath(c)
    if os.path.isfile(c):
      return c
  return None


# =====================================================================
# §9.1: Schema tests (4)
# =====================================================================

def test_response_schema_has_agent_build():
  """Schema entry exists in RESPONSE_V2_SCHEMA."""
  api_schema = _import_api_schema()
  assert "agent_build" in api_schema.RESPONSE_V2_SCHEMA, (
    "RESPONSE_V2_SCHEMA missing 'agent_build' entry")
  print("  PASS: test_response_schema_has_agent_build")


def test_agent_build_has_expected_subfields():
  """All three subfields present with correct types and defaults."""
  api_schema = _import_api_schema()
  spec = api_schema.RESPONSE_V2_SCHEMA["agent_build"]
  assert spec["type"] == dict, "agent_build type should be dict"
  assert spec["default"] == {}, "agent_build default should be {}"
  subfields = spec.get("subfields", {})
  expected = {"version", "defaults_fingerprint", "started_at"}
  assert set(subfields.keys()) == expected, (
    "subfields mismatch: got %s, expected %s" % (
      sorted(subfields.keys()), sorted(expected)))
  # All three subfields are str-typed
  for name in expected:
    assert subfields[name]["type"] == str, (
      "subfield %s should be type str" % name)
  # Defaults match the fallback dict
  assert subfields["version"]["default"] == "unknown"
  assert subfields["defaults_fingerprint"]["default"] == ""
  assert subfields["started_at"]["default"] == ""
  print("  PASS: test_agent_build_has_expected_subfields")


def test_api_version_still_2_0():
  """No protocol-envelope bump.  H2 is response-additive only."""
  api_schema = _import_api_schema()
  assert api_schema.API_VERSION == "2.0", (
    "API_VERSION drifted from 2.0 to %r -- H2 should not bump it"
    % api_schema.API_VERSION)
  print("  PASS: test_api_version_still_2_0")


def test_response_v2_required_unchanged():
  """agent_build is NOT a required field (backwards-compat)."""
  api_schema = _import_api_schema()
  assert "agent_build" not in api_schema.RESPONSE_V2_REQUIRED, (
    "agent_build must not be in RESPONSE_V2_REQUIRED -- old "
    "responses without it must still validate")
  # Existing required fields unchanged
  assert set(api_schema.RESPONSE_V2_REQUIRED) == {
    "api_version", "decision"}, (
    "RESPONSE_V2_REQUIRED changed: %s" % api_schema.RESPONSE_V2_REQUIRED)
  print("  PASS: test_response_v2_required_unchanged")


# =====================================================================
# §9.2: Constructor backwards-compat (3)
# =====================================================================

def test_create_response_signature_unchanged():
  """create_response signature has the same parameters as pre-H2.

  H2 adds no new parameter to any constructor.  The agent_build
  field is added later by the post-processor.
  """
  import inspect
  api_schema = _import_api_schema()
  sig = inspect.signature(api_schema.create_response)
  params = set(sig.parameters.keys())
  # Pre-H2 parameter list (from supplied api_schema.py).
  expected = {
    "program", "command", "reasoning", "strategy",
    "stop", "stop_reason",
    "experiment_type", "workflow_state",
    "warnings", "red_flags",
    "debug_log", "error",
    "server_version", "events",
  }
  assert params == expected, (
    "create_response signature changed.\n"
    "  got:      %s\n  expected: %s\n"
    "  added:    %s\n  removed:  %s" % (
      sorted(params), sorted(expected),
      sorted(params - expected), sorted(expected - params)))
  print("  PASS: test_create_response_signature_unchanged")


def test_create_response_sets_agent_build_to_empty_dict():
  """A constructor-built response has agent_build = {} (not populated).

  This is the apply_response_defaults branch-A behavior: when
  agent_build is missing from the constructor's dict (which it is,
  since the constructor doesn't pass it), apply_response_defaults
  sets it to the schema default ({}).  Subfield defaults are NOT
  applied -- that requires the field to be present-as-partial-dict.

  The post-processor inject_agent_build() then overwrites this {}
  with a populated dict before serialization.

  This test documents the exact behavior so future maintainers don't
  trip on it.
  """
  api_schema = _import_api_schema()
  resp = api_schema.create_response(
    program="phenix.refine",
    command="phenix.refine foo.pdb")
  assert "agent_build" in resp, (
    "create_response should populate agent_build (to {}) via "
    "apply_response_defaults")
  assert resp["agent_build"] == {}, (
    "create_response should set agent_build = {} (branch A); "
    "got %r" % resp["agent_build"])
  print("  PASS: test_create_response_sets_agent_build_to_empty_dict")


def test_apply_response_defaults_branch_A_for_missing_agent_build():
  """apply_response_defaults behavior gotcha: branch A vs branch B.

  The existing logic has two MUTUALLY EXCLUSIVE branches:
    if field not in response:
      response[field] = spec["default"]            # branch A
    elif spec["type"] == dict and "subfields" in spec:
      # ... fill missing subfields                  # branch B

  When agent_build is missing, branch A fires (sets to {}).
  Branch B does NOT fire -- subfield defaults are not applied.

  Subfield defaults only ever fire when agent_build is
  present-as-partial-dict.  In the normal constructor flow,
  agent_build is missing -> branch A -> {}; then the post-processor
  overwrites.

  This test exists so the gotcha doesn't bite a future maintainer
  who expects subfield defaults to apply in all cases.
  """
  api_schema = _import_api_schema()
  # Case 1: agent_build missing -> branch A -> {}, no subfields
  resp = {"api_version": "2.0", "decision": {}}
  result = api_schema.apply_response_defaults(resp)
  assert result["agent_build"] == {}, (
    "Branch A: missing field gets default {}, not subfield-populated")
  # Case 2: agent_build present-as-partial -> branch B -> subfields filled
  resp2 = {
    "api_version": "2.0",
    "decision": {},
    "agent_build": {"version": "test"},  # missing other 2 subfields
  }
  result2 = api_schema.apply_response_defaults(resp2)
  assert result2["agent_build"]["version"] == "test", (
    "Branch B: existing subfield preserved")
  assert result2["agent_build"]["defaults_fingerprint"] == "", (
    "Branch B: missing subfield gets default")
  assert result2["agent_build"]["started_at"] == "", (
    "Branch B: missing subfield gets default")
  print("  PASS: test_apply_response_defaults_branch_A_for_missing_agent_build")


# =====================================================================
# §9.3: Version helpers (4)
# =====================================================================

def test_version_file_exists_at_langchain_root():
  """A VERSION file exists at the langchain/ root."""
  here = os.path.dirname(os.path.abspath(__file__))
  candidates = [
    os.path.join(here, "..", "VERSION"),
    os.path.join(here, "..", "..", "VERSION"),
  ]
  found = any(os.path.isfile(os.path.normpath(c)) for c in candidates)
  assert found, (
    "VERSION file not found.  Expected at langchain/ root.  "
    "Searched: %s" % [os.path.normpath(c) for c in candidates])
  print("  PASS: test_version_file_exists_at_langchain_root")


def test_get_version_returns_string():
  """get_version() returns a non-empty string."""
  v = _import_version_module()
  version = v.get_version()
  assert isinstance(version, str), (
    "get_version should return str, got %r" % type(version))
  assert version, "get_version returned empty string"
  print("  PASS: test_get_version_returns_string")


def test_get_version_returns_unknown_when_file_missing():
  """Missing VERSION file -> 'unknown'."""
  v = _import_version_module()
  original = v._version_file_path
  try:
    v._version_file_path = lambda: "/nonexistent/path/VERSION"
    assert v.get_version() == "unknown"
  finally:
    v._version_file_path = original
  print("  PASS: test_get_version_returns_unknown_when_file_missing")


def test_get_version_returns_unknown_when_file_empty():
  """Empty (or whitespace-only) VERSION file -> 'unknown'."""
  import tempfile
  v = _import_version_module()
  original = v._version_file_path
  empty_path = None
  try:
    fd, empty_path = tempfile.mkstemp(prefix="empty_version_")
    os.close(fd)
    # File exists but empty
    v._version_file_path = lambda p=empty_path: p
    assert v.get_version() == "unknown", (
      "Empty VERSION should return 'unknown'")
    # Whitespace-only
    with open(empty_path, "w") as f:
      f.write("   \n  \t\n")
    assert v.get_version() == "unknown", (
      "Whitespace-only VERSION should return 'unknown'")
  finally:
    v._version_file_path = original
    if empty_path and os.path.exists(empty_path):
      os.unlink(empty_path)
  print("  PASS: test_get_version_returns_unknown_when_file_empty")


def test_get_version_returns_unknown_when_file_undecodable():
  """Binary garbage in VERSION must not raise (never-raises contract).

  Defensive: VERSION is a hand-written ASCII file but a deploy
  process bug or filesystem corruption could leave non-UTF-8 bytes.
  get_version's contract says 'never raises' -- this test verifies
  that contract holds for binary garbage.

  On UTF-8 platforms (Linux/Mac default) the broad except catches
  UnicodeDecodeError and we get "unknown".  On platforms with
  permissive default encodings the bytes may decode-as-garbage
  without raising -- that's also contract-compliant (returns a
  string, didn't raise).  The assertion checks the contract, not
  the specific value.
  """
  import tempfile
  v = _import_version_module()
  original = v._version_file_path
  garbage_path = None
  try:
    fd, garbage_path = tempfile.mkstemp(prefix="garbage_version_")
    os.write(fd, b"\xff\xfe\x00\x01\xbb\xcc some binary garbage")
    os.close(fd)
    v._version_file_path = lambda p=garbage_path: p
    # Must not raise.  Must return a string.
    # (Specific value depends on platform default encoding.)
    result = v.get_version()
    assert isinstance(result, str), (
      "get_version should return str even on binary input, "
      "got %r (type %s)" % (result, type(result).__name__))
  finally:
    v._version_file_path = original
    if garbage_path and os.path.exists(garbage_path):
      os.unlink(garbage_path)
  print("  PASS: test_get_version_returns_unknown_when_file_undecodable")


# =====================================================================
# §9.4: Fingerprint tests (3)
# =====================================================================

def test_fingerprint_format():
  """Fingerprint matches 'sha256:' + 64 hex chars."""
  llm = _import_core_llm()
  fp = llm.compute_defaults_fingerprint()
  assert isinstance(fp, str), "fingerprint must be str"
  assert fp.startswith("sha256:"), (
    "fingerprint must start with 'sha256:', got %r" % fp[:20])
  hex_part = fp[len("sha256:"):]
  assert len(hex_part) == 64, (
    "fingerprint hex part must be 64 chars, got %d" % len(hex_part))
  assert re.match(r"^[0-9a-f]{64}$", hex_part), (
    "fingerprint hex part must be all lowercase hex, got %r" % hex_part)
  print("  PASS: test_fingerprint_format")


def test_fingerprint_deterministic():
  """Two consecutive calls produce the same fingerprint."""
  llm = _import_core_llm()
  fp1 = llm.compute_defaults_fingerprint()
  fp2 = llm.compute_defaults_fingerprint()
  fp3 = llm.compute_defaults_fingerprint()
  assert fp1 == fp2 == fp3, (
    "fingerprint not deterministic: %r %r %r" % (fp1, fp2, fp3))
  print("  PASS: test_fingerprint_deterministic")


def test_fingerprint_excludes_retired_models():
  """Mutating RETIRED_MODELS does not change the fingerprint.

  Gemini Q5 critique: retirement updates are tombstones; they
  must not perturb the active-build fingerprint signal.  Without
  this property, H3 startup canary expectations would need
  updating on every retirement.
  """
  llm = _import_core_llm()
  original = llm.RETIRED_MODELS
  fp_before = llm.compute_defaults_fingerprint()
  try:
    # Synthesize a different RETIRED_MODELS
    llm.RETIRED_MODELS = frozenset([
      "gemini-2.0-flash",
      "synthetic-retired-test-model-1",
      "synthetic-retired-test-model-2",
    ])
    fp_after = llm.compute_defaults_fingerprint()
    assert fp_before == fp_after, (
      "Fingerprint changed when RETIRED_MODELS was mutated.\n"
      "  before: %s\n  after:  %s\n"
      "RETIRED_MODELS must not be part of the fingerprint payload."
      % (fp_before, fp_after))
  finally:
    llm.RETIRED_MODELS = original
  print("  PASS: test_fingerprint_excludes_retired_models")


# =====================================================================
# §9.5: build_info tests (3)
# =====================================================================

def test_get_agent_build_info_returns_three_keys():
  """get_agent_build_info returns exactly three keys."""
  bi = _import_build_info()
  info = bi.get_agent_build_info()
  expected = {"version", "defaults_fingerprint", "started_at"}
  assert set(info.keys()) == expected, (
    "get_agent_build_info keys mismatch.\n"
    "  got:      %s\n  expected: %s" % (
      sorted(info.keys()), sorted(expected)))
  print("  PASS: test_get_agent_build_info_returns_three_keys")


def test_get_agent_build_info_subfields_are_strings():
  """All three subfields are str-valued."""
  bi = _import_build_info()
  info = bi.get_agent_build_info()
  for key in ("version", "defaults_fingerprint", "started_at"):
    assert isinstance(info[key], str), (
      "subfield %s must be str, got %r" % (key, type(info[key])))
  print("  PASS: test_get_agent_build_info_subfields_are_strings")


def test_started_at_format():
  """started_at is strict UTC ISO 8601: YYYY-MM-DDTHH:MM:SSZ.

  No microseconds, no offset.  Pinning the format here keeps it
  stable for log aggregators and uptime calculations.
  """
  bi = _import_build_info()
  info = bi.get_agent_build_info()
  started_at = info["started_at"]
  pattern = r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$"
  assert re.match(pattern, started_at), (
    "started_at does not match strict UTC ISO 8601 "
    "(YYYY-MM-DDTHH:MM:SSZ): got %r" % started_at)
  print("  PASS: test_started_at_format")


# =====================================================================
# §9.6: Injection helper tests (3) -- REV 3
# =====================================================================

def test_inject_agent_build_adds_field_to_empty_response():
  """inject_agent_build({}) adds agent_build with all three subfields."""
  bi = _import_build_info()
  response = {}
  bi.inject_agent_build(response)
  assert "agent_build" in response, (
    "inject_agent_build should add 'agent_build' key")
  info = response["agent_build"]
  assert set(info.keys()) == {
    "version", "defaults_fingerprint", "started_at"}, (
    "agent_build subfields mismatch: %s" % sorted(info.keys()))
  print("  PASS: test_inject_agent_build_adds_field_to_empty_response")


def test_inject_agent_build_overwrites_existing_agent_build():
  """inject_agent_build overwrites a pre-existing agent_build = {}.

  This confirms the apply_response_defaults {} doesn't leak through.
  After the constructor sets agent_build = {} (branch A behavior),
  the post-processor must replace it with a populated dict.
  """
  bi = _import_build_info()
  # Simulate post-constructor state: agent_build is empty dict
  response = {"agent_build": {}}
  bi.inject_agent_build(response)
  assert response["agent_build"] != {}, (
    "inject_agent_build must overwrite empty {}, not preserve it")
  assert "version" in response["agent_build"]
  assert "defaults_fingerprint" in response["agent_build"]
  assert "started_at" in response["agent_build"]
  # Also: a populated dict gets fully replaced (not merged)
  response2 = {"agent_build": {"version": "spoofed", "extra": "junk"}}
  bi.inject_agent_build(response2)
  assert response2["agent_build"]["version"] != "spoofed", (
    "inject should replace, not merge")
  assert "extra" not in response2["agent_build"], (
    "inject should replace, not merge")
  print("  PASS: test_inject_agent_build_overwrites_existing_agent_build")


def test_inject_agent_build_preserves_other_fields():
  """inject_agent_build only touches agent_build; other fields intact."""
  bi = _import_build_info()
  response = {
    "api_version": "2.0",
    "server_version": "phenix-1.21.0",
    "decision": {"program": "phenix.refine", "command": "x"},
    "stop": False,
    "stop_reason": None,
    "metadata": {"warnings": [], "red_flags": []},
    "debug": {"log": ["msg1", "msg2"], "timing_ms": 100},
    "error": None,
    "events": [{"type": "x"}],
  }
  # Snapshot of all keys except agent_build (which doesn't exist yet)
  import copy
  before = copy.deepcopy(response)
  bi.inject_agent_build(response)
  for key in before:
    assert response[key] == before[key], (
      "inject_agent_build clobbered key %r: %r != %r"
      % (key, response[key], before[key]))
  print("  PASS: test_inject_agent_build_preserves_other_fields")


def test_inject_agent_build_never_raises_on_helper_failure():
  """inject_agent_build's never-raises contract survives a broken helper.

  Defense-in-depth: the run_ai_agent.py call site only catches
  ImportError around the import-and-call combo, trusting
  inject_agent_build's "Never raises" docstring promise.  This test
  locks in the contract by monkey-patching get_agent_build_info to
  raise, then asserting inject_agent_build catches it and falls
  back to _FALLBACK_BUILD_INFO.

  Without this test, a future refactor could silently break the
  contract; the next runtime exception would propagate up through
  _build_group_args_response and fail the entire response build.
  """
  bi = _import_build_info()
  original = bi.get_agent_build_info
  try:
    def boom():
      raise RuntimeError("synthetic failure for testing")
    bi.get_agent_build_info = boom
    response = {}
    # Must not raise.
    bi.inject_agent_build(response)
    # Must still produce a valid agent_build (the fallback dict).
    assert "agent_build" in response, (
      "inject_agent_build should fall back to _FALLBACK_BUILD_INFO "
      "when get_agent_build_info raises")
    ab = response["agent_build"]
    assert set(ab.keys()) == {
      "version", "defaults_fingerprint", "started_at"}, (
      "fallback dict should have all three subfields")
    # Fallback values are syntactically valid strings
    assert isinstance(ab["version"], str)
    assert isinstance(ab["defaults_fingerprint"], str)
    assert isinstance(ab["started_at"], str)
  finally:
    bi.get_agent_build_info = original
  print("  PASS: test_inject_agent_build_never_raises_on_helper_failure")


# =====================================================================
# §9.7: server_version compatibility tests (2)
# =====================================================================

def test_server_version_still_string_type():
  """Schema declares server_version as type=str, NOT converted to dict.

  Backwards-compat: client code (any future code) that reads
  server_version expects a string.  H2 must not change this.
  """
  api_schema = _import_api_schema()
  spec = api_schema.RESPONSE_V2_SCHEMA["server_version"]
  assert spec["type"] == str, (
    "server_version type changed from str to %r -- backwards-compat "
    "violation" % spec["type"])
  print("  PASS: test_server_version_still_string_type")


def test_server_version_separate_from_agent_build():
  """server_version and agent_build are independent fields."""
  api_schema = _import_api_schema()
  resp = api_schema.create_response(
    program="phenix.refine",
    command="x",
    server_version="phenix-1.21.0-test")
  assert resp["server_version"] == "phenix-1.21.0-test", (
    "server_version not preserved through create_response")
  # And agent_build is independent (initialized to {} by branch A)
  assert "agent_build" in resp
  assert resp["server_version"] != resp["agent_build"], (
    "fields should be independent values")
  print("  PASS: test_server_version_separate_from_agent_build")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
  # §9.1: Schema (4)
  test_response_schema_has_agent_build()
  test_agent_build_has_expected_subfields()
  test_api_version_still_2_0()
  test_response_v2_required_unchanged()
  # §9.2: Constructor backwards-compat (3)
  test_create_response_signature_unchanged()
  test_create_response_sets_agent_build_to_empty_dict()
  test_apply_response_defaults_branch_A_for_missing_agent_build()
  # §9.3: Version helpers (5)
  test_version_file_exists_at_langchain_root()
  test_get_version_returns_string()
  test_get_version_returns_unknown_when_file_missing()
  test_get_version_returns_unknown_when_file_empty()
  test_get_version_returns_unknown_when_file_undecodable()
  # §9.4: Fingerprint (3)
  test_fingerprint_format()
  test_fingerprint_deterministic()
  test_fingerprint_excludes_retired_models()
  # §9.5: build_info (3)
  test_get_agent_build_info_returns_three_keys()
  test_get_agent_build_info_subfields_are_strings()
  test_started_at_format()
  # §9.6: Injection helper (4)
  test_inject_agent_build_adds_field_to_empty_response()
  test_inject_agent_build_overwrites_existing_agent_build()
  test_inject_agent_build_preserves_other_fields()
  test_inject_agent_build_never_raises_on_helper_failure()
  # §9.7: server_version compat (2)
  test_server_version_still_string_type()
  test_server_version_separate_from_agent_build()


if __name__ == "__main__":
  run_all_tests()
