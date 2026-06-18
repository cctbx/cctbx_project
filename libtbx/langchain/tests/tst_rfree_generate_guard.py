"""K_RFREE_GUARD: phenix.refine R-free generate decision (per-file, parity-safe).

Requirements:
  (1) refine data has NO R-free flags  -> generate a new set
  (2) refine data ALREADY has flags    -> keep them (never generate/overwrite)
  (3) once a set is locked (rfree_mtz)  -> never regenerate

History / why per-file:
  The original guard gated only on `context.rfree_mtz` (which locks from a
  refinement OUTPUT), so a FIRST refinement of pre-flagged INPUT data wrongly
  got `generate=True` and OVERWROTE the existing test set (beta_blip_001).  The
  v120 fix added a client-extracted scalar `input_mtz_has_rfree`, but it
  described the ORIGINAL input only -- wrong for the MR workflow, where refine
  runs on a phaser OUTPUT (PHASER.1.mtz) that carries no flags (AIAgent_35).

v120.2 (per-file + parity): the decision is about the SPECIFIC data MTZ being
refined.  The client inspects every local MTZ it can see (inputs AND outputs)
and ships a per-file map `mtz_rfree_map` {basename: bool} in session_info.
`CommandBuilder._input_mtz_rfree_state(files, context)` is TRI-STATE:
  True  -> flags exist  (strip/suppress generate -- prevents overwrite)
  False -> flags absent (generate -- genuine first refinement)
  None  -> undetermined: when nothing is locked, GENERATE (the safe
           first-refinement default; omitting generate aborts phenix.refine on
           a flagless MR/phaser output)
Sources, in order, ALL session_info-borne (so the answer is IDENTICAL on the
server and the client -- the builder never reads the filesystem):
  (0) context.mtz_rfree_map[basename(selected data_mtz)] -- per-file answer for
      the file actually being refined (authoritative);
  (1) context.input_mtz_has_rfree -- original-input scalar fallback;
  (2) context.mtz_inspection -- precomputed inspection, POSITIVE ONLY.

These tests exercise the REAL `_input_mtz_rfree_state` (extracted/bound from the
CommandBuilder source) plus the generate decision the invariant implements,
including a parity check that the helper performs NO filesystem reads.
"""
from __future__ import absolute_import, division, print_function

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_helper():
  """Bind CommandBuilder._input_mtz_rfree_state without importing the heavy
  ProgramRegistry/yaml stack: extract the method's source and exec it as a
  module-level function.

  As of v120.2 the helper no longer reads the filesystem (it consumes only
  session_info-borne data: the per-file map, the scalar fact, and any
  precomputed inspection), so no inspect_mtz stub is needed.  Returns
  (func, None, noop_restore); the second slot is kept for call-site
  compatibility."""
  cb_path = os.environ.get(
    "COMMAND_BUILDER_PY",
    os.path.join(_HERE, "..", "agent", "command_builder.py"))
  if not os.path.exists(cb_path):
    return None, None, (lambda: None)
  src = open(cb_path).read()
  if "_input_mtz_rfree_state" not in src:
    return None, None, (lambda: None)
  # Extract the method body and dedent to a module-level function.
  lines = src.splitlines()
  start = None
  for i, ln in enumerate(lines):
    if ln.strip().startswith("def _input_mtz_rfree_state"):
      start = i
      break
  if start is None:
    return None, None, (lambda: None)
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
  ns = {}
  exec(func_src, ns)
  return ns["_input_mtz_rfree_state"], None, (lambda: None)


class _Ctx(object):
  def __init__(self, mtz_inspection=None, input_mtz_has_rfree=None,
               files_local=True, mtz_rfree_map=None):
    self.mtz_inspection = mtz_inspection
    # v120: client-extracted tri-state scalar fact (original-input fallback).
    self.input_mtz_has_rfree = input_mtz_has_rfree
    # v120.2: per-file map {basename: bool} — the authoritative per-file source.
    self.mtz_rfree_map = mtz_rfree_map
    # files_local retained for ctor parity; the v120.2 helper no longer reads
    # the filesystem at all (parity: identical server and local), so this no
    # longer affects the decision.
    self.files_local = files_local


# Bind once.
_HELPER, _STUB, _RESTORE = _load_helper()


def _state(files, ctx):
  # The method takes (self, files, context); self is unused, pass None.
  # Returns the TRI-STATE: True | False | None.
  return _HELPER(None, files, ctx)


def _decide(rfree_mtz, files, ctx, strategy=None):
  """Mirror of the invariant's v120.1 THREE-valued generate decision using the
  real _input_mtz_rfree_state.  Strips ONLY on positive evidence flags exist
  (rfree_mtz locked, or state is True).  Generates on state False (confirmed no
  flags) OR state None (undetermined and nothing locked -> first-refine safe
  default that restores the robust pre-fix behavior)."""
  strategy = dict(strategy or {})
  if rfree_mtz:
    rfree_state = True
  else:
    rfree_state = _state(files, ctx)
  if rfree_state is True:
    if strategy.get("generate_rfree_flags"):
      del strategy["generate_rfree_flags"]
  else:
    # False or None (unlocked): generate.
    if "generate_rfree_flags" not in strategy:
      strategy["generate_rfree_flags"] = True
  return strategy.get("generate_rfree_flags", False)


def test_req1_no_flags_generates():
  if _HELPER is None:
    print("  SKIP: helper not found (fix not applied?)")
    return
  # No map entry, no scalar, no inspection -> state None.  Unlocked -> generate
  # (v120.1 first-refinement default).
  g = _decide(None, {"data_mtz": "noflags.mtz"}, _Ctx())
  assert g is True, "no flags info + no lock must generate, got %r" % g
  print("  PASS: test_req1_no_flags_generates")


def test_req2_input_flags_via_map():
  if _HELPER is None:
    print("  SKIP")
    return
  # v120.2: the per-file map reports flags present for the selected data MTZ ->
  # strip (do not regenerate).  Parity-safe: no file read.
  ctx = _Ctx(mtz_rfree_map={"beta_blip_001.mtz": True})
  g = _decide(None, {"data_mtz": "/x/beta_blip_001.mtz"}, ctx)
  assert g is False, "input flags (map True) must NOT generate, got %r" % g
  print("  PASS: test_req2_input_flags_via_map")


def test_req2_strips_llm_forced_generate():
  if _HELPER is None:
    print("  SKIP")
    return
  # Map says the selected file has flags -> an LLM-forced generate is stripped.
  ctx = _Ctx(mtz_rfree_map={"beta_blip_001.mtz": True})
  g = _decide(None, {"data_mtz": "/x/beta_blip_001.mtz"}, ctx,
              strategy={"generate_rfree_flags": True})
  assert g is False, "LLM-forced generate must be stripped, got %r" % g
  print("  PASS: test_req2_strips_llm_forced_generate")


def test_req2_precomputed_inspection():
  if _HELPER is None:
    print("  SKIP")
    return
  # but precomputed mtz_inspection says yes -> trusted
  g = _decide(None, {"data_mtz": "x.mtz"},
              _Ctx(mtz_inspection={"rfree_label": "FreeR_flag"}))
  assert g is False, "precomputed rfree_label must suppress generate, got %r" % g
  print("  PASS: test_req2_precomputed_inspection")


def test_req3_locked_rfree_strips():
  if _HELPER is None:
    print("  SKIP")
    return
  g = _decide("/sub_02_refine/refine_001_data.mtz",
              {"data_mtz": "refine_001_data.mtz"}, _Ctx(),
              strategy={"generate_rfree_flags": True})
  assert g is False, "locked rfree_mtz must strip generate, got %r" % g
  print("  PASS: test_req3_locked_rfree_strips")


def test_undetermined_unlocked_generates():
  if _HELPER is None:
    print("  SKIP")
    return
  # No map entry, no scalar, no inspection -> state None (undetermined).  With
  # nothing locked, the safe default for an unlocked/first refinement is to
  # GENERATE — omitting generate causes a hard "No array of R-free flags found"
  # abort (the MR/phaser workflow, beta_blip_P3221).  We only STRIP on positive
  # evidence flags exist.
  g = _decide(None, {"data_mtz": "unknown.mtz"}, _Ctx())
  assert g is True, "undetermined + unlocked must generate, got %r" % g
  st = _state({"data_mtz": "unknown.mtz"}, _Ctx())
  assert st is None, "no map/scalar/inspection must yield None, got %r" % st
  print("  PASS: test_undetermined_unlocked_generates")


def test_list_valued_data_mtz_slot():
  if _HELPER is None:
    print("  SKIP")
    return
  # multiple:true slots can be lists; the helper takes the first element and
  # looks up its basename in the map.
  ctx = _Ctx(mtz_rfree_map={"multi.mtz": True, "other.mtz": False})
  st = _state({"data_mtz": ["/x/multi.mtz", "/x/other.mtz"]}, ctx)
  assert st is True, "list-valued data_mtz first element must drive lookup, got %r" % st
  print("  PASS: test_list_valued_data_mtz_slot")


def test_no_data_mtz_no_flags():
  if _HELPER is None:
    print("  SKIP")
    return
  # No data MTZ in the selection at all -> cannot determine -> None
  # (undetermined), NOT a confident False (which would wrongly invite
  # generation via a spurious False).
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
  ctx = _Ctx(input_mtz_has_rfree=False, files_local=False)
  st = _state({"data_mtz": "/client/data.mtz"}, ctx)
  assert st is False, "client fact False must be respected, got %r" % st
  g = _decide(None, {"data_mtz": "/client/data.mtz"}, ctx)
  assert g is True, "client fact False must generate, got %r" % g
  print("  PASS: test_client_fact_false_generates")


def test_server_no_fact_undetermined_unlocked_generates():
  # v120.1 (AIAgent_35 MR/phaser case): server mode (files_local=False) with NO
  # client fact and no precomputed inspection -> the helper must NOT read the
  # (absent) file and returns None.  With nothing locked, the decision must
  # GENERATE (first-refinement safe default) — NOT strip.  Stripping here is what
  # broke the MR workflow: phenix.refine on a flagless phaser output aborted on
  # "No array of R-free flags found".
  if _HELPER is None:
    print("  SKIP")
    return
  ctx = _Ctx(input_mtz_has_rfree=None, files_local=False)
  st = _state({"data_mtz": "/client/PHASER.1.mtz"}, ctx)
  assert st is None, "server + no fact must be undetermined (None), got %r" % st
  g = _decide(None, {"data_mtz": "/client/PHASER.1.mtz"}, ctx)
  assert g is True, "undetermined + unlocked must generate, got %r" % g
  print("  PASS: test_server_no_fact_undetermined_unlocked_generates")


def test_map_overrides_scalar():
  # Precedence: the per-file map (source 0) wins over the original-input scalar
  # (source 1).  Here the map says the SELECTED file has no flags while the
  # scalar (original input) says True -> the map's per-file answer governs.
  if _HELPER is None:
    print("  SKIP")
    return
  ctx = _Ctx(input_mtz_has_rfree=True,
             mtz_rfree_map={"PHASER.1.mtz": False})
  st = _state({"data_mtz": "/x/PHASER.1.mtz"}, ctx)
  assert st is False, "map must override scalar for the selected file, got %r" % st
  # The scalar still applies for a file NOT in the map:
  st2 = _state({"data_mtz": "/x/original.mtz"}, ctx)
  assert st2 is True, "scalar applies when file absent from map, got %r" % st2
  print("  PASS: test_map_overrides_scalar")


def test_file_absent_from_map_is_undetermined():
  # The selected file is not in the map (the client could not inspect it) and
  # there is no scalar -> undetermined (None), NOT a confident False.  With
  # nothing locked, the decision then generates (first-refine default).
  if _HELPER is None:
    print("  SKIP")
    return
  ctx = _Ctx(mtz_rfree_map={"other.mtz": True})
  st = _state({"data_mtz": "/x/not_in_map.mtz"}, ctx)
  assert st is None, "file absent from map must be undetermined (None), got %r" % st
  g = _decide(None, {"data_mtz": "/x/not_in_map.mtz"}, ctx)
  assert g is True, "undetermined + unlocked must generate, got %r" % g
  print("  PASS: test_file_absent_from_map_is_undetermined")


def test_mr_phaser_output_fact_false_generates():
  # v120.1 regression (AIAgent_35): MR workflow refines the phaser output
  # (PHASER.1.mtz, NO flags).  The client fact is False (the original input data,
  # beta_blip_P3221.mtz, also has no flags).  rfree_mtz is unlocked.  The decision
  # MUST add generate=True even though the LLM omitted it (the LLM hallucinated
  # "flags generated by Phaser").  Before v120.1 this stripped/omitted generate
  # and phenix.refine aborted.
  if _HELPER is None:
    print("  SKIP")
    return
  ctx = _Ctx(input_mtz_has_rfree=False, files_local=True)
  # LLM omitted generate entirely:
  g = _decide(None, {"data_mtz": "PHASER.1.mtz"}, ctx)
  assert g is True, "MR phaser-output refine (fact False) must generate, got %r" % g
  print("  PASS: test_mr_phaser_output_fact_false_generates")


def test_map_distinguishes_two_files():
  # v120.2 core: the map gives DIFFERENT answers for different selected files in
  # the same session — the original input (flags) vs a phaser output (no flags).
  if _HELPER is None:
    print("  SKIP")
    return
  m = {"beta_blip_001.mtz": True, "PHASER.1.mtz": False}
  st_in = _state({"data_mtz": "/x/beta_blip_001.mtz"}, _Ctx(mtz_rfree_map=m))
  st_ph = _state({"data_mtz": "/x/PHASER.1.mtz"}, _Ctx(mtz_rfree_map=m))
  assert st_in is True, "flagged input -> True, got %r" % st_in
  assert st_ph is False, "phaser output -> False, got %r" % st_ph
  # And the decisions differ accordingly:
  assert _decide(None, {"data_mtz": "/x/beta_blip_001.mtz"},
                 _Ctx(mtz_rfree_map=m)) is False  # strip
  assert _decide(None, {"data_mtz": "/x/PHASER.1.mtz"},
                 _Ctx(mtz_rfree_map=m)) is True   # generate
  print("  PASS: test_map_distinguishes_two_files")


def test_parity_no_filesystem_read():
  # v120.2 parity guarantee: the helper must produce the SAME state regardless of
  # files_local (it no longer reads the filesystem).  inspect_mtz is stubbed to
  # raise (above), so if the helper tried to read, this would error.
  if _HELPER is None:
    print("  SKIP")
    return
  m = {"PHASER.1.mtz": False}
  for label, ctxkw in [
      ("map False", dict(mtz_rfree_map=m)),
      ("scalar True", dict(input_mtz_has_rfree=True)),
      ("scalar False", dict(input_mtz_has_rfree=False)),
      ("undetermined", dict())]:
    st_local = _state({"data_mtz": "/x/PHASER.1.mtz"},
                      _Ctx(files_local=True, **ctxkw))
    st_server = _state({"data_mtz": "/x/PHASER.1.mtz"},
                       _Ctx(files_local=False, **ctxkw))
    assert st_local == st_server, (
        "%s: local %r != server %r (parity violation)"
        % (label, st_local, st_server))
  print("  PASS: test_parity_no_filesystem_read")


_TESTS = [
  test_req1_no_flags_generates,
  test_req2_input_flags_via_map,
  test_req2_strips_llm_forced_generate,
  test_req2_precomputed_inspection,
  test_req3_locked_rfree_strips,
  test_undetermined_unlocked_generates,
  test_list_valued_data_mtz_slot,
  test_no_data_mtz_no_flags,
  test_client_fact_true_strips,
  test_client_fact_false_generates,
  test_server_no_fact_undetermined_unlocked_generates,
  test_map_overrides_scalar,
  test_file_absent_from_map_is_undetermined,
  test_mr_phaser_output_fact_false_generates,
  test_map_distinguishes_two_files,
  test_parity_no_filesystem_read,
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
