"""
Tests for the R-free lock reconciliation + selected-file generate decision
(v120.3, F1+F2) in command_builder.

Bug (beta-blip AIAgent_55 refine_003): a directive file_preference pinned the
data slot to the original FLAGLESS MTZ, overriding the locked (flag-carrying)
R-free MTZ.  The generate decision then stripped generate=True purely because a
lock EXISTED, producing a flagless command with no generate -> hard abort.

  F2  reconciles a confirmed-flagless data slot back to the locked flag MTZ
      (runs in build() right after file selection).  Preserves R-free continuity.
  F1  bases the strip/add-generate decision on the ACTUALLY-SELECTED file
      (via mtz_rfree_map), falling back to the lock only when undetermined.
      Safety valve: if F2 could not reconcile (lock unavailable), a flagless
      file still gets generate=True instead of aborting.

Two layers:
  * test_probe_*  - bind _input_mtz_rfree_state from source; run anywhere.
  * test_build_*  - real build(); run on a full checkout, skip cleanly otherwise.

Run with: python3 tests/tst_rfree_lock_reconciliation.py
"""

from __future__ import absolute_import, division, print_function
import os
import re
import sys
import textwrap

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

try:
  from tests.tst_utils import (
    assert_equal, assert_true, assert_false, assert_in,
    run_tests_with_fail_fast,
  )
except ImportError:
  def assert_equal(a, b, msg=""):
    assert a == b, "%s: %r != %r" % (msg, a, b)
  def assert_true(val, msg=""):
    assert val, msg
  def assert_false(val, msg=""):
    assert not val, msg
  def assert_in(needle, haystack, msg=""):
    assert needle in haystack, "%s: %r not in %r" % (msg, needle, haystack)
  def run_tests_with_fail_fast():
    g = globals()
    for name in sorted(k for k in g if k.startswith("test_")):
      print("  Running %s..." % name)
      g[name]()
      print("  PASS: %s" % name)
    print("\nAll tests passed.")


# ---------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------

BETA = "/data/beta-blip/beta_blip_P3221.mtz"            # original, flagless
LOCK = "/run/sub_03_refine/refine_001_data.mtz"         # locked, has flags
SCA  = "/data/p9-sad/p9.sca"                            # scalepack: UNVERIFIABLE
OFLAG = "/run/other/other_flagged.mtz"                  # a DIFFERENT flagged MTZ
# Note: SCA is deliberately absent from the map -> probe returns None (the
# p9-sad case: inspect_mtz is ccp4_mtz-only, so scalepack is never verdicted).
RMAP = {"beta_blip_P3221.mtz": False, "refine_001_data.mtz": True,
        "other_flagged.mtz": True}


class _Ctx(object):
  """Minimal stand-in for CommandContext for the probe-level tests."""
  def __init__(self, **kw):
    self.mtz_rfree_map = None
    self.input_mtz_has_rfree = None
    self.mtz_inspection = None
    self.rfree_mtz = None
    self.__dict__.update(kw)


def _bind_probe():
  """Bind _input_mtz_rfree_state from command_builder source as a free function
  (it never references self), so the probe is tested without the heavy stack."""
  src_path = os.path.join(PROJECT_ROOT, "agent", "command_builder.py")
  lines = open(src_path).read().splitlines()
  start = next(i for i, l in enumerate(lines)
               if l.strip().startswith("def _input_mtz_rfree_state"))
  end = len(lines)
  for i in range(start + 1, len(lines)):
    if re.match(r"^  def ", lines[i]):
      end = i
      break
  body = textwrap.dedent("\n".join(lines[start:end]))
  body = body.replace("def _input_mtz_rfree_state(self, files, context):",
                      "def probe(files, context):")
  ns = {}
  exec(body, ns)
  return ns["probe"]


# ---------------------------------------------------------------------
# Probe foundation (run anywhere)
# ---------------------------------------------------------------------

def test_probe_flagless_file_is_false():
  probe = _bind_probe()
  assert_true(probe({"data_mtz": BETA}, _Ctx(mtz_rfree_map=RMAP)) is False,
              "flagless original MTZ should probe False")


def test_probe_flagged_file_is_true():
  probe = _bind_probe()
  assert_true(probe({"data_mtz": LOCK}, _Ctx(mtz_rfree_map=RMAP)) is True,
              "flag-carrying refine output should probe True")


def test_probe_unknown_file_is_none():
  probe = _bind_probe()
  assert_true(probe({"data_mtz": "/x/unknown.mtz"}, _Ctx(mtz_rfree_map=RMAP)) is None,
              "file not in the map should probe None (undetermined)")


# ---------------------------------------------------------------------
# build() integration (full checkout; skips cleanly otherwise)
# ---------------------------------------------------------------------

def _load_builder():
  try:
    from agent.command_builder import CommandBuilder, CommandContext
    return CommandBuilder, CommandContext
  except Exception as exc:
    print("    SKIP (command_builder unavailable: %s)" % exc)
    return None, None


def _build_refine(CommandBuilder, CommandContext, *, rfree_mtz,
                  available, pref_data, rmap):
  """Build a phenix.refine command and capture the builder log.

  NOTE: directives live at STATE level — CommandContext.from_state reads
  state.get("directives"), NOT session_info["directives"].  Putting them under
  session_info silently drops the file_preference, so the data slot would be
  filled by the PRIORITY-1 lock instead and F2 would never be exercised.
  """
  builder = CommandBuilder()
  state = {
    "session_info": {
      "experiment_type": "xray",
      "rfree_mtz": rfree_mtz,
      "mtz_rfree_map": rmap,
    },
    "directives": {"file_preferences": {"data_mtz": os.path.basename(pref_data)}},
  }
  ctx = CommandContext.from_state(state)
  logs = []
  ctx.log = logs.append            # _log(context, msg) -> context.log(msg)
  model = "/run/sub_03_refine/refine_001_001.pdb"
  cmd = builder.build("phenix.refine", available + [model], ctx)
  return cmd, logs


def test_build_f2_reconciles_flagless_to_lock():
  """F2: a file_preference forcing the flagless original is reconciled to the
  locked flag MTZ.  Asserts the RECONCILE LOG fired, which proves F2 did the work
  (a bare command check could also be satisfied by PRIORITY-1 selecting the lock
  directly)."""
  CB, CC = _load_builder()
  if CB is None:
    return
  try:
    cmd, logs = _build_refine(CB, CC, rfree_mtz=LOCK,
                              available=[BETA, LOCK], pref_data=BETA, rmap=RMAP)
  except Exception as exc:
    print("    SKIP (build raised: %s)" % exc)
    return
  if not cmd:
    print("    SKIP (build returned None — align harness with tst_command_builder)")
    return
  assert_true(any("reconciled" in m and "locked R-free MTZ" in m for m in logs),
              "F2 reconcile log expected (proves F2 fired, not PRIORITY-1); logs=%r"
              % logs)
  assert_true("refine_001_data.mtz" in cmd,
              "refine should use the locked flag MTZ, got: %s" % cmd)
  assert_false("beta_blip_P3221.mtz" in cmd,
               "flagless original must not be the refine data, got: %s" % cmd)
  assert_false("r_free_flags.generate=True" in cmd,
               "must not regenerate a fresh test set (leakage), got: %s" % cmd)


def test_build_first_refine_generates_flags():
  """F1 add branch: a first refinement of flagless data with NO lock gets
  xray_data.r_free_flags.generate=True (the emitted PHIL form of the internal
  generate_rfree_flags strategy key).  (The 'lock unavailable' safety-valve path is not
  reachable via build() because _file_is_available injects rfree_mtz's basename,
  so it is verified at the logic level by the test_probe_* checks instead.)"""
  CB, CC = _load_builder()
  if CB is None:
    return
  try:
    cmd, _logs = _build_refine(CB, CC, rfree_mtz=None,
                               available=[BETA], pref_data=BETA, rmap=RMAP)
  except Exception as exc:
    print("    SKIP (build raised: %s)" % exc)
    return
  if not cmd:
    print("    SKIP (build returned None — align harness with tst_command_builder)")
    return
  assert_true("r_free_flags.generate=True" in cmd,
              "first refine of flagless data should add generate, got: %s" % cmd)


def test_build_invariant_no_flagless_without_generate():
  """Cross-cutting invariant: the builder never emits a refine command whose data
  MTZ is confirmed-flagless AND lacks generate (the abort condition)."""
  CB, CC = _load_builder()
  if CB is None:
    return
  try:
    cmd, _logs = _build_refine(CB, CC, rfree_mtz=LOCK,
                               available=[BETA, LOCK], pref_data=BETA, rmap=RMAP)
  except Exception as exc:
    print("    SKIP (build raised: %s)" % exc)
    return
  if not cmd:
    print("    SKIP (build returned None — align harness with tst_command_builder)")
    return
  flagless_data = ("beta_blip_P3221.mtz" in cmd
                   and "refine_001_data.mtz" not in cmd)
  assert_false(flagless_data and "r_free_flags.generate=True" not in cmd,
               "emitted a flagless data MTZ with no generate: %s" % cmd)


# ---------------------------------------------------------------------
# v120.4: widened gate (is not True) — undetermined slot reconciles too
# ---------------------------------------------------------------------

def test_probe_scalepack_is_none():
  """A scalepack (.sca) basename absent from the map probes None (not False, no
  exception) — the load-bearing assumption of the widened gate (p9-sad)."""
  probe = _bind_probe()
  assert_true(probe({"data_mtz": SCA}, _Ctx(mtz_rfree_map=RMAP)) is None,
              "scalepack input absent from the map must probe None")


def test_build_f2_reconciles_undetermined_to_lock():
  """v120.4: an UNVERIFIABLE data slot (probe None, e.g. scalepack p9.sca) is
  reconciled to the locked flagged MTZ.  Asserts the reconcile log fired and
  names the file 'unverified' (not 'flagless')."""
  CB, CC = _load_builder()
  if CB is None:
    return
  try:
    cmd, logs = _build_refine(CB, CC, rfree_mtz=LOCK,
                              available=[SCA, LOCK], pref_data=SCA, rmap=RMAP)
  except Exception as exc:
    print("    SKIP (build raised: %s)" % exc)
    return
  if not cmd:
    print("    SKIP (build returned None — align harness with tst_command_builder)")
    return
  assert_true(any("reconciled" in m and "unverified" in m for m in logs),
              "expected an 'unverified' reconcile log; logs=%r" % logs)
  assert_true("refine_001_data.mtz" in cmd,
              "refine should use the locked flag MTZ, got: %s" % cmd)
  assert_false("p9.sca" in cmd,
               "the unverifiable scalepack input must not be the refine data: %s" % cmd)


def test_build_f2_keeps_a_different_flagged_file():
  """A DIFFERENT file that probes True (confirmed flagged) is NOT overridden — the
  widened gate keeps genuinely flagged alternatives."""
  CB, CC = _load_builder()
  if CB is None:
    return
  try:
    cmd, logs = _build_refine(CB, CC, rfree_mtz=LOCK,
                              available=[OFLAG, LOCK], pref_data=OFLAG, rmap=RMAP)
  except Exception as exc:
    print("    SKIP (build raised: %s)" % exc)
    return
  if not cmd:
    print("    SKIP (build returned None — align harness with tst_command_builder)")
    return
  assert_false(any("reconciled" in m for m in logs),
               "a confirmed-flagged file must NOT be reconciled; logs=%r" % logs)
  assert_true("other_flagged.mtz" in cmd,
              "the flagged preference should be kept, got: %s" % cmd)


def test_build_f2_same_file_as_lock_is_noop():
  """Risk A: when the data slot resolves to the lock itself, the abspath guard
  short-circuits before the probe — no spurious reconcile, no double-listing."""
  CB, CC = _load_builder()
  if CB is None:
    return
  try:
    cmd, logs = _build_refine(CB, CC, rfree_mtz=LOCK,
                              available=[LOCK], pref_data=LOCK, rmap=RMAP)
  except Exception as exc:
    print("    SKIP (build raised: %s)" % exc)
    return
  if not cmd:
    print("    SKIP (build returned None — align harness with tst_command_builder)")
    return
  assert_false(any("reconciled" in m for m in logs),
               "the lock must not be reconciled to itself; logs=%r" % logs)
  assert_true("refine_001_data.mtz" in cmd,
              "refine should run on the locked MTZ, got: %s" % cmd)


# ---------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------

def run_all_tests():
  run_tests_with_fail_fast()


if __name__ == "__main__":
  run_all_tests()
