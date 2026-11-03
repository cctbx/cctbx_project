"""Regression test: phenix.autobuild sequence-requirement waiver.

Bug (beta-blip AIAgent_17): phenix.autobuild lists `sequence` in inputs.required.
When the user asked to rebuild a placed model but supplied no sequence file, the
command builder failed with "missing: sequence" — even though autobuild can
rebuild a model "in place" and derive the sequence from the input model.

Fix (command_builder.py, just before the `if missing: return None` check):
when the ONLY missing required slot is `sequence` AND the program is
phenix.autobuild AND a model is available, waive the sequence requirement.  We do
NOT set rebuild_in_place — autobuild decides in-place vs not on its own; we only
stop REQUIRING the sequence in this case.  Sequence stays required when no model
is available and for all other programs; a run that DOES supply a sequence is
unaffected.

The waiver decision is self-contained logic; this test exercises it directly
(the full CommandBuilder pulls in ProgramRegistry + programs.yaml + a live
context, which aren't available in a bare sandbox).  A companion assertion checks
the program_registry placeholder-cleanup contract the waiver relies on (an absent
{sequence} is stripped, not left dangling).  2-space indent.
"""

from __future__ import absolute_import, division, print_function

import re


def _waive_sequence(missing, program, model_available):
  """Mirror of the Fix C waiver block in command_builder._select_files (kept in
  sync).  Mutates and returns the `missing` list."""
  if missing == ["sequence"] and program == "phenix.autobuild":
    if model_available:
      missing.remove("sequence")
  return missing


def _assemble(cmd_template, files, all_inputs):
  """Mirror of program_registry.build_command substitution + the unfilled-
  placeholder cleanup (the contract the waiver relies on)."""
  cmd = cmd_template
  for slot, path in files.items():
    ph = "{%s}" % slot
    if ph in cmd:
      flag = all_inputs.get(slot, {}).get("flag", "")
      cmd = cmd.replace(ph, ("%s%s" % (flag, path)) if flag else path)
  cmd = re.sub(r"\{[a-z_]+\}", "", cmd)   # remove unfilled placeholders
  cmd = " ".join(cmd.split())             # collapse whitespace
  return cmd


_AUTOBUILD_CMD = "phenix.autobuild {data_mtz} {sequence} {model} {map_coeffs}"
_AUTOBUILD_INPUTS = {
  "data_mtz": {"flag": "data="},
  "sequence": {"flag": "seq_file="},
  "model": {"flag": "model="},
  "map_coeffs": {"flag": "map_file="},
}


def test_autobuild_waives_sequence_when_model_available():
  """autobuild + only sequence missing + model available -> waived."""
  missing = _waive_sequence(["sequence"], "phenix.autobuild", True)
  assert missing == [], "sequence should be waived, got %r" % missing


def test_autobuild_keeps_sequence_when_no_model():
  """autobuild + only sequence missing + NO model -> still required (nothing to
  derive the sequence from)."""
  missing = _waive_sequence(["sequence"], "phenix.autobuild", False)
  assert missing == ["sequence"], \
    "without a model the sequence stays required, got %r" % missing


def test_waiver_does_not_fire_when_data_also_missing():
  """If data_mtz is ALSO missing, the build must still fail (waiver requires
  sequence to be the SOLE missing slot)."""
  missing = _waive_sequence(["data_mtz", "sequence"], "phenix.autobuild", True)
  assert "data_mtz" in missing, "data_mtz must remain missing, got %r" % missing
  assert missing == ["data_mtz", "sequence"], \
    "waiver should not fire when more than sequence is missing, got %r" % missing


def test_waiver_only_for_autobuild():
  """Other programs missing a sequence are unaffected by the waiver."""
  for prog in ("phenix.predict_and_build", "phenix.phaser",
               "phenix.autobuild_denmod"):
    missing = _waive_sequence(["sequence"], prog, True)
    assert missing == ["sequence"], \
      "%s should still require sequence, got %r" % (prog, missing)


def test_command_clean_when_sequence_waived():
  """After waiving, the assembled autobuild command has no stray seq_file= and no
  leftover {sequence} placeholder."""
  files = {"data_mtz": "/p/refine.mtz", "model": "/p/refine.pdb"}
  cmd = _assemble(_AUTOBUILD_CMD, files, _AUTOBUILD_INPUTS)
  assert cmd == "phenix.autobuild data=/p/refine.mtz model=/p/refine.pdb", \
    "unexpected command: %r" % cmd
  assert "seq_file=" not in cmd
  assert "{" not in cmd and "}" not in cmd
  assert "  " not in cmd and cmd == cmd.strip()


def test_command_includes_sequence_when_present():
  """Sanity: when a sequence IS supplied the command still includes seq_file=."""
  files = {"data_mtz": "/p/r.mtz", "sequence": "/p/s.fa", "model": "/p/r.pdb"}
  cmd = _assemble(_AUTOBUILD_CMD, files, _AUTOBUILD_INPUTS)
  assert "seq_file=/p/s.fa" in cmd, "unexpected command: %r" % cmd


_TESTS = [
  test_autobuild_waives_sequence_when_model_available,
  test_autobuild_keeps_sequence_when_no_model,
  test_waiver_does_not_fire_when_data_also_missing,
  test_waiver_only_for_autobuild,
  test_command_clean_when_sequence_waived,
  test_command_includes_sequence_when_present,
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
