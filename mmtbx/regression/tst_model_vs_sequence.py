from __future__ import absolute_import, division, print_function
import os
import tempfile
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import run_program
from libtbx.utils import null_out, Sorry
from libtbx.test_utils import approx_equal
from mmtbx.programs import model_vs_sequence

# Minimal PDB with a known sequence:
# chain A, residues 10-28
# (same as in tst_sequence_validation.py so expected output is well-understood).
_pdb_str = """\
ATOM      2  CA  ARG A  10      -6.299  36.344   7.806  1.00 55.20           C
ATOM     25  CA  TYR A  11      -3.391  33.962   7.211  1.00 40.56           C
ATOM     46  CA  ALA A  12      -0.693  34.802   4.693  1.00 67.95           C
ATOM     56  CA  ALA A  13       0.811  31.422   3.858  1.00 57.97           C
ATOM     66  CA  GLY A  14       4.466  31.094   2.905  1.00 49.24           C
ATOM     73  CA  ALA A  15       7.163  28.421   2.671  1.00 54.70           C
ATOM     83  CA  ILE A  16       6.554  24.685   2.957  1.00 51.79           C
ATOM    102  CA  LEU A  17       7.691  23.612   6.406  1.00 42.30           C
ATOM    121  CA  PTY A  18       7.292  19.882   5.861  1.00 36.68           C
ATOM    128  CA  PHE A  19       5.417  16.968   4.327  1.00 44.99           C
ATOM    148  CA  GLY A  20       3.466  14.289   6.150  1.00 41.99           C
ATOM    155  CA  GLY A  21       1.756  11.130   4.965  1.00 35.77           C
ATOM    190  CA  ALA A  24       1.294  19.658   3.683  1.00 47.02           C
ATOM    200  CA  VAL A  24A      2.361  22.009   6.464  1.00 37.13           C
ATOM    216  CA  HIS A  25       2.980  25.633   5.535  1.00 42.52           C
ATOM    234  CA  LEU A  26       4.518  28.425   7.577  1.00 47.63           C
ATOM    253  CA  ALA A  27       2.095  31.320   7.634  1.00 38.61           C
ATOM    263  CA  ARG A  28       1.589  34.719   9.165  1.00 37.04           C
END
"""

# Sequence that partially matches the model
# (produces mismatches/gaps).
_seq_mismatch = ">chain_a\nMTTPSHLSDRYELGEILGFGGMSEVHLARD\n"

# Simple 5-residue PDB and a perfectly matching sequence for 100%-identity test.
_pdb_simple = """\
ATOM      1  CA  GLY A   1       1.000   1.000   1.000  1.00  0.00           C
ATOM      2  CA  ALA A   2       2.000   1.000   1.000  1.00  0.00           C
ATOM      3  CA  VAL A   3       3.000   1.000   1.000  1.00  0.00           C
ATOM      4  CA  LEU A   4       4.000   1.000   1.000  1.00  0.00           C
ATOM      5  CA  ILE A   5       5.000   1.000   1.000  1.00  0.00           C
END
"""
_seq_match = ">chain_a\nGAVLI\n"

# ---------------------------------------------------------------------------

def _write_tmp(suffix, content):
  fd, path = tempfile.mkstemp(suffix=suffix)
  os.close(fd)
  with open(path, "w") as f:
    f.write(content)
  return path

# ---------------------------------------------------------------------------

def exercise_with_mismatches():
  """
  Run the program with a sequence that has mismatches and gaps relative to
  the model; verify chain-level statistics.
  """
  pdb_fn = _write_tmp(".pdb", _pdb_str)
  seq_fn = _write_tmp(".fa", _seq_mismatch)
  try:
    log = StringIO()
    result = run_program(
      program_class=model_vs_sequence.Program,
      args=[pdb_fn, seq_fn],
      logger=log)
    # result is the mmtbx.validation.sequence.validation object
    assert len(result.chains) == 1
    ch = result.chains[0]
    # The sequence has 30 residues; the model covers 17 standard residues
    # starting at position 10 → 9 residues missing at the start.
    assert ch.n_missing_start == 9, ch.n_missing_start
    assert ch.n_gaps > 0, "Expected gaps in alignment"
    assert len(ch.mismatch) > 0, "Expected mismatches"
    out_str = log.getvalue()
    assert "sequence identity" in out_str, out_str
  finally:
    os.remove(pdb_fn)
    os.remove(seq_fn)

# ---------------------------------------------------------------------------

def exercise_perfect_match():
  """
  Run the program with a sequence that matches the model exactly; verify
  100 % identity and no mismatches.
  """
  pdb_fn = _write_tmp(".pdb", _pdb_simple)
  seq_fn = _write_tmp(".fa", _seq_match)
  try:
    result = run_program(
      program_class=model_vs_sequence.Program,
      args=[pdb_fn, seq_fn],
      logger=null_out())
    assert len(result.chains) == 1
    ch = result.chains[0]
    assert approx_equal(ch.identity, 1.0, eps=1e-6), \
      "Expected 100%% identity, got %.4f" % ch.identity
    assert len(ch.mismatch) == 0, "Expected no mismatches"
  finally:
    os.remove(pdb_fn)
    os.remove(seq_fn)

# ---------------------------------------------------------------------------

def exercise_missing_sequence():
  """
  Omitting the sequence file must raise Sorry.
  """
  pdb_fn = _write_tmp(".pdb", _pdb_simple)
  try:
    try:
      run_program(
        program_class=model_vs_sequence.Program,
        args=[pdb_fn],
        logger=null_out())
    except Sorry:
      pass
    else:
      raise AssertionError("Expected Sorry when sequence file is missing")
  finally:
    os.remove(pdb_fn)

# ---------------------------------------------------------------------------

def exercise_missing_model():
  """
  Omitting the model file must raise Sorry.
  """
  seq_fn = _write_tmp(".fa", _seq_match)
  try:
    try:
      run_program(
        program_class=model_vs_sequence.Program,
        args=[seq_fn],
        logger=null_out())
    except Sorry:
      pass
    else:
      raise AssertionError("Expected Sorry when model file is missing")
  finally:
    os.remove(seq_fn)

# ---------------------------------------------------------------------------

def run():
  exercise_with_mismatches()
  exercise_perfect_match()
  exercise_missing_sequence()
  exercise_missing_model()
  print("OK")

# ---------------------------------------------------------------------------

if __name__ == "__main__":
  run()
