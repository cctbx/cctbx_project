"""reduce2 with qmready=True calculates the total charge (as qr.finalise does
with calculate_charge=True) and appends it to the default output file name:
model.pdb -> modelH_charge-2.pdb.
"""
from __future__ import absolute_import, division, print_function
import os
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

# Asp-Glu: Asp -1, Glu -1, N-terminal NH3+ +1, C-terminal COO- -1 -> -2
pdb_str_neg = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A   1       5.000   5.000   5.000  1.00 10.00           N
ATOM      2  CA  ASP A   1       6.458   5.000   5.000  1.00 10.00           C
ATOM      3  C   ASP A   1       7.009   6.420   5.000  1.00 10.00           C
ATOM      4  O   ASP A   1       6.251   7.390   5.000  1.00 10.00           O
ATOM      5  CB  ASP A   1       6.994   4.250   3.780  1.00 10.00           C
ATOM      6  CG  ASP A   1       8.510   4.250   3.780  1.00 10.00           C
ATOM      7  OD1 ASP A   1       9.120   4.900   4.650  1.00 10.00           O
ATOM      8  OD2 ASP A   1       9.080   3.580   2.890  1.00 10.00           O
ATOM      9  N   GLU A   2       8.330   6.540   5.000  1.00 10.00           N
ATOM     10  CA  GLU A   2       8.960   7.850   5.000  1.00 10.00           C
ATOM     11  C   GLU A   2      10.470   7.700   5.000  1.00 10.00           C
ATOM     12  O   GLU A   2      11.000   6.600   5.000  1.00 10.00           O
ATOM     13  CB  GLU A   2       8.520   8.640   6.240  1.00 10.00           C
ATOM     14  CG  GLU A   2       8.950  10.100   6.250  1.00 10.00           C
ATOM     15  CD  GLU A   2       8.500  10.850   7.490  1.00 10.00           C
ATOM     16  OE1 GLU A   2       7.910  10.220   8.390  1.00 10.00           O
ATOM     17  OE2 GLU A   2       8.740  12.070   7.560  1.00 10.00           O
END
"""

# Ala-Ser-Lys: Lys +1, N-terminal NH3+ +1, C-terminal COO- -1 -> +1
pdb_str_pos = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      11.204   2.315  10.312  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.555   2.835  10.528  1.00 20.00           C
ATOM      3  C   ALA A   1      13.588   1.726  10.381  1.00 20.00           C
ATOM      4  O   ALA A   1      13.325   0.567  10.721  1.00 20.00           O
ATOM      5  CB  ALA A   1      12.866   3.965   9.556  1.00 20.00           C
ATOM      6  N   SER A   2      14.767   2.090   9.884  1.00 20.00           N
ATOM      7  CA  SER A   2      15.853   1.130   9.701  1.00 20.00           C
ATOM      8  C   SER A   2      16.968   1.370  10.713  1.00 20.00           C
ATOM      9  O   SER A   2      17.281   2.516  11.041  1.00 20.00           O
ATOM     10  CB  SER A   2      16.411   1.224   8.278  1.00 20.00           C
ATOM     11  OG  SER A   2      15.456   0.826   7.311  1.00 20.00           O
ATOM     12  N   LYS A   3      17.560   0.282  11.199  1.00 20.00           N
ATOM     13  CA  LYS A   3      18.643   0.370  12.172  1.00 20.00           C
ATOM     14  C   LYS A   3      19.945  -0.148  11.570  1.00 20.00           C
ATOM     15  O   LYS A   3      19.964  -1.176  10.892  1.00 20.00           O
ATOM     16  CB  LYS A   3      18.290  -0.424  13.432  1.00 20.00           C
ATOM     17  CG  LYS A   3      19.308  -0.301  14.556  1.00 20.00           C
ATOM     18  CD  LYS A   3      18.913  -1.137  15.762  1.00 20.00           C
ATOM     19  CE  LYS A   3      19.919  -1.002  16.894  1.00 20.00           C
ATOM     20  NZ  LYS A   3      19.531  -1.812  18.059  1.00 20.00           N
ATOM     21  OXT LYS A   3      20.995   0.451  11.763  1.00 20.00           O
END
"""

def run(model_text, prefix, expected_files, in_suffix=".pdb", extra_args=()):
  fn_in = prefix + in_suffix
  with open(fn_in, "w") as f:
    f.write(model_text)
  for fn in expected_files:
    if os.path.exists(fn): os.remove(fn)
  results = run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "overwrite=True"] + list(extra_args))
  for fn in expected_files:
    assert os.path.isfile(fn), fn
  return results

def test_qmready_negative_charge():
  prefix = "tst_reduce2_qmready_charge_1"
  r = run(pdb_str_neg, prefix,
    expected_files=[prefix+"H_charge-2.pdb", prefix+"H_charge-2.txt"],
    extra_args=["qmready=True"])
  assert r.total_charge == -2, r.total_charge

def test_qmready_positive_charge():
  # Same format as qr.finalise: no '+' sign
  prefix = "tst_reduce2_qmready_charge_2"
  r = run(pdb_str_pos, prefix,
    expected_files=[prefix+"H_charge1.pdb", prefix+"H_charge1.txt"],
    extra_args=["qmready=True"])
  assert r.total_charge == 1, r.total_charge

def test_qmready_cif():
  # mmCIF in -> mmCIF out, default name with charge
  import iotbx.pdb
  prefix = "tst_reduce2_qmready_charge_3"
  inp = iotbx.pdb.input(lines=pdb_str_neg.split("\n"), source_info=None)
  cif_str = inp.construct_hierarchy().as_mmcif_string(
    crystal_symmetry=inp.crystal_symmetry())
  r = run(cif_str, prefix,
    expected_files=[prefix+"H_charge-2.cif", prefix+"H_charge-2.txt"],
    in_suffix=".cif", extra_args=["qmready=True"])
  assert r.total_charge == -2, r.total_charge

def test_qmready_explicit_filename():
  prefix = "tst_reduce2_qmready_charge_4"
  fn_out = prefix + "_custom.pdb"
  r = run(pdb_str_neg, prefix, expected_files=[fn_out],
    extra_args=["qmready=True", "output.filename="+fn_out])
  assert r.total_charge == -2, r.total_charge
  assert not os.path.exists(prefix + "_custom_charge-2.pdb")

def test_no_qmready_no_charge():
  prefix = "tst_reduce2_qmready_charge_5"
  r = run(pdb_str_neg, prefix,
    expected_files=[prefix+"H.pdb", prefix+"H.txt"])
  assert r.total_charge is None, r.total_charge
  assert not os.path.exists(prefix+"H_charge-2.pdb")

if __name__ == "__main__":
  test_qmready_negative_charge()
  test_qmready_positive_charge()
  test_qmready_cif()
  test_qmready_explicit_filename()
  test_no_qmready_no_charge()
  print("OK")
