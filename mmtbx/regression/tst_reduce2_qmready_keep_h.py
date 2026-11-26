"""reduce2 with qmready=True keeps the input H by default (keep_existing_H=Auto).
An explicit keep_existing_H=False discards them: the result is the same as for
the input without H.
"""
from __future__ import absolute_import, division, print_function
import os
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
import iotbx.pdb

# Asp-Glu (1ywf A 128-129, centred in the box). Without H: Asp -1, Glu -1,
# NH3+ +1, COO- -1 -> -2
pdb_str = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A   1       9.580   7.432   8.699  1.00 37.74           N
ATOM      2  CA  ASP A   1      10.692   8.369   8.630  1.00 38.33           C
ATOM      3  C   ASP A   1      10.828   9.231   9.884  1.00 37.17           C
ATOM      4  O   ASP A   1      11.938   9.511  10.322  1.00 36.50           O
ATOM      5  CB  ASP A   1      10.625   9.197   7.354  1.00 39.43           C
ATOM      6  CG  ASP A   1      11.281   8.467   6.194  1.00 43.72           C
ATOM      7  OD1 ASP A   1      10.549   7.831   5.406  1.00 49.92           O
ATOM      8  OD2 ASP A   1      12.531   8.382   6.067  1.00 51.10           O
%sATOM      9  N   GLU A   2       9.691   9.552  10.506  1.00 36.16           N
ATOM     10  CA  GLU A   2       9.668  10.262  11.762  1.00 35.61           C
ATOM     11  C   GLU A   2      10.374   9.397  12.850  1.00 34.61           C
ATOM     12  O   GLU A   2      11.228   9.895  13.566  1.00 33.43           O
ATOM     13  CB  GLU A   2       8.215  10.575  12.185  1.00 35.86           C
ATOM     14  CG  GLU A   2       7.535  11.724  11.472  1.00 39.43           C
ATOM     15  CD  GLU A   2       8.237  13.046  11.659  1.00 43.81           C
ATOM     16  OE1 GLU A   2       8.519  13.430  12.812  1.00 45.30           O
ATOM     17  OE2 GLU A   2       8.508  13.706  10.637  1.00 49.24           O
END
"""

# H on the Asp carboxyl: kept, it makes Asp neutral
hd2 = "ATOM      0  HD2 ASP A   1      13.000   9.200   5.900  1.00 51.10           H\n"

def run(model_text, prefix, extra_args=()):
  fn_in = prefix + ".pdb"
  fn_out = prefix + "_out.pdb"
  with open(fn_in, "w") as f:
    f.write(model_text)
  r = run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "overwrite=True", "qmready=True",
          "output.filename=" + fn_out] + list(extra_args))
  assert os.path.isfile(fn_out), fn_out
  with open(fn_out) as f:
    records = [l.rstrip() for l in f if l.startswith(("ATOM", "HETATM"))]
  return r.total_charge, records

def test_default_keeps_input_h():
  charge, records = run(pdb_str % hd2, "tst_reduce2_qmready_keep_h_1")
  assert charge == -1, charge
  assert len([r for r in records if " HD2 ASP" in r]) == 1

def test_explicit_true_keeps_input_h():
  charge, records = run(pdb_str % hd2, "tst_reduce2_qmready_keep_h_2",
    extra_args=["keep_existing_H=True"])
  assert charge == -1, charge
  assert len([r for r in records if " HD2 ASP" in r]) == 1

def test_explicit_false_discards_input_h():
  # the same model as from the input without H
  charge_ref, records_ref = run(pdb_str % "", "tst_reduce2_qmready_keep_h_3")
  assert charge_ref == -2, charge_ref
  charge, records = run(pdb_str % hd2, "tst_reduce2_qmready_keep_h_4",
    extra_args=["keep_existing_H=False"])
  assert charge == -2, charge
  assert not [r for r in records if " HD2 ASP" in r]
  assert records == records_ref

def test_auto_without_qmready_discards_input_h():
  # keep_existing_H=Auto is False without qmready, as the former default
  prefix = "tst_reduce2_qmready_keep_h_5"
  fn_in, fn_out = prefix + ".pdb", prefix + "_out.pdb"
  with open(fn_in, "w") as f:
    f.write(pdb_str % hd2)
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "overwrite=True", "output.filename=" + fn_out])
  atoms = iotbx.pdb.input(fn_out).construct_hierarchy().atoms()
  assert not [a for a in atoms if a.name.strip() == "HD2"
              and a.parent().resname == "ASP"]

if __name__ == "__main__":
  test_default_keeps_input_h()
  test_explicit_true_keeps_input_h()
  test_explicit_false_discards_input_h()
  test_auto_without_qmready_discards_input_h()
  print("OK")
