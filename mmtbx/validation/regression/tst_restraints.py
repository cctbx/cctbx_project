
from __future__ import division
from libtbx.utils import null_out
from libtbx import easy_pickle
from cStringIO import StringIO

def run_validation (pdb_file, ignore_hd=True) :
  from mmtbx.validation import restraints
  import mmtbx.utils
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=[pdb_file],
    master_phil=mmtbx.utils.cmdline_input_phil_str,
    process_pdb_file=True,
    require_data=False,
    out=null_out())
  validation = restraints.combined(
    pdb_hierarchy=cmdline.pdb_hierarchy,
    xray_structure=cmdline.xray_structure,
    geometry_restraints_manager=cmdline.geometry,
    ignore_hd=ignore_hd)
  return validation

def exercise_simple () :
  # extracted from 1lyz, with hydrogens from reduce
  pdb_in = """
ATOM      1  N   LYS A   1       3.296   9.888  10.739  1.00  7.00           N
ATOM      2  CA  LYS A   1       2.439  10.217   9.791  1.00  6.00           C
ATOM      3  C   LYS A   1       2.439  11.997   9.160  1.00  6.00           C
ATOM      4  O   LYS A   1       2.637  12.656  10.107  1.00  8.00           O
ATOM      5  CB  LYS A   1       0.659  10.086   8.844  1.00  6.00           C
ATOM      6  CG  LYS A   1       0.198  10.415   8.086  1.00  6.00           C
ATOM      7  CD  LYS A   1      -1.187  10.086   8.212  1.00  6.00           C
ATOM      8  CE  LYS A   1      -2.175  10.086   7.264  1.00  6.00           C
ATOM      9  NZ  LYS A   1      -3.527   9.869   7.288  1.00  7.00           N
ATOM      0  H1  LYS A   1       3.156   9.045  10.986  1.00  7.00           H
ATOM      0  H2  LYS A   1       4.127   9.972  10.431  1.00  7.00           H
ATOM      0  H3  LYS A   1       3.184  10.425  11.440  1.00  7.00           H
ATOM      0  HA  LYS A   1       2.772   9.314   9.912  1.00  6.00           H
ATOM      0  HB2 LYS A   1       0.584   9.128   8.712  1.00  6.00           H
ATOM      0  HB3 LYS A   1       0.046  10.323   9.557  1.00  6.00           H
ATOM      0  HG2 LYS A   1       0.310  11.376   8.015  1.00  6.00           H
ATOM      0  HG3 LYS A   1       0.563  10.027   7.276  1.00  6.00           H
ATOM      0  HD2 LYS A   1      -1.193   9.186   8.573  1.00  6.00           H
ATOM      0  HD3 LYS A   1      -1.516  10.674   8.910  1.00  6.00           H
ATOM      0  HE2 LYS A   1      -2.097  10.964   6.860  1.00  6.00           H
ATOM      0  HE3 LYS A   1      -1.857   9.444   6.610  1.00  6.00           H
ATOM      0  HZ1 LYS A   1      -3.725   9.170   6.774  1.00  7.00           H
ATOM      0  HZ2 LYS A   1      -3.787   9.706   8.123  1.00  7.00           H
ATOM      0  HZ3 LYS A   1      -3.949  10.590   6.982  1.00  7.00           H
ATOM     10  N   VAL A   2       2.637  12.722   7.707  1.00  7.00           N
ATOM     11  CA  VAL A   2       2.307  14.172   7.580  1.00  6.00           C
ATOM     12  C   VAL A   2       0.857  14.041   6.949  1.00  6.00           C
ATOM     13  O   VAL A   2       0.659  13.843   5.875  1.00  8.00           O
ATOM     14  CB  VAL A   2       3.625  14.172   6.759  1.00  6.00           C
ATOM     15  CG1 VAL A   2       3.494  15.491   6.317  1.00  6.00           C
ATOM     16  CG2 VAL A   2       4.746  13.843   7.580  1.00  6.00           C
ATOM      0  H   VAL A   2       2.920  12.338   6.992  1.00  7.00           H
ATOM      0  HA  VAL A   2       2.195  14.925   8.181  1.00  6.00           H
ATOM      0  HB  VAL A   2       3.767  13.528   6.048  1.00  6.00           H
ATOM      0 HG11 VAL A   2       4.250  15.721   5.755  1.00  6.00           H
ATOM      0 HG12 VAL A   2       2.674  15.582   5.808  1.00  6.00           H
ATOM      0 HG13 VAL A   2       3.467  16.087   7.081  1.00  6.00           H
ATOM      0 HG21 VAL A   2       5.554  13.850   7.043  1.00  6.00           H
ATOM      0 HG22 VAL A   2       4.827  14.495   8.294  1.00  6.00           H
ATOM      0 HG23 VAL A   2       4.620  12.960   7.962  1.00  6.00           H
END
"""
  pdb_file = "tst_validate_restraints_simple.pdb"
  open(pdb_file, "w").write(pdb_in)
  v1 = run_validation(pdb_file, ignore_hd=True)
  out1 = StringIO()
  v1.show(out=out1)
  assert ("""
                       ----------Chiral volumes----------

atoms                   ideal    model    delta   sigma  residual   deviation
 A   1  LYS  CA
 A   1  LYS  N
 A   1  LYS  C
 A   1  LYS  CB          2.56     1.12     1.44  2.00e-01  5.19e+01   7.2*sigma
""" in "\n".join([ l.rstrip() for l in out1.getvalue().splitlines() ]))
  s = easy_pickle.dumps(v1)
  v1p = easy_pickle.loads(s)
  out1p = StringIO()
  v1p.show(out=out1p)
  assert (out1.getvalue() == out1p.getvalue())
  v2 = run_validation(pdb_file, ignore_hd=False)
  out2 = StringIO()
  v2.show(out=out2)
  assert (out2.getvalue() != out1.getvalue())
  assert ("""\
 A   1  LYS  HA        110.00    57.00    53.00  3.00e+00  3.12e+02  17.7*sigma
 A   2  VAL  N
 A   2  VAL  CA
""" in "\n".join([ l.rstrip() for l in out2.getvalue().splitlines() ]))

if (__name__ == "__main__") :
  exercise_simple()
  print "OK"
