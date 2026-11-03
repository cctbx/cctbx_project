"""reduce2 on 6UF7, a cyclic peptide whose N-terminal PRO A 1 N is bonded to
C of the symmetry mate 3666 of DAR A 5 (LINK record; the asymmetric unit
holds half of the macrocycle). Neither end is a terminus: no NH3 at the
PRO N, and with qmready=True no OXT on the DAR C (it stopped with a "Bad
clash" between that OXT and the PRO N of the mate).
"""
import os
from io import StringIO
import iotbx.pdb
import mmtbx.model
from iotbx.cli_parser import run_program
from libtbx.utils import null_out
from mmtbx.programs import reduce2

pdb_str = """\
LINK         C   PRO A   1                 N   DSN A   2     1555   1555  1.35
LINK         N   PRO A   1                 C   DAR A   5     1555   3666  1.33
LINK         C   DSN A   2                 N   DGL A   3     1555   1555  1.34
LINK         C   DGL A   3                 N   ILE A   4     1555   1555  1.34
LINK         C   ILE A   4                 N   DAR A   5     1555   1555  1.34
CRYST1   12.200   10.980   25.540  90.00 114.90  90.00 P 1 21/c 1    4
ATOM      1  N   PRO A   1      -1.949   6.575  14.436  1.00  1.43           N
ATOM      2  CA  PRO A   1      -1.155   6.925  15.610  1.00  1.37           C
ATOM      3  C   PRO A   1      -0.626   5.702  16.352  1.00  1.26           C
ATOM      4  O   PRO A   1      -1.090   4.574  16.222  1.00  1.74           O
ATOM      5  CB  PRO A   1      -2.167   7.666  16.500  1.00  1.81           C
ATOM      6  CG  PRO A   1      -3.488   6.988  16.148  1.00  1.85           C
ATOM      7  CD  PRO A   1      -3.392   6.823  14.635  1.00  1.73           C
ATOM      8  HA  PRO A   1      -0.417   7.520  15.361  1.00  1.52           H
ATOM      9  HB2 PRO A   1      -1.955   7.555  17.440  1.00  2.84           H
ATOM     10  HB3 PRO A   1      -2.191   8.613  16.289  1.00  3.21           H
ATOM     11  HG2 PRO A   1      -3.569   6.128  16.590  1.00  2.21           H
ATOM     12  HG3 PRO A   1      -4.244   7.544  16.393  1.00  2.53           H
ATOM     13  HD2 PRO A   1      -3.924   6.072  14.328  1.00  2.46           H
ATOM     14  HD3 PRO A   1      -3.679   7.627  14.174  1.00  2.24           H
HETATM   15  N   DSN A   2       0.345   6.018  17.235  1.00  1.34           N
HETATM   16  CA  DSN A   2       0.765   5.089  18.267  1.00  1.36           C
HETATM   17  C   DSN A   2       1.979   4.235  17.916  1.00  1.42           C
HETATM   18  O   DSN A   2       2.309   3.320  18.665  1.00  2.10           O
HETATM   19  CB  DSN A   2       1.041   5.830  19.580  1.00  1.67           C
HETATM   20  OG  DSN A   2      -0.162   6.333  20.139  1.00  1.95           O
HETATM   21  H   DSN A   2       0.724   6.788  17.185  1.00  2.84           H
HETATM   22  HA  DSN A   2       0.016   4.478  18.430  1.00  2.16           H
HETATM   23  HB2 DSN A   2       1.653   6.564  19.415  1.00  2.12           H
HETATM   24  HB3 DSN A   2       1.462   5.224  20.210  1.00  1.95           H
HETATM   25  HG  DSN A   2      -0.290   7.096  19.868  1.00  3.01           H
HETATM   26  N   DGL A   3       2.616   4.539  16.778  1.00  1.43           N
HETATM   27  CA  DGL A   3       3.779   3.778  16.343  1.00  1.38           C
HETATM   28  C   DGL A   3       3.485   2.985  15.071  1.00  1.39           C
HETATM   29  O   DGL A   3       4.376   2.353  14.493  1.00  1.65           O
HETATM   30  CB  DGL A   3       5.000   4.691  16.206  1.00  1.53           C
HETATM   31  CG  DGL A   3       5.567   5.049  17.585  1.00  1.62           C
HETATM   32  CD  DGL A   3       6.133   3.798  18.269  1.00  1.60           C
HETATM   33  OE1 DGL A   3       7.043   3.181  17.644  1.00  1.98           O
HETATM   34  OE2 DGL A   3       5.648   3.446  19.364  1.00  2.08           O
HETATM   35  H   DGL A   3       2.334   5.193  16.297  1.00  1.89           H
HETATM   36  HA  DGL A   3       3.980   3.129  17.049  1.00  2.52           H
HETATM   37  HB2 DGL A   3       4.747   5.502  15.739  1.00  2.06           H
HETATM   38  HB3 DGL A   3       5.682   4.242  15.682  1.00  1.92           H
HETATM   39  HG2 DGL A   3       6.269   5.711  17.486  1.00  2.44           H
HETATM   40  HG3 DGL A   3       4.867   5.432  18.136  1.00  2.12           H
ATOM     41  N   ILE A   4       2.217   2.979  14.629  1.00  1.45           N
ATOM     42  CA  ILE A   4       1.872   2.325  13.371  1.00  1.46           C
ATOM     43  C   ILE A   4       2.729   2.885  12.240  1.00  1.41           C
ATOM     44  O   ILE A   4       3.159   2.157  11.342  1.00  1.79           O
ATOM     45  CB  ILE A   4       0.366   2.496  13.090  1.00  1.79           C
ATOM     46  CG1 ILE A   4      -0.492   1.840  14.190  1.00  2.38           C
ATOM     47  CG2 ILE A   4       0.015   1.961  11.700  1.00  2.21           C
ATOM     48  CD1 ILE A   4      -0.322   0.335  14.355  1.00  2.87           C
ATOM     49  H   ILE A   4       1.602   3.365  15.092  1.00  2.53           H
ATOM     50  HA  ILE A   4       2.063   1.368  13.459  1.00  2.47           H
ATOM     51  HB  ILE A   4       0.171   3.456  13.096  1.00  2.56           H
ATOM     52 HG12 ILE A   4      -0.280   2.265  15.036  1.00  3.14           H
ATOM     53 HG13 ILE A   4      -1.425   2.022  13.998  1.00  3.91           H
ATOM     54 HG21AILE A   4       0.567   2.391  11.043  0.80  3.71           H
ATOM     55 HG21BILE A   4       0.206   2.632  11.041  0.20  1.36           H
ATOM     56 HG22AILE A   4       0.165   1.014  11.674  0.80  2.61           H
ATOM     57 HG22BILE A   4       0.536   1.176  11.518  0.20  0.00           H
ATOM     58 HG23AILE A   4      -0.909   2.144  11.511  0.80  4.32           H
ATOM     59 HG23BILE A   4      -0.919   1.739  11.669  0.20  1.54           H
ATOM     60 HD11 ILE A   4      -0.895   0.024  15.059  1.00  3.42           H
ATOM     61 HD12 ILE A   4      -0.556  -0.105  13.534  1.00  3.87           H
ATOM     62 HD13 ILE A   4       0.591   0.138  14.574  1.00  4.55           H
HETATM   63  N   DAR A   5       2.963   4.203  12.252  1.00  1.54           N
HETATM   64  CA  DAR A   5       3.664   4.815  11.132  1.00  1.43           C
HETATM   65  CB  DAR A   5       3.968   6.279  11.451  1.00  1.50           C
HETATM   66  CG  DAR A   5       5.063   6.419  12.504  1.00  1.54           C
HETATM   67  CD  DAR A   5       5.085   7.787  13.159  1.00  1.62           C
HETATM   68  NE  DAR A   5       5.236   8.837  12.148  1.00  1.55           N
HETATM   69  CZ  DAR A   5       4.851  10.088  12.325  1.00  1.43           C
HETATM   70  NH1 DAR A   5       4.279  10.486  13.462  1.00  1.68           N
HETATM   71  NH2 DAR A   5       5.048  10.981  11.340  1.00  1.92           N
HETATM   72  C   DAR A   5       2.787   4.690   9.882  1.00  1.37           C
HETATM   73  O   DAR A   5       1.569   4.860   9.969  1.00  1.79           O
HETATM   74  H   DAR A   5       2.702   4.684  12.916  1.00  3.51           H
HETATM   75  HA  DAR A   5       4.506   4.338  10.981  1.00  1.11           H
HETATM   76  HB2 DAR A   5       3.160   6.711  11.771  1.00  1.68           H
HETATM   77  HB3 DAR A   5       4.246   6.732  10.639  1.00  3.12           H
HETATM   78  HG2 DAR A   5       5.924   6.254  12.088  1.00  2.78           H
HETATM   79  HG3 DAR A   5       4.931   5.744  13.188  1.00  2.23           H
HETATM   80  HD2 DAR A   5       5.822   7.834  13.788  1.00  2.02           H
HETATM   81  HD3 DAR A   5       4.260   7.926  13.650  1.00  2.66           H
HETATM   82  HE  DAR A   5       5.595   8.626  11.396  1.00  2.67           H
HETATM   83 HH11 DAR A   5       4.153   9.923  14.100  1.00  4.76           H
HETATM   84 HH12 DAR A   5       4.035  11.305  13.557  1.00  3.37           H
HETATM   85 HH21 DAR A   5       5.422  10.735  10.605  1.00  3.65           H
HETATM   86 HH22 DAR A   5       4.802  11.798  11.445  1.00  3.71           H
TER      87      DAR A   5
HETATM   88  O   HOH A 101      -0.769   4.650   8.594  1.00  1.96           O
HETATM   89  H1  HOH A 101      -0.869   3.867   8.273  1.00  3.21           H
HETATM   90  H2  HOH A 101       0.023   4.608   8.942  1.00  4.20           H
HETATM   91  O   HOH A 102      -0.298   8.979  19.305  1.00  2.13           O
HETATM   92  H1  HOH A 102      -1.043   9.378  19.147  1.00  5.65           H
HETATM   93  H2  HOH A 102       0.158   9.083  18.567  1.00  3.96           H
HETATM   94  O   HOH A 103       6.832   7.710   9.996  1.00  2.47           O
HETATM   95  H1  HOH A 103       6.302   7.714   9.304  1.00  5.16           H
HETATM   96  H2  HOH A 103       7.531   8.084   9.627  1.00  4.51           H
END
"""

def run(prefix, extra_args=(), pdb=pdb_str):
  """reduce2 on pdb; returns the output hierarchy."""
  with open(prefix + ".pdb", "w") as f:
    f.write(pdb)
  fn_out = prefix + "_out.pdb"
  run_program(program_class=reduce2.Program, logger=StringIO(), args=[
    prefix + ".pdb", "output.filename=" + fn_out,
    "output.description_file_name=" + prefix + ".txt", "overwrite=True"]
    + list(extra_args))
  assert os.path.isfile(fn_out)
  return iotbx.pdb.input(file_name=fn_out).construct_hierarchy()

def names(hierarchy, resseq):
  return sorted(a.name.strip() for a in hierarchy.atoms()
                if a.parent().parent().resseq_as_int() == resseq)

def check(hierarchy):
  pro = names(hierarchy, 1)
  assert not set(pro) & {"H", "H1", "H2", "H3"}, pro
  assert pro == sorted(["N", "CA", "C", "O", "CB", "CG", "CD", "HA", "HB2",
    "HB3", "HG2", "HG3", "HD2", "HD3"]), pro
  dar = names(hierarchy, 5)
  assert "OXT" not in dar, dar
  assert dar.count("H") == 1, dar

def test_add_terminal_oxt():
  # C of DAR 5 is bonded to N of a symmetry mate of PRO 1: no chain end
  model = mmtbx.model.manager(
    model_input=iotbx.pdb.input(source_info=None, lines=pdb_str),
    log=null_out())
  assert model.add_terminal_oxt() == 0

def test_add():
  check(run("tst_reduce2_symmetry_link_1"))

def test_qmready():
  check(run("tst_reduce2_symmetry_link_2", ["qmready=True"]))

def test_own_mate_is_no_bond():
  # half-occupied Gly, N 0.6 A from a 2-fold: its own mate is 1.2 A from the
  # N (disorder across the axis), a terminus with NH3 all the same
  h = run("tst_reduce2_symmetry_link_3", pdb="""\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1 2 1
ATOM      1  N   GLY A   1       0.600  10.000   0.000  0.50 20.00           N
ATOM      2  CA  GLY A   1       2.058  10.000   0.000  0.50 20.00           C
ATOM      3  C   GLY A   1       2.580  11.433   0.000  0.50 20.00           C
ATOM      4  O   GLY A   1       1.797  12.383   0.000  0.50 20.00           O
END
""")
  gly = names(h, 1)
  assert "H" not in gly, gly
  for name in ("H1", "H2", "H3"):
    assert gly.count(name) == 1, gly

def test_own_ca_mate_is_no_bond():
  # Gly along a 2-fold: the mate of its own CA is 1.69 A from the N
  h = run("tst_reduce2_symmetry_link_4", pdb="""\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1 2 1
ATOM      1  N   GLY A   1       0.600  10.000   0.000  0.50 20.00           N
ATOM      2  CA  GLY A   1       0.300  11.426   0.000  0.50 20.00           C
ATOM      3  C   GLY A   1       1.500  12.350   0.000  0.50 20.00           C
ATOM      4  O   GLY A   1       2.640  11.900   0.000  0.50 20.00           O
END
""")
  gly = names(h, 1)
  assert "H" not in gly, gly
  for name in ("H1", "H2", "H3"):
    assert gly.count(name) == 1, gly

def test_closure_with_own_mate():
  # cyclo(Gly-Gly) on a 2-fold, one Gly in the asymmetric unit: N is bonded to
  # C of its own mate, 1.44 A. Neither end is a terminus
  h = run("tst_reduce2_symmetry_link_5", ["qmready=True"], pdb="""\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1 2 1
ATOM      1  N   GLY A   1       1.440  10.000   0.000  1.00 20.00           N
ATOM      2  CA  GLY A   1       0.720  10.000   1.247  1.00 20.00           C
ATOM      3  C   GLY A   1      -0.720  10.000   1.247  1.00 20.00           C
ATOM      4  O   GLY A   1      -1.335  10.000   2.312  1.00 20.00           O
END
""")
  gly = names(h, 1)
  assert not set(gly) & {"H1", "H2", "H3", "OXT"}, gly

if __name__ == "__main__":
  test_add_terminal_oxt()
  test_add()
  test_qmready()
  test_own_mate_is_no_bond()
  test_own_ca_mate_is_no_bond()
  test_closure_with_own_mate()
  print("OK")
