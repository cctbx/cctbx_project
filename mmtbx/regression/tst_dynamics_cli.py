
from __future__ import division
from libtbx.utils import null_out
from cStringIO import StringIO

def exercise_dynamics_command () :
  pdb_in = """\
ATOM      1  N   SER A   1       0.043  -0.056   0.000  1.00  0.00
ATOM      2  CA  SER A   1       1.502  -0.020   0.000  1.00  0.00
ATOM      3  C   SER A   1       2.009   1.402   0.000  1.00  0.00
ATOM      4  O   SER A   1       1.244   2.367   0.000  1.00  0.00
ATOM      5  CB  SER A   1       2.073  -0.829  -1.192  1.00  0.00
ATOM      6  N   ASP A   2       3.321   1.549   0.000  1.00  0.00
ATOM      7  CA  ASP A   2       3.951   2.866   0.000  1.00  0.00
ATOM      8  C   ASP A   2       3.648   3.611  -1.278  1.00  0.00
ATOM      9  O   ASP A   2       3.392   4.855  -1.287  1.00  0.00
ATOM     10  CB  ASP A   2       5.477   2.730   0.242  1.00  0.00
ATOM     11  N   PRO A   3       3.664   2.893  -2.385  1.00  0.00
ATOM     12  CA  PRO A   3       3.393   3.399  -3.800  1.00  0.00
ATOM     13  C   PRO A   3       1.922   3.678  -4.096  1.00  0.00
ATOM     14  O   PRO A   3       1.642   4.698  -4.776  1.00  0.00
ATOM     15  CB  PRO A   3       3.863   2.300  -4.742  1.00  0.00
ATOM     16  N   ALA A   4       0.973   2.874  -3.653  1.00  0.00
ATOM     17  CA  ALA A   4      -0.440   3.120  -3.924  1.00  0.00
ATOM     18  C   ALA A   4      -0.904   4.394  -3.258  1.00  0.00
ATOM     19  O   ALA A   4      -1.686   5.175  -3.813  1.00  0.00
ATOM     20  CB  ALA A   4      -1.226   1.883  -3.459  1.00  0.00
ATOM     21  N   ALA A   5      -0.436   4.616  -2.044  1.00  0.00
ATOM     22  CA  ALA A   5      -0.800   5.806  -1.281  1.00  0.00
ATOM     23  C   ALA A   5      -0.280   7.057  -1.948  1.00  0.00
ATOM     24  O   ALA A   5      -0.943   8.100  -1.993  1.00  0.00
ATOM     25  CB  ALA A   5      -0.262   5.624   0.149  1.00  0.00
ATOM     26  N   LEU A   6       0.930   6.974  -2.469  1.00  0.00
ATOM     27  CA  LEU A   6       1.561   8.105  -3.144  1.00  0.00
ATOM     28  C   LEU A   6       0.807   8.474  -4.399  1.00  0.00
ATOM     29  O   LEU A   6       0.621   9.654  -4.732  1.00  0.00
ATOM     30  CB  LEU A   6       3.050   7.789  -3.463  1.00  0.00
ATOM     31  N   LYS A   7       0.357   7.464  -5.122  1.00  0.00
ATOM     32  CA  LYS A   7      -0.388   7.671  -6.360  1.00  0.00
ATOM     33  C   LYS A   7      -1.709   8.352  -6.090  1.00  0.00
ATOM     34  O   LYS A   7      -2.158   9.234  -6.836  1.00  0.00
ATOM     35  CB  LYS A   7      -0.602   6.306  -7.069  1.00  0.00
ATOM     36  N   ARG A   8      -2.356   7.951  -5.011  1.00  0.00
ATOM     37  CA  ARG A   8      -3.644   8.521  -4.625  1.00  0.00
ATOM     38  C   ARG A   8      -3.504   9.981  -4.267  1.00  0.00
ATOM     39  O   ARG A   8      -4.355  10.812  -4.578  1.00  0.00
ATOM     40  CB  ARG A   8      -4.223   7.699  -3.441  1.00  0.00
ATOM     41  N   ALA A   9      -2.427  10.307  -3.578  1.00  0.00
ATOM     42  CA  ALA A   9      -2.161  11.682  -3.163  1.00  0.00
ATOM     43  C   ALA A   9      -1.934  12.574  -4.359  1.00  0.00
ATOM     44  O   ALA A   9      -2.376  13.729  -4.409  1.00  0.00
ATOM     45  CB  ALA A   9      -0.959  11.654  -2.202  1.00  0.00
ATOM     46  N   ARG A  10      -1.225  12.054  -5.343  1.00  0.00
ATOM     47  CA  ARG A  10      -0.927  12.801  -6.561  1.00  0.00
ATOM     48  C   ARG A  10      -2.188  13.100  -7.335  1.00  0.00
ATOM     49  O   ARG A  10      -2.361  14.174  -7.908  1.00  0.00
ATOM     50  CB  ARG A  10       0.084  11.988  -7.416  1.00  0.00
ATOM     51  N   ASN A  11      -3.083  12.131  -7.381  1.00  0.00
ATOM     52  CA  ASN A  11      -4.349  12.279  -8.093  1.00  0.00
ATOM     53  C   ASN A  11      -5.209  13.341  -7.451  1.00  0.00
ATOM     54  O   ASN A  11      -5.871  14.135  -8.118  1.00  0.00
ATOM     55  CB  ASN A  11      -5.091  10.911  -8.123  1.00  0.00
"""
  open("tst_geo_min_ss_phil.pdb", "w").write(pdb_in)
  from mmtbx.command_line import secondary_structure_restraints
  out = StringIO()
  args = [
    "tst_geo_min_ss_phil.pdb",
    "use_ksdssp=False",
  ]
  secondary_structure_restraints.run(args, out=out, log=null_out())
  # phenix.geometry_minimization
  from mmtbx.command_line import dynamics
  open("tst_geo_min_ss_phil.phil", "w").write(out.getvalue())
  args = [
    "tst_geo_min_ss_phil.pdb",
    "tst_geo_min_ss_phil.phil",
    "secondary_structure_restraints=True",
    "secondary_structure.input.find_automatically=False",
  ]
  dynamics.run(args=args, log=null_out())
  pdb_new = "tst_geo_min_ss_phil_shaken.pdb"
  from iotbx import file_reader
  pdb_in = file_reader.any_file("tst_geo_min_ss_phil.pdb")
  pdb_in.assert_file_type("pdb")
  xrs_in = pdb_in.file_object.xray_structure_simple()
  pdb_out = file_reader.any_file(pdb_new)
  pdb_out.assert_file_type("pdb")
  xrs_out = pdb_out.file_object.xray_structure_simple()
  sites_in = xrs_in.sites_cart()
  sites_out = xrs_out.sites_cart()
  assert (sites_in.rms_difference(sites_out) > 0.3)

if (__name__ == "__main__") :
  exercise_dynamics_command()
  print "OK"
