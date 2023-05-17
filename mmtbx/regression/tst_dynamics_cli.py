
from __future__ import absolute_import, division, print_function
from libtbx.utils import null_out
from six.moves import cStringIO as StringIO

def exercise_dynamics_command():
  pdb_in = """\
CRYST1   21.348   25.003   18.365  90.00  90.00  90.00 P 1
ATOM      1  N   SER A   1      10.914   5.773  13.123  1.00  0.00
ATOM      2  CA  SER A   1      12.373   5.809  13.123  1.00  0.00
ATOM      3  C   SER A   1      12.880   7.231  13.123  1.00  0.00
ATOM      4  O   SER A   1      12.115   8.196  13.123  1.00  0.00
ATOM      5  CB  SER A   1      12.944   5.000  11.931  1.00  0.00
ATOM      6  N   ASP A   2      14.192   7.378  13.123  1.00  0.00
ATOM      7  CA  ASP A   2      14.822   8.695  13.123  1.00  0.00
ATOM      8  C   ASP A   2      14.519   9.440  11.845  1.00  0.00
ATOM      9  O   ASP A   2      14.263  10.684  11.836  1.00  0.00
ATOM     10  CB  ASP A   2      16.348   8.559  13.365  1.00  0.00
ATOM     11  N   PRO A   3      14.535   8.722  10.738  1.00  0.00
ATOM     12  CA  PRO A   3      14.264   9.228   9.323  1.00  0.00
ATOM     13  C   PRO A   3      12.793   9.507   9.027  1.00  0.00
ATOM     14  O   PRO A   3      12.513  10.527   8.347  1.00  0.00
ATOM     15  CB  PRO A   3      14.734   8.129   8.381  1.00  0.00
ATOM     16  N   ALA A   4      11.844   8.703   9.470  1.00  0.00
ATOM     17  CA  ALA A   4      10.431   8.949   9.199  1.00  0.00
ATOM     18  C   ALA A   4       9.967  10.223   9.865  1.00  0.00
ATOM     19  O   ALA A   4       9.185  11.004   9.310  1.00  0.00
ATOM     20  CB  ALA A   4       9.645   7.712   9.664  1.00  0.00
ATOM     21  N   ALA A   5      10.435  10.445  11.079  1.00  0.00
ATOM     22  CA  ALA A   5      10.071  11.635  11.842  1.00  0.00
ATOM     23  C   ALA A   5      10.591  12.886  11.175  1.00  0.00
ATOM     24  O   ALA A   5       9.928  13.929  11.130  1.00  0.00
ATOM     25  CB  ALA A   5      10.609  11.453  13.272  1.00  0.00
ATOM     26  N   LEU A   6      11.801  12.803  10.654  1.00  0.00
ATOM     27  CA  LEU A   6      12.432  13.934   9.979  1.00  0.00
ATOM     28  C   LEU A   6      11.678  14.303   8.724  1.00  0.00
ATOM     29  O   LEU A   6      11.492  15.483   8.391  1.00  0.00
ATOM     30  CB  LEU A   6      13.921  13.618   9.660  1.00  0.00
ATOM     31  N   LYS A   7      11.228  13.293   8.001  1.00  0.00
ATOM     32  CA  LYS A   7      10.483  13.500   6.763  1.00  0.00
ATOM     33  C   LYS A   7       9.162  14.181   7.033  1.00  0.00
ATOM     34  O   LYS A   7       8.713  15.063   6.287  1.00  0.00
ATOM     35  CB  LYS A   7      10.269  12.135   6.054  1.00  0.00
ATOM     36  N   ARG A   8       8.515  13.780   8.112  1.00  0.00
ATOM     37  CA  ARG A   8       7.227  14.350   8.498  1.00  0.00
ATOM     38  C   ARG A   8       7.367  15.810   8.856  1.00  0.00
ATOM     39  O   ARG A   8       6.516  16.641   8.545  1.00  0.00
ATOM     40  CB  ARG A   8       6.648  13.528   9.682  1.00  0.00
ATOM     41  N   ALA A   9       8.444  16.136   9.545  1.00  0.00
ATOM     42  CA  ALA A   9       8.710  17.511   9.960  1.00  0.00
ATOM     43  C   ALA A   9       8.937  18.403   8.764  1.00  0.00
ATOM     44  O   ALA A   9       8.495  19.558   8.714  1.00  0.00
ATOM     45  CB  ALA A   9       9.912  17.483  10.921  1.00  0.00
ATOM     46  N   ARG A  10       9.646  17.883   7.780  1.00  0.00
ATOM     47  CA  ARG A  10       9.944  18.630   6.562  1.00  0.00
ATOM     48  C   ARG A  10       8.683  18.929   5.788  1.00  0.00
ATOM     49  O   ARG A  10       8.510  20.003   5.215  1.00  0.00
ATOM     50  CB  ARG A  10      10.955  17.817   5.707  1.00  0.00
ATOM     51  N   ASN A  11       7.788  17.960   5.742  1.00  0.00
ATOM     52  CA  ASN A  11       6.522  18.108   5.030  1.00  0.00
ATOM     53  C   ASN A  11       5.662  19.170   5.672  1.00  0.00
ATOM     54  O   ASN A  11       5.000  19.964   5.005  1.00  0.00
ATOM     55  CB  ASN A  11       5.780  16.740   5.000  1.00  0.00
TER
END
"""
  open("tst_geo_min_ss_phil.pdb", "w").write(pdb_in)
  from mmtbx.command_line import secondary_structure_restraints
  out = StringIO()
  args = [
    "tst_geo_min_ss_phil.pdb",
  ]
  secondary_structure_restraints.run(args, out=out, log=null_out())
  # phenix.geometry_minimization
  from mmtbx.command_line import dynamics
  open("tst_geo_min_ss_phil.phil", "w").write(out.getvalue())
  args = [
    "tst_geo_min_ss_phil.pdb",
    "tst_geo_min_ss_phil.phil",
    "secondary_structure.enabled=True",
  ]
  dynamics.run(args=args, log=null_out())
  pdb_new = "tst_geo_min_ss_phil_shaken.pdb"
  import iotbx.pdb
  pdb_in = iotbx.pdb.input("tst_geo_min_ss_phil.pdb")
  xrs_in = pdb_in.xray_structure_simple()
  pdb_out = iotbx.pdb.input(pdb_new)
  xrs_out = pdb_out.xray_structure_simple()
  sites_in = xrs_in.sites_cart()
  sites_out = xrs_out.sites_cart()
  assert (sites_in.rms_difference(sites_out) > 0.3)

if (__name__ == "__main__"):
  exercise_dynamics_command()
  print("OK")
