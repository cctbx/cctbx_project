
from __future__ import absolute_import, division, print_function
from mmtbx import disorder
import iotbx.pdb
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out

def exercise():
  pdb_str = """\
MODEL        1
ATOM      1  CA  ALA A  73      40.400   8.490  10.792  0.01 13.69           C
ATOM      2  CA  ARG A  74      36.672   7.864  11.359  0.01 13.55           C
ATOM      3  CA  LYS A  75      37.406   4.181  10.841  0.01 13.28           C
ENDMDL
MODEL        2
ATOM      1  CA  ALA A  73      40.320   8.758  11.103  0.01 13.69           C
ATOM      2  CA  ARG A  74      36.630   7.923  11.568  0.01 13.55           C
ATOM      3  CA  LYS A  75      36.856   4.141  11.040  0.01 13.28           C
ENDMDL
MODEL        3
ATOM      1  CA  ALA A  73      40.192   9.182  11.213  0.01 13.69           C
ATOM      2  CA  ARG A  74      36.644   7.871  11.649  0.01 13.55           C
ATOM      3  CA  LYS A  75      37.564   4.235  11.057  0.01 13.28           C
ENDMDL
MODEL        4
ATOM      1  CA  ALA A  73      39.924   9.073  10.840  0.01 13.69           C
ATOM      2  CA  ARG A  74      36.434   7.623  11.282  0.01 13.55           C
ATOM      3  CA  LYS A  75      38.120   4.233  11.095  0.01 13.28           C
ENDMDL
MODEL        5
ATOM      1  CA  ALA A  73      39.870   8.700  11.195  0.01 13.69           C
ATOM      2  CA  ARG A  74      36.328   7.379  11.083  0.01 13.55           C
ATOM      3  CA  LYS A  75      37.543   3.760  10.681  0.01 13.28           C
ENDMDL
MODEL        6
ATOM      1  CA  ALA A  73      40.236   9.111  11.122  0.01 13.69           C
ATOM      2  CA  ARG A  74      36.653   7.850  11.378  0.01 13.55           C
ATOM      3  CA  LYS A  75      37.252   4.117  10.669  0.01 13.28           C
ENDMDL
MODEL        7
ATOM      1  CA  ALA A  73      40.532   8.599  10.939  0.01 13.69           C
ATOM      2  CA  ARG A  74      36.830   7.868  11.208  0.01 13.55           C
ATOM      3  CA  LYS A  75      37.397   4.094  11.107  0.01 13.28           C
ENDMDL
MODEL        8
ATOM      1  CA  ALA A  73      40.340   9.520  11.303  0.01 13.69           C
ATOM      2  CA  ARG A  74      36.701   8.474  11.515  0.01 13.55           C
ATOM      3  CA  LYS A  75      37.522   4.781  11.633  0.01 13.28           C
ENDMDL
MODEL        9
ATOM      1  CA  ALA A  73      40.296   8.642  10.668  0.01 13.69           C
ATOM      2  CA  ARG A  74      36.832   7.269  11.324  0.01 13.55           C
ATOM      3  CA  LYS A  75      38.541   3.900  11.400  0.01 13.28           C
ENDMDL
MODEL       10
ATOM      1  CA  ALA A  73      39.396   8.850  11.064  0.01 13.69           C
ATOM      2  CA  ARG A  74      35.915   7.291  10.970  0.01 13.55           C
ATOM      3  CA  LYS A  75      37.173   3.719  10.662  0.01 13.28           C
ENDMDL
"""
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str)
  hierarchy = pdb_in.construct_hierarchy()
  disorder.set_ensemble_b_factors_to_xyz_displacement(
    pdb_hierarchy=hierarchy,
    method="rmsf",
    log=null_out())
  model = hierarchy.models()[0]
  assert approx_equal(list(model.atoms().extract_b()),
    [0.480574, 0.478505, 0.616683])
  disorder.set_ensemble_b_factors_to_xyz_displacement(
    pdb_hierarchy=hierarchy,
    method="mcs",
    log=null_out())
  model = hierarchy.models()[0]
  assert approx_equal(list(model.atoms().extract_b()),
    [0.684434, 0.779022, 0.881886])

if (__name__ == "__main__"):
  exercise()
  print("OK")
