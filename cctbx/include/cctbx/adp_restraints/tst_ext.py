from cctbx import adp_restraints
from libtbx.test_utils import approx_equal
from iotbx.shelx import from_ins

def exercise_rigid_bond():
  INSstructure = from_ins.from_ins(file_name = "enk_11i.res")

  """
  p = adp_restraints.rigid_bond_pair(
                    vec3<double> const& site1,
                    vec3<double> const& site2,
                    sym_mat3<double> const& ustar1,
                    sym_mat3<double> const& ustar2,
                    cctbx::uctbx::unit_cell const& uc)
  """

def exercise():
  exercise_rigid_bond()
  print "OK"

if (__name__ == "__main__"):
  exercise()
