from libtbx.test_utils import approx_equal
from iotbx.shelx import from_ins

from cctbx import adp_restraints



def exercise_rigid_bond():
  ins_xray_structure = from_ins.from_ins(file_name = "enk_11i.res")

  ins_xray_structure.show_summary()

  sites_frac = ins_xray_structure.sites_frac()
  ustars = ins_xray_structure.scatterers().extract_u_star()

  site_1 = sites_frac[0]
  site_2 = sites_frac[1]

  ustar_1 = ustars[0]
  ustar_1 = ustars[1]

  """
  p = adp_restraints.rigid_bond_pair(site_1,
                                     site_2,
                                     ustar_1,
                                     ustar_2,
                                     ins_xray_structure.unit_cell())
  print p.delta_z()
  """

def exercise():
  exercise_rigid_bond()
  print "OK"

if (__name__ == "__main__"):
  exercise()
