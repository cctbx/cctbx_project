from libtbx.test_utils import approx_equal
from iotbx.shelx import from_ins
from cctbx.array_family import flex
from cctbx import adp_restraints
import math, os

result = [
  ['C1',   'C2',   0.0039,  0.0162,  0.0123],
  ['C1',   'N1',   0.0002,  0.0129,  0.0131],
  ['C2',   'C1',   0.0039,  0.0123,  0.0162],
  ['C3',   'C4',   0.0001,  0.0147,  0.0146],
  ['C3',   'C8',   0.0024,  0.0078,  0.0102],
  ['C4',   'C3',   0.0001,  0.0146,  0.0147],
  ['C4',   'C5',   0.0013,  0.0156,  0.0144],
  ['C5',   'C4',   0.0013,  0.0144,  0.0156],
  ['C5',   'C6',   0.0012,  0.0109,  0.0121],
  ['C6',   'C5',   0.0012,  0.0121,  0.0109],
  ['C6',   'C7',   0.0002,  0.0171,  0.0169],
  ['C6',   'O1',   0.0008,  0.0132,  0.0140],
  ['C7',   'C6',   0.0002,  0.0169,  0.0171],
  ['C7',   'C8',   0.0004,  0.0165,  0.0161],
  ['C8',   'C3',   0.0024,  0.0102,  0.0078],
  ['C8',   'C7',   0.0004,  0.0161,  0.0165],
  ['C9',   'O2',   0.0017,  0.0106,  0.0123],
  ['C11',  'O3',   0.0007,  0.0151,  0.0145],
  ['C11',  'N3',   0.0009,  0.0207,  0.0198],
  ['C12',  'C13',  0.0006,  0.0114,  0.0119],
  ['C12',  'N3',   0.0040,  0.0193,  0.0153],
  ['C13',  'C12',  0.0006,  0.0119,  0.0114],
  ['C13',  'O4',   0.0001,  0.0128,  0.0130],
  ['C13',  'N4',   0.0009,  0.0110,  0.0119],
  ['C14',  'N4',   0.0006,  0.0090,  0.0096],
  ['C16',  'C17',  0.0017,  0.0168,  0.0186],
  ['C16',  'C21',  0.0023,  0.0205,  0.0183],
  ['C17',  'C16',  0.0017,  0.0186,  0.0168],
  ['C17',  'C18',  0.0063,  0.0178,  0.0241],
  ['C18',  'C17',  0.0063,  0.0241,  0.0178],
  ['C18',  'C19',  0.0049,  0.0358,  0.0309],
  ['C19',  'C18',  0.0049,  0.0309,  0.0358],
  ['C19',  'C20',  0.0012,  0.0207,  0.0196],
  ['C20',  'C19',  0.0012,  0.0196,  0.0207],
  ['C20',  'C21',  0.0006,  0.0163,  0.0157],
  ['C21',  'C16',  0.0023,  0.0183,  0.0205],
  ['C21',  'C20',  0.0006,  0.0157,  0.0163],
  ['C22',  'N5',   0.0015,  0.0098,  0.0083],
  ['C23',  'C24',  0.0002,  0.0072,  0.0073],
  ['C24',  'C23',  0.0002,  0.0073,  0.0072],
  ['C25',  'C27',  0.0001,  0.0075,  0.0076],
  ['C27',  'C25',  0.0001,  0.0076,  0.0075],
  ['C28',  'O6',   0.0023,  0.0192,  0.0169],
  ['C28',  'O7',   0.0001,  0.0120,  0.0119],
  ['O1',   'C6',   0.0008,  0.0140,  0.0132],
  ['O2',   'C9',   0.0017,  0.0123,  0.0106],
  ['O3',   'C11',  0.0007,  0.0145,  0.0151],
  ['O4',   'C13',  0.0001,  0.0130,  0.0128],
  ['O6',   'C28',  0.0023,  0.0169,  0.0192],
  ['O7',   'C28',  0.0001,  0.0119,  0.0120],
  ['N1',   'C1',   0.0002,  0.0131,  0.0129],
  ['N3',   'C11',  0.0009,  0.0198,  0.0207],
  ['N3',   'C12',  0.0040,  0.0153,  0.0193],
  ['N4',   'C13',  0.0009,  0.0119,  0.0110],
  ['N4',   'C14',  0.0006,  0.0096,  0.0090],
  ['N5',   'C22',  0.0015,  0.0083,  0.0098]]

def exercise_rigid_bond_test():
  """
  Results compared with THMA11 (Ver. 20-04-91) - TLS Thermal Motion
  Analysis used as a part of WinGX (WinGX - Crystallographic Program
  System for Windows)
  """
  ins_file = os.path.expandvars('$LIBTBX_DIST_ROOT/regression/pdb/enk_11i.res')
  ins_xray_structure = from_ins.from_ins(file_name = ins_file)
  sites_frac = ins_xray_structure.sites_frac()
  sites_cart = ins_xray_structure.sites_cart()
  ustars = ins_xray_structure.scatterers().extract_u_star()
  scatterers = ins_xray_structure.scatterers()
  j = 0
  for site_cart_1,site_frac_1,ustar_1,scat_1 in zip(sites_cart,sites_frac,ustars,scatterers):
    for site_cart_2,site_frac_2,ustar_2, scat_2 in zip(sites_cart,sites_frac,ustars,scatterers):
      d = math.sqrt(flex.sum(flex.pow2(flex.double(site_cart_1)-\
                                       flex.double(site_cart_2))))
      if(d > 1.1 and d < 1.55):
        p = adp_restraints.rigid_bond_pair(site_frac_1,
                                           site_frac_2,
                                           ustar_1,
                                           ustar_2,
                                           ins_xray_structure.unit_cell())
        if(0):
          print "%4s %4s %7.4f %7.4f %7.4f" % \
                (scat_1.label,scat_2.label,p.delta_z(),p.z_12(),p.z_21())
        r = result[j]
        assert r[0] == scat_1.label
        assert r[1] == scat_2.label
        assert approx_equal(r[2], p.delta_z(), 1.e-4)
        assert approx_equal(r[3], p.z_12(), 1.e-4)
        assert approx_equal(r[4], p.z_21(), 1.e-4)
        j += 1
  assert j == 56

def exercise():
  exercise_rigid_bond_test()
  print "OK"

if (__name__ == "__main__"):
  exercise()
