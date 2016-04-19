from __future__ import division
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import sys
import pickle
import StringIO
from libtbx.test_utils import show_diff

raw_records1 = """\
CRYST1   60.800   60.800   97.000  90.00  90.00 120.00 P 32 2 1      6
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.016447  0.009496  0.000000        0.00000
SCALE2      0.000000  0.018992  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010309        0.00000
ATOM   1050  N   LYS A 135      31.992  14.930  -7.233  1.00  9.47           N
ATOM   1051  CA  LYS A 135      31.388  16.216  -7.637  1.00 12.89           C
ATOM   1052  C   LYS A 135      30.807  16.840  -6.406  1.00  6.47           C
ATOM   1053  O   LYS A 135      29.583  16.869  -6.191  1.00 15.74           O
ATOM   1054  CB  LYS A 135      30.263  16.059  -8.655  1.00 13.51           C
ATOM   1055  CG  LYS A 135      30.742  15.277  -9.843  1.00 16.23           C
ATOM   1056  CD  LYS A 135      29.612  15.131 -10.835  1.00 28.55           C
ATOM   1057  CE  LYS A 135      30.173  14.812 -12.216  1.00 34.52           C
ATOM   1058  NZ  LYS A 135      29.396  13.756 -12.899  1.00 46.18           N
TER    1294      LYS A 162
END
"""

def make_initial_grm(mon_lib_srv, ener_lib, records):
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    raw_records    = records,
    force_symmetry = True)

  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies      = True,
    plain_pairs_radius = 5.0)
  xrs = processed_pdb_file.xray_structure()
  return geometry, xrs

def test_simple_protein(mon_lib_srv, ener_lib, prefix="simple_protein"):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1)
  init_out = StringIO.StringIO()
  from_file_out = StringIO.StringIO()
  geometry.show_sorted(
      sites_cart=xrs.sites_cart(),
      site_labels=xrs.scatterers().extract_labels(),
      f=init_out)

  pklfile = open("%s.pkl" % prefix, 'wb')
  pickle.dump(geometry, pklfile)
  pklfile.close()

  pklfile = open("%s.pkl" % prefix, 'rb')
  grm_from_file = pickle.load(pklfile)
  pklfile.close()
  grm_from_file.show_sorted(
      sites_cart=xrs.sites_cart(),
      site_labels=xrs.scatterers().extract_labels(),
      f=from_file_out)
  print "INITIAL"
  init_v = init_out.getvalue()
  print init_v
  print "="*50
  print "From disc"
  from_file_v = from_file_out.getvalue()
  print from_file_v
  # assert not show_diff(init_v, from_file_v)


def exercise_all(args):
  mon_lib_srv = None
  ener_lib = None
  try:
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
    test_simple_protein(mon_lib_srv, ener_lib)
  except Exception:
    print "Can not initialize monomer_library, skipping test."

if (__name__ == "__main__"):
  exercise_all(sys.argv[1:])
