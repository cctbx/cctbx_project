from __future__ import absolute_import, division, print_function
from mmtbx import monomer_library
import mmtbx.model
from six.moves import cStringIO as StringIO
#import libtbx.load_env
from libtbx.test_utils import show_diff, assert_lines_in_text
from libtbx.utils import null_out
import libtbx.load_env
import iotbx
import os

# from tst_grm_modifications import raw_records4, raw_records9
from tst_grm_modifications import make_initial_grm, show_sorted_geometry

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

def make_grm_via_model(mon_lib_srv, ener_lib, records, params):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=records)
  model = mmtbx.model.manager(
      model_input = pdb_inp,
      log=null_out())
  model.process(pdb_interpretation_params=params, make_restraints=True)
  grm = model.get_restraints_manager().geometry
  xrs = model.get_xray_structure()
  return grm, xrs


def show_sorted_geometry_str(geometry, xrs):
  sio = StringIO()
  geometry.show_sorted(
      sites_cart=xrs.sites_cart(),
      site_labels=xrs.scatterers().extract_labels(),
      f=sio)
  return sio.getvalue()

# SS-bridge with self!
raw_records1 = """\
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2
SCALE1      0.009182  0.000000  0.000000        0.00000
SCALE2      0.000000  0.009182  0.000000        0.00000
SCALE3      0.000000  0.000000  0.009182        0.00000
ATOM      1  N   PRO A   1     -50.110 -12.340  -8.990  1.00 87.52           N
ATOM      2  CA  PRO A   1     -49.518 -12.918 -10.200  1.00 88.41           C
ATOM      3  C   PRO A   1     -50.549 -13.676 -11.025  1.00 97.43           C
ATOM      4  O   PRO A   1     -51.129 -14.656 -10.556  1.00 93.89           O
ATOM      5  CB  PRO A   1     -48.446 -13.865  -9.649  1.00 78.46           C
ATOM      6  CG  PRO A   1     -48.110 -13.314  -8.317  1.00 80.77           C
ATOM      7  CD  PRO A   1     -49.402 -12.783  -7.777  1.00 92.02           C
ATOM      8  N   CYS A   2     -50.777 -13.205 -12.246  1.00 95.50           N
ATOM      9  CA  CYS A   2     -51.774 -13.797 -13.124  1.00103.68           C
ATOM     10  C   CYS A   2     -51.135 -14.327 -14.399  1.00101.77           C
ATOM     11  O   CYS A   2     -50.987 -13.597 -15.382  1.00106.09           O
ATOM     12  CB  CYS A   2     -52.861 -12.775 -13.466  1.00103.91           C
ATOM     13  SG  CYS A   2     -54.064 -13.357 -14.680  1.00106.39           S
ATOM     14  N   LYS A   3     -50.756 -15.600 -14.375  1.00105.98           N
ATOM     15  CA  LYS A   3     -50.135 -16.230 -15.528  1.00115.21           C
ATOM     16  C   LYS A   3     -50.110 -17.740 -15.357  1.00105.26           C
ATOM     17  O   LYS A   3     -50.331 -18.255 -14.264  1.00118.65           O
ATOM     18  CB  LYS A   3     -48.714 -15.704 -15.737  1.00110.03           C
ATOM     19  CG  LYS A   3     -47.784 -15.954 -14.563  1.00100.14           C
ATOM     20  CD  LYS A   3     -46.379 -15.452 -14.853  1.00101.77           C
ATOM     21  CE  LYS A   3     -45.456 -15.678 -13.663  1.00112.21           C
ATOM     22  NZ  LYS A   3     -44.063 -15.229 -13.941  1.00110.68           N
ATOM     23  OXT LYS A   3     -49.867 -18.479 -16.312  1.00 30.00           O
TER
END
"""

# CO is on special position, 3 copies around. Links to NE2 and HOH
raw_records2 = """\
CRYST1   81.430   81.430   33.750  90.00  90.00 120.00 H 3          18
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.012280  0.007090  0.000000        0.00000
SCALE2      0.000000  0.014180  0.000000        0.00000
SCALE3      0.000000  0.000000  0.029630        0.00000
ATOM      1  N   HIS D  10      -0.063   5.481  -5.607  1.00 10.05           N
ATOM      2  CA  HIS D  10       1.347   5.653  -5.385  1.00  9.54           C
ATOM      3  C   HIS D  10       1.835   7.103  -5.261  1.00  8.81           C
ATOM      4  O   HIS D  10       2.731   7.418  -4.480  1.00  8.89           O
ATOM      5  CB  HIS D  10       2.061   4.958  -6.530  1.00 10.08           C
ATOM      6  CG  HIS D  10       1.686   3.518  -6.637  1.00 10.90           C
ATOM      7  ND1 HIS D  10       1.964   2.609  -5.629  1.00 11.69           N
ATOM      8  CD2 HIS D  10       0.995   2.864  -7.574  1.00 11.16           C
ATOM      9  CE1 HIS D  10       1.490   1.429  -5.986  1.00 12.82           C
ATOM     10  NE2 HIS D  10       0.918   1.565  -7.171  1.00 10.61           N
HETATM   11 CO    CO D 101      -0.002   0.003  -8.336  0.33 12.86          CO
HETATM   12  O   HOH D 208      -0.901   1.609  -9.921  1.00 16.66           O
END
"""

raw_records4 = """\
CRYST1   15.775   12.565   13.187  90.00  90.00  90.00 P 1
ATOM      1  N   MET A   1       9.821   6.568   5.000  1.00 66.07           N
ATOM      2  CA  MET A   1       9.946   7.171   6.357  1.00 66.55           C
ATOM      3  C   MET A   1      10.571   6.157   7.305  1.00 64.57           C
ATOM      4  O   MET A   1      10.775   5.000   6.933  1.00 66.25           O
ATOM      5  CB  MET A   1       8.570   7.565   6.897  1.00 69.08           C
ATOM      6  CG  MET A   1       7.723   6.373   7.299  1.00 71.37           C
ATOM      7  SD  MET A   1       6.247   6.862   8.187  1.00 76.22           S
ATOM      8  CE  MET A   1       5.000   6.694   6.892  1.00 74.93           C
END
""".splitlines()

raw_records9 = """\
CRYST1   41.028   41.028  183.010  90.00  90.00  90.00 P 43 21 2
SCALE1      0.024374  0.000000  0.000000        0.00000
SCALE2      0.000000  0.024374  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005464        0.00000
ATOM      1  CG  HIS A 319       2.304  23.849  77.123  1.00 23.70           C
ATOM      2  ND1 HIS A 319       1.668  23.871  75.903  1.00 25.36           N
ATOM      3  CD2 HIS A 319       2.302  25.128  77.578  1.00 21.49           C
ATOM      4  CE1 HIS A 319       1.265  25.110  75.643  1.00 25.28           C
ATOM      5  NE2 HIS A 319       1.643  25.884  76.636  1.00 23.98           N
ATOM      6  HD1 HIS A 319       1.550  23.191  75.390  1.00 30.43           H
ATOM      7  HD2 HIS A 319       2.675  25.435  78.373  1.00 25.78           H
ATOM      8  HE1 HIS A 319       0.795  25.383  74.888  1.00 30.33           H
ATOM      9  HG3 MET A 338      -1.284  24.273  77.766  1.00 36.33           H
ATOM     10  CG  GLU A 362       1.743  29.061  80.665  1.00 32.98           C
ATOM     11  CD  GLU A 362       1.505  28.476  79.273  1.00 26.00           C
ATOM     12  OE1 GLU A 362       0.357  28.085  78.927  1.00 33.20           O
ATOM     13  OE2 GLU A 362       2.511  28.378  78.586  1.00 24.88           O
ATOM     14  HG2 GLU A 362       2.252  29.880  80.555  1.00 39.58           H
ATOM     15  HG3 GLU A 362       2.259  28.414  81.171  1.00 39.58           H
TER
ATOM     16  N   HIS B 304     -20.949  11.990  59.962  1.00 23.21           N
ATOM     17  CA  HIS B 304     -22.165  11.249  60.299  1.00 25.44           C
ATOM     18  C   HIS B 304     -23.349  12.161  60.500  1.00 29.20           C
ATOM     19  O   HIS B 304     -24.477  11.760  60.211  1.00 32.52           O
ATOM     20  CB  HIS B 304     -21.963  10.373  61.537  1.00 28.40           C
ATOM     21  CG  HIS B 304     -20.809   9.431  61.430  1.00 27.87           C
ATOM     22  ND1 HIS B 304     -19.517   9.806  61.733  1.00 36.76           N
ATOM     23  CD2 HIS B 304     -20.737   8.128  61.054  1.00 26.50           C
ATOM     24  CE1 HIS B 304     -18.706   8.776  61.549  1.00 36.85           C
ATOM     25  NE2 HIS B 304     -19.428   7.737  61.176  1.00 27.21           N
ATOM     26  HA  HIS B 304     -22.359  10.653  59.559  1.00 30.53           H
ATOM     27  HB2 HIS B 304     -21.805  10.948  62.303  1.00 34.08           H
ATOM     28  HB3 HIS B 304     -22.764   9.845  61.678  1.00 34.08           H
ATOM     29  HD2 HIS B 304     -21.445   7.599  60.766  1.00 31.79           H
ATOM     30  HE1 HIS B 304     -17.783   8.783  61.663  1.00 44.22           H
ATOM     31  HE2 HIS B 304     -19.128   6.944  61.034  1.00 32.64           H
ATOM     32 HG21 VAL B 315     -20.236   8.723  64.319  1.00 36.99           H
TER
HETATM   33 ZN    ZN A   8       1.797  27.888  76.692  1.00 34.36          Zn
HETATM   34  O   HOH A  57       3.694  28.411  75.816  1.00 39.92           O
HETATM   35  O   HOH A  69       5.784  26.244  75.932  1.00 38.11           O
"""

raw_records91 = """\
CRYST1   41.028   41.028  183.010  90.00  90.00  90.00 P 43 21 2
SCALE1      0.024374  0.000000  0.000000        0.00000
SCALE2      0.000000  0.024374  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005464        0.00000
ATOM     16  N   HIS B 304     -20.949  11.990  59.962  1.00 23.21           N
ATOM     17  CA  HIS B 304     -22.165  11.249  60.299  1.00 25.44           C
ATOM     18  C   HIS B 304     -23.349  12.161  60.500  1.00 29.20           C
ATOM     19  O   HIS B 304     -24.477  11.760  60.211  1.00 32.52           O
ATOM     20  CB  HIS B 304     -21.963  10.373  61.537  1.00 28.40           C
ATOM     21  CG  HIS B 304     -20.809   9.431  61.430  1.00 27.87           C
ATOM     22  ND1 HIS B 304     -19.517   9.806  61.733  1.00 36.76           N
ATOM     23  CD2 HIS B 304     -20.737   8.128  61.054  1.00 26.50           C
ATOM     24  CE1 HIS B 304     -18.706   8.776  61.549  1.00 36.85           C
ATOM     25  NE2 HIS B 304     -19.428   7.737  61.176  1.00 27.21           N
ATOM     26  HA  HIS B 304     -22.359  10.653  59.559  1.00 30.53           H
ATOM     27  HB2 HIS B 304     -21.805  10.948  62.303  1.00 34.08           H
ATOM     28  HB3 HIS B 304     -22.764   9.845  61.678  1.00 34.08           H
ATOM     29  HD2 HIS B 304     -21.445   7.599  60.766  1.00 31.79           H
ATOM     30  HE1 HIS B 304     -17.783   8.783  61.663  1.00 44.22           H
ATOM     31  HE2 HIS B 304     -19.128   6.944  61.034  1.00 32.64           H
ATOM     32 HG21 VAL B 315     -20.236   8.723  64.319  1.00 36.99           H
TER
HETATM   33 ZN    ZN A   8       1.797  27.888  76.692  1.00 34.36          Zn
HETATM   34  O   HOH A  57       3.694  28.411  75.816  1.00 39.92           O
HETATM   35  O   HOH A  69       5.784  26.244  75.932  1.00 38.11           O
"""

raw_records91_noh = """\
CRYST1   41.028   41.028  183.010  90.00  90.00  90.00 P 43 21 2
SCALE1      0.024374  0.000000  0.000000        0.00000
SCALE2      0.000000  0.024374  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005464        0.00000
ATOM     16  N   HIS B 304     -20.949  11.990  59.962  1.00 23.21           N
ATOM     17  CA  HIS B 304     -22.165  11.249  60.299  1.00 25.44           C
ATOM     18  C   HIS B 304     -23.349  12.161  60.500  1.00 29.20           C
ATOM     19  O   HIS B 304     -24.477  11.760  60.211  1.00 32.52           O
ATOM     20  CB  HIS B 304     -21.963  10.373  61.537  1.00 28.40           C
ATOM     21  CG  HIS B 304     -20.809   9.431  61.430  1.00 27.87           C
ATOM     22  ND1 HIS B 304     -19.517   9.806  61.733  1.00 36.76           N
ATOM     23  CD2 HIS B 304     -20.737   8.128  61.054  1.00 26.50           C
ATOM     24  CE1 HIS B 304     -18.706   8.776  61.549  1.00 36.85           C
ATOM     25  NE2 HIS B 304     -19.428   7.737  61.176  1.00 27.21           N
TER
HETATM   33 ZN    ZN A   8       1.797  27.888  76.692  1.00 34.36          Zn
HETATM   34  O   HOH A  57       3.694  28.411  75.816  1.00 39.92           O
HETATM   35  O   HOH A  69       5.784  26.244  75.932  1.00 38.11           O
"""


def exercise_remove_bond_restraint_in_place(mon_lib_srv, ener_lib):
  # removing the bond between N and CA, iseqs being (0,1)
  # making sure nothing else changes
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records4)

  assert geometry.is_bonded_atoms(0,1)
  assert geometry.pair_proxies().bond_proxies.simple.size() == 7
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  assert geometry.pair_proxies().nonbonded_proxies.simple.size() == 9
  assert geometry.pair_proxies().nonbonded_proxies.asu.size() == 0
  show_sorted_geometry(geometry, xrs, 'exercise_remove_bond_restraint_in_place_start.geo')

  # Removing
  geometry.remove_bond_restraints_in_place(bonded_pairs=[(0,1)], sites_cart=xrs.sites_cart())

  show_sorted_geometry(geometry, xrs, 'exercise_remove_bond_restraint_in_place_end.geo')
  final_geo = show_sorted_geometry_str(geometry, xrs)

  # Checking
  assert not geometry.is_bonded_atoms(0,1)
  # 1 less bond:
  assert geometry.pair_proxies().bond_proxies.simple.size() == 6, geometry.pair_proxies().bond_proxies.simple.size()
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0, geometry.pair_proxies().bond_proxies.asu.size()
  assert geometry.pair_proxies().nonbonded_proxies.simple.size() == 12, geometry.pair_proxies().nonbonded_proxies.simple.size()
  assert geometry.pair_proxies().nonbonded_proxies.asu.size() == 0, geometry.pair_proxies().nonbonded_proxies.asu.size()
  # Nonbonded - we expect not only additional N-CA, but also
  # N-C as former 1-3 interaction, N-O and N-CB as former 1-4 interactions now.
  # print(final_geo)
  n_pos = final_geo.find("Nonbonded | unspecified")
  assert n_pos>0
  assert not show_diff(final_geo[n_pos:-2], """\
Nonbonded | unspecified | interactions: 12
Sorted by model distance:
nonbonded pdb=" N   MET A   1 "
          pdb=" CA  MET A   1 "
   model   vdw
   1.490 3.550
nonbonded pdb=" N   MET A   1 "
          pdb=" C   MET A   1 "
   model   vdw
   2.459 3.350
nonbonded pdb=" N   MET A   1 "
          pdb=" CB  MET A   1 "
   model   vdw
   2.481 3.520
nonbonded pdb=" N   MET A   1 "
          pdb=" O   MET A   1 "
   model   vdw
   2.666 3.120
nonbonded pdb=" C   MET A   1 "
          pdb=" CG  MET A   1 "
   model   vdw
   2.856 2.936
nonbonded pdb=" N   MET A   1 "
          pdb=" CG  MET A   1 "
   model   vdw
   3.118 3.520
nonbonded pdb=" O   MET A   1 "
          pdb=" CG  MET A   1 "
   model   vdw
   3.367 3.440
nonbonded pdb=" O   MET A   1 "
          pdb=" CB  MET A   1 "
   model   vdw
   3.383 2.752
nonbonded pdb=" CB  MET A   1 "
          pdb=" CE  MET A   1 "
   model   vdw
   3.675 3.088
nonbonded pdb=" CA  MET A   1 "
          pdb=" SD  MET A   1 "
   model   vdw
   4.138 3.064
nonbonded pdb=" C   MET A   1 "
          pdb=" SD  MET A   1 "
   model   vdw
   4.469 3.630
nonbonded pdb=" N   MET A   1 "
          pdb=" SD  MET A   1 "
   model   vdw
   4.798 3.480""")

def exercise_remove_two_bond_restraints_in_place(mon_lib_srv, ener_lib):
  # removing the bond between N and CA, iseqs being (0,1)
  # making sure nothing else changes
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records4)

  assert geometry.is_bonded_atoms(0,1)
  assert geometry.pair_proxies().bond_proxies.simple.size() == 7
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  assert geometry.pair_proxies().nonbonded_proxies.simple.size() == 9
  assert geometry.pair_proxies().nonbonded_proxies.asu.size() == 0
  show_sorted_geometry(geometry, xrs, 'exercise_remove_two_bond_restraints_in_place_start.geo')

  # Removing
  geometry.remove_bond_restraints_in_place(bonded_pairs=[(0,1), (3,2)], sites_cart=xrs.sites_cart())

  show_sorted_geometry(geometry, xrs, 'exercise_remove_two_bond_restraints_in_place_end.geo')
  final_geo = show_sorted_geometry_str(geometry, xrs)
  # print(final_geo)

  # Checking
  assert not geometry.is_bonded_atoms(0,1)
  # 1 less bond:
  assert geometry.pair_proxies().bond_proxies.simple.size() == 5, geometry.pair_proxies().bond_proxies.simple.size()
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0, geometry.pair_proxies().bond_proxies.asu.size()
  assert geometry.pair_proxies().nonbonded_proxies.simple.size() == 14, geometry.pair_proxies().nonbonded_proxies.simple.size()
  assert geometry.pair_proxies().nonbonded_proxies.asu.size() == 0, geometry.pair_proxies().nonbonded_proxies.asu.size()
  # print(final_geo)
  n_pos = final_geo.find("Nonbonded | unspecified")
  assert n_pos>0
  assert not show_diff(final_geo[n_pos:-2], """\
Nonbonded | unspecified | interactions: 14
Sorted by model distance:
nonbonded pdb=" C   MET A   1 "
          pdb=" O   MET A   1 "
   model   vdw
   1.232 3.270
nonbonded pdb=" N   MET A   1 "
          pdb=" CA  MET A   1 "
   model   vdw
   1.490 3.550
nonbonded pdb=" CA  MET A   1 "
          pdb=" O   MET A   1 "
   model   vdw
   2.394 3.470
nonbonded pdb=" N   MET A   1 "
          pdb=" C   MET A   1 "
   model   vdw
   2.459 3.350
nonbonded pdb=" N   MET A   1 "
          pdb=" CB  MET A   1 "
   model   vdw
   2.481 3.520
nonbonded pdb=" N   MET A   1 "
          pdb=" O   MET A   1 "
   model   vdw
   2.666 3.120
nonbonded pdb=" C   MET A   1 "
          pdb=" CG  MET A   1 "
   model   vdw
   2.856 2.936
nonbonded pdb=" N   MET A   1 "
          pdb=" CG  MET A   1 "
   model   vdw
   3.118 3.520
nonbonded pdb=" O   MET A   1 "
          pdb=" CG  MET A   1 "
   model   vdw
   3.367 3.440
nonbonded pdb=" O   MET A   1 "
          pdb=" CB  MET A   1 "
   model   vdw
   3.383 3.440
nonbonded pdb=" CB  MET A   1 "
          pdb=" CE  MET A   1 "
   model   vdw
   3.675 3.088
nonbonded pdb=" CA  MET A   1 "
          pdb=" SD  MET A   1 "
   model   vdw
   4.138 3.064
nonbonded pdb=" C   MET A   1 "
          pdb=" SD  MET A   1 "
   model   vdw
   4.469 3.630
nonbonded pdb=" N   MET A   1 "
          pdb=" SD  MET A   1 "
   model   vdw
   4.798 3.480""")

def exercise_consistency_between_getting_grm(mon_lib_srv, ener_lib):
  simple_grm, simple_xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records91_noh)
  simple_geo_str = show_sorted_geometry_str(simple_grm, simple_xrs)

  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  # params.pdb_interpretation.restraints_library.mcl=False
  # params.pdb_interpretation.use_neutron_distances=False
  model_grm, model_xrs = make_grm_via_model(mon_lib_srv, ener_lib, raw_records91_noh, params)
  model_geo_str = show_sorted_geometry_str(model_grm, model_xrs)
  # print(model_geo_str)
  assert not show_diff(simple_geo_str, model_geo_str)

def exercise_bond_in_symmetry_grm(mon_lib_srv, ener_lib):
  """ There is symmetry bond in GRM, but we removing non-symmetry one.

  Args:
      mon_lib_srv (_type_): _description_
      ener_lib (_type_): _description_
  """
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  grm, xrs = make_grm_via_model(mon_lib_srv, ener_lib, raw_records91_noh, params)
  simple, asu = grm.get_all_bond_proxies()
  assert (simple.size(), asu.size()) == (10, 2), (simple.size(), asu.size())
  sites_cart = xrs.sites_cart()
  initial_geo = show_sorted_geometry_str(grm, xrs)

  # Now we are removing the bonds one by one
  assert grm.is_bonded_atoms(0,1)
  grm.remove_bond_restraints_in_place(bonded_pairs=[(0,1)], sites_cart=sites_cart)
  assert not grm.is_bonded_atoms(0,1)
  simple, asu = grm.get_all_bond_proxies()
  assert (simple.size(), asu.size()) == (9, 2), (simple.size(), asu.size())

  final_geo = show_sorted_geometry_str(grm, xrs)
  # print(final_geo)
  # STOP()
  n_pos = final_geo.find("Nonbonded | unspecified | interactions: 38")
  assert n_pos>0

  diff_out = StringIO()
  show_diff(final_geo, initial_geo, out=diff_out)
  diff_gv = diff_out.getvalue()
  # print("*"*80)
  # print(diff_gv)
  # STOP()
  assert not show_diff(diff_gv,"""\
---
+++
@@ -1,4 +1,4 @@
-Bond | covalent geometry | restraints: 10
+Bond | covalent geometry | restraints: 9
 Sorted by residual:
 bond pdb=" CA  HIS B 304 "
      pdb=" C   HIS B 304 "
@@ -12,10 +12,6 @@
      pdb=" CE1 HIS B 304 "
   ideal  model  delta    sigma   weight residual
   1.321  1.324 -0.003 1.00e-02 1.00e+04 7.90e-02
-bond pdb=" N   HIS B 304 "
-     pdb=" CA  HIS B 304 "
-  ideal  model  delta    sigma   weight residual
-  1.458  1.463 -0.005 1.90e-02 2.77e+03 7.84e-02
 bond pdb=" CB  HIS B 304 "
      pdb=" CG  HIS B 304 "
   ideal  model  delta    sigma   weight residual
@@ -166,12 +162,24 @@
       pdb=" CE1 HIS B 304 "   -0.009 2.00e-02 2.50e+03
       pdb=" NE2 HIS B 304 "    0.020 2.00e-02 2.50e+03

-Nonbonded | unspecified | interactions: 35
+Nonbonded | unspecified | interactions: 38
 Sorted by model distance:
+nonbonded pdb=" N   HIS B 304 "
+          pdb=" CA  HIS B 304 "
+   model   vdw
+   1.463 3.550
 nonbonded pdb="ZN    ZN A   8 "
           pdb=" O   HOH A  57 "
    model   vdw
    2.154 2.230
+nonbonded pdb=" N   HIS B 304 "
+          pdb=" C   HIS B 304 "
+   model   vdw
+   2.465 3.350
+nonbonded pdb=" N   HIS B 304 "
+          pdb=" CB  HIS B 304 "
+   model   vdw
+   2.475 3.520
 nonbonded pdb=" CG  HIS B 304 "
           pdb="ZN    ZN A   8 "
    model   vdw sym.op.
@@ -179,7 +187,7 @@
 nonbonded pdb=" N   HIS B 304 "
           pdb=" CG  HIS B 304 "
    model   vdw
-   2.953 2.672
+   2.953 3.340
 nonbonded pdb=" O   HOH A  57 "
           pdb=" O   HOH A  69 "
    model   vdw
@@ -215,7 +223,7 @@
 nonbonded pdb=" N   HIS B 304 "
           pdb=" O   HIS B 304 "
    model   vdw
-   3.544 2.496
+   3.544 3.120
 nonbonded pdb=" CB  HIS B 304 "
           pdb=" CE1 HIS B 304 "
    model   vdw

""", strip_trailing_whitespace=True)


def exercise_bond_over_symmetry_2(mon_lib_srv, ener_lib):
  """ There is symmetry bond in GRM, but we remove it.
  Making sure that GRM build without it is identical to the one with it.
  Bond does not include angles.

  Args:
      mon_lib_srv (_type_): _description_
      ener_lib (_type_): _description_
  """
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  grm, xrs = make_grm_via_model(mon_lib_srv, ener_lib, raw_records91_noh, params)
  initial_geo_str = show_sorted_geometry_str(grm, xrs)
  simple, asu = grm.get_all_bond_proxies()
  assert (simple.size(), asu.size()) == (10, 2), (simple.size(), asu.size())
  # print('before', simple.size(), asu.size())
  sites_cart = xrs.sites_cart()

  # Now we are removing the bonds one by one
  assert grm.is_bonded_atoms(9,10)
  grm.remove_bond_restraints_in_place(bonded_pairs=[(9,10)], sites_cart=sites_cart)
  assert not grm.is_bonded_atoms(9,10)
  simple, asu = grm.get_all_bond_proxies()
  assert (simple.size(), asu.size()) == (10, 0), (simple.size(), asu.size())

  # print('after', simple.size(), asu.size())
  after_geo_str = show_sorted_geometry_str(grm, xrs)
  # print(after_geo_str)

  # Here we creating GRM without the Zn -- NE2 bond
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.restraints_library.mcl=False
  nobond_grm, xrs = make_grm_via_model(mon_lib_srv, ener_lib, raw_records91_noh, params)
  simple, asu = nobond_grm.get_all_bond_proxies()
  assert (simple.size(), asu.size()) == (10, 0), (simple.size(), asu.size())
  nobond_geo_str = show_sorted_geometry_str(nobond_grm, xrs)
  assert not show_diff(after_geo_str, nobond_geo_str)

def exercise_bond_with_self(mon_lib_srv, ener_lib):
  """Shell_sym_table has the following in this case:
i_seq: 12
  j_seq: 12
    -x-1,z,y
  """
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  grm, xrs = make_grm_via_model(mon_lib_srv, ener_lib, raw_records1, params)
  sites_cart = xrs.sites_cart()
  initial_geo_str = show_sorted_geometry_str(grm, xrs)
  # with open("exercise_bond_with_self_initial.geo", 'w') as f:
  #   f.write(initial_geo_str)
  # print(initial_geo_str)
  assert_lines_in_text(initial_geo_str,"""\
Bond | Disulphide bridge | restraints: 1
Sorted by residual:
bond pdb=" SG  CYS A   2 "
     pdb=" SG  CYS A   2 "
  ideal  model  delta    sigma   weight residual sym.op.
  2.031  2.028  0.003 2.00e-02 2.50e+03 2.48e-02 -x-1,z,y
""")
  simple, asu = grm.get_all_bond_proxies()
  assert (simple.size(), asu.size()) == (23, 1), (simple.size(), asu.size())
  assert grm.is_bonded_atoms(12,12)
  grm.remove_bond_restraints_in_place(bonded_pairs=[(12,12)], sites_cart=sites_cart)
  assert not grm.is_bonded_atoms(12,12)
  simple, asu = grm.get_all_bond_proxies()
  assert (simple.size(), asu.size()) == (23, 0), (simple.size(), asu.size())
  after_geo_str = show_sorted_geometry_str(grm, xrs)
  # with open("exercise_bond_with_self_after.geo", 'w') as f:
  #   f.write(after_geo_str)
  diff_out = StringIO()
  show_diff(after_geo_str, initial_geo_str, out=diff_out)
  diff_gv = diff_out.getvalue()
  # print("*"*80)
  # print(diff_gv)
  # STOP()
  assert_lines_in_text (diff_gv,"""\
-
-Bond | Disulphide bridge | restraints: 1
-Sorted by residual:
-bond pdb=" SG  CYS A   2 "
-     pdb=" SG  CYS A   2 "
-  ideal  model  delta    sigma   weight residual sym.op.
-  2.031  2.028  0.003 2.00e-02 2.50e+03 2.48e-02 -x-1,z,y
""")
  assert_lines_in_text(diff_gv, """\
-Nonbonded | unspecified | interactions: 80
+Nonbonded | unspecified | interactions: 83
 Sorted by model distance:
+nonbonded pdb=" SG  CYS A   2 "
+          pdb=" SG  CYS A   2 "
+   model   vdw sym.op.
+   2.028 3.760 -x-1,z,y
""")

  assert_lines_in_text(diff_gv, """\
+nonbonded pdb=" CB  CYS A   2 "
+          pdb=" SG  CYS A   2 "
+   model   vdw sym.op.
+   2.753 3.800 -x-1,z,y
""")
  assert_lines_in_text(diff_gv, """\
+nonbonded pdb=" SG  CYS A   2 "
+          pdb=" CB  CYS A   2 "
+   model   vdw sym.op.
+   2.753 3.800 -x-1,z,y
""")

def exercise_on_special_position(mon_lib_srv, ener_lib):
  """
  Looks like special position is not different from no symmetry case.
  """
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.automatic_linking.link_all=True
  grm, xrs = make_grm_via_model(mon_lib_srv, ener_lib, raw_records2, params)
  sites_cart = xrs.sites_cart()
  initial_geo_str = show_sorted_geometry_str(grm, xrs)
  # with open("exercise_bond_with_self_initial.geo", 'w') as f:
  #   f.write(initial_geo_str)
  # print(initial_geo_str)
  simple, asu = grm.get_all_bond_proxies()
  assert (simple.size(), asu.size()) == (12, 0), (simple.size(), asu.size())
  assert grm.is_bonded_atoms(9,10)
  assert grm.is_bonded_atoms(10,11)
  grm.remove_bond_restraints_in_place(bonded_pairs=[(9,10), (10,11)], sites_cart=sites_cart)
  assert not grm.is_bonded_atoms(9,10)
  assert not grm.is_bonded_atoms(10,11)
  after_geo_str = show_sorted_geometry_str(grm, xrs)
  # print(final_geo_str)
  diff_out = StringIO()
  show_diff(after_geo_str, initial_geo_str, out=diff_out)
  diff_gv = diff_out.getvalue()
  # print("*"*80)
  # print(diff_gv)
  # STOP()
  assert_lines_in_text(diff_gv,"""\
---
+++
@@ -40,17 +40,6 @@
      pdb=" CE1 HIS D  10 "
   ideal  model  delta    sigma   weight residual
   1.321  1.321  0.000 1.00e-02 1.00e+04 3.82e-04
-
-Bond | Metal coordination | restraints: 2
-Sorted by residual:
-bond pdb="CO    CO D 101 "
-     pdb=" O   HOH D 208 "
-  ideal  model  delta    sigma   weight residual
-  2.090  2.429 -0.339 1.00e-01 1.00e+02 1.15e+01
-bond pdb=" NE2 HIS D  10 "
-     pdb="CO    CO D 101 "
-  ideal  model  delta    sigma   weight residual
-  2.140  2.155 -0.015 1.00e-01 1.00e+02 2.21e-02

 Bond angle | covalent geometry | restraints: 12
 Sorted by residual:
@@ -166,7 +155,7 @@
       pdb=" CE1 HIS D  10 "   -0.005 2.00e-02 2.50e+03
       pdb=" NE2 HIS D  10 "    0.009 2.00e-02 2.50e+03

-Nonbonded | unspecified | interactions: 81
+Nonbonded | unspecified | interactions: 86
""")

def exercise_remove_ss(mon_lib_srv, ener_lib):
  file_path = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/10mh_helix_sheet_na.pdb",
    test=os.path.isfile)
  if (file_path is None):
    print("Skipping exercise_remove_ss: input file not available")
    return
  with open(file_path, 'r') as f:
    raw_records=f.readlines()

  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.secondary_structure.enabled=True
  params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair[0].restrain_planarity=True
  grm, xrs = make_grm_via_model(mon_lib_srv, ener_lib, raw_records, params)
  sites_cart = xrs.sites_cart()
  initial_geo_str = show_sorted_geometry_str(grm, xrs)
  # with open('start_exercise_remove_ss.geo', 'w') as f:
  #   f.write(initial_geo_str)
  assert_lines_in_text(initial_geo_str, "Bond | Bond-like | restraints: 90")
  assert_lines_in_text(initial_geo_str, "Bond angle | Secondary Structure restraints around h-bond | restraints: 190")
  assert_lines_in_text(initial_geo_str, "Parallelity | Stacking parallelity | restraints: 16")
  assert_lines_in_text(initial_geo_str, "Planarity | Basepair planarity | restraints: 10")

  grm.remove_secondary_structure_restraints(sites_cart=sites_cart)
  end_geo_str = show_sorted_geometry_str(grm, xrs)
  # with open('end_exercise_remove_ss.geo', 'w') as f:
  #   f.write(end_geo_str)
  diff_out = StringIO()
  show_diff(end_geo_str, initial_geo_str, out=diff_out)
  diff_gv = diff_out.getvalue()
  assert_lines_in_text(diff_gv, "-Bond | Bond-like | restraints: 90")
  assert_lines_in_text(diff_gv, "-Bond angle | Secondary Structure restraints around h-bond | restraints: 190")
  assert_lines_in_text(diff_gv, "-Parallelity | Stacking parallelity | restraints: 16")
  assert_lines_in_text(diff_gv, "-Planarity | Basepair planarity | restraints: 10")


def exercise():
  mon_lib_srv = None
  ener_lib = None
  try:
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
  except: # intentional
    print("Can not initialize monomer_library, skipping test.")
  if mon_lib_srv is not None and ener_lib is not None:
    exercise_consistency_between_getting_grm(mon_lib_srv, ener_lib)
    exercise_remove_bond_restraint_in_place(mon_lib_srv, ener_lib)
    exercise_remove_two_bond_restraints_in_place(mon_lib_srv, ener_lib)
    exercise_bond_in_symmetry_grm(mon_lib_srv, ener_lib)
    exercise_bond_over_symmetry_2(mon_lib_srv, ener_lib)
    exercise_bond_with_self(mon_lib_srv, ener_lib)
    exercise_on_special_position(mon_lib_srv, ener_lib)
    exercise_remove_ss(mon_lib_srv, ener_lib)

if (__name__ == "__main__"):
  exercise()
  print("OK")
