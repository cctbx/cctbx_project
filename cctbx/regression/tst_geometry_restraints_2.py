from __future__ import absolute_import, division, print_function
from cctbx import geometry_restraints
from cctbx.geometry_restraints.distance_least_squares \
  import distance_and_repulsion_least_squares
import cctbx.geometry_restraints.manager
from cctbx import crystal
from cctbx.array_family import flex
from six.moves import cStringIO as StringIO
import libtbx.utils
from libtbx.test_utils import approx_equal, show_diff
import libtbx.load_env
import sys, os
from mmtbx.monomer_library import pdb_interpretation
from six.moves import range

# ============================================================
# Edit notes: show_interactions function is obsolete and removed from GRM
# 5/22/2015
# ============================================================

def exercise_with_zeolite(verbose):
  if (not libtbx.env.has_module("iotbx")):
    print("Skipping exercise_with_zeolite(): iotbx not available")
    return
  from iotbx.kriber import strudat
  atlas_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/misc/strudat_zeolite_atlas",
    test=os.path.isfile)
  if (atlas_file is None):
    print("Skipping exercise_with_zeolite(): input file not available")
    return
  with open(atlas_file) as f:
    strudat_contents = strudat.read_all_entries(f)
  strudat_entry = strudat_contents.get("YUG")
  si_structure = strudat_entry.as_xray_structure()
  if (verbose):
    out = sys.stdout
  else:
    out = StringIO()
  drls = distance_and_repulsion_least_squares(
    si_structure=si_structure,
    distance_cutoff=3.5,
    nonbonded_repulsion_function_type="prolsq",
    n_macro_cycles=2,
    out=out)
  #
  nbp = drls.geometry_restraints_manager.pair_proxies().nonbonded_proxies
  assert nbp.n_total() > 50
    # expected is 60, but the exact number depends on the minimizer
  #
  site_labels = drls.minimized_structure.scatterers().extract_labels()
  sites_cart = drls.start_structure.sites_cart()
  pair_proxies = drls.geometry_restraints_manager.pair_proxies()
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    site_labels=site_labels,
    f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert len(out.getvalue().splitlines()) == 48*4+2
  assert out.getvalue().splitlines()[-1].find("remaining") < 0
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    site_labels=site_labels,
    f=out,
    prefix="0^",
    max_items=28)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue().replace("e-00", "e-0"), """\
0^Bond restraints: 48
0^Sorted by residual:
0^bond O3
0^     O4
0^  ideal  model  delta    sigma   weight residual
0^  2.629  2.120  0.509 1.56e+00 4.10e-01 1.06e-01
...
0^bond SI1
0^     SI1
0^  ideal  model  delta    sigma   weight residual sym.op.
0^  3.071  3.216 -0.145 2.08e+00 2.31e-01 4.83e-03 -x+1/2,-y+1/2,-z+1
0^... (remaining 20 not shown)
""",
    selections=[range(6), range(-5,0)])
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted(
    by_value="delta",
    sites_cart=sites_cart,
    site_labels=site_labels,
    f=out,
    prefix="0^",
    max_items=28)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue().replace("e-00", "e-0"), """\
0^Bond restraints: 48
0^Sorted by delta:
0^bond O3
0^     O4
0^  ideal  model  delta    sigma   weight residual
0^  2.629  2.120  0.509 1.56e+00 4.10e-01 1.06e-01
...
0^... (remaining 20 not shown)
""",
    selections=[range(6), [-1]])
  site_labels_long = ["abc"+label+"def" for label in site_labels]
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    site_labels=site_labels_long,
    f=out,
    prefix="^0",
    max_items=28)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue().replace("e-00", "e-0"), """\
^0Bond restraints: 48
^0Sorted by residual:
^0bond abcO3def
^0     abcO4def
^0  ideal  model  delta    sigma   weight residual
^0  2.629  2.120  0.509 1.56e+00 4.10e-01 1.06e-01
...
^0bond abcSI1def
^0     abcSI1def
^0  ideal  model  delta    sigma   weight residual sym.op.
^0  3.071  3.216 -0.145 2.08e+00 2.31e-01 4.83e-03 -x+1/2,-y+1/2,-z+1
^0... (remaining 20 not shown)
""",
    selections=[range(6), range(-5,0)])
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    f=out,
    prefix=".=",
    max_items=28)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue().replace("e-00", "e-0"), """\
.=Bond restraints: 48
.=Sorted by residual:
.=bond 4
.=     5
.=  ideal  model  delta    sigma   weight residual
.=  2.629  2.120  0.509 1.56e+00 4.10e-01 1.06e-01
...
.=bond 0
.=     0
.=  ideal  model  delta    sigma   weight residual sym.op.
.=  3.071  3.216 -0.145 2.08e+00 2.31e-01 4.83e-03 -x+1/2,-y+1/2,-z+1
.=... (remaining 20 not shown)
""",
    selections=[range(6), range(-5,0)])
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    f=out,
    prefix="-+",
    max_items=1)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue().replace("e-00", "e-0"), """\
-+Bond restraints: 48
-+Sorted by residual:
-+bond 4
-+     5
-+  ideal  model  delta    sigma   weight residual
-+  2.629  2.120  0.509 1.56e+00 4.10e-01 1.06e-01
-+... (remaining 47 not shown)
""")
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    f=out,
    prefix="=+",
    max_items=0)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
=+Bond restraints: 48
=+Sorted by residual:
=+... (remaining 48 not shown)
""")
  #
  sites_cart = si_structure.sites_cart()
  site_labels = [sc.label for sc in si_structure.scatterers()]
  asu_mappings = si_structure.asu_mappings(buffer_thickness=3.5)
  for min_cubicle_edge in [0, 5]:
    pair_generator = crystal.neighbors_fast_pair_generator(
      asu_mappings=asu_mappings,
      distance_cutoff=asu_mappings.buffer_thickness(),
      minimal=False,
      min_cubicle_edge=min_cubicle_edge)
    sorted_asu_proxies = geometry_restraints.nonbonded_sorted_asu_proxies(
      asu_mappings=asu_mappings)
    while (not pair_generator.at_end()):
      p = geometry_restraints.nonbonded_asu_proxy(
        pair=next(pair_generator),
        vdw_distance=3)
      sorted_asu_proxies.process(p)
    out = StringIO()
    sorted_asu_proxies.show_sorted(
      by_value="delta",
      sites_cart=sites_cart,
      site_labels=site_labels,
      f=out,
      prefix="d%")
    if (verbose):
      sys.stdout.write(out.getvalue())
    assert not show_diff(out.getvalue(), """\
d%Nonbonded interactions: 7
d%Sorted by model distance:
...
d%nonbonded SI2
d%          SI2
d%   model   vdw sym.op.
d%   3.092 3.000 -x+1,y,-z
...
d%nonbonded SI1
d%          SI1
d%   model   vdw sym.op.
d%   3.216 3.000 -x+1/2,-y+1/2,-z+1
""",
      selections=[range(2), range(10,14), range(26,30)])
    out = StringIO()
    sorted_asu_proxies.show_sorted(
      by_value="delta",
      sites_cart=sites_cart,
      f=out,
      prefix="*j",
      max_items=5)
    if (verbose):
      sys.stdout.write(out.getvalue())
    assert not show_diff(out.getvalue(), """\
*jNonbonded interactions: 7
*jSorted by model distance:
...
*jnonbonded 0
*j          1
*j   model   vdw
*j   3.107 3.000
*jnonbonded 0
*j          0
*j   model   vdw sym.op.
*j   3.130 3.000 -x+1,y,-z+1
*j... (remaining 2 not shown)
""",
      selections=[range(2), range(-9,0)])
    out = StringIO()
    sorted_asu_proxies.show_sorted(
      by_value="delta",
      sites_cart=sites_cart,
      f=out,
      prefix="@r",
      max_items=0)
    if (verbose):
      sys.stdout.write(out.getvalue())
    assert not show_diff(out.getvalue(), """\
@rNonbonded interactions: 7
""")

enk_pdb = """\
CRYST1   10.851   13.095   21.192  90.00  90.00  90.00 P 21 21 21
ATOM      1  CA  TYR A   1       8.787   2.175   5.487  1.00  0.91           C
ATOM      2  CB  TYR A   1       8.968   2.012   6.998  1.00  1.05           C
ATOM      3  CG  TYR A   1       9.527   0.669   7.410  1.00  1.04           C
ATOM      4  CD2 TYR A   1       8.768  -0.222   8.157  1.00  1.37           C
ATOM      5  CE2 TYR A   1       9.275  -1.449   8.537  1.00  1.50           C
ATOM      6  CZ  TYR A   1      10.559  -1.798   8.175  1.00  1.15           C
ATOM      7  CE1 TYR A   1      11.334  -0.929   7.437  1.00  1.37           C
ATOM      8  CD1 TYR A   1      10.818   0.296   7.062  1.00  1.34           C
ATOM      9  C   TYR A   1       7.767   1.172   4.959  1.00  0.85           C
ATOM     10  OH  TYR A   1      11.069  -3.019   8.553  1.00  1.40           O
ATOM     11  O   TYR A   1       6.576   1.472   4.870  1.00  1.01           O
ATOM     12  N   TYR A   1       8.388   3.541   5.166  1.00  1.11           N
ATOM     13  CA  GLY A   2       7.387  -1.050   4.044  1.00  1.14           C
ATOM     14  C   GLY A   2       6.345  -1.569   5.016  1.00  0.97           C
ATOM     15  O   GLY A   2       5.234  -1.918   4.619  1.00  1.15           O
ATOM     16  N   GLY A   2       8.241  -0.020   4.608  1.00  0.93           N
ATOM     17  CA  GLY A   3       5.804  -2.100   7.324  1.00  1.36           C
ATOM     18  C   GLY A   3       4.651  -1.149   7.578  1.00  1.01           C
ATOM     19  O   GLY A   3       3.598  -1.553   8.071  1.00  1.38           O
ATOM     20  N   GLY A   3       6.706  -1.622   6.294  1.00  1.11           N
ATOM     21  CA  PHE A   4       3.819   1.134   7.419  1.00  0.89           C
ATOM     22  CB  PHE A   4       4.397   2.380   8.094  1.00  1.13           C
ATOM     23  CG  PHE A   4       4.930   2.130   9.475  1.00  1.00           C
ATOM     24  CD1 PHE A   4       6.267   1.825   9.673  1.00  1.51           C
ATOM     25  CE1 PHE A   4       6.760   1.595  10.943  1.00  1.80           C
ATOM     26  CZ  PHE A   4       5.916   1.667  12.033  1.00  1.59           C
ATOM     27  CE2 PHE A   4       4.582   1.970  11.850  1.00  1.49           C
ATOM     28  CD2 PHE A   4       4.095   2.199  10.577  1.00  1.24           C
ATOM     29  C   PHE A   4       3.185   1.509   6.084  1.00  0.94           C
ATOM     30  N   PHE A   4       4.852   0.121   7.242  1.00  0.88           N
ATOM     31  O   PHE A   4       2.361   2.421   6.010  1.00  1.47           O
ATOM     32  CA  LEU A   5       3.055   1.059   3.693  1.00  0.87           C
ATOM     33  CB  LEU A   5       3.965   0.435   2.634  1.00  1.13           C
ATOM     34  CG  LEU A   5       3.531   0.603   1.177  1.00  1.16           C
ATOM     35  CD1 LEU A   5       3.411   2.076   0.818  1.00  1.88           C
ATOM     36  CD2 LEU A   5       4.502  -0.103   0.245  1.00  1.67           C
ATOM     37  C   LEU A   5       1.634   0.527   3.541  1.00  0.87           C
ATOM     38  N   LEU A   5       3.576   0.800   5.030  1.00  0.92           N
ATOM     39  OXT LEU A   5       1.246  -0.440   4.196  1.00  1.23           O
ATOM     41  O   HOH B   1      14.655  -4.248   8.995  1.00  1.51           O
ATOM     81  H1  HOH B   1      14.954  -3.400   8.664  1.00  1.13           H
ATOM     82  H2  HOH B   1      13.712  -4.246   8.805  1.00  0.76           H
ATOM     42  O   HOH B   2      12.055  -3.540   8.243  1.00  1.81           O
ATOM     83  H1  HOH B   2      11.456  -4.167   7.841  1.00  0.96           H
ATOM     84  H2  HOH B   2      11.476  -3.010   8.803  1.00  1.17           H
ATOM     43  O   HOH B   3       9.622  -2.103   5.551  1.00  2.24           O
ATOM     85  H1  HOH B   3      10.394  -1.650   5.193  1.00  1.02           H
ATOM     86  H2  HOH B   3       9.588  -2.843   4.937  1.00  1.29           H
END
"""

def exercise_with_pdb(verbose):
  if (not libtbx.env.has_module(name="mmtbx")):
    print("Skipping exercise_with_pdb():", \
      "mmtbx.monomer_library.pdb_interpretation not available")
    return
  if (libtbx.env.find_in_repositories(relative_path="chem_data") is None):
    print("Skipping exercise_with_pdb(): chem_data directory not available")
    return
  if (verbose):
    out = sys.stdout
  else:
    out = StringIO()
  with open("tmp_cctbx_geometry_restraints.pdb", "w") as f:
    f.write(enk_pdb)
  pdb_interpretation_params = pdb_interpretation.master_params.extract()
  pdb_interpretation_params.sort_atoms=False
  processed_pdb_file = pdb_interpretation.run(
    args=["tmp_cctbx_geometry_restraints.pdb"],
    strict_conflict_handling=False,
    params=pdb_interpretation_params,
    log=out,)
  geo = processed_pdb_file.geometry_restraints_manager()
  site_labels = processed_pdb_file.xray_structure().scatterers() \
    .extract_labels()
  #
  assert approx_equal(flex.min(geo.nonbonded_model_distances()), 0.4777342)
  #
  geo._sites_cart_used_for_pair_proxies = None
  #
  sel0 = geo.simple_edge_list()
  assert len(sel0) == 46
  assert sel0[:4] == [(0, 1), (0, 8), (0, 11), (1, 2)]
  assert sel0[-4:] == [(42, 43), (42, 44), (45, 46), (45, 47)]
  geo.bond_params_table[13][14].slack = 0.1
  geo.bond_params_table[28][30].slack = 0.3
  sel = geo.simple_edge_list()
  assert sorted(set(sel0) - set(sel)) == [(13, 14), (28, 30)]
  sel = geo.simple_edge_list(omit_slack_greater_than=0.2)
  assert sorted(set(sel0) - set(sel)) == [(28, 30)]
  #
  d = geo.discard_symmetry(new_unit_cell=(10,10,10,90,90,90))
  assert d.site_symmetry_table.special_position_indices().size()==0
  #
  clusters = geo.rigid_clusters_due_to_dihedrals_and_planes(
    constrain_dihedrals_with_sigma_less_than=10)
  assert sorted([tuple(sorted(c)) for c in clusters]) == [
    (0, 8, 10, 15), (0, 8, 12, 15), (1, 2, 3, 4, 5, 6, 7, 9),
    (5, 6, 7, 9), (12, 13, 14, 19), (12, 13, 16, 19), (16, 17, 18, 29),
    (16, 17, 20, 29), (20, 28, 30, 37), (20, 28, 31, 37),
    (21, 22, 23, 24, 25, 26, 27)]

def exercise_non_crystallographic_conserving_bonds_and_angles():
  sites_cart, geo = geometry_restraints.manager \
    .construct_non_crystallographic_conserving_bonds_and_angles(
      sites_cart=flex.vec3_double([
        (10.949, 12.815, 15.189),
        (10.405, 13.954, 15.917),
        (10.779, 15.262, 15.227),
        ( 9.916, 16.090, 14.936)]),
      edge_list_bonds=[(0, 1), (1, 2), (2, 3)],
      edge_list_angles=[(0, 2), (1, 3)])
  assert approx_equal(sites_cart, [
    (6.033, 5.000, 5.253),
    (5.489, 6.139, 5.981),
    (5.863, 7.447, 5.291),
    (5.000, 8.275, 5.000)])
  assert approx_equal(geo.energies_sites(sites_cart=sites_cart).target, 0)
  sites_cart_noise = flex.vec3_double([ # Just to make all residuals unique,
    (6.043, 5.030, 5.233),              # so that the sorted bond list below
    (5.469, 6.119, 5.941),              # has the same order on all platforms.
    (5.893, 7.487, 5.281),
    (5.040, 8.225, 5.020)])
  sio = StringIO()
  geo.show_sorted(sites_cart=sites_cart_noise, f=sio)
  expected_first_part = """\
Bond | covalent geometry | restraints: 5
Sorted by residual:
bond 2
     3
  ideal  model  delta    sigma   weight residual
  1.231  1.158  0.073 1.00e-01 1.00e+02 5.35e-01
bond 1
     2
  ideal  model  delta    sigma   weight residual
  1.525  1.577 -0.052 1.00e-01 1.00e+02 2.66e-01
bond 1
     3
  ideal  model  delta    sigma   weight residual
  2.401  2.338  0.063 1.41e-01 5.00e+01 1.96e-01
bond 0
     1
  ideal  model  delta    sigma   weight residual
  1.457  1.420  0.037 1.00e-01 1.00e+02 1.37e-01
bond 0
     2
  ideal  model  delta    sigma   weight residual
  2.453  2.462 -0.009 1.41e-01 5.00e+01 3.92e-03

"""
  assert not show_diff(sio.getvalue(), expected_first_part + """\
Nonbonded | unspecified | interactions: 0

""")
  #
  sites_cart, geo = geometry_restraints.manager \
    .construct_non_crystallographic_conserving_bonds_and_angles(
      sites_cart=flex.vec3_double([
        (10.949, 12.815, 15.189),
        (10.405, 13.954, 15.917),
        (10.779, 15.262, 15.227),
        ( 9.916, 16.090, 14.936),
        (10.749, 12.615, 15.389)]),
      edge_list_bonds=[(0, 1), (1, 2), (2, 3)],
      edge_list_angles=[(0, 2), (1, 3)])
  sites_cart_noise.append(sites_cart[-1])
  sio = StringIO()
  geo.show_sorted(sites_cart=sites_cart_noise, f=sio)
  assert not show_diff(sio.getvalue(), expected_first_part + """\
Nonbonded | unspecified | interactions: 2
Sorted by model distance:
nonbonded 0
          4
   model   vdw
   0.306 1.200
nonbonded 1
          4
   model   vdw
   1.274 1.200

""")

def exercise_na_restraints_output_to_geo(verbose=False):
  for dependency in ("chem_data", "ksdssp"):
    if not libtbx.env.has_module(dependency):
      print("Skipping exercise_na_restraints_output_to_geo(): %s not available" %(
        dependency))
      return
  pdb_str_1dpl_cutted="""\
CRYST1   24.627   42.717   46.906  90.00  90.00  90.00 P 21 21 21    8
ATOM    184  P    DG A   9       9.587  13.026  19.037  1.00  6.28           P
ATOM    185  OP1  DG A   9       9.944  14.347  19.602  1.00  8.07           O
ATOM    186  OP2  DG A   9      10.654  12.085  18.639  1.00  8.27           O
ATOM    187  O5'  DG A   9       8.717  12.191  20.048  1.00  5.88           O
ATOM    188  C5'  DG A   9       7.723  12.833  20.854  1.00  5.45           C
ATOM    189  C4'  DG A   9       7.145  11.818  21.807  1.00  5.40           C
ATOM    190  O4'  DG A   9       6.435  10.777  21.087  1.00  5.77           O
ATOM    191  C3'  DG A   9       8.142  11.036  22.648  1.00  5.10           C
ATOM    192  O3'  DG A   9       8.612  11.838  23.723  1.00  5.90           O
ATOM    193  C2'  DG A   9       7.300   9.857  23.068  1.00  5.97           C
ATOM    194  C1'  DG A   9       6.619   9.536  21.805  1.00  5.97           C
ATOM    195  N9   DG A   9       7.390   8.643  20.931  1.00  5.97           N
ATOM    196  C8   DG A   9       8.074   8.881  19.775  1.00  6.62           C
ATOM    197  N7   DG A   9       8.647   7.820  19.249  1.00  6.57           N
ATOM    198  C5   DG A   9       8.308   6.806  20.141  1.00  6.22           C
ATOM    199  C6   DG A   9       8.620   5.431  20.136  1.00  6.03           C
ATOM    200  O6   DG A   9       9.297   4.803  19.296  1.00  7.21           O
ATOM    201  N1   DG A   9       8.101   4.773  21.247  1.00  6.10           N
ATOM    202  C2   DG A   9       7.365   5.351  22.260  1.00  6.24           C
ATOM    203  N2   DG A   9       6.948   4.569  23.241  1.00  7.88           N
ATOM    204  N3   DG A   9       7.051   6.652  22.257  1.00  6.53           N
ATOM    205  C4   DG A   9       7.539   7.295  21.184  1.00  5.69           C
ATOM    206  P    DC A  10      10.081  11.538  24.300  1.00  5.91           P
ATOM    207  OP1  DC A  10      10.273  12.645  25.291  1.00  7.27           O
ATOM    208  OP2  DC A  10      11.063  11.363  23.228  1.00  6.84           O
ATOM    209  O5'  DC A  10       9.953  10.128  25.026  1.00  5.75           O
ATOM    210  C5'  DC A  10       9.077   9.959  26.149  1.00  5.87           C
ATOM    211  C4'  DC A  10       9.188   8.549  26.672  1.00  5.56           C
ATOM    212  O4'  DC A  10       8.708   7.612  25.667  1.00  5.70           O
ATOM    213  C3'  DC A  10      10.580   8.059  27.007  1.00  5.27           C
ATOM    214  O3'  DC A  10      11.010   8.447  28.315  1.00  5.83           O
ATOM    215  C2'  DC A  10      10.422   6.549  26.893  1.00  5.34           C
ATOM    216  C1'  DC A  10       9.436   6.405  25.754  1.00  5.23           C
ATOM    217  N1   DC A  10      10.113   6.168  24.448  1.00  5.30           N
ATOM    218  C2   DC A  10      10.514   4.871  24.152  1.00  5.28           C
ATOM    219  O2   DC A  10      10.283   3.972  25.000  1.00  5.75           O
ATOM    220  N3   DC A  10      11.131   4.627  22.965  1.00  5.65           N
ATOM    221  C4   DC A  10      11.395   5.628  22.138  1.00  5.80           C
ATOM    222  N4   DC A  10      12.034   5.327  21.005  1.00  6.75           N
ATOM    223  C5   DC A  10      11.029   6.970  22.449  1.00  5.99           C
ATOM    224  C6   DC A  10      10.394   7.203  23.612  1.00  5.56           C
ATOM    226  O5'  DG B  11      12.424  -4.393  18.427  1.00 22.70           O
ATOM    227  C5'  DG B  11      12.380  -5.516  19.282  1.00 14.75           C
ATOM    228  C4'  DG B  11      11.969  -5.112  20.676  1.00 10.42           C
ATOM    229  O4'  DG B  11      12.972  -4.192  21.210  1.00 10.51           O
ATOM    230  C3'  DG B  11      10.649  -4.394  20.782  1.00  8.57           C
ATOM    231  O3'  DG B  11       9.618  -5.363  20.846  1.00  8.69           O
ATOM    232  C2'  DG B  11      10.822  -3.597  22.051  1.00  8.63           C
ATOM    233  C1'  DG B  11      12.236  -3.233  21.980  1.00  9.81           C
ATOM    234  N9   DG B  11      12.509  -1.902  21.305  1.00  8.66           N
ATOM    235  C8   DG B  11      13.175  -1.667  20.135  1.00  9.57           C
ATOM    236  N7   DG B  11      13.255  -0.407  19.824  1.00  9.04           N
ATOM    237  C5   DG B  11      12.613   0.235  20.869  1.00  7.63           C
ATOM    238  C6   DG B  11      12.388   1.612  21.119  1.00  7.05           C
ATOM    239  O6   DG B  11      12.723   2.590  20.419  1.00  7.81           O
ATOM    240  N1   DG B  11      11.715   1.819  22.317  1.00  6.27           N
ATOM    241  C2   DG B  11      11.264   0.828  23.159  1.00  6.05           C
ATOM    242  N2   DG B  11      10.611   1.219  24.248  1.00  5.85           N
ATOM    243  N3   DG B  11      11.483  -0.457  22.942  1.00  6.55           N
ATOM    244  C4   DG B  11      12.150  -0.687  21.797  1.00  6.84           C
ATOM    245  P    DC B  12       8.134  -5.009  20.350  1.00  8.13           P
ATOM    246  OP1  DC B  12       7.367  -6.252  20.459  1.00 10.02           O
ATOM    247  OP2  DC B  12       8.172  -4.307  19.052  1.00  9.79           O
ATOM    248  O5'  DC B  12       7.564  -3.912  21.389  1.00  8.18           O
ATOM    249  C5'  DC B  12       7.275  -4.296  22.719  1.00  8.00           C
ATOM    250  C4'  DC B  12       6.856  -3.057  23.487  1.00  8.01           C
ATOM    251  O4'  DC B  12       8.006  -2.146  23.615  1.00  7.35           O
ATOM    252  C3'  DC B  12       5.763  -2.208  22.890  1.00  7.04           C
ATOM    253  O3'  DC B  12       4.456  -2.800  23.100  1.00  9.82           O
ATOM    254  C2'  DC B  12       6.019  -0.916  23.630  1.00  6.50           C
ATOM    255  C1'  DC B  12       7.467  -0.808  23.608  1.00  7.35           C
ATOM    256  N1   DC B  12       8.040  -0.143  22.396  1.00  6.64           N
ATOM    257  C2   DC B  12       8.017   1.257  22.382  1.00  5.68           C
ATOM    258  O2   DC B  12       7.524   1.832  23.357  1.00  6.32           O
ATOM    259  N3   DC B  12       8.543   1.930  21.312  1.00  6.18           N
ATOM    260  C4   DC B  12       9.009   1.236  20.266  1.00  6.48           C
ATOM    261  N4   DC B  12       9.518   1.926  19.243  1.00  7.43           N
ATOM    262  C5   DC B  12       9.012  -0.198  20.248  1.00  6.83           C
ATOM    263  C6   DC B  12       8.502  -0.825  21.311  1.00  6.80           C
  """
  identical_portions = [
  """\
  Histogram of bond lengths:
        1.23 -     1.31: 5
        1.31 -     1.39: 25
        1.39 -     1.46: 27
        1.46 -     1.54: 25
        1.54 -     1.61: 5
  Bond restraints: 87""",
  '''\
  Histogram of bond angle deviations from ideal:
        0.00 -     1.10: 75
        1.10 -     2.19: 30
        2.19 -     3.29: 13
        3.29 -     4.38: 8
        4.38 -     5.48: 4
  Bond angle restraints: 130''',
  ]
  with open("tst_cctbx_geometry_restraints_2_na.pdb", "w") as f:
    f.write(pdb_str_1dpl_cutted)
  out1 = StringIO()
  out2 = StringIO()
  from mmtbx.monomer_library.server import MonomerLibraryServerError
  try:
    processed_pdb_file = pdb_interpretation.run(
      args=["tst_cctbx_geometry_restraints_2_na.pdb"],
      strict_conflict_handling=False,
      log=out1)
  except MonomerLibraryServerError:
    print("Skipping exercise_na_restraints_output_to_geo(): Encountered MonomerLibraryServerError.\n")
    print("Is the CCP4 monomer library installed and made available through environment variables MMTBX_CCP4_MONOMER_LIB or CLIBD_MON?")
    return
  geo1 = processed_pdb_file.geometry_restraints_manager()
  hbp = geo1.get_n_hbond_proxies()
  from mmtbx import monomer_library
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.secondary_structure.enabled=True
  processed_pdb_file = pdb_interpretation.run(
    args=["tst_cctbx_geometry_restraints_2_na.pdb"],
    params=params,
    strict_conflict_handling=False,
    log=out2)
  geo2 = processed_pdb_file.geometry_restraints_manager()
  hbp = geo2.get_n_hbond_proxies()
  v_out1 = out1.getvalue()
  v_out2 = out2.getvalue()
  assert v_out2.find("""\
    Restraints generated for nucleic acids:
      6 hydrogen bonds
      12 hydrogen bond angles
      0 basepair planarities
      2 basepair parallelities
      2 stacking parallelities""") > 0
  for v in [v_out1, v_out2]:
    for portion in identical_portions:
      if not v.find(portion) > 0:
        print("This portion was not found:\n%s\n=====End of portion." % portion)
        assert 0, "the portion above does not match expected portion."
  # check .geo output
  geo_identical_portions = ["Bond | covalent geometry | restraints: 87",
      "Bond angle | covalent geometry | restraints: 130",
      "Dihedral angle | covalent geometry | restraints: 33",
      "Chirality | covalent geometry | restraints: 15",
      "Planarity | covalent geometry | restraints: 4"]
  ss_geo_portions = [
      # "Bond-like restraints: 6",
      "Bond | Bond-like | restraints: 6",
      # 'Secondary Structure restraints around h-bond angle restraints: 12',
      'Bond angle | Secondary Structure restraints around h-bond | restraints: 12',
      # "Stacking parallelity restraints: 2",
      # 'Basepair parallelity restraints: 2',
      'Parallelity | Stacking parallelity | restraints: 2',
      'Parallelity | Basepair parallelity | restraints: 2',
      "Nonbonded | unspecified | interactions: 504"]
  non_ss_geo_portions = [
      #"Bond-like restraints: 0",
      #'Secondary Structure restraints around h-bond angle restraints: 0',
      # "Parallelity restraints: 0", removed because zero
      "Nonbonded | unspecified | interactions: 526"]
  acp = processed_pdb_file.all_chain_proxies
  sites_cart = acp.sites_cart_exact()
  site_labels = [atom.id_str() for atom in acp.pdb_atoms]
  geo_out1 = StringIO()
  geo_out2 = StringIO()
  geo1.show_sorted(sites_cart=sites_cart, site_labels=site_labels, f=geo_out1)
  geo2.show_sorted(sites_cart=sites_cart, site_labels=site_labels, f=geo_out2)
  v_geo_out_noss = geo_out1.getvalue()
  v_geo_out_ss = geo_out2.getvalue()
  for portion in geo_identical_portions+ss_geo_portions:
    assert v_geo_out_ss.find(portion) >= 0, 'did not find %s\n in \n%s\n%s\n' % (
      portion,
      v_geo_out_ss,
      portion)
  for portion in geo_identical_portions+non_ss_geo_portions:
    assert v_geo_out_noss.find(portion) >= 0, 'did not find %s\n in \n%s' % (
      portion,
      v_geo_out_noss)

def exercise_hydrogen_bond_rmsd():
  for dependency in ("chem_data", "ksdssp"):
    if not libtbx.env.has_module(dependency):
      print("Skipping exercise_na_restraints_output_to_geo(): %s not available" %(
        dependency))
      return
  if not os.path.exists('tst_cctbx_geometry_restraints_2_na.pdb'):
    print('exercise_hydrogen_bond_rmsd SKIPPED')
  from libtbx import easy_run
  cmd = 'phenix.geometry_minimization tst_cctbx_geometry_restraints_2_na.pdb'
  print(cmd)
  rc=easy_run.go(cmd)
  for token in ['  87  Z', ' 130  Z']:
    for line in rc.stdout_lines:
      if line.find(token)>-1: break
    else: assert 0, 'token "%s" not found' % token
  cmd = 'phenix.geometry_minimization tst_cctbx_geometry_restraints_2_na.pdb'
  cmd += ' secondary_structure.enabled=True'
  print(cmd)
  rc=easy_run.go(cmd)
  for token in ['  87  Z', ' 130  Z']:
    for line in rc.stdout_lines:
      if line.find(token)>-1: break
    else: assert 0, 'token "%s" not found' % token

def exercise_all(args):
  verbose = "--verbose" in args
  exercise_with_zeolite(verbose=verbose)
  exercise_with_pdb(verbose=verbose)
  exercise_non_crystallographic_conserving_bonds_and_angles()
  exercise_na_restraints_output_to_geo(verbose=verbose)
  exercise_hydrogen_bond_rmsd()
  print(libtbx.utils.format_cpu_times())

if (__name__ == "__main__"):
  exercise_all(sys.argv[1:])
