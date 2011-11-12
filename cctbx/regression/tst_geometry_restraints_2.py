from cctbx import geometry_restraints
from cctbx.geometry_restraints.distance_least_squares \
  import distance_and_repulsion_least_squares
import cctbx.geometry_restraints.manager
from cctbx import crystal
from cctbx.array_family import flex
from cStringIO import StringIO
import libtbx.utils
from libtbx.test_utils import approx_equal, show_diff, blocks_show_diff
import libtbx.load_env
import sys, os

def exercise_with_zeolite(verbose):
  if (not libtbx.env.has_module("iotbx")):
    print "Skipping exercise_with_zeolite(): iotbx not available"
    return
  from iotbx.kriber import strudat
  atlas_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/misc/strudat_zeolite_atlas",
    test=os.path.isfile)
  if (atlas_file is None):
    print "Skipping exercise_with_zeolite(): input file not available"
    return
  strudat_contents = strudat.read_all_entries(open(atlas_file))
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
  out = StringIO()
  drls.geometry_restraints_manager.show_interactions(f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
bond simple: (0, 1)
  distance_model: 3.22437
  distance_ideal: 3.07097
  weight: 0.2308
...
bond asu: (0, 0) -x+1,y,-z+1
  distance_model: 3.4233
  distance_ideal: 3.07097
  weight: 0.2308
""",
    selections=[range(4), range(20,24)])
  nbp = drls.geometry_restraints_manager.pair_proxies().nonbonded_proxies
  assert nbp.n_total() > 50
    # expected is 60, but the exact number depends on the minimizer
  #
  site_labels = drls.minimized_structure.scatterers().extract_labels()
  out = StringIO()
  drls.geometry_restraints_manager.show_interactions(
    site_labels=site_labels, f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
bond simple: (0, 1)
  SI1
  SI2
  distance_model: 3.22437
  distance_ideal: 3.07097
  weight: 0.2308
...
bond asu: (0, 0)
  SI1
  SI1 -x+1,y,-z+1
  distance_model: 3.4233
  distance_ideal: 3.07097
  weight: 0.2308
""",
    selections=[range(6), range(30,36)])
  #
  out = StringIO()
  drls.geometry_restraints_manager._sites_cart_used_for_pair_proxies = None
  drls.geometry_restraints_manager.show_interactions(f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
bond simple: (0, 1)
  distance_ideal: 3.07097
  weight: 0.2308
...
bond asu: (0, 0) -x+1,y,-z+1
  distance_ideal: 3.07097
  weight: 0.2308
""",
    selections=[range(3), range(15,18)])
  #
  site_labels = drls.minimized_structure.scatterers().extract_labels()
  out = StringIO()
  drls.geometry_restraints_manager.show_interactions(
    site_labels=site_labels, f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
bond simple: (0, 1)
  SI1
  SI2
  distance_ideal: 3.07097
  weight: 0.2308
...
bond asu: (0, 0)
  SI1
  SI1 -x+1,y,-z+1
  distance_ideal: 3.07097
  weight: 0.2308
""",
    selections=[range(5), range(25,30)])
  #
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
        pair=pair_generator.next(),
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
    print "Skipping exercise_with_pdb():", \
      "mmtbx.monomer_library.pdb_interpretation not available"
    return
  if (libtbx.env.find_in_repositories(relative_path="chem_data") is None):
    print "Skipping exercise_with_pdb(): chem_data directory not available"
    return
  from mmtbx.monomer_library import pdb_interpretation
  if (verbose):
    out = sys.stdout
  else:
    out = StringIO()
  open("tmp.pdb", "w").write(enk_pdb)
  processed_pdb_file = pdb_interpretation.run(
    args=["tmp.pdb"],
    strict_conflict_handling=False,
    log=out)
  geo = processed_pdb_file.geometry_restraints_manager()
  site_labels = processed_pdb_file.xray_structure().scatterers() \
    .extract_labels()
  #
  out = StringIO()
  geo.show_interactions(f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not blocks_show_diff(out.getvalue(), """\
bond simple: (0, 1)
  distance_model: 1.53051
  distance_ideal: 1.53
  weight: 2500
...
angle: (1, 0, 8)
  angle_model: 110.543
  angle_ideal: 110.1
  weight: 0.277008
...
dihedral: (2, 0, 1, 11)
  angle_model: 173.684
  angle_ideal: -180
  weight: 0.00444444
  periodicity: 3
...
chirality: (0, 1, 8, 11)
  volume_model: -2.50215
  volume_ideal: -2.50318
  both_signs: 0
  weight: 25
...
planarity: (0, 8, 10, 15)
  delta: -0.00012, weight: 2500
  delta:  0.00042, weight: 2500
  delta: -0.00016, weight: 2500
  delta: -0.00014, weight: 2500
...
nonbonded asu: (36, 46) -x+1,y+1/2,-z+1/2
  distance_model: 4.89425
  vdw_distance: 2.95
nonbonded simple: (0, 23)
  distance_model: 4.89852
  vdw_distance: 3.77
""")
  #
  out = StringIO()
  geo.show_interactions(site_labels=site_labels, f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not blocks_show_diff(out.getvalue(), """\
bond simple: (0, 1)
  pdb=" CA  TYR A   1 "
  pdb=" CB  TYR A   1 "
  distance_model: 1.53051
  distance_ideal: 1.53
  weight: 2500
bond simple: (0, 8)
...
angle: (1, 0, 8)
  pdb=" CB  TYR A   1 "
  pdb=" CA  TYR A   1 "
  pdb=" C   TYR A   1 "
  angle_model: 110.543
  angle_ideal: 110.1
  weight: 0.277008
...
dihedral: (2, 0, 1, 11)
  pdb=" CG  TYR A   1 "
  pdb=" CA  TYR A   1 "
  pdb=" CB  TYR A   1 "
  pdb=" N   TYR A   1 "
  angle_model: 173.684
  angle_ideal: -180
  weight: 0.00444444
  periodicity: 3
...
chirality: (0, 1, 8, 11)
  pdb=" CA  TYR A   1 "
  pdb=" CB  TYR A   1 "
  pdb=" C   TYR A   1 "
  pdb=" N   TYR A   1 "
  volume_model: -2.50215
  volume_ideal: -2.50318
  both_signs: 0
  weight: 25
...
planarity: (0, 8, 10, 15)
  pdb=" CA  TYR A   1 " delta: -0.00012, weight: 2500
  pdb=" C   TYR A   1 " delta:  0.00042, weight: 2500
  pdb=" O   TYR A   1 " delta: -0.00016, weight: 2500
  pdb=" N   GLY A   2 " delta: -0.00014, weight: 2500
""")
  #
  assert approx_equal(flex.min(geo.nonbonded_model_distances()), 0.4777342)
  #
  geo._sites_cart_used_for_pair_proxies = None
  out = StringIO()
  geo.show_interactions(f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not blocks_show_diff(out.getvalue(), """\
bond simple: (0, 1)
  distance_ideal: 1.53
  weight: 2500
...
angle: (1, 0, 8)
  angle_ideal: 110.1
  weight: 0.277008
...
dihedral: (2, 0, 1, 11)
  angle_ideal: -180
  weight: 0.00444444
  periodicity: 3
...
chirality: (0, 1, 8, 11)
  volume_ideal: -2.50318
  both_signs: 0
  weight: 25
...
planarity: (0, 8, 10, 15)
  weight: 2500
...
nonbonded simple: (5, 39)
  vdw_distance: 3.26
...
nonbonded asu: (7, 29) x+1,y,z
  vdw_distance: 3.42
""")
  #
  out = StringIO()
  geo.show_interactions(site_labels=site_labels, f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not blocks_show_diff(out.getvalue(), """\
bond simple: (0, 1)
  pdb=" CA  TYR A   1 "
  pdb=" CB  TYR A   1 "
  distance_ideal: 1.53
  weight: 2500
...
angle: (1, 0, 8)
  pdb=" CB  TYR A   1 "
  pdb=" CA  TYR A   1 "
  pdb=" C   TYR A   1 "
  angle_ideal: 110.1
  weight: 0.277008
...
dihedral: (2, 0, 1, 11)
  pdb=" CG  TYR A   1 "
  pdb=" CA  TYR A   1 "
  pdb=" CB  TYR A   1 "
  pdb=" N   TYR A   1 "
  angle_ideal: -180
  weight: 0.00444444
  periodicity: 3
...
chirality: (0, 1, 8, 11)
  pdb=" CA  TYR A   1 "
  pdb=" CB  TYR A   1 "
  pdb=" C   TYR A   1 "
  pdb=" N   TYR A   1 "
  volume_ideal: -2.50318
  both_signs: 0
  weight: 25
...
planarity: (0, 8, 10, 15)
  pdb=" CA  TYR A   1 " weight: 2500
  pdb=" C   TYR A   1 " weight: 2500
  pdb=" O   TYR A   1 " weight: 2500
  pdb=" N   GLY A   2 " weight: 2500
...
nonbonded simple: (5, 39)
  pdb=" CZ  TYR A   1 "
  pdb=" O   HOH B   1 "
  vdw_distance: 3.26
...
nonbonded asu: (7, 29)
  pdb=" CD1 TYR A   1 "
  pdb=" N   PHE A   4 " x+1,y,z
  vdw_distance: 3.42
""")
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
Bond restraints: 5
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
Nonbonded interactions: 0

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
Nonbonded interactions: 2
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

def exercise_all(args):
  verbose = "--verbose" in args
  exercise_with_zeolite(verbose=verbose)
  exercise_with_pdb(verbose=verbose)
  exercise_non_crystallographic_conserving_bonds_and_angles()
  print libtbx.utils.format_cpu_times()

if (__name__ == "__main__"):
  exercise_all(sys.argv[1:])
