from cctbx.geometry_restraints.distance_least_squares \
  import distance_and_repulsion_least_squares
from cctbx.array_family import flex
from iotbx.kriber import strudat
from cStringIO import StringIO
import libtbx.utils
from libtbx.test_utils import approx_equal, show_diff
import libtbx.load_env
import sys, os

def exercise_with_zeolite(verbose):
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
    out=out)
  #
  out = StringIO()
  drls.geometry_restraints_manager.show_interactions(f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue().replace("(2, 1)","(1, 2)"), """\
bond simple: (0, 1)
  distance_model: 3.22437
  distance_ideal: 3.07097
  weight: 0.2308
...
bond asu: (0, 0) -x+1,y,-z+1
  distance_model: 3.4233
  distance_ideal: 3.07097
  weight: 0.2308
...
nonbonded simple: (1, 5)
  distance_model: 3.74919
  vdw_distance: 1
...
nonbonded asu: (1, 2) x,y,z
  distance_model: 4.08434
  vdw_distance: 1
""",
    selections=[range(4), range(20,24), range(264,267), range(-3,0)])
  #
  site_labels = drls.minimized_structure.scatterers().extract_labels()
  out = StringIO()
  drls.geometry_restraints_manager.show_interactions(
    site_labels=site_labels, f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue().replace("(2, 1)","(1, 2)"), """\
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
...
nonbonded simple: (1, 5)
  SI2
  O4
  distance_model: 3.74919
  vdw_distance: 1
...
  distance_model: 4.08434
  vdw_distance: 1
""",
    selections=[range(6), range(30,36), range(408,413), range(-2,0)])
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
...
nonbonded simple: (1, 5)
  vdw_distance: 1
...
nonbonded asu: (7, 4) x+1/2,-y+1/2,z
  vdw_distance: 2
""",
    selections=[range(3), range(15,18), range(144,146), range(-2,0)])
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
...
nonbonded simple: (1, 5)
  SI2
  O4
  vdw_distance: 1
...
nonbonded asu: (7, 4)
  O6
  O3 x+1/2,-y+1/2,z
  vdw_distance: 2
""",
    selections=[range(5), range(25,30), range(240,244), range(-4,0)])
  #
  sites_cart = drls.start_structure.sites_cart()
  pair_proxies = drls.geometry_restraints_manager.pair_proxies()
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted_by_residual(
    sites_cart=sites_cart,
    labels=site_labels,
    f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert len(out.getvalue().splitlines()) == 50
  assert out.getvalue().splitlines()[-1].find("remaining") < 0
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted_by_residual(
    sites_cart=sites_cart,
    labels=site_labels,
    f=out,
    prefix="0^",
    max_lines=28)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue().replace("e-00", "e-0"), """\
0^Bond restraints sorted by residual:
0^  i - j   ideal  model  delta   weight residual sym.op. j
0^O3  - O4  2.629  2.120  0.509 4.10e-01 1.06e-01
...
0^SI1 - SI1 3.071  3.216 -0.145 2.31e-01 4.83e-03 -x+1/2,-y+1/2,-z+1
0^... (remaining 20 not shown)
""",
    selections=[range(3), range(-2,0)])
  site_labels_long = ["abc"+label+"def" for label in site_labels]
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted_by_residual(
    sites_cart=sites_cart,
    labels=site_labels_long,
    f=out,
    prefix="^0",
    max_lines=28)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue().replace("e-00", "e-0"), """\
^0Bond restraints sorted by residual:
^0   atom i - atom j    ideal  model  delta   weight residual sym.op. j
^0abcO3def  - abcO4def  2.629  2.120  0.509 4.10e-01 1.06e-01
...
^0abcSI1def - abcSI1def 3.071  3.216 -0.145 2.31e-01 4.83e-03 -x+1/2,-y+1/2,-z+1
^0... (remaining 20 not shown)
""",
    selections=[range(3), range(-2,0)])
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted_by_residual(
    sites_cart=sites_cart,
    f=out,
    prefix=".=",
    max_lines=28)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue().replace("e-00", "e-0"), """\
.=Bond restraints sorted by residual:
.=ideal  model  delta   weight residual sym.op. j
.=2.629  2.120  0.509 4.10e-01 1.06e-01
...
.=3.071  3.216 -0.145 2.31e-01 4.83e-03 -x+1/2,-y+1/2,-z+1
.=... (remaining 20 not shown)
""",
    selections=[range(3), range(-2,0)])
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted_by_residual(
    sites_cart=sites_cart,
    f=out,
    prefix="-+",
    max_lines=1)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue().replace("e-00", "e-0"), """\
-+Bond restraints sorted by residual:
-+ideal  model  delta   weight residual
-+2.629  2.120  0.509 4.10e-01 1.06e-01
-+... (remaining 47 not shown)
""")
  out = StringIO()
  pair_proxies.bond_proxies.show_sorted_by_residual(
    sites_cart=sites_cart,
    f=out,
    prefix="=+",
    max_lines=0)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
=+... (remaining 48 not shown)
""")
  #
  out = StringIO()
  pair_proxies.nonbonded_proxies.show_sorted_by_model_distance(
    sites_cart=sites_cart,
    labels=site_labels,
    f=out,
    prefix="d%")
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
d%Nonbonded interactions sorted by model distance:
d%  i - j    model   vdw sym.op. j
d%O3  - O3   3.067 2.000 -x+1/2,-y+1/2,-z
...
d%SI2 - O4   3.386 1.000
""",
    selections=[range(3), range(8,9)])
  out = StringIO()
  pair_proxies.nonbonded_proxies.show_sorted_by_model_distance(
    sites_cart=sites_cart,
    labels=site_labels_long,
    f=out,
    prefix="&u",
    max_lines=7)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
&uNonbonded interactions sorted by model distance:
&u   atom i - atom j     model   vdw sym.op. j
&uabcO3def  - abcO3def   3.067 2.000 -x+1/2,-y+1/2,-z
...
&uabcSI2def - abcO4def   3.386 1.000
&u... (remaining 45 not shown)
""",
    selections=[range(3), range(-2,0)])
  out = StringIO()
  pair_proxies.nonbonded_proxies.show_sorted_by_model_distance(
    sites_cart=sites_cart,
    f=out,
    prefix="*j",
    max_lines=7)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
*jNonbonded interactions sorted by model distance:
*j model   vdw sym.op. j
*j 3.067 2.000 -x+1/2,-y+1/2,-z
...
*j 3.386 1.000
*j... (remaining 45 not shown)
""",
    selections=[range(3), range(-2,0)])
  out = StringIO()
  pair_proxies.nonbonded_proxies.show_sorted_by_model_distance(
    sites_cart=sites_cart,
    f=out,
    prefix="@r",
    max_lines=0)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
@r... (remaining 52 not shown)
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
  try: from mmtbx.monomer_library import pdb_interpretation
  except ImportError:
    print "Skipping exercise_with_pdb():", \
      "mmtbx.monomer_library.pdb_interpretation not available"
    return
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
  assert not show_diff(out.getvalue(), """\
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
  vdw_distance: 2.9
nonbonded simple: (0, 23)
  distance_model: 4.89852
  vdw_distance: 3.75
""",
    selections=[range(4), range(184,188), range(404,409), range(459,464),
      range(488,493), range(-6,0)])
  #
  out = StringIO()
  geo.show_interactions(site_labels=site_labels, f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
bond simple: (0, 1)
  " CA  TYR A   1 "
  " CB  TYR A   1 "
  distance_model: 1.53051
  distance_ideal: 1.53
  weight: 2500
bond simple: (0, 8)
...
angle: (1, 0, 8)
  " CB  TYR A   1 "
  " CA  TYR A   1 "
  " C   TYR A   1 "
  angle_model: 110.543
  angle_ideal: 110.1
  weight: 0.277008
...
dihedral: (2, 0, 1, 11)
  " CG  TYR A   1 "
  " CA  TYR A   1 "
  " CB  TYR A   1 "
  " N   TYR A   1 "
  angle_model: 173.684
  angle_ideal: -180
  weight: 0.00444444
  periodicity: 3
...
chirality: (0, 1, 8, 11)
  " CA  TYR A   1 "
  " CB  TYR A   1 "
  " C   TYR A   1 "
  " N   TYR A   1 "
  volume_model: -2.50215
  volume_ideal: -2.50318
  both_signs: 0
  weight: 25
...
planarity: (0, 8, 10, 15)
  " CA  TYR A   1 " delta: -0.00012, weight: 2500
  " C   TYR A   1 " delta:  0.00042, weight: 2500
  " O   TYR A   1 " delta: -0.00016, weight: 2500
  " N   GLY A   2 " delta: -0.00014, weight: 2500
...
nonbonded asu: (36, 46)
  " C   LEU A   5 "
  " H1  HOH B   3 " -x+1,y+1/2,-z+1/2
  distance_model: 4.89425
  vdw_distance: 2.9
nonbonded simple: (0, 23)
  " CA  TYR A   1 "
  " CD1 PHE A   4 "
  distance_model: 4.89852
  vdw_distance: 3.75
""",
    selections=[range(7), range(276,283), range(661,670), range(760,769),
      range(805,810), range(-10,0)])
  #
  assert approx_equal(flex.min(geo.nonbonded_model_distances()), 0.4777342)
  #
  geo._sites_cart_used_for_pair_proxies = None
  out = StringIO()
  geo.show_interactions(f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
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
  weight: 2500
  weight: 2500
  weight: 2500
...
nonbonded simple: (5, 39)
  vdw_distance: 3.22
...
nonbonded asu: (7, 29) x+1,y,z
  vdw_distance: 3.55
""",
    selections=[range(3), range(138,141), range(303,307), range(347,351),
      range(372,377), range(400,402), range(-2,0)])
  #
  out = StringIO()
  geo.show_interactions(site_labels=site_labels, f=out)
  if (verbose):
    sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), """\
bond simple: (0, 1)
  " CA  TYR A   1 "
  " CB  TYR A   1 "
  distance_ideal: 1.53
  weight: 2500
...
angle: (1, 0, 8)
  " CB  TYR A   1 "
  " CA  TYR A   1 "
  " C   TYR A   1 "
  angle_ideal: 110.1
  weight: 0.277008
...
dihedral: (2, 0, 1, 11)
  " CG  TYR A   1 "
  " CA  TYR A   1 "
  " CB  TYR A   1 "
  " N   TYR A   1 "
  angle_ideal: -180
  weight: 0.00444444
  periodicity: 3
...
chirality: (0, 1, 8, 11)
  " CA  TYR A   1 "
  " CB  TYR A   1 "
  " C   TYR A   1 "
  " N   TYR A   1 "
  volume_ideal: -2.50318
  both_signs: 0
  weight: 25
...
planarity: (0, 8, 10, 15)
  " CA  TYR A   1 " weight: 2500
  " C   TYR A   1 " weight: 2500
  " O   TYR A   1 " weight: 2500
  " N   GLY A   2 " weight: 2500
...
nonbonded simple: (5, 39)
  " CZ  TYR A   1 "
  " O   HOH B   1 "
  vdw_distance: 3.22
...
nonbonded asu: (7, 29)
  " CD1 TYR A   1 "
  " N   PHE A   4 " x+1,y,z
  vdw_distance: 3.55
""",
    selections=[range(5), range(230,236), range(560,568), range(648,656),
      range(689,694), range(717,721), range(-4,0)])

def exercise_all(args):
  verbose = "--verbose" in args
  exercise_with_zeolite(verbose=verbose)
  exercise_with_pdb(verbose=verbose)
  print libtbx.utils.format_cpu_times()

if (__name__ == "__main__"):
  exercise_all(sys.argv[1:])
