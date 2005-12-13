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
    relative_path="regression/misc/strudat_zeolite_atlas",
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

def exercise_with_pdb(verbose):
  try: from mmtbx.monomer_library import pdb_interpretation
  except ImportError:
    print "Skipping exercise_with_pdb():", \
      "mmtbx.monomer_library.pdb_interpretation not available"
    return
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="regression/pdb/enk.pdb",
    test=os.path.isfile)
  if (pdb_file is None):
    print "Skipping exercise_with_pdb(): input file not available"
    return
  if (verbose):
    out = sys.stdout
  else:
    out = StringIO()
  processed_pdb_file = pdb_interpretation.run(
    args=[pdb_file],
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
nonbonded asu: (36, 47) -x+1,y+1/2,-z+1/2
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
nonbonded asu: (36, 47)
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
nonbonded simple: (5, 40)
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
nonbonded simple: (5, 40)
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
