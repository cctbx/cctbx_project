from cctbx.geometry_restraints.distance_least_squares \
  import distance_and_repulsion_least_squares
from iotbx.kriber import strudat
from cStringIO import StringIO
import libtbx.utils
from libtbx.test_utils import show_diff
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
...
nonbonded simple: (1, 5)
  SI2
  O4
  distance_model: 3.74919
  vdw_distance: 1
...
nonbonded asu: (1, 2)
  SI2
  O1 x,y,z
  distance_model: 4.08434
  vdw_distance: 1
""",
    selections=[range(6), range(30,36), range(408,413), range(-5,0)])
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
  assert not show_diff(out.getvalue(), """\
bond simple: (0, 1)
  " CA  TYR     1 "
  " CB  TYR     1 "
  distance_model: 1.53051
  distance_ideal: 1.53
  weight: 2500
bond simple: (0, 8)
...
angle: (1, 0, 8)
  " CB  TYR     1 "
  " CA  TYR     1 "
  " C   TYR     1 "
  angle_model: 110.543
  angle_ideal: 110.1
  weight: 0.277008
...
dihedral: (2, 0, 1, 11)
  " CG  TYR     1 "
  " CA  TYR     1 "
  " CB  TYR     1 "
  " N   TYR     1 "
  angle_model: 173.684
  angle_ideal: -180
  weight: 0.00444444
  periodicity: 3
...
chirality: (0, 1, 8, 11)
  " CA  TYR     1 "
  " CB  TYR     1 "
  " C   TYR     1 "
  " N   TYR     1 "
  volume_model: -2.50215
  volume_ideal: -2.50318
  both_signs: 0
  weight: 25
...
planarity: (0, 8, 10, 15)
  " CA  TYR     1 " delta: -0.00012, weight: 2500
  " C   TYR     1 " delta:  0.00042, weight: 2500
  " O   TYR     1 " delta: -0.00016, weight: 2500
  " N   GLY     2 " delta: -0.00014, weight: 2500
...
nonbonded asu: (36, 47)
  " C   LEU     5 "
  " H1  HOH     3 " -x+1,y+1/2,-z+1/2
  distance_model: 4.89425
  vdw_distance: 2.9
nonbonded simple: (0, 23)
  " CA  TYR     1 "
  " CD1 PHE     4 "
  distance_model: 4.89852
  vdw_distance: 3.75
""",
    selections=[range(7), range(276,283), range(661,670), range(760,769),
      range(805,810), range(-10,0)])
  #
  geo._sites_cart_used_for_pair_proxies = None
  out = StringIO()
  geo.show_interactions(f=out)
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
  assert not show_diff(out.getvalue(), """\
bond simple: (0, 1)
  " CA  TYR     1 "
  " CB  TYR     1 "
  distance_ideal: 1.53
  weight: 2500
...
angle: (1, 0, 8)
  " CB  TYR     1 "
  " CA  TYR     1 "
  " C   TYR     1 "
  angle_ideal: 110.1
  weight: 0.277008
...
dihedral: (2, 0, 1, 11)
  " CG  TYR     1 "
  " CA  TYR     1 "
  " CB  TYR     1 "
  " N   TYR     1 "
  angle_ideal: -180
  weight: 0.00444444
  periodicity: 3
...
chirality: (0, 1, 8, 11)
  " CA  TYR     1 "
  " CB  TYR     1 "
  " C   TYR     1 "
  " N   TYR     1 "
  volume_ideal: -2.50318
  both_signs: 0
  weight: 25
...
planarity: (0, 8, 10, 15)
  " CA  TYR     1 " weight: 2500
  " C   TYR     1 " weight: 2500
  " O   TYR     1 " weight: 2500
  " N   GLY     2 " weight: 2500
...
nonbonded simple: (5, 40)
  " CZ  TYR     1 "
  " O   HOH     1 "
  vdw_distance: 3.22
...
nonbonded asu: (7, 29)
  " CD1 TYR     1 "
  " N   PHE     4 " x+1,y,z
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
