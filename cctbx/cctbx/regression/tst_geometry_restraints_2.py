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

def exercise_all(args):
  verbose = "--verbose" in args
  exercise_with_zeolite(verbose=verbose)
  print libtbx.utils.format_cpu_times()

if (__name__ == "__main__"):
  exercise_all(sys.argv[1:])
