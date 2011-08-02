import sys, os
op = os.path

def exercise_image_simple():
  expected_sum_image_pixels_iter = iter(
      (200, 156) + (400, 308)*2
    + (450, 351) + (900, 693)*2
    + ( 91,  69) + (181, 139)*2
    + (160, 124) + (320, 248)*2
    + (246, 188) + (490, 377)*2)
  from cctbx import uctbx
  unit_cell = uctbx.unit_cell((11,12,13,85,95,105))
  from cctbx.array_family import flex
  miller_indices = flex.miller_index([(-1,2,1)])
  from scitbx.math import euler_angles
  crystal_rotation_matrix = euler_angles.xyz_matrix(80,20,30)
  import rstbx.simage
  from libtbx.test_utils import approx_equal, show_diff
  dpx, dpy = 4, 5
  for ewald_proximity,star in [(0.1, " "), (0.5, "*")]:
    image_lines = []
    for point_spread in xrange(1,5+1):
      for spot_intensity_factor in [0.5, 1, None]:
        if (spot_intensity_factor is None):
          spot_intensity_factors = None
        else:
          spot_intensity_factors = flex.double([spot_intensity_factor])
        for apply_proximity_factor in [False, True]:
          if (star == "*"):
            expected_sum_image_pixels = expected_sum_image_pixels_iter.next()
          for code in xrange(16):
            store_miller_index_i_seqs = bool(code & 0x1)
            store_spots = bool(code & 0x2)
            store_signals = bool(code & 0x4)
            set_pixels = bool(code & 0x8)
            image = rstbx.simage.image_simple(
              apply_proximity_factor=apply_proximity_factor,
              store_miller_index_i_seqs=store_miller_index_i_seqs,
              store_spots=store_spots,
              store_signals=store_signals,
              set_pixels=set_pixels).compute(
                unit_cell=unit_cell,
                miller_indices=miller_indices,
                spot_intensity_factors=spot_intensity_factors,
                crystal_rotation_matrix=crystal_rotation_matrix,
                ewald_radius=0.5,
                ewald_proximity=ewald_proximity,
                signal_max=100,
                detector_distance=5,
                detector_size=(10,12),
                detector_pixels=(dpx,dpy),
                point_spread=point_spread,
                gaussian_falloff_scale=4)
            if (store_signals and image.signals.size() == 1):
              partialities = rstbx.simage.image_simple(
                apply_proximity_filter=False,
                apply_proximity_factor=apply_proximity_factor,
                store_signals=True).compute(
                  unit_cell=unit_cell,
                  miller_indices=miller_indices,
                  spot_intensity_factors=None,
                  crystal_rotation_matrix=crystal_rotation_matrix,
                  ewald_radius=0.5,
                  ewald_proximity=ewald_proximity,
                  signal_max=1,
                  detector_distance=5,
                  detector_size=(10,12),
                  detector_pixels=(dpx,dpy),
                  point_spread=point_spread,
                  gaussian_falloff_scale=4).signals
              f = 100
              if (spot_intensity_factor is not None):
                f *= spot_intensity_factor
              assert approx_equal(partialities*f, image.signals)
            if (store_miller_index_i_seqs and star == "*"):
              assert image.miller_index_i_seqs.size() == 1
            else:
              assert image.miller_index_i_seqs.size() == 0
            if (store_spots and star == "*"):
              assert image.spots.size() == 1
            else:
              assert image.spots.size() == 0
            if (store_signals and star == "*"):
              assert image.signals.size() == 1
            else:
              assert image.signals.size() == 0
            if (not set_pixels):
              assert image.pixels.size() == 0
            else:
              assert image.pixels.size() == 20
              sum_image_pixels = flex.sum(image.pixels)
              if (star == "*"):
                assert sum_image_pixels == expected_sum_image_pixels
              else:
                assert sum_image_pixels == 0
      assert image.pixels.all() == (dpx,dpy)
      for i in xrange(dpx):
        line = []
        for j in xrange(dpy):
          if (image.pixels[(i,j)]): c = star
          else: c = " "
          line.append(c)
        image_lines.append("|"+"".join(line)+"|")
      image_lines.append("")
    assert not show_diff("\n".join(image_lines), """\
|     |
|     |
|  ** |
|  ** |

|     |
|  ***|
|  ***|
|  ***|

|     |
|  ** |
| *** |
|  ** |

|  *  |
| ****|
| ****|
| ****|

| *** |
|*****|
|*****|
|*****|
""".replace("*", star))

def exercise_create():
  from rstbx.simage import create
  from libtbx.test_utils import block_show_diff
  from libtbx.str_utils import show_string
  from cStringIO import StringIO
  #
  def check(args, expected_block):
    sio = StringIO()
    work_params = create.process_args(args=args, out=sio)
    assert not block_show_diff(sio.getvalue(), expected_block)
    pixels = create.compute_image(work_params)
    assert pixels.all() == tuple(work_params.detector.pixels)
  #
  check(args=[], expected_block="""\
detector {
  distance = 180
  size = 200 200
  pixels = 1000 1000
}
""")
  #
  relative_path = "phenix_regression/pdb/start.pdb"
  import libtbx.load_env
  full_path = libtbx.env.find_in_repositories(
    relative_path=relative_path, test=op.isfile)
  if (full_path is None):
    print "Skipping some tests due to missing file: %s" % relative_path
  else:
    check(
      args=["pdb_file="+show_string(full_path)],
      expected_block="""\
change_of_basis_op_to_niggli_cell = "a,b,c"
unit_cell = 32.9 32.9 96.1 90 90 120
intensity_symmetry = "P 3 2 1"
lattice_symmetry = "P 6 2 2"
""")

def exercise_explore_completeness():
  import libtbx.load_env
  if (not libtbx.env.has_module("spotfinder")):
    print "Skipping some tests due to missing module: spotfinder"
    return
  from libtbx.test_utils import contains_substring
  from libtbx import easy_run
  def run(args):
    cmd = " ".join(["rstbx.simage.explore_completeness"] + args)
    print cmd
    buf = easy_run.fully_buffered(
      command=cmd, stdout_splitlines=False).raise_if_errors().stdout_buffer
    for key in [
          "Complete with ",
          "Observations per reflection:",
          "  Median: "]:
      assert contains_substring(buf, key)
    return buf
  run(["d_min=10"])
  args = ["d_min=10", "intensity_symmetry=P4", "use_symmetry=True"]
  from libtbx import easy_mp
  if (easy_mp.detect_problem() is None):
    args.append("multiprocessing=True")
  buf = run(args)
  assert contains_substring(buf, 'lattice_symmetry = "P 4 2 2"')

def exercise_solver():
  import libtbx.load_env
  if (not libtbx.env.has_module("spotfinder")):
    print "Skipping some tests due to missing module: spotfinder"
    return
  from libtbx.test_utils import block_show_diff, contains_substring
  from libtbx import easy_run
  def run(args):
    cmd = " ".join(["rstbx.simage.solver"] + args)
    print cmd
    buf = easy_run.fully_buffered(
      command=cmd, stdout_splitlines=False).raise_if_errors().stdout_buffer
    for key in [
          "Final:"]:
      assert contains_substring(buf, key)
    return buf
  buf = run(["d_min=5"])
  assert not block_show_diff(buf, """\
input_im0_i_perm: 0

Correlation of input and estimated I-obs:
  i_perm=0:  1.00000
""")
  buf = run(["d_min=5", "lattice_symmetry=R32:R", "intensity_symmetry=R3:R"])
  assert not block_show_diff(buf, """\
input_im0_i_perm: 1

Correlation of input and estimated I-obs:
  i_perm=0:  0.06799
  i_perm=1:  1.00000
""")
  buf = run(["d_min=5", "lattice_symmetry=R32:R", "intensity_symmetry=P1"])
  assert not block_show_diff(buf, """\
input_im0_i_perm: 5

Correlation of input and estimated I-obs:
  i_perm=0:  0.07524
  i_perm=1: -0.02385
  i_perm=2: -0.04577
  i_perm=3: -0.00099
  i_perm=4:  0.00764
  i_perm=5:  1.00000
""")
  if (not libtbx.env.has_module("labelit")):
    print "Skipping some tests due to missing module: labelit"
  else:
    from libtbx import easy_mp
    mp_problem = easy_mp.detect_problem()
    if (mp_problem is not None):
      print "Skipping some tests:", mp_problem
    else:
      buf = run([
        "d_min=5", "lattice_symmetry=P422", "intensity_symmetry=P4",
        "index_and_integrate=True", "multiprocessing=True"])
      assert contains_substring(buf, "Refined unit cell 9 (")
      assert contains_substring(
        buf, "Correlation of input and estimated I-obs:")
      assert contains_substring(buf, "  Best correlation:  0.999")

def run(args):
  assert len(args) == 0
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  exercise_image_simple()
  exercise_create()
  exercise_explore_completeness()
  exercise_solver()
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
