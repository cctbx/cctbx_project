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

def run(args):
  assert len(args) == 0
  exercise_image_simple()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
