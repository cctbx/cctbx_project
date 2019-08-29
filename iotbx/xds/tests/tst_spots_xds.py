from __future__ import absolute_import, division, print_function

from libtbx.test_utils import open_tmp_file
from libtbx.test_utils import approx_equal
from iotbx.xds import spot_xds


def exercise_spots_xds():
  txt = """\
 1104.20 1290.27 2.20 632. -4 -3 -1
 912.22 1303.37 3.49 346. 0 0 0
 1427.55 1339.93 3.34 259. 7 3 -1
 1380.54 1187.54 3.58 243. 4 3 4
 1222.32 1220.09 3.95 241. -1 0 2
 1491.71 1322.33 2.72 237. 9 4 0
 1053.50 1227.71 2.87 221. -6 -4 1
"""

  f = open_tmp_file(suffix="SPOTS.XDS", mode="w")
  f.write(txt)
  f.close()

  spots_in = spot_xds.reader()
  spots_in.read_file(f.name)
  assert approx_equal(
    spots_in.centroid,
    [[1104.2, 1290.27, 2.2], [912.22, 1303.37, 3.49],
     [1427.55, 1339.93, 3.34], [1380.54, 1187.54, 3.58],
     [1222.32, 1220.09, 3.95], [1491.71, 1322.33, 2.72],
     [1053.5, 1227.71, 2.87]])
  assert approx_equal(
    spots_in.intensity,
    [632.0, 346.0, 259.0, 243.0, 241.0, 237.0, 221.0])
  assert approx_equal(
    spots_in.miller_index,
    [[-4, -3, -1], [0, 0, 0], [7, 3, -1], [4, 3, 4],
     [-1, 0, 2], [9, 4, 0], [-6, -4, 1]])

  spots_out = spot_xds.writer(
    spots_in.centroid, spots_in.intensity, spots_in.miller_index)
  f = open_tmp_file(suffix="SPOTS.XDS", mode="wb")
  f.close()
  spots_out.write_file(filename=f.name)
  spots_in = spot_xds.reader()
  spots_in.read_file(f.name)
  assert approx_equal(spots_in.centroid, spots_out.centroids)
  assert approx_equal(spots_in.intensity, spots_out.intensities)
  assert approx_equal(spots_in.miller_index, spots_out.miller_indices)

  # now without miller indices
  spots_out = spot_xds.writer(spots_in.centroid, spots_in.intensity)
  f = open_tmp_file(suffix="SPOTS.XDS", mode="wb")
  f.close()
  spots_out.write_file(filename=f.name)
  spots_in = spot_xds.reader()
  spots_in.read_file(f.name)
  assert approx_equal(spots_in.centroid, spots_out.centroids)
  assert approx_equal(spots_in.intensity, spots_out.intensities)
  assert len(spots_in.miller_index) == 0


def run():
  exercise_spots_xds()

if __name__ == '__main__':
  run()
  print("OK")
