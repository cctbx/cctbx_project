from cctbx import maptbx
import iotbx.ccp4_map
from libtbx.test_utils import approx_equal
import sys

def run(args):
  for file_name in args:
    print file_name
    m = iotbx.ccp4_map.map_reader(file_name=file_name)
    print "header_min: ", m.header_min
    print "header_max: ", m.header_min
    print "header_mean:", m.header_mean
    print "header_rms: ", m.header_rms
    map_stats = maptbx.statistics(m.data)
    assert approx_equal(map_stats.min(), m.header_min)
    assert approx_equal(map_stats.max(), m.header_max)
    assert approx_equal(map_stats.mean(), m.header_mean)
    if (m.header_rms != 0):
      assert approx_equal(map_stats.sigma(), m.header_rms)
    print m.space_group
    print m.unit_cell
    print m.map_grid
    print m.map_origin
    print m.axes_order
    print m.map_dim

if (__name__ == "__main__"):
  run(sys.argv[1:])
