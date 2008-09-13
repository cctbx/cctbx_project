from cctbx.array_family import flex
import iotbx.ccp4_map
import sys

def run(args):
  for file_name in args:
    print file_name
    iotbx.ccp4_map.test_read(file_name=file_name)
    print
    m = iotbx.ccp4_map.map_reader(file_name=file_name)
    d = m.data
    print d
    print flex.min(d), flex.max(d)
    print m.space_group
    print m.unit_cell
    print m.map_min
    print m.map_max
    print m.map_grid
    print m.map_origin
    print m.axes_order
    print m.map_dim

if (__name__ == "__main__"):
  run(sys.argv[1:])
