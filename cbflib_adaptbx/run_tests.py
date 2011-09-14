def run(args):
  assert len(args) == 0
  import libtbx.load_env
  import os
  op = os.path
  cbf = libtbx.env.under_dist(
    module_name="cbflib",
    path="examples/fit2d_data.cbf")
  assert op.isfile(cbf)
  from cbflib_adaptbx.command_line import dump
  from cStringIO import StringIO
  sio = StringIO()
  dump.process(file_name=cbf, out=sio)
  from libtbx.test_utils import show_diff
  assert not show_diff(sio.getvalue(), """\
File name: %s
Number of blocks: 1
  Block name: image_1
  Number of categories: 12
    Category name: diffrn
    Category name: diffrn_source
    Category name: diffrn_radiation
    Category name: diffrn_radiation_wavelength
    Category name: diffrn_measurement
    Category name: diffrn_detector
    Category name: diffrn_detector_element
    Category name: diffrn_data_frame
    Category name: array_structure_list
    Category name: array_element_size
    Category name: array_intensities
    Category name: array_data

""" % cbf)
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
