# LIBTBX_SET_DISPATCHER_NAME cxi.gain_map
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import sys
import numpy

from libtbx import easy_pickle
from libtbx.option_parser import option_parser
from scitbx.array_family import flex
from xfel.cxi.cspad_ana import cspad_tbx

def run(args):
  print args
  command_line = (option_parser()
                  .option("-o", "--output_filename",
                          action="store",
                          type="string",
                          help="Filename for the output pickle file",
                          default="gain_map.pickle")
                  ).process(args=args)
  output_filename = command_line.options.output_filename
  args = command_line.args
  assert len(args) == 1

  assert args[0].endswith('.npy')
  data = numpy.load(args[0])
  gain_map = convert(data)
  d = cspad_tbx.dpack(data=gain_map)
  easy_pickle.dump(output_filename, d)


def convert(data):

  angles = (90, 90)
  centers = ((92.5, 195.5), (277.5, 195.5))

  det  = numpy.zeros((2 * 185, 2 * 194 + 3))
  mask = [[1], []]

  # Some code below borrowed from function CsPad2x2Image in xfel/cxi/cspad_tbx.py
  print data.shape
  for s in xrange(2):
    # Philip's data comes in at least two formats depending on when it was created
    if data.shape == (2, 185, 388):
      asics = numpy.vsplit(numpy.rot90(data[s, :, :], -1), 2)
    elif data.shape == (370, 388):
      asics = numpy.vsplit(numpy.rot90(data[185*s:185*(s+1), :], -1), 2)
    else:
      raise RuntimeError("Shape of numpy array is %s: I don't know what to do with an array this shape!")
    gap    = numpy.zeros((3, 185), dtype = data.dtype)
    s_data = numpy.vstack((asics[0], gap, asics[1]))

    angle  = angles[s]
    center = centers[s]
    cspad_tbx.rplace(det, s_data, angle, center)
  det = flex.double(det.astype(numpy.float64))
  img_diff = det
  img_sel = (img_diff > 0).as_1d()
  gain_map = flex.double(img_diff.accessor(), 0)
  gain_map.as_1d().set_selected(img_sel.iselection(), 1/img_diff.as_1d().select(img_sel))
  gain_map /= flex.mean(gain_map.as_1d().select(img_sel))

  return gain_map


if __name__ == '__main__':
  run(sys.argv[1:])
