# LIBTBX_SET_DISPATCHER_NAME cxi.gain_map
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import sys
import numpy

from libtbx import easy_pickle
from libtbx.option_parser import option_parser
from scitbx.array_family import flex
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import parse_calib


# Fake objects to emulate the minimal functionality so that we can reuse the
# functions CsPad2x2Image and CsPadDetector in cspad_tbx.py
class fake_cspad_ElementV2(object):
  def __init__(self, data, quad):
    self._data = data
    self._quad = quad

  def data(self):
    return self._data

  def quad(self):
    return self._quad

class fake_config(object):

  def sections(self, i):
    return list(range(8))



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
  if args[0].endswith('.npy'):
    data = numpy.load(args[0])
    det = convert_2x2(data)
  elif args[0].endswith('.txt'):
    raw_data = numpy.loadtxt(args[0])
    assert raw_data.shape == (5920, 388)
    det = convert_detector(raw_data)
  det = flex.double(det.astype(numpy.float64))
  img_diff = det
  img_sel = (img_diff > 0).as_1d()
  gain_map = flex.double(img_diff.accessor(), 0)
  gain_map.as_1d().set_selected(img_sel.iselection(), 1/img_diff.as_1d().select(img_sel))
  gain_map /= flex.mean(gain_map.as_1d().select(img_sel))
  d = cspad_tbx.dpack(data=gain_map)
  easy_pickle.dump(output_filename, d)


def convert_detector(raw_data):
  # https://confluence.slac.stanford.edu/display/PCDS/CSPad+metrology+and+calibration+files%2C+links
  config = fake_config()
  calib_dir = "/reg/d/ana11/cxi/data/CSPAD-metrology/run4/CxiDs1.0:Cspad.0"
  sections = parse_calib.calib2sections(calib_dir)
  data3d = []
  asic_start = 0
  for i_quad in range(4):
    asic_size = 185 * 388
    section_size = asic_size * 2
    quad_start = i_quad * section_size * 4
    quad_asics = []
    for i_2x2 in range(4):
      for i_asic in range(2):
        asic_end = asic_start + 185
        quad_asics.append(raw_data[asic_start:asic_end, :])
        asic_start = asic_end
    quad_data = numpy.dstack(quad_asics)
    quad_data = numpy.rollaxis(quad_data, 2,0)
    data3d.append(fake_cspad_ElementV2(quad_data, i_quad))
  return cspad_tbx.CsPadDetector(data3d, config, sections)


def convert_2x2(data):
  config = fake_config()
  sections = [[parse_calib.Section(90, (185 / 2 + 0,   (2 * 194 + 3) / 2)),
               parse_calib.Section(90, (185 / 2 + 185, (2 * 194 + 3) / 2))]]
  if data.shape == (2, 185, 388):
    data = numpy.rollaxis(data, 0, 3)
  elif data.shape == (370, 388):
    data = numpy.dstack((data[:185, :], data[185:, :]))
  else:
    raise RuntimeError("Shape of numpy array is %s: I don't know what to do with an array this shape!")
  return cspad_tbx.CsPad2x2Image(data, config, sections)


if __name__ == '__main__':
  run(sys.argv[1:])
