from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.gain_map
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import sys,time,math
import numpy

from libtbx import easy_pickle
import libtbx.load_env
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

class fake_env(object):
  # XXX Not tested!

  def __init__(self, config):
    self._config = config

  def getConfig(self, Id, address):
    return self._config

class fake_evt(object):
  # XXX Not tested!

  def __init__(self, data3d):
    self._data3d = data3d

  def getCsPadQuads(self, address, env):
    return self._data3d

  def getTime(self):
    class fakeTime(object):
      def __init__(self):
        t = time.time()
        s = int(math.floor(t))
        self.s = s
        self.n = int(round((t - s) * 1000))

      def seconds(self): return self.s
      def nanoseconds(self): return self.n
    return fakeTime()

def run(args):
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
  elif args[0].endswith('.txt') or args[0].endswith('.gain'):
    raw_data = numpy.loadtxt(args[0])
    assert raw_data.shape in [(5920, 388), (11840, 194)]
    det = convert_detector(raw_data)
  img_diff = det
  img_sel = (img_diff > 0).as_1d()
  gain_map = flex.double(img_diff.accessor(), 0)
  gain_map.as_1d().set_selected(img_sel.iselection(), 1/img_diff.as_1d().select(img_sel))
  gain_map /= flex.mean(gain_map.as_1d().select(img_sel))
  d = cspad_tbx.dpack(data=gain_map)
  easy_pickle.dump(output_filename, d)


def convert_detector(raw_data):
  # https://confluence.slac.stanford.edu/display/PCDS/CSPad+metrology+and+calibration+files%2C+links
  data3d = []
  if raw_data.shape == (5920,388):
    asic_start = 0
    calib_dir = libtbx.env.find_in_repositories("xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0")
    sections = parse_calib.calib2sections(calib_dir)
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
    env = fake_env(fake_config())
    evt = fake_evt(data3d)
    return flex.double(cspad_tbx.CsPadDetector('CxiDs1-0|Cspad-0', evt, env, sections).astype(numpy.float64))
  else:
    asic_start = 0
    use_v1_2_metrology = True
    if use_v1_2_metrology:
      rotations = [
                   3,3,3,3,2,2,2,2,1,1,1,1,2,2,2,2,
                   2,2,2,2,1,1,1,1,0,0,0,0,1,1,1,1,
                   1,1,1,1,0,0,0,0,3,3,3,3,0,0,0,0,
                   0,0,0,0,3,3,3,3,2,2,2,2,3,3,3,3
                  ]
      active_areas = flex.int([
           865, 1121, 1059, 1306, 1062, 1121, 1256, 1306,
           864,  909, 1058, 1094, 1061,  909, 1255, 1094,
          1083, 1534, 1268, 1728, 1083, 1337, 1268, 1531,
           871, 1538, 1056, 1732,  871, 1341, 1056, 1535,
          1495, 1326, 1689, 1511, 1298, 1326, 1492, 1511,
          1496, 1539, 1690, 1724, 1299, 1539, 1493, 1724,
          1482, 1105, 1667, 1299, 1482,  908, 1667, 1102,
          1270, 1107, 1455, 1301, 1270,  910, 1455, 1104,
          1123,  706, 1308,  900, 1123,  509, 1308,  703,
           910,  706, 1095,  900,  910,  509, 1095,  703,
          1535,  498, 1729,  683, 1338,  498, 1532,  683,
          1534,  711, 1728,  896, 1337,  711, 1531,  896,
          1324,   77, 1509,  271, 1324,  274, 1509,  468,
          1537,   75, 1722,  269, 1537,  272, 1722,  466,
          1104,   97, 1298,  282,  907,   97, 1101,  282,
          1105,  310, 1299,  495,  908,  310, 1102,  495,
           706,  457,  900,  642,  509,  457,  703,  642,
           705,  669,  899,  854,  508,  669,  702,  854,
           496,   36,  681,  230,  496,  233,  681,  427,
           709,   38,  894,  232,  709,  235,  894,  429,
            77,  256,  271,  441,  274,  256,  468,  441,
            77,   44,  271,  229,  274,   44,  468,  229,
            98,  467,  283,  661,   98,  664,  283,  858,
           311,  467,  496,  661,  311,  664,  496,  858,
           457,  867,  642, 1061,  457, 1064,  642, 1258,
           670,  865,  855, 1059,  670, 1062,  855, 1256,
            37, 1084,  231, 1269,  234, 1084,  428, 1269,
            37,  871,  231, 1056,  234,  871,  428, 1056,
           256, 1495,  441, 1689,  256, 1298,  441, 1492,
            43, 1497,  228, 1691,   43, 1300,  228, 1494,
           469, 1481,  663, 1666,  666, 1481,  860, 1666,
           467, 1269,  661, 1454,  664, 1269,  858, 1454])
      det = flex.double([0]*(1765*1765))
      det.reshape(flex.grid((1765,1765)))
      for i in xrange(64):
        row = active_areas[i*4]
        col = active_areas[i*4 + 1]
        block = flex.double(raw_data[i * 185:(i+1)*185, :])
        det.matrix_paste_block_in_place(block.matrix_rot90(rotations[i]), row, col)
      return det

    else:
      quad_order = [2,3,0,1]
      for i_quad in range(4):
        asic_size = 185 * 194
        section_size = asic_size * 4
        quad_start = i_quad * section_size * 4
        quad_asics = []
        for i_2x2 in range(4):
          for i_asic in range(2):
            asic_end = asic_start + 185
            a = raw_data[asic_start:asic_end, :]
            asic_start = asic_end

            asic_end = asic_start + 185
            b = raw_data[asic_start:asic_end, :]
            asic_start = asic_end

            quad_asics.append(numpy.concatenate((a,b),axis=1))
        quad_data = numpy.dstack(quad_asics)
        quad_data = numpy.rollaxis(quad_data, 2,0)
        data3d.append(fake_cspad_ElementV2(quad_data, quad_order[i_quad]))

      env = fake_env(fake_config())
      evt = fake_evt(data3d)
      return flex.double(cspad_tbx.CsPadDetector('XppGon-0|Cspad-0', evt, env, sections).astype(numpy.float64))


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
