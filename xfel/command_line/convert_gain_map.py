from __future__ import absolute_import, division, print_function
from six.moves import range
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
from iotbx.detectors.cspad_detector_formats import address_and_timestamp_from_detector_format_version
from xfel.cxi.cspad_ana.cspad_tbx import evt_timestamp

# Fake objects to emulate the minimal functionality so that we can reuse the
# functions CsPad2x2Image and CsPadDetector in cspad_tbx.py
def fake_get_config(address, env):
  return fake_config()

cspad_tbx.getConfig = fake_get_config

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

  def quadMask(self):
    return 15

  def roiMask(self, i):
    return 255

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
                  .option("-f", "--detector_format_version",
                          action="store",
                          type="string",
                          help="Detector format version to use for generating active areas and laying out tiles",
                          default=None)
                  .option("-m", "--optical_metrology_path",
                          action="store",
                          type="string",
                          help="Path to slac optical metrology file. If not set, use Run 4 metrology",
                          default=None)
                  .option("-d", "--distance",
                          action="store",
                          type="int",
                          help="Detector distance put into the gain pickle file. Not needed for processing.",
                          default="0")
                  .option("-w", "--wavelength",
                          action="store",
                          type="float",
                          help="Incident beam wavelength put into the gain pickle file. Not needed for processing.",
                          default="0")
                     ).process(args=args)
  output_filename = command_line.options.output_filename
  detector_format_version = command_line.options.detector_format_version
  if detector_format_version is None or 'XPP' not in detector_format_version:
    beam_center_x = None
    beam_center_y = None
  else:
    beam_center_x = 1765 // 2 * 0.11
    beam_center_y = 1765 // 2 * 0.11
  address, timestamp = address_and_timestamp_from_detector_format_version(detector_format_version)

  # if no detector format version is provided, make sure to write no address to the image pickle
  # but CsPadDetector (called later), needs an address, so give it a fake one
  save_address = address is not None
  if not save_address:
    address = "CxiDs1-0|Cspad-0" # time stamp will still be None
  timestamp = evt_timestamp((timestamp,0))
  args = command_line.args
  assert len(args) == 1
  if args[0].endswith('.npy'):
    data = numpy.load(args[0])
    det, active_areas = convert_2x2(data, detector_format_version, address)
  elif args[0].endswith('.txt') or args[0].endswith('.gain'):
    raw_data = numpy.loadtxt(args[0])
    assert raw_data.shape in [(5920, 388), (11840, 194)]
    det, active_areas = convert_detector(raw_data, detector_format_version, address, command_line.options.optical_metrology_path)
  img_diff = det
  img_sel = (img_diff > 0).as_1d()
  gain_map = flex.double(img_diff.accessor(), 0)
  gain_map.as_1d().set_selected(img_sel.iselection(), 1/img_diff.as_1d().select(img_sel))
  gain_map /= flex.mean(gain_map.as_1d().select(img_sel))

  if not save_address:
    address = None
  d = cspad_tbx.dpack(data=gain_map, address=address, active_areas=active_areas, timestamp=timestamp,
    distance=command_line.options.distance,wavelength=command_line.options.wavelength,
    beam_center_x = beam_center_x, beam_center_y = beam_center_y)
  easy_pickle.dump(output_filename, d)


def convert_detector(raw_data, detector_format_version, address, optical_metrology_path=None):
  # https://confluence.slac.stanford.edu/display/PCDS/CSPad+metrology+and+calibration+files%2C+links
  data3d = []
  if raw_data.shape == (5920,388):
    asic_start = 0
    if optical_metrology_path is None:
      calib_dir = libtbx.env.find_in_repositories("xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0")
      sections = parse_calib.calib2sections(calib_dir)
    else:
      sections = parse_calib.calib2sections(optical_metrology_path)
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
    return flex.double(cspad_tbx.CsPadDetector(address, evt, env, sections).astype(numpy.float64)), None
  else:
    asic_start = 0
    if detector_format_version is not None and 'XPP' in detector_format_version:
      from xfel.cxi.cspad_ana.cspad_tbx import xpp_active_areas
      rotations = xpp_active_areas[detector_format_version]['rotations']
      active_areas = xpp_active_areas[detector_format_version]['active_areas']
      det = flex.double([0]*(1765*1765))
      det.reshape(flex.grid((1765,1765)))
      for i in range(64):
        row = active_areas[i*4]
        col = active_areas[i*4 + 1]
        block = flex.double(raw_data[i * 185:(i+1)*185, :])
        det.matrix_paste_block_in_place(block.matrix_rot90(rotations[i]), row, col)
      return det, active_areas

    else:
      if optical_metrology_path is None:
        calib_dir = libtbx.env.find_in_repositories("xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0")
        sections = parse_calib.calib2sections(calib_dir)
      else:
        sections = parse_calib.calib2sections(optical_metrology_path)
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
        data3d.append(fake_cspad_ElementV2(quad_data, i_quad))

      env = fake_env(fake_config())
      evt = fake_evt(data3d)
      beam_center, active_areas = cspad_tbx.cbcaa(fake_config(),sections)
      return flex.double(cspad_tbx.CsPadDetector(address, evt, env, sections).astype(numpy.float64)), active_areas


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
