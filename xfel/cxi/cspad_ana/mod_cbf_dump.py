# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id: mod_cbf_dump.py 18310 2013-09-28 01:05:21Z hattne $

"""Output image to the file system.
"""


from __future__ import division

__version__ = "$Revision: 18310 $"

from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cftbx.detector.cspad_cbf_tbx import write_cspad_cbf
from parse_calib import calib2sections
import os
import libtbx.load_env

class mod_cbf_dump(common_mode.common_mode_correction):
  """Class for outputting cbf images to the file system within the pyana
  analysis framework.
  """

  def __init__(self,
               address,
               out_dirname,
               out_basename,
               calib_dir=None,
               metrology=None,
               **kwds):
    """The mod_cbf_dump class constructor stores the parameters passed from
    the pyana configuration file in instance variables.

    @param address      Full data source address of the DAQ device
    @param out_dirname  Directory portion of output image pathname
    @param out_basename Filename prefix of output image pathname
    @param calib_dir    Directory with calibration information that common_mode
                        used to do dark subtraction and read the image.  If None,
                        use the Run4 metrology in the xfel source tree.
    @param metrology    Directory with calibration information or cbf file
                        (header only or full file) from which to apply metrology
                        information
    """

    self._calib_dir = cspad_tbx.getOptString(calib_dir)
    if self._calib_dir is None:
      self._calib_dir = libtbx.env.find_in_repositories(
              "xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0")

    super(mod_cbf_dump, self).__init__(address=address, calib_dir=self._calib_dir, **kwds)

    self._basename = cspad_tbx.getOptString(out_basename)
    self._dirname = cspad_tbx.getOptString(out_dirname)
    self._metrology = cspad_tbx.getOptString(metrology)


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.  It
    outputs the detector image associated with the event @p evt to the
    file system.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_cbf_dump, self).event(evt, env)
    if (evt.get('skip_event')):
      return

    # Where the sample-detector distance is not available, set it to
    # zero.
    distance = cspad_tbx.env_distance(self.address, env, self._detz_offset)
    if distance is None:
      distance = 0

    # See r17537 of mod_average.py.
    device = cspad_tbx.address_split(self.address)[2]
    if device == 'Cspad':
      pixel_size = cspad_tbx.pixel_size
      saturated_value = cspad_tbx.dynamic_range
    elif device == 'marccd':
      pixel_size = 0.079346
      saturated_value = 2**16 - 1

    beam_center_x=pixel_size * self.beam_center[0]
    beam_center_y=pixel_size * self.beam_center[1]

    tiles = {}
    data = self.cspad_img

    sections = calib2sections(cspad_tbx.getOptString(self._calib_dir))
    for p in xrange(len(sections)):
      for s in xrange(len(sections[p])):

        # Pull the sensor block from the image, and rotate it back to
        # the "lying down" convention.
        c = sections[p][s].corners_asic()
        k = (int(round(-sections[p][s].angle / 90.0)) + 1) % 4
        for a in xrange(2):
          asic = data.matrix_copy_block(
            i_row=c[a][0],
            i_column=c[a][1],
            n_rows=c[a][2] - c[a][0],
            n_columns=c[a][3] - c[a][1])
          tiles[(0, p, s, a)] = asic.matrix_rot90(k)

    if os.path.isdir(self._metrology):
      metrostyle = 'calibdir'
      from xfel.cftbx.detector.metrology2phil import metrology2phil
      from iotbx import phil
      metro = metrology2phil(self._metrology,False)

      args = [
        "beam_center=(%f,%f)"%(beam_center_x, beam_center_y),
        "timestamp=%s"%        self.timestamp,
        ]

      for arg in args:
        metro = metro.fetch(sources=[phil.parse(arg)])

    else:
      assert os.path.isfile(self._metrology)
      metrostyle = 'cbf'

      from xfel.cftbx.detector.cspad_cbf_tbx import cbf_file_to_basis_dict
      metro = cbf_file_to_basis_dict(self._metrology)

    t = self.timestamp
    s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]

    write_cspad_cbf(tiles, metro, metrostyle, self.timestamp,
      os.path.join(self._dirname, self._basename + s + ".cbf"), self.wavelength, distance)
