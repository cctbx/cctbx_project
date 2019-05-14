# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

"""Output image to the file system.
"""


from __future__ import absolute_import, division, print_function

__version__ = "$Revision$"

from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import rayonix_tbx


class mod_dump(common_mode.common_mode_correction):
  """Class for outputting images to the file system within the pyana
  analysis framework.  XXX This should eventually deprecate the
  'write_dict' dispatch from mod_hitfind.
  """

  def __init__(self, address, out_dirname, out_basename, out_format="pickle", **kwds):
    """The mod_dump class constructor stores the parameters passed from
    the pyana configuration file in instance variables.

    @param address      Full data source address of the DAQ device
    @param out_dirname  Directory portion of output image pathname
    @param out_basename Filename prefix of output image pathname
    @param out_format   Output the data as pickle or TIFF
    """

    super(mod_dump, self).__init__(address=address, **kwds)

    self._basename = cspad_tbx.getOptString(out_basename)
    self._dirname = cspad_tbx.getOptString(out_dirname)
    self._format = cspad_tbx.getOptString(out_format)


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.  It
    outputs the detector image associated with the event @p evt to the
    file system.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_dump, self).event(evt, env)
    if (evt.get('skip_event')):
      return

    if self.cspad_img is None:
      print("No image to save for %s"%self.timestamp)
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
      saturated_value = cspad_tbx.cspad_saturated_value
      output_filename = self._basename
    elif device == 'Rayonix':
      pixel_size = rayonix_tbx.get_rayonix_pixel_size(self.bin_size)
      saturated_value = rayonix_tbx.rayonix_saturated_value
      output_filename = self._basename
    elif device == 'marccd':
      if distance == 0:
        distance = evt.get('marccd_distance')
      pixel_size = 0.079346
      saturated_value = 2**16 - 1
      output_filename = self._basename + evt.get(str, 'mccd_name') + "_"

    d = cspad_tbx.dpack(
      active_areas=self.active_areas,
      address=self.address,
      beam_center_x=pixel_size * self.beam_center[0],
      beam_center_y=pixel_size * self.beam_center[1],
      data=self.cspad_img.iround(), # XXX ouch!
      distance=distance,
      pixel_size=pixel_size,
      saturated_value=saturated_value,
      timestamp=self.timestamp,
      wavelength=self.wavelength)
    if self._format == "pickle":
      cspad_tbx.dwritef(d, self._dirname, output_filename)
    elif self._format == "tiff":
      cspad_tbx.write_tiff(d, self._dirname, output_filename)
    output_filename = None
