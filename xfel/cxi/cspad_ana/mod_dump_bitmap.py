# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

"""Output image to the file system.
"""


from __future__ import absolute_import, division, print_function

__version__ = "$Revision$"

import logging
import os
from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx


class mod_dump_bitmap(common_mode.common_mode_correction):
  """Class for outputting images to the file system within the pyana
  analysis framework.  XXX This should eventually deprecate the
  'write_dict' dispatch from mod_hitfind.
  """

  def __init__(self, address, out_dirname, out_basename,
               binning=1, brightness=1.0, color_scheme=0,
               format='png', **kwds):
    """The mod_dump_bitmap class constructor stores the parameters passed from
    the pyana configuration file in instance variables.

    @param address      Full data source address of the DAQ device
    @param out_dirname  Directory portion of output image pathname
    @param out_basename Filename prefix of output image pathname

    """

    #define COLOR_GRAY 0
    #define COLOR_RAINBOW 1
    #define COLOR_HEAT 2
    #define COLOR_INVERT 3

    super(mod_dump_bitmap, self).__init__(address=address, **kwds)

    self._basename = cspad_tbx.getOptString(out_basename)
    self._dirname = cspad_tbx.getOptString(out_dirname)
    self._binning = cspad_tbx.getOptInteger(binning)
    self._brightness = cspad_tbx.getOptFloat(brightness)
    self._color_scheme = cspad_tbx.getOptInteger(color_scheme)
    self._format = cspad_tbx.getOptString(format)
    self._ext = self._format.lower()
    self._logger = logging.getLogger(self.__class__.__name__)
    if (not os.path.isdir(self._dirname)):
      os.makedirs(self._dirname)


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.  It
    outputs the detector image associated with the event @p evt to the
    file system.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_dump_bitmap, self).event(evt, env)
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
      saturated_value = cspad_tbx.cspad_saturated_value
    elif device == 'marccd':
      pixel_size = 0.079346
      saturated_value = 2**16 - 1

    from iotbx.detectors import FlexImage_d as FlexImage
    vendortype = device
    saturation = 65535
    flex_img = FlexImage(
      rawdata=self.cspad_img,
      binning=self._binning,
      vendortype=vendortype,
      brightness=self._brightness,
      saturation=saturated_value)

    flex_img.setWindow(0, 0, 1)
    flex_img.adjust(color_scheme=self._color_scheme)
    flex_img.prep_string()
    try:
      import PIL.Image as Image
    except ImportError:
      import Image
    # XXX is size//self._binning safe here?
    try:
      pil_img = Image.fromstring(
        'RGB', (flex_img.size2()//self._binning,
                flex_img.size1()//self._binning),
        flex_img.export_string)
    except NotImplementedError:
      pil_img = Image.frombytes(
        'RGB', (flex_img.size2()//self._binning,
                flex_img.size1()//self._binning),
        flex_img.export_string)

    # The output path should not contain any funny characters which may
    # not work in all environments.  This constructs a sequence number a
    # la evt_seqno() from the dictionary's timestamp.
    t = self.timestamp
    s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]

    path = os.path.join(
      self._dirname, self._basename + s + '.' + self._ext)

    self._logger.info("Exporting %s" %path)
    tmp_stream = open(path, 'wb')
    pil_img.save(tmp_stream, format=self._format)
    tmp_stream.close()
