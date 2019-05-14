# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$

"""The mod_hdf5 module writes all events to a single HDF5
(hierarchical data format, version 5) file.  XXX Not sure how this
module will behave under multiprocessing.
"""
from __future__ import absolute_import, division, print_function

__version__ = "$Revision$"

import h5py

from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx


class mod_hdf5(common_mode.common_mode_correction):
  def __init__(self, address, path, **kwds):
    """The constructor stores the path the output HDF5 file.  All
    intermediate directories must exist.

    @param address Address string XXX Que?!
    @param path    Path to output HDF5 file
    """

    super(mod_hdf5, self).__init__(address=address, **kwds)
    self._file = h5py.File(cspad_tbx.getOptString(path), 'w')


  def __del__(self):
    """The destructor closes the HDF5 file opened by __init__().
    """

    self._file.close()


  def event(self, evt, env):
    """The event() function creates a HDF5 group for the event, unless
    it contains a "skip_event" object with value @c True.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_hdf5, self).event(evt, env)
    if (evt.get('skip_event')):
      return

    # If no detector distance is available set it to NaN, since
    # Python's None is not permitted in HDF5
    distance = cspad_tbx.env_distance(self.address, env, self._detz_offset)
    if distance is None:
      distance = float('nan')

    cspad_tbx.hdf5pack(
      hdf5_file=self._file,
      active_areas=self.active_areas,
      address=self.address,
      attenuation=self.sifoil,
      beam_center_x=cspad_tbx.pixel_size * self.beam_center[0],
      beam_center_y=cspad_tbx.pixel_size * self.beam_center[1],
      ccd_image_saturation=cspad_tbx.cspad_saturated_value,
      data=self.cspad_img,
      distance=distance,
      pixel_size=cspad_tbx.pixel_size,
      pulse_length=self.pulse_length,
      saturated_value=cspad_tbx.cspad_saturated_value,
      timestamp=self.timestamp,
      wavelength=self.wavelength,
      xtal_target=repr(None))
