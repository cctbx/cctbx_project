# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#

"""Create the CSPAD image dict and put it in the event
"""
from __future__ import absolute_import, division, print_function

from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import rayonix_tbx
from xfel.cxi.cspad_ana import skip_event_flag

class mod_image_dict(common_mode.common_mode_correction):
  """Class for saving the CSPAD image dict in the psana event
  """

  def __init__(self,
               address,
               out_key  = "cctbx.xfel.image_dict",
               **kwds):
    """The mod_image_dict class constructor stores the parameters passed
    from the psana configuration file in instance variables.  All
    parameters, except @p address are optional, and hence need not be
    defined in psana.cfg.

    @param address      Full data source address of the DAQ device
    @param out_key      The image dict will be saved in the event with this
                        name
    """

    super(mod_image_dict, self).__init__(address=address, **kwds)

    self.m_out_key = cspad_tbx.getOptString(out_key)


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_image_dict, self).beginjob(evt, env)

  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_image_dict, self).event(evt, env)
    if (evt.get("skip_event")):
      return

    if self.cspad_img is None:
      return

    # This module only applies to detectors for which a distance is
    # available.
    distance = cspad_tbx.env_distance(self.address, env, self._detz_offset)
    if distance is None:
      self.nfail += 1
      self.logger.warning("event(): no distance, shot skipped")
      evt.put(skip_event_flag(), "skip_event")
      return

    device = cspad_tbx.address_split(self.address)[2]

    self.logger.info("Subprocess %02d: process image #%05d @ %s" %
                     (env.subprocess(), self.nshots, self.timestamp))

    # See r17537 of mod_average.py.
    if device == 'Cspad':
      pixel_size = cspad_tbx.pixel_size
      saturated_value = cspad_tbx.cspad_saturated_value
    elif device == 'Rayonix':
      pixel_size = rayonix_tbx.get_rayonix_pixel_size(self.bin_size)
      saturated_value = rayonix_tbx.rayonix_saturated_value
    elif device == 'marccd':
      pixel_size = evt.get("marccd_pixel_size")
      saturated_value = evt.get("marccd_saturated_value")
      if distance == 0:
        distance = evt.get("marccd_distance")

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

    evt.put(d, self.m_out_key)

    # Diagnostic message emitted only when all the processing is done.
    if (env.subprocess() >= 0):
      self.logger.info("Subprocess %02d: accepted #%05d @ %s" %
                       (env.subprocess(), self.nshots, self.timestamp))
    else:
      self.logger.info("Accepted #%05d @ %s" %
                       (self.nshots, self.timestamp))

  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

  def endjob(self, obj1, obj2=None):
    """The endjob() function logs the number of processed shots.

    @param evt Event object (psana only)
    @param env Environment object
    """

    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2

    super(mod_image_dict, self).endjob(env)
    if (env.subprocess() >= 0):
      self.logger.info("Subprocess %02d: processed %d shots" %
                       (env.subprocess(), self.nshots))
    else:
      self.logger.info("Processed %d shots" % self.nshots)
