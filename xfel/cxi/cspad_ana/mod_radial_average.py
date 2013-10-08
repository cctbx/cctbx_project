# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id: mod_radial_average.py

"""Calculate the radial average of the images in a stream, for each image.
"""
from __future__ import division

__version__ = "$Revision$"

from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx

class mod_radial_average(common_mode.common_mode_correction):
  """Class for calulating radial average of images
  """

  def __init__(self,
               address,
               xtal_target            = None,
               **kwds):
    """The mod_radial_average class constructor stores the parameters passed
    from the pyana configuration file in instance variables.  All
    parameters, except @p address are optional, and hence need not be
    defined in pyana.cfg.

    @param address      Address string XXX Que?!
    @param xtal_target  Phil file with target paramters, including metrology corrections
    """

    super(mod_radial_average, self).__init__(address=address, **kwds)
    self.m_xtal_target          = cspad_tbx.getOptString(xtal_target)


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_radial_average, self).beginjob(evt, env)

  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.
    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_radial_average, self).event(evt, env)
    if (evt.get("skip_event")):
      return

    # This module only applies to detectors for which a distance is
    # available.
    distance = cspad_tbx.env_distance(self.address, env, self._detz_offset)
    if distance is None:
      self.nfail += 1
      self.logger.warning("event(): no distance, shot skipped")
      evt.put(True, "skip_event")
      return

    # See r17537 of mod_average.py.
    device = cspad_tbx.address_split(self.address)[2]
    if device == 'Cspad':
      pixel_size = cspad_tbx.pixel_size
      saturated_value = cspad_tbx.dynamic_range
    elif device == 'marccd':
      pixel_size = 0.079346
      saturated_value = 2**16 - 1

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
      wavelength=self.wavelength,
      xtal_target=self.m_xtal_target)

    from xfel.command_line.radial_average import run
    args = [
      "file_path=XTC stream",
      "xfel_target=%s"%self.m_xtal_target,
      "verbose=False"
    ]

    run(args, d)

  def endjob(self, env):
    """The endjob() function logs the number of processed shots.

    @param env Environment object
    """

    super(mod_radial_average, self).endjob(env)
    if (env.subprocess() >= 0):
      self.logger.info("Subprocess %02d: processed %d shots" %
                       (env.subprocess(), self.nshots))
    else:
      self.logger.info("Processed %d shots" % self.nshots)
