# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id: mod_radial_average.py

"""Calculate the radial average of the images in a stream, for each image.
"""
from __future__ import absolute_import, division, print_function

__version__ = "$Revision$"

from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import skip_event_flag
from xfel.cxi.cspad_ana import rayonix_tbx
import os

class mod_radial_average(common_mode.common_mode_correction):
  """Class for calulating radial average of images
  """

  def __init__(self,
               address,
               out_dirname  = None,
               out_basename = None,
               xtal_target  = None,
               two_theta_low = None,
               two_theta_high = None,
               **kwds):
    """The mod_radial_average class constructor stores the parameters passed
    from the pyana configuration file in instance variables.  All
    parameters, except @p address are optional, and hence need not be
    defined in pyana.cfg.

    @param address      Address string XXX Que?!
    @param out_dirname  Optional directory portion of output average pathname
    @param out_basename Optional filename prefix of output average pathname
    @param xtal_target  Phil file with target paramters, including metrology corrections
    @param two_theta_low Optional two theta value of interest
    @param two_theta_high Optional two theta value of interest
    """

    super(mod_radial_average, self).__init__(address=address, **kwds)
    self.m_xtal_target          = cspad_tbx.getOptString(xtal_target)

    self._basename = cspad_tbx.getOptString(out_basename)
    self._dirname = cspad_tbx.getOptString(out_dirname)
    self._two_theta_low = cspad_tbx.getOptFloat(two_theta_low)
    self._two_theta_high = cspad_tbx.getOptFloat(two_theta_high)

    if self._dirname is not None or self._basename is not None:
      assert self._dirname is not None and self._basename is not None

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
      evt.put(skip_event_flag(), "skip_event")
      return

    # See r17537 of mod_average.py.
    device = cspad_tbx.address_split(self.address)[2]
    if device == 'Cspad':
      pixel_size = cspad_tbx.pixel_size
      saturated_value = cspad_tbx.cspad_saturated_value
    elif device == 'marccd':
      pixel_size = 0.079346
      saturated_value = 2**16 - 1
    elif device == 'Rayonix':
      pixel_size = rayonix_tbx.get_rayonix_pixel_size(self.bin_size)
      saturated_value = rayonix_tbx.rayonix_saturated_value

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

    t = self.timestamp
    s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]

    if self._dirname is not None:
      dest_path = os.path.join(self._dirname, self._basename + s + ".txt")
      args.append("output_file=%s"%dest_path)

    self.logger.info("Calculating radial average for image %s"%s)
    xvals, results = run(args, d)

    evt.put(xvals, "cctbx.xfel.radial_average.xvals")
    evt.put(results, "cctbx.xfel.radial_average.results")

    def get_closest_idx(data, val):
      from scitbx.array_family import flex
      deltas = flex.abs(data - val)
      return flex.first_index(deltas, flex.min(deltas))

    if self._two_theta_low is not None:
      i_low = results[get_closest_idx(xvals, self._two_theta_low)]
      evt.put(i_low, "cctbx.xfel.radial_average.two_theta_low")

    if self._two_theta_high is not None:
      i_high = results[get_closest_idx(xvals, self._two_theta_high)]
      evt.put(i_high, "cctbx.xfel.radial_average.two_theta_high")

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

    super(mod_radial_average, self).endjob(env)
    if (env.subprocess() >= 0):
      self.logger.info("Subprocess %02d: processed %d shots" %
                       (env.subprocess(), self.nshots))
    else:
      self.logger.info("Processed %d shots" % self.nshots)
