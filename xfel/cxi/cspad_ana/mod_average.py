# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$
"""First- and second-order statistics for CS-PAD images

The mod_average user analysis module computes the (arithmetic) mean and the
standard deviation from the images in an XTC stream.  On successful
completion, the mean and standard deviation images are written to disk
as pickled dictionaries.
"""
from __future__ import absolute_import, division, print_function

__version__ = "$Revision$"

from xfel.cxi.cspad_ana import average_tbx
from xfel.cxi.cspad_ana import cspad_tbx


class mod_average(average_tbx.average_mixin):
  """Class for generating first- and second-order statistics within
  the pyana framework

  XXX Maybe this module should be renamed to mod_stat12, mod_sstat or
  some such?
  """


  def __init__(self,
               address,
               max_out=None,
               mean_out=None,
               std_out=None,
               **kwds):
    """The mod_average class constructor stores the parameters passed
    from the pyana configuration file in instance variables.  All
    parameters, except @p address are optional, and hence need not be
    defined in pyana.cfg.

    @param address Full data source address of the DAQ device
    """

    super(mod_average, self).__init__(address=address, **kwds)

    self._max_out = cspad_tbx.getOptString(max_out)
    self._mean_out = cspad_tbx.getOptString(mean_out)
    self._std_out = cspad_tbx.getOptString(std_out)

    # XXX Ugly hack here instead of modifying avg_tbx.py
    if max_out is not None:
      self._have_max = True
    if mean_out is not None:
      self._have_mean = True
    if std_out is not None:
      self._have_std = True


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_average, self).beginjob(evt, env)

    # XXX Potential problem if one uses stream substitution, which may
    # change between the time of substitution, and the time the stream
    # is written.  Since files are written in endjob() (which doesn't
    # trigger on an event) substitutions cannot be made there.  Just
    # note the caveats here, perhaps?

#    print "*** P _mean_out : ", self._mean_out
#    print "*** P _max_out  : ", self._max_out
#    print "*** P _std_out  : ", self._std_out

    if self._mean_out is not None:
      self._mean_out = cspad_tbx.pathsubst(self._mean_out, evt, env)

    if self._max_out is not None:
      self._max_out = cspad_tbx.pathsubst(self._max_out, evt, env)

    if self._std_out is not None:
      self._std_out = cspad_tbx.pathsubst(self._std_out, evt, env)

#    print "*** S _mean_out : ", self._mean_out
#    print "*** S _max_out  : ", self._max_out
#    print "*** S _std_out  : ", self._std_out


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_average, self).event(evt, env)
    if evt.get('skip_event'):
      return
    self.logger.info("Shot number %i"  % self._nmemb)

  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

  def endjob(self, obj1, obj2=None):
    """The endjob() function writes the mean and standard deviation images
    to disk.

    @param evt Event object (psana only)
    @param env Environment object
    """
    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2

    stats = super(mod_average, self).endjob(env)
    if stats is None:
      return

    device = cspad_tbx.address_split(self.address)[2]
    if device == 'Andor':
      beam_center = (0, 0) # XXX Fiction!
      pixel_size = 13.5e-3 # XXX Should not be hardcoded here!
      saturated_value = 10000
    elif device == 'Cspad' or device == 'Cspad2x2':
      beam_center = self.beam_center
      pixel_size = cspad_tbx.pixel_size
      saturated_value = cspad_tbx.cspad_saturated_value
    elif device == 'marccd':
      beam_center = tuple(t // 2 for t in d['mean_img'].focus())
      pixel_size = 0.079346
      saturated_value = 2**16 - 1

    if stats['nmemb'] > 0:
      if self.avg_dirname  is not None or \
         self.avg_basename is not None or \
         self._mean_out    is not None:
        d = cspad_tbx.dpack(
          active_areas=self.active_areas,
          address=self.address,
          beam_center_x=pixel_size * beam_center[0],
          beam_center_y=pixel_size * beam_center[1],
          data=stats['mean_img'],
          distance=stats['distance'],
          pixel_size=pixel_size,
          saturated_value=saturated_value,
          timestamp=cspad_tbx.evt_timestamp(stats['time']),
          wavelength=stats['wavelength'])
        if self._mean_out is not None:
          p = cspad_tbx.dwritef2(d, self._mean_out)
        else:
          p = cspad_tbx.dwritef(d, self.avg_dirname, self.avg_basename)
        self.logger.info("Average written to %s" % p)

      if self.stddev_dirname  is not None or \
         self.stddev_basename is not None or \
         self._std_out    is not None:
        d = cspad_tbx.dpack(
          active_areas=self.active_areas,
          address=self.address,
          beam_center_x=pixel_size * beam_center[0],
          beam_center_y=pixel_size * beam_center[1],
          data=stats['std_img'],
          distance=stats['distance'],
          pixel_size=pixel_size,
          saturated_value=saturated_value,
          timestamp=cspad_tbx.evt_timestamp(stats['time']),
          wavelength=stats['wavelength'])
        if self._std_out is not None:
          p = cspad_tbx.dwritef2(d, self._std_out)
        else:
          p = cspad_tbx.dwritef(d, self.stddev_dirname, self.stddev_basename)
        self.logger.info("Standard deviation written to %s" % p)

      if self.max_dirname  is not None or \
         self.max_basename is not None or \
         self._max_out    is not None:
        d = cspad_tbx.dpack(
          active_areas=self.active_areas,
          address=self.address,
          beam_center_x=pixel_size * beam_center[0],
          beam_center_y=pixel_size * beam_center[1],
          data=stats['max_img'],
          distance=stats['distance'],
          pixel_size=pixel_size,
          saturated_value=saturated_value,
          timestamp=cspad_tbx.evt_timestamp(stats['time']),
          wavelength=stats['wavelength'])
        if self._max_out is not None:
          p = cspad_tbx.dwritef2(d, self._max_out)
        else:
          p = cspad_tbx.dwritef(d, self.max_dirname, self.max_basename)
        self.logger.info("Max written to %s" % p)

    if stats['nfail'] == 0:
      self.logger.info("%d images processed" % stats['nmemb'])
    else:
      self.logger.warning(
        "%d images processed, %d failed" % (stats['nmemb'], stats['nfail']))
