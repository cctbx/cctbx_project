# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

from __future__ import absolute_import, division, print_function

import math
import multiprocessing
import numpy

from libtbx import easy_pickle
from scitbx.array_family import flex
import scitbx.math
from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import skip_event_flag

class average_mixin(common_mode.common_mode_correction):
  def __init__(self,
               address,
               avg_dirname=None,
               avg_basename=None,
               stddev_dirname=None,
               stddev_basename=None,
               max_dirname=None,
               max_basename=None,
               background_path=None,
               flags=None,
               hot_threshold=None,
               gain_threshold=None,
               noise_threshold=7,
               elastic_threshold=9,
               symnoise_threshold=4,
               **kwds):
    """
    @param address         Full data source address of the DAQ device
    @param avg_dirname     Directory portion of output average image
                           XXX mean
    @param avg_basename    Filename prefix of output average image XXX
                           mean
    @param flags inactive:  Eliminate the inactive pixels
                 noelastic: Eliminate elastic scattering
                 nohot:     Eliminate the hot pixels
                 nonoise:   Eliminate noisy pixels
                 symnoise:  Symmetrically eliminate noisy pixels
    @param stddev_dirname  Directory portion of output standard
                           deviation image XXX std
    @param stddev_basename Filename prefix of output standard
                           deviation image XXX std
    @param max_dirname     Directory portion of output maximum
                           projection image
    @param max_basename    Filename prefix of output maximum
                           projection image
     """

    super(average_mixin, self).__init__(
      address=address,
      **kwds
    )
    self.roi = None
    self.avg_basename = cspad_tbx.getOptString(avg_basename)
    self.avg_dirname = cspad_tbx.getOptString(avg_dirname)
    self.detector = cspad_tbx.address_split(address)[0]
    self.flags = cspad_tbx.getOptStrings(flags, default = [])
    self.stddev_basename = cspad_tbx.getOptString(stddev_basename)
    self.stddev_dirname = cspad_tbx.getOptString(stddev_dirname)
    self.max_basename = cspad_tbx.getOptString(max_basename)
    self.max_dirname = cspad_tbx.getOptString(max_dirname)
    self.background_path = cspad_tbx.getOptString(background_path)
    self.hot_threshold = cspad_tbx.getOptFloat(hot_threshold)
    self.gain_threshold = cspad_tbx.getOptFloat(gain_threshold)
    self.noise_threshold = cspad_tbx.getOptFloat(noise_threshold)
    self.elastic_threshold = cspad_tbx.getOptFloat(elastic_threshold)
    self.symnoise_threshold = cspad_tbx.getOptFloat(symnoise_threshold)

    if background_path is not None:
      background_dict = easy_pickle.load(background_path)
      self.background_img = background_dict['DATA']

    self._have_max = self.max_basename is not None or \
                     self.max_dirname is not None
    self._have_mean = self.avg_basename is not None or \
                      self.avg_dirname is not None
    self._have_std = self.stddev_basename is not None or \
                     self.stddev_dirname is not None

    # Start a server process which holds a set of Python objects that
    # other processes can manipulate using proxies.  The queues will
    # be used in endjob() to pass images between the worker processes,
    # and the lock will ensure the transfer is treated as a critical
    # section.  There is therefore the risk of a hang if the queues
    # cannot hold all the data one process will supply before another
    # empties it.
    #
    # In an attempt to alleviate this issue, separate queues are used
    # for the potentially big images.  The hope is to prevent
    # producers from blocking while consumers are locked out by using
    # more buffers.
    mgr = multiprocessing.Manager()
    self._lock = mgr.Lock()
    self._metadata = mgr.dict()
    self._queue_max = mgr.Queue()
    self._queue_sum = mgr.Queue()
    self._queue_ssq = mgr.Queue()


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(average_mixin, self).beginjob(evt, env)

    if self.dark_img is not None and self.hot_threshold is not None:
      self.hot_threshold *= flex.median(self.dark_img.as_1d())
      self.logger.info("HOT THRESHOLD: %.2f" % self.hot_threshold)
      self.logger.info("Number of pixels above hot threshold: %i" % \
                       (self.dark_img > self.hot_threshold).count(True))

    self._nfail = 0
    self._nmemb = 0

    # The time_base metadata item is a two-long array of seconds and
    # milliseconds, and it must be recorded only once.  It indicates
    # the base time, which is subtracted from all per-shot times.
    # Assuming an average run length of ten minutes, five minutes past
    # the start of a run is a good base time.
    self._lock.acquire()
    if 'time_base' not in self._metadata:
      self._metadata['time_base'] = (cspad_tbx.evt_time(evt)[0] + 5 * 60, 500)
    self._lock.release()


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(average_mixin, self).event(evt, env)
    if evt.get('skip_event'):
      return

    # Get the distance for the detectors that should have it, and set
    # it to NaN for those that should not.
    if self.detector == 'CxiDs1' or \
       self.detector == 'CxiDs2' or \
       self.detector == 'CxiDsd' or \
       self.detector == 'XppGon':
      distance = cspad_tbx.env_distance(self.address, env, self._detz_offset)
      if distance is None:
        self._nfail += 1
        self.logger.warning("event(): no distance, shot skipped")
        evt.put(skip_event_flag(), 'skip_event')
        return
    else:
      distance = float('nan')

    if ("skew" in self.flags):
      # Take out inactive pixels
      if self.roi is not None:
        pixels = self.cspad_img[self.roi[2]:self.roi[3], self.roi[0]:self.roi[1]]
        dark_mask = self.dark_mask[self.roi[2]:self.roi[3], self.roi[0]:self.roi[1]]
        pixels = pixels.as_1d().select(dark_mask.as_1d())
      else:
        pixels = self.cspad_img.as_1d().select(self.dark_mask.as_1d()).as_double()
      stats = scitbx.math.basic_statistics(pixels.as_double())
      #stats.show()
      self.logger.info("skew: %.3f" %stats.skew)
      self.logger.info("kurtosis: %.3f" %stats.kurtosis)
      if 0:
        from matplotlib import pyplot
        hist_min, hist_max = flex.min(flex_cspad_img.as_double()), flex.max(flex_cspad_img.as_double())
        print(hist_min, hist_max)
        n_slots = 100
        n, bins, patches = pyplot.hist(flex_cspad_img.as_1d().as_numpy_array(), bins=n_slots, range=(hist_min, hist_max))
        pyplot.show()

      # XXX This skew threshold probably needs fine-tuning
      skew_threshold = 0.35
      if stats.skew < skew_threshold:
        self._nfail += 1
        self.logger.warning("event(): skew < %f, shot skipped" % skew_threshold)
        evt.put(skip_event_flag(), 'skip_event')
        return
      #self.cspad_img *= stats.skew

    if ("inactive" in self.flags):
      self.cspad_img.set_selected(self.dark_stddev <= 0, 0)

    if ("noelastic" in self.flags):
      ELASTIC_THRESHOLD = self.elastic_threshold
      self.cspad_img.set_selected(self.cspad_img > ELASTIC_THRESHOLD, 0)

    if self.hot_threshold is not None:
      HOT_THRESHOLD = self.hot_threshold
      self.cspad_img.set_selected(self.dark_img > HOT_THRESHOLD, 0)

    if self.gain_map is not None and self.gain_threshold is not None:
      # XXX comparing each pixel to a moving average would probably be better
      # since the gain should vary approximately smoothly over different areas
      # of the detector
      GAIN_THRESHOLD = self.gain_threshold
      #self.logger.debug(
        #"rejecting: %i" %(self.gain_map > GAIN_THRESHOLD).count(True))
      self.cspad_img.set_selected(self.gain_map > GAIN_THRESHOLD, 0)

    if ("nonoise" in self.flags):
      NOISE_THRESHOLD = self.noise_threshold
      self.cspad_img.set_selected(self.cspad_img < NOISE_THRESHOLD, 0)

    if ("sigma_scaling" in self.flags):
      self.do_sigma_scaling()

    if ("symnoise" in self.flags):
      SYMNOISE_THRESHOLD = self.symnoise_threshold
      self.cspad_img.set_selected((-SYMNOISE_THRESHOLD < self.cspad_img) &
                                  ( self.cspad_img  < SYMNOISE_THRESHOLD), 0)

    if ("output" in self.flags):
      try:
        from six.moves import cPickle as pickle
      except ImportError:
        import pickle
      import os
      if (not os.path.isdir(self.pickle_dirname)):
        os.makedirs(self.pickle_dirname)
      flexdata = flex.int(self.cspad_img.astype(numpy.int32))
      d = cspad_tbx.dpack(
        address=self.address,
        data=flexdata,
        timestamp=cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt))
      )
      G = open(os.path.join(".",self.pickle_dirname)+"/"+self.pickle_basename,
               "ab")
      pickle.dump(d,G,pickle.HIGHEST_PROTOCOL)
      G.close()

    if self.photon_threshold is not None and self.two_photon_threshold is not None:
      self.do_photon_counting()

    if self.background_path is not None:
      self.cspad_img -= self.background_img


    # t and self._sum_time are a two-long arrays of seconds and
    # milliseconds which hold time with respect to the base time.
    t = [t1 - t2 for (t1, t2) in zip(cspad_tbx.evt_time(evt),
                                     self._metadata['time_base'])]
    if self._nmemb == 0:
      # The peers metadata item is a bit field where a bit is set if
      # the partial sum from the corresponding worker process is
      # pending.  If this is the first frame a worker process sees,
      # set its corresponding bit in the bit field since it will
      # contribute a partial sum.
      if env.subprocess() >= 0:
        self._lock.acquire()
        if 'peers' in self._metadata:
          self._metadata['peers'] |= (1 << env.subprocess())
        else:
          self._metadata['peers'] = (1 << env.subprocess())
        self._lock.release()

      self._sum_distance = distance
      self._sum_time = (t[0], t[1])
      self._sum_wavelength = self.wavelength

      if self._have_max:
        self._max_img = self.cspad_img.deep_copy()
      if self._have_mean:
        self._sum_img = self.cspad_img.deep_copy()
      if self._have_std:
        self._ssq_img = flex.pow2(self.cspad_img)

    else:
      self._sum_distance += distance
      self._sum_time = (self._sum_time[0] + t[0], self._sum_time[1] + t[1])
      self._sum_wavelength += self.wavelength

      if self._have_max:
        sel = (self.cspad_img > self._max_img).as_1d()
        self._max_img.as_1d().set_selected(
          sel, self.cspad_img.as_1d().select(sel))
      if self._have_mean:
        self._sum_img += self.cspad_img
      if self._have_std:
        self._ssq_img += flex.pow2(self.cspad_img)

    self._nmemb += 1


  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

  def endjob(self, obj1, obj2=None):
    """The endjob() function finalises the mean and standard deviation
    images.  The distance and wavelength in all images is actually the
    mean distance and wavelength, since standard deviations or maximum
    values of those quantities do not make much sense in
    visualisation.

    @param evt Event object (psana only)
    @param env Environment object
    @return    A dictionary object with accumulated statistics or @c
               none if the contribution from the worker process is
               accounted for elsewhere
    """

    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2
    from Queue import Empty

    super(average_mixin, self).endjob(env)

    # This entire function is protected by self._lock to guard against
    # race conditions.
    self._lock.acquire()

    # Attempt to get all the information from the shared objects,
    # without blocking.
    try:
      queue_max = self._queue_max.get_nowait()
      queue_sum = self._queue_sum.get_nowait()
      queue_ssq = self._queue_ssq.get_nowait()

      queue_distance = self._metadata['distance']
      queue_nfail = self._metadata['nfail']
      queue_nmemb = self._metadata['nmemb']
      queue_time = self._metadata['time']
      queue_wavelength = self._metadata['wavelength']
    except (Empty, KeyError):
      pass

    # If a complete set of items could be retrieved from the shared
    # objects, add them to this process's partial sums.  If only a
    # subset of the expected items could be retrieved from the shared
    # objects, log an error and proceed.
    items = [not self._have_max or 'queue_max' in locals(),
             not self._have_mean or 'queue_sum' in locals(),
             not self._have_std or 'queue_ssq' in locals(),
             'queue_distance' in locals(),
             'queue_nfail' in locals(),
             'queue_nmemb' in locals(),
             'queue_time' in locals(),
             'queue_wavelength' in locals()]
    if items.count(False) == 0:
      if self._have_max:
        if hasattr(self, '_max_img'):
          sel = (queue_max > self._max_img).as_1d()
          self._max_img.as_1d().set_selected(sel, queue_max.as_1d().select(sel))
        else:
          self._max_img = queue_max
      if self._have_mean:
        self._sum_img = getattr(self, '_sum_img', 0) + queue_sum
      if self._have_std:
        self._ssq_img = getattr(self, '_ssq_img', 0) + queue_ssq

      self._sum_distance = getattr(self, '_sum_distance', 0) + queue_distance
      self._nfail = getattr(self, '_nfail', 0) + queue_nfail
      self._nmemb = getattr(self, '_nmemb', 0) +  queue_nmemb
      self._sum_time = (getattr(self, '_sum_time', (0, 0))[0] + queue_time[0],
                        getattr(self, '_sum_time', (0, 0))[1] + queue_time[1])
      self._sum_wavelength = getattr(self, '_sum_wavelength', 0) + \
                             queue_wavelength

    elif items.count(True) > 0:
      self.logger.error("Queue holds incomplete set of data items")

    # Clear the bit field for the worker process.
    if env.subprocess() >= 0:
      self._metadata['peers'] &= ~(1 << env.subprocess())

    if self._metadata.get('peers', -1) > 0:
      # There are other processes left.  Place the accumulated sums
      # back in the shared objects.  If this ever blocks, the buffer
      # size of the queues will probably have to be increased.
      if self._nmemb > 0:
        self._metadata['distance'] = self._sum_distance
        self._metadata['nfail'] = self._nfail
        self._metadata['nmemb'] = self._nmemb
        self._metadata['time'] = self._sum_time
        self._metadata['wavelength'] = self._sum_wavelength

        if self._have_max:
          self._queue_max.put(self._max_img)
        if self._have_mean:
          self._queue_sum.put(self._sum_img)
        if self._have_std:
          self._queue_ssq.put(self._ssq_img)

      self._lock.release()
      return None

    # This is the last worker process, all others must have
    # contributed their partial sums.  Finalise the max, mean, and
    # standard deviation images if requested.
    d = {'nfail': self._nfail,
         'nmemb': self._nmemb}
    if self._nmemb > 0:
      d['distance'] = self._sum_distance / self._nmemb
      d['time'] = (
        self._metadata['time_base'][0] +
        int(round(self._sum_time[0] / self._nmemb)),
        self._metadata['time_base'][1] +
        int(round(self._sum_time[1] / self._nmemb)))
      d['wavelength'] = self._sum_wavelength / self._nmemb

      if self._have_max:
        d['max_img'] = self._max_img
      if self._have_mean or self._have_std:
        mean_img = self._sum_img / self._nmemb
      if self._have_mean:
        d['mean_img'] = mean_img
      if self._have_std:
        # Accumulating floating-point numbers introduces errors,
        # which may cause negative variances.  Since a two-pass
        # approach is unacceptable, the standard deviation is
        # clamped at zero.
        d['std_img'] = self._ssq_img - self._sum_img * mean_img
        d['std_img'].set_selected(d['std_img'] < 0, 0)
        if self._nmemb == 1:
          d['std_img'] = flex.sqrt(d['std_img'])
        else:
          d['std_img'] = flex.sqrt(d['std_img'] / (self._nmemb - 1))

    self._lock.release()
    return d
