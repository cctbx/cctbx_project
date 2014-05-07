# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id

from __future__ import division

import copy, multiprocessing, os

from scitbx.array_family import flex
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana.mod_cspad_cbf import mod_cspad_cbf

class mod_average_cbf(mod_cspad_cbf):
  def __init__(self,
               address,
               avg_dirname=None,
               avg_basename=None,
               stddev_dirname=None,
               stddev_basename=None,
               max_dirname=None,
               max_basename=None,
               **kwds):
    """
    @param address         Full data source address of the DAQ device
    @param avg_dirname     Directory portion of output average image
    @param avg_basename    Filename prefix of output average image
    @param stddev_dirname  Directory portion of output standard
                           deviation image
    @param stddev_basename Filename prefix of output standard
                           deviation image
    @param max_dirname     Directory portion of output maximum
                           projection image
    @param max_basename    Filename prefix of output maximum
                           projection image
     """

    super(mod_average_cbf, self).__init__(
      address=address,
      **kwds
    )
    self.avg_basename = cspad_tbx.getOptString(avg_basename)
    self.avg_dirname = cspad_tbx.getOptString(avg_dirname)
    self.detector = cspad_tbx.address_split(address)[0]
    self.stddev_basename = cspad_tbx.getOptString(stddev_basename)
    self.stddev_dirname = cspad_tbx.getOptString(stddev_dirname)
    self.max_basename = cspad_tbx.getOptString(max_basename)
    self.max_dirname = cspad_tbx.getOptString(max_dirname)

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

    super(mod_average_cbf, self).beginjob(evt, env)

    self._nfail = 0
    self._nmemb = 0

    # The time_base metadata item is a two-long array of seconds and
    # milliseconds, and it must be recorded only once.  It indicates
    # the base time, which is subtracted from all per-shot times.
    # Assuming an average run length of ten minutes, five minutes past
    # the start of a run is a good base time.
    self._lock.acquire()
    if 'time_base' not in self._metadata.keys():
      self._metadata['time_base'] = (cspad_tbx.evt_time(evt)[0] + 5 * 60, 500)
    self._lock.release()


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_average_cbf, self).event(evt, env)
    if evt.get('skip_event'):
      return

    # Get the distance for the detectors that should have it, and set
    # it to NaN for those that should not.
    if self.detector == 'CxiDs1' or \
       self.detector == 'CxiDsd' or \
       self.detector == 'XppGon':
      distance = cspad_tbx.env_distance(self.address, env, self._detz_offset)
      if distance is None:
        self._nfail += 1
        self.logger.warning("event(): no distance, shot skipped")
        evt.put(True, 'skip_event')
        return
    else:
      distance = float('nan')

    # t and self._sum_time are a two-long arrays of seconds and
    # milliseconds which hold time with respect to the base time.
    t = [t1 - t2 for (t1, t2) in zip(cspad_tbx.evt_time(evt),
                                     self._metadata['time_base'])]

    det = self.cspad_img.get_detector()
    data = [self.cspad_img.get_raw_data(i) for i in xrange(len(det))]
    if self._nmemb == 0:
      # The peers metadata item is a bit field where a bit is set if
      # the partial sum from the corresponding worker process is
      # pending.  If this is the first frame a worker process sees,
      # set its corresponding bit in the bit field since it will
      # contribute a partial sum.
      if env.subprocess() >= 0:
        self._lock.acquire()
        if 'peers' in self._metadata.keys():
          self._metadata['peers'] |= (1 << env.subprocess())
        else:
          self._metadata['peers'] = (1 << env.subprocess())
        self._lock.release()

      self._sum_distance = distance
      self._sum_time = (t[0], t[1])
      self._sum_wavelength = self.wavelength

      if self._have_max:
        self._max_img = copy.deepcopy(data)
      if self._have_mean:
        self._sum_img = copy.deepcopy(data)
      if self._have_std:
        self._ssq_img = [flex.pow2(d) for d in data]

    else:
      self._sum_distance += distance
      self._sum_time = (self._sum_time[0] + t[0], self._sum_time[1] + t[1])
      self._sum_wavelength += self.wavelength

      if self._have_max:
        sel = [(d > max_d).as_1d() for d, max_d in zip(data, self._max_img)]
        for d, max_d, s in zip(data, self._max_img, sel):
          max_d.as_1d().set_selected(s, d.as_1d().select(s))
      if self._have_mean:
        self._sum_img = [sum_d + d for sum_d, d in zip(self._sum_img, data)]
      if self._have_std:
        self._ssq_img = [ssq_d + flex.pow2(d) for ssq_d, d in zip(self._ssq_img, data)]

    self._nmemb += 1


  def endjob(self, env):
    """The endjob() function finalises the mean and standard deviation
    images.  The distance and wavelength in all images is actually the
    mean distance and wavelength, since standard deviations or maximum
    values of those quantities do not make much sense in
    visualisation.

    @param env Environment object
    @return    A dictionary object with accumulated statistics or @c
               none if the contribution from the worker process is
               accounted for elsewhere
    """

    from Queue import Empty

    super(mod_average_cbf, self).endjob(env)

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
    except Empty, KeyError:
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
          sel = [(queue_max_d > max_d).as_1d() for queue_max_d, max_d in zip(queue_max, self._max_img)]
          for max_d, s, queue_max_d in zip(self._max_img, sel, queue_max):
            max_d.as_1d().set_selected(s, queue_max_d.as_1d().select(s))
        else:
          self._max_img = queue_max
      if self._have_mean:
        self._sum_img = [sum_d + queue_sum_d for sum_d, queue_sum_d in zip(getattr(self, '_sum_img', [0]*len(queue_sum)), queue_sum)]
      if self._have_std:
        self._ssq_img = [ssq_d + queue_ssq_d for ssq_d, queue_ssq_d in zip(getattr(self, '_ssq_img', [0]*len(queue_ssq)), queue_ssq)]

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
      return

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
      d['time'] = cspad_tbx.evt_timestamp(d['time'])
      d['wavelength'] = self._sum_wavelength / self._nmemb

      if self._have_max:
        d['max_img'] = self._max_img
      if self._have_mean or self._have_std:
        mean_img = [sum_d.as_double() / self._nmemb for sum_d in self._sum_img]
      if self._have_mean:
        d['mean_img'] = mean_img
        #d['mean_img'] = [mean_d.iround() for mean_d in mean_img]
      if self._have_std:
        # Accumulating floating-point numbers introduces errors,
        # which may cause negative variances.  Since a two-pass
        # approach is unacceptable, the standard deviation is
        # clamped at zero.
        d['std_img'] = [ssq_d.as_double() - sum_d.as_double() * mean_d for ssq_d, sum_d, mean_d in zip(self._ssq_img, self._sum_img, mean_img)]

        for stddev_d in d['std_img']:
          stddev_d.set_selected(stddev_d < 0, 0)

        if self._nmemb == 1:
          d['std_img'] = [flex.sqrt(stddev_d) for stddev_d in d['std_img']]
        else:
          d['std_img'] = [flex.sqrt(stddev_d / (self._nmemb - 1)) for stddev_d in d['std_img']]

    self._lock.release()

    #beam_center = self.beam_center
    if d['nmemb'] > 0:
      from xfel.cftbx.detector.cspad_cbf_tbx import cbf_file_to_basis_dict, write_cspad_cbf
      metro = cbf_file_to_basis_dict(self.metrology)
      detector = self.base_dxtbx.get_detector()

      def make_tiles(data, detector):
        """
        Assemble a tiles dictionary as required by write_cspad_cbf
        Assumes the order in the data array matches the order of the enumerated detector panels
        """
        assert len(data) == 64
        tiles = {}

        tile_id = 0
        for p_id, p in enumerate(detector.hierarchy()):
          for s_id, s in enumerate(p):
            for a_id, a in enumerate(s):
              data[tile_id].resize(flex.grid(a.get_image_size()[1], a.get_image_size()[0]))
              tiles[(0, p_id, s_id, a_id)] = data[tile_id]
              tile_id += 1
        return tiles

      def make_path(dirname, basename):
        """
        Check that the directories are available, then assemble the path
        """
        if basename is None:
          basename = ""
        if dirname is None:
          dirname = "."
        if not os.path.isdir(dirname):
          os.makedirs(dirname)

        return os.path.join(dirname, basename + '.cbf')

      if self.avg_dirname  is not None or \
        self.avg_basename is not None:
        tiles = make_tiles(d['mean_img'], detector)
        path = make_path(self.avg_dirname, self.avg_basename)
        write_cspad_cbf(tiles, metro, 'cbf', d['time'], path, d['wavelength'], d['distance'])

        self.logger.info("Average written to %s" % path)

      if self.stddev_dirname  is not None or \
         self.stddev_basename is not None:
        tiles = make_tiles(d['std_img'], detector)
        path = make_path(self.stddev_dirname, self.stddev_basename)
        write_cspad_cbf(tiles, metro, 'cbf', d['time'], path, d['wavelength'], d['distance'])
        self.logger.info("Standard deviation written to %s" % path)

      if self.max_dirname  is not None or \
         self.max_basename is not None:
        tiles = make_tiles(d['max_img'], detector)
        path = make_path(self.max_dirname, self.max_basename)
        write_cspad_cbf(tiles, metro, 'cbf', d['time'], path, d['wavelength'], d['distance'])
        self.logger.info("Max written to %s" % path)

    if d['nfail'] == 0:
      self.logger.info("%d images processed" % d['nmemb'])
    else:
      self.logger.warning(
        "%d images processed, %d failed" % (d['nmemb'], d['nfail']))

