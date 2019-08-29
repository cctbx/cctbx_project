# -*- mode: python; indent-tabs-mode: nil; python-indent: 2; tab-width: 8 -*-
#
# $Id$

"""XXX"""

from __future__ import absolute_import, division, print_function

import os
import logging

from . import cspad_tbx


class mod_pulnix_dump(object):
  """Class for dumping simple intensity images into an uncompressed
  tarball.  At 120 Hz it is inconvenient to dump every single image to
  its own file.
  """

  def __init__(self,
               address,
               out_dirname=None,
               out_basename=None):
    """The mod_pulnix_dump class constructor stores the parameters passed
    from the pyana configuration file in instance variables.

    @param address      Address string XXX Que?!
    @param out_dirname  Directory portion of output image pathname
    @param out_basename Filename prefix of output image pathname
    """

    self._logger = logging.getLogger(self.__class__.__name__)
    self._logger.setLevel(logging.WARNING)

    self._address = cspad_tbx.getOptString(address)
    self._out_basename = cspad_tbx.getOptString(out_basename)
    self._out_dirname = cspad_tbx.getOptString(out_dirname)


  def beginjob(self, evt, env):
    """The beginjob() function is called at an XTC configure transition.
    It opens per-process tar files to avoid clobbering while
    multiprocessing.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    from collections import deque
    from tarfile import open as opentar

    # The maxlen parameter specifies for how many shots the direction
    # of the principal component must be stable before it can be
    # considered a clear view of the capillary.
    self._angles = deque(maxlen=10)

    if not os.path.isdir(self._out_dirname):
      os.makedirs(self._out_dirname)

    if env.subprocess() >= 0:
      basename = '%sp%02d-r%04d.tar' % (
        self._out_basename, env.subprocess(), evt.run())
    else:
      basename = '%sr%04d.tar' % (
        self._out_basename, evt.run())

    self._tar = opentar(os.path.join(self._out_dirname, basename), 'w')


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    from math import atan2
    from numpy import cov, mean, roll, std
    from numpy.linalg import eig
    from scipy.misc import toimage
    from tempfile import mkstemp

    value = evt.getFrameValue(self._address)
    if value is None:
      self._logger.info("event(): no image data, shot skipped")
      return

    seqno = cspad_tbx.evt_seqno(evt)
    if seqno is None:
      self._logger.warning("event(): no sequence number, shot skipped")
      return

    # Heuristic thresholds to determine whether image is blank or not.
    # XXX Better moved to a separate filtering module?

    # Normalize the image, compute its principal components, and sort
    # them by descending modulus of the corresponding eigenvalue.
    img_norm = (value.data() - mean(value.data())) / std(value.data())
    (w, v) = eig(cov(img_norm))
    if abs(w[0]) < abs(w[1]):
      w = roll(w, 1)
      v = roll(v, 1, 1)

    # For a clear view of the capillary the direction of the principal
    # compenents should be stable over time.  A standard deviation of
    # 0.1 radians corresponds to 5.7 degrees.  Three standard
    # deviations around the mean account for almost 100% of the area
    # under a Gaussian curve.
    #
    # XXX This works perfectly for r0047 and r0210 from cxi78513, but
    # has problems with r0090, r0098, r0099, r0109, r0197, r0267,
    # r0273, r0274, r0281, r0287, and probably others.
    visible = False
    self._angles.append(atan2(v[0, 0], v[0, 1]))
    if len(self._angles) == self._angles.maxlen:
      a_mean = mean(self._angles)
      a_std = std(self._angles)
      if a_std < 0.1 and abs(self._angles[-1] - a_mean) < 3 * a_std:
        visible = True

    if not visible:
      self._logger.info("event(): image blank, shot skipped")
      return

    # Write image to temporary image file.  Add image file to tar
    # archive, and remove temporary file.
    (tmp_fd, tmp_path) = mkstemp()
    tmp_stream = os.fdopen(tmp_fd, 'wb')
    img = toimage(value.data())
    img.save(tmp_stream, format='PNG')
    tmp_stream.close()

    arcname = '%s%s.png' % (self._out_basename, seqno)
    self._tar.add(tmp_path, arcname=arcname)
    os.unlink(tmp_path)

    self._logger.info("event(): archived %s" % arcname)

  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

  def endjob(self, obj1, obj2=None):
    """The endjob() function closes the tar file opened in beginjob().

    @param evt Event object (psana only)
    @param env Environment object
    """

    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2

    self._tar.close()
