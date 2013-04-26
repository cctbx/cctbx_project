# -*- mode: python; indent-tabs-mode: nil; python-indent: 2; tab-width: 8 -*-
#
# $Id$

"""XXX"""

from __future__ import division

import os
import tarfile
import scipy.misc

import cspad_tbx


class mod_pulnix_dump(object):
  """Class for dumping simple intensity images into a tarball.  At 120
  Hz it is inconvenient to dump every single image to its own file.
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

    self._address = cspad_tbx.getOptString(address)
    self._out_basename = cspad_tbx.getOptString(out_basename)
    self._out_dirname = cspad_tbx.getOptString(out_dirname)


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    if not os.path.isdir(self._out_dirname):
      os.makedirs(self._out_dirname)

    # Open per-process tar files to avoid clobbering while
    # multiprocessing.
    if env.subprocess() >= 0:
      basename = '%sp%02d-r%04d.tar' % (
        self._out_basename, env.subprocess(), evt.run())
    else:
      basename = '%sr%04d.tar' % (
        self._out_basename, evt.run())

    self._tar = tarfile.open(os.path.join(self._out_dirname, basename), 'w')


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    seqno = cspad_tbx.evt_seqno(evt)
    if seqno is None:
      #print "Shot skipped"
      return

    value = evt.getFrameValue(self._address)
    if value is None:
      #print "Shot skipped"
      return

    # Write image to temporary image file.  Add image file to tar
    # archive, and remove temporary file.
    path = os.path.join(
      self._out_dirname, '%s%s.png' % (self._out_basename, seqno))
    im = scipy.misc.toimage(value.data())
    im.save(path)

    self._tar.add(path, arcname=os.path.basename(path))
    os.remove(path)
    #print "*** Saved image to", path


  def endjob(self, env):
    """The endjob() function closes the tar file opened in beginjob().

    @param env Environment object
    """

    self._tar.close()
