# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id: mod_cspad_cbf_dump.py

"""Output image to the file system.
"""


from __future__ import division

__version__ = ""

from  xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana.mod_cspad_cbf import mod_cspad_cbf
import os, pycbf

class mod_cspad_cbf_dump(mod_cspad_cbf):
  """Class for outputting cbf images to the file system within the pyana
  analysis framework.
  """

  def __init__(self,
               address,
               out_dirname,
               out_basename,
               **kwds):
    """The mod_cspad_cbf_dump class constructor stores the parameters passed from
    the pyana configuration file in instance variables.

    @param address      Full data source address of the DAQ device
    @param out_dirname  Directory portion of output image pathname
    @param out_basename Filename prefix of output image pathname
    """

    super(mod_cspad_cbf_dump, self).__init__(address=address, **kwds)

    self._basename = cspad_tbx.getOptString(out_basename)
    self._dirname = cspad_tbx.getOptString(out_dirname)


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.  It
    outputs the detector image associated with the event @p evt to the
    file system.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_cspad_cbf_dump, self).event(evt, env)
    if (evt.get('skip_event')):
      return

    self.cspad_img.sync_detector_to_cbf()

    t = self.timestamp
    s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]

    dest_path = os.path.join(self._dirname, self._basename + s + ".cbf")

    self.cspad_img._cbf_handle.write_widefile(dest_path, pycbf.CBF,\
      pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K, 0)

    self.logger.info("Wrote %s"%os.path.basename(dest_path))
