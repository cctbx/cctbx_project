# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id: mod_cbf_index.py

"""Index an in-memory cbf image
"""


from __future__ import division
from dials.algorithms.peak_finding.spotfinder_factory import SpotFinderFactory
from dials.util.options import SystemConfig
from dxtbx.imageset import MemImageSet

__version__ = ""

from xfel.cxi.cspad_ana.mod_cspad_cbf import mod_cspad_cbf
from xfel.cxi.cspad_ana import cspad_tbx

import os
from libtbx.utils import Sorry

class mod_cbf_index(mod_cspad_cbf):
  """Class for indexing a cbf image
  """

  def __init__(self,
               out_dirname,
               out_basename,
               target_phil = None,
               **kwds):
    """The mod_cbf_index class constructor stores the parameters passed from
    the pyana configuration file in instance variables.

    @param out_dirname  Directory portion of output image pathname
    @param out_basename Filename prefix of output image pathname
    """

    super(mod_cbf_index, self).__init__(**kwds)

    self._basename = cspad_tbx.getOptString(out_basename)
    self._dirname = cspad_tbx.getOptString(out_dirname)

    # get default parameters
    sysconfig = SystemConfig()
    params = sysconfig.config()

    # load custom parameters
    target_path = cspad_tbx.getOptString(target_phil)
    if target_path is None:
      self.dials_phil = params.extract()
    else:
      if not os.path.isfile(target_path):
        raise Sorry("Target not found: " + target_path)
      from iotbx.phil import parse
      source = parse(file_name=target_path,process_includes=True)
      self.dials_phil = params.fetch(source = source).extract()

    # create the spot finder
    self.spotfinder = SpotFinderFactory.from_parameters(self.dials_phil)

  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.  It
    outputs the detector image associated with the event @p evt to the
    file system.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_cbf_index, self).event(evt, env)
    if (evt.get('skip_event')):
      return

    # spotfind
    imgset = MemImageSet([self.cspad_img])
    reflections = self.spotfinder(imgset)
    if len(reflections) < 16:
      self.logger.info("Not enough spots into index")
      evt.put(True, "skip_event")
      return

    # index
    from dials_regression.indexing_test_data.i04_weak_data.run_indexing_api import run as run_index
    try:
      run_index([imgset,reflections])
      self.logger.info("Successfully indexed")
    except Exception, e:
      self.logger.info("Couldn't index, " + e.message)
      evt.put(True, "skip_event")
      return

