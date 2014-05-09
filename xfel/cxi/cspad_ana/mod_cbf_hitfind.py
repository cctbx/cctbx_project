# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id

"""Spotfind on  an in-memory cbf image and call it a hit or miss
depending on input parameters
"""


from __future__ import division
from dials.algorithms.peak_finding.spotfinder_factory import SpotFinderFactory
from dials.framework.registry import Registry
from dxtbx.imageset import MemImageSet
from dxtbx.datablock import DataBlockFactory

__version__ = ""

from xfel.cxi.cspad_ana.mod_cspad_cbf import mod_cspad_cbf
from xfel.cxi.cspad_ana import cspad_tbx

import os
from libtbx.utils import Sorry

class mod_cbf_hitfind(mod_cspad_cbf):
  """Class for hitfinding on a cbf image
  """

  def __init__(self,
               out_dirname,
               out_basename,
               target_phil = None,
               **kwds):
    """The mod_cbf_hitfind class constructor stores the parameters passed from
    the pyana configuration file in instance variables.

    @param out_dirname  Directory portion of output image pathname
    @param out_basename Filename prefix of output image pathname
    @param target_phil  File with spotfinding parameters
    """

    super(mod_cbf_hitfind, self).__init__(**kwds)

    self._basename = cspad_tbx.getOptString(out_basename)
    self._dirname = cspad_tbx.getOptString(out_dirname)

    # get default parameters
    config = Registry().config()

    # load custom parameters
    target_path = cspad_tbx.getOptString(target_phil)
    if target_path is not None:
      if not os.path.isfile(target_path):
        raise Sorry("Target not found: " + target_path)
      from iotbx.phil import parse
      args = []
      args.append(parse(file_name=target_path,process_includes=True))
      config.fetch(sources = args)

    self.dials_phil = config.params()

    # create the spot finder
    self.spotfinder = SpotFinderFactory.from_parameters(self.dials_phil)

  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.  It
    outputs the detector image associated with the event @p evt to the
    file system.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_cbf_hitfind, self).event(evt, env)
    if (evt.get('skip_event')):
      return

    # spotfind
    imgset = MemImageSet([self.cspad_img])
    datablock = DataBlockFactory.from_imageset(imgset)[0]
    reflections = self.spotfinder(datablock)

    if len(reflections) < self.dials_phil.refinement.reflections.minimum_number_of_reflections:
      self.logger.info("Not enough spots to index")
      evt.put(True, "skip_event")
      return
