# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#

"""The mod_event_code module creates an "alist" for cxi.merge based on what events integrated
(using DIALS or LABELIT backends), given a set of psana evr event codes.
"""
from __future__ import absolute_import, division, print_function

import logging, os

from xfel.cxi.cspad_ana import cspad_tbx
import psana

class mod_event_code(object):
  def __init__(
    self,
    integration_dirname,
    out_dirname,
    event_codes,
    alist_names):
    """
    @param integration_dirname directory with integration pickle files
    @param out_dirname directory for alists
    @event_codes space-seperated list of evr event codes
    @alist_names corresponding list of names to give the alists
    """

    self.logger = logging.getLogger(self.__class__.__name__)
    self.logger.setLevel(logging.INFO)

    self.integration_dirname  = cspad_tbx.getOptString(integration_dirname)
    self.out_dirname          = cspad_tbx.getOptString(out_dirname)
    self.event_codes          = cspad_tbx.getOptStrings(event_codes)
    self.alist_names          = cspad_tbx.getOptStrings(alist_names)

    self.event_codes = [int(s) for s in self.event_codes]

    # Try to guess multiprocessing method
    if 'SGE_TASK_ID'    in os.environ and \
       'SGE_TASK_FIRST' in os.environ and \
       'SGE_TASK_LAST'  in os.environ:
      if 'SGE_STEP_SIZE' in os.environ:
        assert int(os.environ['SGE_STEP_SIZE']) == 1
      if os.environ['SGE_TASK_ID'] == 'undefined' or os.environ['SGE_TASK_ID'] == 'undefined' or os.environ['SGE_TASK_ID'] == 'undefined':
        self.rank = 0
        self.size = 1
      else:
        self.rank = int(os.environ['SGE_TASK_ID']) - int(os.environ['SGE_TASK_FIRST'])
        self.size = int(os.environ['SGE_TASK_LAST']) - int(os.environ['SGE_TASK_FIRST']) + 1
    else:
      try:
        from mpi4py import MPI
      except ImportError:
        self.rank = 0
        self.size = 1
      else:
        comm = MPI.COMM_WORLD
        self.rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
        self.size = comm.Get_size() # size: number of processes running in this job

    # Save a dicitonary of timestamps that satisfy any of the event codes given
    self.timestamps_d = {}
    for alist in self.alist_names:
      self.timestamps_d[alist] = []

  def __del__(self):
    logging.shutdown()

  def beginjob(self, evt, env):
    pass

  def event(self, evt, env):
    """The event() function puts a "skip_event" object with value @c
    True into the event if the shot is to be skipped.

    Read the evr codes for this event and save the timestamp if it has an evr code we are
    looking for.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    if (evt.get("skip_event")):
      return

    evr = psana.Detector('evr0')
    codes = evr.eventCodes(evt)
    timestamp = cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt)) # human readable format
    logging.info("mod_event_code, rank: %d, timestamp: %s, code_list: %s"%(self.rank, timestamp, ",".join(["%d"%c for c in codes])))

    # Save this timestamp if it has an event_code w are looking for
    for alist, code in zip(self.alist_names, self.event_codes):
      if code in codes:
        self.timestamps_d[alist].append(timestamp)

  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

  def endjob(self, obj1, obj2=None):
    """
    Write the alist files, seperated by rank. We do it in the endjob method so this module can be called
    before or after indexing and integration occurs during processing.

    @param evt Event object (psana only)
    @param env Environment object
    """

    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2

    for alist in self.alist_names:
      alist_f = open(os.path.join(self.out_dirname, "alist_%s_rank_%d.txt"%(alist, self.rank)), 'w')
      for t in self.timestamps_d[alist]:
        s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]

        if os.path.exists(os.path.join(self.integration_dirname, "int-0-%s.pickle"%s)): # DIALS output
          alist_f.write(os.path.join(self.integration_dirname, "int-0-%s.pickle\n"%s))
        elif os.path.exists(os.path.join(self.integration_dirname, "int-%s_00000.pickle"%t)): # LABELIT output
          alist_f.write(os.path.join(self.integration_dirname, "int-%s_00000.pickle\n"%t))
      alist_f.close()
