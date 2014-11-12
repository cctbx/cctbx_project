from __future__ import division
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.xtc_process
#
from psana import *
import os, sys
from libtbx.utils import Sorry
from libtbx.phil import parse

phil_scope = parse('''
  dispatch {
    max_events = None
      .type = int
      .help = "If not specified, process all events. Otherwise, only process this many"
  }
  input {
    cfg = None
      .type = str
      .help = "Path to psana config file"
    experiment = None
      .type = str
      .help = "Experiment identifier, e.g. cxi84914"
    run_num = None
      .type = int
      .help = "Run number or run range to process"
  }
''', process_includes=True)

class Script(object):
  """ Script to process XFEL data at LCLS """
  def __init__(self):
    pass

  def run(self):
    """ Process all images assigned to this thread """
    user_phil = []
    for arg in sys.argv[1:]:
      if (os.path.isfile(arg)):
        user_phil.append(parse(file_name=arg))
      else:
        try:
          user_phil.append(parse(arg))
        except RuntimeError, e:
          raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

    params = phil_scope.fetch(sources=user_phil).extract()

    assert params.input.cfg is not None
    assert params.input.experiment is not None
    assert params.input.run_num is not None

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
    size = comm.Get_size() # size: number of processes running in this job

    # set up psana
    setConfigFile(params.input.cfg)
    dataset_name = "exp=%s:run=%s:idx"%(params.input.experiment,params.input.run_num)
    ds = DataSource(dataset_name)

    # set this to sys.maxint to analyze all events
    if params.dispatch.max_events is None:
      max_events = sys.maxint
    else:
      max_events = params.dispatch.max_events

    for run in ds.runs():
      # list of all events
      times = run.times()
      nevents = min(len(times),max_events)
      mylength = nevents//size # easy but sloppy. lose few events at end of run.
      # chop the list into pieces, depending on rank
      mytimes= times[rank*mylength:(rank+1)*mylength]
      for i in range(mylength):
        evt = run.event(mytimes[i])
        id = evt.get(EventId)
        print "Event #",i," has id:",id

        # nop since the module does all the work

if __name__ == "__main__":
  script = Script()
  script.run()
