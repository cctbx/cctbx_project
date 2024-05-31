from __future__ import absolute_import, division, print_function
from six.moves import range
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.xtc_process
#
from psana import *
import os, sys, socket
from libtbx.utils import Sorry
from libtbx.phil import parse
from xfel.cxi.cspad_ana import cspad_tbx

help_str = """
Use cxi.xtc_process to analyze XTC streams from LCLS using psana modules
specified in a config file.

You will likely find cxi.mpi_submit more useful generally, but use this prog-
ram to test your config files for a few events.

Should be ran from your myrelease directory. Be sure to run sit_setup first.

Example usages:

cxi.xtc_process input.cfg=cxi49812/thermo.cfg input.experiment=cxi49812 \\
  input.run_num=25

This will use one process on the current node to analyze every event from run
25 of experiment cxi49812 using the modules specfied in cxi49812/thermo.cfg.

mpirun -n 16 cxi.xtc_process input.cfg=cxi49812/thermo.cfg \\
  input.experiment=cxi49812 input.run_num=25 dispatch.max_events=1000

As above, but use 16 processes on the current node and only process 1000 events.

bsub -a mympi -n 100 -o out.log -q psanaq cxi.xtc_process \\
  input.cfg=cxi49812/thermo.cfg input.experiment=cxi49812 input.run_num=25

Submit a processing job to the psana queue using 100 cores and put the
resultant log in out.log.

"""

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
    stream = None
      .type = int
      .expert_level = 2
      .help = Stream number to read from. Usually not necessary as psana will read the data \
              from all streams by default
    use_ffb = False
      .type = bool
      .help = "Run on the ffb if possible. Only for active users!"
    xtc_dir = None
      .type = str
      .help = "Optional path to data directory if it's non-standard. Only needed if xtc"
      .help = "streams are not in the standard location for your PSDM installation."
  }
  output {
    output_dir = "."
      .type = str
      .help = "Directory for output files"
  }
  mp {
    method = *mpi sge
      .type = choice
      .help = "Muliprocessing method"
    mpi {
      method = *client_server striping
        .type = choice
        .help = Method of serving data to child processes in MPI. client_server:    \
                use one process as a server that sends timestamps to each process.  \
                All processes will stay busy at all times at the cost of MPI send/  \
                recieve overhead. striping: each process uses its rank to determine \
                which events to process. Some processes will finish early and go    \
                idle, but no MPI overhead is incurred.
    }
  }
  debug {
    write_debug_files = False
      .type = bool
      .help = "If true, will write out a tiny diagnostic file for each event before"
      .help = "and after the event is processed"
    use_debug_files = False
      .type = bool
      .help = "If true, will look for debug diagnostic files in the output_dir and"
      .help = "only process events that crashed previously"
    event_timestamp = None
      .type = str
      .help = "If set to be a timestamp, will only process the event that matches it"
  }
''', process_includes=True)

class Script(object):
  """ Script to process XFEL data at LCLS """
  def __init__(self):
    pass

  def run(self):
    """ Process all images assigned to this thread """
    if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv or "-c" in sys.argv:
      print(help_str)
      print("Showing phil parameters:")
      print(phil_scope.as_str(attributes_level = 2))
      return

    user_phil = []
    for arg in sys.argv[1:]:
      if (os.path.isfile(arg)):
        user_phil.append(parse(file_name=arg))
      else:
        try:
          user_phil.append(parse(arg))
        except RuntimeError as e:
          raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

    params = phil_scope.fetch(sources=user_phil).extract()
    self.params = params

    assert params.input.cfg is not None
    assert params.input.experiment is not None
    assert params.input.run_num is not None

    print("Processing run %d of experiment %s using config file %s"%(params.input.run_num, params.input.experiment, params.input.cfg))

    if params.mp.method == "mpi":
      from libtbx.mpi4py import MPI
      comm = MPI.COMM_WORLD
      rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
      size = comm.Get_size() # size: number of processes running in this job
    elif params.mp.method == "sge" and \
        'SGE_TASK_ID'    in os.environ and \
        'SGE_TASK_FIRST' in os.environ and \
        'SGE_TASK_LAST'  in os.environ:
      if 'SGE_STEP_SIZE' in os.environ:
        assert int(os.environ['SGE_STEP_SIZE']) == 1
      if os.environ['SGE_TASK_ID'] == 'undefined' or os.environ['SGE_TASK_ID'] == 'undefined' or os.environ['SGE_TASK_ID'] == 'undefined':
        rank = 0
        size = 1
      else:
        rank = int(os.environ['SGE_TASK_ID']) - int(os.environ['SGE_TASK_FIRST'])
        size = int(os.environ['SGE_TASK_LAST']) - int(os.environ['SGE_TASK_FIRST']) + 1
    else:
      rank = 0
      size = 1

    # set up psana
    setConfigFile(params.input.cfg)
    dataset_name = "exp=%s:run=%s:idx"%(params.input.experiment,params.input.run_num)
    if params.input.xtc_dir is not None:
      if params.input.use_ffb:
        raise Sorry("Cannot specify the xtc_dir and use SLAC's ffb system")
      dataset_name += ":dir=%s"%params.input.xtc_dir
    elif params.input.use_ffb:
      # as ffb is only at SLAC, ok to hardcode /reg/d here
      dataset_name += ":dir=/reg/d/ffb/%s/%s/xtc"%(params.input.experiment[0:3],params.input.experiment)
    if params.input.stream is not None:
      dataset_name += ":stream=%d"%params.input.stream
    ds = DataSource(dataset_name)

    # set this to sys.maxint to analyze all events
    if params.dispatch.max_events is None:
      max_events = sys.maxsize
    else:
      max_events = params.dispatch.max_events

    for run in ds.runs():
      # list of all events
      times = run.times()
      nevents = min(len(times),max_events)

      if params.mp.method == "mpi" and size > 2 and params.mp.mpi.method == 'client_server':
        # use a client/server approach to be sure every process is busy as much as possible
        # only do this if there are more than 2 processes, as one process will be a server
        print("Using MPI client server")
        if rank == 0:
          # server process
          for t in times[:nevents]:
            # a client process will indicate it's ready by sending its rank
            rankreq = comm.recv(source=MPI.ANY_SOURCE)
            comm.send(t,dest=rankreq)
          # send a stop command to each process
          for rankreq in range(size-1):
            rankreq = comm.recv(source=MPI.ANY_SOURCE)
            comm.send('endrun',dest=rankreq)
        else:
          # client process
          while True:
            # inform the server this process is ready for an event
            comm.send(rank,dest=0)
            evttime = comm.recv(source=0)
            if evttime == 'endrun': break
            self.process_event(run, evttime)
      else:
        # chop the list into pieces, depending on rank.  This assigns each process
        # events such that the get every Nth event where N is the number of processes
        print("Striping events")
        mytimes = [times[i] for i in range(nevents) if (i+rank)%size == 0]

        for i in range(len(mytimes)):
          self.process_event(run, mytimes[i])
      run.end()
    ds.end()

  def process_event(self, run, timestamp):
    """
    Process a single event from a run
    @param run psana run object
    @param timestamp psana timestamp object
    """

    ts = cspad_tbx.evt_timestamp((timestamp.seconds(),timestamp.nanoseconds()/1e6))
    if self.params.debug.event_timestamp is not None and self.params.debug.event_timestamp != ts:
      return

    ts_path = os.path.join(self.params.output.output_dir, "debug-" + ts + ".txt")

    if self.params.debug.use_debug_files:
      if not os.path.exists(ts_path):
        print("Skipping event %s: no debug file found"%ts)
        return

      f = open(ts_path, "r")
      if len(f.readlines()) > 1:
        print("Skipping event %s: processed successfully previously"%ts)
        return
      f.close()
      print("Accepted", ts)

    if self.params.debug.write_debug_files:
      f = open(ts_path, "w")
      f.write("%s about to process %s\n"%(socket.gethostname(), ts))
      f.close()

    # load the event.  This will cause any modules to be run.
    evt = run.event(timestamp)

    if self.params.debug.write_debug_files:
      f = open(ts_path, "a")
      f.write("done\n")
      f.close()

    id = evt.get(EventId)

    # nop since the module does all the work.

if __name__ == "__main__":
  script = Script()
  script.run()
