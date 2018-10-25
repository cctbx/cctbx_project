from __future__ import absolute_import, division, print_function

if __name__ == "__main__":
  import argparse

  parser = argparse.ArgumentParser(
    description = "Server process handling far end of remote execution",
    usage = "%(prog)s job-factory queue-factory (do not run directly!)",
    )
  parser.add_argument(
    "job",
    help = "Job factory pickle string",
    )
  parser.add_argument(
    "queue",
    help = "Job factory pickle string",
    )
  parser.add_argument( "--folder", default = ".", help = "Process workdir" )
  parser.add_argument(
    "--polltime",
    type = int,
    default = 0.01,
    help = "Poll intervall",
    )

  params = parser.parse_args()

  import os
  import sys
  os.chdir( params.folder )
  sys.path.append( os.path.abspath( params.folder ) )

  from libtbx.queuing_system_utils import remote

  jfact = remote.argument_to_object( arg = params.job )
  qfact = remote.argument_to_object( arg = params.queue )

  from libtbx.queuing_system_utils import scheduling

  manager = scheduling.Scheduler(
    handler = scheduling.Unlimited(
      factory = lambda:
        scheduling.ExecutionUnit(
          factory = jfact,
          processor = scheduling.RetrieveProcessor( queue = qfact(), timeout = 1 ),
          ),
      ),
    polling_interval = params.polltime,
    )

  server = remote.SchedulerServer(
    instream = sys.stdin,
    outstream = sys.stdout,
    manager = manager,
    )

  # Redirect standard filehandles, so that the communication stream is intact
  sys.stdout = sys.stderr

  server.serve()

