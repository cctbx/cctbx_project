# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.stripe_experiment_wrapper

import sys
from dials.command_line.combine_experiments import run as combine_run
from dials.command_line.refine import run as refine_run
from xfel.command_line.recompute_mosaicity import run as recompute_mosaicity_run
from xfel.merging.command_line.mpi_integrate import run as reintegration_run
from xfel.command_line.frame_extractor import run as frame_extractor_run
from libtbx.phil import parse
from libtbx.utils import Sorry

"""
Thin wrapper for time depending ensemble refinement from Brewster 2018

Executes the following programs in order:
dials.combine_experiments
dials.refine
cctbx.xfel.recompute_mosaicity
cctbx.xfel.mpi_integrate
cctbx.xfel.frame_extractor

Also manages MPI, by only using mpi for mpi_integrate, using rank 0 for
the other programs
"""

phil_scope = parse("""
combine_experiments_phil = None
  .type = str
refinement_phil = None
  .type = str
recompute_mosaicity_phil = None
  .type = str
reintegration_phil = None
  .type = str
postprocessing_phil = None
  .type = str
""")

def run(args):
  user_phil = []
  for arg in args:
    try:
      user_phil.append(parse(arg))
    except Exception as e:
      raise Sorry("Unrecognized argument: %s" % arg)

  user_scope, unused = phil_scope.fetch(sources=user_phil, track_unused_definitions=True)
  if any(unused):
    msg = "\n".join([str(loc) for loc in unused])
    raise Sorry("Unrecognized argument(s): " + msg)

  params = user_scope.extract()

  try:
    from libtbx.mpi4py import MPI
  except ImportError:
    rank = 0
    size = 1
  else:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()  # each process in MPI has a unique id, 0-indexed
    size = comm.Get_size()  # size: number of processes running in this job

  if rank == 0:
    combine_run([params.combine_experiments_phil])
    refine_run([params.refinement_phil])
    recompute_mosaicity_run([params.recompute_mosaicity_phil])
  
  if params.reintegration_phil:
    reintegration_run([params.reintegration_phil])

  if rank == 0 and params.postprocessing_phil:
    frame_extractor_run([params.postprocessing_phil])

if __name__ == "__main__":
  run(sys.argv[1:])
