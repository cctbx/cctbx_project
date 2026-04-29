from __future__ import absolute_import, division, print_function
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.ensemble_refinement_pipeline
#

help_message = """
This is an MPI enabled pipeline to run several dials and cctbx.xfel commands in row.
Specifically dials.combine_experiments, dials.refine, cctbx.xfel.recompute_mosaicity,
and optionally cctbx.xfel.mpi_integrate
"""

from libtbx.phil import parse
from dials.util import show_mail_on_error
from libtbx.mpi4py import MPI, mpi_abort_on_exception

ensemble_refinement_pipline_str = """
combine_experiments_phil = None
  .type = path
  .help = Path to the phil file for dials.combine_experiments
refine_phil = None
  .type = path
  .help = Path to the phil file for dials.refine
recompute_mosaicity_phil = None
  .type = path
  .help = Path to the phil file for cctbx.xfel.recompute_mosaicity
integration_phil = None
  .type = path
  .help = Path to the phil file for cctbx.xfel.mpi_integrate
"""

phil_scope = parse(ensemble_refinement_pipline_str)

class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import ArgumentParser
    import libtbx.load_env

    # Create the option parser
    usage = "usage: %s combine_experiments_phil=combine.phil refine_phil=refine.phil recompute_mosaicity_phil=recompute_mosaicity.phil integration_phil=integration.phil" % libtbx.env.dispatcher_name
    self.parser = ArgumentParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message)

  @mpi_abort_on_exception
  def run(self):
    ''' Parse the options. '''
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
    size = comm.Get_size() # size: number of processes running in this job

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True if rank == 0 else False)

    if rank == 0:
      try:
        if params.combine_experiments_phil:
          from dials.command_line.combine_experiments import run
          run(args=[params.combine_experiments_phil])
        if params.refine_phil:
          from dials.command_line.refine import run
          run(args=[params.refine_phil])
        if params.recompute_mosaicity_phil:
          from xfel.command_line.recompute_mosaicity import Script as RecomputeScript
          script = RecomputeScript()
          script.run_with_preparsed(*script.parser.parse_args([params.recompute_mosaicity_phil], show_diff_phil=True))
      except Exception as e:
        if hasattr(e, 'message'):
          print(e.message)
        else:
          print(str(e))
        status = False
      else:
        status = True
    else:
      status = None
    status = comm.bcast(status, root=0)
    if not status:
      print("Rank %d shutting down due to job failure"%rank)
      return

    if params.integration_phil:
      import xfel.merging.command_line.mpi_integrate # reconfigure phils for integration
      from xfel.merging.command_line.merge import Script as IntegrateScript
      IntegrateScript().run(args=[params.integration_phil])

if __name__ == '__main__':
  with show_mail_on_error():
    script = Script()
    script.run()
