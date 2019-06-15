from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.merge
import sys
from xfel.merging.application.mpi_helper import mpi_helper
from xfel.merging.application.mpi_logger import mpi_logger
from six.moves import cStringIO as StringIO

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    self.mpi_helper = mpi_helper()
    self.mpi_logger = mpi_logger()

  def __del__(self):
    self.mpi_helper.finalize()

  def parse_input(self):
    '''Parse input at rank 0 and broadcast the input parameters and options to all ranks'''

    if self.mpi_helper.rank == 0:
      from xfel.merging.application.phil.phil import phil_scope
      help_message = '''Merge xfel data.'''

      # The script usage
      import libtbx.load_env
      self.usage = "usage: %s [options] [param.phil] " % libtbx.env.dispatcher_name
      self.parser = None

      '''Initialize the script.'''
      from dials.util.options import OptionParser
      # Create the parser
      self.parser = OptionParser(
        usage=self.usage,
        phil=phil_scope,
        epilog=help_message)

      # Parse the command line. quick_parse is required for MPI compatibility
      try:
        bkp = sys.stdout
        sys.stdout = out = StringIO()
        params, options = self.parser.parse_args(show_diff_phil=True,quick_parse=True)
        self.mpi_logger.log(out.getvalue())
        if self.mpi_helper.rank == 0:
          self.mpi_logger.main_log(out.getvalue())
      finally:
        sys.stdout = bkp

      # prepare for transmitting input parameters to all ranks
      self.mpi_logger.log("Broadcasting input parameters...")
      transmitted = dict(params = params, options = options)
    else:
      transmitted = None

    # broadcast parameters and options to all ranks
    self.mpi_logger.log_step_time("BROADCAST_INPUT_PARAMS")

    transmitted = self.mpi_helper.comm.bcast(transmitted, root = 0)

    self.params = transmitted['params']
    self.options = transmitted['options']

    self.mpi_logger.set_log_file_paths(self.params)

    self.mpi_logger.log("Received input parameters and options")
    self.mpi_logger.log_step_time("BROADCAST_INPUT_PARAMS", True)

  def run(self):

    import datetime
    time_now = datetime.datetime.now()

    self.mpi_logger.log(str(time_now))
    if self.mpi_helper.rank == 0:
      self.mpi_logger.main_log(str(time_now))


    self.mpi_logger.log_step_time("TOTAL")

    self.mpi_logger.log_step_time("PARSE_INPUT_PARAMS")
    self.parse_input()
    self.mpi_logger.log_step_time("PARSE_INPUT_PARAMS", True)

    # Create the workers using the factories
    from xfel.merging import application
    import importlib

    workers = []
    for step in ['input',
                 'model scaling', # the full miller set is based on the target unit cell
                 'modify', # polarization correction, etc.
                 'edit',   # add asu HKL column, remove unnecessary columns from reflection table
                 'filter', # unit cell, I/Sigma
                 'errors pre_merge', # e.g. ha14
                 'scale',
                 'postrefine',
                 'statistics unit_cell', # if required, saves the average unit cell to the phil parameters
                 'statistics beam', # saves the average wavelength to the phil parameters
                 'model statistics', # if required, the full miller set is based on the average unit cell
                 'statistics experiment_resolution',
                 'group', # MPI-alltoall: this must be done before any analysis or merging that requires all measurements of an HKL
                 'errors post_merge', # e.g. errors_from_sample_residuals
                 'statistics intensity',
                 'merge', # merge HKL intensities, MPI-gather all HKLs at rank 0, output "odd", "even" and "all" HKLs as mtz files
                 'statistics intensity cxi', # follows the merge step and uses cxi_cc code ported from cxi-xmerge
                 ]:

      step_factory_name = step
      step_additional_info = []

      step_info = step.split(' ')
      assert len(step_info) > 0
      if len(step_info) > 1:
        step_factory_name = step_info[0]
        step_additional_info = step_info[1:]

      factory = importlib.import_module('xfel.merging.application.' + step_factory_name + '.factory')
      workers.extend(factory.factory.from_parameters(self.params, step_additional_info))

    # Perform phil validation up front
    for worker in workers:
      worker.validate()

    # Do the work
    experiments = reflections = None
    step = 0
    while(workers):
      worker = workers.pop(0)

      # Log worker name, i.e. execution step name
      step += 1
      if step > 1:
        self.mpi_logger.log('')
      step_desc = "STEP %d: %s"%(step, worker)
      self.mpi_logger.log(step_desc)

      if self.mpi_helper.rank == 0:
        if step > 1:
          self.mpi_logger.main_log('')
        self.mpi_logger.main_log(step_desc)

      # Execute worker
      experiments, reflections = worker.run(experiments, reflections)

    self.mpi_logger.log_step_time("TOTAL", True)

if __name__ == '__main__':
  script = Script()

  result = script.run()
