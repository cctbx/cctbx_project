from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.merge

import os, sys

# vvv Temporary fix for https://github.com/dials/dials/issues/1998
# Remove this when https://github.com/pydata/numexpr/pull/400 and
# https://github.com/dials/dials/pull/2000 have been merged and released
if os.environ.get("CCTBX_NO_UUID", None) is not None:
  def mp_system(): return sys.platform.capitalize()
  def mp_machine(): return 'x86_64'
  import platform
  platform.system = mp_system
  platform.machine = mp_machine
# ^^^

from xfel.merging.application.mpi_helper import mpi_helper
from xfel.merging.application.mpi_logger import mpi_logger
from libtbx.mpi4py import mpi_abort_on_exception

# Note, when modifying this list, be sure to modify the README in xfel/merging

default_steps = [
  'input',
  'balance', # balance input load
  'model_scaling', # build full Miller list, model intensities, and resolution binner - for scaling and post-refinement
  'modify', # apply polarization correction, etc.
  'filter', # reject whole experiments or individual reflections
  'scale',
  'postrefine',
  'statistics_unitcell', # if required, save the average unit cell to the phil parameters
  'statistics_beam', # save the average wavelength to the phil parameters
  'model_statistics', # build full Miller list, model intensities, and resolution binner - for statistics. May use average unit cell.
  'statistics_resolution', # calculate resolution statistics for experiments
  'group', # re-distribute reflections over the ranks, so that all measurements of every HKL are gathered at the same rank
  'errors_merge', # correct errors using a per-HKL algorithm, e.g. errors_from_sample_residuals
  'statistics_intensity', # calculate resolution statistics for intensities
  'merge', # merge HKL intensities, MPI-gather all HKLs at rank 0, output "odd", "even" and "all" HKLs as mtz files
  'statistics_intensity_cxi', # run cxi_cc code ported from cxi-xmerge
]

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    self.mpi_helper = mpi_helper()
    self.mpi_logger = mpi_logger()

  def __del__(self):
    self.mpi_helper.finalize()

  def parse_input(self, args=None):
    '''Parse input at rank 0 and broadcast the input parameters and options to all ranks'''

    if self.mpi_helper.rank == 0:
      from xfel.merging.application.phil.phil import phil_scope
      help_message = '''Merge xfel data.'''

      # The script usage
      import libtbx.load_env
      self.usage = "usage: %s [options] [param.phil] " % libtbx.env.dispatcher_name
      self.parser = None

      '''Initialize the script.'''
      from dials.util.options import ArgumentParser
      # Create the parser
      self.parser = ArgumentParser(
        usage=self.usage,
        phil=phil_scope,
        epilog=help_message)

      # Parse the command line. quick_parse is required for MPI compatibility
      params, options = self.parser.parse_args(args, show_diff_phil=True,quick_parse=True)

      # Log the modified phil parameters
      diff_phil_str = self.parser.diff_phil.as_str()
      if diff_phil_str != "":
        self.mpi_logger.main_log("The following parameters have been modified:\n%s"%diff_phil_str)

      # prepare for transmitting input parameters to all ranks
      transmitted = dict(params = params, options = options)

      # make the output folders
      try:
        os.mkdir(params.output.output_dir)
      except FileExistsError:
        pass

    else:
      transmitted = None

    # broadcast parameters and options to all ranks
    self.mpi_logger.log("Broadcasting input parameters...")
    self.mpi_logger.log_step_time("BROADCAST_INPUT_PARAMS")

    transmitted = self.mpi_helper.comm.bcast(transmitted, root = 0)

    self.params = transmitted['params']
    self.options = transmitted['options']

    self.mpi_logger.set_log_file_paths(self.params)

    self.mpi_logger.log("Received input parameters and options")
    self.mpi_logger.log_step_time("BROADCAST_INPUT_PARAMS", True)

  @mpi_abort_on_exception
  def run(self, args=None):
    import datetime
    time_now = datetime.datetime.now()

    self.mpi_logger.log(str(time_now))
    if self.mpi_helper.rank == 0:
      self.mpi_logger.main_log(str(time_now))

    self.mpi_logger.log_step_time("TOTAL")

    self.mpi_logger.log_step_time("PARSE_INPUT_PARAMS")
    self.parse_input(args)
    self.mpi_logger.log_step_time("PARSE_INPUT_PARAMS", True)

    if self.params.mp.debug.cProfile:
      import cProfile
      pr = cProfile.Profile()
      pr.enable()

    # Create the workers using the factories
    self.mpi_logger.log_step_time("CREATE_WORKERS")
    from xfel.merging import application
    import importlib, copy

    self._resolve_persistent_columns()

    workers = []
    self.params.dispatch.step_list = self.params.dispatch.step_list or default_steps
    for step in self.params.dispatch.step_list:
      step_factory_name = step
      step_additional_info = []

      step_info = step.split('_')
      assert len(step_info) > 0
      if len(step_info) > 1:
        step_factory_name = step_info[0]
        step_additional_info = step_info[1:]

      try:
        factory = importlib.import_module('xfel.merging.application.' + step_factory_name + '.factory')
      except ModuleNotFoundError:
        custom_worker_path = os.environ.get('XFEL_CUSTOM_WORKER_PATH')
        if custom_worker_path is None: raise
        # remember the system path so the custom worker can temporarily modify it
        sys_path = copy.deepcopy(sys.path)
        pathstr = os.path.join(
            custom_worker_path, step_factory_name, 'factory.py'
        )

        modulename = 'xfel.merging.application.' + step_factory_name + '.factory'
        spec = importlib.util.spec_from_file_location(modulename, pathstr)
        factory = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(factory)
        # reset the path
        sys.path = sys_path

      workers.extend(factory.factory.from_parameters(self.params, step_additional_info, mpi_helper=self.mpi_helper, mpi_logger=self.mpi_logger))

    # Perform phil validation up front
    for worker in workers:
      worker.validate()
    self.mpi_logger.log_step_time("CREATE_WORKERS", True)

    # Do the work
    experiments = reflections = None
    step = 0
    while(workers):
      worker = workers.pop(0)
      self.mpi_logger.log_step_time("STEP_" + worker.__repr__())
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
      self.mpi_logger.log_step_time("STEP_" + worker.__repr__(), True)
      if experiments:
        self.mpi_logger.log("Ending step with %d experiments"%len(experiments))

    if self.params.output.save_experiments_and_reflections:
      if self.mpi_helper.size == 1:
        filename_suffix = ""
      else:
        filename_suffix = "_%06d"%self.mpi_helper.rank

      if len(reflections):
        assert 'id' in reflections
        reflections.as_file(os.path.join(self.params.output.output_dir, "%s%s.refl"%(self.params.output.prefix, filename_suffix)))
      if experiments:
        experiments.as_file(os.path.join(self.params.output.output_dir, "%s%s.expt"%(self.params.output.prefix, filename_suffix)))

    self.mpi_logger.log_step_time("TOTAL", True)

    if self.params.mp.debug.cProfile:
      pr.disable()
      pr.dump_stats(os.path.join(self.params.output.output_dir, "cpu_%s_%d.prof"%(self.params.output.prefix, self.mpi_helper.rank)))

  def _resolve_persistent_columns(self):
    if self.params.output.expanded_bookkeeping:
      if self.params.input.persistent_refl_cols is None:
        self.params.input.persistent_refl_cols = []
      keysCreatedByMerge = ["input_refl_index", "original_id", "file_list_mapping", "is_odd_experiment"]
      for key in keysCreatedByMerge:
        if key not in self.params.input.persistent_refl_cols:
          self.params.input.persistent_refl_cols.append(key)
    if "correlation_after_post" not in self.params.input.persistent_refl_cols:
      self.params.input.persistent_refl_cols.append("correlation_after_post")
    if "correlation" not in self.params.input.persistent_refl_cols:
      self.params.input.persistent_refl_cols.append("correlation")

if __name__ == '__main__':
  script = Script()

  result = script.run()
