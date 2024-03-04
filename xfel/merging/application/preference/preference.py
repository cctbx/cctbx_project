from __future__ import absolute_import, division, print_function

from functools import wraps

import numpy as np

from xfel.merging.application.worker import worker
from xfel.util.preference import ascii_plot, DirectSpaceBases, \
  find_preferential_distribution, space_group_auto


def return_input_expts_refls_if_error(run_method):
  @wraps(run_method)
  def run_method_wrapper(worker_instance, experiments, reflections):
    try:
      return run_method(worker_instance, experiments, reflections)
    except Exception:  # noqa - this is meant to catch everything
      worker_instance.logger.log('Error when running worker, returning input')
      return experiments, reflections
  return run_method_wrapper


class PreferenceWorker(worker):
  """Read lattice vectors, apply symmetry, look for preferential orientation"""

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    kwargs = dict(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)
    super(PreferenceWorker, self).__init__(**kwargs)

  def __repr__(self):
    return 'Evaluate preferential orientation'

  @return_input_expts_refls_if_error
  def run(self, experiments, reflections):
    sg, message_ = space_group_auto(experiments, comm=self.mpi_helper.comm)
    if self.mpi_helper.rank == 0:
      for message_line in message_.split('\n'):
        self.logger.main_log(message_line)
      self.logger.main_log('')
    abc_stack = DirectSpaceBases.from_expts(experiments, sg)
    abc_stack = abc_stack.symmetrize(sg.build_derived_point_group())
    abc_stacks = self.mpi_helper.comm.gather(abc_stack)
    if self.mpi_helper.comm.rank != 0:
      return experiments, reflections
    abc_stack = DirectSpaceBases(np.concatenate(abc_stacks, axis=0))
    distributions = find_preferential_distribution(abc_stack, sg)
    self.logger.main_log(distributions.table)
    self.logger.main_log('')
    pref_direction, pref_distribution = distributions.best
    self.logger.main_log('Vector distribution for zone axes family '
                         + pref_direction)
    self.logger.main_log(ascii_plot(pref_distribution.vectors))
    self.logger.main_log('')
    return experiments, reflections
