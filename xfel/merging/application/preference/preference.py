from __future__ import absolute_import, division, print_function

from functools import wraps

import numpy as np

from xfel.merging.application.worker import worker
from xfel.util.preference import DirectSpaceBases, \
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
      self.logger.log(message_)
      self.logger.log('')
    abc_stack = DirectSpaceBases.from_expts(experiments, sg)
    abc_stack = abc_stack.symmetrize(sg.build_derived_point_group())
    abc_stacks = self.mpi_helper.comm.gather(abc_stack)
    if self.mpi_helper.comm.rank != 0:
      return experiments, reflections
    abc_stack = DirectSpaceBases(np.concatenate(abc_stacks, axis=0))
    distributions = find_preferential_distribution(abc_stack, sg)
    self.logger.log('')
    self.logger.log(distributions.table)
    self.logger.log('')
    return experiments, reflections
