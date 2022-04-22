from __future__ import division
from xfel.merging.application.worker import worker
from rstbx.symmetry.constraints import parameter_reduction
from cctbx import sgtbx

def reindex_experiments(experiments, cb_op, space_group):
  for experiment in experiments:
    experiment.crystal.set_space_group(sgtbx.space_group("P 1"))
    experiment.crystal = experiment.crystal.change_basis(cb_op)
    experiment.crystal.set_space_group(space_group)
    S = parameter_reduction.symmetrize_reduce_enlarge(
        experiment.crystal.get_space_group()
    )
    S.set_orientation(experiment.crystal.get_B())
    S.symmetrize()
    experiment.crystal.set_B(S.orientation.reciprocal_matrix())
  return experiments

class reindex_to_abc(worker):
  """
  Reindex according to an a,b,c basis tranformation
  """

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(reindex_to_abc, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Reindex according to an a,b,c basis tranformation'

  def run(self, experiments, reflections):
    if len(experiments) == 0:
      return experiments, reflections
    cb_op_str = self.params.modify.reindex_to_abc.change_of_basis_op
    change_of_basis_op = sgtbx.change_of_basis_op(cb_op_str)

    space_group = self.params.modify.reindex_to_abc.space_group
    assert space_group is not None
    space_group = space_group.group()

    experiments = reindex_experiments(
        experiments, change_of_basis_op, space_group=space_group
    )

    miller_indices = reflections["miller_index"]
    miller_indices_reindexed = change_of_basis_op.apply(miller_indices)
    reflections["miller_index"] = miller_indices_reindexed

    return experiments, reflections
