from __future__ import division
from xfel.merging.application.worker import worker
from rstbx.symmetry.constraints import parameter_reduction
from cctbx import sgtbx
from cctbx import miller
import copy

def reindex_experiments(experiments, cb_op, space_group):
  for experiment in experiments:
    # debug print("BEFORE",experiment.crystal.get_unit_cell())
    experiment.crystal.set_space_group(sgtbx.space_group("P 1"))
    # On the experiment side we must apply the change of basis
    experiment.crystal = experiment.crystal.change_basis(cb_op)
    # Reset the space group
    experiment.crystal.set_space_group(space_group)
    # Re-symmetrize the unit cell parameters with the new space group constraints
    S = parameter_reduction.symmetrize_reduce_enlarge(
        experiment.crystal.get_space_group()
    )
    S.set_orientation(experiment.crystal.get_B())
    S.symmetrize()
    # And reset the orientation matrix with the re-symmetrized value
    experiment.crystal.set_B(S.orientation.reciprocal_matrix())
    # debug print("AFTER",experiment.crystal.get_unit_cell(),"\n")
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
    from xfel.merging.application.utils.data_counter import data_counter
    if self.mpi_helper.rank == 0: self.logger.main_log("Data count before reindexing")
    data_counter(self.params).count(experiments, reflections)

    # Empty ranks must not return early: data_counter.count() below is an
    # MPI collective and every rank has to participate or the job hangs.
    if len(experiments) > 0:
      cb_op_str = self.params.modify.reindex_to_abc.change_of_basis_op
      change_of_basis_op = sgtbx.change_of_basis_op(cb_op_str)

      space_group = self.params.modify.reindex_to_abc.space_group
      assert space_group is not None
      reindex_target_space_group_type = space_group.type()
      space_group = space_group.group()

      # On the reflection side we must reindex the original miller index
      experiments = reindex_experiments(
          experiments, change_of_basis_op, space_group=space_group
      )

      miller_indices = reflections["miller_index"]
      miller_indices_reindexed = change_of_basis_op.apply(miller_indices)
      reflections["miller_index"] = miller_indices_reindexed

      # And recalculate the new asymmetric unit indices consistent with the new space group
      reflections['miller_index_asymmetric'] = copy.deepcopy(reflections['miller_index'])
      miller.map_to_asu(reindex_target_space_group_type,
                        not self.params.merging.merge_anomalous,
                        reflections['miller_index_asymmetric'])

    if self.mpi_helper.rank == 0: self.logger.main_log("Data count after reindexing")
    data_counter(self.params).count(experiments, reflections)
    return experiments, reflections
