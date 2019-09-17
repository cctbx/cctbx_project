from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from cctbx.crystal import symmetry
from libtbx import Auto

class experiment_filter(worker):
  '''Reject experiments based on various criteria'''

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(experiment_filter, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Filter experiments'

  def check_unit_cell(self, experiment):

    experiment_unit_cell = experiment.crystal.get_unit_cell()

    is_ok = experiment_unit_cell.is_similar_to(self.params.filter.unit_cell.value.target_unit_cell,
                                               self.params.filter.unit_cell.value.relative_length_tolerance,
                                               self.params.filter.unit_cell.value.absolute_angle_tolerance)
    return is_ok

  def check_space_group(self, experiment):
    # build patterson group from the target space group
    target_unit_cell = self.params.filter.unit_cell.value.target_unit_cell
    target_space_group_info = self.params.filter.unit_cell.value.target_space_group
    target_symmetry = symmetry(unit_cell=target_unit_cell, space_group_info=target_space_group_info)
    target_space_group = target_symmetry.space_group()
    target_patterson_group_sn = target_space_group.build_derived_patterson_group().info().symbol_and_number()

    # build patterson group from the experiment space group
    experiment_space_group = experiment.crystal.get_space_group()
    experiment_patterson_group_sn = experiment_space_group.build_derived_patterson_group().info().symbol_and_number()

    is_ok = (target_patterson_group_sn == experiment_patterson_group_sn)

    return is_ok

  @staticmethod
  def remove_experiments(experiments, reflections, experiment_ids_to_remove):
    '''Remove specified experiments from the experiment list. Remove corresponding reflections from the reflection table'''

    new_experiments = ExperimentList()
    new_reflections = flex.reflection_table()

    for experiment in experiments:
      if experiment.identifier in experiment_ids_to_remove:
        continue
      new_experiments.append(experiment)
      refls = reflections.select(reflections['exp_id'] == experiment.identifier)
      new_reflections.extend(refls)

    return new_experiments, new_reflections

  def run(self, experiments, reflections):
    if 'unit_cell' not in self.params.filter.algorithm: # so far only "unit_cell" algorithm is supported
      return experiments, reflections

    self.logger.log_step_time("FILTER_EXPERIMENTS")

    # If the filter unit cell and/or space group params are Auto, use the corresponding scaling targets.
    if self.params.filter.unit_cell.value.target_unit_cell == Auto:
      self.params.filter.unit_cell.value.target_unit_cell = self.params.scaling.unit_cell
    if self.params.filter.unit_cell.value.target_space_group == Auto:
      self.params.filter.unit_cell.value.target_space_group = self.params.scaling.space_group

    self.logger.log("Using filter target unit cell: %s"%str(self.params.filter.unit_cell.value.target_unit_cell))
    self.logger.log("Using filter target space group: %s"%str(self.params.filter.unit_cell.value.target_space_group))

    experiment_ids_to_remove = []
    removed_for_unit_cell = 0
    removed_for_space_group = 0
    for experiment in experiments:
      if not self.check_space_group(experiment):
        experiment_ids_to_remove.append(experiment.identifier)
        removed_for_space_group += 1
      elif not self.check_unit_cell(experiment):
        experiment_ids_to_remove.append(experiment.identifier)
        removed_for_unit_cell += 1

    new_experiments, new_reflections = experiment_filter.remove_experiments(experiments, reflections, experiment_ids_to_remove)

    removed_reflections = len(reflections) - len(new_reflections)
    assert removed_for_space_group + removed_for_unit_cell == len(experiments) - len(new_experiments)

    self.logger.log("Experiments rejected because of unit cell dimensions: %d"%removed_for_unit_cell)
    self.logger.log("Experiments rejected because of space group %d"%removed_for_space_group)
    self.logger.log("Reflections rejected because of rejected experiments: %d"%removed_reflections)

    # MPI-reduce total counts
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    total_removed_for_unit_cell  = comm.reduce(removed_for_unit_cell, MPI.SUM, 0)
    total_removed_for_space_group  = comm.reduce(removed_for_space_group, MPI.SUM, 0)
    total_reflections_removed  = comm.reduce(removed_reflections, MPI.SUM, 0)

    # rank 0: log total counts
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Total experiments rejected because of unit cell dimensions: %d"%total_removed_for_unit_cell)
      self.logger.main_log("Total experiments rejected because of space group %d"%total_removed_for_space_group)
      self.logger.main_log("Total reflections rejected because of rejected experiments %d"%total_reflections_removed)

    self.logger.log_step_time("FILTER_EXPERIMENTS", True)

    # Do we have any data left?
    from xfel.merging.application.utils.data_counter import data_counter
    data_counter(self.params).count(new_experiments, new_reflections)

    return new_experiments, new_reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(experiment_filter)
