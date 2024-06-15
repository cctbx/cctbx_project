from __future__ import absolute_import, division, print_function

from xfel.merging.application.worker import worker
class simple_cell_listing(worker):
  '''A class that reads the experiments and writes the unit cell / space group list in text format.'''

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(simple_cell_listing, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Output unit cells in tdata format'

  def run(self, all_experiments, all_reflections):
    if self.mpi_helper.rank == 0:
      self.logger.log("Writing out a list of all unit cells")
  # START OUTPUT ALL UNIT CELLS
    all_results = []
    for experiment_id, experiment in enumerate(all_experiments):
          if experiment.identifier is None or len(experiment.identifier) == 0:
            experiment.identifier = create_experiment_identifier(experiment, experiments_filename, experiment_id)

          uc = experiment.crystal.get_unit_cell()
          sg = "".join(experiment.crystal.get_space_group().type().lookup_symbol().split())

          uc_ps = uc.parameters()
          result = "%f %f %f %f %f %f %s"%(uc_ps[0], uc_ps[1], uc_ps[2], uc_ps[3], uc_ps[4], uc_ps[5], sg)
          all_results.append(result)

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    self.global_results = comm.reduce(all_results, MPI.SUM, 0)
    if self.mpi_helper.rank == 0:
      file_cells = open("%s.tdata"%(self.params.tdata.output_path), "w")
      for result in self.global_results:
        line =  str(result) + "\n"
        file_cells.write(line)
      file_cells.close()
      self.logger.main_log("output a list of %d unit cells"%len(self.global_results))
  # END OUTPUT ALL UNIT CELLS

    return all_experiments, all_reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(simple_cell_listing)
