from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dxtbx.model.experiment_list import ExperimentList
from dials.array_family import flex

class lunus(worker):
  """
  Calls into the Lunus library to do diffuse scatter integration

  See DOI: 10.1007/978-1-59745-483-4_17
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(lunus, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Compute diffuse scatter using Lunus'

  def filter_by_n_laticces(self, experiments, reflections, n = 1):
    # filter out experiments with more than one lattice
    image_list = {}
    for expt in experiments:
      assert len(expt.imageset.paths()) == 1
      path = expt.imageset.paths()[0]
      if path not in image_list:
        image_list[path] = {}
      assert len(expt.imageset.indices()) == 1
      index = expt.imageset.indices()[0]
      if index in image_list[path]:
        image_list[path][index] += 1
      else:
        image_list[path][index] = 1
    all_image_lists = self.mpi_helper.comm.gather(image_list)
    if self.mpi_helper.rank == 0:
      all_image_list = {}
      for ilist in all_image_lists:
        for path in ilist:
          if path not in all_image_list:
            all_image_list[path] = {}
          for index in ilist[path]:
            if index in all_image_list[path]:
              all_image_list[path][index] += ilist[path][index]
            else:
              all_image_list[path][index] = ilist[path][index]
    else:
      all_image_list = None
    image_list = self.mpi_helper.comm.bcast(all_image_list, root=0)

    filtered_expts = ExperimentList()
    refls_sel = flex.bool(len(reflections), False)
    for expt_id, expt in enumerate(experiments):
      path = expt.imageset.paths()[0]
      index = expt.imageset.indices()[0]
      if image_list[path][index] == n:
        filtered_expts.append(expt)
        refls_sel |= reflections['exp_id'] == expt.identifier
    self.logger.log("Filtered out %d experiments with more than %d lattices out of %d"%((len(experiments)-len(filtered_expts), n, len(experiments))))
    return filtered_expts, reflections.select(refls_sel)

  def run(self, experiments, reflections):

    self.logger.log_step_time("LUNUS")

    experiments, reflections = self.filter_by_n_laticces(experiments, reflections)

    self.logger.log_step_time("LUNUS", True)

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(lunus)
