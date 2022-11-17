from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.imageset import ImageSetFactory
from dials.array_family import flex
from lunus.command_line.process import get_experiment_params, get_experiment_xvectors
import lunus as lunus_processor
from scitbx import matrix
import numpy as np

class lunus(worker):
  """
  Calls into the Lunus library to do diffuse scatter integration

  See DOI: 10.1007/978-1-59745-483-4_17
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(lunus, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

    self.current_path = None
    self.current_imageset = None
#    mpi_init()

  def __repr__(self):
    return 'Process diffuse scattering using Lunus, yielding a 3D dataset'

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
        refls_sel |= reflections['id'] == expt_id
    filtered_refls = reflections.select(refls_sel)
    filtered_refls.reset_ids()
    self.logger.log("Filtered out %d experiments with more than %d lattices out of %d"%((len(experiments)-len(filtered_expts), n, len(experiments))))
    return filtered_expts, filtered_refls

  def run(self, experiments, reflections):

    self.logger.log_step_time("LUNUS")

    experiments, reflections = self.filter_by_n_laticces(experiments, reflections)

    if self.mpi_helper.rank == 0:
      self.reference_experiments = experiments[0:1]
      self.logger.log("LUNUS: Using reference experiment %s with image %s %d" % (self.reference_experiments[0].identifier,self.reference_experiments[0].imageset.paths()[0],self.reference_experiments[0].imageset.indices()[0]))
    else:
      self.reference_experiments = None

    self.reference_experiments = self.mpi_helper.comm.bcast(self.reference_experiments, root=0)

    reference_experiment_params = get_experiment_params(self.reference_experiments[0:1])
    self.processor = lunus_processor.Process(len(reference_experiment_params))

    if self.mpi_helper.rank == 0:
      with open(self.params.lunus.deck_file) as f:
        self.deck = f.read()
    else:
      self.deck = None

    self.deck = self.mpi_helper.comm.bcast(self.deck, root=0)

    deck_and_extras = self.deck+reference_experiment_params[0]
    self.processor.LunusSetparamslt(deck_and_extras)

    self.lunus_integrate(self.reference_experiments, is_reference = True)

    for i in range(len(experiments)):
      self.lunus_integrate(experiments[i:i+1])

    self.finalize()

    self.logger.log_step_time("LUNUS", True)

    return experiments, reflections

  def lunus_integrate(self, experiments, is_reference = False):
    assert len(experiments) == 1

    experiment_params = get_experiment_params(experiments)
    p = self.processor

#    self.logger.log("LUNUS_INTEGRATE: Passed %s %d" % (experiments[0].imageset.paths()[0],experiments[0].imageset.indices()[0]))

    if self.current_path != experiments[0].imageset.paths()[0]:
      self.current_imageset = ImageSetFactory.make_imageset(experiments[0].imageset.paths())
    idx = experiments[0].imageset.indices()[0]
    experiments[0].imageset = self.current_imageset[idx:idx+1]

    self.logger.log("LUNUS_INTEGRATE: Processing image %s %d" % (experiments[0].imageset.paths()[0],experiments[0].imageset.indices()[0]))

    data = experiments[0].imageset[0]
    if not isinstance(data, tuple):
      data = data,
    for panel_idx, panel in enumerate(data):
      self.processor.set_image(panel_idx, panel)
#      self.logger.log("LUNUS_INTEGRATE: file %s panel_idx %d panel[0:10] = %s " % (experiments[0].imageset.paths()[0],panel_idx,str(list(panel[0:10]))))

    for pidx in range(len(experiment_params)):
      deck_and_extras = self.deck+experiment_params[pidx]
      p.LunusSetparamsim(pidx,deck_and_extras)

#    self.logger.log("LUNUS: Processing image")

    if is_reference:
      x = get_experiment_xvectors(experiments)
      for pidx in range(len(x)):
        p.set_xvectors(pidx,x[pidx])

      # We need an amatrix for the next call, to set up the lattice size
      if self.params.merging.set_average_unit_cell:
        assert 'average_unit_cell' in (self.params.statistics).__dict__
        uc = self.params.statistics.__phil_get__('average_unit_cell')
      else:
        uc = self.params.scaling.unit_cell

      assert uc is not None, "Lunus needs a target unit cell"
      A_matrix = matrix.sqr(uc.orthogonalization_matrix())
      At_flex = A_matrix.transpose().as_flex_double_matrix()

      p.set_amatrix(At_flex)

      p.LunusProcimlt(0)
    else:
      crystal = experiments[0].crystal
      A_matrix = matrix.sqr(crystal.get_A()).inverse()
      At_flex = A_matrix.transpose().as_flex_double_matrix()

      p.set_amatrix(At_flex)

      p.LunusProcimlt(1)

  def finalize(self):
    self.mpi_helper.comm.Barrier() # Need to synchronize at this point so that all the server/client ranks finish

    self.logger.log("LUNUS: STARTING FINALIZE, rank %d, size %d"%(self.mpi_helper.rank, self.mpi_helper.size))

    p = self.processor
    p = self.mpi_reduce_p(p)

    if self.mpi_helper.rank == 0:
      p.divide_by_counts()
      p.write_as_hkl('results.hkl')
      p.write_as_vtk('results.vtk')

  def mpi_reduce_p(self,p, root=0):
    l = p.get_lattice().as_numpy_array()
    c = p.get_counts().as_numpy_array()

    if self.mpi_helper.rank == root:
      lt = np.zeros_like(l)
      ct = np.zeros_like(c)
    else:
      lt = None
      ct = None

    self.mpi_helper.comm.Reduce(l,lt,op=self.mpi_helper.MPI.SUM,root=root)
    self.mpi_helper.comm.Reduce(c,ct,op=self.mpi_helper.MPI.SUM,root=root)

    if self.mpi_helper.rank == root:
      p.set_lattice(flex.double(lt))
      p.set_counts(flex.int(ct))

    return p

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(lunus)
