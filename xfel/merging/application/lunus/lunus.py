from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.imageset import ImageSetFactory
from dials.array_family import flex
from lunus.command_line.process import get_experiment_params, get_experiment_xvectors
import lunus as lunus_processor
from scitbx import matrix

def mpi_enabled():
  return ('MPI_LOCALRANKID' in os.environ.keys() or 'OMPI_COMM_WORLD_RANK' in os.environ.keys())
#  return False

def mpi_init():
  global mpi_comm
  global MPI
  from mpi4py import MPI as MPI
  mpi_comm = MPI.COMM_WORLD

def get_mpi_rank():
  return mpi_comm.Get_rank() if mpi_enabled() else 0


def get_mpi_size():
  return mpi_comm.Get_size() if mpi_enabled() else 1

def mpi_send(d,dest):
  if mpi_enabled():
    mpi_comm.Send(d,dest=dest)

def mpi_recv(d,source):
  if mpi_enabled():
    mpi_comm.Recv(d,source=source)

def mpi_bcast(d,root=0):
  if mpi_enabled():
    db = mpi_comm.bcast(d,root=root)
  else:
    db = d

  return db

def mpi_barrier():
  if mpi_enabled():
    mpi_comm.Barrier()

def mpi_sum(d, dt, root=0):
  if mpi_enabled():
    mpi_comm.Reduce(d,dt,op=MPI.SUM,root=root)
  else:
    dt = d

def mpi_reduce_p(p, root=0):
  if mpi_enabled():
#    if get_mpi_rank() == 0:
#      print("LUNUS.PROCESS: Convertinf flex arrays to numpy arrays")
#      sys.stdout.flush()

    l = p.get_lattice().as_numpy_array()
    c = p.get_counts().as_numpy_array()

    if get_mpi_rank() == root:
      lt = np.zeros_like(l)
      ct = np.zeros_like(c)
    else:
      lt = None
      ct = None


    mpi_comm.Reduce(l,lt,op=MPI.SUM,root=root)
    mpi_comm.Reduce(c,ct,op=MPI.SUM,root=root)

    if get_mpi_rank() == root:
      #logger.info("LUNUS.PROCESS: Converting numpy arrays to flex arrays")
      p.set_lattice(flex.double(lt))
      p.set_counts(flex.int(ct))

  return p



class lunus(worker):
  """
  Calls into the Lunus library to do diffuse scatter integration

  See DOI: 10.1007/978-1-59745-483-4_17
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(lunus, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

    self.current_path = None
    self.current_imageset = None

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

    if self.mpi_helper.rank == 0:
      self.reference_experiments = experiments[0:1]
    else:
      self.reference_experiments = None

    self.reference_experiments = self.mpi_helper.comm.bcast(self.reference_experiments, root=0)

    reference_experiment_params = get_experiment_params(experiments[0:1])
    self.processor = lunus_processor.Process(len(reference_experiment_params))
    with open(self.params.lunus.deck_file) as f:
      self.deck = f.read()

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

    if self.current_path != experiments[0].imageset.paths()[0]:
      self.current_imageset = ImageSetFactory.make_imageset(experiments[0].imageset.paths())
    idx = experiments[0].imageset.indices()[0]
    experiments[0].imageset = self.current_imageset[idx:idx+1]

    data = self.reference_experiments[0].imageset[0]
    if not isinstance(data, tuple):
      data = data,
    for panel_idx, panel in enumerate(data):
      self.processor.set_image(panel_idx, panel)

    for pidx in range(len(experiment_params)):
      deck_and_extras = self.deck+experiment_params[pidx]
      p.LunusSetparamsim(pidx,deck_and_extras)

    self.logger.log("LUNUS: Processing image")

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

      self.logger.log("LUNUS: Entering lprocimlt()")
      p.LunusProcimlt(0)
      self.logger.log("LUNUS: Done with lprocimlt()")
    else:
      crystal = experiments[0].crystal
      A_matrix = matrix.sqr(crystal.get_A()).inverse()
      At_flex = A_matrix.transpose().as_flex_double_matrix()

      p.set_amatrix(At_flex)

      p.LunusProcimlt(1)

  def finalize(self):
    mpi_barrier() # Need to synchronize at this point so that all the server/client ranks finish

    self.logger.log("STARTING FINALIZE, rank %d, size %d"%(get_mpi_rank(), get_mpi_size()))

    p = self.processor
    p = mpi_reduce_p(p)

    if get_mpi_rank() == 0:
      p.divide_by_counts()
      p.write_as_hkl('results.hkl')
      p.write_as_vtk('results.vtk')

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(lunus)
