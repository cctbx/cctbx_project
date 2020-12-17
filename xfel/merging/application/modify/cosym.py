from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from xfel.merging.application.utils.memory_usage import get_memory_usage
from cctbx.array_family import flex
from cctbx import sgtbx

class cosym(worker):
  """
  Resolve indexing ambiguity using dials.cosym
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(cosym, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Resolve indexing ambiguity using dials.cosym'

  def run(self, experiments, reflections):

    self.logger.log_step_time("COSYM")
    if False:
      import pickle
      with open("special.pickle","wb") as F:
        pickle.dump((experiments,reflections), F)

      return experiments,reflections
    if True:
      import pickle
      with open("special.pickle","rb") as F:
        experiments,reflections = pickle.load(F)

    all_sampling_experiments = experiments
    all_sampling_reflections = reflections
    # because cosym has a problem with hashed identifiers, use simple experiment identifiers
    from dxtbx.model.experiment_list import ExperimentList
    sampling_experiments_for_cosym = ExperimentList()
    sampling_reflections_for_cosym = [] # is a list of flex.reflection_table

    def task_a():
      # add an anchor
      if self.params.modify.cosym.anchor:
        from xfel.merging.application.model.crystal_model import crystal_model
        XM = crystal_model(params = self.params, purpose="cosym")
        model_intensities = XM.run([],[])
        from dxtbx.model import Experiment, Crystal
        from scitbx.matrix import sqr
        O = sqr(model_intensities.unit_cell().orthogonalization_matrix()).transpose().elems
        real_a = (O[0],O[1],O[2])
        real_b = (O[3],O[4],O[5])
        real_c = (O[6],O[7],O[8])
        nc = Crystal(real_a,real_b,real_c, model_intensities.space_group())
        sampling_experiments_for_cosym.append(Experiment(crystal=nc)) # prepends the reference model to the cosym E-list
        from dials.array_family import flex

        exp_reflections = flex.reflection_table()
        exp_reflections['intensity.sum.value'] = model_intensities.data()
        exp_reflections['intensity.sum.variance'] = flex.pow(model_intensities.sigmas(),2)
        exp_reflections['miller_index'] = model_intensities.indices()
        exp_reflections['miller_index_asymmetric'] = model_intensities.indices()

        # prepare individual reflection tables for each experiment

        simple_experiment_id = len(sampling_experiments_for_cosym) - 1
        #experiment.identifier = "%d"%simple_experiment_id
        sampling_experiments_for_cosym[-1].identifier = "%d"%simple_experiment_id
        # experiment identifier must be a string according to *.h file
        # the identifier is changed on the _for_cosym Experiment list, not the master experiments for through analysis

        exp_reflections['id'] = flex.int(len(exp_reflections), simple_experiment_id)
        # register the integer id as a new column in the per-experiment reflection table

        exp_reflections.experiment_identifiers()[simple_experiment_id] = sampling_experiments_for_cosym[-1].identifier
        #apparently the reflection table holds a map from integer id (reflection table) to string id (experiment)

        sampling_reflections_for_cosym.append(exp_reflections)


    if self.mpi_helper.rank == 0:
      task_a()

    def task_1():

      for experiment in all_sampling_experiments:
        sampling_experiments_for_cosym.append(experiment)

        exp_reflections = all_sampling_reflections.select(all_sampling_reflections['exp_id'] == experiment.identifier)
        # prepare individual reflection tables for each experiment

        simple_experiment_id = len(sampling_experiments_for_cosym) - 1
        #experiment.identifier = "%d"%simple_experiment_id
        sampling_experiments_for_cosym[-1].identifier = "%d"%simple_experiment_id
        # experiment identifier must be a string according to *.h file
        # the identifier is changed on the _for_cosym Experiment list, not the master experiments for through analysis

        exp_reflections['id'] = flex.int(len(exp_reflections), simple_experiment_id)
        # register the integer id as a new column in the per-experiment reflection table

        exp_reflections.experiment_identifiers()[simple_experiment_id] = sampling_experiments_for_cosym[-1].identifier
        #apparently the reflection table holds a map from integer id (reflection table) to string id (experiment)

        sampling_reflections_for_cosym.append(exp_reflections)

      from dials.command_line import cosym as cosym_module
      cosym_module.logger = self.logger

      from dials.command_line.cosym import cosym
      COSYM = cosym(sampling_experiments_for_cosym, sampling_reflections_for_cosym, params=self.params.modify.cosym)
      return COSYM
    COSYM = task_1()

    rank_N_refl=flex.double([r.size() for r in COSYM.reflections])
    message = """Task 1. Prepare the data for cosym
    change_of_basis_ops_to_minimum_cell
    eliminate_sys_absent
    transform models into Miller arrays, putting data in primitive triclinic reduced cell
    There are %d experiments with %d reflections, averaging %.1f reflections/experiment"""%(
      len(COSYM.experiments), flex.sum(rank_N_refl), flex.mean(rank_N_refl))
    self.logger.log(message)
    if self.mpi_helper.rank == 0:
      self.logger.main_log(message)

    COSYM.run()
    exit()
    if self.params.modify.cosym.dataframe:
      from collections import OrderedDict
      assert len(experiments) == len(COSYM._experiments)
      keyval = [("experiment",[]),("reindex_op", [])]
      raw = OrderedDict(keyval)
      for sidx in range(len(experiments)):
        raw["experiment"].append(experiments[sidx].identifier)
        raw["reindex_op"].append(sgtbx.change_of_basis_op(COSYM.cosym_analysis.reindexing_ops[sidx][0]).as_hkl())
      keys = list(raw.keys())
      print (raw)
      from pandas import DataFrame as df
      data = df(raw)
      data.to_pickle(path = "cosym_myfile") # XXX

    if self.mpi_helper.rank == 0:
      self.logger.main_log("Task 2. Analyze the correlation coefficient matrix")

    self.logger.log_step_time("COSYM", True)
    self.logger.log("Memory usage: %d MB"%get_memory_usage())

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(reindex_to_reference)
