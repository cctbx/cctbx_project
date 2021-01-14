from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from xfel.merging.application.utils.memory_usage import get_memory_usage
from cctbx.array_family import flex
from cctbx import sgtbx
from cctbx.sgtbx import change_of_basis_op

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
    if False:
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
      self.uuid_cache = []
      if self.mpi_helper.size == 1: # simple case, one rank
       for experiment in all_sampling_experiments:
        sampling_experiments_for_cosym.append(experiment)
        self.uuid_cache.append(experiment.identifier)

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
      else: # complex case, overlap tranches for mutual coset determination
        self.mpi_helper.MPI.COMM_WORLD.barrier()
        from xfel.merging.application.modify.token_passing_left_right import token_passing_left_right
        values = token_passing_left_right((experiments,reflections))
        for tranch_experiments, tranch_reflections in values:
          for experiment in tranch_experiments:
            sampling_experiments_for_cosym.append(experiment)
            self.uuid_cache.append(experiment.identifier)

            exp_reflections = tranch_reflections.select(tranch_reflections['exp_id'] == experiment.identifier)
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

      from xfel.merging.application.modify.aux_cosym import dials_cl_cosym_subclass as dials_cl_cosym_wrapper
      COSYM = dials_cl_cosym_wrapper(
                sampling_experiments_for_cosym, sampling_reflections_for_cosym,
                params=self.params.modify.cosym)
      return COSYM
    COSYM = task_1()

    # runtime code specialization, replace Gildea algorithm with Paley
    from dials.algorithms.symmetry.cosym.target import Target
    from xfel.merging.test import Reproducer
    Target._compute_rij_wij = Reproducer.compute_rij_wij_cplusplus # fastest implementation so far
    #Target._compute_rij_wij = Reproducer.compute_rij_wij_Gildea # not so fast

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

    if self.params.modify.cosym.dataframe:
      from collections import OrderedDict
      #assert len(sampling_experiments_for_cosym) + 1 anchor if present == len(COSYM._experiments)
      keyval = [("experiment", []), ("reindex_op", []), ("coset", [])]
      raw = OrderedDict(keyval)
      print("Rank",self.mpi_helper.rank,"experiments:",len(sampling_experiments_for_cosym))

      for sidx in range(len(self.uuid_cache)):
        raw["experiment"].append(self.uuid_cache[sidx])
      #for sidx in range(len(experiments)):
      #  raw["experiment"].append(experiments[sidx].identifier)

        sidx_plus = sidx
        if self.mpi_helper.rank == 0 and self.params.modify.cosym.anchor:# add 1 if COSYM had an anchor
          sidx_plus += 1
          if sidx == 0:
            reindex_op = COSYM.cb_op_to_minimum[0].inverse() * \
                     sgtbx.change_of_basis_op(COSYM.cosym_analysis.reindexing_ops[0][0]) * \
                     COSYM.cb_op_to_minimum[0]
            print("The consensus for the anchor is",reindex_op.as_hkl())

        minimum_to_input = COSYM.cb_op_to_minimum[sidx_plus].inverse()
        reindex_op = minimum_to_input * \
                     sgtbx.change_of_basis_op(COSYM.cosym_analysis.reindexing_ops[sidx_plus][0]) * \
                     COSYM.cb_op_to_minimum[sidx_plus]

        # Keep this block even though not currently used; need for future assertions:
        LG = COSYM.cosym_analysis.target._lattice_group
        LGINP = LG.change_basis(COSYM.cosym_analysis.cb_op_inp_min.inverse()).change_basis(minimum_to_input)
        SG = COSYM.cosym_analysis.space_groups[sidx_plus]
        SGINP = SG.change_basis(COSYM.cosym_analysis.cb_op_inp_min.inverse()).change_basis(minimum_to_input)
        CO = sgtbx.cosets.left_decomposition(LGINP, SGINP)
        # if sidx_plus==10: CO.show() # for debugging
        partitions = CO.partitions
        this_reindex_op = reindex_op.as_hkl()
        this_coset = None
        for p_no, partition in enumerate(partitions):
          partition_ops = [change_of_basis_op(ip).as_hkl() for ip in partition]
          if this_reindex_op in partition_ops:
            this_coset = p_no; break
        assert this_coset is not None
        raw["coset"].append(this_coset)
        raw["reindex_op"].append(this_reindex_op)

      keys = list(raw.keys())
      from pandas import DataFrame as df
      data = df(raw)

      # report back to rank==0 and reconcile all coset assignments
      reports = self.mpi_helper.MPI.COMM_WORLD.gather((data, CO),root=0)
      if self.mpi_helper.rank == 0:
        from xfel.merging.application.modify.df_cosym import reconcile_cosym_reports
        REC = reconcile_cosym_reports(reports)
        results = REC.simple_merge(voting_method="consensus")
        results.to_pickle(path = "cosym_myfile")



    if self.mpi_helper.rank == 0:
      self.logger.main_log("Task 2. Analyze the correlation coefficient matrix")



    self.logger.log_step_time("COSYM", True)
    self.logger.log("Memory usage: %d MB"%get_memory_usage())

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(reindex_to_reference)
