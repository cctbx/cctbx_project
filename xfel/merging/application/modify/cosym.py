from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from xfel.merging.application.utils.memory_usage import get_memory_usage
from cctbx.array_family import flex
from cctbx import sgtbx
from cctbx.sgtbx import change_of_basis_op, rt_mx
from dxtbx.model.crystal import MosaicCrystalSauter2014
from dxtbx.model.experiment_list import ExperimentList
#from libtbx.development.timers import Profiler, Timer


class cosym(worker):
  """
  Resolve indexing ambiguity using dials.cosym
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(cosym, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Resolve indexing ambiguity using dials.cosym'


  @staticmethod
  def experiment_id_detail(expts, refls, exp_reflections): # function arguments are modified by the function
        simple_experiment_id = len(expts) - 1
        #experiment.identifier = "%d"%simple_experiment_id
        expts[-1].identifier = "%d"%simple_experiment_id
        # experiment identifier must be a string according to *.h file
        # the identifier is changed on the _for_cosym Experiment list, not the master experiments for through analysis

        exp_reflections['id'] = flex.int(len(exp_reflections), simple_experiment_id)
        # register the integer id as a new column in the per-experiment reflection table

        exp_reflections.experiment_identifiers()[simple_experiment_id] = expts[-1].identifier
        #apparently the reflection table holds a map from integer id (reflection table) to string id (experiment)

        refls.append(exp_reflections)

  @staticmethod
  def task_a(params):
      # add an anchor
      sampling_experiments_for_cosym = ExperimentList()
      sampling_reflections_for_cosym = []
      if params.modify.cosym.anchor:
        from xfel.merging.application.model.crystal_model import crystal_model
        #P = Timer("construct the anchor reference model")
        XM = crystal_model(params = params, purpose="cosym")
        model_intensities = XM.run([],[])
        #del P
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
        exp_reflections['flags'] = flex.size_t(model_intensities.size(), flex.reflection_table.flags.integrated_sum)

        # prepare individual reflection tables for each experiment
        cosym.experiment_id_detail(sampling_experiments_for_cosym, sampling_reflections_for_cosym, exp_reflections)
      return sampling_experiments_for_cosym, sampling_reflections_for_cosym

  @staticmethod
  def task_c(params, mpi_helper, logger, tokens,
      sampling_experiments_for_cosym, sampling_reflections_for_cosym, uuid_starting=[], communicator_size=1, do_plot=False):
      # Purpose: assemble a composite tranch from src inputs, and instantiate the COSYM class
      uuid_cache = uuid_starting

      for tranch_experiments, tranch_reflections in tokens:
          for expt_id, experiment in enumerate(tranch_experiments):
            sampling_experiments_for_cosym.append(experiment)
            uuid_cache.append(experiment.identifier)

            exp_reflections = tranch_reflections.select(tranch_reflections['id'] == expt_id)
            # prepare individual reflection tables for each experiment

            cosym.experiment_id_detail(sampling_experiments_for_cosym, sampling_reflections_for_cosym, exp_reflections)

      from dials.command_line import cosym as cosym_module
      cosym_module.logger = logger

      i_plot = mpi_helper.rank
      from xfel.merging.application.modify.aux_cosym import dials_cl_cosym_subclass as dials_cl_cosym_wrapper
      COSYM = dials_cl_cosym_wrapper(
                sampling_experiments_for_cosym, sampling_reflections_for_cosym,
                uuid_cache, params=params.modify.cosym,
                output_dir=params.output.output_dir, do_plot=do_plot,
                i_plot=i_plot)
      return COSYM


  def run(self, input_experiments, input_reflections):
    from collections import OrderedDict
    if self.mpi_helper.rank == 0:
      print("Starting cosym worker")
      #Overall = Profiler("Cosym total time")

    #  Evenly distribute all experiments from mpi_helper ranks
    reports = self.mpi_helper.comm.gather((len(input_experiments)),root=0) # report from all ranks on experiment count
    if self.mpi_helper.rank == 0:
      from xfel.merging.application.modify.token_passing_left_right import construct_src_to_dst_plan
      plan = construct_src_to_dst_plan(flex.int(reports), self.params.modify.cosym.tranch_size, self.mpi_helper.comm)
    else:
      plan = 0
    plan = self.mpi_helper.comm.bcast(plan, root = 0)
    dst_offset = 1 if self.mpi_helper.size>1 else 0 # decision whether to reserve rank 0 for parallel anchor determination
                                                    # FIXME XXX probably need to look at plan size to decide dst_offset or not
    from xfel.merging.application.modify.token_passing_left_right import apply_all_to_all
    tokens = apply_all_to_all(plan=plan, dst_offset=dst_offset,
                   value=(input_experiments, input_reflections), comm = self.mpi_helper.comm)

    if self.params.modify.cosym.anchor:
      if self.mpi_helper.rank == 0:
        MIN_ANCHOR = 20
        from xfel.merging.application.modify.token_passing_left_right import construct_anchor_src_to_dst_plan
        anchor_plan = construct_anchor_src_to_dst_plan(MIN_ANCHOR, flex.int(reports), self.params.modify.cosym.tranch_size, self.mpi_helper.comm)
      else:
        anchor_plan = 0
      anchor_plan = self.mpi_helper.comm.bcast(anchor_plan, root = 0)
    self.logger.log_step_time("COSYM")

    if self.params.modify.cosym.plot.interactive:
      self.params.modify.cosym.plot.filename = None

    has_tokens = len(tokens) > 0
    all_has_tokens = self.mpi_helper.comm.allgather(has_tokens)
    ranks_with_tokens = [i for (i, val) in enumerate(all_has_tokens) if val]
    ranks_to_plot = ranks_with_tokens[:self.params.modify.cosym.plot.n_max]
    do_plot = (self.params.modify.cosym.plot.do_plot
        and self.mpi_helper.rank in ranks_to_plot)

    if len(tokens) > 0: # Only select ranks that have been assigned tranch data, for mutual coset determination
      # because cosym has a problem with hashed identifiers, use simple experiment identifiers
      sampling_experiments_for_cosym = ExperimentList()
      sampling_reflections_for_cosym = [] # is a list of flex.reflection_table
      COSYM = self.task_c(self.params, self.mpi_helper, self.logger, tokens,
        sampling_experiments_for_cosym, sampling_reflections_for_cosym,
        communicator_size=self.mpi_helper.size, do_plot=do_plot)
      self.uuid_cache = COSYM.uuid_cache # reformed uuid list after n_refls filter


      rank_N_refl=flex.double([r.size() for r in COSYM.reflections])
      message = """Task 1. Prepare the data for cosym
    change_of_basis_ops_to_minimum_cell
    eliminate_sys_absent
    transform models into Miller arrays, putting data in primitive triclinic reduced cell
    There are %d experiments with %d reflections, averaging %.1f reflections/experiment"""%(
      len(COSYM.experiments), flex.sum(rank_N_refl), flex.mean(rank_N_refl))
      self.logger.log(message)
      if self.mpi_helper.rank == 1: print(message) #; P = Timer("COSYM.run")
      COSYM.run()
      #if self.mpi_helper.rank == 1: del P

      keyval = [("experiment", []), ("reindex_op", []), ("coset", [])]
      raw = OrderedDict(keyval)

      if self.mpi_helper.rank == 0: print("Rank",self.mpi_helper.rank,"experiments:",len(sampling_experiments_for_cosym))

      for sidx in range(len(self.uuid_cache)):
        raw["experiment"].append(self.uuid_cache[sidx])

        sidx_plus = sidx

        try:
          minimum_to_input = COSYM.cb_op_to_minimum[sidx_plus].inverse()
        except Exception as e:
          print ("raising",e,sidx_plus, len(COSYM.cb_op_to_minimum))
          raise e

        reindex_op = minimum_to_input * \
                     sgtbx.change_of_basis_op(COSYM.cosym_analysis.reindexing_ops[sidx_plus]) * \
                     COSYM.cb_op_to_minimum[sidx_plus]

        # Keep this block even though not currently used; need for future assertions:
        LG = COSYM.cosym_analysis.target._lattice_group
        LGINP = LG.change_basis(COSYM.cosym_analysis.cb_op_inp_min.inverse()).change_basis(minimum_to_input)
        SG = COSYM.cosym_analysis.input_space_group
        SGINP = SG.change_basis(COSYM.cosym_analysis.cb_op_inp_min.inverse()).change_basis(minimum_to_input)
        CO = sgtbx.cosets.left_decomposition(LGINP, SGINP)
        partitions = CO.partitions
        this_reindex_op = reindex_op.as_hkl()
        this_coset = None
        for p_no, partition in enumerate(partitions):
          # Note: For centered cells it appears the partitions may have a translational part.
          # We eliminate it by round-tripping the rt_mx to its rotational part and back.
          partition_ops = [change_of_basis_op(rt_mx(ip.r())).as_hkl() for ip in partition]
          if this_reindex_op in partition_ops:
            this_coset = p_no; break
        assert this_coset is not None
        raw["coset"].append(this_coset)
        raw["reindex_op"].append(this_reindex_op)

      keys = list(raw.keys())
      from pandas import DataFrame as df
      data = df(raw)
      # major assumption is that all the coset decompositions "CO" are the same.  NOT sure if a test is needed.
      reports = self.mpi_helper.comm.gather((data, CO),root=0)
    else:
      reports = self.mpi_helper.comm.gather(None,root=0)
    if self.mpi_helper.rank == 0:
      # report back to rank==0 and reconcile all coset assignments
      while None in reports:
        reports.pop(reports.index(None))
      # global CO
      global_coset_decomposition = reports[0][1] # again, assuming here they are all the same XXX
    else:
      global_coset_decomposition = 0
    global_coset_decomposition = self.mpi_helper.comm.bcast(global_coset_decomposition, root=0)
    partitions = global_coset_decomposition.partitions
    self.mpi_helper.comm.barrier()
    # end of distributed embedding

    if self.params.modify.cosym.anchor:
        anchor_tokens = apply_all_to_all(plan=anchor_plan, dst_offset=0,
        value=(input_experiments, input_reflections), comm = self.mpi_helper.comm)

    if self.mpi_helper.rank == 0:
        from xfel.merging.application.modify.df_cosym import reconcile_cosym_reports
        REC = reconcile_cosym_reports(reports)
        results = REC.composite_tranch_merge(voting_method="consensus")

        # at this point we have the opportunity to reconcile the results with an anchor
        # recycle the data structures for anchor determination
        if self.params.modify.cosym.anchor:
          sampling_experiments_for_cosym, sampling_reflections_for_cosym = self.task_a(self.params)
          ANCHOR = self.task_c(self.params, self.mpi_helper, self.logger, anchor_tokens,
            sampling_experiments_for_cosym, sampling_reflections_for_cosym,
            uuid_starting=["anchor structure"], communicator_size=1) # only run on the rank==0 tranch.
          self.uuid_cache = ANCHOR.uuid_cache # reformed uuid list after n_refls filter
          #P = Timer("ANCHOR.run")
          ANCHOR.run() # Future redesign XXX FIXME do this in rank 0 in parallel with distributed composite tranches
          #del P

          keyval = [("experiment", []), ("coset", [])]
          raw = OrderedDict(keyval)
          print("Anchor","experiments:",len(sampling_experiments_for_cosym))

          anchor_op = ANCHOR.cb_op_to_minimum[0].inverse() * \
                     sgtbx.change_of_basis_op(ANCHOR.cosym_analysis.reindexing_ops[0]) * \
                     ANCHOR.cb_op_to_minimum[0]
          anchor_coset = None
          for p_no, partition in enumerate(partitions):
              partition_ops = [change_of_basis_op(rt_mx(ip.r())).as_hkl() for ip in partition]
              if anchor_op.as_hkl() in partition_ops:
                anchor_coset = p_no; break
          assert anchor_coset is not None
          print("The consensus for the anchor is",anchor_op.as_hkl()," anchor coset", anchor_coset)

          raw["experiment"].append("anchor structure"); raw["coset"].append(anchor_coset)
          for sidx in range(1,len(self.uuid_cache)):
            raw["experiment"].append(self.uuid_cache[sidx])

            sidx_plus = sidx

            minimum_to_input = ANCHOR.cb_op_to_minimum[sidx_plus].inverse()
            reindex_op = minimum_to_input * \
                     sgtbx.change_of_basis_op(ANCHOR.cosym_analysis.reindexing_ops[sidx_plus]) * \
                     ANCHOR.cb_op_to_minimum[sidx_plus]
            this_reindex_op = reindex_op.as_hkl()
            this_coset = None
            for p_no, partition in enumerate(partitions):
              partition_ops = [change_of_basis_op(rt_mx(ip.r())).as_hkl() for ip in partition]
              if this_reindex_op in partition_ops:
                this_coset = p_no; break
            assert this_coset is not None
            raw["coset"].append(this_coset)

          from pandas import DataFrame as df
          anchor_data = df(raw)
          REC.reconcile_with_anchor(results, anchor_data, anchor_op)
          # no need for return value; results dataframe is modified in place

        if self.params.modify.cosym.dataframe:
          import os
          results.to_pickle(path = os.path.join(self.params.output.output_dir,self.params.modify.cosym.dataframe))
        transmitted = results
    else:
        transmitted = 0
    self.mpi_helper.comm.barrier()
    transmitted = self.mpi_helper.comm.bcast(transmitted, root = 0)
    # "transmitted" holds the global coset assignments

    #subselect expt and refl on the successful coset assignments
    # output:  experiments-->result_experiments_for_cosym; reflections-->reflections (modified in place)
    result_experiments_for_cosym = ExperimentList()
    good_refls = flex.bool(len(input_reflections), False)
    good_expt_id = list(transmitted["experiment"])
    good_coset = list(transmitted["coset"]) # would like to understand how to use pandas rather than Python list
    for iexpt in range(len(input_experiments)):
        iexpt_id = input_experiments[iexpt].identifier
        keepit = iexpt_id in good_expt_id
        if keepit:
          this_coset = good_coset[ good_expt_id.index(iexpt_id) ]
          this_cb_op = change_of_basis_op(global_coset_decomposition.partitions[this_coset][0])
          accepted_expt = input_experiments[iexpt]
          if this_coset > 0:
            accepted_expt.crystal = MosaicCrystalSauter2014( accepted_expt.crystal.change_basis(this_cb_op) )
                                    # need to use wrapper because of cctbx/dxtbx#5
          result_experiments_for_cosym.append(accepted_expt)
          good_refls |= input_reflections["id"] == iexpt
    selected_reflections = input_reflections.select(good_refls)
    selected_reflections.reset_ids()
    self.mpi_helper.comm.barrier()

    # still have to reindex the reflection table, but try to do it efficiently
    from xfel.merging.application.modify.reindex_cosym import reindex_refl_by_coset
    if (len(result_experiments_for_cosym) > 0):
      reindex_refl_by_coset(refl = selected_reflections,
                          data = transmitted,
                          symms=[E.crystal.get_crystal_symmetry() for E in result_experiments_for_cosym],
                          uuids=[E.identifier for E in result_experiments_for_cosym],
                          co=global_coset_decomposition,
                          anomalous_flag = self.params.merging.merge_anomalous==False,
                          verbose=False)
    # this should have re-indexed the refls in place, no need for return value

    self.mpi_helper.comm.barrier()
    # Note:  this handles the simple case of lattice ambiguity (P63 in P/mmm lattice group)
    # in this use case we assume all inputs and outputs are in P63.
    # more complex use cases would have to reset the space group in the crystal, and recalculate
    # the ASU "miller_indicies" in the reflections table.

    self.logger.log_step_time("COSYM", True)
    self.logger.log("Memory usage: %d MB"%get_memory_usage())

    from xfel.merging.application.utils.data_counter import data_counter
    data_counter(self.params).count(result_experiments_for_cosym, selected_reflections)
    return result_experiments_for_cosym, selected_reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(reindex_to_reference)
