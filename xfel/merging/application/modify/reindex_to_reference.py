from __future__ import absolute_import, division, print_function
from dials.array_family import flex
from xfel.merging.application.worker import worker
from xfel.merging.application.scale.experiment_scaler import experiment_scaler
from xfel.merging.application.utils.memory_usage import get_memory_usage
from mmtbx.scaling.twin_analyses import twin_laws
from cctbx import sgtbx, miller
from cctbx.crystal import symmetry
from dxtbx.model.crystal import MosaicCrystalSauter2014

class reindex_to_reference(worker):
  """
  Resolve indexing ambiguity usin twin laws and correlation to a reference set
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(reindex_to_reference, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Resolve indexing ambiguity by comparison to a reference'

  def process_dataframe(self, raw, params):
    #two purposes: 1) write table to disk, 2) make plot
    import pandas as pd
    pd.set_option('display.max_rows', 300)
    pd.set_option('display.max_columns', 8)
    pd.set_option("display.width", None)
    from pandas import DataFrame as df
    data = df(raw)
    print (data)
    symops = list(data.keys())[1:]
    data["max"] = data[symops].max(axis=1)
    data["min"] = data[symops].min(axis=1)
    data['diff'] = data["max"] - data["min"]
    data["reindex_op"] = data[symops].idxmax(axis=1)

    # Now collect dataframes from all ranks and merge
    reports = self.mpi_helper.comm.gather(data,root=0)
    if self.mpi_helper.rank == 0:

      result = pd.concat(reports)

      # Best guess coset decomposition.  Works fine for P6 use case.
      change_of_basis_op_to_niggli_cell = \
      self.params.scaling.i_model.change_of_basis_op_to_niggli_cell()
      minimum_cell_symmetry = self.params.scaling.i_model.crystal_symmetry().change_basis(
      cb_op=change_of_basis_op_to_niggli_cell)
      lattice_group = sgtbx.lattice_symmetry.group(
      reduced_cell=minimum_cell_symmetry.unit_cell())
      Laue_group = self.params.scaling.i_model.space_group().build_derived_laue_group()
      LATG = lattice_group.change_basis(change_of_basis_op_to_niggli_cell.inverse())
      LAUG = Laue_group.build_derived_acentric_group()
      CO = sgtbx.cosets.left_decomposition(LATG,LAUG)
      partitions = CO.partitions

      # Initialize a coset column using best estimate
      working_coset = []
      working_reindex_op = list(result["reindex_op"])
      for iexpt in range(len(working_reindex_op)):
            this_coset = None
            for p_no, partition in enumerate(partitions):
              partition_ops = [sgtbx.change_of_basis_op(ip).as_hkl() for ip in partition]
              if working_reindex_op[iexpt] in partition_ops:
                this_coset = p_no; break
            assert this_coset is not None
            working_coset.append(this_coset)
      result["coset"]=working_coset
      print (result)

      import os
      result.to_pickle(path = os.path.join(params.output.output_dir, params.modify.reindex_to_reference.dataframe))
      #from matplotlib import pyplot as plt
      #plt.hist(data["min"], bins=24, range=[-0.5,1], )
      #plt.hist(data["max"], bins=24, range=[-0.5,1], alpha=0.75)
      #plt.show()
      # later change this plot to produce PDF file

  def run(self, experiments, reflections):

    self.logger.log_step_time("REINDEX")

    # Get list of twinning operators for this space group
    operators = twin_laws(self.params.scaling.i_model).operators
    if not operators:
      self.logger.log("No indexing ambiguity. Skipping this step.")
      return experiments, reflections
    self.logger.log("Resolving indexing ambiguity using operators h,k,l, %s"%", ".join( \
      [op.operator.r().as_hkl() for op in operators]))
    if self.params.modify.reindex_to_reference.dataframe:
      from collections import OrderedDict
      keyval = [("experiment",[]),("h,k,l", [])]
      for op in operators:
        keyval.append((op.operator.r().as_hkl(), []))
      raw = OrderedDict(keyval)
      keys = list(raw.keys())

    operators = [sgtbx.change_of_basis_op(op.operator.r().as_hkl()) for op in operators]

    result = flex.reflection_table()
    scaler = experiment_scaler(self.params, self.mpi_helper, self.logger)
    model_intensities = self.params.scaling.i_model
    target_symm = symmetry(unit_cell = self.params.scaling.unit_cell, space_group_info = self.params.scaling.space_group)

    def get_correlation(cb_op=None):
      """ Helper function to get CC to the reference given an operator """
      # Build a miller array for the experiment reflections
      exp_miller_indices = miller.set(target_symm, exp_reflections['miller_index_asymmetric'], True)
      exp_intensities = miller.array(exp_miller_indices, exp_reflections['intensity.sum.value'],
                                     flex.sqrt(exp_reflections['intensity.sum.variance']))
      if cb_op:
        exp_intensities = exp_intensities.change_basis(cb_op).map_to_asu()

      # Extract an array of HKLs from the model to match the experiment HKLs
      matching_indices = miller.match_multi_indices(miller_indices_unique = model_intensities.indices(), miller_indices = exp_intensities.indices())

      # Least squares
      scaling_result = scaler.fit_experiment_to_reference(model_intensities, exp_intensities, matching_indices)
      return scaling_result.correlation if scaling_result.correlation is not None else -1

    # Test each experiment to see if an operator gives a better CC to the reference, and if it does, apply it
    for expt_id, experiment in enumerate(experiments):
      exp_reflections = reflections.select(reflections['exp_id'] == experiment.identifier)
      all_correlations = []
      best_correlation = get_correlation()
      all_correlations.append(best_correlation)
      best_op = None
      for cb_op in operators:
        test_correlation = get_correlation(cb_op)
        all_correlations.append(test_correlation)
        if test_correlation > best_correlation:
          best_correlation = test_correlation
          best_op = cb_op
      if best_op:
        exp_miller_indices = miller.set(target_symm, exp_reflections['miller_index'], True).change_basis(best_op)
        exp_reflections['miller_index_asymmetric'] = exp_miller_indices.map_to_asu().indices()
        exp_reflections['miller_index'] = exp_miller_indices.indices()
        experiment.crystal = MosaicCrystalSauter2014(experiment.crystal.change_basis(best_op)) # need to use wrapper because of cctbx/dxtbx#5
      result.extend(exp_reflections)

      self.logger.log("Expt %d, reindexing op correlations: %s"%(expt_id, ", ".join(["%6.3f"%c for c in all_correlations])))
      if self.params.modify.reindex_to_reference.dataframe:
        raw["experiment"].append(experiment.identifier)
        for ix in range(len(all_correlations)):
          raw[keys[ix+1]].append(all_correlations[ix])

    if self.params.modify.reindex_to_reference.dataframe:
      self.process_dataframe(raw, self.params)

    self.logger.log_step_time("REINDEX", True)
    self.logger.log("Memory usage: %d MB"%get_memory_usage())

    from xfel.merging.application.utils.data_counter import data_counter
    data_counter(self.params).count(experiments, reflections)
    return experiments, result

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(reindex_to_reference)
