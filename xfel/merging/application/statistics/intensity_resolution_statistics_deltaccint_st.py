from __future__ import absolute_import, division, print_function
from six.moves import range
from xfel.merging.application.worker import worker
from dials.array_family import flex
from xfel.merging.application.reflection_table_utils import reflection_table_utils
from xfel.merging.application.statistics.intensity_resolution_statistics import intensity_resolution_statistics
import numpy as np
from dxtbx import flumpy
import io

class intensity_resolution_statistics_deltaccint_st(intensity_resolution_statistics):
  '''Calculates hkl intensity statistics for resolution bins'''

  def __repr__(self):
    return 'Intensity resolution statistics deltaccint st'

  def run(self, experiments, reflections):
    self.logger.log_step_time("INTENSITY_STATISTICS_DELTACCINT_ST")

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI

    assert experiments == None, "Must be run after the group worker"
    # need at least 3 reflections to keep multiplicty > 2 when removing an image
    filtered = flex.reflection_table()
    min_mult = max(3, self.params.merging.minimum_multiplicity)
    for refls in reflection_table_utils.get_next_hkl_reflection_table(reflections=reflections):
      if len(set(refls['id'])) >= min_mult:
        filtered.extend(refls)

    # Get list of all experiment identifiers
    expt_ids = comm.gather(filtered.experiment_identifiers().values(), 0)
    if self.mpi_helper.rank == 0:
      all_expt_ids = flex.std_string()
      for expt_ids_ in expt_ids:
        all_expt_ids.extend(expt_ids_)
      expt_ids = list(set(all_expt_ids))
    else:
      expt_ids = None
    all_expt_ids = comm.bcast(expt_ids, 0)
    all_expts_map = dict(zip(all_expt_ids, range(len(all_expt_ids))))

    if self.mpi_helper.rank == 0:
      self.logger.main_log("N experiments to run delta cc st on: %d"%len(all_expt_ids))

    # We need to compute
    # 1) variance of the average intensities -> compute averages of all intensities and then compute variance per bin
    # 2) average variance of the intensities -> compute variances of all intensities and then compute average per bin

    resolution_binner = self.params.statistics.resolution_binner
    hkl_resolution_bins = self.params.statistics.hkl_resolution_bins

    hkl_set = [hkl for hkl in list(set(filtered['miller_index_asymmetric'])) if hkl in hkl_resolution_bins]
    n_hkl = len(hkl_set)
    hkl_map = dict(zip(hkl_set, range(n_hkl)))

    expt_map = filtered.experiment_identifiers()

#for i, refls in enumerate(reflection_table_utils.get_next_hkl_reflection_table(reflections=filtered)):
#  if i == 103: break
#  if len(refls) != len(set(refls['id'])):
#    print(i, len(refls), len(set(refls['id'])), len(refls) == len(set(refls['id'])))



    # N expts (all ranks) x N hkl (this rank)
    sums     = np.zeros((len(all_expt_ids), n_hkl), float)
    n        = np.zeros((len(all_expt_ids), n_hkl), int)
    merged   = np.zeros((len(all_expt_ids), n_hkl), float)
    variance = np.zeros((len(all_expt_ids), n_hkl), float)

    for refls in reflection_table_utils.get_next_hkl_reflection_table(reflections=filtered):
      hkl = refls['miller_index_asymmetric'][0]
      if hkl not in hkl_map: continue
      hkl_idx = hkl_map[hkl]
      self.logger.log("%s %d"%(str(hkl), len(refls)))
      intensity = flumpy.to_numpy(refls['intensity.sum.value'])

      # set the sum and n for all expts for this hkl
      sums[:,hkl_idx] = np.sum(intensity)
      n   [:,hkl_idx] = len(refls)

      # For each reflection, subtract it off the experiment it came from, leaving it for the rest of the experiments
      for i in range(len(refls)):
        expt_idx = all_expts_map[expt_map[refls['id'][i]]]
        sums[expt_idx,hkl_idx] -= intensity[i]
        n   [expt_idx,hkl_idx] -= 1
      merged[:,hkl_idx] = sums[:,hkl_idx]/n[:,hkl_idx]

      # compute variance for each hkl (less the intensity from each experiment)
      diff_sq_ = np.zeros(len(all_expt_ids))
      for expt_idx in range(len(all_expt_ids)):
        mean_intensity = merged[expt_idx,hkl_idx]
        diff_sq_[expt_idx] = np.sum((intensity-mean_intensity)**2)

      for i in range(len(refls)):
        expt_idx = all_expts_map[expt_map[refls['id'][i]]]
        mean_intensity = merged[expt_idx,hkl_idx]
        diff_sq_[expt_idx] -= (intensity[i] - mean_intensity)**2
        
      variance[:,hkl_idx] = diff_sq_ / (n[:,hkl_idx]-1)

    # N expts (all ranks) x N bins
    n_bins = resolution_binner.n_bins_all()
    all_i_sums   = np.zeros((len(all_expt_ids), n_bins), float) # sum of the averaged intensities
    all_i_n      = np.zeros((len(all_expt_ids), n_bins), int)   # count of the averaged intensities
    all_var_sums = np.zeros((len(all_expt_ids), n_bins), float) # sums of the variances of the intensities

    # Sum up and reduce the bins
    for hkl in hkl_set:
      bin_idx = hkl_resolution_bins[hkl]
      all_i_sums  [:,bin_idx] += merged[:,hkl_map[hkl]]
      all_i_n     [:,bin_idx] += 1
      all_var_sums[:,bin_idx] += variance[:,hkl_map[hkl]]

    total_i_sums   = comm.reduce(all_i_sums,   MPI.SUM, 0)
    total_i_n      = comm.reduce(all_i_n,      MPI.SUM, 0)
    total_var_sums = comm.reduce(all_var_sums, MPI.SUM, 0)

    # Broadcast the average intensities
    if self.mpi_helper.rank == 0:
      total_i_average_ = total_i_sums / total_i_n
    else:
      total_i_average_ = None
    total_i_average = comm.bcast(total_i_average_, 0)

    # Compute the variance of the average intensities
    # First, the numerator, the difference between each hkl and the average for that hkl's bin (ommiting each experiment once)
    diff_sq = np.zeros(merged.shape, float)
    for hkl in hkl_set:
      diff_sq[:,hkl_map[hkl]] = (merged[:,hkl_map[hkl]] - total_i_average[:,hkl_resolution_bins[hkl]]) ** 2

    # N expts (all ranks) x N bins
    diff_sq_sum = np.zeros(all_i_sums.shape, float)
    diff_sq_n   = np.zeros(all_i_sums.shape, int)

    # Complete the numerator for this rank
    for hkl in hkl_set:
      diff_sq_sum[:,hkl_resolution_bins[hkl]] += diff_sq[:,hkl_map[hkl]]
      diff_sq_n  [:,hkl_resolution_bins[hkl]] += 1

    total_diff_sq_sum = comm.reduce(diff_sq_sum, MPI.SUM, 0)
    total_diff_sq_n   = comm.reduce(diff_sq_n,   MPI.SUM, 0)

    # Report
    if self.mpi_helper.rank == 0:
      self.logger.main_log("st")
      sigma_sq_y = total_diff_sq_sum / (total_diff_sq_n-1)
      sigma_sq_e = total_var_sums / total_i_n
      deltaccint_st = (sigma_sq_y - (0.5 * sigma_sq_e)) / (sigma_sq_y + (0.5 * sigma_sq_e))
      for i in range(10):
        self.logger.main_log("Delta CC1/2 sigma-tau, expt %d"%i)
        with io.StringIO() as f:
          resolution_binner.show_data(deltaccint_st[i] * 100, data_fmt="%6.2f", f=f)
          self.logger.main_log(f.getvalue())
      from IPython import embed; embed()

    self.logger.log_step_time("INTENSITY_STATISTICS_DELTACCINT_ST", True)

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(intensity_resolution_statistics_deltaccint)
