from __future__ import absolute_import, division, print_function
from six.moves import range
from dxtbx import flumpy
from dials.array_family import flex
import itertools
import numpy as np
from scitbx.math import five_number_summary
from xfel.merging.application.reflection_table_utils import reflection_table_utils
from xfel.merging.application.worker import worker


class deltaccint(worker):
  '''Calculates ΔCC½ using the σ-τ method from Assmann 2016'''

  def __repr__(self):
    return 'ΔCC½ statistics (σ-τ method)'

  def run(self, experiments, reflections):
    self.logger.log_step_time("STATISTICS_DELTA_CCINT")

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI

    assert experiments is None, "Must be run after the group worker"
    # need at least 3 reflections to keep multiplicty > 2 when removing an image
    min_mult = max(3, self.params.merging.minimum_multiplicity)
    hkl_resolution_bins = self.params.statistics.hkl_resolution_bins

    if self.mpi_helper.rank == 0:
      self.logger.main_log("Beginning ΔCC½ analysis (σ-τ method from Assmann 2016)")
      self.logger.main_log("Removing reflections with less than %d measurements"%min_mult)
      self.logger.main_log("")

    # Iterate: after each removal pass, recompute ΔCC½ on the reduced set so
    # that a second (or third) poison lattice hidden behind a larger outlier can
    # emerge and be caught.  Stop when a pass removes nothing or iqr_ratio is
    # not configured.
    for iteration in range(10):
      if self.mpi_helper.rank == 0:
        self.logger.main_log("ΔCC½ pass %d"%(iteration+1))
        self.logger.main_log("")

      # A rank may hold zero reflections after the group worker; guard the hkl
      # iterator, which raises on an empty table, so every rank still reaches the
      # collective calls below (otherwise the empty rank aborts and the others hang).
      filtered = flex.reflection_table()
      if len(reflections) > 0:
        for refls in reflection_table_utils.get_next_hkl_reflection_table(reflections=reflections):
          if len(set(refls['id'])) >= min_mult:
            filtered.extend(refls)

      expt_map = filtered.experiment_identifiers()

      all_expt_ids = sorted(set(itertools.chain.from_iterable(comm.allgather(expt_map.values()))))
      all_expts_map = {v: k for k, v in enumerate(all_expt_ids)}

      if self.mpi_helper.rank == 0:
        self.logger.main_log("N experiments after filtering: %d"%len(all_expt_ids))
        self.logger.main_log("")

      # We need to compute
      # 1) variance of the average intensities -> compute averages of all intensities and then compute variance per bin
      # 2) average variance of the intensities -> compute variances of all intensities and then compute average per bin
      # In reality for 2), we compute the average of the standard error of the mean squared instead of the variance (semsq)

      # filtered has no columns when this rank kept nothing, so guard the lookup
      if len(filtered) > 0:
        hkl_set = [hkl for hkl in set(filtered['miller_index_asymmetric']) if hkl in hkl_resolution_bins]
      else:
        hkl_set = []
      n_hkl = len(hkl_set)
      hkl_map = {v: k for k, v in enumerate(hkl_set)}

      # N expts (all ranks) x N hkl (this rank)
      sums   = np.zeros((len(all_expt_ids), n_hkl), float)
      n      = np.zeros((len(all_expt_ids), n_hkl), int)
      merged = np.zeros((len(all_expt_ids), n_hkl), float)
      semsq  = np.zeros((len(all_expt_ids), n_hkl), float)

      # empty on ranks that kept nothing; skip the iterator that raises on empty
      filtered_hkl_tables = reflection_table_utils.get_next_hkl_reflection_table(reflections=filtered) if len(filtered) > 0 else []
      for refls in filtered_hkl_tables:
        hkl = refls['miller_index_asymmetric'][0]
        if hkl not in hkl_map: continue
        hkl_idx = hkl_map[hkl]
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

        # compute semsq for each hkl (less the intensity from each experiment).
        # diff_sq_ holds, per experiment, the sum of squared deviations of the
        # measurements remaining once that experiment is removed. Track which
        # experiments contributed with a separate boolean mask: a sum of squares
        # can legitimately be zero (or dip slightly negative from floating-point
        # cancellation), so its sign must not double as an "unset" sentinel.
        diff_sq_ = np.zeros(len(all_expt_ids), dtype=float)
        contributed = np.zeros(len(all_expt_ids), dtype=bool)

        for i in range(len(refls)):
          expt_idx = all_expts_map[expt_map[refls['id'][i]]]
          mean_intensity_modified = merged[expt_idx,hkl_idx]
          # First measurement from this experiment: start from the full sum of
          # squared deviations (about the leave-one-out mean) and remove this term.
          # Further measurements (a single image measuring this hkl more than once)
          # just remove their term, excluding all of the experiment's measurements.
          if not contributed[expt_idx]:
            diff_sq_modified = np.sum((intensity-mean_intensity_modified)**2)
            diff_sq_[expt_idx] = diff_sq_modified - (intensity[i] - mean_intensity_modified)**2
            contributed[expt_idx] = True
          else:
            diff_sq_[expt_idx] -= (intensity[i] - mean_intensity_modified)**2

        # Experiments that never measured this hkl keep the full sum of squared
        # deviations about the unmodified mean (their leave-one-out set is unchanged).
        mean_intensity_unmodified = np.mean(intensity)
        diff_sq_unmodified = np.sum((intensity-mean_intensity_unmodified)**2)
        diff_sq_[~contributed] = diff_sq_unmodified
        # clamp tiny negatives from floating-point cancellation (a SS is >= 0)
        np.clip(diff_sq_, 0, None, out=diff_sq_)
        semsq[:,hkl_idx] = diff_sq_ / (n[:,hkl_idx]-1) / n[:,hkl_idx]

      # N expts (all ranks)
      all_i_sums     = np.zeros(len(all_expt_ids), float) # sum of the averaged intensities
      all_i_n        = np.zeros(len(all_expt_ids), int)   # count of the averaged intensities
      all_semsq_sums = np.zeros(len(all_expt_ids), float) # sums of the semsq of the intensities

      # Sum up and reduce the bins
      for hkl in hkl_set:
        all_i_sums     += merged[:,hkl_map[hkl]]
        all_i_n        += 1
        all_semsq_sums += semsq[:,hkl_map[hkl]]

      # Broadcast the semsq and average intensities
      total_i_sums     = comm.allreduce(all_i_sums, op=MPI.SUM)
      total_i_n        = comm.allreduce(all_i_n,    op=MPI.SUM)
      total_semsq_sums = comm.reduce(all_semsq_sums,  op=MPI.SUM)
      total_i_average  = total_i_sums / total_i_n

      # Compute the variance of the average intensities
      # First, the numerator, the difference between each hkl and the average for that hkl's bin (ommiting each experiment once)
      diff_sq = np.zeros(merged.shape, float)
      for hkl in hkl_set:
        diff_sq[:,hkl_map[hkl]] = (merged[:,hkl_map[hkl]] - total_i_average) ** 2

      # N expts (all ranks) x N bins
      diff_sq_sum = np.zeros(all_i_sums.shape, float)

      # Complete the numerator for this rank
      for hkl in hkl_set:
        diff_sq_sum += diff_sq[:,hkl_map[hkl]]

      total_diff_sq_sum = comm.reduce(diff_sq_sum, MPI.SUM, 0)

      # Report
      if self.mpi_helper.rank == 0:
        sigma_sq_y = total_diff_sq_sum / (total_i_n-1) # variance of the average intensities
        sigma_sq_e = 2 * total_semsq_sums / total_i_n    # average semsq of the intensities
        deltaccint_st = (sigma_sq_y - (0.5 * sigma_sq_e)) / (sigma_sq_y + (0.5 * sigma_sq_e))

        data = flex.double(deltaccint_st) * 100
        sorted_data = data.select(flex.sort_permutation(data))

        mini, q1, med, q3, maxi = five_number_summary(data)
        self.logger.main_log("Five number summary")
        self.logger.main_log("% 8.4f%% min"%mini)
        self.logger.main_log("% 8.4f%% q1"%q1)
        self.logger.main_log("% 8.4f%% median"%med)
        self.logger.main_log("% 8.4f%% q3"%q3)
        self.logger.main_log("% 8.4f%% max"%maxi)
        self.logger.main_log("")

        if self.params.statistics.deltaccint.verbose:
          self.logger.main_log("Showing ΔCC½ for all lattices")
          for e, identifier in enumerate(all_expt_ids):
            self.logger.main_log("%s %f"%(identifier, data[e]))

        n_worst = min(len(data), 30)
        worst = sorted_data[-n_worst:]
        iqr = q3-q1

        self.logger.main_log("Showing ΔCC½ of the worst %d lattices and IQR ratio needed to remove them"%n_worst)
        self.logger.main_log(" ΔCC½ (%) IQR ratio Lattices removed")
        for i in range(len(worst)):
          iqr_ratio_needed = (worst[i]-med)/iqr if iqr > 0 else float('nan')
          self.logger.main_log("% 8.4f % 10.1f %d"%(worst[i], iqr_ratio_needed, n_worst-i))
        self.logger.main_log("")

        if self.params.statistics.deltaccint.iqr_ratio:
          iqr_ratio = self.params.statistics.deltaccint.iqr_ratio
          sd = float(np.std(deltaccint_st, ddof=1)) * 100 if len(deltaccint_st) >= 2 else 0.0
          if iqr >= sd * 0.1:
            # Normal case: IQR is a meaningful fraction of the spread (for a
            # Gaussian, IQR/SD ≈ 1.35, well above the 0.1 threshold).
            cut = q3 + iqr * iqr_ratio
          else:
            # Degenerate distribution: IQR is zero or negligibly small relative
            # to the overall spread. This happens when one or more outliers with
            # enormous intensities compress the bulk into a near-identical cluster,
            # making IQR ≈ 0 while SD is large. Using IQR would place the cutoff
            # at ~q3, sweeping up the bulk. Fall back to SD, which is inflated by
            # the very outliers we want to remove and places the cut well above the
            # bulk spike.
            cut = q3 + iqr_ratio * sd
            self.logger.main_log("IQR of ΔCC½ is negligible relative to spread (degenerate distribution); using standard-deviation spread instead")
          sel = data > cut
          worst_expts_ = flex.std_string(all_expt_ids).select(sel)
          self.logger.main_log("Removing %d experiments out of %d using ΔCC½ cutoff %.4f%%"%(len(worst_expts_), len(all_expt_ids), cut))
        else:
          worst_expts_ = None

      # Broadcast the worst experiments to cut
      else:
        worst_expts_ = None

      if self.params.statistics.deltaccint.iqr_ratio:
        worst_expts = comm.bcast(worst_expts_, 0)
        self.logger.log("Starting number of reflections: %d"%len(reflections))
        self.logger.log("Reflections after filtering by minimum multiplicity of %d: %d"%(min_mult, len(filtered)))
        if len(worst_expts) == 0:
          break  # converged: nothing removed this pass
        reflections.remove_on_experiment_identifiers(worst_expts)
        reflections.reset_ids()
        self.logger.log("Reflections after filtering by ΔCC½: %d"%len(reflections))
      else:
        break  # iqr_ratio not configured; report only, no removal

    self.logger.log_step_time("STATISTICS_DELTA_CCINT", True)

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(deltaccint)
