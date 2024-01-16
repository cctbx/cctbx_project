from __future__ import absolute_import, division, print_function
from dials.array_family import flex
from xfel.merging.application.worker import worker
from xfel.merging.application.reflection_table_utils import reflection_table_utils
import math
import numpy as np
import os
import sys
from scipy.special import gamma
from scipy.special import polygamma
import scipy.stats

number_of_intensity_bins = 100


class error_modifier_ev11(worker):
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(error_modifier_ev11, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)
    self.algorithm = self.params.merging.error.ev11.algorithm
    self.likelihood = self.params.merging.error.ev11.likelihood
    self.limit_differences = self.params.merging.error.ev11.limit_differences
    self.n_coefs = self.params.merging.error.ev11.n_degrees + 1
    self.n_max_differences = self.params.merging.error.ev11.n_max_differences
    self.overall_scaling_method = self.params.merging.error.ev11.overall_scaling_method
    self.tuning_param = self.params.merging.error.ev11.tuning_param
    self.tuning_param_opt = self.params.merging.error.ev11.tuning_param_opt
    if self.params.merging.error.ev11.cc_after_pr:
      self.cc_key = 'correlation_after_post'
    else:
      self.cc_key = 'correlation'

  def __repr__(self):
    return 'Adjust intensity errors -- ev11'

  def run(self, experiments, reflections):
    '''Modify intensity errors according to EV11 -- Brewster2019 / Mittan-Moreau 202X'''
    assert self.params.merging.error.model == "ev11"
    self.logger.log_step_time("ERROR_MODIFIER_EV11")
    self.logger.log("Modifying intensity errors -- ev11 method (starting with %d reflections)"%(len(reflections)))
    reflections = self.modify_errors(reflections)
    self.logger.log_step_time("ERROR_MODIFIER_EV11", True)
    return experiments, reflections

  def setup_work_arrays(self, reflections):
    '''Select multiply-measured HKLs.
    Calculate and cache reflection deltas, deltas squared, and HKL means for every reflection
    The mean_without_val for a reflection is the mean intensity of its hkl mates not including itself
    '''
    self.deltas     = flex.double()
    self.work_table = flex.reflection_table()
    I_wt            = flex.double()
    delta_all       = flex.double()
    delta_sq        = flex.double()
    mean            = flex.double() # mean = <I'_hj>
    biased_mean     = flex.double() # biased_mean = <I_h>, so dont leave out any reflection
    var             = flex.double()
    all_biased_mean = flex.double()
    correlation     = flex.double()
    # Need miller indices to scale by pairwise differences
    miller_index_asymmetric = flex.miller_index()

    for refls in reflection_table_utils.get_next_hkl_reflection_table(reflections):
      number_of_measurements = refls.size()
      # if the returned "refls" list is empty, it's the end of the input "reflections" list
      if number_of_measurements == 0:
        break
      refls_biased_mean = flex.double(len(refls), flex.mean(refls['intensity.sum.value']))
      all_biased_mean.extend(refls_biased_mean)
      if number_of_measurements > self.params.merging.minimum_multiplicity:
        I_wt.extend(refls['intensity.sum.value'])
        nn_factor_sqrt = math.sqrt((number_of_measurements - 1) / number_of_measurements)
        i_sum = flex.double(number_of_measurements, flex.sum(refls['intensity.sum.value']))
        i_sum_minus_val = i_sum - refls['intensity.sum.value']
        mean_without_val = i_sum_minus_val / (number_of_measurements - 1)
        delta = nn_factor_sqrt * (refls['intensity.sum.value'] - mean_without_val)
        delta_all.extend(delta)
        # Please be careful about where to put the var
        self.deltas.extend(delta/flex.sqrt(refls['intensity.sum.variance']))
        delta_sq.extend(delta**2)
        mean.extend(mean_without_val)
        biased_mean.extend(refls_biased_mean)
        var.extend(refls['intensity.sum.variance'])
        if self.algorithm == "ev11_mll":
          correlation.extend(refls[self.cc_key])
        else:
          correlation.extend(flex.double(number_of_measurements, 0))
        miller_index_asymmetric.extend(refls['miller_index_asymmetric'])

    self.work_table[self.cc_key]    = correlation
    self.work_table["delta_sq"]     = delta_sq
    self.work_table["delta_all"]    = delta_all
    self.work_table["mean"]         = mean
    self.work_table["biased_mean"]  = biased_mean
    self.work_table["intensity.sum.variance"] = var
    self.work_table["I"]            = I_wt
    self.work_table["miller_index_asymmetric"] = miller_index_asymmetric
    reflections['biased_mean']      = all_biased_mean

    self.logger.log("Number of work reflections selected: %d"%self.deltas.size())
    return reflections

  def calculate_delta_statistics(self):
    '''Calculate min, max, mean, and stddev for the normalized deltas'''
    delta_min = flex.min(self.deltas) if self.deltas.size() > 0 else float('inf')
    delta_max = flex.max(self.deltas) if self.deltas.size() > 0 else float('-inf')
    delta_sum = flex.sum(self.deltas) if self.deltas.size() > 0 else 0.0

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI

    # global min and max
    self.global_delta_min = comm.allreduce(delta_min, MPI.MIN)
    self.global_delta_max = comm.allreduce(delta_max, MPI.MAX)

    # global mean
    self.global_delta_count = comm.allreduce(self.deltas.size(), MPI.SUM)
    if self.global_delta_count < 20:
      raise ValueError("Too few reflections available for ev11 algorithm")
    global_delta_sum = comm.allreduce(delta_sum, MPI.SUM)
    self.global_delta_mean = global_delta_sum / self.global_delta_count

    # global standard deviation
    array_of_global_delta_means = flex.double(self.deltas.size(), self.global_delta_mean)
    array_of_diffs = self.deltas - array_of_global_delta_means
    array_of_square_diffs = array_of_diffs * array_of_diffs
    sum_of_square_diffs = flex.sum(array_of_square_diffs)
    global_sum_of_square_diffs = comm.allreduce(sum_of_square_diffs, MPI.SUM)
    self.global_delta_stddev = math.sqrt(global_sum_of_square_diffs / (self.global_delta_count - 1))
    if self.mpi_helper.rank == 0:
      self.logger.main_log(
        "Global delta statistics (count,min,max,mean,stddev): ("
        + f"{self.global_delta_count:d},"
        + f"{self.global_delta_min:f},"
        + f"{self.global_delta_max:f},"
        + f"{self.global_delta_mean:f},"
        + f"{self.global_delta_stddev:f})"
      )

  def calculate_delta_bin_limits(self):
    '''Divide the delta (min,max) range into "number of ranks" bins. For a balanced rank load,
       bin limits should be chosen so that the bins are equally populated by the deltas.
       Assuming the normal distribution of deltas, we use the probability density function
       for the bin calculations.'''
    from scipy.stats import norm
    cdf_min = norm.cdf(self.global_delta_min, loc=self.global_delta_mean, scale=self.global_delta_stddev)
    cdf_max = norm.cdf(self.global_delta_max, loc=self.global_delta_mean, scale=self.global_delta_stddev)
    self.delta_bin_limits = flex.double()
    for val in np.linspace(cdf_min, cdf_max, self.mpi_helper.size + 1):
      self.delta_bin_limits.append(norm.ppf(val, loc=self.global_delta_mean, scale=self.global_delta_stddev))
    # To fool-proof the binning, set the first and last bin limits to infinity
    self.delta_bin_limits[0]                                = float('-inf')
    self.delta_bin_limits[self.delta_bin_limits.size() - 1] = float('inf')

  def distribute_deltas_over_bins(self):
    '''Have each rank distribute its deltas over the global delta bins'''
    self.logger.log("Delta count: %d"%self.deltas.size())
    self.delta_bins = []
    for bin_begin in range(self.delta_bin_limits.size() - 1):
      test_1 = self.deltas >= flex.double(self.deltas.size(), self.delta_bin_limits[bin_begin])
      test_2 = self.deltas < flex.double(self.deltas.size(), self.delta_bin_limits[bin_begin + 1])
      d = self.deltas.select(test_1 & test_2)
      self.delta_bins.append(d)

    total_deltas_distributed = 0
    for delta_bin in self.delta_bins:
      total_deltas_distributed += delta_bin.size()
    self.logger.log("Total deltas distributed: %d"%total_deltas_distributed)

  def distribute_deltas_over_ranks(self):
    '''Use alltoall to accumulate all deltas of one delta bin at a single rank'''
    new_delta_bins = self.mpi_helper.comm.alltoall(self.delta_bins)

    self.deltas = flex.double()
    for delta_bin in new_delta_bins:
      self.deltas.extend(delta_bin)

    self.deltas = flex.sorted(self.deltas)

    self.logger.log("New deltas count: %d"%self.deltas.size())

  def calculate_delta_rankits(self):
    '''Implement expression (12) of Brewster2019'''
    # Get the base global index for this rank's deltas.
    # Example: if rank 0 has 10 deltas, the first delta on rank 1 will be the 10th global delta.
    from scitbx.math import distributions
    delta_count_per_rank = self.mpi_helper.comm.allreduce([self.deltas.size()])
    base_delta_index = sum(delta_count_per_rank[0:self.mpi_helper.rank])
    self.logger.log("Delta base index: %d"%base_delta_index)

    norm = distributions.normal_distribution()

    a = 3./8. if self.global_delta_count < 10. else 0.5

    self.rankits = flex.double()
    for i in range(self.deltas.size()):
      global_delta_index = base_delta_index + i
      rankit = norm.quantile((global_delta_index+1-a)/(self.global_delta_count+1-(2*a)))
      self.rankits.append(rankit)

  def calculate_initial_ev11_parameters(self):
    '''Do a global LS fit of deltas to rankits. Work only in the [0.5,0.5] range of rankits'''
    sum_xx = 0
    sum_yy = 0
    sum_xy = 0
    sum_x  = 0
    sum_y  = 0
    count = 0
    for delta,rankit in zip(self.deltas, self.rankits):
      if rankit >= -0.5 and rankit <= 0.5:
        sum_xx += delta ** 2
        sum_yy += rankit ** 2
        sum_xy += delta * rankit
        sum_x  += delta
        sum_y  += rankit
        count += 1

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    global_sum_xx = comm.reduce(sum_xx, MPI.SUM, root =0)
    global_sum_yy = comm.reduce(sum_yy, MPI.SUM, root =0)
    global_sum_xy = comm.reduce(sum_xy, MPI.SUM, root =0)
    global_sum_x  = comm.reduce(sum_x,  MPI.SUM, root =0)
    global_sum_y  = comm.reduce(sum_y,  MPI.SUM, root =0)
    global_count  = comm.reduce(count,  MPI.SUM, root =0)

    if self.mpi_helper.rank == 0:
      slope = 0
      offset = 0
      corr = 0
      try:
        DELTA = global_count * global_sum_xx - global_sum_x**2 # see p. 105 in Bevington & Robinson
        #assert abs(DELTA) > sys.float_info.epsilon, "Cannot initialize ev11 parameters"
        slope = (global_count * global_sum_xy - global_sum_x * global_sum_y) / DELTA
        offset = (global_sum_xx * global_sum_y - global_sum_x * global_sum_xy) / DELTA
      except ZeroDivisionError:
        pass
      self.logger.main_log("SLOPE: %f; OFFSET: %f"%(slope,offset))

      # Calculate initial EV11 parameters
      self.sfac = 1/slope
<<<<<<< HEAD
      self.sb = math.sqrt(offset)
      self.sadd = [0 for i in range(self.n_coefs)]
      self.sadd[0] = offset
=======
      self.sadd = offset
      self.sb = math.sqrt(self.sadd) if self.sadd > 0 else 0
>>>>>>> master

      '''
      if True:
        from matplotlib import pyplot as plt
        f = plt.figure(0)
        lim = -5, 5
        x = np.linspace(lim[0],lim[1],100) # 100 linearly spaced numbers
        y = slope * x + offset
        plt.plot(self.deltas, self.rankits, '-')

        #plt.plot(x,y)
        plt.title("CC: %.3f Slope: %.3f Offset: %.3f"%(corr, slope, offset))
        plt.xlabel("Sorted data")
        plt.ylabel("Rankits")
        plt.xlim(lim); plt.ylim(lim)
        plt.axes().set_aspect('equal')

        f = plt.figure(1)
        h = flex.histogram(self.deltas, n_slots=100, data_min = lim[0], data_max = lim[1])
        stats = flex.mean_and_variance(self.deltas)
        plt.plot(h.slot_centers().as_numpy_array(), h.slots().as_numpy_array(), '-')
        plt.xlim(lim)
        plt.xlabel("Sorted data")
        plt.ylabel("Count")
        plt.title("Normalized data mean: %.3f +/- %.3f"%(stats.mean(), stats.unweighted_sample_standard_deviation()))

        #plt.show()

        if True:
          plt.ion()
          plt.pause(0.05)
      '''

      initial_params = (self.sfac, self.sadd, self.sb)
    else:
      initial_params = None

    initial_params = self.mpi_helper.comm.bcast(initial_params, root=0)
    self.sfac = initial_params[0]
    self.sadd = initial_params[1]
    self.sb   = initial_params[2]

  def calculate_intensity_bin_limits(self):
    if self.algorithm == 'ev11':
      self.calculate_intensity_bin_limits_ev11()
    elif self.algorithm == 'ev11_mll':
      self.calculate_intensity_bin_limits_mll()

  def calculate_intensity_bin_limits_mll(self):
    '''Calculate the minimum and maximum values of the mean intensities for each HKL
    This is specific to the updated maximum log-likelihood error model (Mittan-Moreau 202X).'''
    senddata = self.work_table['biased_mean'].as_numpy_array()
    sendcounts = np.array(self.mpi_helper.comm.gather(len(senddata), root=0))
    if self.mpi_helper.rank == 0:
      all_mean_intensities = np.empty(sum(sendcounts), dtype=senddata.dtype)
    else:
      all_mean_intensities = None
    self.mpi_helper.comm.Gatherv(senddata, (all_mean_intensities, sendcounts), root=0)
    if self.mpi_helper.rank == 0:
      all_mean_intensities = np.sort(all_mean_intensities)
      lower_percentile = 0.005
      upper_percentile = 0.995
      n = all_mean_intensities.size
      lower = all_mean_intensities[int(lower_percentile * n)]
      upper = all_mean_intensities[int(upper_percentile * n)]
      self.intensity_bin_limits = np.linspace(lower, upper, number_of_intensity_bins + 1)

    else:
      self.intensity_bin_limits = None
    self.intensity_bin_limits = self.mpi_helper.comm.bcast(self.intensity_bin_limits, root=0)

  def calculate_intensity_bin_limits_ev11(self):
    '''Calculate the minimum and maximum values of the mean intensities for each HKL
    This is bin limit calculation in the original Ev11 implementation (Brewster 2019)'''
    count = self.work_table.size()
    mean_intensity_min = flex.min(self.work_table['biased_mean']) if count > 0 else float('inf')
    mean_intensity_max = flex.max(self.work_table['biased_mean']) if count > 0 else float('-inf')
    if count > 0:
      self.logger.log("Using %d multiply-measured HKLs; mean intensity (min,max): (%f,%f)"%(count, mean_intensity_min, mean_intensity_max))
    else:
      self.logger.log("No multiply-measured HKLs available")

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    global_mean_intensity_min = comm.allreduce(mean_intensity_min, MPI.MIN)
    global_mean_intensity_max = comm.allreduce(mean_intensity_max, MPI.MAX)
    self.logger.log("Global mean intensity (min,max): (%f,%f)"%(global_mean_intensity_min, global_mean_intensity_max))

    self.intensity_bin_limits = np.linspace(global_mean_intensity_min, global_mean_intensity_max, number_of_intensity_bins + 1)
    self.intensity_bin_limits[0] = float('-inf')
    self.intensity_bin_limits[len(self.intensity_bin_limits) - 1] = float('inf')

  def distribute_reflections_over_intensity_bins(self):
    self.intensity_bins = []
    count = self.work_table.size()
    for bin_begin in range(number_of_intensity_bins):
      self.intensity_bins.append(flex.reflection_table())
      test_1 = self.work_table['biased_mean'] >= flex.double(count, self.intensity_bin_limits[bin_begin])
      test_2 = self.work_table['biased_mean'] < flex.double(count, self.intensity_bin_limits[bin_begin + 1])
      sel = self.work_table.select(test_1 & test_2)
      self.intensity_bins[bin_begin].extend(sel)

    # for debugging
    number_of_refls_distributed = 0
    for intensity_bin in self.intensity_bins:
      number_of_refls_distributed += intensity_bin.size()
    self.logger.log("Distributed over intensity bins %d out of %d reflections"%(number_of_refls_distributed, count))

  def setup_pairwise_differences(self):
    rng = np.random.default_rng(seed=self.mpi_helper.rank)
    self.pairwise_differences = [flex.reflection_table() for i in range(number_of_intensity_bins)]
    self.n_refls_in_bin = [0 for i in range(number_of_intensity_bins)]
    for bin_index, bin_reflections in enumerate(self.intensity_bins):
      n_refls_in_bin_rank = 0
      if len(bin_reflections) > 0:
        pairwise_differences = flex.double()
        biased_mean = flex.double()
        var_cs_i = flex.double()
        var_cs_j = flex.double()
        correlation_i = flex.double()
        correlation_j = flex.double()
        miller_index_asymmetric = flex.miller_index()

        for same_reflections in reflection_table_utils.get_next_hkl_reflection_table(bin_reflections):
          number_of_reflections = same_reflections.size()
          n_refls_in_bin_rank += number_of_reflections
          if number_of_reflections > 1:
            I = same_reflections['I'].as_numpy_array()
            var_cs = same_reflections['intensity.sum.variance'].as_numpy_array()
            correlation = same_reflections[self.cc_key].as_numpy_array()

            indices = np.triu_indices(n=I.size, k=1)
            N = indices[0].size
            if self.limit_differences and N > self.n_max_differences:
              subset_indices = rng.choice(N, size=self.n_max_differences, replace=False, shuffle=True)
              indices = (indices[0][subset_indices], indices[1][subset_indices])
              N = self.n_max_differences
            differences = np.abs((I[np.newaxis] - I[:, np.newaxis])[indices])
            v_cs_i = var_cs[indices[0]]
            v_cs_j = var_cs[indices[1]]
            cc_i = correlation[indices[0]]
            cc_j = correlation[indices[1]]

            pairwise_differences.extend(flex.double(differences))
            biased_mean.extend(flex.double(N, same_reflections['biased_mean'][0]))
            var_cs_i.extend(flex.double(v_cs_i))
            var_cs_j.extend(flex.double(v_cs_j))
            correlation_i.extend(flex.double(cc_i))
            correlation_j.extend(flex.double(cc_j))
            miller_index_asymmetric.extend(flex.miller_index(
              N, same_reflections['miller_index_asymmetric'][0]
              ))
        self.pairwise_differences[bin_index]['pairwise_differences'] = pairwise_differences
        self.pairwise_differences[bin_index]['biased_mean'] = biased_mean
        self.pairwise_differences[bin_index]['var_cs_i'] = var_cs_i
        self.pairwise_differences[bin_index]['var_cs_j'] = var_cs_j
        self.pairwise_differences[bin_index]['correlation_i'] = correlation_i
        self.pairwise_differences[bin_index]['correlation_j'] = correlation_j
        self.pairwise_differences[bin_index]['miller_index_asymmetric'] = miller_index_asymmetric
      self.n_refls_in_bin[bin_index] = self.mpi_helper.comm.reduce(
          n_refls_in_bin_rank,
          op=self.mpi_helper.MPI.SUM, root=0
          )

  def run_minimizer(self):
    from scitbx import lbfgsb

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    size = self.mpi_helper.size
    if self.tuning_param_opt:
      param_offset = 3
      param_shift = 1
      self.x = flex.double([self.tuning_param, self.sfac, self.sb, *self.sadd])
    else:
      param_offset = 2
      param_shift = 0
      self.x = flex.double([self.sfac, self.sb, *self.sadd])
    self.n = self.n_coefs + param_offset
    if self.mpi_helper.rank == 0:
      self.logger.main_log(
        'Initial Parameter Estimates = '
        + f'sfac: {self.sfac} '
        + f'sb: {self.sb} '
        + f'sadd: {self.sadd[0]} '
        + f'nu: {self.tuning_param} '
      )
    l = flex.double(self.n, 1e-8)
    u = flex.double(self.n, 0)
    if self.tuning_param_opt:
      # normalization for the truncated t-distribution is numerically unstable for nu < 2
      l[0] = 2
      if self.x[0] < 2:
        self.x[0] = 2
    for degree_index in range(self.n_coefs):
      l[degree_index + param_offset] = -1000
    if self.mpi_helper.rank == 0:
      self.minimizer = lbfgsb.minimizer(
        n = self.n,
        l = l,
        u = u,
        nbd = flex.int(self.n, 1),
      )
    #self.verify_derivatives()
    while True:
      self.compute_functional_and_gradients()
      status=-1
      if self.mpi_helper.rank == 0:
        if self.minimizer.process(self.x, self.f, self.g):
          if self.tuning_param_opt:
            self.tuning_param = self.x[0]
            tuning_param = f'{self.tuning_param:0.3f}'
          self.sfac = self.x[0 + param_shift]
          self.sb   = self.x[1 + param_shift]
          self.sadd = self.x[2 + param_shift:]
          sfac = f'{self.sfac:0.3f}'
          sb = f'{self.sb:0.6f}'
          sadd = [f'{self.sadd[i]:0.3f}' for i in range(self.n_coefs)]
          log_out = 'intermediate minimization results = '\
            + f'functional: {self.f:.2f} '\
            + f'sfac: {sfac} '\
            + f'sb: {sb} '\
            + f'sadd: {sadd} '
          if self.tuning_param_opt:
            log_out += f'nu: {tuning_param}'
          self.logger.main_log(log_out)
          status=1
        elif self.minimizer.is_terminated():
          status=0

      comm.barrier()
      status=comm.bcast(status, root=0)
      if status==1:
        self.tuning_param=comm.bcast(self.tuning_param, root=0)
        self.sfac=comm.bcast(self.sfac, root=0)
        self.sb=comm.bcast(self.sb, root=0)
        self.sadd=comm.bcast(self.sadd, root=0)
        pass
      if status==0:
        break

    if self.mpi_helper.rank == 0:
      tuning_param = f'{self.tuning_param:0.3f}'
      sfac = f'{self.sfac:0.3f}'
      sb = f'{self.sb:0.3f}'
      sadd = [f'{self.sadd[i]:0.3f}' for i in range(self.n_coefs)]
      log_out = 'FINAL EV11 VALUES = '\
        + f'functional: {self.f:.2f} '\
        + f'sfac: {sfac} '\
        + f'sb: {sb} '\
        + f'sadd: {sadd} '
      if self.tuning_param_opt:
        log_out += f'nu: {tuning_param}'
      self.logger.main_log(log_out)
    if self.algorithm == 'ev11_mll':
      self.scale_sfac()
    if self.params.merging.error.ev11.do_diagnostics:
      #self.verify_derivatives()
      self.plot_diagnostics()

  def scale_sfac(self):
    scale = self._get_scale_pairwise_difference()
    if self.mpi_helper.rank == 0:
      if self.overall_scaling_method == 'pairwise_differences':
        sfac_old = self.sfac
        self.sfac = self.sfac * scale
        log_out = 'SCALING SFAC  '\
          + f'Unscaled standard deviation: {scale:0.3f} '\
          + f'Optimized sfac: {sfac_old:0.3f} '\
          + f'Scaled sfac: {self.sfac:0.3f} '
      else:
        log_out = f'Unscaled standard deviation: {scale:0.3f} '
      self.logger.main_log(log_out)
    if self.overall_scaling_method == 'pairwise_differences':
      self.sfac = self.mpi_helper.comm.bcast(self.sfac, root=0)

  def _get_scale_pairwise_difference(self, return_in_bins=False):
    def min_fun_t(c, df):
      term0 = scipy.stats.t.ppf(3/4, df=df)
      term1 = scipy.stats.t.cdf(term0 + 1/c, df=df)
      term2 = scipy.stats.t.cdf(term0 - 1/c, df=df)
      term3 = 1/2
      diff = term1 - term2 - term3
      return diff**2

    if self.mpi_helper.rank == 0:
      if self.likelihood == 'normal':
        conversion_factor = 1.1926
      elif self.likelihood == 't-dist':
        results = scipy.optimize.minimize_scalar(
          min_fun_t,
          bounds=(0.1, 2),
          args=(self.tuning_param)
          )
        conversion_factor = results.x
    else:
      conversion_factor = None
    conversion_factor = self.mpi_helper.comm.bcast(conversion_factor)

    median_differences = []
    if return_in_bins:
      median_differences_bin = np.zeros(number_of_intensity_bins)
    for bin_index, bin_reflections in enumerate(self.intensity_bins):
      if len(bin_reflections) > 0:
        for same_reflections in reflection_table_utils.get_next_hkl_reflection_table(bin_reflections):
          number_of_reflections = same_reflections.size()
          if number_of_reflections > 1:
            I = same_reflections['I'].as_numpy_array()
            var_ev11 = self._get_var_ev11(
              same_reflections['intensity.sum.variance'],
              same_reflections['biased_mean'],
              same_reflections[self.cc_key]
              )
            std_estimate = np.median(np.sqrt(var_ev11.as_numpy_array()))
            differences = np.abs(I[np.newaxis] - I[:, np.newaxis])
            differences[np.diag_indices(number_of_reflections)] = np.nan
            median_differences.append(
              np.median(np.nanmedian(differences, axis=1)) / std_estimate
              )

      if return_in_bins:
        all_median_differences_bin = self.mpi_helper.comm.gather(median_differences, root=0)
        median_differences = []
        if self.mpi_helper.rank == 0:
          all_median_differences_bin = np.concatenate(all_median_differences_bin)
          median_differences_bin[bin_index] = np.median(all_median_differences_bin) * conversion_factor

    if return_in_bins:
      return median_differences_bin
    else:
      all_median_differences = self.mpi_helper.comm.gather(median_differences, root=0)
      if self.mpi_helper.rank == 0:
        all_median_differences = np.concatenate(all_median_differences)
        scale = np.median(all_median_differences) * conversion_factor
        return float(scale)
      else:
        return None

  def _get_scale_pairwise_difference_nd(self, return_in_bins=False):
    def min_fun_t(c, df):
      term0 = scipy.stats.t.ppf(3/4, df=df)
      term1 = scipy.stats.t.cdf(term0 + 1/c, df=df)
      term2 = scipy.stats.t.cdf(term0 - 1/c, df=df)
      term3 = 1/2
      diff = term1 - term2 - term3
      return diff**2

    if self.mpi_helper.rank == 0:
      if self.likelihood == 'normal':
        conversion_factor = 1.1926
      elif self.likelihood == 't-dist':
        results = scipy.optimize.minimize_scalar(
          min_fun_t,
          bounds=(0.1, 2),
          args=(self.tuning_param)
          )
        conversion_factor = results.x
    else:
      conversion_factor = None
    conversion_factor = self.mpi_helper.comm.bcast(conversion_factor)

    median_differences = []
    if return_in_bins:
      median_differences_bin = np.zeros(number_of_intensity_bins)
    for bin_index, bin_reflections in enumerate(self.intensity_bins):
      if len(bin_reflections) > 0:
        for same_reflections in reflection_table_utils.get_next_hkl_reflection_table(bin_reflections):
          number_of_reflections = same_reflections.size()
          if number_of_reflections > 0:
            var_ev11 = self._get_var_ev11(
              same_reflections['intensity.sum.variance'],
              same_reflections['biased_mean'],
              same_reflections[self.cc_key]
              )
            normalized_deviations = (same_reflections['delta_all'] / flex.sqrt(var_ev11)).as_numpy_array()
            differences = normalized_deviations[np.newaxis] - normalized_deviations[:, np.newaxis]
            median_differences.append(np.median(np.abs(differences), axis=1))
      if return_in_bins:
        if len(median_differences) > 0:
          median_differences = np.concatenate(median_differences)
        all_median_differences_bin = self.mpi_helper.comm.gather(median_differences, root=0)
        median_differences = []
        if self.mpi_helper.rank == 0:
          all_median_differences_bin = np.concatenate(all_median_differences_bin)
          median_differences_bin[bin_index] = np.median(all_median_differences_bin) * conversion_factor

    if return_in_bins:
      return median_differences_bin
    else:
      all_median_differences = self.mpi_helper.comm.gather(np.concatenate(median_differences), root=0)
      if self.mpi_helper.rank == 0:
        all_median_differences = np.concatenate(all_median_differences)
        scale = np.median(all_median_differences) * conversion_factor
        return float(scale)
      else:
        return None

  def compute_functional_and_gradients(self):
    if self.algorithm == 'ev11':
      self.calculate_functional_ev11()
    elif self.algorithm == 'ev11_mll':
      self.calculate_functional_mll()
    self.f = self.functional
    if self.tuning_param_opt:
      self.g = flex.double([self.der_wrt_nu, self.der_wrt_sfac, self.der_wrt_sb, *self.der_wrt_sadd])
    else:
      self.g = flex.double([self.der_wrt_sfac, self.der_wrt_sb, *self.der_wrt_sadd])
    return self.f, self.g

  def verify_derivatives(self):
    shift = 0.000001
    import copy
    if self.tuning_param_opt:
      self.n = self.n_coefs + 3
    else:
      self.n = self.n_coefs + 2
    if self.algorithm == 'ev11':
      tf_func = self.calculate_functional_ev11
    elif self.algorithm == 'ev11_mll':
      tf_func = self.calculate_functional_mll

    tf_func()
    TF = copy.copy(self.functional)
    der_wrt_sfac = copy.copy(self.der_wrt_sfac)
    der_wrt_sb = copy.copy(self.der_wrt_sb)
    der_wrt_sadd = copy.copy(self.der_wrt_sadd)
    sfac = copy.copy(self.sfac)
    sb = copy.copy(self.sb)
    sadd = copy.copy(self.sadd)

    # Tuning param
    if self.tuning_param_opt:
      der_wrt_nu = copy.copy(self.der_wrt_nu)
      tuning_param = copy.copy(self.tuning_param)
      self.tuning_param = tuning_param * (1 + shift)
      tf_func()
      TF_p = copy.copy(self.functional)
      self.tuning_param = tuning_param * (1 - shift)
      tf_func()
      TF_m = copy.copy(self.functional)
      check_der_wrt_nu = (TF_p - TF_m) / (2 * shift * tuning_param)
      self.tuning_param = tuning_param
      if self.mpi_helper.rank == 0:
        print(f'der_wrt_nu numerical: {check_der_wrt_nu} analytical {der_wrt_nu}')

    # sfac
    self.sfac = sfac * (1 + shift)
    tf_func()
    TF_p = copy.copy(self.functional)
    self.sfac = sfac * (1 - shift)
    tf_func()
    TF_m = copy.copy(self.functional)
    check_der_wrt_sfac = (TF_p - TF_m) / (2 * shift * sfac)
    self.sfac = sfac
    if self.mpi_helper.rank == 0:
      print(f'der_wrt_sfac numerical: {check_der_wrt_sfac} analytical {der_wrt_sfac}')

    # sb
    self.sb = sb * (1 + shift)
    tf_func()
    TF_p = copy.copy(self.functional)
    self.sb = sb * (1 - shift)
    tf_func()
    TF_m = copy.copy(self.functional)
    check_der_wrt_sb = (TF_p - TF_m) / (2 * shift * sb)
    self.sb = sb
    if self.mpi_helper.rank == 0:
      print(f'der_wrt_sb numerical: {check_der_wrt_sb} analytical {der_wrt_sb}')

    # sadd:
    for degree_index in range(self.n_coefs):
      if sadd[degree_index] == 0:
        self.sadd[degree_index] = shift
      else:
        self.sadd[degree_index] = sadd[degree_index] * (1 + shift)
      tf_func()
      TF_p = copy.copy(self.functional)
      if sadd[degree_index] == 0:
        self.sadd[degree_index] = -shift
      else:
        self.sadd[degree_index] = sadd[degree_index] * (1 - shift)
      tf_func()
      TF_m = copy.copy(self.functional)
      if sadd[degree_index] == 0:
        check_der_wrt_sadd = (TF_p - TF_m) / (2 * shift)
      else:
        check_der_wrt_sadd = (TF_p - TF_m) / (2 * shift * sadd[degree_index])
      self.sadd[degree_index] = sadd[degree_index]
      if self.mpi_helper.rank == 0:
        print(f'der_wrt_sadd - degree {degree_index} numerical: {check_der_wrt_sadd} analytical {der_wrt_sadd[degree_index]}')
    return None

  def plot_diagnostics(self):
    def get_rankits(n, down_sample, distribution):
      prob_level = (np.arange(1, n+1) - 0.5) / n
      if distribution == 'normal':
        return scipy.stats.norm.ppf(prob_level[::down_sample])
      elif distribution == 'half normal':
        return scipy.stats.halfnorm.ppf(prob_level[::down_sample])
      elif distribution == 't-dist':
        return scipy.stats.t.ppf(prob_level[::down_sample], df=self.tuning_param)
      elif distribution == 'half t-dist':
        prob_level = (prob_level + 1) / 2
        return scipy.stats.t.ppf(prob_level[::down_sample], df=self.tuning_param)

    cc_rank = self.work_table[self.cc_key].as_numpy_array()
    cc_all = self.mpi_helper.comm.gather(cc_rank, root=0)

    output_nd_keys = ['deviations', 'biased_mean', 'var_ev11', 'I', 'intensity.sum.variance']
    output_nd_rank = {}
    for key in output_nd_keys:
      output_nd_rank[key] = flex.double()
    for reflections in self.intensity_bins:
      number_of_reflections_in_bin = reflections.size()
      if number_of_reflections_in_bin > 0:
        var_ev11 = self._get_var_ev11(
          reflections['intensity.sum.variance'],
          reflections['biased_mean'],
          reflections[self.cc_key]
          )
        output_nd_rank['deviations'].extend(reflections['delta_all'])
        output_nd_rank['biased_mean'].extend(reflections['biased_mean'])
        output_nd_rank['var_ev11'].extend(var_ev11)
        output_nd_rank['I'].extend(reflections['I'])
        output_nd_rank['intensity.sum.variance'].extend(reflections['intensity.sum.variance'])
    output_nd_all = {}
    for key in output_nd_keys:
      output_nd_all[key] = self.mpi_helper.comm.gather(output_nd_rank[key].as_numpy_array(), root=0)

    output_pd_keys = ['differences', 'var_ev11_i', 'var_ev11_j', 'biased_mean']
    output_pd_rank = {}
    for key in output_pd_keys:
      output_pd_rank[key] = flex.double()
    if self.algorithm == 'ev11':
      self.setup_pairwise_differences()
    for differences in self.pairwise_differences:
      number_of_differences_in_bin = differences.size()
      if number_of_differences_in_bin > 0:
        var_i = self._get_var_ev11(
          differences['var_cs_i'],
          differences['biased_mean'],
          differences['correlation_i'],
          )
        var_j = self._get_var_ev11(
          differences['var_cs_j'],
          differences['biased_mean'],
          differences['correlation_j'],
          )
        output_pd_rank['differences'].extend(differences['pairwise_differences'])
        output_pd_rank['var_ev11_i'].extend(var_i)
        output_pd_rank['var_ev11_j'].extend(var_j)
        output_pd_rank['biased_mean'].extend(differences['biased_mean'])
    output_pd_all = {}
    for key in output_pd_keys:
      output_pd_all[key] = self.mpi_helper.comm.gather(output_pd_rank[key].as_numpy_array(), root=0)

    sigma_pairwise_differences_nd = self._get_scale_pairwise_difference_nd(return_in_bins=True)
    sigma_pairwise_differences = self._get_scale_pairwise_difference(return_in_bins=True)
    if self.mpi_helper.rank == 0:
      import matplotlib.pyplot as plt
      import pandas as pd
      for key in output_nd_keys:
        output_nd_all[key] = np.concatenate(output_nd_all[key])
      for key in output_pd_keys:
        output_pd_all[key] = np.concatenate(output_pd_all[key])
      lim = 5
      downsample = 10
      normalized_deviations = output_nd_all['deviations'] / np.sqrt(output_nd_all['var_ev11'])
      sorted_normalized_deviations = np.sort(normalized_deviations)
      normalized_deviations_bins = np.linspace(-lim, lim, 101)
      normalized_deviations_db = normalized_deviations_bins[1] - normalized_deviations_bins[0]
      normalized_deviations_centers = (normalized_deviations_bins[1:] + normalized_deviations_bins[:-1]) / 2
      normalized_deviations_hist, _ = np.histogram(
        sorted_normalized_deviations, bins=normalized_deviations_bins, density=True
        )
      nd_rankits_normal = get_rankits(sorted_normalized_deviations.size, downsample, 'normal')
      nd_rankits_t_dist = get_rankits(sorted_normalized_deviations.size, downsample, 't-dist')

      var = output_pd_all['var_ev11_i'] + output_pd_all['var_ev11_j']
      normalized_pd = output_pd_all['differences'] / np.sqrt(var)
      sorted_normalized_pd = np.sort(normalized_pd)
      normalized_pd_bins = np.linspace(0, lim, 101)
      normalized_pd_db = normalized_pd_bins[1] - normalized_pd_bins[0]
      normalized_pd_centers = (normalized_pd_bins[1:] + normalized_pd_bins[:-1]) / 2
      normalized_pd_hist, _ = np.histogram(
        sorted_normalized_pd, bins=normalized_pd_bins, density=True
        )
      npd_rankits_half_normal = get_rankits(sorted_normalized_pd.size, downsample, 'half normal')
      npd_rankits_half_t_dist = get_rankits(sorted_normalized_pd.size, downsample, 'half t-dist')

      intensity_centers = (self.intensity_bin_limits[1:] + self.intensity_bin_limits[:-1]) / 2
      sigma_bin = np.zeros((number_of_intensity_bins, 2, 3))
      sigma_bin[:, 0, 1] = sigma_pairwise_differences_nd
      sigma_bin[:, 1, 1] = sigma_pairwise_differences
      sigma_bin[:, 0, 2] = intensity_centers
      sigma_bin[:, 1, 2] = intensity_centers
      for index in range(number_of_intensity_bins):
          indices = np.logical_and(
            output_nd_all['biased_mean'] >= self.intensity_bin_limits[index],
            output_nd_all['biased_mean'] < self.intensity_bin_limits[index + 1]
            )
          if np.count_nonzero(indices) > 0:
            nd = np.sort(normalized_deviations[indices])
            sigma_bin[index, 0, 0] = nd.std()
          else:
            sigma_bin[index, 0, :] = np.nan

          indices = np.logical_and(
            output_pd_all['biased_mean'] >= self.intensity_bin_limits[index],
            output_pd_all['biased_mean'] < self.intensity_bin_limits[index + 1]
            )
          if np.count_nonzero(indices) > 0:
            differences = np.sort(normalized_pd[indices])
            sigma_bin[index, 1, 0] = differences.std()
          else:
            sigma_bin[index, 1, :] = np.nan

      np.save(
        os.path.join(self.params.output.output_dir, self.params.output.prefix + f'_sigma_bin.npy'),
        sigma_bin
        )
      # Normalized Deviations Plots
      fig, axes = plt.subplots(2, 3, figsize=(8, 5))
      I_scale = 100000
      grey1 = np.array([99, 102, 106]) / 255
      grey2 = np.array([177, 179, 179]) / 255
      axes[0, 0].bar(
        normalized_deviations_centers, normalized_deviations_hist,
        width=normalized_deviations_db, label='$\delta_{hkl}$'
        )
      axes[0, 0].plot(
        normalized_deviations_centers,
        scipy.stats.norm.pdf(normalized_deviations_centers),
        color=grey1, label='Normal'
        )
      if self.likelihood == 't-dist':
        axes[0, 0].plot(
          normalized_deviations_centers,
          scipy.stats.t.pdf(normalized_deviations_centers, df=self.tuning_param),
          color=grey2, label=f't-dist\n$\\nu: ${self.tuning_param:0.1f}'
          )

      axes[1, 0].bar(
        normalized_pd_centers, normalized_pd_hist,
        width=normalized_pd_db, label='$DIFF_{hkl}$'
        )
      axes[1, 0].plot(
        normalized_pd_centers,
        scipy.stats.halfnorm.pdf(normalized_pd_centers),
        color=grey1, label='Normal'
        )
      if self.likelihood == 't-dist':
        axes[1, 0].plot(
          normalized_pd_centers,
          2*scipy.stats.t.pdf(normalized_pd_centers, df=self.tuning_param),
          color=grey2, label=f't-dist\n$\\nu: ${self.tuning_param:0.1f}'
          )

      axes[0, 0].legend(frameon=False, fontsize=8, handlelength=1)
      axes[0, 0].set_ylabel('Distribution of $\delta_{hbk}$')
      axes[0, 0].set_xlabel('Normalized Deviations ($\delta_{hbk}$)')
      axes[0, 0].set_xlim([-4.5, 4.5])
      axes[0, 0].set_xticks([-4, -2, 0, 2, 4])
      axes[1, 0].set_ylabel('Distribution of $DIFF_{hbk}$')
      axes[1, 0].set_xlabel('Normalized PD ($DIFF_{hbk}$)')
      axes[1, 0].set_xlim([0, 4.5])
      axes[1, 0].set_xticks([0, 1, 2, 3, 4])

      axes[0, 1].plot([-lim, lim], [-lim, lim], color=[0, 0, 0], linewidth=1, linestyle=':')
      axes[0, 1].plot(
        sorted_normalized_deviations[::downsample],
        nd_rankits_normal,
        color=grey1
        )
      if self.likelihood == 't-dist':
        axes[0, 1].plot(
          sorted_normalized_deviations[::downsample],
          nd_rankits_t_dist,
          color=grey2
          )

      axes[1, 1].plot([0, lim], [0, lim], color=[0, 0, 0], linewidth=1, linestyle=':')
      axes[1, 1].plot(
        sorted_normalized_pd[::downsample],
        npd_rankits_half_normal,
        color=grey1
        )
      if self.likelihood == 't-dist':
        axes[1, 1].plot(
          sorted_normalized_pd[::downsample],
          npd_rankits_half_t_dist,
          color=grey2
          )

      axes[0, 1].set_ylim([-lim, lim])
      axes[0, 1].set_ylabel('Rankits')
      axes[0, 1].set_xlabel('Sorted Normalized Deviations ($\delta_{hbk}$)')
      axes[0, 1].set_box_aspect(1)
      axes[0, 1].set_xticks([-4, -2, 0, 2, 4])
      axes[0, 1].set_yticks([-4, -2, 0, 2, 4])
      axes[0, 1].set_xlim([-4.5, 4.5])
      axes[0, 1].set_ylim([-4.5, 4.5])

      axes[1, 1].set_ylim([0, lim])
      axes[1, 1].set_ylabel('Rankits')
      axes[1, 1].set_xlabel('Sorted Normalized PD ($DIFF_{hbk}$)')
      axes[1, 1].set_box_aspect(1)
      axes[1, 1].set_xticks([0, 1, 2, 3, 4])
      axes[1, 1].set_yticks([0, 1, 2, 3, 4])
      axes[1, 1].set_xlim([0, 4.5])
      axes[1, 1].set_ylim([0, 4.5])

      x = intensity_centers / I_scale
      axes[0, 2].plot([x[0], x[-1]], [1, 1], color=[0, 0, 0])
      axes[0, 2].plot(x, sigma_bin[:, 0, 0], label='STD')
      axes[0, 2].plot(x, sigma_bin[:, 0, 1], label='Pairwise')
      axes[0, 2].set_ylabel('Standard Deviation of $\delta_{hbk}$')
      axes[0, 2].set_xlabel('Mean Intensity X 100,000')
      axes[0, 2].legend(frameon=False, labelspacing=0.1, handlelength=1)

      axes[1, 2].plot([x[0], x[-1]], [1, 1], color=[0, 0, 0])
      if self.likelihood == 't-dist':
        v = self.tuning_param
        expected_sigma = np.sqrt(v / (v - 2))
        axes[0, 2].plot([x[0], x[-1]], [expected_sigma, expected_sigma], color=grey2)
        axes[1, 2].plot([x[0], x[-1]], [expected_sigma, expected_sigma], color=grey2)

      axes[1, 2].plot(x, sigma_bin[:, 1, 0], label='STD')
      axes[1, 2].plot(x, sigma_bin[:, 1, 1], label='Pairwise')
      axes[1, 2].set_ylabel('Standard Deviation of $DIFF_{hbk}$')
      axes[1, 2].set_xlabel('Mean Intensity X 100,000')

      fig.tight_layout()
      fig.savefig(os.path.join(
        self.params.output.output_dir,
        self.params.output.prefix + '_NormalizedDeviations.png'
        ))
      plt.close()

      df = pd.DataFrame(output_nd_all)
      df.to_parquet(os.path.join(
        self.params.output.output_dir,
        self.params.output.prefix + f'_deltas.parquet'
        ))

      df = pd.DataFrame(output_pd_all)
      df.to_parquet(os.path.join(
        self.params.output.output_dir,
        self.params.output.prefix + f'_differences.parquet'
        ))

      # CC & sadd plots #
      cc_all = np.concatenate(cc_all)
      cc_unique = np.unique(cc_all)
      bins = np.linspace(cc_unique.min(), cc_unique.max(), 101)
      dbin = bins[1] - bins[0]
      centers = (bins[1:] + bins[:-1]) / 2
      hist_all, _ = np.histogram(cc_unique, bins=bins)

      hist_color = np.array([0, 49, 60]) / 256
      line_color = np.array([213, 120, 0]) / 256
      sadd, _ = self._get_sadd(flex.double(centers))
      fig, axes_hist = plt.subplots(1, 1, figsize=(5, 3))
      axes_sadd = axes_hist.twinx()
      axes_hist.bar(centers, hist_all / 1000, width=dbin, color=hist_color)
      axes_sadd.plot(centers, (self.sfac * sadd)**2, color=line_color)
      axes_hist.set_xlabel('Correlation Coefficient')
      axes_hist.set_ylabel('Lattices (x1,000)')
      axes_sadd.set_ylabel('$s_{\mathrm{fac}}^2 \\times s_{\mathrm{add}}^2$')
      fig.tight_layout()
      fig.savefig(os.path.join(
        self.params.output.output_dir,
        self.params.output.prefix + '_sadd.png'
        ))
      plt.close()

      # I/sigma plots
      fig, axes = plt.subplots(1, 2, figsize=(7, 3), sharex=True, sharey=True)
      x = output_nd_all['I']
      y0 = output_nd_all['I'] / np.sqrt(output_nd_all['intensity.sum.variance'])
      y1 = output_nd_all['I'] / np.sqrt(output_nd_all['var_ev11'])
      I_s_asy =  1 / (self.sfac * flex.abs(sadd))
      for index in range(sadd.size()):
        axes[1].plot(
          [x.min(), x.max()], [I_s_asy[index], I_s_asy[index]],
          alpha=hist_all[index] / hist_all.max(),
          color=[0.8, 0, 0],
          linewidth=0.2
          )
      axes[0].plot(x[::10], y0[::10], marker='.', linestyle='none', markersize=0.1, alpha=0.5, color=[0, 0, 0])
      axes[1].plot(x[::10], y1[::10], marker='.', linestyle='none', markersize=0.1, alpha=0.5, color=[0, 0, 0])
      axes[0].set_xscale('log')
      axes[0].set_yscale('log')
      axes[0].set_xlabel('Intensity')
      axes[1].set_xlabel('Intensity')
      axes[0].set_ylabel('Intensity / $\sigma$')
      axes[0].set_title('Only Counting Statistics Error')
      axes[1].set_title('Ev11 Error')
      fig.tight_layout()
      fig.savefig(os.path.join(
        self.params.output.output_dir,
        self.params.output.prefix + '_I_sigma_plots.png'
        ))
      plt.close()
    return None

  def _loss_function_gaus(self, differences, var_i, var_j):
    var = var_i + var_j
    dvar_dvar_x = 1
    z = differences / flex.sqrt(var)
    dz_dvar = -differences / (2 * var**(3/2))
    L1 = 1/2*flex.log(var)
    dL1_dvar = 1/2 * 1/var
    L2 = 1/2 * z**2
    dL2_dz = z
    L = L1 + L2
    dL_dvar = dL1_dvar + dL2_dz * dz_dvar
    dL_dvar_i = dL_dvar * dvar_dvar_x
    dL_dvar_j = dL_dvar * dvar_dvar_x
    return L, dL_dvar_i, dL_dvar_j

  def _loss_function_t(self, differences_flex, var_i_flex, var_j_flex):
    v = self.tuning_param
    differences = differences_flex.as_numpy_array()
    var = (var_i_flex + var_j_flex).as_numpy_array()
    dvar_dvar_x = 1
    z = differences / np.sqrt(var)
    dz_dvar = -differences / (2 * var**(3/2))
    arg = 1 + 1/v * z**2
    darg_dz = 2*z/v

    L1 = 1/2 * np.log(var)
    dL1_dvar = 1/2 * 1/var
    L2 = (v+1)/2 * np.log(arg)
    dL2_darg = (v+1)/2 * 1/arg
    dL2_dvar = dL2_darg * darg_dz * dz_dvar
    L = L1 + L2
    dL_dvar_i = (dL1_dvar + dL2_dvar) * dvar_dvar_x
    dL_dvar_j = (dL1_dvar + dL2_dvar) * dvar_dvar_x
    return flex.double(L), flex.double(dL_dvar_i), flex.double(dL_dvar_j)

  def _loss_function_t_opt(self, differences_flex, var_i_flex, var_j_flex):
    v = self.tuning_param
    differences = differences_flex.as_numpy_array()
    var = (var_i_flex + var_j_flex).as_numpy_array()
    dvar_dvar_x = 1
    z = differences / np.sqrt(var)
    dz_dvar = -differences / (2 * var**(3/2))
    arg = 1 + 1/v * z**2
    darg_dz = 2*z/v
    darg_dv = -z**2 / v**2
    darg_dvar = darg_dz * dz_dvar

    L0 = -np.log(gamma((v+1)/2))
    dL0_dv = -polygamma(0, (v+1)/2) * 1/2

    L1 = 1/2 * np.log(np.pi)

    L2 = 1/2 * np.log(v)
    dL2_dv = 1 / (2*v)

    L3 = np.log(gamma(v/2))
    dL3_dv = polygamma(0, v/2) * 1/2

    L4 = 1/2 * np.log(var)
    dL4_dvar = 1/2 * 1/var

    L5 = (v+1)/2 * np.log(arg)
    dL5_dvar = (v+1)/2 * 1/arg * darg_dvar
    dL5_dv = 1/2 * np.log(arg) + (v+1)/2 * 1/arg * darg_dv

    L = L0 + L1 + L2 + L3 + L4 + L5
    dL_dvar = dL4_dvar + dL5_dvar
    dL_dv = dL0_dv + dL2_dv + dL3_dv + dL5_dv
    return flex.double(L), flex.double(dL_dvar), flex.double(dL_dvar), flex.double(dL_dv)

  def _loss_function_original(self, delta_sq, var):
    # This just sums the delta squares
    L = delta_sq / var
    dL_dvar = -delta_sq / var**2
    return L, dL_dvar

  def _get_sadd(self, correlation):
    sadd = flex.double(len(correlation), 0)
    dsadd_dsaddi = [flex.double(len(correlation), 0) for i in range(self.n_coefs)]
    for degree_index in range(self.n_coefs):
      sadd += self.sadd[degree_index] * correlation**degree_index
      dsadd_dsaddi[degree_index] = correlation**degree_index
    return sadd, dsadd_dsaddi

  def _get_var_ev11(self, counting_err, biased_mean, correlation, return_der=False):
    sadd, dsadd_dsaddi = self._get_sadd(correlation)
    var = self.sfac**2 * (counting_err + self.sb**2 * biased_mean + sadd**2 * biased_mean**2)
    if return_der:
      dvar_dsfac = 2 * self.sfac * (counting_err + self.sb**2 * biased_mean + sadd**2 * biased_mean**2)
      dvar_dsb = self.sfac**2 * 2 * self.sb * biased_mean
      dvar_dsadd = self.sfac**2 * 2 * sadd * biased_mean**2
      return var, dvar_dsfac, dvar_dsb, dvar_dsadd, dsadd_dsaddi
    else:
      return var

  def calculate_functional_mll(self):
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    L        = 0
    dL_dsfac = 0
    dL_dsb   = 0
    dL_dsadd = [0 for i in range(self.n_coefs)]
    if self.tuning_param_opt:
      dL_dnu   = 0

    for bin_index, differences in enumerate(self.pairwise_differences):
      number_of_differences_in_bin = differences.size()
      global_number_of_differences_in_bin = comm.reduce(number_of_differences_in_bin, MPI.SUM, root=0)
      L_in_bin = flex.double(number_of_differences_in_bin, 0)
      dL_dsfac_in_bin = flex.double(number_of_differences_in_bin, 0)
      dL_dsb_in_bin = flex.double(number_of_differences_in_bin, 0)
      dL_dsadd_in_bin = [flex.double(number_of_differences_in_bin, 0) for i in range(self.n_coefs)]
      if self.tuning_param_opt:
        dL_dnu_in_bin = flex.double(number_of_differences_in_bin, 0)

      if number_of_differences_in_bin > 0:
        var_i, dvar_i_dsfac, dvar_i_dsb, dvar_i_dsadd, dsadd_i_dsaddi = self._get_var_ev11(
          differences['var_cs_i'],
          differences['biased_mean'],
          differences['correlation_i'],
          return_der=True
          )
        var_j, dvar_j_dsfac, dvar_j_dsb, dvar_j_dsadd, dsadd_j_dsaddi = self._get_var_ev11(
          differences['var_cs_j'],
          differences['biased_mean'],
          differences['correlation_j'],
          return_der=True
          )

        if self.likelihood == 'normal':
          L_in_bin, dL_dvar_i, dL_dvar_j = self._loss_function_gaus(
            differences['pairwise_differences'], var_i, var_j
            )
        elif self.likelihood == 't-dist':
          if self.tuning_param_opt:
            L_in_bin, dL_dvar_i, dL_dvar_j, dL_dnu_in_bin = self._loss_function_t_opt(
              differences['pairwise_differences'], var_i, var_j
              )
          else:
            L_in_bin, dL_dvar_i, dL_dvar_j = self._loss_function_t(
              differences['pairwise_differences'], var_i, var_j
              )

        dL_dsfac_in_bin = dL_dvar_i * dvar_i_dsfac + dL_dvar_j * dvar_j_dsfac
        dL_dsb_in_bin = dL_dvar_i * dvar_i_dsb + dL_dvar_j * dvar_j_dsb
        for degree_index in range(self.n_coefs):
          dL_dsadd_in_bin[degree_index] = \
            dL_dvar_i * dvar_i_dsadd * dsadd_i_dsaddi[degree_index] \
            + dL_dvar_j * dvar_j_dsadd * dsadd_j_dsaddi[degree_index]

      global_sum_of_L_in_bin = comm.reduce(flex.sum(L_in_bin), MPI.SUM, root=0)
      global_sum_of_dL_dsfac_in_bin = comm.reduce(flex.sum(dL_dsfac_in_bin), MPI.SUM, root=0)
      global_sum_of_dL_dsb_in_bin = comm.reduce(flex.sum(dL_dsb_in_bin), MPI.SUM, root=0)
      global_sum_of_dL_dsadd_in_bin = [None for i in range(self.n_coefs)]
      for degree_index in range(self.n_coefs):
        global_sum_of_dL_dsadd_in_bin[degree_index] = \
          comm.reduce(flex.sum(dL_dsadd_in_bin[degree_index]), MPI.SUM, root=0)
      if self.tuning_param_opt:
        global_sum_of_dL_dnu_in_bin = comm.reduce(flex.sum(dL_dnu_in_bin), MPI.SUM, root=0)

      if self.mpi_helper.rank == 0 and global_number_of_differences_in_bin > 0:
        #global_weight_for_bin = 1 / math.sqrt(global_number_of_differences_in_bin)
        global_weight_for_bin = math.sqrt(self.n_refls_in_bin[bin_index]) / global_number_of_differences_in_bin
        L += global_weight_for_bin * global_sum_of_L_in_bin
        dL_dsfac += global_weight_for_bin * global_sum_of_dL_dsfac_in_bin
        dL_dsb += global_weight_for_bin * global_sum_of_dL_dsb_in_bin
        for degree_index in range(self.n_coefs):
          dL_dsadd[degree_index] += global_weight_for_bin * global_sum_of_dL_dsadd_in_bin[degree_index]
        if self.tuning_param_opt:
          dL_dnu += global_weight_for_bin * global_sum_of_dL_dnu_in_bin

    # Broadcast these derivates and functional values to all ranks
    self.functional = comm.bcast(L, root=0)
    self.der_wrt_sfac = comm.bcast(dL_dsfac, root=0)
    self.der_wrt_sb = comm.bcast(dL_dsb, root=0)
    self.der_wrt_sadd = comm.bcast(dL_dsadd, root=0)
    if self.tuning_param_opt:
      self.der_wrt_nu = comm.bcast(dL_dnu, root=0)

  def calculate_functional_ev11(self):
    # Results of calculation (on rank 0):
    TF        = 0
    dTF_dsfac = 0
    dTF_dsb   = 0
    dTF_dsadd = [0 for i in range(self.n_coefs)]
    if self.tuning_param_opt:
      dTF_dnu   = 0

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI

    # these reflections come from self.work_table
    # self.work_table['delta_sq'] = [(n-1)/n * (I_bk - <I'_hk>)]^2
    #   delta_bk^2 = self.work_table['delta_sq'] / var_ev11 (eq 13 from Brewster 2019)
    # sum_of_deviations_squared_in_bin = sum( (1 - delta_bk^2)^2 )
    for reflections in self.intensity_bins:
      number_of_reflections_in_bin = reflections.size()
      global_number_of_reflections_in_bin = comm.reduce(number_of_reflections_in_bin, MPI.SUM, root=0)
      TF_in_bin = flex.double(number_of_reflections_in_bin, 0)
      dTF_dsfac_in_bin = flex.double(number_of_reflections_in_bin, 0)
      dTF_dsb_in_bin = flex.double(number_of_reflections_in_bin, 0)
      dTF_dsadd_in_bin = [flex.double(number_of_reflections_in_bin, 0) for i in range(self.n_coefs)]
      if self.tuning_param_opt:
        dTF_dnu_in_bin = flex.double(number_of_reflections_in_bin, 0)

      if number_of_reflections_in_bin > 0:
        var, dvar_dsfac_in_bin, dvar_dsb_in_bin, dvar_dsadd_in_bin, dsadd_dsaddi = \
          self._get_var_ev11(
          reflections['intensity.sum.variance'],
          reflections['biased_mean'],
          reflections[self.cc_key],
          return_der=True
          )
        TF_in_bin, dL_dvar = self._loss_function_original(reflections['delta_sq'], var)
        dTF_dsfac_in_bin = dL_dvar * dvar_dsfac_in_bin
        dTF_dsb_in_bin = dL_dvar * dvar_dsb_in_bin
        for degree_index in range(self.n_coefs):
          dTF_dsadd_in_bin[degree_index] = dL_dvar * dvar_dsadd_in_bin * dsadd_dsaddi[degree_index]

      # TF_in_bin = delta_sq
      # dTF_d*_in_bin = ddelta_sq_d*
      global_sum_of_delta_sq_in_bin = comm.reduce(flex.sum(TF_in_bin), MPI.SUM, root=0)
      global_sum_of_ddelta_sq_dsfac_in_bin = comm.reduce(flex.sum(dTF_dsfac_in_bin), MPI.SUM, root=0)
      global_sum_of_ddelta_sq_dsb_in_bin = comm.reduce(flex.sum(dTF_dsb_in_bin), MPI.SUM, root=0)
      global_sum_of_ddelta_sq_dsadd_in_bin = [None for i in range(self.n_coefs)]
      for degree_index in range(self.n_coefs):
        global_sum_of_ddelta_sq_dsadd_in_bin[degree_index] = \
          comm.reduce(flex.sum(dTF_dsadd_in_bin[degree_index]), MPI.SUM, root=0)
      if self.mpi_helper.rank == 0 and global_number_of_reflections_in_bin > 0:
        global_weight_for_bin = math.sqrt(global_number_of_reflections_in_bin)

        std_delta_in_bin = (global_sum_of_delta_sq_in_bin / global_number_of_reflections_in_bin)**(1/2)
        prefactor = 1/2 \
          * (global_sum_of_delta_sq_in_bin / global_number_of_reflections_in_bin)**(-1/2) \
          / global_number_of_reflections_in_bin
        dstd_delta_in_bin_dsfac = prefactor * global_sum_of_ddelta_sq_dsfac_in_bin
        dstd_delta_in_bin_dsb = prefactor * global_sum_of_ddelta_sq_dsb_in_bin
        dstd_delta_in_bin_dsadd = [None for i in range(self.n_coefs)]
        for degree_index in range(self.n_coefs):
          dstd_delta_in_bin_dsadd[degree_index] = prefactor * global_sum_of_ddelta_sq_dsadd_in_bin[degree_index]

        TF += global_weight_for_bin * (1 - std_delta_in_bin)**2
        prefactor = -2 * global_weight_for_bin * (1 - std_delta_in_bin)
        dTF_dsfac += prefactor * dstd_delta_in_bin_dsfac
        dTF_dsb += prefactor * dstd_delta_in_bin_dsb
        for degree_index in range(self.n_coefs):
          dTF_dsadd[degree_index] += prefactor * dstd_delta_in_bin_dsadd[degree_index]

    # Broadcast these derivates and functional values to all ranks
    self.functional = comm.bcast(TF, root=0)
    self.der_wrt_sfac = comm.bcast(dTF_dsfac, root=0)
    self.der_wrt_sb = comm.bcast(dTF_dsb, root=0)
    self.der_wrt_sadd = comm.bcast(dTF_dsadd, root=0)
    if self.tuning_param_opt:
      self.der_wrt_nu = comm.bcast(dTF_dnu, root=0)

  def modify_errors(self, reflections):
    # First set up a reflection table to do work downstream.
    # Needed to calculate delta_sq and bin reflections
    reflections = self.setup_work_arrays(reflections)
    # Make sure every rank knows the global mean/stdev for deltas and use them to get the bin limits
    self.calculate_delta_statistics()
    self.calculate_delta_bin_limits()
    # assign deltas for each reflection to appropriate bin
    self.distribute_deltas_over_bins()
    # Each rank gets its own bin. Make sure all deltas in that bin are on that rank and sorted.
    self.distribute_deltas_over_ranks()
    # calculate rankits, each rank does its own rankits calculation
    self.calculate_delta_rankits()
    #   ev11 params using slope and offset of fit to rankits
    self.calculate_initial_ev11_parameters()
    # Now moving to intensities, find the bin limits using global min/max of the means of each reflection
    self.calculate_intensity_bin_limits()
    # Once bin limits are determined, assign intensities on each rank to appropriate bin limits
    self.distribute_reflections_over_intensity_bins()
    if self.algorithm == 'ev11_mll':
      self.setup_pairwise_differences()
    # Run LBFGSB minimizer
    #  -- only rank0 does minimization but gradients/functionals are calculated using all rank
    self.run_minimizer()
    # Finally update the variances of each reflection as per Eq (10) in Brewster et. al (2019)
    reflections['intensity.sum.variance'] = self._get_var_ev11(
      reflections['intensity.sum.variance'],
      reflections['biased_mean'],
      reflections[self.cc_key]
      )
    del reflections['biased_mean']
    return reflections


if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(error_modifier)
