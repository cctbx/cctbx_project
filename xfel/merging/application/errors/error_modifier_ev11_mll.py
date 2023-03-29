"""
ToDo:
  add use_groups to merging params
  self.setup_work_arrays
    use a median as an estimator for the biased_mean and mean_without_value
How does this code work
* Indicates modifications
  1: run
    calls self.modify_errors
  2: self.modify_errors
    does the steps
    *1: self.setup_work_arrays(reflections)
      Modified to determine if groups of ev11 will be used and calculate the total number of groups

      Create work_table and add biased_mean to reflections
      self.deltas                      <- eq11 deltas
      self.work_table["delta_sq"]    = delta_sq <- eq13 delta without dividing by var
      self.work_table["mean"]        = mean <- mean of Ihkl not including itself
      self.work_table["biased_mean"] = biased_mean <- mean intensity of I_hkl (including itself)
      self.work_table["var"]         = var <- intensity.sum.variance
      reflections['biased_mean'] = all_biased_mean

    *2: self.calculate_delta_statistics()
      count, min, max, mean, std of all deltas (eq11)
      edited a really long line
    *3: self.calculate_delta_bin_limits()
      Removed redundant numpy import
    4: self.distribute_deltas_over_bins()
    5: self.distribute_deltas_over_ranks()
    *6: self.calculate_delta_rankits()
      does function title
      Removed redundant numpy import

    *7: self.calculate_initial_ev11_parameters()
      These should be the same for all groups.
      Make these parameters lists of

    *8: self.calculate_intensity_bin_limits()
      Moved so it was ordered in the file in the order it is called by self.modify_errors
      calculates the min / max of the biased means and makes bins

    9: self.distribute_reflections_over_intensity_bins()
      does as says

    *10: self.run_minimizer()
      extensive modified for additional params

    *11: modifying the errors in self.modify_errors

*get_overall_correlation_flex
  removed redundant math import

"""
from __future__ import absolute_import, division, print_function
from dials.array_family import flex
from xfel.merging.application.worker import worker
from xfel.merging.application.reflection_table_utils import reflection_table_utils
import math
import numpy as np
import os
import sys


number_of_intensity_bins = 100

class error_modifier_ev11_mll(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(error_modifier_ev11_mll, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)
    self.n_degrees = self.params.merging.error.ev11_mll.n_degrees
    self.tuning_param = self.params.merging.error.ev11_mll.tuning_param

  def __repr__(self):
    return 'Adjust intensity errors -- ev11_mll'

  def run(self, experiments, reflections):
    '''Modify intensity errors according to EV11 -- Brewster2019'''
    assert self.params.merging.error.model == "ev11_mll"
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

    for refls in reflection_table_utils.get_next_hkl_reflection_table(reflections):
      number_of_measurements = refls.size()
      if number_of_measurements == 0: # if the returned "refls" list is empty, it's the end of the input "reflections" list
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
        self.deltas.extend(delta/flex.sqrt(refls['intensity.sum.variance'])) # Please be careful about where to put the var
        delta_sq.extend(delta**2)
        mean.extend(mean_without_val)
        biased_mean.extend(refls_biased_mean)
        var.extend(refls['intensity.sum.variance'])
        correlation.extend(refls['correlation'])

    self.work_table["correlation"] = correlation
    self.work_table["delta_sq"]    = delta_sq
    self.work_table["delta_all"]   = delta_all
    self.work_table["mean"]        = mean
    self.work_table["biased_mean"] = biased_mean
    self.work_table["var"]         = var
    self.work_table["I"]           = I_wt
    reflections['biased_mean']     = all_biased_mean

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
    '''Divide the delta (min,max) range into "number of ranks" bins. For a balanced rank load, bin limits should be
       chosen so that the bins are equally populated by the deltas. Assuming the normal distribution of deltas,
       we use the probability density function for the bin calculations.'''
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

  def get_overall_correlation_flex(self, data_a, data_b):
    # This function does not get called within error_modifier
    """
    Correlate any two sets of data.
    @param data_a - references
    @param data_b - observations
    @return tuple containing correlation coefficent, slope and offset.
    """
    assert len(data_a) == len(data_b)
    corr = 0
    slope = 0
    offset = 0
    try:
      sum_xx = 0
      sum_xy = 0
      sum_yy = 0
      sum_x  = 0
      sum_y  = 0
      N      = 0
      for i in range(len(data_a)):
        I_r       = data_a[i]
        I_o       = data_b[i]
        N      += 1
        sum_xx += I_r**2
        sum_yy += I_o**2
        sum_xy += I_r * I_o
        sum_x  += I_r
        sum_y  += I_o
      slope = (N * sum_xy - sum_x * sum_y) / (N * sum_xx - sum_x**2)
      offset = (sum_xx * sum_y - sum_x * sum_xy) / (N * sum_xx - sum_x**2)
      corr  = (N * sum_xy - sum_x * sum_y) / (math.sqrt(N * sum_xx - sum_x**2) *
                 math.sqrt(N * sum_yy - sum_y**2))
    except ZeroDivisionError:
      pass

    return corr, slope, offset

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
      self.sb = math.sqrt(offset)
      self.sadd = [0 for i in range(self.n_degrees)]
      self.sadd[0] = offset

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
    '''Calculate the minimum and maximum values of the mean intensities for each HKL'''
    count = self.work_table.size()
    mean_intensity_min = flex.min(self.work_table['biased_mean']) if count > 0 else float('inf')
    mean_intensity_max = flex.max(self.work_table['biased_mean']) if count > 0 else float('-inf')
    if count > 0:
      self.logger.log(
        f"Using {count} multiply-measured HKLs; "
        + f"mean intensity (min,max): ({mean_intensity_min},{mean_intensity_max})"
      )
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

  def run_minimizer(self):
    from scitbx import lbfgsb

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    size = self.mpi_helper.size
    self.n = self.n_degrees + 2
    # The * combines the grouped params into one group
    self.x = flex.double([self.sfac, self.sb, *self.sadd])
    if self.mpi_helper.rank == 0:
      self.logger.main_log(
        'Initial Parameter Estimates = '
        + f'sfac: {self.sfac} '
        + f'sb: {self.sb} '
        + f'sadd: {self.sadd[0]}'
      )
    l = flex.double(self.n, 1e-8)
    for degree_index in range(self.n_degrees):
      l[degree_index + 2] = -10
    if self.mpi_helper.rank == 0:
      self.minimizer = lbfgsb.minimizer(
        n = self.n,
        l = l,
        u = flex.double(self.n, 0),
        nbd = flex.int(self.n, 1),
      )
    delta_print_index = 0
    while True:
      self.compute_functional_and_gradients()
      #self.get_deltas(delta_print_index)
      delta_print_index += 1
      status=-1
      if self.mpi_helper.rank == 0:
        if self.minimizer.process(self.x, self.f, self.g):
          self.sfac = self.x[0]
          self.sb   = self.x[1]
          self.sadd = self.x[2:]
          sfac = f'{self.sfac:0.3f}'
          sb = f'{self.sb:0.6f}'
          sadd = [f'{self.sadd[i]:0.3f}' for i in range(self.n_degrees)]
          self.logger.main_log(
            'intermediate minimization results = '
            + f'functional: {self.f:.2f} '
            + f'sfac: {sfac} '
            + f'sb: {sb} '
            + f'sadd: {sadd}'
          )
          status=1
        elif self.minimizer.is_terminated():
          status=0

      comm.barrier()
      status=comm.bcast(status, root=0)
      if status==1:
        self.sfac=comm.bcast(self.sfac, root=0)
        self.sb=comm.bcast(self.sb, root=0)
        self.sadd=comm.bcast(self.sadd, root=0)
        pass
      if status==0:
        break

    if self.mpi_helper.rank == 0:
      sfac = f'{self.sfac:0.3f}'
      sb = f'{self.sb:0.3f}'
      sadd = [f'{self.sadd[i]:0.3f}' for i in range(self.n_degrees)]
      self.logger.main_log(
        'FINAL EV11 VALUES = '
        + f'functional: {self.f:.2f} '
        + f'sfac: {sfac} '
        + f'sb: {sb} '
        + f'sadd: {sadd}'
      )
      self.plot_sadd()
    if self.params.merging.error.ev11_mll.verify_derivatives:
      self.verify_derivatives()

  def plot_sadd(self):
    import matplotlib.pyplot as plt
    import os
    c = flex.double(np.linspace(0, 1, 101))
    sadd, _ = self._get_sadd(c)
    fig, axes = plt.subplots(1, 1, figsize=(6, 4))
    axes.plot(c, sadd)
    axes.set_xlabel('correlation')
    axes.set_ylabel('$s_{add}$')
    fig.tight_layout()
    fig.savefig(os.path.join(
      self.params.output.output_dir,
      self.params.output.prefix + '_sadd.png'
      ))
    plt.close()
    return None

  def compute_functional_and_gradients(self):
    self.calculate_functional_ev11()
    self.f = self.functional
    self.g = flex.double([self.der_wrt_sfac, self.der_wrt_sb, *self.der_wrt_sadd])
    return self.f, self.g

  def verify_derivatives(self):
    shift = 0.0001
    import copy
    self.n = self.n_degrees + 2
    self.calculate_functional_ev11()
    TF = copy.copy(self.functional)
    der_wrt_sfac = copy.copy(self.der_wrt_sfac)
    der_wrt_sb = copy.copy(self.der_wrt_sb)
    der_wrt_sadd = copy.copy(self.der_wrt_sadd)
    sfac = copy.copy(self.sfac)
    sb = copy.copy(self.sb)
    sadd = copy.copy(self.sadd)

    # sfac
    self.sfac = sfac * (1 + shift)
    self.calculate_functional_ev11()
    TF_p = copy.copy(self.functional)
    self.sfac = sfac * (1 - shift)
    self.calculate_functional_ev11()
    TF_m = copy.copy(self.functional)
    check_der_wrt_sfac = (TF_p - TF_m) / (2 * shift * sfac)
    self.sfac = sfac
    if self.mpi_helper.rank == 0:
      print(f'der_wrt_sfac numerical: {check_der_wrt_sfac} analytical {der_wrt_sfac}')

    # sb
    self.sb = sb * (1 + shift)
    self.calculate_functional_ev11()
    TF_p = copy.copy(self.functional)
    self.sb = sb * (1 - shift)
    self.calculate_functional_ev11()
    TF_m = copy.copy(self.functional)
    check_der_wrt_sb = (TF_p - TF_m) / (2 * shift * sb)
    self.sb = sb
    if self.mpi_helper.rank == 0:
      print(f'der_wrt_sb numerical: {check_der_wrt_sb} analytical {der_wrt_sb}')

    # sadd:
    for degree_index in range(self.n_degrees):
      if sadd[degree_index] == 0:
        self.sadd[degree_index] = shift
      else:
        self.sadd[degree_index] = sadd[degree_index] * (1 + shift)
      self.calculate_functional_ev11()
      TF_p = copy.copy(self.functional)
      if sadd[degree_index] == 0:
        self.sadd[degree_index] = -shift
      else:
        self.sadd[degree_index] = sadd[degree_index] * (1 - shift)
      self.calculate_functional_ev11()
      TF_m = copy.copy(self.functional)
      if sadd[degree_index] == 0:
        check_der_wrt_sadd = (TF_p - TF_m) / (2 * shift)
      else:
        check_der_wrt_sadd = (TF_p - TF_m) / (2 * shift * sadd[degree_index])
      self.sadd[degree_index] = sadd[degree_index]
      if self.mpi_helper.rank == 0:
        print(f'der_wrt_sadd - degree {degree_index} numerical: {check_der_wrt_sadd} analytical {der_wrt_sadd[degree_index]}')
    return None

  def get_deltas(self, delta_print_index):
    deltas = flex.double()
    deltas_no_ev11 = flex.double()
    deltas_outlier = flex.double()
    I = flex.double()
    for reflections in self.intensity_bins:
      number_of_reflections_in_bin = reflections.size()
      if number_of_reflections_in_bin > 0:
        sfac = self.sfac
        sb = self.sb
        v = reflections['var']
        mi = reflections['biased_mean']
        sadd, _ = self._get_sadd(reflections['correlation'])
        var_ev11 = sfac**2 * (v + sb**2 * mi + sadd**2 * mi**2)

        deltas.extend(reflections['delta_all'] / flex.sqrt(var_ev11))
        deltas_no_ev11.extend(reflections['delta_all'] / flex.sqrt(v))
        deltas_outlier.extend(reflections['delta_all'] / flex.sqrt(mi))
        I.extend(reflections['I'])

    output_data_rank = np.column_stack((
        deltas.as_numpy_array(),
        I.as_numpy_array(),
        deltas_no_ev11.as_numpy_array(),
        deltas_outlier.as_numpy_array(),
      ))

    output_data_all = self.mpi_helper.comm.gather(output_data_rank, root=0)
    if self.mpi_helper.rank == 0:
      for rank_index in range(self.mpi_helper.comm.Get_size()):
        if rank_index == 0:
          output_data_to_save = output_data_all[rank_index]
        else:
          output_data_to_save = np.concatenate((output_data_to_save, output_data_all[rank_index]))
      save_to_name = os.path.join(
        self.params.output.output_dir,
        self.params.output.prefix + f'_deltas_{delta_print_index:03d}.npy'
        )
      np.save(save_to_name, output_data_to_save)
    return None

  def _loss_function_gaus(self, delta_sq, var):
    s = self.tuning_param
    delta_sq = delta_sq.as_numpy_array()
    var = var.as_numpy_array()
    L = 1 - np.exp(-1/2 * (1 - delta_sq / var)**2 / s**2)
    dL_dvar = (1 - delta_sq / var) / s**2 * delta_sq / var**2 * np.exp(-1/2 * (1 - delta_sq / var)**2 / s**2)
    L = flex.double(L)
    dL_dvar = flex.double(dL_dvar)
    return L, dL_dvar

  def _loss_function_t(self, delta_sq, var):
    v = self.tuning_param
    delta_sq = delta_sq.as_numpy_array()
    var = var.as_numpy_array()
    arg = 1 + (1 - delta_sq / var)**2 / v
    L = (v+1)/2 * np.log(arg)
    dL_dvar = (v+1)/2 * 1/arg * 2*(1-delta_sq/var)/v * delta_sq/var**2
    L = flex.double(L)
    dL_dvar = flex.double(dL_dvar)
    return L, dL_dvar

  def _get_sadd(self, correlation):
    sadd = flex.double(len(correlation), 0)
    dsadd_dsaddi = [flex.double(len(correlation), 0) for i in range(self.n_degrees)]
    for degree_index in range(self.n_degrees):
      sadd += self.sadd[degree_index] * correlation**degree_index
      dsadd_dsaddi[degree_index] = correlation**degree_index
    return sadd, dsadd_dsaddi

  def calculate_functional_ev11(self):
    # Results of calculation (on rank 0):
    TF        = 0
    dTF_dsfac = 0
    dTF_dsb   = 0
    dTF_dsadd = [0 for i in range(self.n_degrees)]

    if self.mpi_helper.rank == 0:
      multiplier = 0
      for degree_index in range(2, self.n_degrees):
        TF += (multiplier * degree_index * self.sadd[degree_index])**2
        dTF_dsadd[degree_index] += 2 * multiplier**2 * degree_index**2 * self.sadd[degree_index]

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
      dTF_dsadd_in_bin = [flex.double(number_of_reflections_in_bin, 0) for i in range(self.n_degrees)]

      if number_of_reflections_in_bin > 0:
        sfac = self.sfac
        sb = self.sb
        v = reflections['var']
        mi = reflections['biased_mean']
        sadd, dsadd_dsaddi = self._get_sadd(reflections['correlation'])
        var = sfac**2 * (v + sb**2 * mi + sadd**2 * mi**2)

        #L, dL_dvar = self._loss_function_gaus(reflections['delta_sq'], var)
        TF_in_bin, dL_dvar = self._loss_function_t(reflections['delta_sq'], var)

        dvar_dsfac_in_bin = 2 * sfac * (v + sb**2 * mi + sadd**2 * mi**2)
        dvar_dsb_in_bin = sfac**2 * 2 * sb * mi
        dvar_dsadd_in_bin = sfac**2 * 2 * sadd * mi**2

        dTF_dsfac_in_bin = dL_dvar * dvar_dsfac_in_bin
        dTF_dsb_in_bin = dL_dvar * dvar_dsb_in_bin
        for degree_index in range(self.n_degrees):
          dTF_dsadd_in_bin[degree_index] = dL_dvar * dvar_dsadd_in_bin * dsadd_dsaddi[degree_index]

      global_sum_of_TF_in_bin = comm.reduce(flex.sum(TF_in_bin),  MPI.SUM, root=0)
      global_sum_of_dTF_dsfac_in_bin = comm.reduce(flex.sum(dTF_dsfac_in_bin), MPI.SUM, root=0)
      global_sum_of_dTF_dsb_in_bin = comm.reduce(flex.sum(dTF_dsb_in_bin), MPI.SUM, root=0)
      global_sum_of_dTF_dsadd_in_bin = [None for i in range(self.n_degrees)]
      for degree_index in range(self.n_degrees):
        global_sum_of_dTF_dsadd_in_bin[degree_index] = comm.reduce(flex.sum(dTF_dsadd_in_bin[degree_index]), MPI.SUM, root=0)

      if self.mpi_helper.rank == 0 and global_number_of_reflections_in_bin > 0:
        global_weight_for_bin = 1 / math.sqrt(global_number_of_reflections_in_bin)
        TF += global_weight_for_bin * global_sum_of_TF_in_bin
        dTF_dsfac += global_weight_for_bin * global_sum_of_dTF_dsfac_in_bin
        dTF_dsb += global_weight_for_bin * global_sum_of_dTF_dsb_in_bin
        for degree_index in range(self.n_degrees):
          dTF_dsadd[degree_index] += global_weight_for_bin * global_sum_of_dTF_dsadd_in_bin[degree_index]
    # Broadcast these derivates and functional values to all ranks
    self.functional = comm.bcast(TF, root=0)
    self.der_wrt_sfac = comm.bcast(dTF_dsfac, root=0)
    self.der_wrt_sb = comm.bcast(dTF_dsb, root=0)
    self.der_wrt_sadd = comm.bcast(dTF_dsadd, root=0)

  def modify_errors(self, reflections):
    # First set up a reflection table to do work downstream. Needed to calculate delta_sq and bin reflections
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
    # Run LBFGSB minimizer -- only rank0 does minimization but gradients/functionals are calculated using all rank
    self.run_minimizer()
    # Finally update the variances of each reflection as per Eq (10) in Brewster et. al (2019)

    sadd, _ = self._get_sadd(reflections['correlation'])
    reflections['intensity.sum.variance'] = \
      self.sfac**2 * ( \
        reflections['intensity.sum.variance'] \
        + self.sb**2 * reflections['biased_mean']
        + sadd**2 * reflections['biased_mean']**2
      )

    del reflections['biased_mean']

    return reflections


if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(error_modifier)
