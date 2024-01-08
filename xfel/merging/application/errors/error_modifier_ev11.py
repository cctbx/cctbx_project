from __future__ import absolute_import, division, print_function
from dials.array_family import flex
from xfel.merging.application.worker import worker
from xfel.merging.application.reflection_table_utils import reflection_table_utils
import math
import numpy as np
import sys

number_of_intensity_bins = 100

class error_modifier_ev11(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(error_modifier_ev11, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Adjust intensity errors -- ev11'

  def run(self, experiments, reflections):
    '''Modify intensity errors according to EV11 -- Brewster2019'''
    assert self.params.merging.error.model == "ev11"

    self.logger.log_step_time("ERROR_MODIFIER_EV11")
    self.logger.log("Modifying intensity errors -- ev11 method (starting with %d reflections)"%(len(reflections)))
    reflections = self.modify_errors(reflections)
    self.logger.log_step_time("ERROR_MODIFIER_EV11", True)

    return experiments, reflections

  def calculate_intensity_bin_limits(self):
    '''Calculate the minimum and maximum values of the mean intensities for each HKL'''
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

  def setup_work_arrays(self, reflections):
    '''Select multiply-measured HKLs. Calculate and cache reflection deltas, deltas squared, and HKL means for every reflection'''
    self.deltas     = flex.double()
    self.work_table = flex.reflection_table()
    delta_sq        = flex.double()
    mean            = flex.double() # mean = <I'_hj>
    biased_mean     = flex.double() # biased_mean = <I_h>, so dont leave out any reflection
    var             = flex.double()
    all_biased_mean = flex.double()

    for refls in reflection_table_utils.get_next_hkl_reflection_table(reflections):
      number_of_measurements = refls.size()
      if number_of_measurements == 0: # if the returned "refls" list is empty, it's the end of the input "reflections" list
        break
      refls_biased_mean = flex.double(len(refls), flex.mean(refls['intensity.sum.value']))
      all_biased_mean.extend(refls_biased_mean)

      if number_of_measurements > self.params.merging.minimum_multiplicity:
        nn_factor_sqrt = math.sqrt((number_of_measurements - 1) / number_of_measurements)
        i_sum = flex.double(number_of_measurements, flex.sum(refls['intensity.sum.value']))
        i_sum_minus_val = i_sum - refls['intensity.sum.value']
        mean_without_val = i_sum_minus_val/(number_of_measurements-1)
        delta = nn_factor_sqrt * (refls['intensity.sum.value'] - mean_without_val)
        self.deltas.extend(delta/flex.sqrt(refls['intensity.sum.variance'])) # Please be careful about where to put the var
        delta_sq.extend(delta**2)
        mean.extend(mean_without_val)
        biased_mean.extend(refls_biased_mean)
        var.extend(refls['intensity.sum.variance'])

    self.work_table["delta_sq"]    = delta_sq
    self.work_table["mean"]        = mean
    self.work_table["biased_mean"] = biased_mean
    self.work_table["var"]         = var
    reflections['biased_mean'] = all_biased_mean
    self.logger.log("Number of work reflections selected: %d"%self.deltas.size())
    return reflections

  def calculate_functional_ev11(self):
    # Results of calculation (on rank 0):
    func           = 0
    der_wrt_sfac   = 0
    der_wrt_sb     = 0
    der_wrt_sadd   = 0

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI

    for reflections in self.intensity_bins:
      number_of_reflections_in_bin = reflections.size()

      sum_of_delta_squared_in_bin = flex.double(number_of_reflections_in_bin, 0)
      sum_of_der_wrt_sfac_in_bin  = flex.double(number_of_reflections_in_bin, 0)
      sum_of_der_wrt_sb_in_bin    = flex.double(number_of_reflections_in_bin, 0)
      sum_of_der_wrt_sadd_in_bin  = flex.double(number_of_reflections_in_bin, 0)

      if number_of_reflections_in_bin > 0:
        variance = reflections['var']
        mean_intensity = reflections['biased_mean'] # the mean intensity of the sample of HKLs which includes this reflection

        var_ev11               = self.sfac**2*(variance + self.sb**2*mean_intensity + self.sadd**2*mean_intensity**2)
        var_ev11_der_over_sfac = 2*self.sfac * (variance + self.sb**2 * mean_intensity + self.sadd**2 * mean_intensity**2)
        var_ev11_der_over_sb   = self.sfac**2 * 2*self.sb   * mean_intensity
        var_ev11_der_over_sadd = self.sfac**2 * 2*self.sadd * mean_intensity**2

        sum_of_delta_squared_in_bin += reflections['delta_sq'] / var_ev11

        sum_of_der_wrt_sfac_in_bin  -= reflections['delta_sq'] / var_ev11**2 * var_ev11_der_over_sfac
        sum_of_der_wrt_sb_in_bin    -= reflections['delta_sq'] / var_ev11**2 * var_ev11_der_over_sb
        sum_of_der_wrt_sadd_in_bin  -= reflections['delta_sq'] / var_ev11**2 * var_ev11_der_over_sadd

      global_number_of_reflections_in_bin   = comm.reduce(number_of_reflections_in_bin, MPI.SUM, root=0)
      global_sum_of_delta_squared_in_bin    = comm.reduce(flex.sum(sum_of_delta_squared_in_bin),  MPI.SUM, root=0)
      global_sum_of_der_wrt_sfac_in_bin     = comm.reduce(flex.sum(sum_of_der_wrt_sfac_in_bin),   MPI.SUM, root=0)
      global_sum_of_der_wrt_sb_in_bin       = comm.reduce(flex.sum(sum_of_der_wrt_sb_in_bin),     MPI.SUM, root=0)
      global_sum_of_der_wrt_sadd_in_bin     = comm.reduce(flex.sum(sum_of_der_wrt_sadd_in_bin),   MPI.SUM, root=0)

      if self.mpi_helper.rank == 0:
        if global_number_of_reflections_in_bin > 0 and global_sum_of_delta_squared_in_bin > 0:

          global_weight_for_bin = math.sqrt(global_number_of_reflections_in_bin)

          func += global_weight_for_bin * (1 - math.sqrt(global_sum_of_delta_squared_in_bin/global_number_of_reflections_in_bin))**2

          #if global_sum_of_delta_squared_in_bin == 0:
          #  from IPython import embed;embed()

          der_temp = global_weight_for_bin * (1 / math.sqrt(global_sum_of_delta_squared_in_bin/global_number_of_reflections_in_bin) - 1)  / global_number_of_reflections_in_bin

          der_wrt_sfac  -= der_temp * global_sum_of_der_wrt_sfac_in_bin
          der_wrt_sb    -= der_temp * global_sum_of_der_wrt_sb_in_bin
          der_wrt_sadd  -= der_temp * global_sum_of_der_wrt_sadd_in_bin


    #if self.mpi_helper.rank==0:
    #  self.functional = functional
    #  self.der_wrt_sfac = der_wrt_sfac
    #  self.der_wrt_sb = der_wrt_sb
    #  self.der_wrt_sadd = der_wrt_sadd

    # Broadcast these derivates and functional values to all ranks
    self.functional = comm.bcast(func, root=0)
    self.der_wrt_sfac = comm.bcast(der_wrt_sfac, root=0)
    self.der_wrt_sb = comm.bcast(der_wrt_sb, root=0)
    self.der_wrt_sadd = comm.bcast(der_wrt_sadd, root=0)

    #return (functional, der_wrt_sfac, der_wrt_sb, der_wrt_sadd)

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
      self.logger.main_log("Global delta statistics (count,min,max,mean,stddev): (%d,%f,%f,%f,%f)"%(self.global_delta_count, self.global_delta_min, self.global_delta_max, self.global_delta_mean, self.global_delta_stddev))

  def calculate_delta_bin_limits(self):
    '''Divide the delta (min,max) range into "number of ranks" bins. For a balanced rank load, bin limits should be
       chosen so that the bins are equally populated by the deltas. Assuming the normal distribution of deltas,
       we use the probability density function for the bin calculations.'''
    from scipy.stats import norm
    import numpy as np
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
    # Get the base global index for this rank's deltas. Example: if rank 0 has 10 deltas, the first delta on rank 1 will be the 10th global delta.
    delta_count_per_rank = self.mpi_helper.comm.allreduce([self.deltas.size()])
    base_delta_index = sum(delta_count_per_rank[0:self.mpi_helper.rank])
    self.logger.log("Delta base index: %d"%base_delta_index)

    from scitbx.math import distributions
    import numpy as np
    norm = distributions.normal_distribution()

    a = 3./8. if self.global_delta_count < 10. else 0.5

    self.rankits = flex.double()
    for i in range(self.deltas.size()):
      global_delta_index = base_delta_index + i
      rankit = norm.quantile((global_delta_index+1-a)/(self.global_delta_count+1-(2*a)))
      self.rankits.append(rankit)

  def get_overall_correlation_flex(self, data_a, data_b) :
    """
    Correlate any two sets of data.
    @param data_a - references
    @param data_b - observations
    @return tuple containing correlation coefficent, slope and offset.
    """
    import math

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
      self.sadd = offset
      self.sb = math.sqrt(self.sadd) if self.sadd > 0 else 0

      '''
      if True:
        from matplotlib import pyplot as plt
        import numpy as np
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
    self.sfac  = initial_params[0]
    self.sadd  = initial_params[1]
    self.sb    = initial_params[2]

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
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    size = self.mpi_helper.size
    self.n = 3
    self.x = flex.double([self.sfac, self.sb, self.sadd])
    self.logger.main_log('Initial Parameter Estimates = sdfac: %.2f  sdb: %.2f  sdadd: %.2f'%(self.sfac, self.sb, self.sadd))
    if True:
      from scitbx import lbfgsb
      l = flex.double(self.n, 1e-8)

      if len(l) > 3:
        for p in range(7,len(l)):
          l[p] = 1e-15 # g*

      if self.mpi_helper.rank == 0:
        self.minimizer = lbfgsb.minimizer(
          n = self.n,
          l = l,
          u = flex.double(self.n, 0),
          nbd = flex.int(self.n, 1),
        )
      while True:
        self.compute_functional_and_gradients()
        status=-1
        if self.mpi_helper.rank == 0:
          if self.minimizer.process(self.x, self.f, self.g):
            self.logger.main_log('intermediate minimization results = functional: %.2f  sdfac: %.2f sdb: %.2f sdadd: %.2f' %(self.f,self.x[0],self.x[1], self.x[2]))
            status=1
            self.sfac = self.x[0]
            self.sb = self.x[1]
            self.sadd = self.x[2]
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
      self.logger.main_log('FINAL SDFAC VALUES = functional: %.2f  sdfac: %.2f sdb: %.2f sdadd: %.2f' %(self.f,self.x[0],self.x[1], self.x[2]))

  def compute_functional_and_gradients(self):
    self.calculate_functional_ev11()
    self.f = self.functional
    self.g = flex.double([self.der_wrt_sfac, self.der_wrt_sb, self.der_wrt_sadd])
    return self.f, self.g

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
    # initial ev11 params using slope and offset of fit to rankits
    self.calculate_initial_ev11_parameters()
    # Now moving to intensities, find the bin limits using global min/max of the means of each reflection
    self.calculate_intensity_bin_limits()
    # Once bin limits are determined, assign intensities on each rank to appropriate bin limits
    self.distribute_reflections_over_intensity_bins()
    # Run LBFGSB minimizer -- only rank0 does minimization but gradients/functionals are calculated using all rank
    self.run_minimizer()
    # Finally update the variances of each reflection as per Eq (10) in Brewster et. al (2019)
    reflections['intensity.sum.variance'] = (self.sfac**2)*(reflections['intensity.sum.variance'] +
                                                            self.sb*self.sb*reflections['biased_mean'] +
                                                            self.sadd*self.sadd*reflections['biased_mean']**2)

    del reflections['biased_mean']

    return reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(error_modifier)
