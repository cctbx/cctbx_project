from __future__ import absolute_import, division, print_function
from six.moves import range
from xfel.merging.application.worker import worker
from dials.array_family import flex
import math
from libtbx import adopt_init_args
from libtbx.str_utils import format_value
from libtbx import table_utils
from xfel.merging.application.reflection_table_utils import reflection_table_utils
from six.moves import cStringIO as StringIO
from cctbx.crystal import symmetry
from cctbx import miller
from scitbx.math import basic_statistics

class intensity_table_bin(object):
  '''Storage class for parameters of a resolution bin used by the hkl intensity statistics table'''
  def __init__(self,
               i_bin=None,
               d_range=None,
               d_min=None,
               redundancy_asu=None,
               redundancy_obs=None,
               redundancy_to_edge=None,
               absent=None,
               complete_tag=None,
               completeness=None,
               measurements=None,
               multiply_measured_asu=None,
               predictions=None,
               mean_I=None,
               mean_I_sigI=None,
               sigmaa=None,
               unmerged_meanIsig=None,
               unmerged_stddevIsig=None,
               unmerged_skewIsig=None):
    adopt_init_args(self, locals())

class intensity_table(object):
  '''Represents a table of hkl intensity statistics for resolution bins'''
  def __init__(self):
    self.table = []
    self.cumulative_theor_asu_count = 0
    self.cumulative_observed_asu_count = 0
    self.cumulative_observed_count = 0
    self.cumulative_multiply_observed_asu_count = 0
    self.cumulative_I = 0.0
    self.cumulative_Isigma = 0.0

  def get_table_text(self):
    '''Produce formatted table text'''
    table_header = ["","","","","<asu","<obs","<pred","","","","Merged","Merged","Unmerged","Unmerged","Unmerged"]
    table_header2 = ["Bin","Resolution Range","Completeness","%","multi>","multi>","multi>",
                      "n_meas", "asu_m_meas", "n_pred","<I>","<I/sig(I)>", "mean(I/sig(I))", "stddev(I/sig(I))", "skew(I/sig(I))"]

    use_preds = False # TODO?

    include_columns = [True, True, True, True, True, True, use_preds, True, True, use_preds, True, True, True, True, True]

    table_data = []
    table_data.append(table_header)
    table_data.append(table_header2)

    for bin in self.table:
      table_row = []
      table_row.append("%3d" % bin.i_bin)
      table_row.append("%-13s" % bin.d_range)
      table_row.append("%13s" % bin.complete_tag)
      table_row.append("%5.2f" % (100*bin.completeness))
      table_row.append("%6.2f" % bin.redundancy_asu)
      table_row.append("%6.2f" % bin.redundancy_obs)
      table_row.append("%6.2f" % (0)) # if redundancy_to_edge is None else bin.redundancy_to_edge))
      table_row.append("%6d" % bin.measurements)
      table_row.append("%6d" % bin.multiply_measured_asu)
      table_row.append("%6d" % (0)) # if redundancy_to_edge is None else flex.sum(bin.predictions)))
      table_row.append("%8.0f" % bin.mean_I)
      table_row.append("%8.3f" % bin.mean_I_sigI)
      table_row.append("%8.3f" % bin.unmerged_meanIsig)
      table_row.append("%8.3f" % bin.unmerged_stddevIsig)
      table_row.append("%8.3f" % bin.unmerged_skewIsig)
      table_data.append(table_row)

    if len(table_data) <= 2:
      return ("Intensity statistics table could not be constructed -- no bins accepted.")

    table_data.append([""] * len(table_header))

    total_completeness = 0
    total_asu_multiplicity = 0
    if self.cumulative_theor_asu_count > 0:
      total_completeness = 100 * (self.cumulative_observed_asu_count / self.cumulative_theor_asu_count)
      total_asu_multiplicity = self.cumulative_observed_count / self.cumulative_theor_asu_count

    total_obs_multiplicity = 0
    if self.cumulative_observed_asu_count:
      total_obs_multiplicity = self.cumulative_observed_count / self.cumulative_observed_asu_count

    total_mean_I = 0
    total_mean_Isigma = 0
    total_unmerged_mean_Isigma = 0
    total_unmerged_stddev_Isigma = 0
    total_unmerged_skew_Isigma = 0
    if self.cumulative_multiply_observed_asu_count > 0:
      total_mean_I        = self.cumulative_I / self.cumulative_multiply_observed_asu_count
      total_mean_Isigma   = self.cumulative_Isigma / self.cumulative_multiply_observed_asu_count

    table_data.append(  [
        format_value("%3s",   "All"),
        format_value("%-13s", "                 "),
        format_value("%13s",  "[%d/%d]"%(self.cumulative_observed_asu_count, self.cumulative_theor_asu_count)),
        format_value("%5.2f", total_completeness),
        format_value("%6.2f", total_asu_multiplicity),
        format_value("%6.2f", total_obs_multiplicity),
        format_value("%6.2f", (0)), # if redundancy_to_edge is None else cumulative_n_pred/cumulative_theor)),
        format_value("%6d",   self.cumulative_observed_count),
        format_value("%6d",   self.cumulative_multiply_observed_asu_count),
        format_value("%6d",   (0)), #if redundancy_to_edge is None else flex.sum(redundancy_to_edge))),
        format_value("%8.0f", total_mean_I),
        format_value("%8.3f", total_mean_Isigma),
        format_value("%8.3f", total_unmerged_mean_Isigma),
        format_value("%8.3f", total_unmerged_stddev_Isigma),
        format_value("%8.3f", total_unmerged_skew_Isigma),
    ])
    table_data = table_utils.manage_columns(table_data, include_columns)

    return table_utils.format(table_data, has_header = 2, justify ='center', delim = " ")

class cross_correlation_resolution_bin(object):
  '''Storage class for parameters of a cross-correlation resolution bin'''
  def __init__(self,
               i_bin=None,
               theor_asu_count=0,
               observed_matching_asu_count=0,
               cross_correlation=None):
    adopt_init_args(self, locals())

class cross_correlation_table(object):
  '''Represents a table of cross-correlations for resolution bins'''
  def __init__(self):
    self.table = []
    self.cumulative_observed_matching_asu_count = 0
    self.cumulative_theor_asu_count = 0
    self.cumulative_cross_correlation = 0.0

class intensity_resolution_statistics(worker):
  '''Calculates hkl intensity statistics for resolution bins'''

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(intensity_resolution_statistics, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Intensity resolution statistics'

  def run(self, experiments, reflections):
    self.last_bin_incomplete = False
    self.suggested_resolution_scalar = -1.0

    reflections_before_stats_count = reflections.size()
    total_reflections_before_stats_count = self.mpi_helper.sum(reflections_before_stats_count)
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Total reflections before doing intensity statistics: %d"%(total_reflections_before_stats_count))

    self.logger.log_step_time("INTENSITY_STATISTICS")

    title = "\n                     Intensity Statistics (odd accepted experiments)\n"
    self.logger.log(title, rank_prepend=False)
    if self.mpi_helper.rank == 0:
      self.logger.main_log(title)
    reflections_odd = reflection_table_utils.select_odd_experiment_reflections(reflections)
    self.run_detail(reflections_odd)

    title = "\n                     Intensity Statistics (even accepted experiments)\n"
    self.logger.log(title, rank_prepend=False)
    if self.mpi_helper.rank == 0:
      self.logger.main_log(title)
    reflections_even = reflection_table_utils.select_even_experiment_reflections(reflections)
    self.run_detail(reflections_even)
    assert len(reflections_even) + len(reflections_odd) == len(reflections)

    title = "\n                     Intensity Statistics (all accepted experiments)\n"
    self.logger.log(title, rank_prepend=False)
    if self.mpi_helper.rank == 0:
      self.logger.main_log(title)
    self.run_detail(reflections)

    title = "\n                     CC 1/2, CC ISO\n"
    self.logger.log(title, rank_prepend=False)
    self.calculate_cc_int(reflections_odd, reflections_even)
    self.calculate_cc_iso(reflections)

    if self.mpi_helper.rank == 0:
      for entry in [(self.Total_CC_OneHalf_Table, "\n                     CC 1/2\n"),
                    (self.Total_CC_Iso_Table,     "\n                     CC ISO\n")]:
        Table = entry[0]
        Title = entry[1]
        if Table is not None:
          self.logger.main_log(Title)
          for row in Table.table:
            self.logger.main_log("%d/%d\t\t%f"%(row.observed_matching_asu_count, row.theor_asu_count, row.cross_correlation))
          self.logger.main_log("--------------------------------------------------------")
          self.logger.main_log("%d/%d\t\t%f\n"%(Table.cumulative_observed_matching_asu_count,
                                                     Table.cumulative_theor_asu_count,
                                                     Table.cumulative_cross_correlation))

    self.logger.log_step_time("INTENSITY_STATISTICS", True)

    return experiments, reflections

  def calculate_cc_iso(self, reflections):
    if self.params.statistics.cciso.mtz_file == None:
      self.Total_CC_Iso_Table = None
      return

    reflections_merged = reflection_table_utils.merge_reflections(reflections, self.params.merging.minimum_multiplicity)

    # Create target symmetry
    if self.params.merging.set_average_unit_cell:
      assert 'average_unit_cell' in (self.params.statistics).__dict__
      unit_cell = self.params.statistics.__phil_get__('average_unit_cell')
    else:
      unit_cell = self.params.scaling.unit_cell
    target_symm = symmetry(unit_cell=unit_cell, space_group_info = self.params.scaling.space_group)

    # Build miller array for experimental data
    miller_indices = miller.set(target_symm, reflections_merged['miller_index'], True)
    exp_intensities = miller.array(miller_indices,
                                   reflections_merged['intensity'],
                                   flex.double(reflections_merged.size(), 1.0))

    self.Total_CC_Iso_Table = self.calculate_cross_correlation(self.params.scaling.i_model, exp_intensities)

  def calculate_cc_int(self, odd_reflections, even_reflections):
    odd_reflections_merged = reflection_table_utils.merge_reflections(odd_reflections, self.params.merging.minimum_multiplicity)
    even_reflections_merged = reflection_table_utils.merge_reflections(even_reflections, self.params.merging.minimum_multiplicity)

    # Create target symmetry
    if self.params.merging.set_average_unit_cell:
      assert 'average_unit_cell' in (self.params.statistics).__dict__
      unit_cell = self.params.statistics.__phil_get__('average_unit_cell')
    else:
      unit_cell = self.params.scaling.unit_cell
    target_symm = symmetry(unit_cell=unit_cell, space_group_info = self.params.scaling.space_group)

    # Build miller arrays
    miller_indices_odd = miller.set(target_symm, odd_reflections_merged['miller_index'], True)
    intensities_odd = miller.array(miller_indices_odd,
                                   odd_reflections_merged['intensity'],
                                   flex.double(odd_reflections_merged.size(), 1.0))

    miller_indices_even = miller.set(target_symm, even_reflections_merged['miller_index'], True)
    intensities_even = miller.array(miller_indices_even, even_reflections_merged['intensity'],
                                    flex.double(even_reflections_merged.size(), 1.0))

    # Calculate crosss-correlation
    self.Total_CC_OneHalf_Table = self.calculate_cross_correlation(intensities_odd, intensities_even)

  def calculate_cross_correlation(self, miller_array_1, miller_array_2):
    # Get pre-created resolution binning objects from the parameters
    self.resolution_binner = self.params.statistics.resolution_binner
    self.hkl_resolution_bins = self.params.statistics.hkl_resolution_bins

    # How many bins do we have?
    n_bins = self.resolution_binner.n_bins_all() # (self.params.statistics.n_bins + 2), 2 - to account for the hkls outside of the binner resolution range

    # To enable MPI all-rank reduction, every rank must initialize statistics array(s), even if the rank doesn't have any reflections.
    self.cc_N         = flex.int(n_bins, 0)
    self.cc_sum_xx    = flex.double(n_bins, 0.0)
    self.cc_sum_xy    = flex.double(n_bins, 0.0)
    self.cc_sum_yy    = flex.double(n_bins, 0.0)
    self.cc_sum_x     = flex.double(n_bins, 0.0)
    self.cc_sum_y     = flex.double(n_bins, 0.0)

    # Find matching indices in the two data sets
    matching_indices = miller.match_multi_indices(miller_indices_unique = miller_array_1.indices(),
                                                  miller_indices = miller_array_2.indices())

    # Perform binned summations for all components of the cross-correlation formula
    for pair in matching_indices.pairs():

      hkl = miller_array_1.indices()[pair[0]]
      assert hkl == miller_array_2.indices()[pair[1]]

      if hkl in self.hkl_resolution_bins:
        i_bin = self.hkl_resolution_bins[hkl]

        I_x = miller_array_1.data()[pair[0]]
        I_y = miller_array_2.data()[pair[1]]

        self.cc_N[i_bin]        += 1
        self.cc_sum_xx[i_bin]   += I_x**2
        self.cc_sum_yy[i_bin]   += I_y**2
        self.cc_sum_xy[i_bin]   += I_x * I_y
        self.cc_sum_x[i_bin]    += I_x
        self.cc_sum_y[i_bin]    += I_y

    # Accumulate binned counts (cc_N) and sums (cc_sum) from all ranks
    all_ranks_cc_N          = self.mpi_helper.cumulative_flex(self.cc_N,      flex.int)
    all_ranks_cc_sum_xx     = self.mpi_helper.cumulative_flex(self.cc_sum_xx, flex.double)
    all_ranks_cc_sum_yy     = self.mpi_helper.cumulative_flex(self.cc_sum_yy, flex.double)
    all_ranks_cc_sum_xy     = self.mpi_helper.cumulative_flex(self.cc_sum_xy, flex.double)
    all_ranks_cc_sum_x      = self.mpi_helper.cumulative_flex(self.cc_sum_x,  flex.double)
    all_ranks_cc_sum_y      = self.mpi_helper.cumulative_flex(self.cc_sum_y,  flex.double)

    # Reduce all binned counts (cc_N) and sums (cc_sum) from all ranks
    if self.mpi_helper.rank == 0:
      return self.build_cross_correlation_table(
                                                all_ranks_cc_N,
                                                all_ranks_cc_sum_xx,
                                                all_ranks_cc_sum_yy,
                                                all_ranks_cc_sum_xy,
                                                all_ranks_cc_sum_x,
                                                all_ranks_cc_sum_y)
    else:
      return None

  def run_detail(self, reflections):
    # Get pre-created resolution binning objects from the parameters
    self.resolution_binner = self.params.statistics.resolution_binner
    self.hkl_resolution_bins = self.params.statistics.hkl_resolution_bins

    # How many bins do we have?
    n_bins = self.resolution_binner.n_bins_all() # (self.params.statistics.n_bins + 2), 2 - to account for the hkls outside of the binner resolution range

    # To enable MPI all-rank reduction, every rank must initialize statistics array(s), even if the rank doesn't have any reflections.
    self.I_sum      = flex.double(n_bins, 0.0) # a sum of the weighted mean intensities from all asu HKLs
    self.Isig_sum   = flex.double(n_bins, 0.0) # a sum of I/sigma from all asu HKLs
    self.Isig_list  = [flex.double() for _ in range(n_bins)]
    self.n_sum      = flex.int(n_bins, 0) # number of theoretically prediced asu hkls
    self.m_sum      = flex.int(n_bins, 0) # number of observed asu hkls
    self.mm_sum     = flex.int(n_bins, 0) # number of observed asu hkls with multiplicity > 1

    # Calculate, format and output statistics for each rank
    self.logger.log("Calculating intensity statistics...")
    self.calculate_intensity_statistics(reflections)
    Intensity_Table = self.build_intensity_table(
                                                I_sum = self.I_sum,
                                                Isig_sum = self.Isig_sum,
                                                n_sum = self.n_sum,
                                                m_sum = self.m_sum,
                                                mm_sum = self.mm_sum,
                                                unmerged_meanIsig = None,
                                                unmerged_stddevIsig = None,
                                                unmerged_skewIsig = None)
    if self.params.output.log_level == 0:
      self.logger.log(Intensity_Table.get_table_text(), rank_prepend=False)

    # Accumulate statistics from all ranks
    all_ranks_I_sum       = self.mpi_helper.cumulative_flex(self.I_sum, flex.double)
    all_ranks_Isig_sum    = self.mpi_helper.cumulative_flex(self.Isig_sum, flex.double)
    all_ranks_n_sum       = self.mpi_helper.cumulative_flex(self.n_sum, flex.int)
    all_ranks_m_sum       = self.mpi_helper.cumulative_flex(self.m_sum, flex.int)
    all_ranks_mm_sum      = self.mpi_helper.cumulative_flex(self.mm_sum, flex.int)

    all_ranks_unmerged_meanIsig = []
    all_ranks_unmerged_stddevIsig = []
    all_ranks_unmerged_skewIsig = []
    for bin_id in range(n_bins):
      all_ranks_isigi_list = self.mpi_helper.comm.gather(self.Isig_list[bin_id], 0)
      if self.mpi_helper.rank == 0:
        all_isigi = flex.double()
        for ranklist in all_ranks_isigi_list:
          all_isigi.extend(ranklist)
        stats = basic_statistics(all_isigi)
        all_ranks_unmerged_meanIsig.append(stats.mean)
        all_ranks_unmerged_stddevIsig.append(stats.bias_corrected_standard_deviation)
        all_ranks_unmerged_skewIsig.append(stats.skew)

    # Calculate, format and output all-rank total statistics
    if self.mpi_helper.rank == 0:
      Intensity_Table = self.build_intensity_table(
                                                  I_sum      = all_ranks_I_sum,
                                                  Isig_sum   = all_ranks_Isig_sum,
                                                  n_sum      = all_ranks_n_sum,
                                                  m_sum      = all_ranks_m_sum,
                                                  mm_sum     = all_ranks_mm_sum,
                                                  unmerged_meanIsig   = all_ranks_unmerged_meanIsig,
                                                  unmerged_stddevIsig = all_ranks_unmerged_stddevIsig,
                                                  unmerged_skewIsig   = all_ranks_unmerged_skewIsig)
      self.logger.main_log(Intensity_Table.get_table_text())
      if self.last_bin_incomplete:
        self.logger.main_log("Warning: the last resolution shell is incomplete. If your data was integrated to that resolution,\nconsider using scaling.resolution_scalar=%f or lower."%self.suggested_resolution_scalar)

  def calculate_intensity_statistics(self, reflections):
    '''Calculate statistics for hkl intensities distributed over resolution bins'''
    # The count of all reflections in the input reflection list, not restricted by any resolution limits, multiplicity thresholds, etc.
    zero_intensity_count_all      = (reflections['intensity.sum.value'] == 0.0).count(True)
    positive_intensity_count_all  = (reflections['intensity.sum.value'] >  0.0).count(True)
    negative_intensity_count_all  = (reflections['intensity.sum.value'] <  0.0).count(True)

    # The count of reflections filtered out from the input list by resolution limits
    zero_intensity_count_resolution_limited = 0
    positive_intensity_count_resolution_limited = 0
    negative_intensity_count_resolution_limited = 0

    # a list of reflections _actually_ used for the intenisty statistics calculations
    #used_reflections = flex.reflection_table()

    # How many bins do we have?
    n_bins = self.resolution_binner.n_bins_all() # (self.params.statistics.n_bins + 2), 2 - to account for the HKLs outside of the binner resolution range

    # Calculate auxiliary per-bin sums needed for the intensity statistics
    for refls in reflection_table_utils.get_next_hkl_reflection_table(reflections=reflections):
      if refls.size() == 0:
        break # unless the input "reflections" list is empty, generated "refls" lists cannot be empty

      hkl = refls[0]['miller_index_asymmetric']
      if hkl in self.hkl_resolution_bins:

        i_bin = self.hkl_resolution_bins[hkl]

        if i_bin > 0 and i_bin < n_bins - 1:
          zero_intensity_count_resolution_limited     += (refls['intensity.sum.value'] == 0.0).count(True)
          positive_intensity_count_resolution_limited += (refls['intensity.sum.value'] >  0.0).count(True)
          negative_intensity_count_resolution_limited += (refls['intensity.sum.value'] <  0.0).count(True)

          refls = refls.select(refls['intensity.sum.variance'] > 0.0)

          multiplicity = refls.size()

          if multiplicity >= self.params.merging.minimum_multiplicity:
            self.n_sum[i_bin] += 1
            self.m_sum[i_bin] += multiplicity
            self.mm_sum[i_bin] += 1
            weighted_intensity_array = refls['intensity.sum.value'] / refls['intensity.sum.variance']
            weights_array = flex.double(refls.size(), 1.0) / refls['intensity.sum.variance']
            self.I_sum[i_bin]     += flex.sum(weighted_intensity_array) / flex.sum(weights_array)
            self.Isig_sum[i_bin]  += flex.sum(weighted_intensity_array) / math.sqrt(flex.sum(weights_array))
            self.Isig_list[i_bin].extend(refls['intensity.sum.value.unmodified']/flex.sqrt(refls['intensity.sum.variance.unmodified']))
            #used_reflections.extend(refls)

    self.logger.log_step_time("INTENSITY_HISTOGRAM")

    # Accumulate intensities, which were used in the above statistics table, from all ranks
    #all_used_intensities = self.mpi_helper.extend_flex(used_reflections['intensity.sum.value'], flex.double)

    total_zero_intensity_count_all = self.mpi_helper.sum(zero_intensity_count_all)
    total_zero_intensity_count_resolution_limited = self.mpi_helper.sum(zero_intensity_count_resolution_limited)

    total_positive_intensity_count_all = self.mpi_helper.sum(positive_intensity_count_all)
    total_positive_intensity_count_resolution_limited = self.mpi_helper.sum(positive_intensity_count_resolution_limited)

    total_negative_intensity_count_all = self.mpi_helper.sum(negative_intensity_count_all)
    total_negative_intensity_count_resolution_limited = self.mpi_helper.sum(negative_intensity_count_resolution_limited)

    # Build a histogram of all intensities
    if self.mpi_helper.rank == 0:

      self.logger.main_log("Total reflections (I == 0.0): \t\t%d"%(total_zero_intensity_count_all))
      self.logger.main_log("Total reflections (I > 0.0): \t\t%d"%(total_positive_intensity_count_all))
      self.logger.main_log("Total reflections (I < 0.0): \t\t%d\n"%(total_negative_intensity_count_all))

      self.logger.main_log("Resolution-limited reflections (I == 0.0): \t\t%d"%(total_zero_intensity_count_resolution_limited))
      self.logger.main_log("Resolution-limited reflections (I > 0.0): \t\t%d"%(total_positive_intensity_count_resolution_limited))
      self.logger.main_log("Resolution-limited reflections (I < 0.0): \t\t%d\n"%(total_negative_intensity_count_resolution_limited))

      #self.histogram(all_used_intensities)
    self.logger.log_step_time("INTENSITY_HISTOGRAM", True)

  def build_intensity_table(self,
                            n_sum,
                            m_sum,
                            mm_sum,
                            I_sum,
                            Isig_sum,
                            unmerged_meanIsig,
                            unmerged_stddevIsig,
                            unmerged_skewIsig):
    '''Produce a table with hkl intensity statistics for resolution bins'''
    Intensity_Table = intensity_table()

    for i_bin in self.resolution_binner.range_used():
      theor_asu_count             = self.resolution_binner.counts()[i_bin]
      observed_asu_count          = n_sum[i_bin]
      observed_count              = m_sum[i_bin]
      multiply_observed_asu_count = mm_sum[i_bin]
      intensity_sum               = I_sum[i_bin]
      intensity_to_sigma_sum      = Isig_sum[i_bin]

      if observed_asu_count > 0:
        mean_I = mean_Isig = 0
        if multiply_observed_asu_count > 0:
          mean_I      = intensity_sum / multiply_observed_asu_count
          mean_Isig   = intensity_to_sigma_sum / multiply_observed_asu_count

        res_bin = intensity_table_bin(i_bin                 = i_bin,
                                  d_range                   = self.resolution_binner.bin_legend(i_bin=i_bin, show_bin_number=False, show_counts=False),
                                  d_min                     = self.resolution_binner.bin_d_min(i_bin),
                                  redundancy_asu            = observed_count / theor_asu_count,
                                  redundancy_obs            = observed_count / observed_asu_count,
                                  redundancy_to_edge        = 0, # TODO
                                  complete_tag              = "[%d/%d]" % (observed_asu_count, theor_asu_count),
                                  completeness              = observed_asu_count / theor_asu_count,
                                  measurements              = observed_count,
                                  multiply_measured_asu     = multiply_observed_asu_count,
                                  predictions               = None, # TODO
                                  mean_I                    = mean_I,
                                  mean_I_sigI               = mean_Isig,
                                  unmerged_meanIsig         = unmerged_meanIsig[i_bin] if unmerged_meanIsig is not None else 0,
                                  unmerged_stddevIsig       = unmerged_stddevIsig[i_bin] if unmerged_stddevIsig is not None else 0,
                                  unmerged_skewIsig         = unmerged_skewIsig[i_bin] if unmerged_skewIsig is not None else 0)

        Intensity_Table.table.append(res_bin)

      Intensity_Table.cumulative_observed_asu_count                   += observed_asu_count
      Intensity_Table.cumulative_observed_count                       += observed_count
      Intensity_Table.cumulative_theor_asu_count                      += theor_asu_count
      Intensity_Table.cumulative_multiply_observed_asu_count          += multiply_observed_asu_count
      Intensity_Table.cumulative_I                                    += intensity_sum
      Intensity_Table.cumulative_Isigma                               += intensity_to_sigma_sum

    '''
    # TODO
    if redundancy_to_edge is not None:
      cumulative_n_pred += flex.sum(sel_redundancy_pred)
      cumulative_pred   += redundancy_to_edge
    '''
    return Intensity_Table

  def histogram(self, data):
    from matplotlib import pyplot as plt
    nslots = 100
    histogram = flex.histogram(
                               data=data,
                               n_slots=nslots)
    out = StringIO()
    histogram.show(f=out, prefix="    ", format_cutoffs="%6.2f")
    self.logger.main_log(out.getvalue() + '\n' + "Total: %d"%data.size() + '\n')

    if False:
      fig = plt.figure()
      plt.bar(histogram.slot_centers(), histogram.slots(), align="center", width=histogram.slot_width())
      plt.show()

  def build_cross_correlation_table(self,
                                    count_array,
                                    sum_xx_array,
                                    sum_yy_array,
                                    sum_xy_array,
                                    sum_x_array,
                                    sum_y_array):

    Cross_Correlation_Table = cross_correlation_table()

    cumulative_observed_matching_asu_count = 0
    cumulative_theor_asu_count = 0
    cumulative_sum_xx  = 0.0
    cumulative_sum_yy  = 0.0
    cumulative_sum_xy  = 0.0
    cumulative_sum_x   = 0.0
    cumulative_sum_y   = 0.0

    for i_bin in self.resolution_binner.range_used():
      count   = count_array[i_bin]
      sum_xx  = sum_xx_array[i_bin]
      sum_yy  = sum_yy_array[i_bin]
      sum_xy  = sum_xy_array[i_bin]
      sum_x   = sum_x_array[i_bin]
      sum_y   = sum_y_array[i_bin]

      cross_correlation = self.cross_correlation_formula(count,
                                                         sum_xx,
                                                         sum_yy,
                                                         sum_xy,
                                                         sum_x,
                                                         sum_y);

      Cross_Correlation_Table.table.append(
                    cross_correlation_resolution_bin(
                            i_bin = i_bin,
                            theor_asu_count = self.resolution_binner.counts()[i_bin],
                            observed_matching_asu_count = count,
                            cross_correlation = cross_correlation))

      cumulative_observed_matching_asu_count += count
      cumulative_theor_asu_count += self.resolution_binner.counts()[i_bin]
      cumulative_sum_xx  += sum_xx
      cumulative_sum_yy  += sum_yy
      cumulative_sum_xy  += sum_xy
      cumulative_sum_x   += sum_x
      cumulative_sum_y   += sum_y

    Cross_Correlation_Table.cumulative_observed_matching_asu_count = cumulative_observed_matching_asu_count
    Cross_Correlation_Table.cumulative_theor_asu_count = cumulative_theor_asu_count
    Cross_Correlation_Table.cumulative_cross_correlation = self.cross_correlation_formula(cumulative_observed_matching_asu_count,
                                                                                          cumulative_sum_xx,
                                                                                          cumulative_sum_yy,
                                                                                          cumulative_sum_xy,
                                                                                          cumulative_sum_x,
                                                                                          cumulative_sum_y)
    return Cross_Correlation_Table


  def cross_correlation_formula(self,
                                count,
                                sum_xx,
                                sum_yy,
                                sum_xy,
                                sum_x,
                                sum_y):

    numerator = (count * sum_xy - sum_x * sum_y)
    denominator = (math.sqrt(count * sum_xx - sum_x**2) * math.sqrt(count * sum_yy - sum_y**2))
    cross_correlation = 0.0
    if denominator != 0.0:
      cross_correlation = numerator / denominator

    return cross_correlation

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(intensity_resolution_statistics)
