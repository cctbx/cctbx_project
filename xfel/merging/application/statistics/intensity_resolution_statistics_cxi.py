from __future__ import absolute_import, division, print_function
from six.moves import range
from xfel.merging.application.worker import worker
from dials.array_family import flex
import math
from libtbx.str_utils import format_value
from libtbx import table_utils
from cctbx.crystal import symmetry
from cctbx.miller import binned_data
import os
from iotbx import mtz

class intensity_resolution_statistics_cxi(worker):
  '''Calculates hkl intensity statistics for resolution bins using adapted cxi-xmerge code'''

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(intensity_resolution_statistics_cxi, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Intensity resolution statistics (cxi-xmerge method)'

  def run(self, experiments, reflections):
    self.logger.log_step_time("INTENSITY_STATISTICS CXI")

    reflections_before_stats_count = reflections.size()
    total_reflections_before_stats_count = self.mpi_helper.sum(reflections_before_stats_count)
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Total reflections before doing intensity statistics (cxi-xmerge method): %d"%(total_reflections_before_stats_count))

    if self.mpi_helper.rank == 0:
      self.run_cc()

    self.logger.log_step_time("INTENSITY_STATISTICS CXI", True)

    return experiments, reflections

  def run_cc(self):
    uniform, selected_uniform, have_iso_ref = self.load_cc_data()

    include_negatives = True
    if have_iso_ref:
      slope, offset, corr_iso, N_iso = self.correlation(
        selected_uniform[1], selected_uniform[0], include_negatives)
      self.logger.main_log("C.C. iso is %.1f%% on %d indices"%(100 * corr_iso, N_iso))

    slope, offset, corr_int, N_int = self.correlation(selected_uniform[2], selected_uniform[3], include_negatives)
    self.logger.main_log("C.C. int is %.1f%% on %d indices"%(100.*corr_int, N_int))

    if have_iso_ref:
      binned_cc_ref, binned_cc_ref_N = self.binned_correlation(
        selected_uniform[1], selected_uniform[0], include_negatives)
      #binned_cc_ref.show(f=output)

      ref_scale = self.scale_factor(
        selected_uniform[1], selected_uniform[0],
        weights=flex.pow(selected_uniform[1].sigmas(), -2),
        use_binning=True)
      #ref_scale.show(f=output)

      ref_riso = self.r1_factor(
        selected_uniform[1], selected_uniform[0],
        scale_factor=ref_scale, use_binning=True)
      #ref_riso.show(f=output)

      ref_scale_all = self.scale_factor(
        selected_uniform[1], selected_uniform[0],
        weights=flex.pow(selected_uniform[1].sigmas(), -2))

      ref_riso_all = self.r1_factor(
        selected_uniform[1], selected_uniform[0],
        scale_factor=ref_scale_all)

    binned_cc_int,binned_cc_int_N = self.binned_correlation(
      #selected_uniform[2], selected_uniform[3], params.include_negatives)
      selected_uniform[2], selected_uniform[3], True)
    #binned_cc_int.show(f=output)

    oe_scale = self.scale_factor(selected_uniform[2],selected_uniform[3],
      weights = flex.pow(selected_uniform[2].sigmas(),-2)
              + flex.pow(selected_uniform[3].sigmas(),-2),
      use_binning=True)
    #oe_scale.show(f=output)

    oe_rint = self.r1_factor(selected_uniform[2],selected_uniform[3],
                         scale_factor = oe_scale, use_binning=True)
    #oe_rint.show(f=output)

    oe_rsplit = self.r_split(selected_uniform[2], selected_uniform[3],
                        use_binning=True)

    oe_scale_all = self.scale_factor(selected_uniform[2],selected_uniform[3],
      weights = flex.pow(selected_uniform[2].sigmas(),-2)
              + flex.pow(selected_uniform[3].sigmas(),-2),)

    oe_rint_all = self.r1_factor(selected_uniform[2],selected_uniform[3],
                         scale_factor = oe_scale_all)
    oe_rsplit_all = self.r_split(selected_uniform[2], selected_uniform[3])
    if have_iso_ref:
      self.logger.main_log("R factors Riso = %.1f%%, Rint = %.1f%%"%(100.*ref_riso_all, 100.*oe_rint_all))
    else:
      self.logger.main_log("R factor Rint = %.1f%%"%(100.*oe_rint_all))

    split_sigma_data = self.split_sigma_test(selected_uniform[2],selected_uniform[3],
                                        scale=oe_scale,use_binning=True,show_plot=False)
    split_sigma_data_all = self.split_sigma_test(selected_uniform[2],selected_uniform[3],
                                            scale=oe_scale_all,use_binning=False,show_plot=False)

    self.logger.main_log('')
    if self.params.scaling.model_reindex_op == "h,k,l":
      self.logger.main_log("Table of Scaling Results:")
    else:
      self.logger.main_log("Table of Scaling Results with Model Reindexing as %s:"%reindexing_op)

    from libtbx import table_utils
    table_header = ["","","","CC"," N","CC"," N","R","R","R","Scale","Scale","SpSig"]
    table_header2 = ["Bin","Resolution Range","Completeness","int","int","iso","iso","int","split","iso","int","iso","Test"]
    table_data = []
    table_data.append(table_header)
    table_data.append(table_header2)

    items = binned_cc_int.binner.range_used()

    # XXX Make it clear what the completeness here actually is!
    cumulative_counts_given = 0
    cumulative_counts_complete = 0
    for bin in items:
      table_row = []
      table_row.append("%3d"%bin)
      table_row.append("%-13s"%binned_cc_int.binner.bin_legend(i_bin=bin,show_bin_number=False,show_bin_range=False,
                                                   show_d_range=True, show_counts=False))
      table_row.append("%13s"%binned_cc_int.binner.bin_legend(i_bin=bin,show_bin_number=False,show_bin_range=False,
                                                   show_d_range=False, show_counts=True))
      cumulative_counts_given += binned_cc_int.binner._counts_given[bin]
      cumulative_counts_complete += binned_cc_int.binner._counts_complete[bin]
      table_row.append("%.1f%%"%(100.*binned_cc_int.data[bin]))
      table_row.append("%7d"%(binned_cc_int_N.data[bin]))

      if have_iso_ref and binned_cc_ref.data[bin] is not None:
        table_row.append("%.1f%%" % (100 * binned_cc_ref.data[bin]))
      else:
        table_row.append("--")

      if have_iso_ref and binned_cc_ref_N.data[bin] is not None:
        table_row.append("%6d" % (binned_cc_ref_N.data[bin]))
      else:
        table_row.append("--")

      if oe_rint.data[bin] is not None:
        table_row.append("%.1f%%"%(100.*oe_rint.data[bin]))
      else:
        table_row.append("--")

      if oe_rsplit.data[bin] is not None:
        table_row.append("%.1f%%" % (100 * oe_rsplit.data[bin]))
      else:
        table_row.append("--")

      if have_iso_ref and ref_riso.data[bin] is not None:
        table_row.append("%.1f%%" % (100 * ref_riso.data[bin]))
      else:
        table_row.append("--")

      if oe_scale.data[bin] is not None:
        table_row.append("%.3f"%oe_scale.data[bin])
      else:
        table_row.append("--")

      if have_iso_ref and ref_scale.data[bin] is not None:
        table_row.append("%.3f" % ref_scale.data[bin])
      else:
        table_row.append("--")

      if split_sigma_data.data[bin] is not None:
        table_row.append("%.4f" % split_sigma_data.data[bin])
      else:
        table_row.append("--")

      table_data.append(table_row)
    table_data.append([""]*len(table_header))

    table_row = [format_value("%3s",   "All"),
                 format_value("%-13s", "                 "),
                 format_value("%13s",  "[%d/%d]"%(cumulative_counts_given,
                                                  cumulative_counts_complete)),
                 format_value("%.1f%%", 100 * corr_int),
                 format_value("%7d", N_int)]

    if have_iso_ref:
      table_row.extend((format_value("%.1f%%", 100 * corr_iso),
                        format_value("%6d", N_iso)))
    else:
      table_row.extend(("--", "--"))

    table_row.extend((format_value("%.1f%%", 100 * oe_rint_all),
                      format_value("%.1f%%", 100 * oe_rsplit_all)))
    if have_iso_ref:
      table_row.append(format_value("%.1f%%", 100 * ref_riso_all))
    else:
      table_row.append("--")

    table_row.append(format_value("%.3f", oe_scale_all))
    if have_iso_ref:
      table_row.append(format_value("%.3f", ref_scale_all))
    else:
      table_row.append("--")

    if split_sigma_data_all is not None:
      table_row.append("%.1f" % split_sigma_data_all)
    else:
      table_row.append("--")

    table_data.append(table_row)

    self.logger.main_log(' ')
    self.logger.main_log(table_utils.format(table_data,has_header=2,justify='center',delim=" "))
    self.logger.main_log("""CCint is the CC-1/2 defined by Diederichs; correlation between odd/even images.
    Similarly, Scale int and R int are the scaling factor and scaling R factor between odd/even images.
    "iso" columns compare the whole XFEL dataset to the isomorphous reference.""")

    self.logger.main_log("Niso: result vs. reference common set")

    if have_iso_ref:
      assert N_iso == flex.sum(flex.double([x for x in binned_cc_ref_N.data if x is not None]))
    assert N_int == flex.sum(flex.double([x for x in binned_cc_int_N.data if x is not None]))

    # TODO: how is plotting handled in the new phil design?
    '''
    if params.scaling.show_plots:
      from matplotlib import pyplot as plt
      plt.plot(flex.log(selected_uniform[-2].data()),
               flex.log(selected_uniform[-1].data()), 'r.')
      plt.show()
      if have_iso_ref:
        plt.plot(flex.log(selected_uniform[0].data()),
                 flex.log(selected_uniform[1].data()), 'r.')
        plt.show()
    '''
    self.logger.main_log(' ')

  def scale_factor(self, this, other, weights=None, cutoff_factor=None,
                     use_binning=False):
      """
      The analytical expression for the least squares scale factor.

      K = sum(w * yo * yc) / sum(w * yc^2)

      If the optional cutoff_factor argument is provided, only the reflections
      whose magnitudes are greater than cutoff_factor * max(yo) will be included
      in the calculation.
      """
      assert not use_binning or this.binner() is not None
      if use_binning: assert cutoff_factor is None
      assert other.size() == this.data().size()
      if not use_binning:
        if this.data().size() == 0: return None
        obs = this.data()
        calc = other.data()
        if cutoff_factor is not None:
          assert cutoff_factor < 1
          sel = obs >= flex.max(this.data()) * cutoff_factor
          obs = obs.select(sel)
          calc = calc.select(sel)
          if weights is not None:
            weights = weights.select(sel)
        if weights is None:
          return flex.sum(obs*calc) / flex.sum(flex.pow2(calc))
        else:
          return flex.sum(weights * obs * calc) \
               / flex.sum(weights * flex.pow2(calc))
      results = []
      for i_bin in this.binner().range_all():
        sel = this.binner().selection(i_bin)
        weights_sel = None
        if weights is not None:
          weights_sel = weights.select(sel)
        results.append(
          self.scale_factor(this.select(sel),other.select(sel), weights_sel))
      return binned_data(binner=this.binner(), data=results, data_fmt="%7.4f")

  def split_sigma_test(self, this, other, scale, use_binning=False, show_plot=False):
    """
    Calculates the split sigma ratio test by Peter Zwart:
    ssr = sum( (Iah-Ibh)^2 ) / sum( sigma_ah^2 + sigma_bh^2)

    where Iah and Ibh are merged intensities for a given hkl from two halves of
    a dataset (a and b). Likewise for sigma_ah and sigma_bh.

    ssr (split sigma ratio) should approximately equal 1 if the errors are correctly estimated.
    """

    assert other.size() == this.data().size()
    assert (this.indices() == other.indices()).all_eq(True)
    assert not use_binning or this.binner() is not None

    if use_binning:
      results = []
      for i_bin in this.binner().range_all():
        sel = this.binner().selection(i_bin)
        i_this = this.select(sel)
        i_other = other.select(sel)
        scale_rel = scale.data[i_bin]
        if i_this.size() == 0:
          results.append(None)
        else:
          results.append(self.split_sigma_test(i_this,i_other,scale=scale_rel,show_plot=show_plot))
      return binned_data(binner=this.binner(), data=results, data_fmt="%7.4f")

    a_data = this.data(); b_data = scale * other.data()
    a_sigmas = this.sigmas(); b_sigmas = scale * other.sigmas()

    if show_plot:
      """
      # Diagnostic use of the (I - <I>) / sigma distribution, should have mean=0, std=1
      a_variance = a_sigmas * a_sigmas
      b_variance = b_sigmas * b_sigmas
      mean_num = (a_data/ (a_variance) ) + (b_data/ (b_variance) )
      mean_den = (1./ (a_variance) ) + (1./ (b_variance) )
      mean_values = mean_num / mean_den

      delta_I_a = a_data - mean_values
      normal_a = delta_I_a / (a_sigmas)
      stats_a = flex.mean_and_variance(normal_a)
      print "\nA mean %7.4f std %7.4f"%(stats_a.mean(),stats_a.unweighted_sample_standard_deviation())
      order_a = flex.sort_permutation(normal_a)

      delta_I_b = b_data - mean_values
      normal_b = delta_I_b / (b_sigmas)
      stats_b = flex.mean_and_variance(normal_b)
      print "B mean %7.4f std %7.4f"%(stats_b.mean(),stats_b.unweighted_sample_standard_deviation())
      order_b = flex.sort_permutation(normal_b)
      # plots for debugging
      from matplotlib import pyplot as plt
      plt.plot(range(len(order_a)),normal_a.select(order_a),"b.")
      plt.plot(range(len(order_b)),normal_b.select(order_b),"r.")
      plt.show()
      """
      from cctbx.examples.merging.sigma_correction import ccp4_model
      Correction = ccp4_model()
      Correction.plots(a_data, b_data, a_sigmas, b_sigmas)
      #a_new_variance,b_new_variance = Correction.optimize(a_data, b_data, a_sigmas, b_sigmas)
      #Correction.plots(a_data, b_data, flex.sqrt(a_new_variance), flex.sqrt(b_new_variance))

    n = flex.pow(a_data - b_data,2)
    d = flex.pow(a_sigmas,2)+flex.pow(b_sigmas,2)

    return flex.sum(n)/flex.sum(d)

  def r1_factor(self, this, other, scale_factor=None, assume_index_matching=False,
                  use_binning=False):
      r"""Get the R1 factor according to this formula

      .. math::
         R1 = \dfrac{\sum{||F| - k|F'||}}{\sum{|F|}}

      where F is this.data() and F' is other.data() and
      k is the factor to put F' on the same scale as F"""
      assert not use_binning or this.binner() is not None
      assert other.indices().size() == this.indices().size()
      if not use_binning:
        if this.data().size() == 0: return None
        if (assume_index_matching):
          o, c = this, other
        else:
          o, c = this.common_sets(other=other, assert_no_singles=True)
        o  = flex.abs(o.data())
        c = flex.abs(c.data())
        if (scale_factor is None):
          den = flex.sum(c * c)
          if (den != 0):
            c *= (flex.sum(o * c) / den)
        elif (scale_factor is not None):
          c *= scale_factor
        return flex.sum(flex.abs(o - c)) / flex.sum(o)
      results = []
      for i_bin in this.binner().range_all():
        sel = this.binner().selection(i_bin)
        results.append(self.r1_factor(this.select(sel),
          other.select(sel), scale_factor.data[i_bin], assume_index_matching))
      return binned_data(binner=this.binner(), data=results, data_fmt="%7.4f")

  def r_split(self, this, other, assume_index_matching=False, use_binning=False):
      # Used in Boutet et al. (2012), which credit it to Owen et al
      # (2006).  See also R_mrgd_I in Diederichs & Karplus (1997)?
      # Barends cites Collaborative Computational Project Number 4. The
      # CCP4 suite: programs for protein crystallography. Acta
      # Crystallogr. Sect. D-Biol. Crystallogr. 50, 760-763 (1994) and
      # White, T. A. et al. CrystFEL: a software suite for snapshot
      # serial crystallography. J. Appl. Cryst. 45, 335-341 (2012).
      if not use_binning:
        assert other.indices().size() == this.indices().size()
        if this.data().size() == 0:
          return None

        if assume_index_matching:
          (o, c) = (this, other)
        else:
          (o, c) = this.common_sets(other=other, assert_no_singles=True)

        # The case where the denominator is less or equal to zero is
        # pathological and should never arise in practice.
        den = flex.sum(flex.abs(o.data() + c.data()))
        assert den > 0
        return math.sqrt(2) * flex.sum(flex.abs(o.data() - c.data())) / den

      assert this.binner is not None
      results = []
      for i_bin in this.binner().range_all():
        sel = this.binner().selection(i_bin)
        results.append(self.r_split(this.select(sel), other.select(sel),
          assume_index_matching=assume_index_matching,
          use_binning=False))
      return binned_data(binner=this.binner(), data=results, data_fmt='%7.4f')

  def binned_correlation(self, this, other, include_negatives=False):
      results = []
      bin_count = []
      for i_bin in this.binner().range_all():
        sel = this.binner().selection(i_bin)
        if sel.count(True)==0:
          results.append(0.)
          bin_count.append(0.)
          continue
        result_tuple = self.correlation(this.select(sel), other.select(sel), include_negatives)
        results.append(result_tuple[2])
        bin_count.append(result_tuple[3])
        # plots for debugging
        #from matplotlib import pyplot as plt
        #plt.plot(flex.log(this.select(sel).data()),flex.log(other.select(sel).data()),"b.")
        #plt.show()

      return binned_data(binner=this.binner(), data=results, data_fmt="%7.4f"),\
             binned_data(binner=this.binner(), data=bin_count, data_fmt="%7d")

  def correlation(self, this, other, include_negatives=False):
      N = 0
      sum_xx = 0
      sum_xy = 0
      sum_yy = 0
      sum_x = 0
      sum_y = 0
      for idx in range(this.indices().size()):

        assert this.indices()[idx]==other.indices()[idx]
        I_r = other.data()[idx]
        I_o = this.data()[idx]
        #assert I_r >= 0. or I_o >= 0.
        #why does this go from 81% to 33% when uniform selection is made?

        if not include_negatives and (I_r < 0. or I_o < 0.): continue
        #print "%15s %15s %10.0f %10.0f"%(
      #this.indices()[idx], this.indices()[idx],
      #this.data()[idx], other.data()[idx],)
        N      += 1
        sum_xx += I_r**2
        sum_yy += I_o**2
        sum_xy += I_r * I_o
        sum_x  += I_r
        sum_y  += I_o
      # Linearly fit I_r to I_o, i.e. find slope and offset such that
      # I_o = slope * I_r + offset, optimal in a least-squares sense.
      if N < 2:  return 0,0,0,N
      slope = (N * sum_xy - sum_x * sum_y) / (N * sum_xx - sum_x**2)
      offset = (sum_xx * sum_y - sum_x * sum_xy) / (N * sum_xx - sum_x**2)
      corr  = (N * sum_xy - sum_x * sum_y) / (math.sqrt(N * sum_xx - sum_x**2) *
               math.sqrt(N * sum_yy - sum_y**2))
      return slope,offset,corr,N

  def load_mpi_merge_data(self, filename_postfix):
    filepath = os.path.join(self.params.output.output_dir, "%s_%s.mtz"%(self.params.output.prefix, filename_postfix))
    data = mtz.object(filepath)
    return data

  def load_cc_data(self):
    if self.params.scaling.model_reindex_op != "h,k,l":
      self.logger.main_log('''Recalculating after reindexing the new data with %s
       (it is necessary to pick which indexing choice gives the sensible CC iso):'''%self.params.scaling.model_reindex_op)

    try:
      model_file_path = self.params.statistics.cciso.mtz_file
      assert model_file_path.endswith(("mtz", "sf.cif"))
      # support both old-style *.mtz and structure factor *-sf.cif
      from iotbx import reflection_file_reader
      data_SR = reflection_file_reader.any_reflection_file(
                file_name = model_file_path)
      have_iso_ref = True
    except (RuntimeError, AttributeError): # Attribute error if model_file_path is None
      data_SR = None
      have_iso_ref = False

    data_d0 = self.load_mpi_merge_data('all')
    data_d1 = self.load_mpi_merge_data('odd')
    data_d2 = self.load_mpi_merge_data('even')

    uniform = []
    for idx, item in enumerate([data_SR, data_d0, data_d1, data_d2]):
      if not have_iso_ref and idx == 0:
        uniform.append(None)
        continue
      #item.show_summary()
      self.logger.main_log("-------------------------------")
      for array in item.as_miller_arrays():
         this_label = array.info().label_string().lower()
         self.logger.main_log(this_label + ' ' + self.params.statistics.cciso.mtz_column_F)
         if this_label.find("fobs") >= 0:
           self.logger.main_log(this_label + ' ' + str(array.observation_type()))
           uniform.append(array.as_intensity_array())
           break
         if this_label.find("iobs") >= 0:
           """This test is added for the use case of unmerged anomalous data without
              an isomorphous reference (in other words, mark1).  Without this,
              the unmerged reflections are not picked up by the cc comparison.
              Indicates that this section probably has to be reanalyzed and redesigned.
           """
           self.logger.main_log(this_label + ' ' + str(array.observation_type()))
           uniform.append(array.as_intensity_array())
           break
         if this_label.find("imean") >= 0:
           self.logger.main_log(this_label + ' ' + str(array.observation_type()))
           uniform.append(array.as_intensity_array())
           break
         if self.params.statistics.cciso.mtz_column_F.lower() in this_label:
           self.logger.main_log(this_label + ' ' + str(array.observation_type()))
           uniform.append(array.as_intensity_array())
           break

    # If necesssary, generate Bijvoet mates for the isomorphous reference.
    if have_iso_ref \
       and not self.params.merging.merge_anomalous \
       and not uniform[0].anomalous_flag():
        uniform[0] = uniform[0].generate_bijvoet_mates()

    # TODO: Shouldn't reindexing be done as part of filtering?

    #for x in [1,2,3]:
    #  # reindex the experimental data
    #  uniform[x] = uniform[x].change_basis(reindexing_op).map_to_asu()

    d_max_min = uniform[1].d_max_min()
    if have_iso_ref:
      sgi = uniform[0].space_group_info()
    else:
      sgi = self.params.scaling.space_group

    unit_cell = uniform[1].unit_cell()
    unit_cell_formatted = "(%.6f, %.6f, %.6f, %.3f, %.3f, %.3f)"\
                      %(unit_cell.parameters()[0], unit_cell.parameters()[1], unit_cell.parameters()[2], \
                        unit_cell.parameters()[3], unit_cell.parameters()[4], unit_cell.parameters()[5])
    self.logger.main_log("\nUnit cell: %s\n"%unit_cell_formatted)

    for x in [0,1,2,3]:
      if not have_iso_ref and x == 0:
        continue
      self.logger.main_log("%6d indices:"%uniform[x].size() + ' ' + {0:"Reference intensities",
                       1:"Merged structure factors",
                       2:"Semi-dataset 1",
                       3:"Semi-dataset 2"}[x])

      uniform[x] = uniform[x].customized_copy(
        crystal_symmetry = symmetry(unit_cell=uniform[1].unit_cell(), space_group_info=sgi),
        ).resolution_filter(d_min = d_max_min[1], d_max = self.params.merging.d_max,
        ).complete_array(d_min = d_max_min[1], d_max = self.params.merging.d_max).map_to_asu()

    self.logger.main_log("%6d indices: An asymmetric unit in the resolution interval %.2f - %.2f Angstrom"%(
       uniform[1].size(), d_max_min[0], d_max_min[1]))

    if have_iso_ref:
      uniform[0] = uniform[0].common_set(uniform[1])
      assert len(uniform[0].indices()) == len(uniform[1].indices())

    uniform[2] = uniform[2].common_set(uniform[1])
    uniform[3] = uniform[3].common_set(uniform[1])

    self.logger.main_log("-------------------------------")

    NBIN = self.params.statistics.n_bins

    for x in [0, 1, 2, 3]:
      if not have_iso_ref and x == 0:
        continue
      uniform[x].setup_binner(n_bins=NBIN)

    for x in range(len(uniform[1].indices())):
      if have_iso_ref:
        assert uniform[0].indices()[x] == uniform[1].indices()[x]
      assert uniform[1].indices()[x] == uniform[2].indices()[x]
      assert uniform[2].indices()[x] == uniform[3].indices()[x]

    selected_uniform = []
    if have_iso_ref:
      # quickly circumvent the odd case where the reference intensities have no sigmas
      # (which is the case for model F's)
      if uniform[0].sigmas() is None:
        uniform[0].set_sigmas( uniform[0].data() )
      uniformA = (uniform[0].sigmas() > 0.).__and__(uniform[1].sigmas() > 0.)

      for x in [0, 1]:
        selected_uniform.append(uniform[x].select(uniformA))
        selected_uniform[x].setup_binner(
          d_max=(self.params.merging.d_max or 100000), d_min=self.params.merging.d_min, n_bins=NBIN)

    else:
      selected_uniform = [None, None]

    uniformB = (uniform[2].sigmas() > 0.).__and__(uniform[3].sigmas() > 0.)

    for x in [2, 3]:
      selected_uniform.append(uniform[x].select(uniformB))
      selected_uniform[x].setup_binner(
        d_max=(self.params.merging.d_max or 100000), d_min=self.params.merging.d_min, n_bins=NBIN)

    return uniform, selected_uniform, have_iso_ref

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(intensity_resolution_statistics_cxi)
