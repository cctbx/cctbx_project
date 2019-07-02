# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

from __future__ import absolute_import, division, print_function
from six.moves import range

import math
from libtbx.str_utils import format_value
from iotbx import mtz
from cctbx import crystal
from scitbx.array_family import flex
from cctbx.miller import binned_data

# A possible development path 9/11/2015 NKS
# can we first filter on the outliers to avoid outliers in the I - <I> plot and check that in
# then figure out what our target function is: delta or SpSig?
# figure out what happens to these statistics when confronted with known normal distributions
# thus settling the issues of factors of two and root(two)
# then figure out some empirical correction that fits the plot.  Doesn't have to be sdfac/sdb/adadd just has to be smooth.

def scale_factor(self, other, weights=None, cutoff_factor=None,
                   use_binning=False):
    """
    The analytical expression for the least squares scale factor.

    K = sum(w * yo * yc) / sum(w * yc^2)

    If the optional cutoff_factor argument is provided, only the reflections
    whose magnitudes are greater than cutoff_factor * max(yo) will be included
    in the calculation.
    """
    assert not use_binning or self.binner() is not None
    if use_binning: assert cutoff_factor is None
    assert other.size() == self.data().size()
    if not use_binning:
      if self.data().size() == 0: return None
      obs = self.data()
      calc = other.data()
      if cutoff_factor is not None:
        assert cutoff_factor < 1
        sel = obs >= flex.max(self.data()) * cutoff_factor
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
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      weights_sel = None
      if weights is not None:
        weights_sel = weights.select(sel)
      results.append(
        scale_factor(self.select(sel),other.select(sel), weights_sel))
    return binned_data(binner=self.binner(), data=results, data_fmt="%7.4f")

def split_sigma_test(self, other, scale, use_binning=False, show_plot=False):
  """
  Calculates the split sigma ratio test by Peter Zwart:
  ssr = sum( (Iah-Ibh)^2 ) / sum( sigma_ah^2 + sigma_bh^2)

  where Iah and Ibh are merged intensities for a given hkl from two halves of
  a dataset (a and b). Likewise for sigma_ah and sigma_bh.

  ssr (split sigma ratio) should approximately equal 1 if the errors are correctly estimated.
  """

  assert other.size() == self.data().size()
  assert (self.indices() == other.indices()).all_eq(True)
  assert not use_binning or self.binner() is not None

  if use_binning:
    results = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      i_self = self.select(sel)
      i_other = other.select(sel)
      scale_rel = scale.data[i_bin]
      if i_self.size() == 0:
        results.append(None)
      else:
        results.append(split_sigma_test(i_self,i_other,scale=scale_rel,show_plot=show_plot))
    return binned_data(binner=self.binner(), data=results, data_fmt="%7.4f")

  a_data = self.data(); b_data = scale * other.data()
  a_sigmas = self.sigmas(); b_sigmas = scale * other.sigmas()

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

def r1_factor(self, other, scale_factor=None, assume_index_matching=False,
                use_binning=False):
    """Get the R1 factor according to this formula

    .. math::
       R1 = \dfrac{\sum{||F| - k|F'||}}{\sum{|F|}}

    where F is self.data() and F' is other.data() and
    k is the factor to put F' on the same scale as F"""
    assert not use_binning or self.binner() is not None
    assert other.indices().size() == self.indices().size()
    if not use_binning:
      if self.data().size() == 0: return None
      if (assume_index_matching):
        o, c = self, other
      else:
        o, c = self.common_sets(other=other, assert_no_singles=True)
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
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      results.append(r1_factor(self.select(sel),
        other.select(sel), scale_factor.data[i_bin], assume_index_matching))
    return binned_data(binner=self.binner(), data=results, data_fmt="%7.4f")

def r_split(self, other, assume_index_matching=False, use_binning=False):
    # Used in Boutet et al. (2012), which credit it to Owen et al
    # (2006).  See also R_mrgd_I in Diederichs & Karplus (1997)?
    # Barends cites Collaborative Computational Project Number 4. The
    # CCP4 suite: programs for protein crystallography. Acta
    # Crystallogr. Sect. D-Biol. Crystallogr. 50, 760-763 (1994) and
    # White, T. A. et al. CrystFEL: a software suite for snapshot
    # serial crystallography. J. Appl. Cryst. 45, 335–341 (2012).

    if not use_binning:
      assert other.indices().size() == self.indices().size()
      if self.data().size() == 0:
        return None

      if assume_index_matching:
        (o, c) = (self, other)
      else:
        (o, c) = self.common_sets(other=other, assert_no_singles=True)

      # The case where the denominator is less or equal to zero is
      # pathological and should never arise in practice.
      den = flex.sum(flex.abs(o.data() + c.data()))
      assert den > 0
      return math.sqrt(2) * flex.sum(flex.abs(o.data() - c.data())) / den

    assert self.binner is not None
    results = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      results.append(r_split(self.select(sel), other.select(sel),
        assume_index_matching=assume_index_matching,
        use_binning=False))
    return binned_data(binner=self.binner(), data=results, data_fmt='%7.4f')

def binned_correlation(self,other,include_negatives=False):
    results = []
    bin_count = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      if sel.count(True)==0:
        results.append(0.)
        bin_count.append(0.)
        continue
      result_tuple = correlation(self.select(sel),other.select(sel),include_negatives)
      results.append(result_tuple[2])
      bin_count.append(result_tuple[3])
      # plots for debugging
      #from matplotlib import pyplot as plt
      #plt.plot(flex.log(self.select(sel).data()),flex.log(other.select(sel).data()),"b.")
      #plt.show()

    return binned_data(binner=self.binner(), data=results, data_fmt="%7.4f"),\
           binned_data(binner=self.binner(), data=bin_count, data_fmt="%7d")

def correlation(self,other, include_negatives=False):
    N = 0
    sum_xx = 0
    sum_xy = 0
    sum_yy = 0
    sum_x = 0
    sum_y = 0
    for idx in range(self.indices().size()):

      assert self.indices()[idx]==other.indices()[idx]
      I_r = other.data()[idx]
      I_o = self.data()[idx]
      #assert I_r >= 0. or I_o >= 0.
      #why does this go from 81% to 33% when uniform selection is made?

      if not include_negatives and (I_r < 0. or I_o < 0.): continue
      #print "%15s %15s %10.0f %10.0f"%(
    #self.indices()[idx], self.indices()[idx],
    #self.data()[idx], other.data()[idx],)
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

def load_cc_data(params,reindexing_op,output):
  if reindexing_op is not "h,k,l":
    print("""Recalculating after reindexing the new data with %s
     (it is necessary to pick which indexing choice gives the sensible CC iso):"""%reindexing_op)
  try:
    data_SR = mtz.object(params.scaling.mtz_file)
    have_iso_ref = True
  except RuntimeError:
    data_SR = None
    have_iso_ref = False
  data_d0 = mtz.object(params.output.prefix+"_s0_"+params.scaling.algorithm+".mtz")
  data_d1 = mtz.object(params.output.prefix+"_s1_"+params.scaling.algorithm+".mtz")
  data_d2 = mtz.object(params.output.prefix+"_s2_"+params.scaling.algorithm+".mtz")

  uniform = []
  for idx,item in enumerate([data_SR,data_d0,data_d1,data_d2]):
    if not have_iso_ref and idx == 0:
      uniform.append(None)
      continue
    #item.show_summary()
    print("-------------------------------", file=output)
    for array in item.as_miller_arrays():
       this_label = array.info().label_string().lower()
       print(this_label, params.scaling.mtz_column_F)
       if this_label.find("fobs")>=0:
         print(this_label,array.observation_type(), file=output)
         uniform.append(array.as_intensity_array())
         break
       if this_label.find("iobs")>=0:
         """This test is added for the use case of unmerged anomalous data without
            an isomorphous reference (in other words, mark1).  Without this,
            the unmerged reflections are not picked up by the cc comparison.
            Indicates that this section probably has to be reanalyzed and redesigned.
         """
         print(this_label,array.observation_type(), file=output)
         uniform.append(array.as_intensity_array())
         break
       if this_label.find("imean")>=0:
         print(this_label,array.observation_type(), file=output)
         uniform.append(array.as_intensity_array())
         break
       if this_label.find(params.scaling.mtz_column_F)==0:
         print(this_label,array.observation_type(), file=output)
         uniform.append(array.as_intensity_array())
         break

  # If necesssary, generate Bijvoet mates for the isomorphous
  # reference.
  if have_iso_ref \
     and not params.merge_anomalous \
     and not uniform[0].anomalous_flag():
      uniform[0] = uniform[0].generate_bijvoet_mates()

  for x in [1,2,3]:
    # reindex the experimental data
    uniform[x] = uniform[x].change_basis(reindexing_op).map_to_asu()

  d_max_min = uniform[1].d_max_min()
  if have_iso_ref:
    sgi = uniform[0].space_group_info()
  else:
    sgi = params.target_space_group

  for x in [0,1,2,3]:
    if not have_iso_ref and x == 0:
      continue
    print("%6d indices:"%uniform[x].size(),{0:"Reference intensities",
                     1:"Merged structure factors",
                     2:"Semi-dataset 1",
                     3:"Semi-dataset 2"}[x], file=output)
    uniform[x] = uniform[x].customized_copy(
      crystal_symmetry = crystal.symmetry(unit_cell=uniform[1].unit_cell(),
                                          space_group_info=sgi),
      ).resolution_filter(d_min=uniform[1].d_min(),d_max=params.d_max,
      ).complete_array(d_min=uniform[1].d_min(),d_max=params.d_max).map_to_asu()
  print("%6d indices: An asymmetric unit in the resolution interval %.2f - %.2f Angstrom"%(
     uniform[1].size(),d_max_min[0],uniform[1].d_min()), file=output)

  if have_iso_ref:
    uniform[0] = uniform[0].common_set(uniform[1])
    assert len(uniform[0].indices()) == len(uniform[1].indices())
  uniform[2] = uniform[2].common_set(uniform[1])
  uniform[3] = uniform[3].common_set(uniform[1])
  print("-------------------------------", file=output)
  NBIN = params.output.n_bins
  for x in [0, 1, 2, 3]:
    if not have_iso_ref and x == 0:
      continue
    uniform[x].setup_binner(n_bins=NBIN)

  for x in range(len(uniform[1].indices())):
    if have_iso_ref:
      assert uniform[0].indices()[x] == uniform[1].indices()[x]
    assert uniform[1].indices()[x] == uniform[2].indices()[x]
    assert uniform[2].indices()[x] == uniform[3].indices()[x]

  cutoff = math.exp(params.scaling.log_cutoff or -100.)

  selected_uniform = []
  if have_iso_ref:
    if params.include_negatives:
      uniformA = (uniform[0].sigmas() > 0.).__and__(uniform[1].sigmas() > 0.)
    elif params.scaling.log_cutoff is None:
      uniformA = (uniform[0].data() > 0.).__and__(uniform[1].data() > 0.)
    else:
      uniformA = (uniform[0].data() > cutoff).__and__(uniform[1].data() > cutoff)
    for x in [0, 1]:
      selected_uniform.append(uniform[x].select(uniformA))
      selected_uniform[x].setup_binner(
        d_max=(params.d_max or 100000), d_min=params.d_min, n_bins=NBIN)

  else:
    selected_uniform = [None, None]

  if params.include_negatives:
      uniformB = (uniform[2].sigmas() > 0.).__and__(uniform[3].sigmas() > 0.)
  elif params.scaling.log_cutoff is None:
      uniformB = (uniform[2].data() > 0.).__and__(uniform[3].data() > 0.)
  else:
      uniformB = (uniform[2].data() > cutoff).__and__(uniform[3].data() > cutoff)

  for x in [2, 3]:
    selected_uniform.append(uniform[x].select(uniformB))
    selected_uniform[x].setup_binner(
      d_max=(params.d_max or 100000), d_min=params.d_min, n_bins=NBIN)

  return uniform, selected_uniform, have_iso_ref

def run_cc(params,reindexing_op,output):
  uniform, selected_uniform, have_iso_ref = load_cc_data(params, reindexing_op, output)
  NBIN = params.output.n_bins

  if have_iso_ref:
    slope, offset, corr_iso, N_iso = correlation(
      selected_uniform[1], selected_uniform[0], params.include_negatives)
    print("C.C. iso is %.1f%% on %d indices" % (
      100 * corr_iso, N_iso), file=output)

  slope,offset,corr_int,N_int = correlation(selected_uniform[2],selected_uniform[3], params.include_negatives)
  print("C.C. int is %.1f%% on %d indices"%(100.*corr_int, N_int), file=output)

  if have_iso_ref:
    binned_cc_ref, binned_cc_ref_N = binned_correlation(
      selected_uniform[1], selected_uniform[0], params.include_negatives)
    #binned_cc_ref.show(f=output)

    ref_scale = scale_factor(
      selected_uniform[1], selected_uniform[0],
      weights=flex.pow(selected_uniform[1].sigmas(), -2),
      use_binning=True)
    #ref_scale.show(f=output)

    ref_riso = r1_factor(
      selected_uniform[1], selected_uniform[0],
      scale_factor=ref_scale, use_binning=True)
    #ref_riso.show(f=output)

    ref_scale_all = scale_factor(
      selected_uniform[1], selected_uniform[0],
      weights=flex.pow(selected_uniform[1].sigmas(), -2))

    ref_riso_all = r1_factor(
      selected_uniform[1], selected_uniform[0],
      scale_factor=ref_scale_all)

  binned_cc_int,binned_cc_int_N = binned_correlation(
    selected_uniform[2], selected_uniform[3], params.include_negatives)
  #binned_cc_int.show(f=output)

  oe_scale = scale_factor(selected_uniform[2],selected_uniform[3],
    weights = flex.pow(selected_uniform[2].sigmas(),-2)
            + flex.pow(selected_uniform[3].sigmas(),-2),
    use_binning=True)
  #oe_scale.show(f=output)

  oe_rint = r1_factor(selected_uniform[2],selected_uniform[3],
                       scale_factor = oe_scale, use_binning=True)
  #oe_rint.show(f=output)

  oe_rsplit = r_split(selected_uniform[2], selected_uniform[3],
                      use_binning=True)

  oe_scale_all = scale_factor(selected_uniform[2],selected_uniform[3],
    weights = flex.pow(selected_uniform[2].sigmas(),-2)
            + flex.pow(selected_uniform[3].sigmas(),-2),)

  oe_rint_all = r1_factor(selected_uniform[2],selected_uniform[3],
                       scale_factor = oe_scale_all)
  oe_rsplit_all = r_split(selected_uniform[2], selected_uniform[3])
  if have_iso_ref:
    print("R factors Riso = %.1f%%, Rint = %.1f%%"%(100.*ref_riso_all, 100.*oe_rint_all), file=output)
  else:
    print("R factor Rint = %.1f%%"%(100.*oe_rint_all), file=output)

  split_sigma_data = split_sigma_test(selected_uniform[2],selected_uniform[3],
                                      scale=oe_scale,use_binning=True,show_plot=False)
  split_sigma_data_all = split_sigma_test(selected_uniform[2],selected_uniform[3],
                                          scale=oe_scale_all,use_binning=False,show_plot=False)

  print(file=output)
  if reindexing_op == "h,k,l":
    print("Table of Scaling Results:", file=output)
  else:
    print("Table of Scaling Results Reindexing as %s:"%reindexing_op, file=output)

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

  print(file=output)
  print(table_utils.format(table_data,has_header=2,justify='center',delim=" "), file=output)
  print("""CCint is the CC-1/2 defined by Diederichs; correlation between odd/even images.
  Similarly, Scale int and R int are the scaling factor and scaling R factor between odd/even images.
  "iso" columns compare the whole XFEL dataset to the isomorphous reference.""", file=output)

  print("""Niso: result vs. reference common set""", end=' ', file=output)
  if params.include_negatives:
    print("""including negative merged intensities (set by phil parameter).""", file=output)
  elif params.scaling.log_cutoff is None:
    print(file=output)
  else:
    print("""with intensites < %7.2g filtered out (controlled by
    scaling.log_cutoff phil parameter set to %5.1f)"""%(math.exp(params.scaling.log_cutoff),
    params.scaling.log_cutoff), file=output)

  if have_iso_ref:
    assert N_iso == flex.sum(flex.double([x for x in binned_cc_ref_N.data if x is not None]))
  assert N_int == flex.sum(flex.double([x for x in binned_cc_int_N.data if x is not None]))

  if params.scaling.show_plots:
    from matplotlib import pyplot as plt
    plt.plot(flex.log(selected_uniform[-2].data()),
             flex.log(selected_uniform[-1].data()), 'r.')
    plt.show()
    if have_iso_ref:
      plt.plot(flex.log(selected_uniform[0].data()),
               flex.log(selected_uniform[1].data()), 'r.')
      plt.show()
  print(file=output)
