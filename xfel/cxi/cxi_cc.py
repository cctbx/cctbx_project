import math
from libtbx.str_utils import format_value
from iotbx import mtz
from cctbx import crystal
from scitbx.array_family import flex
from cctbx.miller import binned_data

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

def binned_correlation(self,other):
    results = []
    bin_count = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      if sel.count(True)==0:
        results.append(0.)
        bin_count.append(0.)
        continue
      result_tuple = correlation(self.select(sel),other.select(sel))
      results.append(result_tuple[2])
      bin_count.append(result_tuple[3])
      # plots for debugging
      #from matplotlib import pyplot as plt
      #plt.plot(flex.log(self.select(sel).data()),flex.log(other.select(sel).data()),"b.")
      #plt.show()

    return binned_data(binner=self.binner(), data=results, data_fmt="%7.4f"),\
           binned_data(binner=self.binner(), data=bin_count, data_fmt="%7d")

def correlation(self,other):
    N = 0
    sum_xx = 0
    sum_xy = 0
    sum_yy = 0
    sum_x = 0
    sum_y = 0
    for idx in xrange(self.indices().size()):

      assert self.indices()[idx]==other.indices()[idx]
      I_r = other.data()[idx]
      I_o = self.data()[idx]
      #assert I_r >= 0. or I_o >= 0.
      #why does this go from 81% to 33% when uniform selection is made?

      if I_r < 0. or I_o < 0.: continue
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
    slope = (N * sum_xy - sum_x * sum_y) / (N * sum_xx - sum_x**2)
    offset = (sum_xx * sum_y - sum_x * sum_xy) / (N * sum_xx - sum_x**2)
    corr  = (N * sum_xy - sum_x * sum_y) / (math.sqrt(N * sum_xx - sum_x**2) *
             math.sqrt(N * sum_yy - sum_y**2))
    return slope,offset,corr,N

def run_cc(params,output):
  data_SR = mtz.object(params.scaling.mtz_file)
  data_d0 = mtz.object(params.output.prefix+"_s0_"+params.scaling.algorithm+".mtz")
  data_d1 = mtz.object(params.output.prefix+"_s1_"+params.scaling.algorithm+".mtz")
  data_d2 = mtz.object(params.output.prefix+"_s2_"+params.scaling.algorithm+".mtz")

  uniform = []
  for idx,item in enumerate([data_SR,data_d0,data_d1,data_d2]):
    #item.show_summary()
    print >>output, "-------------------------------"
    for array in item.as_miller_arrays():
       this_label = array.info().label_string().lower()
       if this_label.find("fobs")>=0:
         print >>output, this_label,array.observation_type()
         uniform.append(array.as_intensity_array())
       if this_label.find("imean")>=0:
         print >>output, this_label,array.observation_type()
         uniform.append(array.as_intensity_array())
  reserve_highres = params.d_min
  for x in [0,1,2,3]:
    print >>output, uniform[x].size()
    uniform[x] = uniform[x].customized_copy(
      crystal_symmetry = crystal.symmetry(unit_cell=uniform[1].unit_cell(),
                                          space_group_info=uniform[0].space_group_info()),
      ).resolution_filter(d_min=uniform[1].d_min()
      ).complete_array(d_min=uniform[1].d_min()).map_to_asu()
    print >>output, uniform[x].size()

  uniform[2] = uniform[2].common_set(uniform[1])
  uniform[0] = uniform[0].common_set(uniform[1])
  uniform[3] = uniform[3].common_set(uniform[1])
  print >>output, "-------------------------------"
  NBIN = params.output.n_bins
  uniform[0].setup_binner(n_bins=NBIN)
  for x in [1,2,3]: uniform[x].setup_binner(n_bins=NBIN)
  for x in xrange(len(uniform[0].indices())):
    assert uniform[1].indices()[x] == uniform[0].indices()[x]
    assert uniform[2].indices()[x] == uniform[0].indices()[x]
    assert uniform[3].indices()[x] == uniform[2].indices()[x]
    continue
    print >>output, "%15s %15s %15s %15s %10.0f %10.0f %10.0f %10.0f"%(
      uniform[0].indices()[x], uniform[1].indices()[x],
      uniform[2].indices()[x], uniform[3].indices()[x],
      uniform[0].data()[x], uniform[1].data()[x],
      uniform[2].data()[x], uniform[3].data()[x],)


  ref_scale = scale_factor(uniform[1],uniform[0])
  oe_scale = scale_factor(uniform[2],uniform[3])

  cutoff = math.exp(params.scaling.log_cutoff or -100.)
  uniformA = (uniform[0].data()>=0.).__and__(uniform[1].data()>=0.).__and__(uniform[0].data() > cutoff).__and__(uniform[1].data() > cutoff)

  slope,offset,corr_iso,N_iso = correlation(uniform[1].select(uniformA),uniform[0].select(uniformA))
  print >>output,"C.C. iso is %.1f%% on %d indices"%(100.*corr_iso, N_iso)

  uniformB = (uniform[2].data()>=0.).__and__(uniform[3].data()>=0.).__and__(uniform[2].data()>cutoff).__and__(uniform[3].data() > cutoff)

  slope,offset,corr_int,N_int = correlation(uniform[2].select(uniformB),uniform[3].select(uniformB))
  print >>output, "C.C. int is %.1f%% on %d indices"%(100.*corr_int, N_int)

  selected_uniform = []
  for x in [0,1]: selected_uniform.append(uniform[x].select(uniformA))
  for x in [2,3]: selected_uniform.append(uniform[x].select(uniformB))
  for x in [0,1,2,3]: selected_uniform[x].setup_binner(d_max=100000, d_min=reserve_highres, n_bins=NBIN)

  binned_cc_ref,binned_cc_ref_N = binned_correlation(selected_uniform[1],selected_uniform[0])
  #binned_cc_ref.show(f=output)

  binned_cc_int,binned_cc_int_N = binned_correlation(selected_uniform[2],selected_uniform[3])
  #binned_cc_int.show(f=output)

  ref_scale = scale_factor(selected_uniform[1],selected_uniform[0],
    weights = flex.pow(selected_uniform[1].sigmas(),-2),
    use_binning=True)
  oe_scale = scale_factor(selected_uniform[2],selected_uniform[3],
    weights = flex.pow(selected_uniform[2].sigmas(),-2)
            + flex.pow(selected_uniform[3].sigmas(),-2),
    use_binning=True)


  #ref_scale.show(f=output)
  #oe_scale.show(f=output)

  ref_riso = r1_factor(selected_uniform[1],selected_uniform[0],
                       scale_factor = ref_scale, use_binning=True)
  #ref_riso.show(f=output)
  oe_rint = r1_factor(selected_uniform[2],selected_uniform[3],
                       scale_factor = oe_scale, use_binning=True)
  #oe_rint.show(f=output)


  ref_scale_all = scale_factor(selected_uniform[1],selected_uniform[0],
    weights = flex.pow(selected_uniform[1].sigmas(),-2),)
  oe_scale_all = scale_factor(selected_uniform[2],selected_uniform[3],
    weights = flex.pow(selected_uniform[2].sigmas(),-2)
            + flex.pow(selected_uniform[3].sigmas(),-2),)
  print >>output, "Scale factors iso/int", ref_scale_all,oe_scale_all

  ref_riso_all = r1_factor(selected_uniform[1],selected_uniform[0],
                       scale_factor = ref_scale_all)
  oe_rint_all = r1_factor(selected_uniform[2],selected_uniform[3],
                       scale_factor = oe_scale_all)
  print >>output, "R factors Riso = %.1f%%, Rint = %.1f%%"%(100.*ref_riso_all, 100.*oe_rint_all)

  print >>output
  print >> output, "Table of Scaling Results:"
  from libtbx import table_utils
  table_header = ["","","","CC","","CC","","R","R","Scale","Scale"]
  table_header2 = ["Bin","Resolution Range","Completeness","int","N","iso","N","int","iso","int","iso"]
  table_data = []
  table_data.append(table_header)
  table_data.append(table_header2)

  items = binned_cc_ref.binner.range_used()

  for bin in items:
    table_row = []
    table_row.append("%3d"%bin)
    table_row.append("%-13s"%binned_cc_ref.binner.bin_legend(i_bin=bin,show_bin_number=False,show_bin_range=False,
                                                 show_d_range=True, show_counts=False))
    table_row.append("%13s"%binned_cc_ref.binner.bin_legend(i_bin=bin,show_bin_number=False,show_bin_range=False,
                                                 show_d_range=False, show_counts=True))
    table_row.append("%.1f%%"%(100.*binned_cc_int.data[bin]))
    table_row.append("%7d"%(binned_cc_int_N.data[bin]))
    table_row.append("%.1f%%"%(100.*binned_cc_ref.data[bin]))
    table_row.append("%7d"%(binned_cc_ref_N.data[bin]))
    if oe_rint.data[bin] is not None:
      table_row.append("%.1f%%"%(100.*oe_rint.data[bin]))
    else:
      table_row.append("--")
    if ref_riso.data[bin] is not None:
      table_row.append("%.1f%%"%(100.*ref_riso.data[bin]))
    else:
      table_row.append("--")
    if oe_scale.data[bin] is not None:
      table_row.append("%.3f"%oe_scale.data[bin])
    else:
      table_row.append("--")
    if ref_scale.data[bin] is not None:
      table_row.append("%.3f"%ref_scale.data[bin])
    else:
      table_row.append("--")
    table_data.append(table_row)
  table_data.append([""]*len(table_header))
  table_data.append(  [
      format_value("%3s",   "All"),
      format_value("%-13s", "                 "),
      format_value("%13s",  ""),
      format_value("%.1f%%", 100.*corr_int),
      format_value("%7d", N_int),
      format_value("%.1f%%", 100.*corr_iso),
      format_value("%7d", N_iso),
      format_value("%.1f%%", 100.*oe_rint_all),
      format_value("%.1f%%", 100.*ref_riso_all),
      format_value("%.3f", oe_scale_all),
      format_value("%.3f", ref_scale_all),
  ])

  print >>output
  print >>output,table_utils.format(table_data,has_header=2,justify='center',delim=" ")
  print >>output,"""CCint is the CC-1/2 defined by Diedrichs; correlation between odd/even images.
  Similarly, Scale int and R int are the scaling factor and scaling R factor between odd/even images.
  "iso" columns compare the whole XFEL dataset to the isomorphous reference."""


  if params.scaling.show_plots:
    from matplotlib import pyplot as plt
    plt.plot(flex.log(uniform[2].select(uniformB).data()),flex.log(uniform[3].select(uniformB).data()),"r.")
    plt.show()
    plt.plot(flex.log(uniform[0].select(uniformA).data()),flex.log(uniform[1].select(uniformA).data()),"r.")
    plt.show()
  print >>output
