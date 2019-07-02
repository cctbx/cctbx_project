from __future__ import absolute_import, division, print_function
from six.moves import range
from cctbx.array_family import flex

class show_observations:
  def __init__(self,obs,unobstructed,params,out=None, n_bins=12):
    if out==None:
      import sys
      out = sys.stdout
    from libtbx.str_utils import format_value
    self.params = params
    self.out = out
    self.obs = obs

    obs.setup_binner(n_bins = n_bins)
    attenuated = obs.select(~unobstructed)
    unattenuated = obs.select(unobstructed)

    attenuated.use_binning_of(obs)
    unattenuated.use_binning_of(obs)

    self.result = []
    counts_given = obs.binner().counts_given()
    counts_complete = obs.binner().counts_complete()

    counts_given_attenuated = attenuated.binner().counts_given()
    counts_given_unattenuated = unattenuated.binner().counts_given()

    for i_bin in obs.binner().range_used():
      sel_w = obs.binner().selection(i_bin)
      sel_fo_all = obs.select(sel_w)
      d_max_,d_min_ = sel_fo_all.d_max_min()
      d_range = obs.binner().bin_legend(
        i_bin=i_bin, show_bin_number=False, show_counts=True)
      sel_data = obs.select(sel_w).data()
      sel_sig = obs.select(sel_w).sigmas()

      sel_unatten_w = unattenuated.binner().selection(i_bin)
      sel_unatten_data = unattenuated.select(sel_unatten_w).data()
      sel_unatten_sig = unattenuated.select(sel_unatten_w).sigmas()

      sel_atten_w = attenuated.binner().selection(i_bin)
      sel_atten_data = attenuated.select(sel_atten_w).data()
      sel_atten_sig = attenuated.select(sel_atten_w).sigmas()

      if len(sel_unatten_data)>0:
        unatten_mean_I = flex.mean(sel_unatten_data)
        unatten_mean_I_sigI  = flex.mean(sel_unatten_data/sel_unatten_sig)
      else:
        unatten_mean_I = 0
        unatten_mean_I_sigI  = 0

      if len(sel_atten_data)>0:
        atten_mean_I = flex.mean(sel_atten_data)
        atten_mean_I_sigI  = flex.mean(sel_atten_data/sel_atten_sig)
      else:
        atten_mean_I = 0
        atten_mean_I_sigI  = 0

      if(sel_data.size() > 0):
        bin = resolution_bin(
          i_bin        = i_bin,
          d_range      = d_range,
          mean_I       = flex.mean(sel_data),
          n_work       = sel_data.size(),
          mean_I_sigI  = flex.mean(sel_data/sel_sig),
          d_max_min    = (d_max_, d_min_),
          completeness = (counts_given[i_bin], counts_complete[i_bin]),
          given_unatten = counts_given_unattenuated[i_bin],
          unatten_mean_I       = unatten_mean_I,
          unatten_mean_I_sigI  = unatten_mean_I_sigI,
          given_atten  = counts_given_attenuated[i_bin],
          atten_mean_I       = atten_mean_I,
          atten_mean_I_sigI  = atten_mean_I_sigI,
        )
        self.result.append(bin)
    self.set_limits(unobstructed)
    print("\n Bin  Resolution Range  Compl.         <I>     <I/sig(I)> Unobstructed <I>     <I/sig(I)> Obstructed <I>     <I/sig(I)>", file=out)
    for bin in self.result:
      fmt = " %s %s    %s  %s%s   %s   %s  %s%s   %s  %s   %s%s"
      print(fmt%(
        format_value("%3d",   bin.i_bin),
        format_value("%-17s", bin.d_range),
        format_value("%8.1f", bin.mean_I),
        format_value("%8.2f", bin.mean_I_sigI),
        format_value("%1s",   getattr(bin,"limit"," ")),
        format_value("%6d",   bin.given_unatten),
        format_value("%8.1f", bin.unatten_mean_I),
        format_value("%8.2f", bin.unatten_mean_I_sigI),
        format_value("%1s",   getattr(bin,"unatten_limit"," ")),
        format_value("%6d",   bin.given_atten),
        format_value("%8.1f", bin.atten_mean_I),
        format_value("%8.2f", bin.atten_mean_I_sigI),
        format_value("%1s",   getattr(bin,"atten_limit"," ")),
        ), file=out)

  def set_limits(self, unobstructed):
    acceptable_resolution_bins = [
      bin.mean_I_sigI > self.params.significance_filter.sigma for bin in self.result]
    acceptable_nested_bin_sequences = [i for i in range(len(acceptable_resolution_bins))
                                       if False not in acceptable_resolution_bins[:i+1]]

    N_acceptable_bins = max(acceptable_nested_bin_sequences)
    self.result[N_acceptable_bins].limit="*"

    unatten_acceptable = N_acceptable_bins
    for x in range(N_acceptable_bins,len(self.result)):
      if self.result[x].unatten_mean_I_sigI > self.params.significance_filter.sigma:
        unatten_acceptable = x
      else:
        break

    self.result[unatten_acceptable].unatten_limit="*"

    atten_acceptable = N_acceptable_bins
    for x in range(N_acceptable_bins,1,-1):
      if self.result[x].atten_mean_I_sigI < self.params.significance_filter.sigma:
        atten_acceptable = x - 1

    self.result[atten_acceptable].atten_limit="*"

    # Now compute the appropriate selections
    unattenuated_res_limit = float(self.result[unatten_acceptable].d_range.split()[2])
    attenuated_res_limit = float(self.result[atten_acceptable].d_range.split()[2])
    print("New combination resolution filter at %7.2f and %7.2f"%(unattenuated_res_limit,
    attenuated_res_limit), file=self.out)

    unattenuated_res_selection = self.obs.resolution_filter_selection(
            d_min=unattenuated_res_limit)
    attenuated_res_selection = self.obs.resolution_filter_selection(
            d_min=attenuated_res_limit)

    self.master_selection = (unobstructed & unattenuated_res_selection) | (~unobstructed & attenuated_res_selection)


class resolution_bin(object):
  def __init__(self,
               **kwargs):
    for key in kwargs.keys():
      assert not hasattr(self.__dict__, key)
      self.__dict__[key] = kwargs[key]
    #similar funtionality: from libtbx import adopt_init_args
    # XXX not sure if this introduces a memory leak.  Keep an eye on it.
