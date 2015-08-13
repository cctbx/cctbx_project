from __future__ import division
import sys
from cctbx.array_family import flex
from libtbx.str_utils import format_value
from iotbx import data_plots
from libtbx import adopt_init_args

def show_overall_observations(Fit_I,Fit_I_stddev,model_I,I_visited,ordered,frames,sim,
  out = None,title = None,work_params = None):

  # at minimum we need the Miller set of merged intensities and sigmas
  model_subset = ordered[0:len(Fit_I)]
  uc = ordered.unit_cell()
  model_subset.show_summary()

  d_min = flex.min(uc.d(model_subset.indices())) # extent of data
  n_bins = min( len(Fit_I)//20, 15 )

  #obs, redundancy, summed_wt_I, summed_weight, ISIGI, n_bins=15, out=None, title=None, work_params=None):
  if out is None:
    out = sys.stdout
  model_subset.setup_binner(d_max=100000, d_min=d_min, n_bins=n_bins)
  result = []
  multiplicity = flex.int(len(Fit_I))
  for iraw in xrange(len(sim.miller)):
    multiplicity[sim.miller[iraw]] += 1
  cumulative_unique = 0
  cumulative_meas   = 0
  cumulative_theor  = 0
  cumulative_In     = 0
  cumulative_I      = 0.0
  cumulative_Isigma = 0.0

  for i_bin in model_subset.binner().range_used():
    sel_w = model_subset.binner().selection(i_bin)
    sel_fo_all = model_subset.select(sel_w)
    d_range = model_subset.binner().bin_legend(
      i_bin=i_bin, show_bin_number=False, show_counts=False)

    sel_multiplicity = multiplicity.select(sel_w)
    sel_absent = sel_multiplicity.count(0)
    n_present = sel_multiplicity.size() - sel_absent
    sel_complete_tag = "[%d/%d]" % (n_present, sel_multiplicity.size())
    sel_measurements = flex.sum(sel_multiplicity)

    # Alternatively, redundancy (or multiplicity) is calculated as the
    # average number of observations for the observed
    # reflections--missing reflections do not affect the redundancy
    # adversely, and the reported value becomes
    # completeness-independent.
    val_multiplicity_obs = 0
    if n_present > 0:
      val_multiplicity_obs = flex.sum(sel_multiplicity) / n_present


    # Per-bin sum of I and I/sig(I).  For any reflection, the weight
    # of the merged intensity must be positive for this to make sense.
    sel_o = (sel_w & (Fit_I_stddev > 0.))
    selected_intensity = Fit_I.select(sel_o)
    selected_stddev    = Fit_I_stddev.select(sel_o)
    I_sum = flex.sum(selected_intensity)
    assert selected_stddev.count(0.) == 0
    I_sigI_sum = flex.sum(selected_intensity / selected_stddev)
    I_n = sel_o.count(True)
    '''
    # Per-bin sum of I and I/sig(I) for each observation.
    # R-merge statistics have been removed because
    #  >> R-merge is defined on whole structure factor intensities, either
    #     full observations or summed partials from the rotation method.
    #     For XFEL data all reflections are assumed to be partial; no
    #     method exists now to convert partials to fulls.
    if work_params.plot_single_index_histograms: import numpy as np
    for i in obs.binner().array_indices(i_bin) :
      index = obs.indices()[i]
      if (index in ISIGI) :
        # Compute m, the "merged" intensity, as the average intensity
        # of all observations of the reflection with the given index.
        N = 0
        m = 0
        for t in ISIGI[index] :
          N += 1
          m += t[0]
        if work_params is not None and \
           (work_params.plot_single_index_histograms and \
            N >= 30 and \
            work_params.data_subset == 0):
          print "Miller %20s n-obs=%4d  sum-I=%10.0f"%(index, N, m)
          plot_n_bins = N//10
          hist,bins = np.histogram([t[0] for t in ISIGI[index]],bins=25)
          width = 0.7*(bins[1]-bins[0])
          center = (bins[:-1]+bins[1:])/2
          import matplotlib.pyplot as plt
          plt.bar(center, hist, align="center", width=width)
          plt.show()
    '''
    if sel_measurements > 0:
      mean_I = mean_I_sigI = 0
      if I_n > 0:
        mean_I = I_sum / I_n
        mean_I_sigI = I_sigI_sum / I_n
      bin = resolution_bin(
        i_bin=i_bin,
        d_range=d_range,
        d_min=model_subset.binner().bin_d_min(i_bin),
        redundancy_asu=flex.mean(sel_multiplicity.as_double()),
        redundancy_obs=val_multiplicity_obs,
        complete_tag=sel_complete_tag,
        completeness=n_present / sel_multiplicity.size(),
        measurements=sel_measurements,
        mean_I=mean_I,
        mean_I_sigI=mean_I_sigI
        )
      result.append(bin)
    cumulative_unique += n_present
    cumulative_meas   += sel_measurements
    cumulative_theor  += sel_multiplicity.size()
    cumulative_In     += I_n
    cumulative_I      += I_sum
    cumulative_Isigma += I_sigI_sum

  if (title is not None) :
    print >> out, title
  from libtbx import table_utils
  table_header = ["","","","<asu","<obs","","","",""]
  table_header2 = ["Bin","Resolution Range","Completeness","multi>","multi>","n_meas","<I>","<I/sig(I)>"]
  table_data = []
  table_data.append(table_header)
  table_data.append(table_header2)
  for bin in result:
    table_row = []
    table_row.append("%3d" % bin.i_bin)
    table_row.append("%-13s" % bin.d_range)
    table_row.append("%13s" % bin.complete_tag)
    table_row.append("%6.2f" % bin.redundancy_asu)
    table_row.append("%6.2f" % bin.redundancy_obs)
    table_row.append("%6d" % bin.measurements)
    table_row.append("%8.3f" % bin.mean_I)
    table_row.append("%8.3f" % bin.mean_I_sigI)
    table_data.append(table_row)
  table_data.append([""]*len(table_header))
  table_data.append(  [
      format_value("%3s",   "All"),
      format_value("%-13s", "                 "),
      format_value("%13s",  "[%d/%d]"%(cumulative_unique,cumulative_theor)),
      format_value("%6.2f", cumulative_meas/cumulative_theor),
      format_value("%6.2f", cumulative_meas/cumulative_unique),
      format_value("%6d",   cumulative_meas),
      format_value("%8.3f", cumulative_I/cumulative_In),
      format_value("%8.3f", cumulative_Isigma/cumulative_In),
  ])

  print
  print >>out,table_utils.format(table_data,has_header=2,justify='center',delim=" ")

  # XXX generate table object for displaying plots
  if (title is None) :
    title = "Data statistics by resolution"
  table = data_plots.table_data(
    title=title,
    x_is_inverse_d_min=True,
    force_exact_x_labels=True)
  table.add_column(
    column=[1 / bin.d_min**2 for bin in result],
    column_name="d_min",
    column_label="Resolution")
  table.add_column(
    column=[bin.redundancy_asu for bin in result],
    column_name="redundancy",
    column_label="Redundancy")
  table.add_column(
    column=[bin.completeness for bin in result],
    column_name="completeness",
    column_label="Completeness")
  table.add_column(
    column=[bin.mean_I_sigI for bin in result],
    column_name="mean_i_over_sigI",
    column_label="<I/sig(I)>")
  table.add_graph(
    name="Redundancy vs. resolution",
    type="GRAPH",
    columns=[0,1])
  table.add_graph(
    name="Completeness vs. resolution",
    type="GRAPH",
    columns=[0,2])
  table.add_graph(
    name="<I/sig(I)> vs. resolution",
    type="GRAPH",
    columns=[0,3])
  return table,n_bins,d_min

class resolution_bin(object):
  def __init__(self,
               i_bin=None,
               d_range=None,
               d_min=None,
               redundancy_asu=None,
               redundancy_obs=None,
               absent=None,
               complete_tag=None,
               completeness=None,
               measurements=None,
               mean_I=None,
               mean_I_sigI=None,
               sigmaa=None):
    adopt_init_args(self, locals())
