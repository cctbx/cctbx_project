from __future__ import division
from libtbx.str_utils import make_sub_header, format_value
from libtbx.utils import Sorry, null_out
from libtbx import group_args, Auto, slots_getstate_setstate
from math import sqrt
import cStringIO
import sys

citations_str = """\
  Diederichs K & Karplus PA (1997) Nature Structural Biology 4:269-275
    (with erratum in: Nat Struct Biol 1997 Jul;4(7):592)
  Weiss MS (2001) J Appl Cryst 34:130-135.
  Karplus PA & Diederichs K (2012) Science 336:1030-3."""

sigma_filtering_phil_str = """
sigma_filtering = *auto xds scala scalepack
  .type = choice
  .short_caption = Sigma(I) filtering convention
  .help = Determines how data are filtered by SigmaI and I/SigmaI.  XDS \
    discards reflections whose intensity after merging is less than -3*sigma, \
    Scalepack uses the same cutoff before merging, and SCALA does not do any \
    filtering.  Reflections with negative SigmaI will always be discarded.
"""

merging_params_str = """
high_resolution = None
  .type = float
  .input_size = 64
low_resolution = None
  .type = float
  .input_size = 64
n_bins = 10
  .type = int
  .short_caption = Number of resolution bins
  .input_size = 64
  .style = spinner
anomalous = False
  .type = bool
  .short_caption = Keep anomalous pairs separate in merging statistics
%s
""" % sigma_filtering_phil_str

class model_based_arrays (object) :
  """
  Container for observed and calculated intensities, along with the selections
  for work and free sets; these should be provided by mmtbx.f_model.  It is
  assumed (or hoped) that the resolution range of these arrays will be
  the same as that of the unmerged data, but the current implementation does
  not force this.
  """
  def __init__ (self, f_obs, i_obs, i_calc, work_sel, free_sel) :
    assert (i_obs.data().size() == i_calc.data().size() ==
            work_sel.data().size() == free_sel.data().size())
    assert (len(f_obs.data()) <= len(i_obs.data()))
    self.f_obs = f_obs
    self.i_obs = i_obs.common_set(other=self.f_obs)
    self.i_calc = i_calc.common_set(other=self.f_obs)
    self.work_sel = work_sel.common_set(other=self.f_obs)
    self.free_sel = free_sel.common_set(other=self.f_obs)

  def cc_work_and_free (self, other) :
    """
    Given a unique array of arbitrary resolution range, extract the equivalent
    reflections from the observed and calculated intensities, and calculate
    CC and R-factor for work and free sets.  Currently, these statistics will
    be None if there are no matching reflections.
    """
    assert (self.i_obs.is_similar_symmetry(other))
    i_obs_sel = self.i_obs.common_set(other=other)
    f_obs_sel = self.f_obs.common_set(other=other)
    i_calc_sel = self.i_calc.common_set(other=other)
    work_sel = self.work_sel.common_set(other=other)
    free_sel = self.free_sel.common_set(other=other)
    if (len(i_obs_sel.data()) == 0) : # XXX should this raise an error?
      return [None] * 4
    i_obs_work = i_obs_sel.select(work_sel.data())
    i_calc_work = i_calc_sel.select(work_sel.data())
    i_obs_free = i_obs_sel.select(free_sel.data())
    i_calc_free = i_calc_sel.select(free_sel.data())
    f_obs_work = f_obs_sel.select(work_sel.data())
    f_obs_free = f_obs_sel.select(free_sel.data())
    if (len(f_obs_work.data()) > 0) and (len(f_obs_free.data()) > 0) :
      from scitbx.array_family import flex
      cc_work = flex.linear_correlation(i_obs_work.data(),
        i_calc_work.data()).coefficient()
      cc_free = flex.linear_correlation(i_obs_free.data(),
        i_calc_free.data()).coefficient()
      r_work = f_obs_work.r1_factor(i_calc_work.f_sq_as_f())
      r_free = f_obs_free.r1_factor(i_calc_free.f_sq_as_f())
      return cc_work, cc_free, r_work, r_free
    return [None] * 4

def get_filtering_convention (i_obs, sigma_filtering=Auto) :
  info = i_obs.info()
  if (sigma_filtering in [Auto, "auto"]) :
    if (info.source_type == "xds_ascii") :
      sigma_filtering = "xds"
    elif (info.source_type == "ccp4_mtz") :
      sigma_filtering = "scala"
    elif (info.source_type == "scalepack_no_merge_original_index") :
      sigma_filtering = "scalepack"
    else : # XXX default to the most conservative method
      sigma_filtering = "scala"
  return sigma_filtering

class filter_intensities_by_sigma (object) :
  """
  Wrapper for filtering intensities based on one of several different
  conventions:

    - in XDS, reflections where I < -3*sigmaI after merging are deleted from
      both the merged and unmerged arrays
    - in Scalepack, the filtering is done before merging
    - SCALA and AIMLESS do not do any filtering

  note that ctruncate and cctbx.french_wilson (any others?) do their own
  filtering, e.g. discarding I < -4*sigma in cctbx.french_wilson.
  """
  def __init__ (self, array, sigma_filtering=Auto) :
    sigma_filtering = get_filtering_convention(array, sigma_filtering)
    assert (sigma_filtering in ["scala","scalepack","xds", None])
    self.n_rejected_before_merge = self.n_rejected_after_merge = 0
    merge = array.merge_equivalents(use_internal_variance=False)
    array_merged = merge.array()
    reject_sel = None
    self.observed_criterion_sigma_I = None
    if (sigma_filtering == "xds") :
      self.observed_criterion_sigma_I = -3
      reject_sel = (array_merged.data() < -3*array_merged.sigmas())
      self.n_rejected_after_merge = reject_sel.count(True)
      bad_data = array_merged.select(reject_sel)
      array = array.delete_indices(other=bad_data)
      # and merge again...
      merge = array.merge_equivalents(use_internal_variance=False)
      array_merged = merge.array()
    elif (sigma_filtering == "scalepack") :
      self.observed_criterion_sigma_I = -3
      reject_sel = (array.data() < -3* array.sigmas())
      self.n_rejected_before_merge = reject_sel.count(True)
      array = array.select(~reject_sel)
      merge = array.merge_equivalents(use_internal_variance=False)
      array_merged = merge.array()
    elif (sigma_filtering == "scala") or (sigma_filtering is None) :
      pass
    else :
      raise ValueError("Unrecognized sigmaI filtering convention '%s'." %
        sigma_filtering)
    self.array = array
    self.merge = merge
    self.array_merged = array_merged

class merging_stats (object) :
  """
  Calculate standard merging statistics for (scaled) unmerged data.  Usually
  these statistics will consider I(+) and I(-) as observations of the same
  reflection, but these can be kept separate instead if desired.

  Reflections with negative sigmas will be discarded, and depending on the
  program we're trying to mimic, excessively negative intensities.
  """
  def __init__ (self,
      array,
      model_arrays=None,
      anomalous=False,
      debug=None,
      sigma_filtering="scala") :
    import cctbx.miller
    from scitbx.array_family import flex
    assert (array.sigmas() is not None)
    array = array.eliminate_sys_absent()
    array = array.customized_copy(anomalous_flag=anomalous).map_to_asu()
    non_negative_sel = array.sigmas() >= 0
    self.n_neg_sigmas = non_negative_sel.count(False)
    array = array.select(non_negative_sel)
    filter = filter_intensities_by_sigma(
      array=array,
      sigma_filtering=sigma_filtering)
    self.d_max, self.d_min = array.d_max_min()
    self.observed_criterion_sigma_I = filter.observed_criterion_sigma_I
    array = filter.array
    merge = filter.merge
    array_merged = filter.array_merged
    self.n_rejected_before_merge = filter.n_rejected_before_merge
    self.n_rejected_after_merge = filter.n_rejected_after_merge
    self.n_obs = array.indices().size()
    self.n_uniq = array_merged.indices().size()
    complete_set = array_merged.complete_set().resolution_filter(
      d_min=self.d_min, d_max=self.d_max)
    n_expected = len(complete_set.indices())
    if (n_expected == 0) :
      raise RuntimeError(("No reflections within specified resolution range "+
        "(%g - %g)") % (self.d_max, self.d_min))
    self.completeness = min(self.n_uniq / n_expected, 1.)
    self.anom_completeness = None
    # TODO also calculate when anomalous=False, since it is customary to
    # calculate merging statistics with F+ and F- treated as redundant
    # observations even when we're going to keep them separate.
    if (anomalous) :
      self.anom_completeness = array_merged.anomalous_completeness()
    redundancies = merge.redundancies().data()
    self.redundancies = {}
    for x in sorted(set(redundancies)) :
      self.redundancies[x] = redundancies.count(x)
    self.mean_redundancy = flex.mean(redundancies.as_double())
    self.i_mean = flex.mean(array_merged.data())
    self.sigi_mean = flex.mean(array_merged.sigmas())
    nonzero_array = array_merged.select(array_merged.sigmas() > 0)
    i_over_sigma = nonzero_array.data() / nonzero_array.sigmas()
    self.i_over_sigma_mean = flex.mean(i_over_sigma)
    self.i_mean_over_sigi_mean = self.i_mean/self.sigi_mean
    self.r_merge = merge.r_merge()
    self.r_meas = merge.r_meas()
    self.r_pim = merge.r_pim()
    # XXX Pure-Python reference implementation
    if (debug) :
      from libtbx.test_utils import approx_equal
      r_merge_num = r_meas_num = r_pim_num = r_merge_den = 0
      indices = array_merged.indices()
      data = array_merged.data()
      for hkl, i_mean in zip(indices, data) :
        sele = (array.indices() == hkl)
        hkl_array = array.select(sele)
        n_hkl = hkl_array.indices().size()
        if (n_hkl > 1) :
          sum_num = sum_den = 0
          for i_obs in hkl_array.data() :
            sum_num += abs(i_obs - i_mean)
            sum_den += i_obs
          r_merge_num += sum_num
          r_meas_num += sqrt(n_hkl/(n_hkl-1.)) * sum_num
          r_pim_num += sqrt(1./(n_hkl-1)) * sum_num
          r_merge_den += sum_den
      assert (approx_equal(self.r_merge, r_merge_num / r_merge_den))
      assert (approx_equal(self.r_meas, r_meas_num / r_merge_den))
      assert (approx_equal(self.r_pim, r_pim_num / r_merge_den))
    #---
    self.cc_one_half = cctbx.miller.compute_cc_one_half(
      unmerged=array)
    if (self.cc_one_half == 0) :
      self.cc_star = 0
    else :
      mult = 1.
      if (self.cc_one_half < 0) :
        mult = -1.
      self.cc_star = mult * sqrt((2*abs(self.cc_one_half)) /
                                 (1 + self.cc_one_half))
    self.cc_work = self.cc_free = self.r_work = self.r_free = None
    if (model_arrays is not None) :
      self.cc_work, self.cc_free, self.r_work, self.r_free = \
        model_arrays.cc_work_and_free(array_merged)

  def format (self) :
    return "%6.2f %6.2f %6d %6d   %5.2f %6.2f  %8.1f  %6.1f  %5.3f  %5.3f  %5.3f  %5.3f" % (
      self.d_max, self.d_min,
      self.n_obs, self.n_uniq,
      self.mean_redundancy, self.completeness*100,
      self.i_mean, self.i_over_sigma_mean,
      self.r_merge, self.r_meas, self.r_pim,
      self.cc_one_half)

  def format_for_model_cc (self) :
    return "%6.2f  %6.2f  %6d  %6.2f  %6.2f  %5.3f  %5.3f   %s   %s  %s  %s"%(
      self.d_max, self.d_min, self.n_uniq,
      self.completeness*100, self.i_over_sigma_mean,
      self.cc_one_half, self.cc_star,
      format_value("%5.3f", self.cc_work), format_value("%5.3f", self.cc_free),
      format_value("%5.3f", self.r_work), format_value("%5.3f", self.r_free))

  def format_for_gui (self) :
    return [ "%.2f - %.2f" % (self.d_max, self.d_min),
             str(self.n_obs),
             str(self.n_uniq),
             "%.1f" % self.mean_redundancy,
             "%.1f %%" % (self.completeness * 100),
             "%.1f" % self.i_over_sigma_mean,
             "%.3f" % self.r_merge,
             "%.3f" % self.r_meas,
             "%.3f" % self.r_pim,
             "%.3f" % self.cc_one_half ]

  def format_for_cc_star_gui (self) :
      return [ "%.2f - %.2f" % (self.d_max, self.d_min),
             str(self.n_uniq),
             "%.1f %%" % (self.completeness * 100),
             "%.1f" % self.i_over_sigma_mean,
             "%.3f" % self.cc_one_half,
             "%.3f" % self.cc_star,
              format_value("%5.3f", self.cc_work),
              format_value("%5.3f", self.cc_free),
              format_value("%5.3f", self.r_work),
              format_value("%5.3f", self.r_free) ]

  def table_data (self) :
    table = [(1/self.d_min**2), self.n_obs, self.n_uniq, self.mean_redundancy,
            self.completeness*100, self.i_mean, self.i_over_sigma_mean,
            self.r_merge, self.r_meas, self.r_pim, self.cc_one_half]
    if (self.cc_work is not None) :
      table.extend([self.cc_star, self.cc_work, self.cc_free, self.r_work,
        self.r_free])
    return table

  def show_summary (self, out=sys.stdout) :
    print >> out, "Resolution: %.2f - %.2f" % (self.d_max, self.d_min)
    print >> out, "Observations: %d" % self.n_obs
    print >> out, "Unique reflections: %d" % self.n_uniq
    print >> out, "Redundancy: %.1f" % self.mean_redundancy
    print >> out, "Completeness: %.2f%%" % (self.completeness*100)
    print >> out, "Mean intensity: %.1f" % self.i_mean
    print >> out, "Mean I/sigma(I): %.1f" % self.i_over_sigma_mean
    # negative sigmas are rejected before merging
    if (self.n_neg_sigmas > 0) :
      print >> out, "SigI < 0 (rejected): %d observations" % self.n_neg_sigmas
    # excessively negative intensities can be rejected either before or after
    # merging, depending on convention used
    if (self.n_rejected_before_merge > 0) :
      print >> out, "I < -3*SigI (rejected): %d observations" % \
        self.n_rejected_before_merge
    if (self.n_rejected_after_merge > 0) :
      print >> out, "I < -3*SigI (rejected): %d reflections" % \
        self.n_rejected_after_merge
    print >> out, "R-merge: %5.3f" % self.r_merge
    print >> out, "R-meas:  %5.3f" % self.r_meas
    print >> out, "R-pim:   %5.3f" % self.r_pim

class dataset_statistics (object) :
  """
  Container for overall and by-shell merging statistics, plus a table_data
  object suitable for displaying graphs (or outputting loggraph format).
  """
  def __init__ (self,
      i_obs,
      crystal_symmetry=None,
      d_min=None,
      d_max=None,
      anomalous=False,
      n_bins=10,
      debug=False,
      file_name=None,
      model_arrays=None,
      sigma_filtering=Auto,
      d_min_tolerance=1.e-6,
      estimate_cutoffs=False, # XXX slow
      log=None) :
    self.file_name = file_name
    if (log is None) : log = null_out()
    from iotbx import data_plots
    assert (i_obs.sigmas() is not None)
    info = i_obs.info()
    sigma_filtering = get_filtering_convention(i_obs, sigma_filtering)
    if (crystal_symmetry is None) :
      assert (i_obs.space_group() is not None)
      crystal_symmetry = i_obs
    i_obs = i_obs.customized_copy(
      crystal_symmetry=crystal_symmetry).set_info(info)
    if (i_obs.is_unique_set_under_symmetry()) :
      raise Sorry(("The data in %s are already merged.  Only unmerged (but "+
        "scaled) data may be used in this program.")%
        i_obs.info().label_string())
    d_min_cutoff = d_min
    d_max_cutoff = d_max
    if (d_min is not None) :
      d_min_cutoff *= (1-d_min_tolerance)
    if (d_max is not None) :
      d_max_cutoff *= 1+d_min_tolerance
    i_obs = i_obs.resolution_filter(
      d_min=d_min_cutoff,
      d_max=d_max_cutoff).set_info(info)
    i_obs.show_summary(f=log)
    self.anom_extra = ""
    if (not anomalous) :
      i_obs = i_obs.customized_copy(anomalous_flag=False).set_info(info)
      self.anom_extra = " (non-anomalous)"
    i_obs.setup_binner(n_bins=n_bins)
    merge = i_obs.merge_equivalents(use_internal_variance=False)
    self.overall = merging_stats(i_obs,
      model_arrays=model_arrays,
      anomalous=anomalous,
      debug=debug,
      sigma_filtering=sigma_filtering)
    self.bins = []
    title = "Intensity merging statistics"
    column_labels = ["1/d**2","N(obs)","N(unique)","Redundancy","Completeness",
        "Mean(I)", "Mean(I/sigma)", "R-merge", "R-meas", "R-pim", "CC1/2"]
    graph_names = ["Reflection counts", "Redundancy", "Completeness",
        "Mean(I)", "Mean(I/sigma)", "R-factors", "CC1/2"]
    graph_columns = [[0,1,2],[0,3],[0,4],[0,5],[0,6],[0,7,8,9],[0,10]]
    #--- CC* mode
    if (model_arrays is not None) :
      title = "Model quality and intensity merging statistics"
      column_labels.extend(["CC*", "CC(work)", "CC(free)", "R-work", "R-free"])
      graph_names.extend(["CC*", "Model R-factors"])
      graph_columns.extend([[0,11,12,13],[0,14,15]])
    #---
    self.table = data_plots.table_data(
      title=title,
      column_labels=column_labels,
      graph_names=graph_names,
      graph_columns=graph_columns,
      x_is_inverse_d_min=True,
      force_exact_x_labels=True)
    last_bin = None
    for bin in i_obs.binner().range_used() :
      sele_unmerged = i_obs.binner().selection(bin)
      bin_stats = merging_stats(i_obs.select(sele_unmerged),
        model_arrays=model_arrays,
        anomalous=anomalous,
        debug=debug,
        sigma_filtering=sigma_filtering)
      self.bins.append(bin_stats)
      self.table.add_row(bin_stats.table_data())
    self.cutoffs = None
    if (estimate_cutoffs) :
      self.cutoffs = estimate_crude_resolution_cutoffs(i_obs=i_obs)

  def get_estimated_cutoffs (self) :
    return getattr(self, "cutoffs", None)

  def show_estimated_cutoffs (self, out=sys.stdout, prefix="") :
    cutoffs = self.get_estimated_cutoffs()
    if (cutoffs is None) : return
    print >> out, ""
    print >> out, ""
    cutoffs.show(out=out, prefix=prefix)
    print >> out, ""
    print >> out, "NOTE: we do not endorse the use of any of these cutoffs!"
    print >> out, "It is preferrable that the resolution limit be determined"
    print >> out, "by the information content of the data, which can only"
    print >> out, "be judged by refinement."

  def show_loggraph (self, out=None) :
    if (out is None) : out = sys.stdout
    print >> out, ""
    print >> out, self.table.format_loggraph()
    print >> out, ""

  def show (self, out=None, header=True) :
    if (out is None) : out = sys.stdout
    if (header) :
      make_sub_header("Merging statistics", out=out)
    self.overall.show_summary(out)
    print >> out, ""
    print >> out, "Redundancies%s:" % self.anom_extra
    n_obs = sorted(self.overall.redundancies.keys())
    for x in n_obs :
      print >> out, "  %d : %d" % (x, self.overall.redundancies[x])
    print >> out, ""
    print >> out, """\
  Statistics by resolution bin:
 d_max  d_min   #obs  #uniq   mult.  %comp       <I>  <I/sI>  r_mrg r_meas  r_pim  cc1/2"""
    for bin_stats in self.bins :
      print >> out, bin_stats.format()
    print >> out, self.overall.format()

  def show_cc_star (self, out=None) :
    make_sub_header("CC* and related statistics", out=out)
    print >> out, """\
 d_max   d_min  n_uniq  compl. <I/sI>  cc_1/2    cc* cc_work cc_free r_work r_free"""
    for k, bin in enumerate(self.bins) :
      print >> out, bin.format_for_model_cc()
    print >> out, self.overall.format_for_model_cc()

  def extract_outer_shell_stats (self) :
    """
    For compatibility with iotbx.logfiles (which should probably now be
    deprecated) and phenix.table_one
    """
    shell = self.bins[-1]
    return group_args(
      d_max_min=(shell.d_max, shell.d_min),
      n_refl=shell.n_uniq,
      n_refl_all=shell.n_obs,
      completeness=shell.completeness,
      multiplicity=shell.mean_redundancy, # XXX bad
      r_sym=shell.r_merge,
      r_meas=shell.r_meas,
      cc_one_half=shell.cc_one_half,
      cc_star=shell.cc_star,
      i_over_sigma=shell.i_over_sigma_mean)

  def as_cif_block(self, cif_block=None):
    import iotbx.cif.model
    if cif_block is None:
      cif_block = iotbx.cif.model.block()

    observed_criterion_sigma_I = self.overall.observed_criterion_sigma_I
    if observed_criterion_sigma_I is None:
      observed_criterion_sigma_I = "?"

    cif_block["_reflns.d_resolution_low"] = self.overall.d_max
    cif_block["_reflns.d_resolution_high"] = self.overall.d_min
    cif_block["_reflns.percent_possible_obs"] = self.overall.completeness * 100
    cif_block["_reflns.pdbx_number_measured_all"] = self.overall.n_obs
    cif_block["_reflns.number_obs"] = self.overall.n_uniq
    cif_block["_reflns.pdbx_redundancy"] = self.overall.mean_redundancy
    cif_block["_reflns.phenix_mean_I"] = self.overall.i_mean
    cif_block["_reflns.pdbx_netI_over_sigmaI"] = self.overall.i_over_sigma_mean
    cif_block["_reflns.pdbx_Rmerge_I_obs"] = self.overall.r_merge
    cif_block["_reflns.pdbx_Rrim_I_obs"] = self.overall.r_meas
    cif_block["_reflns.pdbx_Rpim_I_obs"] = self.overall.r_pim
    cif_block["_reflns.phenix_cc_star"] = self.overall.cc_star
    cif_block["_reflns.phenix_cc_1/2"] = self.overall.cc_one_half
    cif_block["_reflns.observed_criterion_sigma_I"] = observed_criterion_sigma_I
    cif_block["_reflns.observed_criterion_sigma_F"] = "?"

    reflns_shell_loop = iotbx.cif.model.loop(header=(
      "_reflns_shell.d_res_high",
      "_reflns_shell.d_res_low",
      "_reflns_shell.number_measured_obs",
      "_reflns_shell.number_unique_obs",
      "_reflns_shell.pdbx_redundancy",
      "_reflns_shell.percent_possible_obs",
      "_reflns_shell.phenix_mean_I",
      "_reflns_shell.pdbx_netI_over_sigmaI_obs",
      "_reflns_shell.meanI_over_sigI_obs",
      "_reflns_shell.Rmerge_I_obs",
      "_reflns_shell.pdbx_Rrim_I_obs",
      "_reflns_shell.pdbx_Rpim_I_obs",
      "_reflns_shell.phenix_cc_star",
      "_reflns_shell.phenix_cc_1/2",
    ))
    for bin_stats in self.bins:
      reflns_shell_loop.add_row((
        bin_stats.d_min,
        bin_stats.d_max,
        bin_stats.n_obs,
        bin_stats.n_uniq,
        bin_stats.mean_redundancy,
        bin_stats.completeness*100,
        bin_stats.i_mean,
        bin_stats.i_over_sigma_mean,
        bin_stats.i_mean_over_sigi_mean,
        bin_stats.r_merge,
        bin_stats.r_meas,
        bin_stats.r_pim,
        bin_stats.cc_star,
        bin_stats.cc_one_half))
    cif_block.add_loop(reflns_shell_loop)
    return cif_block

  def as_remark_200 (self, wavelength=None) :
    from libtbx.test_utils import approx_equal
    synchrotron = wl = "NULL"
    if (wavelength is not None) :
      out = cStringIO.StringIO()
      # XXX somewhat risky...
      if (not approx_equal(wavelength, 1.5418, eps=0.01, out=out) and
          not approx_equal(wavelength, 0.7107, eps=0.01, out=out)) :
        synchrotron = "Y"
      else :
        synchrotron = "N"
      wl = "%.4f" % wavelength
    lines = []
    lines.append("")
    lines.append("EXPERIMENTAL DETAILS")
    lines.append(" EXPERIMENT TYPE                : X-RAY DIFFRACTION")
    lines.append(" DATE OF DATA COLLECTION        : NULL")
    lines.append(" TEMPERATURE           (KELVIN) : NULL")
    lines.append(" PH                             : NULL")
    lines.append(" NUMBER OF CRYSTALS USED        : NULL")
    lines.append("")
    lines.append(" SYNCHROTRON              (Y/N) : NULL")
    lines.append(" RADIATION SOURCE               : NULL")
    lines.append(" BEAMLINE                       : NULL")
    lines.append(" X-RAY GENERATOR MODEL          : NULL")
    lines.append(" MONOCHROMATIC OR LAUE    (M/L) : M")
    lines.append(" WAVELENGTH OR RANGE        (A) : %s" % wl)
    lines.append(" MONOCHROMATOR                  : NULL")
    lines.append(" OPTICS                         : NULL")
    lines.append("")
    lines.append(" DETECTOR TYPE                  : NULL")
    lines.append(" DETECTOR MANUFACTURER          : NULL")
    lines.append(" INTENSITY-INTEGRATION SOFTWARE : NULL")
    lines.append(" DATA SCALING SOFTWARE          : NULL")
    lines.append("")
    lines.append("OVERALL.")
    comp_overall = format_value("%.1f", self.overall.completeness * 100)
    mult_overall = format_value("%.1f", self.overall.mean_redundancy)
    rmerg_overall = format_value("%.5f", self.overall.r_merge)
    s2n_overall = format_value("%.4f", self.overall.i_over_sigma_mean)
    lines.append(" COMPLETENESS FOR RANGE     (%%) : %s" % comp_overall)
    lines.append(" DATA REDUNDANCY                : %s" % mult_overall)
    lines.append(" R MERGE                    (I) : %s" % rmerg_overall)
    lines.append(" R SYM                      (I) : NULL")
    lines.append(" <I/SIGMA(I)> FOR THE DATA SET  : %s" % s2n_overall)
    lines.append("")
    lines.append("IN THE HIGHEST RESOLUTION SHELL.")
    bin_stats = self.bins[-1]
    d_max = format_value("%.2f", bin_stats.d_max)
    d_min = format_value("%.2f", bin_stats.d_min)
    comp_lastbin = format_value("%.1f", bin_stats.completeness * 100)
    mult_lastbin = format_value("%.1f", bin_stats.mean_redundancy)
    rmerg_lastbin = format_value("%.5f", bin_stats.r_merge)
    s2n_lastbin = format_value("%.4f", bin_stats.i_over_sigma_mean)
    lines.append(" HIGHEST RESOLUTION SHELL, RANGE HIGH (A) : %s" % d_min)
    lines.append(" HIGHEST RESOLUTION SHELL, RANGE LOW  (A) : %s" % d_max)
    lines.append(" COMPLETENESS FOR SHELL     (%%) : %s" % comp_lastbin)
    lines.append(" DATA REDUNDANCY IN SHELL       : %s" % mult_lastbin)
    lines.append(" R MERGE FOR SHELL          (I) : %s" % rmerg_lastbin)
    lines.append(" R SYM FOR SHELL            (I) : NULL")
    lines.append(" <I/SIGMA(I)> FOR SHELL         : %s" % s2n_lastbin)
    lines.append("")
    remark_lines = [ "REMARK 200 %s" % line for line in lines ]
    return "\n".join(remark_lines)

  def show_model_vs_data (self, out=None, prefix="") :
    assert (self.overall.cc_work is not None)
    if (out is None) : out = sys.stdout
    outer_shell = self.bins[-1]
    print >> out, prefix + "Merging statistics and CC*:"
    print >> out, prefix + "  Resolution      : %.3f - %.3f (%.3f - %.3f)" % (
      self.overall.d_max, self.overall.d_min, outer_shell.d_max,
      outer_shell.d_min)
    print >> out, prefix + "  Mean(I/sigmaI)  : %6.3f (%.3f)" % (
      self.overall.i_over_sigma_mean, outer_shell.i_over_sigma_mean)
    print >> out, prefix + "  Redundancy      :  %4.2f  (%.2f)" % (
      self.overall.mean_redundancy, outer_shell.mean_redundancy)
    print >> out, prefix + "  R-merge         :  %5.3f (%.3f)" % (
      self.overall.r_merge, outer_shell.r_merge)
    print >> out, prefix + "  R-meas          :  %5.3f (%.3f)" % (
      self.overall.r_meas, outer_shell.r_meas)
    print >> out, prefix + "  R-pim           :  %5.3f (%.3f)" % (
      self.overall.r_pim, outer_shell.r_pim)
    print >> out, prefix + "  CC1/2           :  %5.3f (%.3f)" % (
      self.overall.cc_one_half, outer_shell.cc_one_half)
    print >> out, prefix + "  CC*             :  %5.3f (%.3f)" % (
      self.overall.cc_star, outer_shell.cc_star)
    print >> out, prefix + "  CC(work)        :  %6.4f (%.4f)" % (
      self.overall.cc_work, outer_shell.cc_work)
    if (self.overall.cc_free is not None) :
      print >> out, prefix + "  CC(free)        :  %6.4f (%.4f)" % (
        self.overall.cc_free, outer_shell.cc_free)
    else :
      print >> out, prefix + "  CC(free)        :  not available"

def select_data (file_name, data_labels, log=None,
    assume_shelx_observation_type_is=None) :
  if (log is None) : log = null_out()
  from iotbx import reflection_file_reader
  hkl_in = reflection_file_reader.any_reflection_file(file_name)
  print >> log, "Format:", hkl_in.file_type()
  miller_arrays = hkl_in.as_miller_arrays(merge_equivalents=False,
    assume_shelx_observation_type_is=assume_shelx_observation_type_is)
  if ((hkl_in.file_type() == "shelx_hklf") and (not "=" in file_name)) :
    print >> log, "WARNING: SHELX file is assumed to contain intensities"
  i_obs = None
  all_i_obs = []
  for array in miller_arrays :
    labels = array.info().label_string()
    if (labels == data_labels) :
      i_obs = array
      break
    elif (array.is_xray_intensity_array()) :
      all_i_obs.append(array)
  if (i_obs is None) :
    if (len(all_i_obs) == 0) :
      raise Sorry("No intensities found in %s." % file_name)
    elif (len(all_i_obs) > 1) :
      raise Sorry("Multiple intensity arrays - please specify one:\n%s" %
        "\n".join(["  labels=%s"%a.info().label_string() for a in all_i_obs]))
    else :
      i_obs = all_i_obs[0]
  if (not i_obs.is_xray_intensity_array()) :
    raise Sorry("%s is not an intensity array." % i_obs.info().label_string())
  return i_obs

class estimate_crude_resolution_cutoffs (slots_getstate_setstate) :
  """
  Uses incredibly simplistic criteria to determine the approximate
  resolution limit of the data based on the merging statistics (using much
  smaller bins than normal).  Not really appropriate for routine use, but
  useful for the pedantic and just-curious.
  """
  cutoffs_attr = ["i_over_sigma_cut", "r_merge_cut", "r_meas_cut",
    "completeness_cut_conservative", "completeness_cut_permissive",
    "cc_one_half_cut",]
  __slots__ = ["n_bins", "i_over_sigma_min", "r_merge_max", "r_meas_max",
    "completeness_min_conservative", "completeness_min_permissive",
    "cc_one_half_min", "d_min_overall",] + cutoffs_attr
  def __init__ (self,
      i_obs,
      crystal_symmetry=None,
      n_bins=100,
      i_over_sigma_min=2.0,
      r_merge_max=0.5,
      r_meas_max=0.5,
      completeness_min_conservative=0.9,
      completeness_min_permissive=0.5,
      cc_one_half_min=0.5) :
    self.n_bins = n_bins
    self.i_over_sigma_min = i_over_sigma_min
    self.r_merge_max = r_merge_max
    self.r_meas_max = r_meas_max
    self.completeness_min_conservative = completeness_min_conservative
    self.completeness_min_permissive = completeness_min_permissive
    self.cc_one_half_min = cc_one_half_min
    for attr in self.cutoffs_attr :
      setattr(self, attr, None)
    stats = dataset_statistics(
      i_obs=i_obs,
      crystal_symmetry=crystal_symmetry)
    self.d_min_overall = stats.overall.d_min
    d_min_last = float("inf")
    for bin in stats.bins :
      if (self.i_over_sigma_cut is None) :
        if (bin.i_over_sigma_mean < self.i_over_sigma_min) :
          self.i_over_sigma_cut = d_min_last
      if (self.cc_one_half_cut is None) :
        if (bin.cc_one_half < self.cc_one_half_min) :
          self.cc_one_half_cut = d_min_last
      if (self.r_merge_cut is None) :
        if (bin.r_merge > self.r_merge_max) :
          self.r_merge_cut = d_min_last
      if (self.r_meas_cut is None) :
        if (bin.r_meas > self.r_meas_max) :
          self.r_meas_cut = d_min_last
      if (self.completeness_cut_conservative is None) :
        if (bin.completeness < completeness_min_conservative) :
          self.completeness_cut_conservative = d_min_last
      if (self.completeness_cut_permissive is None) :
        if (bin.completeness < completeness_min_permissive) :
          self.completeness_cut_permissive = d_min_last
      d_min_last = bin.d_min

  def show (self, out=sys.stdout, prefix="") :
    def format_d_min (value) :
      if (value is None) :
        return "(use all data)" #% self.d_min_overall
      return "%7.3f" % value
    print >> out, prefix + "Crude resolution cutoff estimates:"
    print >> out, prefix + "  resolution of all data          : %7.3f" % \
      self.d_min_overall
    if (self.cc_one_half_min is not None) :
      print >> out, prefix + "  based on CC(1/2) > %-6g       : %s" % \
        (self.cc_one_half_min, format_d_min(self.cc_one_half_cut))
    if (self.i_over_sigma_min is not None) :
      print >> out, prefix + "  based on mean(I/sigma) > %-6g : %s" % \
        (self.i_over_sigma_min, format_d_min(self.i_over_sigma_cut))
    if (self.r_merge_max is not None) :
      print >> out, prefix + "  based on R-merge < %-6g       : %s" % \
        (self.r_merge_max, format_d_min(self.r_merge_cut))
    if (self.r_meas_max is not None) :
      print >> out, prefix + "  based on R-meas < %-6g        : %s" % \
        (self.r_meas_max, format_d_min(self.r_meas_cut))
    if (self.completeness_min_conservative is not None) :
      print >> out, prefix + "  based on completeness > %-6g  : %s" % \
        (self.completeness_min_conservative,
         format_d_min(self.completeness_cut_conservative))
    if (self.completeness_min_permissive is not None) :
      print >> out, prefix + "  based on completeness > %-6g  : %s" % \
        (self.completeness_min_permissive,
         format_d_min(self.completeness_cut_permissive))
