# LIBTBX_SET_DISPATCHER_NAME phenix.merging_statistics

from __future__ import division
from libtbx.str_utils import make_sub_header, format_value
from libtbx.utils import Sorry, Usage, null_out
from libtbx import runtime_utils
from libtbx import group_args
from math import sqrt
import sys

citations_str = """\
  Diederichs K & Karplus PA (1997) Nature Structural Biology 4:269-275
    (with erratum in: Nat Struct Biol 1997 Jul;4(7):592)
  Weiss MS (2001) J Appl Cryst 34:130-135.
  Karplus PA & Diederichs K (2012) Science 336:1030-3."""

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
"""

master_phil = """
file_name = None
  .type = path
  .short_caption = Unmerged data
  .style = file_type:hkl OnChange:extract_unmerged_intensities bold
labels = None
  .type = str
  .input_size = 200
  .style = renderer:draw_unmerged_intensities_widget
space_group = None
  .type = space_group
  .input_size = 120
unit_cell = None
  .type = unit_cell
symmetry_file = None
  .type = path
  .style = file_type:pdb,hkl OnChange:extract_symmetry
%s
debug = False
  .type = bool
loggraph = False
  .type = bool
include scope libtbx.phil.interface.tracking_params
""" % merging_params_str

class merging_stats (object) :
  """
  Calculate standard merging statistics for (scaled) unmerged data.  Usually
  these statistics will consider I(+) and I(-) as observations of the same
  reflection, but these can be kept separate instead if desired.

  Reflections with negative sigmas will be discarded, and also, per the
  recommendation of Kay Diederichs, reflections where I < -3 * sigmaI.
  """
  def __init__ (self, array, anomalous=False, debug=None) :
    import cctbx.miller
    from scitbx.array_family import flex
    assert (array.sigmas() is not None)
    array = array.customized_copy(anomalous_flag=anomalous).map_to_asu()
    non_negative_sel = array.sigmas() >= 0
    self.n_neg_sigmas = non_negative_sel.count(False)
    array = array.select(non_negative_sel)
    merge = array.merge_equivalents()
    array_merged = merge.array()
    reject_sel = (array_merged.data() < -3*array_merged.sigmas())
    self.n_rejected = reject_sel.count(True)
    array_merged = array_merged.select(~reject_sel)
    self.d_max, self.d_min = array.d_max_min()
    self.n_obs = array.indices().size()
    self.n_uniq = array_merged.indices().size()
    complete_set = array_merged.complete_set().resolution_filter(
      d_min=self.d_min, d_max=self.d_max)
    n_expected = len(complete_set.indices())
    if (n_expected == 0) :
      raise RuntimeError(("No reflections within specified resolution range "+
        "(%g - %g)") % (self.d_max, self.d_min))
    self.completeness = min(self.n_uniq / n_expected, 1.)
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

  def format (self) :
    return "%6.2f  %6.2f %6d %6d   %5.2f %6.2f  %8.1f  %6.1f  %5.3f  %5.3f  %5.3f  %5.3f" % (
      self.d_max, self.d_min,
      self.n_obs, self.n_uniq,
      self.mean_redundancy, self.completeness*100,
      self.i_mean, self.i_over_sigma_mean,
      self.r_merge, self.r_meas, self.r_pim,
      self.cc_one_half)

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

  def table_data (self) :
    return [(1/self.d_min**2), self.n_obs, self.n_uniq, self.mean_redundancy,
            self.completeness*100, self.i_mean, self.i_over_sigma_mean,
            self.r_merge, self.r_meas, self.r_pim, self.cc_one_half]

  def show_summary (self, out=sys.stdout) :
    print >> out, "Resolution: %.2f - %.2f" % (self.d_max, self.d_min)
    print >> out, "Observations: %d" % self.n_obs
    print >> out, "Unique reflections: %d" % self.n_uniq
    print >> out, "Redundancy: %.1f" % self.mean_redundancy
    print >> out, "Completeness: %.2f%%" % (self.completeness*100)
    print >> out, "Mean intensity: %.1f" % self.i_mean
    print >> out, "Mean I/sigma(I): %.1f" % self.i_over_sigma_mean
    if (self.n_neg_sigmas > 0) :
      print >> out, "SigI < 0 (rejected): %d reflections" % self.n_neg_sigmas
    if (self.n_rejected > 0) :
      print >> out, "I < -3*SigI (rejected): %d reflections" % self.n_rejected
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
      log=None) :
    self.file_name = file_name
    if (log is None) : log = null_out()
    from iotbx import data_plots
    assert (i_obs.sigmas() is not None)
    info = i_obs.info()
    if (crystal_symmetry is None) :
      assert (i_obs.space_group() is not None)
      crystal_symmetry = i_obs
    i_obs = i_obs.customized_copy(
      crystal_symmetry=crystal_symmetry).set_info(info)
    if (i_obs.is_unique_set_under_symmetry()) :
      raise Sorry(("The data in %s are already merged.  Only unmerged (but "+
        "scaled) data may be used in this program.")%
        i_obs.info().label_string())
    i_obs = i_obs.resolution_filter(
      d_min=d_min,
      d_max=d_max).set_info(info)
    i_obs.show_summary(f=log)
    self.anom_extra = ""
    if (not anomalous) :
      i_obs = i_obs.customized_copy(anomalous_flag=False)
      self.anom_extra = " (non-anomalous)"
    i_obs.setup_binner(n_bins=n_bins)
    merge = i_obs.merge_equivalents()
    self.overall = merging_stats(i_obs, anomalous=anomalous, debug=debug)
    self.bins = []
    self.table = data_plots.table_data(
      title="Intensity merging statistics",
      column_labels=["1/d**2","N(obs)","N(unique)","Redundancy","Completeness",
        "Mean(I)", "Mean(I/sigma)", "R-merge", "R-meas", "R-pim", "CC1/2"],
      graph_names=["Reflection counts", "Redundancy", "Completeness",
        "Mean(I)", "Mean(I/sigma)", "R-factors", "CC1/2"],
      graph_columns=[[0,1,2],[0,3],[0,4],[0,5],[0,6],[0,7,8,9],[0,10]],
      x_is_inverse_d_min=True,
      force_exact_x_labels=True)
    last_bin = None
    for bin in i_obs.binner().range_used() :
      sele_unmerged = i_obs.binner().selection(bin)
      bin_stats = merging_stats(i_obs.select(sele_unmerged),
        anomalous=anomalous,
        debug=debug)
      self.bins.append(bin_stats)
      self.table.add_row(bin_stats.table_data())

  def show_loggraph (self, out=None) :
    if (out is None) : out = sys.stdout
    print >> out, ""
    print >> out, self.table.format_loggraph()
    print >> out, ""

  def show (self, out=None) :
    if (out is None) : out = sys.stdout
    make_sub_header("Merging statistics", out=out)
    self.overall.show_summary(out)
    print >> out, ""
    print >> out, "Redundancies%s:" % self.anom_extra
    n_obs = sorted(self.overall.redundancies.keys())
    for x in n_obs :
      print "  %d : %d" % (x, self.overall.redundancies[x])
    print >> out, ""
    print >> out, """\
  Statistics by resolution bin:
 d_min   d_max   #obs  #uniq   mult.  %comp       <I>  <I/sI>  r_mrg r_meas  r_pim  cc1/2"""
    for bin_stats in self.bins :
      print >> out, bin_stats.format()
    print >> out, self.overall.format()

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
      i_over_sigma=shell.i_over_sigma_mean)

def select_data (file_name, data_labels, out) :
  from iotbx import reflection_file_reader
  hkl_in = reflection_file_reader.any_reflection_file(file_name)
  print >> out, "Format:", hkl_in.file_type()
  miller_arrays = hkl_in.as_miller_arrays(merge_equivalents=False)
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

def run (args, out=None) :
  if (out is None) : out = sys.stdout
  import iotbx.phil
  master_params = iotbx.phil.parse(master_phil, process_includes=True)
  if (len(args) == 0) :
    raise Usage("""\
phenix.merging_statistics [data_file] [options...]

Calculate merging statistics for non-unique data, including R-merge, R-meas,
R-pim, and redundancy.  Any format supported by Phenix is allowed, including
MTZ, unmerged Scalepack, or XDS/XSCALE (and possibly others).  Data should
already be on a common scale, but with individual observations unmerged.
%s

Full parameters:
%s
""" % (citations_str, master_params.as_str(prefix="  ")))
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_params,
    reflection_file_def="file_name",
    pdb_file_def="symmetry_file")
  params = cmdline.work.extract()
  i_obs = select_data(
    file_name=params.file_name,
    data_labels=params.labels,
    out=out)
  symm = None
  if (params.symmetry_file is not None) :
    from iotbx import crystal_symmetry_from_any
    symm = crystal_symmetry_from_any.extract_from(
      file_name=params.symmetry_file)
    if (symm is None) :
      raise Sorry("No symmetry records found in %s." % params.symmetry_file)
  if (symm is None) :
    sg = i_obs.space_group()
    if (sg is None) :
      if (params.space_group is not None) :
        sg = params.space_group.group()
      else :
        raise Sorry("Missing space group information.")
    uc = i_obs.unit_cell()
    if (uc is None) :
      if (params.unit_cell is not None) :
        uc = params.unit_cell
      else :
        raise Sorry("Missing unit cell information.")
    from cctbx import crystal
    symm = crystal.symmetry(
      space_group=sg,
      unit_cell=uc)
  if (i_obs.sigmas() is None) :
    raise Sorry("Sigma(I) values required for this application.")
  result = dataset_statistics(
    i_obs=i_obs,
    crystal_symmetry=symm,
    d_min=params.high_resolution,
    d_max=params.low_resolution,
    n_bins=params.n_bins,
    anomalous=params.anomalous,
    debug=params.debug,
    file_name=params.file_name,
    log=out)
  result.show(out=out)
  if (params.loggraph) :
    result.show_loggraph()
  print >> out, "References:"
  print >> out, citations_str
  return result

#-----------------------------------------------------------------------
# Phenix GUI stuff
def validate_params (params) :
  if (params.file_name is None) :
    raise Sorry("No data file specified!")
  elif (params.labels is None) :
    raise Sorry("No data labels selected!")
  return True

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    return run(args=list(self.args), out=sys.stdout)

def finish_job (result) :
  stats = []
  if (result is not None) :
    stats = [
      ("High resolution", format_value("%.3g", result.overall.d_min)),
      ("Redundancy", format_value("%.1f", result.overall.mean_redundancy)),
      ("R-meas", format_value("%.3g", result.overall.r_meas)),
      ("R-meas (high-res)", format_value("%.3g", result.bins[-1].r_meas)),
      ("<I/sigma>", format_value("%.2g", result.overall.i_over_sigma_mean)),
      ("<I/sigma> (high-res)", format_value("%.2g",
        result.bins[-1].i_over_sigma_mean)),
      ("Completeness", format_value("%.1f%%", result.overall.completeness*100)),
      ("Completeness (high-res)", format_value("%.1f%%",
        result.bins[-1].completeness*100)),
    ]
  return ([], stats)

if (__name__ == "__main__") :
  run(sys.argv[1:])
