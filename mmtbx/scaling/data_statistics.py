
"""
Collect multiple analyses of experimental data quality, including
signal-to-noise ratio, completeness, ice rings and other suspicious
outliers, anomalous measurability, and Wilson plot.
"""

from __future__ import absolute_import, division, print_function
from mmtbx.scaling import absolute_scaling
from mmtbx import scaling
from iotbx import data_plots
from cctbx import miller
from cctbx import adptbx
from cctbx.array_family import flex
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
from scitbx.math import erf
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
from libtbx import table_utils
import math
from six.moves import zip
from six.moves import range


class i_sigi_completeness_stats(scaling.xtriage_analysis):
  """
  Collects resolution-dependent statistics on I/sigma expressed as percentage
  of reflections above specified cutoffs.
  """
  def __init__(self,
               miller_array,
               n_bins=15,
               isigi_cut=3.0,
               completeness_cut=0.85,
               resolution_at_least=3.5,
               completeness_as_non_anomalous=None):
    miller_array = miller_array.deep_copy().set_observation_type(miller_array)
    self.amplitudes_flag = False
    #check to see if we have intensities
    if miller_array.is_real_array():
      if miller_array.is_xray_amplitude_array():
        self.amplitudes_flag = True
        miller_array = miller_array.f_as_f_sq()
    assert miller_array.is_xray_intensity_array()
    #make sure we have sigmas
    assert miller_array.sigmas() is not None
    # select things with sigma larger then zero
    #miller_array = miller_array.select( miller_array.sigmas() > 0 )
    miller_array.setup_binner(n_bins=n_bins)
    self.resolution_bins = list(miller_array.binner().limits())
    self.overall = miller_array.completeness(use_binning=False,
       as_non_anomalous_array = completeness_as_non_anomalous)
    self.overall_binned = miller_array.completeness(use_binning=True,
       as_non_anomalous_array = completeness_as_non_anomalous).data[1:-1]
    self.completeness_bins = []
    def cut_completeness(cut_value):
      tmp_miller = miller_array.select(
        miller_array.data() > cut_value*miller_array.sigmas() )
      tmp_miller.use_binning_of(miller_array)
      completeness = tmp_miller.completeness(use_binning=True,
               return_fail=0.0,
               as_non_anomalous_array = completeness_as_non_anomalous)
      return completeness.data
    # create table
    table_data = []
    cut_list=[1,2,3,5,10,15]
    for cut_value in cut_list:
      self.completeness_bins.append( cut_completeness(cut_value) )
    legend = ["Res. Range", "I/sigI>1 ", "I/sigI>2 ", "I/sigI>3 ",
              "I/sigI>5 ", "I/sigI>10", "I/sigI>15" ]
    self.table = data_plots.table_data(
      title="Completeness and data strength",
      column_labels=["Low res.", "High res."] + legend[1:],
      column_formats=["%4.2f", "%4.2f" ] + [ "%5.1f" ] * 6,
      graph_names=["I/sigI by shell"],
      graph_labels=[("High resolution of shell", "% of total")],
      graph_columns=[list(range(1,8))],
      x_is_inverse_d_min=True,
      first_two_columns_are_resolution=True)
    for ii in range(1,len(self.resolution_bins)-1):
      a = self.resolution_bins[ii-1]
      b = self.resolution_bins[ii]
      row = [ a, b ]
      for jj in  self.completeness_bins:
        row.append(100.0*jj[ii])
      self.table.add_row(row)
    self.isigi_cut=isigi_cut
    self.completeness_cut=completeness_cut
    self.resolution_at_least=resolution_at_least
    self.resolution_cut = 4.0
    comp_data = cut_completeness(isigi_cut)
    reso=4.0
    for ii in range(1,len(self.resolution_bins)-1):
      a = self.resolution_bins[ii-1]
      b = self.resolution_bins[ii]
      if b < resolution_at_least:
        if comp_data[ii] > completeness_cut:
          self.resolution_cut = b**-0.5

  def _show_impl(self, out):
    out.show_sub_header("Completeness at I/sigma cutoffs")
    out.show_text("""
 The following table lists the completeness in various resolution ranges,
 after applying a I/sigI cut. Miller indices for which individual I/sigI
 values are larger than the value specified in the top row of the table, are
 retained, while other intensities are discarded. The resulting completeness
 profiles are an indication of the strength of the data.
""")
    if (self.amplitudes_flag):
      out.warn("""\
Please be aware that the input data were given as amplitudes and squared for
the purposes of this analysis, therefore the numbers displayed here are less
reliable than the values calculated from the raw intensities.
""")
    out.show_table(self.table, indent=2)
    out.show_plot(self.table)
    message = """\
  The completeness of data for which I/sig(I)>%.2f, exceeds %.0f %%
  for resolution ranges lower than %.2fA.
""" % (self.isigi_cut, self.completeness_cut*100, self.resolution_cut)

    if self.resolution_cut < self.resolution_at_least:
      message += """\
  The data are cut at this resolution for the potential twin tests and
  intensity statistics.
"""
    else:
      message += """\
  As we do not want to throw away too much data, the resolution for
  analyzing the intensity statistics will be limited to %3.2fA.
""" % self.resolution_at_least
    out.show_text(message)


# XXX where is this used?
class completeness_enforcement(object):
  def __init__(self,
               miller_array,
               minimum_completeness=0.75,
               completeness_as_non_anomalous=None):

    self.miller_array = miller_array.deep_copy()
    self.miller_array.setup_binner_d_star_sq_step(auto_binning=True)
    completeness = self.miller_array.completeness(use_binning=True,
      return_fail=1.0,
      as_non_anomalous_array = completeness_as_non_anomalous)
    selection_array = flex.bool( self.miller_array.indices().size(), False )
    #In this pass, we make sure that we have reasonable completeness
    for bin in completeness.binner.range_all():
      selection = completeness.binner.selection(bin).iselection()
      if completeness.data[bin] >= minimum_completeness:
        # the completeness is okai
        # use these indices please
        selection_array = selection_array.set_selected( selection, True )
    # now select the indices please
    self.new_miller = miller_array.select( selection_array )

class analyze_resolution_limits(scaling.xtriage_analysis):
  """
  Check for elliptical truncation, which may be applied to the data by some
  processing software (or as a post-processing step).  As a general rule this
  is not recommended since phenix.refine and related programs will handle
  anisotropy automatically, and users tend to apply it blindly (and even
  deposit the modified data).
  """
  def __init__(self, miller_set, d_min_max_delta=0.25):
    # XXX very important - for high-symmetry space groups we need to examine
    # all possible positive indices
    tmp_miller = miller_set.expand_to_p1()
    if (miller_set.unit_cell().parameters()[3:] != (90,90,90)):
    #if (miller_set.space_group().n_smx() > 2):
      all_pos_neg_indices = flex.miller_index()
      for h,k,l in tmp_miller.indices():
        if (h >= 0 and k >= 0 and l >= 0) or (h <= 0 and k <= 0 and l <= 0):
          all_pos_neg_indices.append((h,k,l))
      if isinstance(tmp_miller, miller.array):
        tmp_miller = tmp_miller.set(indices=all_pos_neg_indices)
      else :
        tmp_miller = tmp_miller.customized_copy(indices=all_pos_neg_indices)
    else :
      tmp_miller = tmp_miller.expand_to_p1()
    self.d_min_overall = tmp_miller.d_min()
    # FIXME this still fails for some space groups - need to test on entire
    # PDB
    assert approx_equal(self.d_min_overall, miller_set.d_min(), eps=0.1)
    self.d_min_a, self.d_min_b, self.d_min_c = \
      tmp_miller.d_min_along_a_b_c_star()
    self.d_min_max_delta = d_min_max_delta

  def max_d_min_delta(self):
    """
    Return the maximum difference in d_min along any two axes.
    """
    limits = [self.d_min_a, self.d_min_b, self.d_min_c]
    max_delta = 0
    for i in range(2):
      for j in range(i,3):
        max_delta = max(max_delta, abs(limits[i] - limits[j]))
    return max_delta

  def is_elliptically_truncated(self, d_min_max_delta=None):
    if (d_min_max_delta is None):
      d_min_max_delta = self.d_min_max_delta
    return self.max_d_min_delta() > d_min_max_delta

  def _show_impl(self, out):
    out.show_sub_header("Analysis of resolution limits")
    out.show_text("""\
Your data have been examined to determine the resolution limits of the data
along the reciprocal space axes (a*, b*, and c*).  These are expected to vary
slightly depending on unit cell parameters and overall resolution, but should
never be significantly different for complete data.  (This is distinct from the
amount of anisotropy present in the data, which changes the effective
resolution but does not actually exclude reflections.)""")
    max_delta = self.max_d_min_delta()
    out.show_preformatted_text("""
    overall d_min                = %.3f
    d_min along a*               = %.3f
    d_min along b*               = %.3f
    d_min along c*               = %.3f
    max. difference between axes = %.3f
""" % (self.d_min_overall, self.d_min_a,
       self.d_min_b, self.d_min_c, max_delta))
    if (max_delta > self.d_min_max_delta):
      out.show_text("""\
The resolution limit appears to be direction-dependent, which may indicate
issues with the data collection geometry, processing errors, or that elliptical
truncation has been applied.  We do not recommend using elliptically truncated
data, as anisotropy is handled automatically by Phaser, phenix.refine, and
related programs, and discarding large numbers of weak reflections may result
in increased map bias and/or artifacts.  You should always deposit the original,
uncorrected reflections in the PDB, not the truncated data.""")
    else :
      out.show_text("""Resolution limits are within expected tolerances.""")


class log_binned_completeness(scaling.xtriage_analysis):
  """
  Table of completeness using log-scale resolution binning.
  """
  def __init__(self, miller_array,
      n_reflections_in_lowest_resolution_bin=100,
      max_number_of_bins=30,
      min_reflections_in_bin=50,
      completeness_as_non_anomalous=None):
    rows = []
    bins = miller_array.log_binning(
      n_reflections_in_lowest_resolution_bin=\
        n_reflections_in_lowest_resolution_bin,
      max_number_of_bins=max_number_of_bins,
      min_reflections_in_bin=min_reflections_in_bin)
    for selection in bins :
      bin_array = miller_array.select(selection)
      d_max, d_min = bin_array.d_max_min()
      n_refl = bin_array.size()
      n_refl_expect = bin_array.complete_set(d_max=d_max).size()
      completeness = bin_array.completeness(d_max=d_max,
          as_non_anomalous_array = completeness_as_non_anomalous)
      rows.append(("%.4f - %.4f" % (d_max, d_min), "%d/%d" % (n_refl,
        n_refl_expect), "%.1f%%" % (completeness*100)))
    self.table = table_utils.simple_table(
      column_headers=["Resolution", "Reflections", "Completeness"],
      table_rows=rows)

  def _show_impl(self, out):
    out.show_sub_header("Completeness (log-binning)")
    out.show_text("""\
The table below presents an alternative overview of data completeness, using
the entire resolution range but on a logarithmic scale.  This is more sensitive
to missing low-resolution data (and is complementary to the separate table
showing low-resolution completeness only).""")
    out.show_table(self.table)

#-----------------------------------------------------------------------
# OUTLIER FILTERING
class possible_outliers(scaling.xtriage_analysis):
  """
  Flag specific reflections with suspicious intensities.  Inspired by:
  Read, Acta Cryst. (1999). D55, 1759-1764
  """
  def __init__(self,
               miller_array,
               prob_cut_ex=1.0E-1,
               prob_cut_wil=1.0E-6):
    if miller_array.observation_type() is None:
      raise Sorry("Unknown observation type")
    if not miller_array.is_xray_intensity_array():
      miller_array = miller_array.f_as_f_sq()
    else:
      miller_array = miller_array.deep_copy()
    assert miller_array.is_xray_intensity_array()
    work_array = miller_array.deep_copy()
    normalizer = absolute_scaling.kernel_normalisation(
      work_array, auto_kernel=True)
    work_array = work_array.array(
      data=normalizer.normalised_miller.data()
      /work_array.epsilons().data().as_double())
    self.n_refl = work_array.size()
    # array processing
    centric_cut = work_array.select_centric().set_observation_type(work_array)
    acentric_cut = work_array.select_acentric().set_observation_type(
      work_array)
    assert centric_cut.data().all_ge(0)
    p_acentric_single = 1.0 - (1.0 - flex.exp(-acentric_cut.data() ))
    p_centric_single = 1.0 - erf(flex.sqrt(centric_cut.data()/2.0) )
    n_centric = p_centric_single.size()
    n_acentric = p_acentric_single.size()
    extreme_acentric = 1.0 -  \
       flex.pow(1.0 - flex.exp(-acentric_cut.data() ),float(n_acentric))
    extreme_centric = 1.0 - \
       flex.pow(erf(flex.sqrt(centric_cut.data()/2.0) ),float(n_centric))

    ## combine both the wilson and extreme value cut-off values
    acentric_outlier = (extreme_acentric < prob_cut_ex) or (
     p_acentric_single < prob_cut_wil)
    centric_outlier = (extreme_centric < prob_cut_ex) or (
     p_centric_single  < prob_cut_wil)

    ## acentrics
    self.acentric_outlier_miller = acentric_cut.indices().select(
      acentric_outlier)
    self.acentric_outlier_e_vals = acentric_cut.data().select(
      acentric_outlier)
    self.acentric_outlier_e_vals = flex.sqrt(self.acentric_outlier_e_vals)
    self.acentric_outlier_d_spacings = acentric_cut.d_spacings().data()\
                                       .select(acentric_outlier)
    self.acentric_outlier_p_val = p_acentric_single.select(acentric_outlier)
    self.acentric_outlier_extreme_val = extreme_acentric.select(
      acentric_outlier)
    self.acentric_outliers_table = data_plots.table_data(
      title="Acentric reflections",
      column_labels=["d_spacing", "H K L", "|E|", "p(wilson)", "p(extreme)"],
      column_formats=["%8.3f","%5i,%5i,%5i", "%6.2f", "%9.2e", "%10.2e"],
      graph_names=["Possible acentric outliers"],
      graph_columns=[[0,1,2,3,4]])
    for d,hkl,e,p,extr in zip(self.acentric_outlier_d_spacings,
                              self.acentric_outlier_miller,
                              self.acentric_outlier_e_vals,
                              self.acentric_outlier_p_val,
                              self.acentric_outlier_extreme_val):
      self.acentric_outliers_table.add_row([d, hkl, e, p, extr])

    ## centrics
    self.centric_outlier_miller = centric_cut.indices().select(
      centric_outlier)
    self.centric_outlier_e_vals = centric_cut.data().select(
      centric_outlier)
    self.centric_outlier_e_vals = flex.sqrt(self.centric_outlier_e_vals)
    self.centric_outlier_d_spacings = centric_cut.d_spacings().data()\
                                       .select(centric_outlier)
    self.centric_outlier_p_val = p_centric_single.select(centric_outlier)
    self.centric_outlier_extreme_val = extreme_centric.select(centric_outlier)
    self.centric_outliers_table = data_plots.table_data(
      title="Centric reflections",
      column_labels=["d_spacing", "H K L", "|E|", "p(wilson)", "p(extreme)"],
      column_formats=["%8.3f","%5i,%5i,%5i", "%6.2f", "%9.2e", "%10.2e"],
      graph_names=["Possible centric outliers"],
      graph_columns=[[0,1,2,3,4]])
    for d,hkl,e,p,extr in zip(self.centric_outlier_d_spacings,
                              self.centric_outlier_miller,
                              self.centric_outlier_e_vals,
                              self.centric_outlier_p_val,
                              self.centric_outlier_extreme_val):
      self.centric_outliers_table.add_row([d, hkl, e, p, extr])

  def fraction_outliers(self):
    return self.n_outliers() / self.n_refl

  def n_outliers(self):
    n_acentric = self.acentric_outlier_miller.size()
    n_centric = self.centric_outlier_miller.size()
    return n_acentric + n_centric

  def remove_outliers(self, miller_array):
    ## remove the outliers please
    centric_matches = miller.match_indices( miller_array.indices(),
                                     self.centric_outlier_miller )
    miller_array = miller_array.select( centric_matches.single_selection(0) )
    acentric_matches = miller.match_indices( miller_array.indices(),
                                      self.acentric_outlier_miller )
    miller_array = miller_array.select( acentric_matches.single_selection(0) )
    return miller_array

  def _show_impl(self, out):
    out.show_sub_header("Possible outliers")
    out.show("  Inspired by: Read, Acta Cryst. (1999). D55, 1759-1764\n")
    out.show("Acentric reflections:")
    if self.acentric_outlier_miller.size() ==0:
      out.show("            None\n")
    else:
      out.show_table(self.acentric_outliers_table, indent=2)
      out.show_preformatted_text("""\n
 p(wilson)  : 1-(1-exp[-|E|^2])
 p(extreme) : 1-(1-exp[-|E|^2])^(n_acentrics)
""")
      out.show("""
 p(wilson) is the probability that an E-value of the specified value would be
 observed if it were selected at random the given data set. p(extreme) is the
 probability that the largest |E| value is larger or equal than the observed
 largest |E| value.

 Both measures can be used for outlier detection. p(extreme) takes into
 account the size of the dataset.
""")

    out.show("Centric reflections:")
    if self.centric_outlier_miller.size() ==0:
      out.show("            None\n")
    else :
      out.show_table(self.centric_outliers_table, indent=2)
      out.show_preformatted_text("""\n
 p(wilson)  : 1-(erf[|E|/sqrt(2)])
 p(extreme) : 1-(erf[|E|/sqrt(2)])^(n_acentrics)
""")
      out.show("""
 p(wilson) is the probability that an E-value of the specified
 value would be observed when it would selected at random from
 the given data set.
 p(extreme) is the probability that the largest |E| value is
 larger or equal than the observed largest |E| value.

 Both measures can be used for outlier detection. p(extreme)
 takes into account the size of the dataset.
""")

  def summarize_issues(self):
    if (self.fraction_outliers() > 0.001):
      return [(1, "There are a large number of outliers in the data.",
        "Possible outliers")]
    else :
      return [(0, "The fraction of outliers in the data is less than 0.1%.",
        "Possible outliers")]

#-----------------------------------------------------------------------
# ICE RINGS
class ice_ring_checker(scaling.xtriage_analysis):
  """
  Check intensity and completeness statistics in specific resolution ranges
  known to have strong diffraction when crystalline ice is present.
  """
  def __init__(self,
               bin_centers,
               completeness_data,
               z_scores_data,
               completeness_abnormality_level=4.0,
               intensity_level=0.1,
               z_score_limit=10):
    self.completeness_abnormality_level = completeness_abnormality_level
    self.intensity_level = intensity_level
    self.z_score_limit = z_score_limit
    self.ice_d_spacings=flex.double(
      [3.897,3.669,3.441,2.671,2.249,
       2.072,1.948,1.918,1.883,1.721])
    self.ice_rel_intens=flex.double(
      [1.000, 0.750, 0.530, 0.170, 0.390,
       0.300, 0.040, 0.180, 0.030, 0.020])
    self.ice_ring_bin_location=\
      [None, None, None, None, None,
       None, None, None, None, None]
    self.mean_comp=None
    self.mean_z_score=None
    tmp_low_d_star_sq=bin_centers[0]
    tmp_high_d_star_sq = bin_centers[bin_centers.size()-1]
    tmp_step = (bin_centers[1]-bin_centers[0])
    count=0
    ## array of weights masking ice ring sensitive areas
    weights = flex.double( bin_centers.size(), 1)
    for ice_ring in self.ice_d_spacings:
      tmp_ice_ring = 1.0/(ice_ring**2.0)
      tmp_ice_ring_bin = tmp_ice_ring - tmp_low_d_star_sq
      tmp_ice_ring_bin = (tmp_ice_ring_bin-tmp_step/2.0)/tmp_step
      tmp_ice_ring_bin = int(tmp_ice_ring_bin+0.5)+1
      if tmp_ice_ring_bin < 0 or tmp_ice_ring_bin > bin_centers.size() - 1:
        tmp_ice_ring_bin = None
      self.ice_ring_bin_location[count] = tmp_ice_ring_bin
      count+=1
      if tmp_ice_ring_bin is not None:
        weights[tmp_ice_ring_bin]=0.0
        ## also ignore flanking bins for safety
        if (tmp_ice_ring_bin-1) >= 0:
          weights[tmp_ice_ring_bin-1]=0.0
        if (tmp_ice_ring_bin+1) <=bin_centers.size()-1:
          weights[tmp_ice_ring_bin+1]=0.0
    mean_z_score = flex.sum(  weights*z_scores_data )
    mean_z_score /= flex.sum(weights)
    std_z_score = flex.sum(  weights*z_scores_data*z_scores_data )
    std_z_score /= flex.sum(weights)
    std_z_score = math.sqrt(std_z_score-mean_z_score*mean_z_score)
    mean_comp = flex.sum(  weights*completeness_data )
    mean_comp /= flex.sum(weights)
    std_comp = flex.sum(  weights*completeness_data*completeness_data)
    std_comp /= flex.sum(weights)
    std_comp = math.sqrt(std_comp-mean_comp*mean_comp)
    self.mean_comp=mean_comp
    self.mean_z_score= mean_z_score
    self.std_comp=std_comp
    self.std_z_score= std_z_score
    ## This array has z-score like features
    ## to detect ice rings (were are looking for spikes)
    self.abnormality_intensity = flex.double(10,0.0)
    self.value_intensity = flex.double(10,0.0)
    ## This array looks at the completeness
    ## and checks for sudden 'dips'
    self.abnormality_completeness = flex.double(10,0.0)
    self.value_completeness = flex.double(10,0.0)
    for ii in range(10):
      if self.ice_ring_bin_location[ii] is not None:
        ## checking is there is any weird, out of order z-score
        self.value_intensity[ii] = z_scores_data[
          self.ice_ring_bin_location[ii]]
        self.abnormality_intensity[ii]=z_scores_data[
          self.ice_ring_bin_location[ii]]
        self.abnormality_intensity[ii]-=mean_z_score
        self.abnormality_intensity[ii]/=(std_z_score+1.0e-6)
        ## checking for sudden dips in completeness
        self.value_completeness[ii]=completeness_data[
          self.ice_ring_bin_location[ii]]
        self.abnormality_completeness[ii]=completeness_data[
          self.ice_ring_bin_location[ii]]
        self.abnormality_completeness[ii]-=mean_comp
        self.abnormality_completeness[ii]/=(std_comp+1.0e-6)
    self.warnings=0 ## Number of ice ring related warnings
    icy_shells = []
    for ii in range(10):
      if self.ice_ring_bin_location[ii] is not None:
        icy_shells.append([ "%9.3f" % self.ice_d_spacings[ii],
                            "%10.3f" % abs(self.ice_rel_intens[ii]),
                            "%7.2f" % abs(self.value_intensity[ii]),
                            "%7.2f" % abs(self.value_completeness[ii]), ])
    self.table = table_utils.simple_table(
      table_rows=icy_shells,
      column_headers=["d_spacing", "Expected rel. I", "Data Z-score", "Completeness"])
    self.warnings = 0 ## Number of ice ring related warnings
    cutoff = self.completeness_abnormality_level
    for ii in range(10):
      if ((self.ice_ring_bin_location[ii] is not None) and
          (self.ice_rel_intens[ii] > self.intensity_level)):
        abnormality_completeness = abs(self.abnormality_completeness[ii])
        if ((abnormality_completeness >= cutoff) or
            (self.value_intensity[ii] > z_score_limit) or
            (abs(self.abnormality_intensity[ii]) >= cutoff)):
          self.warnings += 1

  def _show_impl(self, out):
    out.show_sub_header("Ice ring related problems")
    out.show("""\
 The following statistics were obtained from ice-ring insensitive resolution
 ranges:
""")
    out.show_preformatted_text("""\
    mean bin z_score      : %4.2f
        ( rms deviation   : %4.2f )
    mean bin completeness : %4.2f
        ( rms deviation   : %4.2f )
""" % (self.mean_z_score, self.std_z_score, self.mean_comp, self.std_comp))
    out.show("""\
 The following table shows the Wilson plot Z-scores and completeness for
 observed data in ice-ring sensitive areas.  The expected relative intensity
 is the theoretical intensity of crystalline ice at the given resolution.
 Large z-scores and high completeness in these resolution ranges might
 be a reason to re-assess your data processsing if ice rings were present.
""")
    out.show_table(self.table, indent=2)
    out.show("""\
 Abnormalities in mean intensity or completeness at resolution ranges with a
 relative ice ring intensity lower than %3.2f will be ignored.""" %
      (self.intensity_level))
    problems_detected = False
    cutoff = self.completeness_abnormality_level
    for ii in range(10):
      if ((self.ice_ring_bin_location[ii] is not None) and
          (self.ice_rel_intens[ii] > self.intensity_level)):
        abnormality_completeness = abs(self.abnormality_completeness[ii])
        if (abnormality_completeness >= cutoff):
          problems_detected = True
          out.show("""\
 At %3.2f A there is a lower completeness than expected from the rest of the
 data set.""" % self.ice_d_spacings[ii])
          if (abs(self.abnormality_intensity[ii]) >= cutoff):
            out.show("""\
 At the same resolution range, the expected mean intensity does not behave as
 it should.""")
          else :
            out.show("""\
 Even though the completeness is lower than expected, the mean intensity is
 still reasonable at this resolution.""")
        if ((abs(self.abnormality_intensity[ii]) >= cutoff) or
            (abs(self.value_intensity[ii])>getattr(self,"z_score_limit",10))):
          problems_detected = True
          if (abs(self.abnormality_completeness[ii])<=cutoff):
            out.show("""\
 At %3.2f A the z-score is more than %3.2f times the standard deviation of
 all z-scores, while at the same time, the completeness does not go down.
""" %(self.ice_d_spacings[ii], cutoff))

    if (not problems_detected):
      out.show("""\
 No ice ring related problems detected.
 If ice rings were present, the data does not look worse at ice ring related
 d_spacings as compared to the rest of the data set.""")
    if self.warnings == 1:
      out.show("""\
 As there was only 1 ice-ring related warning, it is not clear whether or not
 ice ring related features are really present.""")
    elif (self.warnings >= 2):
      out.show("""\
 There were %d ice ring related warnings.  This could indicate the presence of
 ice rings.""" % self.warnings)

  def summarize_issues(self):
    if (self.warnings > 0):
      return [(1, "The data appear to have one or more ice rings.",
        "Ice ring related problems")]
    else :
      return [(0, "Ice rings do not appear to be present.",
        "Ice ring related problems")]

#-----------------------------------------------------------------------
# ANOMALOUS MEASURABILITY
class analyze_measurability(scaling.xtriage_analysis):
  def __init__(self,
               d_star_sq,
               smooth_approx,
               meas_data,
               miller_array=None,
               low_level_cut=0.03,
               high_level_cut=0.06):
    low_level_range = smooth_approx > low_level_cut
    high_level_range = smooth_approx > high_level_cut
    tmp_low = d_star_sq.select(low_level_range)
    tmp_high = d_star_sq.select(high_level_range)
    self.low_d_cut = None
    self.high_d_cut = None
    ref_marks = [[None, None]]
    if tmp_low.size()>0:
      self.low_d_cut = flex.max(tmp_low) ** (-0.5)
      ref_marks[0][0] = flex.max(tmp_low)
    if tmp_high.size()>0:
      self.high_d_cut = flex.max(tmp_high) ** (-0.5)
      ref_marks[0][1] = flex.max(tmp_high)
    self.meas_table=None
    if miller_array is not None:
      work_array = miller_array.deep_copy()
      work_array.setup_binner(n_bins=10)
      self.meas_table = work_array.measurability(
        use_binning=True).as_simple_table(data_label="Measurability")
    self.table = data_plots.table_data(
      title="Measurability of Anomalous signal",
      column_labels=["1/resol**2", "Measurability", "Smooth approximation"],
      graph_names=["Anomalous measurability"],
      graph_labels=[("Resolution", "Measurability")],
      graph_columns=[[0,1,2]],
      data=[list(d_star_sq), list(meas_data), list(smooth_approx)],
      x_is_inverse_d_min=True,
      reference_marks=ref_marks,
      force_exact_x_labels=True)

  def _show_impl(self, out):
    out.show_sub_header("Measurability of anomalous signal")
    message = None
    unknown_cutoffs = ([self.low_d_cut,self.high_d_cut]).count(None)
    if unknown_cutoffs == 0 :
      message = None
      if self.low_d_cut ==  self.high_d_cut :
        message = """\
 The full resolution range seems to contain a useful amount of anomalous
 signal. Depending on your specific substructure, you could use all the data
 available for the location of the heavy atoms, or cut the resolution
 to speed up the search."""
      elif self.low_d_cut < self.high_d_cut:
        message = """\
 The anomalous signal seems to extend to about %3.1f A (or to %3.1f A, from a
 more optimistic point of view).  The quoted resolution limits can be used as a
 guideline to decide where to cut the resolution for phenix.hyss.""" % \
          (self.high_d_cut, self.low_d_cut)
        if self.high_d_cut < 3.0:
          message += """
 Depending however on the size and nature of your substructure you could cut
 the data at an even lower resolution to speed up the search."""
        elif self.high_d_cut > 4.5:
          message += """
 As the anomalous signal is not very strong in this dataset substructure
 solution via SAD might prove to be a challenge.  Especially if only low
 resolution reflections are used, the resulting substructures could contain a
 significant amount of false positives."""
    elif unknown_cutoffs == 2 :
      message = """\
 There seems to be no real significant anomalous differences in this dataset."""
    elif unknown_cutoffs == 1 :
      if self.high_d_cut is None:
        message="""\
 The anomalous signal seems to extend to about %3.1f A, but doesn't seem very
 strong. The quoted resolution limits can be used as a guideline to decide
 where to cut the resolution for phenix.hyss."""%(self.low_d_cut)
      else:
         message = """\
 This should not have happend: please contact software authors at
 help@phenix-online.org."""
    else :
      message = """\
 This should not have happend: please contact software authors at
 help@phenix-online.org."""
    if message is not None :
      out.show(message)
    # more details
    if (self.meas_table is not None):
      out.show("\nTable of measurability as a function of resolution:\n")
      out.show_table(self.meas_table, indent=2)
      out.show("""
 The measurability is defined as the fraction of Bijvoet related intensity
 differences for which the following relationship holds:""")
      out.show_preformatted_text("""
    |delta_I|/sigma_delta_I > 3.0
    min[I(+)/sigma_I(+), I(-)/sigma_I(-)] > 3.0
""")
      out.show("""\
  The measurability provides an intuitive feeling of the quality of the data,
  as it is related to the number of reliable Bijvoet differences.  When the
  data are processed properly and the standard deviations have been estimated
  accurately, values larger than 0.05 are encouraging.  Note that this
  analysis relies on the correctness of the estimated standard deviations of
  the intensities.
""")
      out.show_plot(self.table)

# TODO anomalous completeness

#-----------------------------------------------------------------------
# WRAPPER CLASSES

class data_strength_and_completeness(scaling.xtriage_analysis):
  """
  Collect basic info about overall completeness and signal-to-noise ratios,
  independent of scaling.
  """
  def __init__(self,
      miller_array,
      isigi_cut=3.0,
      completeness_cut=0.85,
      completeness_as_non_anomalous=None):
    # Make a deep copy not to upset or change things in
    # the binner that might be present
    tmp_miller = miller_array.deep_copy()
    tmp_miller.setup_binner_d_star_sq_step(auto_binning=True)
    self.d_star_sq_ori = tmp_miller.binner().bin_centers(2)
    # Signal-to-noise and completeness
    self.i_sig_i = None
    self.i_sig_i_table = None
    self.overall_i_sig_i = None
    if tmp_miller.sigmas() is not None:
      self.overall_i_sig_i = tmp_miller.i_over_sig_i(use_binning=False,
        return_fail=0)

      i_over_sigma = tmp_miller.i_over_sig_i(use_binning=True,
        return_fail=0)
      self.i_sig_i = i_over_sigma.data[1:len(i_over_sigma.data)-1]
      self.i_sig_i =  flex.double( self.i_sig_i )
      self.i_sig_i_table = data_plots.table_data(
        title="Signal to noise (<I/sigma_I>)",
        column_labels=["1/resol**2", "<I/sigma_I>"],
        graph_names=["Signal to noise"],
        graph_labels=[("Resolution", "<I/sigma_I>")],
        graph_columns=[[0,1]],
        data=[self.d_star_sq_ori,self.i_sig_i],
        x_is_inverse_d_min=True)
    #
    self.data_strength = None
    if (miller_array.sigmas() is not None):
      self.data_strength = i_sigi_completeness_stats(
        miller_array,
        isigi_cut=isigi_cut,
        completeness_cut=completeness_cut,
        completeness_as_non_anomalous=completeness_as_non_anomalous)
    self.completeness_overall = miller_array.completeness(
       as_non_anomalous_array = completeness_as_non_anomalous)
    # low-resolution completeness
    tmp_miller_lowres = miller_array.resolution_filter(d_min=5.0)
    self.low_resolution_completeness = None
    self.low_resolution_completeness_overall = None
    self.low_res_table = None
    if (tmp_miller_lowres.indices().size()>0):
      tmp_miller_lowres.setup_binner(n_bins=10)
      self.low_resolution_completeness_overall = \
        tmp_miller_lowres.completeness(use_binning=False,
          as_non_anomalous_array = completeness_as_non_anomalous)
      low_resolution_completeness = tmp_miller_lowres.completeness(
        use_binning=True,
        as_non_anomalous_array = completeness_as_non_anomalous)
      self.low_resolution_completeness = \
        low_resolution_completeness.as_simple_table(data_label="Completeness")
      binner = low_resolution_completeness.binner
      d_star_sq_ori = []
      comp = []
      for i_bin in binner.range_used():
        d_min = binner.bin_d_min(i_bin)
        d_star_sq_ori.append(1/(d_min**2))
        frac = low_resolution_completeness.data[i_bin]
        if (frac is None):
          frac = 0
        comp.append(frac*100)
      self.low_res_table = data_plots.table_data(
        title="Low-resolution completeness",
        column_labels=["Max. resolution", "Completeness"],
        graph_names=["Low-resolution completeness"],
        graph_labels=[("High resolution of shell", "% of total")],
        graph_columns=[[0,1]],
        data=[d_star_sq_ori, comp],
        x_is_inverse_d_min=True,
        force_exact_x_labels=True)
    try :
      self.d_min_directional = analyze_resolution_limits(miller_array)
    except AssertionError : # FIXME
      self.d_min_directional = None
    self.log_binned = log_binned_completeness(miller_array,
      completeness_as_non_anomalous=
        completeness_as_non_anomalous)

  def _show_impl(self, out):
    out.show_header("Data strength and completeness")

    if (hasattr(self, 'overall_i_sig_i') and
        (self.overall_i_sig_i is not None)):
      out.show("Overall <I/sigma> for this dataset is %7.1f" %(
        self.overall_i_sig_i))
    if (self.data_strength is not None):
      self.data_strength.show(out)
    if self.low_resolution_completeness is not None:
      out.show_sub_header("Low resolution completeness analyses")
      out.show("""\
The following table shows the completeness of the data to 5.0 A.  Poor
low-resolution completeness often leads to map distortions and other
difficulties, and is typically caused by problems with the crystal orientation
during data collection, overexposure of frames, interference with the beamstop,
or omission of reflections by data-processing software.""")
      out.show_table(self.low_resolution_completeness, indent=2)
    self.log_binned.show(out)
    if (self.d_min_directional is not None):
      self.d_min_directional.show(out)

  def high_resolution_for_twin_tests(self):
    if (self.data_strength is None):
      return None
    elif (self.data_strength.resolution_cut >
        self.data_strength.resolution_at_least):
      return self.data_strength.resolution_at_least
    else:
      return self.data_strength.resolution_cut

  def i_over_sigma_outer_shell(self):
    if (self.i_sig_i is not None):
      return self.i_sig_i[-1]
    return None

  def summarize_issues(self):
    issues = []
    if (self.d_min_directional is not None):
      if self.d_min_directional.is_elliptically_truncated():
        issues.append((1,"The data appear to have been elliptically truncated.",
          "Analysis of resolution limits"))
      else :
        issues.append((0,
          "The resolution cutoff appears to be similar in all directions.",
          "Analysis of resolution limits"))
    if self.low_resolution_completeness_overall and 0 < self.low_resolution_completeness_overall < 0.75:
      issues.append((2, "The overall completeness in low-resolution shells "+
        "is less than 90%.", "Low resolution completeness analyses"))
    elif self.low_resolution_completeness_overall and 0.75 <= self.low_resolution_completeness_overall < 0.9:
      issues.append((1, "The overall completeness in low-resolution shells "+
        "is less than 90%.", "Low resolution completeness analyses"))
    else :
      issues.append((0, "The overall completeness in low-resolution shells "+
        "is at least 90%.", "Low resolution completeness analyses"))
    return issues

class anomalous(scaling.xtriage_analysis):
  def __init__(self, miller_array, merging_stats=None,
      plan_sad_experiment_stats=None):
    assert miller_array.anomalous_flag()
    tmp_miller = miller_array.deep_copy()
    tmp_miller.setup_binner_d_star_sq_step(auto_binning=True)
    self.d_star_sq_ori = tmp_miller.binner().bin_centers(2)
    self.cc_anom_table = None
    if (merging_stats is not None):
      self.cc_anom_table = merging_stats.cc_anom_table
    self.plan_sad_experiment_stats=None
    if (plan_sad_experiment_stats is not None):
      self.plan_sad_experiment_stats=plan_sad_experiment_stats
    self.measurability = None
    # Get measurability if data is anomalous
    if (tmp_miller.sigmas() is not None):
      measurability = tmp_miller.measurability(use_binning=True,
        return_fail=0)
      meas_data = flex.double( measurability.data[1:len(
        measurability.data)-1])
      # make a smooth approximation to the measurability please
      smooth_meas_approx = chebyshev_lsq_fit.chebyshev_lsq_fit(
        int(self.d_star_sq_ori.size()/10) +3,
        self.d_star_sq_ori,
        meas_data )
      smooth_meas_approx = chebyshev_polynome(
        int(self.d_star_sq_ori.size()/10) +3,
        flex.min(self.d_star_sq_ori),
        flex.max(self.d_star_sq_ori),
        smooth_meas_approx.coefs)
      meas_smooth = smooth_meas_approx.f(self.d_star_sq_ori)
      self.measurability = analyze_measurability(
        d_star_sq=self.d_star_sq_ori,
        smooth_approx=meas_smooth,
        meas_data=meas_data,
        miller_array=tmp_miller)

  def _show_impl(self, out):
    out.show_header("Anomalous signal")
    if (self.cc_anom_table is not None):
      out.show_sub_header("Half-dataset anomalous correlation")
      out.show_text("""\
 This statistic (also called CC_anom) is calculated from unmerged data, which
 have been split into two half-datasets of approximately the same size and
 merged to obtain unique (anomalous) intensities, and their anomalous
 differences compared.  Resolution shells with CC_anom above 0.3 contain
 substantial anomalous signal useful for heavy atom location.""")
      if (out.gui_output):
        out.show_plot(self.cc_anom_table)
      else :
        out.show_table(self.cc_anom_table)
    if (self.measurability is not None):
      self.measurability.show(out)

    if (hasattr(self, 'plan_sad_experiment_stats') and
        (self.plan_sad_experiment_stats is not None)):
      out.show_header("Analysis of probability of finding %d sites with this anomalous data" %(
        self.plan_sad_experiment_stats.sites))
      if (out.gui_output):
        self.plan_sad_experiment_stats.show_in_wxgui(out=out)
      else:
        self.plan_sad_experiment_stats.show_summary(out=out)

  def summarize_issues(self):
    '''
    Traffic light
    '''
    issues = list()
    if (hasattr(self, 'plan_sad_experiment_stats') and
        (self.plan_sad_experiment_stats is not None)):
      p_substr = self.plan_sad_experiment_stats.get_p_substr()
      # round to nearest 5% when below 97%, otherwise round to nearest 1%
      if (p_substr < 97):
        p_substr = int(5.0*round(p_substr/5.0))
      else:
        p_substr = int(round(p_substr))
      sites = self.plan_sad_experiment_stats.sites
      message = 'The probability of finding %i sites with this data is around %i%%.' %\
                (sites, p_substr)
      section_name = 'Probability of finding sites and expected FOM'
      hi_limit = 85.0
      lo_limit = 50.0
      if (p_substr >= hi_limit):
        issues.append((0, message, section_name))
      elif ( (p_substr >= lo_limit) and (p_substr < hi_limit) ):
        issues.append((1, message, section_name))
      elif (p_substr < lo_limit):
        issues.append((2, message, section_name))
    return issues

class wilson_scaling(scaling.xtriage_analysis):
  """
  Calculates isotropic and anisotropic scale factors, Wilson plot, and various
  derived analyses such as ice rings and outliers.
  """
  def __init__(self,
               miller_array,
               n_residues,
               remove_aniso_final_b="eigen_min",
               use_b_iso=None,
               n_copies_solc=1,
               n_bases=0,
               z_score_cut=4.5,
               completeness_as_non_anomalous=None):
    assert (n_bases+n_residues != 0) and (n_copies_solc != 0)
    info = miller_array.info()
    tmp_miller = miller_array.deep_copy()
    tmp_miller.setup_binner_d_star_sq_step(auto_binning=True)
    self.d_star_sq_ori = tmp_miller.binner().bin_centers(2)
    if n_residues is None:
      n_residues = 0
    if n_bases is None:
      n_bases = 0
    if n_bases+n_residues==0:
      raise Sorry("No scatterers available")

    #-------------------------------------------------------------------
    # WILSON SCALING
    #
    # Isotropic
    order_z = miller_array.space_group().order_z()
    iso_scale_and_b = absolute_scaling.ml_iso_absolute_scaling(
      miller_array = miller_array,
      n_residues = n_residues * order_z * n_copies_solc,
      n_bases=n_bases * order_z * n_copies_solc)
    self.iso_scale_and_b = iso_scale_and_b
    ## Store the b and scale values from isotropic ML scaling
    self.iso_p_scale = iso_scale_and_b.p_scale
    self.iso_b_wilson =  iso_scale_and_b.b_wilson
    scat_info = self.iso_scale_and_b.scat_info

    ## Anisotropic ml wilson scaling
    order_z = miller_array.space_group().order_z()
    aniso_scale_and_b = absolute_scaling.ml_aniso_absolute_scaling(
      miller_array = miller_array,
      n_residues = n_residues * order_z * n_copies_solc,
      n_bases = n_bases * order_z * n_copies_solc)
    b_cart = aniso_scale_and_b.b_cart
    self.aniso_scale_and_b = aniso_scale_and_b

    self.aniso_p_scale = aniso_scale_and_b.p_scale
    self.aniso_u_star  = aniso_scale_and_b.u_star
    self.aniso_b_cart  = aniso_scale_and_b.b_cart
    # XXX: for GUI
    self.overall_b_cart = getattr(aniso_scale_and_b, "overall_b_cart", None)

    #-------------------------------------------------------------------
    # ANISOTROPY CORRECTION
    b_cart_observed = aniso_scale_and_b.b_cart
    b_trace_average = (b_cart_observed[0]+
                       b_cart_observed[1]+
                       b_cart_observed[2])/3.0
    b_trace_min = b_cart_observed[0]
    if  b_cart_observed[1] <b_trace_min: b_trace_min=b_cart_observed[1]
    if  b_cart_observed[2] <b_trace_min: b_trace_min=b_cart_observed[2]
    if (remove_aniso_final_b == "eigen_min"):
      b_use=aniso_scale_and_b.eigen_values[2]
    elif (remove_aniso_final_b == "eigen_mean"):
      b_use=flex.mean(aniso_scale_and_b.eigen_values)
    elif (remove_aniso_final_b == "user_b_iso"):
      assert remove_aniso_b_iso is not None
      b_use = use_b_iso
    else:
      b_use = 30

    b_cart_aniso_removed = [ -b_use, -b_use, -b_use, 0, 0, 0 ]
    u_star_aniso_removed = adptbx.u_cart_as_u_star(
      miller_array.unit_cell(),
      adptbx.b_as_u( b_cart_aniso_removed  ) )
    ## I do things in two steps, but can easely be done in 1 step
    ## just for clarity, thats all.
    self.no_aniso_array = absolute_scaling.anisotropic_correction(
      miller_array,0.0,aniso_scale_and_b.u_star )
    self.no_aniso_array = absolute_scaling.anisotropic_correction(
      self.no_aniso_array,0.0,u_star_aniso_removed)
    self.no_aniso_array = self.no_aniso_array.set_observation_type(
      miller_array )

    ## Make normalised structure factors please
    sel_big = self.no_aniso_array.data() > 1.e+50
    self.no_aniso_array = self.no_aniso_array.array(
      data = self.no_aniso_array.data().set_selected(sel_big, 0))
    self.no_aniso_array = self.no_aniso_array.set_observation_type(
      miller_array )
    normalistion = absolute_scaling.kernel_normalisation(
      self.no_aniso_array,auto_kernel=True)
    self.normalised_miller = normalistion.normalised_miller.deep_copy()

    # First we have to apply a resolution cut to make sure the reslution limts
    # match those of the empirical gamma array
    absolute_miller = miller_array.resolution_filter(
      d_max = math.sqrt(1.0/0.008),
      d_min = math.sqrt(1.0/0.69))
    #absolute_miller.as_mtz_dataset(column_root_label="F").mtz_object().write("start.mtz")
    ## anisotropy correction, bring data to absolute scale and B-value zero
    absolute_miller = absolute_scaling.anisotropic_correction(
      cache_0=absolute_miller,
      p_scale=self.aniso_p_scale,
      u_star=self.aniso_u_star)
    #absolute_miller.as_mtz_dataset(column_root_label="F").mtz_object().write("end.mtz")
    absolute_miller.set_observation_type( miller_array )
    # Now do some binning ala Popov&Bourenkov
    absolute_miller.setup_binner_d_star_sq_step(auto_binning=True)
    d_star_sq = absolute_miller.binner().bin_centers(2)
    d_star_sq[d_star_sq.size()-1] = 1.0/(
      flex.min( absolute_miller.d_spacings().data())
      *flex.min( absolute_miller.d_spacings().data()) )
    # Binning
    mean_observed_intensity = absolute_miller\
                              .mean_of_intensity_divided_by_epsilon(
                                 use_binning=True, return_fail=0.0)
    completeness = absolute_miller.completeness(use_binning=True,
        return_fail=1.0,
        as_non_anomalous_array = completeness_as_non_anomalous)
    # Recompute the scattering info summats please
    scat_info.scat_data(d_star_sq)
    # Please compute normalised structure factors
    normalisation = absolute_scaling.kernel_normalisation(
      absolute_miller,auto_kernel=True)
    normalised_miller = normalisation.normalised_miller.deep_copy()
    # set up a binner for this array as well please
    normalised_miller.setup_binner_d_star_sq_step(auto_binning=True)
    # Make a deep copy not to upset or change things in
    # the binner that might be present
    tmp_miller = miller_array.deep_copy()
    tmp_miller.setup_binner_d_star_sq_step(auto_binning=True)
    self.d_star_sq_ori = tmp_miller.binner().bin_centers(2)

    # Set up some arrays for plotting and analyses purposes
    self.d_star_sq = d_star_sq
    self.mean_I_normalisation = flex.exp(normalisation.normalizer.f(
      self.d_star_sq))
    self.mean_I_obs_data = flex.double(
      mean_observed_intensity.data[1:len(mean_observed_intensity.data)-1])
    # expected intensity
    theory = absolute_scaling.expected_intensity(
      scat_info,d_star_sq)
    self.mean_I_obs_theory = theory.mean_intensity
    # add standard deviations of experimental part
    # assuming wilson statistics for simplicity
    counts = flex.double( normalised_miller.binner().counts_given())
    counts = counts[1:counts.size()-1]
    self.mean_I_obs_sigma = self.mean_I_obs_data*self.mean_I_obs_data/(counts+1e-6)
    self.mean_I_obs_sigma+=theory.sigma_intensity*theory.sigma_intensity
    self.mean_I_obs_sigma=flex.sqrt(self.mean_I_obs_sigma)

    # z scores and completeness
    self.z_scores = flex.abs( self.mean_I_obs_data - self.mean_I_obs_theory )/\
                    self.mean_I_obs_sigma
    self.completeness = flex.double(completeness.data[
      1:len(mean_observed_intensity.data)-1])

    self.outliers = possible_outliers(absolute_miller)
    self.miller_array_filtered = self.outliers.remove_outliers(absolute_miller)
    self.ice_rings = None
    if (self.d_star_sq.size() > 1) and (flex.min( self.d_star_sq ) > 0.01):
      self.ice_rings = ice_ring_checker(
        bin_centers=self.d_star_sq,
        completeness_data=self.completeness,
        z_scores_data=self.z_scores)

    self.wilson_table = data_plots.table_data(
      title="Intensity plots",
      column_labels=["1/resol**2", "<I> smooth approximation",
                    "<I> via binning", "<I> expected"],
      graph_names=["Intensity plots"],
      graph_labels=[("Resolution", "<I>")],
      graph_columns=[[0,1,2,3]],
      data=[list(self.d_star_sq), list(self.mean_I_normalisation),
            list(self.mean_I_obs_data), list(self.mean_I_obs_theory)],
      x_is_inverse_d_min=True,
      force_exact_x_labels=True)
    # collect suspicious resolution shells
    worrisome = self.z_scores > z_score_cut
    self.n_worrisome = worrisome.count(True)
    self.outlier_shell_table = data_plots.table_data(
      title="Mean intensity by shell (outliers)",
      column_labels=["d_spacing", "z_score", "completeness", "<Iobs>/<Iexp>"],
      column_formats=["%9.3f", "%7.2f", "%7.2f", "%10.3f"],
      graph_names=["Worrisome resolution shells"],
      graph_columns=[[0,1,2,3]])
    for ii in range(self.d_star_sq.size()):
      if worrisome[ii]:
        d_space = self.d_star_sq[ii]**(-0.5)
        z_score = self.z_scores[ii]
        comp =  self.completeness[ii]
        ratio = self.mean_I_obs_data[ii] / self.mean_I_obs_theory[ii]
        self.outlier_shell_table.add_row([d_space,z_score,comp,ratio])
    ## z scores and completeness
    self.zscore_table = data_plots.table_data(
      title="Z scores and completeness",
      column_labels=["1/resol**2", "Z_score", "Fractional completeness"],
      graph_names=["Data sanity and completeness check"],
      graph_labels=[("Resolution", "Z score or fractional completeness")],
      graph_columns=[[0,1,2]],
      data=[list(self.d_star_sq), list(self.z_scores), list(self.completeness)],
      x_is_inverse_d_min=True)

  def show_worrisome_shells(self, out):
    out.show_sub_header("Mean intensity analyses")
    out.show("""\
 Inspired by: Morris et al. (2004). J. Synch. Rad.11, 56-59.
 The following resolution shells are worrisome:""")
    if (self.n_worrisome > 0):
      out.show_table(self.outlier_shell_table, indent=2)
      out.show("""\
 Possible reasons for the presence of the reported unexpected low or elevated
 mean intensity in a given resolution bin are :
   - missing overloaded or weak reflections
   - suboptimal data processing
   - satellite (ice) crystals
   - NCS
   - translational pseudo symmetry (detected elsewhere)
   - outliers (detected elsewhere)
   - ice rings (detected elsewhere)
   - other problems
 Note that the presence of abnormalities in a certain region of reciprocal
 space might confuse the data validation algorithm throughout a large region
 of reciprocal space, even though the data are acceptable in those areas.

""")
    else :
      out.show(" *** None ***")

  def _show_impl(self, out):
    out.show_header("Absolute scaling and Wilson analysis")
    self.iso_scale_and_b.show(out=out)
    self.aniso_scale_and_b.show(out=out)
    out.show_sub_header("Wilson plot")
    # FIXME get feedback from TPTB
    out.show_text("""\
 The Wilson plot shows the falloff in intensity as a function in resolution;
 this is used to calculate the overall B-factor ("Wilson B-factor") for the
 data shown above.  The expected plot is calculated based on analysis of
 macromolecule structures in the PDB, and the distinctive appearance is due to
 the non-random arrangement of atoms in the crystal.  Some variation is
 natural, but major deviations from the expected plot may indicate pathological
 data (including ice rings, detector problems, or processing errors).""")
    out.show_plot(self.wilson_table)
    # XXX is this really necessary?
    #out.show_plot(self.zscore_table)
    self.show_worrisome_shells(out)
    self.outliers.show(out)
    if self.ice_rings is not None:
      self.ice_rings.show(out)

  def summarize_issues(self):
    issues = []
    if self.ice_rings is not None:
      issues.extend(self.ice_rings.summarize_issues())
    issues.extend(self.outliers.summarize_issues())
    issues.extend(self.iso_scale_and_b.summarize_issues())
    issues.extend(self.aniso_scale_and_b.summarize_issues())
    return issues

  # Objects of this class are relatively bulky due to the storage of multiple
  # derived Miller arrays that may be used elsewhere.  Since we do not need
  # these arrays for simply displaying the results (in the Phenix GUI or
  # elsewhere), they are deleted prior to pickling to reduce the amount of
  # data that needs to be transfered or saved.  It is not necessary to
  # implement __setstate__, since we are still just pickling self.__dict__.
  def __getstate__(self):
    """
    Pickling function with storage efficiency optimizations.
    """
    self.miller_array_filtered = None
    self.no_aniso_array = None
    self.normalised_miller = None
    self.d_star_sq_ori = None
    self.iso_scale_and_b.work_array = None
    self.aniso_scale_and_b.work_array = None
    self.iso_scale_and_b.scat_info = None
    self.aniso_scale_and_b.scat_info = None
    return self.__dict__
