from __future__ import division
from cctbx import miller
from cctbx.array_family import flex
from libtbx.utils import Sorry
from iotbx import data_plots
from mmtbx.scaling import absolute_scaling
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
from scitbx.math import erf
from libtbx import table_utils
import sys
import math


class i_sigi_completeness_stats(object):
  def __init__(self,
               miller,
               n_bins=15,
               isigi_cut=3.0,
               completeness_cut=0.85,
               resolution_at_least=3.5):
    self.miller = miller.deep_copy().set_observation_type(miller)
    #check to see if we have intensities
    if self.miller.is_real_array():
      if self.miller.is_xray_amplitude_array():
        self.miller=self.miller.f_as_f_sq()
    assert self.miller.is_xray_intensity_array()
    #make sure we have sigmas
    assert self.miller.sigmas() is not None
    # select things with sigma larger then zero
    #self.miller = self.miller.select( self.miller.sigmas() > 0 )
    self.miller.setup_binner(n_bins=n_bins)

    self.resolution_bins = list(self.miller.binner().limits())
    self.completeness_bins = []

    self.table = None
    self.make_table()

    self.isigi_cut=isigi_cut
    self.completeness_cut=completeness_cut
    self.resolution_at_least=resolution_at_least
    self.resolution_cut = self.guess_resolution_cut(
      isigi_cut=isigi_cut,
      completeness_cut=completeness_cut,
      resolution_at_least=resolution_at_least)


  def cut_completeness(self,cut_value):
    tmp_miller = self.miller.select(
      self.miller.data() > cut_value*self.miller.sigmas() )
    tmp_miller.use_binning_of(self.miller)
    completeness = tmp_miller.completeness(use_binning=True,
                                           return_fail=0.0)
    return completeness.data


  def guess_resolution_cut(self,
                           isigi_cut=3.0,
                           completeness_cut=0.85,
                           resolution_at_least=3.5):
    comp_data = self.cut_completeness(isigi_cut)
    reso=4.0
    for ii in xrange(1,len(self.resolution_bins)-1):
      a = self.resolution_bins[ii-1]
      b = self.resolution_bins[ii]
      if b < resolution_at_least:
        if comp_data[ii]>completeness_cut:
          reso = b**-0.5
    return reso

  def make_table(self):
    table_data = []
    cut_list=[1,2,3,5,10,15]
    for cut_value in cut_list:
      self.completeness_bins.append( self.cut_completeness(cut_value) )

    legend = ("Res. Range", "I/sigI>1 ",
                            "I/sigI>2 ",
                            "I/sigI>3 ",
                            "I/sigI>5 ",
                            "I/sigI>10",
                            "I/sigI>15" )
    # XXX: for GUI
    self.table_for_gui = data_plots.table_data(
      title="Completeness and data strength",
      column_labels=["Max. resolution"] + list(legend)[1:],
      graph_names=["I/sigI by shell"],
      graph_labels=[("High resolution of shell", "% of total")],
      graph_columns=[list(range(7))],
      x_is_inverse_d_min=True)
    for ii in xrange(1,len(self.resolution_bins)-1):
      row = []
      raw_row = []
      a = self.resolution_bins[ii-1]
      b = self.resolution_bins[ii]
      limsa =("%4.2f"%(a**-0.5)).rjust(5)
      limsb =("%4.2f"%(b**-0.5)).rjust(5)

      lims = limsa+" -"+limsb
      row.append( lims )
      raw_row.append(b)
      for jj in  self.completeness_bins:
        raw_row.append(100.0*jj[ii])
        tmp="%3.1f%s"%(100.0*jj[ii],"%")
        row.append( tmp )
      table_data.append( row )
      self.table_for_gui.add_row(raw_row)

    self.table = table_utils.format([legend]+table_data,
                                    has_header=True,
                                    separate_rows=False,
                                    prefix='| ',
                                    postfix=' |')

  def show(self, out=None):
    print >> out, "Completeness and data strength analyses"
    self.completeness_info = """
  The following table lists the completeness in various resolution
  ranges, after applying a I/sigI cut. Miller indices for which
  individual I/sigI values are larger than the value specified in
  the top row of the table, are retained, while other intensities
  are discarded. The resulting completeness profiles are an indication
  of the strength of the data.

"""
    print >> out, self.completeness_info
    if out is None:
      out = sys.stdout
    print >> out, self.table
    print >> out
    message = """\
  The completeness of data for which I/sig(I)>%3.2f, exceeds %3.0f%s
  for resolution ranges lower than %3.2fA.
""" % (self.isigi_cut, self.completeness_cut*100, "%", self.resolution_cut)

    if self.resolution_cut < self.resolution_at_least:
      message += """\
  The data are cut at this resolution for the potential twin tests
  and intensity statistics.
"""
    else:
      message += """\
  As we do not want to throw away too much data, the resolution for
  analyzing the intensity statistics will be limited to %3.2fA.
""" % self.resolution_at_least
    print >> out, message
    self.completeness_info += "\n" + message


class completeness_enforcement(object):
  def __init__(self,
               miller_array,
               minimum_completeness=0.75):

    self.miller_array = miller_array.deep_copy()
    self.miller_array.setup_binner_d_star_sq_step(auto_binning=True)
    completeness = self.miller_array.completeness(use_binning=True,
                                                  return_fail=1.0)

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






class possible_outliers(object):
  def __init__(self,
               miller_array,
               prob_cut_ex=1.0E-1,
               prob_cut_wil=1.0E-6):

    if miller_array.observation_type() is None:
      raise Sorry("Unknown observation type")

    if not miller_array.is_xray_intensity_array():
      work_array = miller_array.f_as_f_sq()
    else:
      work_array = miller_array.deep_copy()


    work_array = miller_array.deep_copy()
    normalizer = absolute_scaling.kernel_normalisation(
      work_array, auto_kernel=True)
    work_array = work_array.array(
      data=normalizer.normalised_miller.data()
      /work_array.epsilons().data().as_double())

    self.centric_cut = work_array.select_centric().set_observation_type(
      work_array)
    self.acentric_cut = work_array.select_acentric().set_observation_type(
      work_array)

    self.p_acentric_single = 1.0 - (1.0 - flex.exp(-self.acentric_cut.data() ))
    self.p_centric_single = 1.0 - erf(flex.sqrt(self.centric_cut.data()/2.0) )

    n_centric = self.p_centric_single.size()
    n_acentric = self.p_acentric_single.size()

    self.extreme_acentric = 1.0 -  \
       flex.pow(1.0 - flex.exp(-self.acentric_cut.data() ),float(n_acentric))

    self.extreme_centric = 1.0 - \
       flex.pow(erf(flex.sqrt(self.centric_cut.data()/2.0) ),float(n_centric))

    ## combine both the wilson and extreme value cut-off values
    acentric_outlier = (self.extreme_acentric < prob_cut_ex) or (
     self.p_acentric_single < prob_cut_wil)
    centric_outlier = (self.extreme_centric < prob_cut_ex) or (
     self.p_centric_single  < prob_cut_wil)


    ## acentrics
    self.acentric_outlier_miller = self.acentric_cut.indices().select(
      acentric_outlier)
    self.acentric_outlier_e_vals = self.acentric_cut.data().select(
      acentric_outlier)
    self.acentric_outlier_e_vals = flex.sqrt(self.acentric_outlier_e_vals)
    self.acentric_outlier_d_spacings = self.acentric_cut.d_spacings().data()\
                                       .select(acentric_outlier)
    self.acentric_outlier_p_val = self.p_acentric_single.select(
      acentric_outlier)
    self.acentric_outlier_extreme_val = self.extreme_acentric.select(
      acentric_outlier)

    ## centrics
    self.centric_outlier_miller = self.centric_cut.indices().select(
      centric_outlier)
    self.centric_outlier_e_vals = self.centric_cut.data().select(
      centric_outlier)
    self.centric_outlier_e_vals = flex.sqrt(self.centric_outlier_e_vals)
    self.centric_outlier_d_spacings = self.centric_cut.d_spacings().data()\
                                       .select(centric_outlier)
    self.centric_outlier_p_val = self.p_centric_single.select(
      centric_outlier)
    self.centric_outlier_extreme_val = self.extreme_centric.select(
      centric_outlier)

  def remove_outliers(self, miller_array):
    ## remove the outliers please
    centric_matches = miller.match_indices( miller_array.indices(),
                                     self.centric_outlier_miller )
    miller_array = miller_array.select( centric_matches.single_selection(0) )

    acentric_matches = miller.match_indices( miller_array.indices(),
                                      self.acentric_outlier_miller )
    miller_array = miller_array.select( acentric_matches.single_selection(0) )
    return miller_array


  def show(self,out=None):
    if out is None:
      out=sys.stdout

    print >> out, "Possible outliers "
    print >> out, "  Inspired by: Read, Acta Cryst. (1999). D55, 1759-1764"
    print >> out
    print >> out, " Acentric reflections:"
    print >> out
    self.acentric_outliers_table = data_plots.table_data(
      title="Acentric reflections",
      column_labels=["d_spacing", "H K L", "|E|", "p(wilson)", "p(extreme)"],
      graph_names=["Possible acentric outliers"],
      graph_columns=[[0,1,2,3,4]])
    if self.acentric_outlier_miller.size() ==0:
      print >> out, "            None "
      print >> out
    else:
      print >> out, "-----------------------------------------------------------------"
      print >> out,"| d_space |      H     K     L |  |E|  | p(wilson) | p(extreme) |"
      print >> out,"-----------------------------------------------------------------"
      for d,hkl,e,p,extr in zip(self.acentric_outlier_d_spacings,
                                self.acentric_outlier_miller,
                                self.acentric_outlier_e_vals,
                                self.acentric_outlier_p_val,
                                self.acentric_outlier_extreme_val):
        self.acentric_outliers_table.add_row([d, hkl, e, p, extr])
        h = hkl[0]
        k = hkl[1]
        l = hkl[2]
        hkl = str(hkl)
        print >> out,"|%8.3f |  %5i,%5i,%5i |%6.2f | %9.2e | %10.2e |"%(d,h,k,l,e,p,extr)
      print >> out,"-----------------------------------------------------------------"
      print >> out
      print >> out, " p(wilson)  : 1-(1-exp[-|E|^2]) "
      print >> out, " p(extreme) : 1-(1-exp[-|E|^2])^(n_acentrics)"
      print >> out, " p(wilson) is the probability that an E-value of the specified"
      print >> out, " value would be observed if it were selected at random"
      print >> out, " the given data set."
      print >> out, " p(extreme) is the probability that the largest |E| value is "
      print >> out, " larger or equal than the observed largest |E| value."
      print >> out
      print >> out, " Both measures can be used for outlier detection. p(extreme)"
      print >> out, " takes into account the size of the dataset."
      print >> out
      print >> out

    print >> out, " Centric reflections:"
    print >> out
    self.centric_outliers_table = data_plots.table_data(
      title="Centric reflections",
      column_labels=["d_spacing", "H K L", "|E|", "p(wilson)", "p(extreme)"],
      graph_names=["Possible centric outliers"],
      graph_columns=[[0,1,2,3,4]])
    if self.centric_outlier_miller.size() ==0:
      print >> out, "            None "
      print >> out
    else:
      print >> out, "-----------------------------------------------------------------"
      print >> out,"| d_space |      H     K     L |  |E|  | p(wilson) | p(extreme) |"
      print >> out,"-----------------------------------------------------------------"
      for d,hkl,e,p,extr in zip(self.centric_outlier_d_spacings,
                                self.centric_outlier_miller,
                                self.centric_outlier_e_vals,
                                self.centric_outlier_p_val,
                                self.centric_outlier_extreme_val):
        self.centric_outliers_table.add_row([d, hkl, e, p, extr])
        h = hkl[0]
        k = hkl[1]
        l = hkl[2]
        hkl = str(hkl)
        print >> out,"|%8.3f |  %5i,%5i,%5i |%6.2f | %9.2e | %10.2e |"%(d,h,k,l,e,p,extr)
      print >> out,"-----------------------------------------------------------------"
      print >> out
      print >> out, " p(wilson)  : 1-(erf[|E|/sqrt(2)]) "
      print >> out, " p(extreme) : 1-(erf[|E|/sqrt(2)])^(n_acentrics)"
      print >> out, " p(wilson) is the probability that an E-value of the specified"
      print >> out, " value would be observed when it would selected at random from"
      print >> out, " the given data set."
      print >> out, " p(extreme) is the probability that the largest |E| value is "
      print >> out, " larger or equal than the observed largest |E| value."
      print >> out
      print >> out, " Both measures can be used for outlier detection. p(extreme)"
      print >> out, " takes into account the size of the dataset."


class ice_ring_checker(object):
  def __init__(self,
               bin_centers,
               completeness_data,
               z_scores_data):
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
      if tmp_ice_ring_bin < 0:
        tmp_ice_ring_bin = None
      if tmp_ice_ring_bin > bin_centers.size()-1:
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
        ## checking is there is any weerd, out of order z-score
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


  def show(self, out=None, verbose=0,level=4.0,intensity_level=0.1):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out, "Ice ring related problems"
    print >> out
    print >> out, " The following statistics were obtained from ice-ring "
    print >> out, " insensitive resolution ranges "
    print >> out, "  mean bin z_score      : %4.2f"%(self.mean_z_score)
    print >> out, "      ( rms deviation   : %4.2f )"%(self.std_z_score)
    print >> out, "  mean bin completeness : %4.2f"%(self.mean_comp)
    print >> out, "     ( rms deviation   : %4.2f )"%(self.std_comp)
    print >> out
    print >> out, " The following table shows the z-scores "
    print >> out, " and completeness in ice-ring sensitive areas."
    print >> out, " Large z-scores and high completeness in these "
    print >> out, " resolution ranges might be a reason to re-assess"
    print >> out, " your data processsing if ice rings were present."
    print >> out
    print >> out, "------------------------------------------------"
    print >> out, "| d_spacing | z_score | compl. | Rel. Ice int. |"
    print >> out, "------------------------------------------------"
    self.icy_shells = []
    for ii in range(10):
      if self.ice_ring_bin_location[ii] is not None:

        print >> \
          out, "|%9.3f  |%7.2f  |%7.2f |%10.3f     |"%(
          (self.ice_d_spacings[ii]),
          abs(self.value_intensity[ii]),
          abs(self.value_completeness[ii]),
          abs(self.ice_rel_intens[ii]))
        self.icy_shells.append([ self.ice_d_spacings[ii],
                                 abs(self.value_intensity[ii]),
                                 abs(self.value_completeness[ii]),
                                 abs(self.ice_rel_intens[ii]) ])
    print >> out, "------------------------------------------------"
    print >> out
    print >> out, " Abnormalities in mean intensity or completeness at"
    print >> out, " resolution ranges with a relative ice ring intensity"
    print >> out, " lower than %3.2f will be ignored."%(intensity_level)
    print >> out
    self.comments = False
    self.message = ""
    for ii in range(10):
      if (self.ice_ring_bin_location[ii] is not None) \
       and (self.ice_rel_intens[ii]>intensity_level):
        if (abs(self.abnormality_completeness[ii])>=level):
          self.comments = True
          self.warnings+=1
          self.message += """\
 At %3.2f A there is an lower occupancy
 than expected from the rest of the data set.""" % self.ice_d_spacings[ii]
          if (abs(self.abnormality_intensity[ii])>=level):
            self.message += """\
 At the same resolution range, the expected
 mean intensity does not behave as it should.
"""
          elif (abs(self.abnormality_intensity[ii])<level):
            self.message += """\
 Even though the completeness is lower as expected,
 the mean intensity is still reasonable at this resolution.
"""
        if (abs(self.abnormality_intensity[ii])>=level):
          self.comments=True
          self.warnings+=1
          if (abs(self.abnormality_completeness[ii])<=level):
            self.message += """\
 At %3.2f A the z-score is more than %3.2f times the standard
 deviation of all z-scores, while at the same time, the
 occupancy does not go down.
""" %(self.ice_d_spacings[ii],level)

    if not self.comments:
      self.message += """\
 No ice ring related problems detected.
 If ice rings were present, the data does not look
 worse at ice ring related d_spacings as compared
 to the rest of the data set."""
    if self.warnings==1:
      self.message += """\
 As there was only 1 ice-ring related warning, it is not
 clear whether or not ice ring related features are
 really present.
"""
    elif self.warnings>=2:
      self.message += """\
 There were %2.0f ice ring related warnings
 This could indicate the presence of ice rings.
""" % self.warnings

    if self.message != "" :
      print >> out, self.message
    print >> out
    print >> out
    print >> out
    print >> out

class analyze_measurability(object):
  def __init__(self,
               d_star_sq,
               smooth_approx,
               miller_array=None,
               low_level_cut=0.03,
               high_level_cut=0.06):
    low_level_range = smooth_approx > low_level_cut
    high_level_range = smooth_approx > high_level_cut
    tmp_low = d_star_sq.select(low_level_range)
    tmp_high = d_star_sq.select(high_level_range)

    self.low_d_cut = None
    self.high_d_cut = None
    if tmp_low.size()>0:
      self.low_d_cut = flex.max(tmp_low)
      self.low_d_cut = self.low_d_cut**(-0.5)
    if tmp_high.size()>0:
      self.high_d_cut = flex.max(tmp_high)
      self.high_d_cut = self.high_d_cut**(-0.5)
    self.meas_table=None

    if miller_array is not None:
      work_array = miller_array.deep_copy()
      work_array.setup_binner(n_bins=10)
      self.meas_table = work_array.measurability(use_binning=True)

  def show(self, out=None):
    if out is None:
      out = sys.stdout
    print >> out, "Analyses of anomalous differences"
    print >> out
    if self.meas_table is not None:
      print >> out, """\
  Table of measurability as a function of resolution

  The measurability is defined as the fraction of
  Bijvoet related intensity differences for which
  |delta_I|/sigma_delta_I > 3.0
  min[I(+)/sigma_I(+), I(-)/sigma_I(-)] > 3.0
  holds.
  The measurability provides an intuitive feeling
  of the quality of the data, as it is related to the
  number of reliable Bijvoet differences.
  When the data are processed properly and the standard
  deviations have been estimated accurately, values larger
  than 0.05 are encouraging.
  Note that this analyses relies on the correctness of the estimated
  standard deviations of the intensities.

"""
      self.meas_table.show(f=out)
      print >> out

    # XXX: refactored for GUI
    message = None
    unknown_cutoffs = ([self.low_d_cut,self.high_d_cut]).count(None)
    if unknown_cutoffs == 0 :
      message = None
      if self.low_d_cut ==  self.high_d_cut :
        message = """\
 The full resolution range seems to contain a useful
 ammount of anomalous signal. Depending on your
 specific substructure, you could use all the data available
 for the location of the heavy atoms, or cut the resolution
 to speed up the search."""
      elif self.low_d_cut < self.high_d_cut:
        message = """\
 The anomalous signal seems to extend to about %3.1f A
 (or to %3.1f A, from a more optimistic point of view).
 The quoted resolution limits can be used as a guideline
 to decide where to cut the resolution for phenix.hyss.""" % \
          (self.high_d_cut, self.low_d_cut)
        if self.high_d_cut < 3.0:
          message += """\
 Depending however on the size and nature of your substructure
 you could cut the data at an even lower resolution to speed up
 the search."""
        elif self.high_d_cut > 4.5:
          message += """\
 As the anomalous signal is not very strong in this dataset
 substructure solution via SAD might prove to be a challenge.
 Especially if only low resolution reflections are used,
 the resulting substructures could contain a significant amount
 of false positives."""
    elif unknown_cutoffs == 2 :
      message = """\
 There seems to be no real significant anomalous differences
 in this dataset."""
    elif unknown_cutoffs == 1 :
      if self.high_d_cut is None:
        message="""\
The anomalous signal seems to extend to about %3.1f A, but doesn't seem very strong.
The quoted resolution limits can be used as a guideline
to decide where to cut the resolution for phenix.hyss."""%(self.low_d_cut)
      else:
         message = """\
 This should not have happend: please contact software authors
 at bugs@phenix-online.org."""
    else :
      message = """\
 This should not have happend: please contact software authors
 at bugs@phenix-online.org.

"""
    if message is not None :
      print >> out, message
    self.message = message



class i_over_sigma_and_completeness(object):
  def __init__(self,
               miller_array,
               n_bins_table=15):
    self.miller_array = miller_array.deep_copy()
    self.resolution_bins =None
    self.completeness_array = []
    self.i_sigi_array = [1.0,2.0,3.0,5.0,10.0,15.0]
    # make sure we actually have sigmas
    if self.miller_array.sigmas() is not None:
      self.miller_array = self.miller_array.select(
        self.miller_array.sigmas()> 0 )
      # I would like to have intensities
      if not self.miller_array.is_xray_intensity_array():
        if self.miller_array.is_real_array():
          self.miller_array = self.miller_array.f_as_f_sq()

      self.miller_array.setup_binner(n_bins=n_bins_table)

      for i_i_sigi in xrange(len(self.i_sigi_array)):
        # please select all reflection for which i sigi is larger than the given number
        tmp_cut_off = self.i_sigi_array[ i_i_sigi ]
        tmp_miller_array = self.miller_array.select(
          (self.miller_array.data() > tmp_cut_off*self.miller_array.sigmas())
        )
        tmp_miller_array.use_binning_of( self.miller_array )
        # now please get the completeness
        compl =  tmp_miller_array.completeness(use_binning=True,
                                               return_fail=0.0)
        tmp_data = flex.double(compl.data)

        tmp_bins = compl.binner.limits()
        self.completeness_array.append( tmp_data )
        self.resolution_bins = list(compl.binner.limits())

      self.table=None
      self.make_table()

  def make_table(self):

    table_data = []
    legend = ("Res. Range", "I/sigI>1 ",
                            "I/sigI>2 ",
                            "I/sigI>3 ",
                            "I/sigI>5 ",
                            "I/sigI>10",
                            "I/sigI>15" )

    for ii in xrange(1,len(self.resolution_bins)-1):
      row = []
      a = self.resolution_bins[ii-1]
      b = self.resolution_bins[ii]
      limsa =("%4.2f"%(a**-0.5)).rjust(5)
      limsb =("%4.2f"%(b**-0.5)).rjust(5)

      lims = limsa+" -"+limsb
      row.append( lims )
      for jj in  self.completeness_array:
        tmp="%3.1f%s"%(100.0*jj[ii],"%")
        row.append( tmp )
      table_data.append( row )

    self.table = table_utils.format([legend]+table_data,
                                    has_header=True,
                                    separate_rows=False,
                                    prefix='| ',
                                    postfix=' |')


  def show(self, out=None):
    print >> out
    print >> out
    print >> out, "Completeness and data strength analyses "
    print >> out
    print >> out, "  The following table lists the completeness in various resolution"
    print >> out, "  ranges, after applying a I/sigI cut. Miller indices for which"
    print >> out, "  individual I/sigI values are larger than the value specified in"
    print >> out, "  the top row of the table, are retained, while other intensities"
    print >> out, "  are discarded. The resulting completeness profiles are an indication"
    print >> out, "  of the strength of the data."
    print >> out
    if out is None:
      out = sys.stdout
    print >> out, self.table









class basic_intensity_statistics:
  def __init__(self,
               miller_array,
               p_scale,
               u_star_tensor,
               scat_info,
               out=None,
               out_plot=None,
               verbose=0):

    if out is None:
      out=sys.stdout

    ## First we have to apply a resolution cut
    ## to make sure the reslution limts match those
    ## of the empirical gamma array
    absolute_miller = miller_array.resolution_filter(
      d_max = math.sqrt(1.0/0.008),
      d_min = math.sqrt(1.0/0.69))

    ## anisotropy correction, binrg data to absolute scale and B-value zero
    absolute_miller = absolute_scaling.anisotropic_correction(
      absolute_miller,
      p_scale,
      u_star_tensor)
    absolute_miller.set_observation_type( miller_array )
    ## Now do some binning ala Popov&Bourenkov
    absolute_miller.setup_binner_d_star_sq_step(auto_binning=True,
                                                )
    d_star_sq = absolute_miller.binner().bin_centers(2)


    d_star_sq[d_star_sq.size()-1] = 1.0/(
      flex.min( absolute_miller.d_spacings().data())
      *flex.min( absolute_miller.d_spacings().data()) )

    ## Binning
    mean_observed_intensity = absolute_miller\
                              .mean_of_intensity_divided_by_epsilon(
                                 use_binning=True, return_fail=0.0)
    completeness = absolute_miller.completeness(use_binning=True,
                                                return_fail=1.0)
    ## Recompute the scattering info summats please
    scat_info.scat_data(d_star_sq)

    ## Please compute normalised structure factors
    normalisation = absolute_scaling.kernel_normalisation(
      absolute_miller,auto_kernel=True)

    normalised_miller = normalisation.normalised_miller.deep_copy()
    ## set up a binner for this array as well please
    normalised_miller.setup_binner_d_star_sq_step(auto_binning=True)

    ## Make a deep copy not to upset or change things in
    ## the binner that might be present
    tmp_miller = miller_array.deep_copy()
    tmp_miller.setup_binner_d_star_sq_step(auto_binning=True)
    self.d_star_sq_ori = tmp_miller.binner().bin_centers(2)
    self.meas_data=None
    ## Get measurability if data is anomalous
    if tmp_miller.anomalous_flag():
      if tmp_miller.sigmas() is not None:
        measurability = tmp_miller.measurability(use_binning=True,return_fail=0)
        self.meas_data = flex.double( measurability.data[1:len(
          measurability.data)-1])
        ## make a smooth approximation to the measurability please
        smooth_meas_approx = chebyshev_lsq_fit.chebyshev_lsq_fit(
          int(self.d_star_sq_ori.size()/10) +3,
          self.d_star_sq_ori,
          self.meas_data )
        smooth_meas_approx = chebyshev_polynome(
          int(self.d_star_sq_ori.size()/10) +3,
          flex.min(self.d_star_sq_ori),
          flex.max(self.d_star_sq_ori),
          smooth_meas_approx.coefs)
        self.meas_smooth = smooth_meas_approx.f(self.d_star_sq_ori)

        self.meas_anal =analyze_measurability(
          self.d_star_sq_ori,
          self.meas_smooth,
          tmp_miller)


    ## Set up some arrays for plotting and analyses purposes
    self.d_star_sq = d_star_sq

    self.mean_I_normalisation = flex.exp(normalisation.normalizer.f(self.d_star_sq))
    self.mean_I_obs_data = flex.double(
      mean_observed_intensity.data[1:len(mean_observed_intensity.data)-1])

    ## expected intensity
    theory = absolute_scaling.expected_intensity(
      scat_info,d_star_sq)
    self.mean_I_obs_theory = theory.mean_intensity
    ## add standard deviations of experimental part
    ## assuming wilson statistics for simplicity
    counts = flex.double( normalised_miller.binner().counts_given())
    counts = counts[1:counts.size()-1]
    self.mean_I_obs_sigma=self.mean_I_obs_data*self.mean_I_obs_data/(counts+1e-6)
    self.mean_I_obs_sigma+=theory.sigma_intensity*theory.sigma_intensity
    self.mean_I_obs_sigma=flex.sqrt(self.mean_I_obs_sigma)

    ## z scores and completeness
    self.z_scores = flex.abs( self.mean_I_obs_data - self.mean_I_obs_theory )/\
                    self.mean_I_obs_sigma
    self.completeness = flex.double(completeness.data[
      1:len(mean_observed_intensity.data)-1])

    ## I over sigma
    self.i_sig_i = None
    if tmp_miller.sigmas() is not None:
      i_over_sigma = tmp_miller.i_over_sig_i(use_binning=True,return_fail=0)
      self.i_sig_i = i_over_sigma.data[1:len(i_over_sigma.data)-1]
      self.i_sig_i =  flex.double( self.i_sig_i )



    self.low_reso_completeness=None
    self.low_reso_meas=None
    self.low_reso_strong_ano=None
    self.suggested_reso_for_hyss=None







    ## For the log file, not plots
    ##
    ##
    tmp_miller = miller_array.resolution_filter(d_min=5.0)
    self.low_resolution_completeness = None
    if (tmp_miller.indices().size()>0):
      tmp_miller.setup_binner(n_bins=10)
      self.low_resolution_completeness = tmp_miller.completeness(
        use_binning=True)




    self.outlier = possible_outliers(absolute_miller)
    self.new_miller = self.outlier.remove_outliers(absolute_miller)

    self.ijsco = None
    if(self.d_star_sq.size()>1):
      if flex.min( self.d_star_sq ) > 0.01:
        self.ijsco = ice_ring_checker(self.d_star_sq,
                                      self.completeness,
                                      self.z_scores)
    else:
      print >> out, " *** WARNING *** "*5
      print >> out, "   Skip ice_ring_checker"
      print >> out, " *** WARNING *** "*5
    self.show(out,out_plot)

  def show(self,out=None, out_plot=None,z_level=4.5):
    if out is None:
      out = sys.stdout

    ## First talk about completeness
    print >> out
    print >> out
    if self.low_resolution_completeness is not None:
      print >> out, "Low resolution completeness analyses "
      print >> out
      print >> out, " The following table shows the completeness"
      print >> out, " of the data to 5 Angstrom."
      self.low_resolution_completeness.show(f=out)
      print >> out
      print >> out
      print >> out
    ## Mean intensity analyses
    print >> out, "Mean intensity analyses "
    print >> out, " Analyses of the mean intensity. "
    print >> out, " Inspired by: Morris et al. (2004). J. Synch. Rad.11, 56-59."
    print >> out, " The following resolution shells are worrisome: "
    worrisome = self.z_scores > z_level
    print >> out, "------------------------------------------------"
    print >> out, "| d_spacing | z_score | compl. | <Iobs>/<Iexp> |"
    print >> out, "------------------------------------------------"
    self.shell_table = data_plots.table_data(
      title="Mean intensity by shell (outliers)",
      column_labels=["d_spacing", "z_score", "completeness", "<Iobs>/<Iexp>"],
      graph_names=["Worrisome resolution shells"],
      graph_columns=[[0,1,2,3]])
    for ii in range(self.d_star_sq.size()):
      if  worrisome[ii]:
        d_space = self.d_star_sq[ii]**(-0.5)
        z_score = self.z_scores[ii]
        comp =  self.completeness[ii]
        ratio = self.mean_I_obs_data[ii]/self.mean_I_obs_theory[ii]
        print >> out, "|%9.3f  |%7.2f  |%7.2f |%10.3f     |"%(
          d_space,z_score,comp,ratio)
        self.shell_table.add_row([d_space,z_score,comp,ratio])

    if (worrisome).count(True) == 0:
      print >> out, "     None"
    print >> out, "------------------------------------------------"
    print >> out
    if (worrisome).count(True) > 0:
      self.shell_info = """\
 Possible reasons for the presence of the reported
 unexpected low or elevated mean intensity in
 a given resolution bin are :
 - missing overloaded or weak reflections
 - suboptimal data processing
 - satellite (ice) crystals
 - NCS
 - translational pseudo symmetry (detected elsewhere)
 - outliers (detected elsewhere)
 - ice rings (detected elsewhere)
 - other problems
 Note that the presence of abnormalities
 in a certain region of reciprocal space might
 confuse the data validation algorithm throughout
 a large region of reciprocal space, even though
 the data are acceptable in those areas.


"""
      print >> out, self.shell_info

    ## outlier analyses
    self.outlier.show(out)
    ## ice ring results are reported here
    if self.ijsco is not None:
      self.ijsco.show(out)

    ## say something about the anomalous data is it is available
    if self.meas_data is not None:
      self.meas_anal.show(out)

    if out_plot is not None:
      ## intensity plots
      wilson_plot = data_plots.plot_data(
        plot_title='Intensity plots',
        x_label='1/resol^2',
        y_label='<I>',
        x_data=self.d_star_sq,
        y_data=self.mean_I_normalisation,
        y_legend = '<I> smooth approximation',
        comments = 'Intensity plots'
        )
      wilson_plot.add_data(y_data=self.mean_I_obs_data,
                           y_legend='<I> via binning')
      wilson_plot.add_data(y_data=self.mean_I_obs_theory,
                           y_legend='<I> expected')
      data_plots.plot_data_loggraph(wilson_plot, out_plot)
      # XXX: Nat's newer table class (for GUI)
      self.wilson_table = data_plots.table_data(
        title="Intensity plots",
        column_labels=["1/resol**2", "<I> smooth approximation",
                      "<I> via binning", "<I> expected"],
        graph_names=["Intensity plots"],
        graph_labels=[("Resolution", "<I>")],
        graph_columns=[[0,1,2,3]],
        data=[self.d_star_sq, self.mean_I_normalisation, self.mean_I_obs_data,
              self.mean_I_obs_theory],
        x_is_inverse_d_min=True)

      ## z scores and completeness
      z_scores_and_completeness = data_plots.plot_data(
        plot_title='Z scores and completeness',
        x_label='1/resol^2',
        y_label='Z_score or fractional completeness',
        x_data=self.d_star_sq,
        y_data=self.z_scores,
        y_legend='Z score',
        comments='Data sanity and completeness check')
      z_scores_and_completeness.add_data(
        y_data=self.completeness,
        y_legend='Completeness')
      data_plots.plot_data_loggraph(z_scores_and_completeness, out_plot)
      # XXX: for GUI
      self.zscore_table = data_plots.table_data(
        title="Z scores and completeness",
        column_labels=["1/resol**2", "Z_score", "Fractional completeness"],
        graph_names=["Data sanity and completeness check"],
        graph_labels=[("Resolution", "Z score or fractional completeness")],
        graph_columns=[[0,1,2]],
        data=[self.d_star_sq, self.z_scores, self.completeness],
        x_is_inverse_d_min=True)

      ## measurability data is anomalous
      if self.meas_data is not None:
        meas_plots = data_plots.plot_data(
          plot_title='Measurability of Anomalous signal',
          x_data = self.d_star_sq_ori,
          y_data = self.meas_data,
          x_label='1/resol^2',
          y_label='Measurability',
          y_legend='Observed anomalous measurability',
          comments='Anomalous measurability')
        meas_plots.domain_flag='N'
        meas_plots.add_data(y_data=self.meas_smooth,
                            y_legend='smooth approximation')
        data_plots.plot_data_loggraph( meas_plots,   out_plot)
        # XXX: for GUI
        self.meas_table = data_plots.table_data(
          title="Measurability of Anomalous signal",
          column_labels=["1/resol**2", "Measurability", "Smooth approximation"],
          graph_names=["Anomalous measurability"],
          graph_labels=[("Resolution", "Measurability")],
          graph_columns=[[0,1,2]],
          data=[self.d_star_sq_ori, self.meas_data, self.meas_smooth],
          x_is_inverse_d_min=True)

      ## I over sigma I
      if self.i_sig_i is not None:
        i_sig_i = data_plots.plot_data(
          plot_title = '<I/sigma_I>',
          x_data = self.d_star_sq_ori,
          y_data = self.i_sig_i,
          x_label = '1/resol^2',
          y_label = '<I/sigma_I>',
          y_legend = '<I/sigma_I>; Signal to noise',
          comments = 'Signal to noise')
        i_sig_i.domain_flag='N'
        data_plots.plot_data_loggraph( i_sig_i,   out_plot)
        # XXX: for GUI
        self.i_sig_i_table = data_plots.table_data(
          title="Signal to noise (<I/sigma_I>)",
          column_labels=["1/resol**2", "<I/sigma_I>"],
          graph_names=["Signal to noise"],
          graph_labels=[("Resolution", "<I/sigma_I>")],
          graph_columns=[[0,1]],
          data=[self.d_star_sq_ori,self.i_sig_i],
          x_is_inverse_d_min=True)

      if (self.low_resolution_completeness is not None) :
        binner = self.low_resolution_completeness.binner
        d_star_sq_ori = []
        comp = []
        for i_bin in binner.range_used() :
          d_min = binner.bin_d_min(i_bin)
          d_star_sq_ori.append(1/(d_min**2))
          comp.append(self.low_resolution_completeness.data[i_bin]*100)
        self.low_res_table = data_plots.table_data(
          title="Low-resolution completeness",
          column_labels=["Max. resolution", "Completeness"],
          graph_names=["Low-resolution completeness"],
          graph_labels=[("High resolution of shell", "% of total")],
          graph_columns=[[0,1]],
          data=[d_star_sq_ori, comp],
          x_is_inverse_d_min=True,
          force_exact_x_labels=True)
