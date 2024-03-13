
"""
Main program driver for Xtriage.
"""

from __future__ import absolute_import, division, print_function
from mmtbx.scaling import data_statistics
from mmtbx.scaling import relative_wilson
from mmtbx.scaling import pair_analyses
from mmtbx.scaling import twin_analyses
from mmtbx.scaling import matthews
import mmtbx.scaling
import mmtbx.utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import iotbx.merging_statistics
import iotbx.bioinformatics
from iotbx import pdb
import iotbx.phil
from cctbx.array_family import flex
from cctbx import crystal
from libtbx.utils import Sorry, multi_out
from libtbx.str_utils import StringIO, wordwrap
from libtbx.utils import null_out
from libtbx import runtime_utils
import libtbx.callbacks # import dependency
from six.moves import cStringIO as StringIO
import operator
import os
import sys

MIN_ACENTRICS = 25
CC_ONE_HALF_OUTER_POOR = 0.3
CC_ONE_HALF_OUTER_MED = 0.7
CC_ONE_HALF_OUTER_HIGH = 0.9

master_params = iotbx.phil.parse("""\
scaling {
input {
  asu_contents
    .help = "Defines the ASU contents"
    .short_caption = ASU contents
  {
    sequence_file = None
      .type = path
      .help = File containing protein or nucleic acid sequences.  Values for \
        n_residues and n_bases will be extracted automatically if this is \
        provided.
      .style = file_type:seq input_file seq_file
    n_residues=None
      .type=float
      .help="Number of residues in structural unit"
      .short_caption = Number of residues
    n_bases=None
      .type=float
      .help="Number of nucleotides in structural unit"
      .short_caption = Number of nucleotides
    n_sites=5
      .type=int
      .help="Number of atoms in anomalous substructure"
      .short_caption = Atoms in anomalous substructure
    n_copies_per_asu=None
      .type=float
      .help="Number of copies per ASU. If not specified, Matthews analyses is performed"
      .short_caption = Number of copies in ASU
  }

  xray_data
   .help="Defines xray data"
    .style = auto_align
  {
    file_name=None
      .type=path
      .help="File name with data"
      .short_caption = Reflections
      .style = bold noauto file_type:hkl input_file process_hkl \
        update_d_max_min child:fobs:obs_labels child:fcalc:calc_labels \
        child:space_group:space_group child:unit_cell:unit_cell \
        child:d_min:high_resolution child:d_max:low_resolution
    obs_labels=None
      .type=strings
      .help="Labels for observed data"
      .short_caption = Data labels
      .input_size = 160
      .style = bold renderer:draw_xtriage_hkl_label_widget \
        OnChange:auto_update_label_choice parent:file_name:file_name \
        child:d_min:high_resolution child:d_max:low_resolution
    calc_labels=None
      .type=strings
      .help="Lables for calculated data"
      .short_caption = Fcalc labels
      .input_size = 160
      .style = bold renderer:draw_xtriage_hkl_label_widget
    unit_cell=None
      .type=unit_cell
      .help="Unit cell parameters"
      .style = bold
    space_group=None
      .type=space_group
      .help="space group"
      .style = bold
    high_resolution=None
      .type=float
      .help="High resolution limit"
      .input_size = 64
      .style = bold renderer:draw_resolution_widget
    low_resolution=None
      .type=float
      .help="Low resolution limit"
      .input_size = 64
      .style = bold renderer:draw_resolution_widget
    skip_sanity_checks = False
      .type = bool
      .help = Disable guards against corrupted input data.
      .expert_level = 3
    completeness_as_non_anomalous = True
      .type = bool
      .help = Report completeness of non-anomalous version of arrays
      .short_caption = Completeness as non-anomalous
    reference
      .help = "A reference data set. For the investigation of possible reindexing options"
      .short_caption = Reference data
      .style = noauto
    {
      data
        .help="Defines an x-ray dataset"
        .short_caption = X-ray data
        .style = box auto_align
      {
        file_name=None
          .type = path
          .help = "File name"
          .short_caption = Reference x-ray file
          .style = bold file_type:hkl process_hkl input_file update_d_max_min \
                   child:fobs:labels child:space_group:space_group \
                   child:unit_cell:unit_cell
        labels=None
          .type=strings
          .help="Labels"
          .input_size = 160
          .short_caption = Data labels
          .style = bold renderer:draw_xtriage_hkl_label_widget \
                   parent:file_name:file_name
        unit_cell=None
          .type=unit_cell
          .help=""Unit cell parameters"
        space_group=None
          .type=space_group
          .help="Space group"
      }
      structure{
         file_name=None
         .type=path
         .help="Filename of reference PDB file"
         .short_caption = Reference PDB file
         .style = file_type:pdb noauto input_file OnChange:extract_symmetry
      }
    }
  }


  parameters
   .help="Basic settings"
  {
    reporting
      .help="Some output issues"
    {
      verbose=1
        .type=int
        .help="Verbosity"
      log=logfile.log
        .type=str
        .help="Logfile"
        .style = hidden
      loggraphs = False
        .type = bool
    }

    merging
      .short_caption = Merging options
      .style = menu_item
    {
      n_bins = 10
        .type = int
        .short_caption = Number of bins for merging
      skip_merging = False
        .type = bool
        .expert_level = 1
    }

    misc_twin_parameters
      .help="Various settings for twinning or symmetry tests"
      .short_caption = Other settings
      .style = menu_item auto_align
    {
      missing_symmetry
         .help = "Settings for missing symmetry tests"
      {
        sigma_inflation = 1.25
         .type=float
         .help="Standard deviations of intensities can be increased to make point group determination more reliable."
       }

       twinning_with_ncs
       .help="Analysing the possibility of an NCS operator parallel to a twin law."
       {
         perform_analyses = False
         .type = bool
         .help = "Determines whether or not this analyses is carried out."
         .short_caption = Analyze twinning with NCS
         n_bins = 7
         .type = int
         .help = "Number of bins used in NCS analyses."
         .short_caption = Number of bins
       }

       twin_test_cuts
       .help = "Various cuts used in determining resolution limit for data used in intensity statistics "
        .short_caption = Twin test cutoffs
       .style = box auto_align
       {
         low_resolution=10.0
         .type=float
         .help="Low resolution"
         high_resolution=None
         .type=float
         .help = "High resolution"
         isigi_cut=3.0
         .type=float
         .help="I/sigI ratio used in completeness cut "
         completeness_cut=0.85
         .type=float
         .help='''Data is cut at resolution where intensities with
            I/sigI greater than isigi_cut are more than
            completeness_cut complete'''
       }
       l_test_dhkl = None
          .type = ints(size=3)
          .short_caption = Reciprocal-space vector for L-test
          .help="Offsets for hkl in applying L-test. Defaults to 2,2,2 but should be changed for commensurate modulation (tNCS) involving shifts other than 1/2 along x,y,z"
          .expert_level = 2
       apply_basic_filters_prior_to_twin_analysis=True
         .type=bool
         .help="Keep data cutoffs from the basic_analyses module (I/sigma,Wilson scaling,Anisotropy) when twin stats are computed."

     }
   }
   optional {
     include scope mmtbx.scaling.massage_twin_detwin_data.output_params_str
     include scope mmtbx.scaling.massage_twin_detwin_data.master_params
   }

   expert_level=1
    .type=int
    .expert_level=10
    .help="Expert level"
}
gui
  .help = GUI-specific parameters, not applicable to command-line version.
{
  result_file = None
    .type = path
    .help = Pickled result file for Phenix GUI
  output_dir = None
    .type = path
    .style = output_dir
    .short_caption = Output directory
  include scope libtbx.phil.interface.tracking_params
}
}
""", process_includes=True)


def make_big_header(text, out):
  if hasattr(out, 'show_big_header'):
    out.show_big_header(text)
  else:
    from libtbx.str_utils import make_big_header as mbh
    mbh(text, out)


def print_banner(appl, out=None):
  if out is None:
    out=sys.stdout
  hashes = "#############################################################"
  def print_centered(s):
    b = max(0, len(hashes) - len(s) - 4)
    l = int(b // 2)
    r = b - l
    print("##%s%s%s##" % (" "*l, s, " "*r), file=out)
  print(hashes, file=out)
  print_centered(appl)
  print_centered("")
  print_centered("P.H. Zwart, R.W. Grosse-Kunstleve & P.D. Adams")
  print_centered("")
  print(hashes, file=out)


def print_help(appl):
  print("""
--------------------------------------------------------------------------------
Usage: %(appl)s file_name=myfile.sca <options>

%(appl)s performs a variety of twinning and related test on the given
x-ray data set. See CCP4 newletter number 42, July 2005 for more information.
(http://www.ccp4.ac.uk/newsletters/newsletter42/newsletter42.pdf)


The program options are summarized below

1. scope: asu_contents
   keys: * n_residues :: Number of residues per monomer/unit
         * n_bases :: Number of nucleotides per monomer/unit
         * n_copies_per_asu :: Number of copies in the ASU.
         * n_sites :: Number of atoms in anomalous sub-structure

   These keywords control the determination of the absolute scale.
   If the number of residues/bases is not specified, a solvent content of 50%%
   is assumed.  n_sites is used in analysis of probability of
   finding the sub-structure


2a.scope: xray_data
   keys: * file_name :: file name with xray data.
         * obs_labels :: labels for observed data is format is mtz of XPLOR/CNS
         * calc_labels :: optional; labels for calculated data
         * unit_cell :: overrides unit cell in reflection file (if present)
         * space_group :: overrides space group in reflection file (if prersent)
         * high_resolution :: High resolution limit of the data
         * low_resolution :: Low resolution limit of the data

   Note that the matching of specified and present labels involves a sub-string
   matching algorithm. See 'Example usage'.


2b.scope: xray_data.reference : A reference data set or structure.
   keys:  data.file_name :: file name for xray data
          structure.file_name :: file name of a PDB.
                                 Specification of a reference structure
                                 triggers RvsR cacluations.

3. scope: parameters.misc_twin_parameters.missing_symmetry
   keys: * sigma_inflation :: Sigma inflation value in scoring function.

   sigma_intensity_in_scoring_function = sigma_intensity*sigma_inflation

   Larger values of sigma inflation will tend to result in the point group
   selection algorithm to favour higher symmetry. If data is processed
   reasonably, the default should be fine.

4. scope: parameters.misc_twin_parameters.twinning_with_ncs
   keys: * perform_test :: can be set to True or False
         * n_bins :: Number of bins in determination of D_ncs

   The perform_test is by default set to False. Setting it to True triggers the
   determination of the twin fraction while taking into account NCS parallel to
   the twin axis.


5. scope: parameters.misc_twin_parameters.twin_test_cuts
   keys: * high_resolution : high resolution for twin tests
         * low_resolution: low resolution for twin tests
         * isigi_cut: I/sig(I) threshold in automatic determination
                      of high resolution limit
         * completeness_cut: completeness threshold in automatic
                             determination of high resolution limit

   The automatic determination of the resolution limit for the twinning test
   is determined on the basis of the completeness after removing intensities for
   which I/sigI<isigi_cut. The lowest limit obtain in this way is 3.5A.
   The value determined by the automatic procedure can be overruled by
   specification of the high_resolution keyword. The low resolution is set to
   10A by default.


6. scope: parameters.reporting
   keys: * verbose :: verbosity level.
         * log :: log file name

7. scope: xray_data
   keys: * file_name :: file name with xray data.
         * obs_labels :: labels for observed data is format is mtz of XPLOR/CNS
         * calc_labels :: optional; labels for calculated data
         * unit_cell :: overrides unit cell in reflection file (if present)
         * space_group :: overrides space group in reflection file (if prersent)
         * high_resolution :: High resolution limit of the data
         * low_resolution :: Low resolution limit of the data

   Note that the matching of specified and present labels involves a sub-string
   matching algorithm. See 'Example usage'.


Example usage:

  The commands
    %(appl)s xray_data.file_name=my_refl.mtz
    %(appl)s my_refl.mtz
  are equivalent.

  The commands
    %(appl)s my_refl.mtz obs=F_HG,SIGF_HG data.high=2.0 log=log.log perform=True
    %(appl)s my_refl.mtz obs=HG data.high=2.0 log=log.log perform=True
  are equivalent if the substring 'HG' is unique in the mtz file. If the labels contain character such as
  '(', use quoation marks: obs='F_HG(+)'.


--------------------------------------------------------------------------------

""" % vars())

class merging_statistics(mmtbx.scaling.xtriage_analysis,
                          iotbx.merging_statistics.dataset_statistics):
  """
  Subclass of iotbx merging statistics class to override the show() method
  and use the Xtriage output style.
  """
  def _show_impl(self, out):
    # overall statistics
    out.show_header("Statistics for data merging")
    out.show_sub_header("Overall statistics")
    out.show_text("""\
Because unmerged intensities were supplied as input, Xtriage has calculated
statistics for merging to assess the usefulness and internal consistency of
the data (both overall and as a function of resolution).  The merged data will
be used for all subsequent steps.  Although some of the other analyses
performed later are partially redundant with these, the statistics
displayed here are the most useful for calculating the information content
and expectations for resolution.
""")
    out.show_text("""\
Note that completeness statistics are calculated with Friedel mates (I+ and I-)
treated as redundant observations; completeness for anomalous data may be
lower than shown here.
""")
    tmp_out = StringIO()
    self.overall.show_summary(out=tmp_out, prefix="    ")
    out.show_lines(tmp_out.getvalue())
    is_p1 = self.crystal_symmetry.space_group_info().type().number() == 1
    if (self.overall.r_merge > 0.2):
      if (not is_p1):
        out.warn("""\
The overall R-merge is greater than 0.2, suggesting that the choice of space
group may be incorrect.  We suggest that you try merging in lower symmetry
to see if this improves the statistics.""")
    elif (self.cc_one_half_outer > CC_ONE_HALF_OUTER_HIGH):
      out.warn("""\
The correlation of half-datasets (CC1/2)in the outer resolution shell suggests
that the useful resolution range extends significantly beyond the cutoff used
in the input data.""")
    elif (self.cc_one_half_outer > CC_ONE_HALF_OUTER_MED):
      out.warn("""\
The correlation of half-datasets (CC1/2)in the outer resolution shell suggests
that the useful resolution range may extend beyond the cutoff used in the input
data.""")
    elif (self.cc_one_half_outer < CC_ONE_HALF_OUTER_POOR):
      out.warn("""\
The correlation of half-datasets (CC1/2)in the outer resolution shell suggests
that the useful resolution range may be lower than that of the input data.""")
    # TODO more warnings?  need to figure out appropriate cutoffs...
    out.show_sub_header("Multiplicity, completeness, and signal")
    out.show_table(self.signal_table, plot_button=True)
    out.show_sub_header("Dataset consistency")
    out.show_table(self.quality_table, plot_button=True)
    out.show_lines("References:\n"+iotbx.merging_statistics.citations_str)

  @property
  def cc_one_half_outer(self):
    return self.bins[-1].cc_one_half

  def summarize_issues(self):
    if (self.overall.r_merge > 0.2):
      is_p1 = self.crystal_symmetry.space_group_info().type().number() == 1
      if (not is_p1):
        return [ (2, "The merging statistics indicate that the data may "+
          "be assigned to the wrong space group.", "Dataset consistency") ]
      else :
        return [ (2, "The merging statistics are unusually poor, but the "+
          "space group is P1, suggesting problems with data processing.",
          "Dataset consistency") ]
    if (self.cc_one_half_outer > CC_ONE_HALF_OUTER_MED):
      return [ (1, "The resolution of the data may be useful to higher "+
        "resolution than the given resolution.", "Dataset consistency") ]
    elif (self.cc_one_half_outer < CC_ONE_HALF_OUTER_POOR):
      return [ (1, "The resolution of the data may be lower than the given "+
        "resolution.", "Dataset consistency") ]
    else :
      return [ (0, "The resolution limit of the data seems appropriate.",
        "Dataset consistency") ]
    return []

class data_summary(mmtbx.scaling.xtriage_analysis):
  """
  Basic info about the input data (somewhat redundant at the moment).
  """
  def __init__(self, miller_array, was_merged=False):
    self.info = miller_array.info()
    self.space_group = miller_array.space_group_info()
    self.unit_cell = miller_array.unit_cell().parameters()
    self.obs_type = miller_array.observation_type()
    self.was_merged = was_merged
    self.n_indices = self.n_indices_merged = miller_array.size()
    self.anomalous_flag = miller_array.anomalous_flag()
    self.completeness = self.completeness_merged = miller_array.completeness(
      as_non_anomalous_array=False)
    self.d_max_min = miller_array.d_max_min()
    self.anomalous_completeness = None
    non_anom = miller_array.average_bijvoet_mates()
    if (self.anomalous_flag):
      self.n_indices_merged = non_anom.size()
      self.completeness_merged = non_anom.completeness(
        as_non_anomalous_array=False)
      self.anomalous_completeness = miller_array.anomalous_completeness()

  def _show_impl(self, out):
    out.show_header("Input data")
    out.show_sub_header("Summary")
    if (self.was_merged):
      out.show_text("""\
The original dataset contained unmerged intensities; the statistics below
are for the merged data.""")
    source = getattr(self.info, "source", None)
    label_string = getattr(self.info, "label_string", lambda: None)
    if (source is not None):
      if (os.path.isfile(source)):
        source = os.path.basename(source)
    lines = [
      ("File name:", str(source)),
      ("Data labels:", str(label_string())),
      ("Space group:", str(self.space_group)),
      ("Unit cell:", "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g" % self.unit_cell),
      ("Data type:", str(self.obs_type)),
      ("Resolution:", "%g - %g" % self.d_max_min),
      ("Anomalous:", str(self.anomalous_flag)),
    ]
    if (self.anomalous_flag):
      lines.extend([
        ("Number of reflections (non-anomalous):", str(self.n_indices_merged)),
        ("Completeness (non-anomalous):", "%.2f%%" %
          (self.completeness_merged*100)),
        ("Number of reflections (all):", str(self.n_indices)),
        ("Completeness (all):", "%.2f%%" % (self.completeness*100)),
        ("Anomalous completeness:", "%.2f%%" %
          (self.anomalous_completeness*100)),
      ])
    else :
      lines.extend([
        ("Number of reflections:", str(self.n_indices)),
        ("Completeness:", "%.2f%%" % (self.completeness*100)),
      ])
    out.show_text_columns(lines, indent=2)
    # description of completeness values
    if (self.anomalous_flag):
      out.show_text("""
  Completeness (non-anomalous) is the completeness of the data after merging
  Friedel pairs.

  Completeness (all) is the completeness of the data before merging Friedel
  pairs.

  Completeness (anomalous) is the completeness of the anomalous data. The
  anomalous completeness is calcluated by dividing the number of measured
  acentric, Bijvoet mates by the total possible number of acentric indices.

  Completeness (non-anomalous) should be used to determine if there is
  sufficient data for non-anomalous purposes (refinement, model-building).
  A value greater than 90% is generally desired, while a value less than
  75% is considered poor. Values in between will provide less than optimal
  results.

  Completeness (anomalous) should be used to determine if there is sufficient
  data for anomalous purposes (phasing, finding anomalous atoms, refinement
  of anomalous occupancies or scattering factors). A value greater than 80%
  is generally desired, while a value less than 50% is considered poor. Values
  in between will provide less than optimal results.
  """)
    else:
      out.show_text("""
  Completeness should be used to determine if there is sufficient data for
  refinement and/or model-building. A value greater than 90% is generally
  desired, while a value less than 75% is considered poor. Values in between
  will provide less than optimal results.
  """)

  def summarize_issues(self):
    issues = []
    if (self.anomalous_flag):
      message = 'The non-anomalous completeness is %.2f%%.'%(100.0*self.completeness_merged)
      if (self.completeness_merged < 0.75):
        issues.append((2, message, 'Summary'))
      elif (self.completeness_merged < 0.90):
        issues.append((1, message, 'Summary'))
      else :
        issues.append((0, message, 'Summary'))
      message = 'The anomalous completeness is %.2f%%.'%(100.0*self.anomalous_completeness)
      if (self.anomalous_completeness < 0.5):
        issues.append((1.5, message, 'Summary'))
      elif (0.5 <= self.anomalous_completeness < 0.8):
        issues.append((1, message, 'Summary'))
      else:
        issues.append((0,message,'Summary'))
    else:
      message = 'The completeness is %.2f%%.'%(100.0*self.completeness)
      if (self.completeness < 0.75):
        issues.append((2, message, 'Summary'))
      elif (self.completeness < 0.90):
        issues.append((1, message, 'Summary'))
      else :
        issues.append((0, message, 'Summary'))
    return issues

class xtriage_analyses(mmtbx.scaling.xtriage_analysis):
  """
  Run all Xtriage analyses for experimental data, with optional Fcalc or
  reference datasets.

  :param miller_obs: array of observed data, should be intensity or amplitude
  :param miller_calc: array of calculated data
  :param miller_ref: array with 'reference' data, for instance a data set with
                     an alternative indexing scheme
  :param text_out: A filehandle or other object with a write method
  :param params: An extracted PHIL parameter block, derived from master_params
  """
  new_format = True # XXX used by Phenix GUI
  def __init__(self,
      # array of observed data, should be intensity or amplitude
      miller_obs,
      miller_calc  = None,
      miller_ref   = None,
      params       = None,
      text_out     = None,
      unmerged_obs = None,
      log_file_name = None):
    self.log_file_name = log_file_name
    assert (miller_obs is not None)


    if miller_obs.is_unmerged_intensity_array():
      unmerged_obs = miller_obs.deep_copy()
    def process_input_array(miller_array):
      if (miller_array is None) : return None
      info = miller_array.info()
      miller_array = miller_array.merge_equivalents().array()
      return miller_array.remove_systematic_absences().set_info(info)
    miller_obs = process_input_array(miller_obs)
    miller_original = miller_obs.deep_copy()
    miller_calc = process_input_array(miller_calc)
    miller_ref = process_input_array(miller_ref)
    self.data_summary = data_summary(miller_obs,
      was_merged=(unmerged_obs is not None))
    if text_out == "silent": # suppress all output
      text_out = null_out()
    if params is None:         # if nothing specified, defaults are used
      params = master_params.fetch(sources=[]).extract()

    self.completeness_as_non_anomalous=\
       params.scaling.input.xray_data.completeness_as_non_anomalous
    if (not params.scaling.input.xray_data.skip_sanity_checks):
      check_for_pathological_input_data(miller_obs,
        completeness_as_non_anomalous=
          self.completeness_as_non_anomalous)

    ###
    self.merging_stats = None
    if ((unmerged_obs is not None) and
        (unmerged_obs.is_xray_intensity_array())):
      try :
        self.merging_stats = merging_statistics(
          i_obs=unmerged_obs,
          crystal_symmetry=miller_obs.crystal_symmetry(),
          d_min=params.scaling.input.xray_data.high_resolution,
          d_max=params.scaling.input.xray_data.low_resolution,
          n_bins=params.scaling.input.parameters.merging.n_bins,
          log=text_out)
        self.merging_stats.show(out=text_out)
        print("", file=text_out)
        print("References:", file=text_out)
        print(iotbx.merging_statistics.citations_str, file=text_out)
        print("", file=text_out)
      except Exception as e :
        print("WARNING: calculation of merging statistics failed", file=text_out)
        print("  error: %s" % str(e), file=text_out)

    self.plan_sad_experiment_stats= None
    if ((unmerged_obs is not None) and
        (unmerged_obs.is_xray_intensity_array())):
      try :
        from phenix.command_line.anomalous_signal import anomalous_signal
        self.plan_sad_experiment_stats= anomalous_signal(
          i_obs=unmerged_obs,
          d_min=params.scaling.input.xray_data.high_resolution,
          d_max=params.scaling.input.xray_data.low_resolution,
          sites=params.scaling.input.asu_contents.n_sites,
          skip_scaling=True,
          out=null_out()).analysis
      except Exception as e :
        print("WARNING: calculation of anomalous statistics failed", file=text_out)
        print("  error: %s\n" % str(e), file=text_out)
    ###
    make_big_header("Basic statistics", out=text_out)
    n_copies_solc = 1.0
    nres_known = False
    if (params.scaling.input.asu_contents.n_residues is not None or
        params.scaling.input.asu_contents.n_bases is not None):
      nres_known = True
      if (params.scaling.input.asu_contents.sequence_file is not None):
        print("  warning: ignoring sequence file", file=text_out)
    elif (params.scaling.input.asu_contents.sequence_file is not None):
      print("  determining composition from sequence file %s" % \
        params.scaling.input.asu_contents.sequence_file, file=text_out)
      seq_comp = iotbx.bioinformatics.composition_from_sequence_file(
        file_name=params.scaling.input.asu_contents.sequence_file,
        log=text_out)
      if (seq_comp is not None):
        params.scaling.input.asu_contents.n_residues = seq_comp.n_residues
        params.scaling.input.asu_contents.n_bases = seq_comp.n_bases
        nres_known = True
    matthews_results = matthews.matthews_rupp(
      crystal_symmetry = miller_obs,
      n_residues = params.scaling.input.asu_contents.n_residues,
      n_bases=params.scaling.input.asu_contents.n_bases,
      out=text_out)
    self.matthews = matthews_results
    self.matthews.show(out=text_out)
    params.scaling.input.asu_contents.n_residues = matthews_results.n_residues
    params.scaling.input.asu_contents.n_bases = matthews_results.n_bases
    n_copies_solc = matthews_results.n_copies
    if params.scaling.input.asu_contents.n_copies_per_asu is not None:
      n_copies_solc = params.scaling.input.asu_contents.n_copies_per_asu
      print("Number of copies per asymmetric unit provided", file=text_out)
      print(" Will use user specified value of ", n_copies_solc, file=text_out)
    else:
      params.scaling.input.asu_contents.n_copies_per_asu = n_copies_solc

    # Signal-to-noise and completeness
    self.data_strength_and_completeness = \
      data_statistics.data_strength_and_completeness(
        miller_array=miller_obs,completeness_as_non_anomalous=\
           self.completeness_as_non_anomalous
        ).show(out=text_out)
    # completeness summary
    self.data_summary.show(out=text_out)
    # Anomalous signal
    self.anomalous_info = None
    if miller_obs.anomalous_flag():
      self.anomalous_info = data_statistics.anomalous(
        miller_array=miller_obs,
        merging_stats=self.merging_stats,
        plan_sad_experiment_stats=self.plan_sad_experiment_stats).show(
          out=text_out)
    # Wilson statistics
    self.wilson_scaling = data_statistics.wilson_scaling(
      miller_array=miller_obs,
      n_residues=params.scaling.input.asu_contents.n_residues,
      n_bases=params.scaling.input.asu_contents.n_bases,
      n_copies_solc=params.scaling.input.asu_contents.n_copies_per_asu,
      completeness_as_non_anomalous=\
           self.completeness_as_non_anomalous)
    self.wilson_scaling.show(out=text_out)
    self.relative_wilson = None
    # XXX The resolution filter isn't perfect - sometimes the subsequent steps
    # in relative_wilson.py result in an effective resolution worse than 4.0.
    # An example in the PDB is 2bf1:FOBS,SIGFOBS (d_min=3.98).
    # This will now raise Sorry instead and the program will be allowed to
    # continue, but should probably be fixed.  Should we simply be using a more
    # stringent resolution cutoff?
    if (miller_calc is not None) and (miller_calc.d_min() < 4.0):
      try :
        self.relative_wilson = relative_wilson.relative_wilson(
          miller_obs=miller_obs,
          miller_calc=miller_calc).summary()
      except Sorry as e :
        print("*** Error calculating relative Wilson plot - skipping.", file=text_out)
        print(str(e), file=text_out)
        print("", file=text_out)
      except Exception as e :
        print("RELATIVE WILSON ERROR")
        raise
        print("", file=text_out)
      else : # FIXME
        pass # self.relative_wilson.show(out=text_out)
    text_out.flush()
    # outliers are removed, make a new copy
    twin_params = params.scaling.input.parameters.misc_twin_parameters
    if (twin_params.twin_test_cuts.high_resolution is None):
      twin_params.twin_test_cuts.high_resolution = \
        self.data_strength_and_completeness.high_resolution_for_twin_tests()
    if twin_params.apply_basic_filters_prior_to_twin_analysis:
      new_miller = self.wilson_scaling.miller_array_filtered
      miller_obs = new_miller.deep_copy()
      normalised_array = self.wilson_scaling.normalised_miller.deep_copy()

    print("", file=text_out)
    #Do the twinning analyses
    ## resolution check
    if (twin_params.twin_test_cuts.high_resolution is not None
        and flex.min(miller_obs.d_spacings().data())
        > twin_params.twin_test_cuts.high_resolution):
      params.scaling.input.xray_data.high_resolution = flex.min(
        miller_obs.d_spacings().data())

    default_high_reso_limit_wilson_ratio = \
      twin_params.twin_test_cuts.high_resolution
    if default_high_reso_limit_wilson_ratio is None:
      default_high_reso_limit_wilson_ratio = 0.0

    d_star_sq_high_limit = default_high_reso_limit_wilson_ratio
    d_star_sq_high_limit = 1.0/((d_star_sq_high_limit+1e-6)**2.0)
    default_low_reso_limit_wilson_ratio = \
      twin_params.twin_test_cuts.low_resolution

    d_star_sq_low_limit = default_low_reso_limit_wilson_ratio
    d_star_sq_low_limit = 1.0/((d_star_sq_low_limit+1e-6)**2.0)
    self.twin_results = None
    acentrics = miller_obs.select_acentric()
    n_acentrics = acentrics.size()
    if (n_acentrics > 0):
      make_big_header("Twinning and symmetry analyses", out=text_out)
      self.twin_results = twin_analyses.twin_analyses(
        miller_array=miller_obs,
        d_star_sq_low_limit=d_star_sq_low_limit,
        d_star_sq_high_limit=d_star_sq_high_limit,
        d_hkl_for_l_test=twin_params.l_test_dhkl,
        normalise=True,
        out=text_out,
        miller_calc=miller_calc,
        additional_parameters=twin_params,
        original_data=miller_original,
        completeness_as_non_anomalous=\
           self.completeness_as_non_anomalous)
      self.twin_results.show(text_out)
      text_out.flush()
    else :
      assert miller_obs.space_group().is_centric()
      print("", file=text_out)
      print("Centric space group - skipping twin analyses.", file=text_out)
      print("", file=text_out)

    if miller_ref is not None:
      self.reference_analyses = pair_analyses.reindexing(
        miller_ref,
        miller_obs,
        file_name=params.scaling.input.xray_data.file_name)

  def _show_impl(self, out):
    self.data_summary.show(out)
    if (self.merging_stats is not None):
      self.merging_stats.show(out)
    if isinstance(out, mmtbx.scaling.printed_output):
      make_big_header("Basic analyses", out=out)
    self.matthews.show(out)
    self.data_strength_and_completeness.show(out)
    self.wilson_scaling.show(out)
    if (self.relative_wilson is not None):
      self.relative_wilson.show(out)
    if (self.anomalous_info is not None):
      self.anomalous_info.show(out)
    if (self.twin_results is not None):
      if isinstance(out, mmtbx.scaling.printed_output):
        make_big_header("Twinning and symmetry", out=out)
      self.twin_results.show(out)
    summary = self.summarize_issues()
    summary.show(out)

  def matthews_n_copies(self):
    """
    Convenience method for retrieving the number of copies.
    """
    return self.matthews.n_copies

  def resolution_cut(self):
    """
    Convenience method for retrieving a conservative resolution cutoff.
    """
    ds = self.data_strength_and_completeness
    return getattr(ds.data_strength, "completeness_cut", None)

  def is_twinned(self):
    """
    Convenience method for indicating whether the data are likely twinned.
    """
    if (self.twin_results is not None):
      return self.twin_results.twin_summary.has_twinning()
    return False

  def resolution_limit_of_anomalous_signal(self):
    """
    Convenience method for retrieving the recommended resolution cutoff for
    anomalous substructures search.  Used in AutoSol.
    """
    return getattr(self.anomalous_info, "low_d_cut", None)

  @property
  def overall_i_sig_i(self):
    return self.data_strength_and_completeness.overall_i_sig_i

  @property
  def i_over_sigma_outer_shell(self):
    return self.data_strength_and_completeness.i_over_sigma_outer_shell()

  @property
  def iso_b_wilson(self):
    """Convenience method for isotropic Wilson B-factor"""
    return self.wilson_scaling.iso_b_wilson

  @property
  def low_d_cut(self):
    """Shortcut to resolution_limit_of_anomalous_signal()."""
    return self.resolution_limit_of_anomalous_signal()

  @property
  def aniso_b_min(self):
    """
    Convenience method for retrieving the minimum anisotropic B_cart tensor.
    Used in AutoSol.
    """
    b_cart = self.wilson_scaling.aniso_scale_and_b.b_cart
    return min(b_cart[0:3])

  @property
  def aniso_range_of_b(self):
    """
    Convenience method for retrieving the range of anisotropic B_cart tensors.
    Used in AutoSol.
    """
    b_cart = self.wilson_scaling.aniso_scale_and_b.b_cart
    return max(b_cart[0:3]) - min(b_cart[0:3])

  @property
  def aniso_b_ratio(self):
    """
    Ratio of the maximum difference between anisotropic B_cart tensors to the
    mean of the tensors.  Used in PDB validation server.
    """
    b_cart = self.wilson_scaling.aniso_scale_and_b.b_cart
    return (max(b_cart[0:3]) - min(b_cart[0:3])) / (sum(b_cart[0:3]) / 3)

  @property
  def number_of_wilson_outliers(self):
    """
    Number of centric and acentric outliers flagged by Wilson plot analysis.
    Used in PDB validation server.
    """
    return self.wilson_scaling.outliers.n_outliers()

  @property
  def l_test_mean_l(self):
    """
    <|L|> from the L test for abnormal intensity distributions.  Used in PDB
    validation server.
    """
    return self.twin_results.l_test.mean_l

  @property
  def l_test_mean_l_squared(self):
    """
    <L^2> from the L test for abnormal intensity distributions.  Used in PDB
    validation server.
    """
    return self.twin_results.l_test.mean_l2

  # FIXME this uses the Britton test, but the PDB validation server appears to
  # use the H test.  Which is correct?
  @property
  def max_estimated_twin_fraction(self):
    """
    Estimated twin fraction from the most worrysome twin law.  Used by PDB
    validation server.
    """
    return self.twin_results.twin_summary.max_twin_fraction()

  @property
  def patterson_verdict(self):
    """
    Plain-English explanation of Patterson analysis for TNCS detection.  Used
    by PDB validation server.
    """
    return self.twin_results.twin_summary.patterson_verdict()

  def estimate_d_min(self, **kwds):
    """
    Suggest resolution cutoffs based on selected statistics (if merging was
    included in analyses).  See
    :py:class:`iotbx.merging_statistics.dataset_statistics` for underlying
    function documentation.
    """
    if (self.merging_stats is not None):
      return self.merging_stats.estimate_d_min(**kwds)
    return None

  def summarize_issues(self):
    issues = []
    if (self.twin_results is not None):
      issues.extend(self.twin_results.twin_summary.summarize_issues())
    issues.extend(self.wilson_scaling.summarize_issues())
    if (self.merging_stats is not None):
      issues.extend(self.merging_stats.summarize_issues())
    issues.extend(self.data_strength_and_completeness.summarize_issues())
    issues.extend(self.data_summary.summarize_issues())
    if (self.anomalous_info is not None):
      if (hasattr(self.anomalous_info, 'summarize_issues')):
        issues.extend(self.anomalous_info.summarize_issues())
    return summary(issues)

class summary(mmtbx.scaling.xtriage_analysis):
  def __init__(self, issues, sort=True):
    self._issues = issues
    if (sort):
      self._issues.sort(key=operator.itemgetter(0), reverse=True)

  @property
  def n_problems(self):
    return len(self._issues)

  def _show_impl(self, out):
    out.show_header("Xtriage summary")
    some_problems = False
    have_problems = False
    if (self.n_problems > 0):
      if hasattr(out, "show_issues") : # XXX GUI hack
        for severity, message, linkto in self._issues :
          if (severity == 2):
            have_problems = True
            break
          elif ( (severity < 2) and (severity >= 1) ):
            some_problems = True
            break
        out.show_issues(self._issues)
      else :
        for severity, message, linkto in self._issues :
          if (severity == 2):
            have_problems = True
            out.warn(message)
          elif ( (severity < 2) and (severity >= 1) ):
            some_problems = True
            out.warn(message)
          else :
            out.show(wordwrap(message, max_chars=78))
    if (have_problems):
      out.show("""
Please inspect all individual results closely, as it is difficult to
automatically detect all issues.""")
    elif (some_problems):
      out.show("""
There are some aspects of the data that are of concern, please check the
individual results.""")
    else:
      out.show("""
No obvious problems were found with this dataset.  However, we recommend that
you inspect the individual results closely, as it is difficult to automatically
detect all issues.""")

def check_for_pathological_input_data(miller_array,
   completeness_as_non_anomalous=None):
  acentrics = miller_array.select_acentric()
  n_acentrics = acentrics.size()
  if (n_acentrics == 0) and (not miller_array.space_group().is_centric()):
    raise Sorry(("No acentric reflections present in the data.  Since the "+
      "space group %s is not centric, this is probably an error in "+
      "the input file.  Please check your data processing and/or file "+
      "conversion.") % (str(miller_array.space_group_info())))
  elif (n_acentrics < MIN_ACENTRICS):
    raise Sorry(("Only %d acentric reflections are present in these data; "+
      "this usually indicates a corrupted input file.") % n_acentrics)
  cmpl = miller_array.completeness(
    as_non_anomalous_array=completeness_as_non_anomalous)
  if (cmpl < 0.1):
    raise Sorry(("Data are unusually incomplete (%.1f%%); this usually "+
      "indicates a corrupted input file or data processing errors.") %
      (cmpl * 100))
  min_mi, max_mi = miller_array.min_max_indices()
  if ((min_mi[0] == max_mi[0]) or (min_mi[1] == max_mi[1]) or
      (min_mi[2] == max_mi[2])):
    raise Sorry("These data appear to be limited to a single plane of "+
      "reciprocal space (max/min Miller indices: %s, %s)." % (min_mi, max_mi))

def run(args, command_name="phenix.xtriage", return_result=False,
    out=None, data_file_name=None):
  if (out is None):
    out = sys.stdout
  log = multi_out()
  log.register(label="stdout", file_object=out)
  string_buffer = StringIO()
  string_buffer_plots = StringIO()
  log.register(label="log_buffer", file_object=string_buffer)

  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_params,
    pdb_file_def="scaling.input.xray_data.reference.structure.file_name",
    seq_file_def="scaling.input.asu_contents.sequence_file",
    reflection_file_def="scaling.input.xray_data.file_name",
    unit_cell_def="scaling.input.xray_data.unit_cell",
    space_group_def="scaling.input.xray_data.space_group",
    usage_string="phenix.xtriage [options] reflection_file parameters [...]")
  effective_params = cmdline.work
  params = effective_params.extract()
  if (params.scaling.input.xray_data.file_name is None):
    raise Sorry("No reflection file in input.")
  reference_params = params.scaling.input.xray_data.reference
  verbose = params.scaling.input.parameters.reporting.verbose
  scope =  params.scaling.input.xray_data

  ## Check for number of residues

  if params.scaling.input.asu_contents.n_residues is None:
    print("##-------------------------------------------##", file=log)
    print("## WARNING:                                  ##", file=log)
    print("## Number of residues unspecified            ##", file=log)
    print("##-------------------------------------------##", file=log)

  if params.scaling.input.xray_data.file_name is None:
    raise Sorry("No reflection file defined")

  ## Check for unit cell and spacegroup definitions
  crystal_symmetry = crystal_symmetry_from_any.extract_from(
      file_name=params.scaling.input.xray_data.file_name)
  if (crystal_symmetry is not None):
    space_group = crystal_symmetry.space_group_info()
    unit_cell = crystal_symmetry.unit_cell()
  else :
    space_group = params.scaling.input.xray_data.space_group
    unit_cell = params.scaling.input.xray_data.unit_cell
  if (None in [crystal_symmetry, space_group, unit_cell]):
    print("Cell and/or symmetry not specified in reflection file", file=log)
    if (reference_params.structure.file_name is not None):
      print("Using reference PDB file", file=out)
      crystal_symmetry = crystal_symmetry_from_any.extract_from(
        file_name=reference_params.structure.file_name)
    if (crystal_symmetry is None and
        (params.scaling.input.xray_data.reference.data.file_name is not None)):
      crystal_symmetry = crystal_symmetry_from_any.extract_from(
        file_name=params.scaling.input.xray_data.reference.data.file_name)
  if crystal_symmetry is None:
    if params.scaling.input.xray_data.space_group is None:
      raise Sorry("""No space group info available.
Use keyword 'xray_data.space_group' to specify space group""" )

    if (params.scaling.input.xray_data.unit_cell is None):
      raise Sorry("""
No unit cell info available.
Use keyword 'xray_data.unit_cell' to specify unit_cell
                  """ )
  #provisions for nomerge original index
  else:
    if crystal_symmetry.unit_cell() is None:
      if params.scaling.input.xray_data.unit_cell is None:
        raise Sorry("""No unit cell info available.
  Use keyword 'xray_data.unit_cell' to specify unit_cell""" )
      else: pass
    if (crystal_symmetry.space_group() is None):
      if params.scaling.input.xray_data.space_group is None:
        raise Sorry("""No space group info available.
  Use keyword 'xray_data.space_group' to specify space_group""" )

  if (params.scaling.input.xray_data.unit_cell is None):
    params.scaling.input.xray_data.unit_cell = crystal_symmetry.unit_cell()
  if (params.scaling.input.xray_data.space_group is None):
    params.scaling.input.xray_data.space_group = \
      crystal_symmetry.space_group_info()

  new_params = master_params.format(python_object=params)
  if (params.scaling.input.parameters.reporting.verbose>0):
    print(file=log)
    print(file=log)
    print("Effective parameters: ", file=log)
    print("#phil __ON__", file=log)
    new_params.show(out=log,expert_level=params.scaling.input.expert_level)
    print("#phil __END__", file=log)
  crystal_symmetry = crystal.symmetry(
    unit_cell = params.scaling.input.xray_data.unit_cell,
    space_group_symbol = str(params.scaling.input.xray_data.space_group) )

  ## Please check if we have a acentric space group
  if crystal_symmetry.space_group().is_centric():
    libtbx.warn(("The specificed space group (%s) is centric; Xtriage will "+
      "still run, but many analyses will be skipped.") %
      str(crystal_symmetry.space_group_info()))

  ## Now it time to read in reflection files somehow
  ## We do this via a reflection_file_server

  xray_data_server =  reflection_file_utils.reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    force_symmetry = True,
    reflection_files=[])

  miller_array = unmerged_array = None
  miller_array = xray_data_server.get_xray_data(
    file_name = params.scaling.input.xray_data.file_name,
    labels = params.scaling.input.xray_data.obs_labels,
    ignore_all_zeros = True,
    parameter_scope = 'scaling.input.xray_data',
    parameter_name = 'obs_labels'
  )

  # If the input data were unmerged intensities, read them in again without
  # merging - this array will be used to calcluate merging statistics
  info = miller_array.info()
  if (info.merged) and (miller_array.is_xray_intensity_array()):
    hkl_in_raw = reflection_file_reader.any_reflection_file(
      file_name=params.scaling.input.xray_data.file_name)
    assert (hkl_in_raw.file_type() is not None)
    raw_arrays = hkl_in_raw.as_miller_arrays(
      crystal_symmetry=miller_array,
      merge_equivalents=False)
    for array in raw_arrays :
      if (array.info().labels == info.labels):
        if (array.is_unmerged_intensity_array()):
          unmerged_array = array
          print("", file=log)
          print("Also reading unmerged data as %s" % \
            array.info().label_string(), file=log)
          print("", file=log)
          break

  if not miller_array.is_real_array():
    miller_array = abs( miller_array )
    from cctbx.xray import observation_types
    miller_array = miller_array.set_observation_type(
      observation_types.amplitude() )
  miller_array = miller_array.map_to_asu()

  if miller_array.observation_type() is None:
    raise Sorry("Observation type of data unkown. Please check input reflection file")

  miller_array = miller_array.select(
    miller_array.indices() != (0,0,0))
  if (miller_array.is_complex_array()):
    miller_array = abs(miller_array)

  # first do a low reso cutn if applicable
  if params.scaling.input.xray_data.low_resolution is not None:
    miller_array = miller_array.resolution_filter(
      d_max=params.scaling.input.xray_data.low_resolution)
  if params.scaling.input.xray_data.high_resolution is not None:
    miller_array = miller_array.resolution_filter(
      d_min=params.scaling.input.xray_data.high_resolution)

  # make sure sigmas are okai, otherwise, cut them
  if (not miller_array.sigmas_are_sensible()):
    #clearly there is something wrong with the sigmas
    #forget about them I would say
    miller_array = miller_array.customized_copy(
      indices=miller_array.indices(),
      data=miller_array.data(),
      sigmas=None ).set_observation_type( miller_array )

  miller_array = miller_array.eliminate_sys_absent(
    integral_only=True, log=log)

  ## Check if Fcalc label is available
  f_calc_miller = None
  f_calc_miller_complex = None
  reference_structure = None
  reference_params = params.scaling.input.xray_data.reference
  if (params.scaling.input.xray_data.calc_labels is not None):
    f_calc_miller = xray_data_server.get_amplitudes(
      file_name = params.scaling.input.xray_data.file_name,
      labels = params.scaling.input.xray_data.calc_labels,
      convert_to_amplitudes_if_necessary=False,
      parameter_scope = 'scaling.input.xray_data',
      parameter_name = 'calc_labels')
    if not f_calc_miller.is_real_array():
      f_calc_miller = f_calc_miller.customized_copy(
        data = flex.abs( f_calc_miller.data() ) ).set_observation_type(
          f_calc_miller)
    if (f_calc_miller.is_xray_intensity_array()):
      print("Converting %s to amplitudes" % \
        (params.scaling.input.xray_data.calc_labels), file=log)
      f_calc_miller = f_calc_miller.f_sq_as_f()
    f_calc_miller = f_calc_miller.eliminate_sys_absent(integral_only=True,
      log=log)
    if (miller_array.anomalous_flag()):
      if (not f_calc_miller.anomalous_flag()):
        f_calc_miller = f_calc_miller.generate_bijvoet_mates()
    elif (f_calc_miller.anomalous_flag()):
      f_calc_miller = f_calc_miller.average_bijvoet_mates()
    f_calc_miller.set_observation_type_xray_amplitude()
  elif (reference_params.structure.file_name is not None):
    if (not os.path.isfile(reference_params.structure.file_name)):
      raise Sorry("Can't open reference structure - not a valid file.")
    assert f_calc_miller is None
    pdb_in = iotbx.pdb.input(
      file_name=reference_params.structure.file_name,
      raise_sorry_if_format_error=True)
    pdb_in.atoms().set_chemical_element_simple_if_necessary()
    reference_structure = pdb_in.xray_structure_simple(
        crystal_symmetry = miller_array.crystal_symmetry())
    miller_tmp = miller_array
    if (not miller_tmp.is_unique_set_under_symmetry()):
      miller_tmp = miller_tmp.merge_equivalents().array()
    tmp_obj = mmtbx.utils.fmodel_from_xray_structure(
      xray_structure = reference_structure,
      f_obs = miller_tmp)
    f_calc_miller_complex = tmp_obj.f_model
    f_calc_miller = abs( tmp_obj.f_model ).eliminate_sys_absent(
      integral_only=True,
      log=log).set_observation_type_xray_amplitude()
  twin_results = None

  #-------------------------------------------------------------------
  # REFERENCE DATA
  reference_array = None
  if params.scaling.input.xray_data.reference.data.file_name is not None:
    user_cell=None
    user_group=None
    if params.scaling.input.xray_data.reference.data.unit_cell is not None:
      user_cell = params.scaling.input.xray_data.reference.data.unit_cell
    if params.scaling.input.xray_data.reference.data.space_group is not None:
      user_group = params.scaling.input.xray_data.reference.data.space_group
    reference_symmetry = None
    if (user_cell is not None) and (user_group is not None):
      reference_symmetry= crystal.symmetry(
        unit_cell=user_cell,
        space_group=user_group.group())
    if reference_symmetry is None:
      reference_symmetry = crystal_symmetry_from_any.extract_from(
        file_name=params.scaling.input.xray_data.reference.data.file_name)
    if reference_symmetry is None:
      print("No reference unit cell and space group could be deduced", file=log)
      raise Sorry(
        "Please provide unit cell and space group for reference data")
    xray_data_server =  reflection_file_utils.reflection_file_server(
      crystal_symmetry = reference_symmetry,
      force_symmetry = True,
      reflection_files=[])
    reference_array =  None
    reference_array = xray_data_server.get_xray_data(
      file_name = params.scaling.input.xray_data.reference.data.file_name,
      labels = params.scaling.input.xray_data.reference.data.labels,
      ignore_all_zeros = True,
      parameter_scope = 'scaling.input.xray_data.reference.data',
      parameter_name = 'labels')
    info = reference_array.info()
    reference_array = reference_array.map_to_asu()
    reference_array = reference_array.select(
      reference_array.indices() != (0,0,0))
    reference_array = reference_array.select(reference_array.data() > 0)
    if (reference_array.is_xray_intensity_array()):
      reference_array = reference_array.f_sq_as_f()
    elif (reference_array.is_complex_array()):
      reference_array = abs(reference_array)
    reference_array = reference_array.eliminate_sys_absent(
      integral_only=True, log=log)

  # make sure we hold on to the 'raw' data for later usage is desired
  raw_data = miller_array.deep_copy()
  xtriage_results = None
  assert miller_array.is_real_array()
  print(file=log)
  print(file=log)
  print("Symmetry, cell and reflection file content summary", file=log)
  print(file=log)
  miller_array.set_info(info=info)
  miller_array.show_comprehensive_summary(f=log)

  minimal_pass = True
  twin_pass = True
  reference_pass = False

  if params.scaling.input.xray_data.reference.data.file_name is not None:
    reference_pass = True

  xtriage_results = xtriage_analyses(
    miller_obs    = miller_array,
    miller_calc   = f_calc_miller,
    miller_ref    = reference_array,
    params        = params,
    text_out      = log,
    unmerged_obs  = unmerged_array,
    log_file_name = params.scaling.input.parameters.reporting.log)

  if(params.scaling.input.parameters.reporting.log is not None):
    with open( params.scaling.input.parameters.reporting.log  ,'w') as output_file:
      output_file.write(string_buffer.getvalue())

  if (params.scaling.input.optional.hklout is not None):
    # FIXME DEPRECATED, replace with mmtbx.command_line.massage_data
    from mmtbx.scaling import massage_twin_detwin_data
    massaged_obs = massage_twin_detwin_data.massage_data(
      miller_array=raw_data,
      parameters=params.scaling.input.optional,
      out=log)
    params2 = params.scaling.input.optional
    massaged_obs.write_data(
      file_name=params2.hklout,
      output_type=params2.hklout_type,
      label_extension=params2.label_extension)


  if (data_file_name is not None):
    from libtbx import easy_pickle
    easy_pickle.dump(data_file_name, raw_data)
  if (params.scaling.input.parameters.reporting.loggraphs):
    with open(params.scaling.input.parameters.reporting.log, 'a') as output_file:
      graph_out = mmtbx.scaling.loggraph_output(output_file)
      xtriage_results.show(out=graph_out)
  return xtriage_results

def change_symmetry(miller_array, space_group_symbol, file_name=None,
    log=sys.stdout):
  """
  Encapsulates all operations required to convert the original data to a
  different symmetry as suggested by Xtriage.
  """
  miller_array = miller_array.change_symmetry(
    space_group_symbol=space_group_symbol,
    log=log)
  if (file_name is not None):
    column_root_label = None
    if (miller_array.is_xray_amplitude_array()):
      column_root_label = "F"
    elif (miller_array.is_xray_intensity_array()):
      column_root_label = "I"
    if (column_root_label is None):
      raise RuntimeError("Only amplitudes and intensites supported.")
    miller_array.as_mtz_dataset(
      column_root_label=column_root_label).mtz_object().write(file_name)
  return miller_array

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=list(self.args),
      data_file_name="xtriage_data.pkl")

def finish_job(result):
  output_files = []
  stats = []
  if (result is not None):
    if (result.log_file_name is not None):
      output_files.append((result.log_file_name, "Log file and graphs"))
    if (result.low_d_cut is not None):
      stats.append(("Limit of useful anomalous signal",
        "%.2f" % result.low_d_cut))
    if (result.iso_b_wilson is not None):
      stats.append(("Wilson B", "%.2f" % result.iso_b_wilson))
  return (output_files, stats)

def validate_params(params, callback=None):
  if (params.scaling.input.xray_data.file_name is None):
    raise Sorry("You must supply a reflection file first!")
  if (params.scaling.input.xray_data.obs_labels is None):
    raise Sorry("Please select labels for input data.")
  elif (params.scaling.input.xray_data.obs_labels ==
        params.scaling.input.xray_data.calc_labels):
    raise Sorry("You may not select the same array labels for both the "+
      "experimental data and F(calc).  (It is okay to run Xtriage without "+
      "F(calc) data if it is not available - most analyses do not require "+
      "this information.)")
  d_min = params.scaling.input.xray_data.high_resolution
  d_max = params.scaling.input.xray_data.low_resolution
  if (d_min is not None) and (d_max is not None):
    if d_min > d_max :
      raise Sorry(("The specified high-resolution cutoff (%.3f) is greater "+
        "than the low-resolution cutoff (%.3f)!") % (d_min, d_max))
  if ((params.scaling.input.xray_data.space_group is None) or
      (params.scaling.input.xray_data.unit_cell is None)):
    raise Sorry("Missing or incomplete symmetry information.")
  space_group = params.scaling.input.xray_data.space_group.group()
  unit_cell = params.scaling.input.xray_data.unit_cell
  if (not space_group.is_compatible_unit_cell(unit_cell)):
    raise Sorry("Unit cell parameters are not consistent with the "+
        "currently set space group.  Please make sure that the symmetry "+
        "information is entered correctly.")
  ref_structure = params.scaling.input.xray_data.reference.structure.file_name
  if (ref_structure is not None) and (not os.path.isfile(ref_structure)):
    raise Sorry(("A path was defined for the reference PDB file (%s), but it "+
      "could not be recognized as a file or does not exist.") % ref_structure)
  if (params.scaling.input.asu_contents.sequence_file is not None):
    # XXX this won't actually crash (xtriage will just ignore the sequence
    # file) but it might be better to avoid user confusion
    if ((params.scaling.input.asu_contents.n_residues is not None) or
        (params.scaling.input.asu_contents.n_bases is not None)):
      raise Sorry("Please leave the number of residues and bases blank when "+
        "providing a sequence file.")
  if (params.scaling.gui.output_dir is None):
    raise Sorry("Please specify an output directory.")
  twin_params = params.scaling.input.optional.symmetry.twinning_parameters
  if (twin_params.fraction is not None) and (twin_params.fraction >= 0.5):
    raise Sorry("The twin fraction (if defined) must be less than 0.5.")
  return True

########################################################################
# XXX BACKWARDS COMPATIBILITY
class xtriage_summary(object):
  """
  Old result class, minus initialization.  Provides backwards compatibility
  with pickle files from Phenix 1.9 and earlier.
  """
  new_format = False
  def get_relative_wilson(self):
    if hasattr(self, "rel_wilson_caption"):
      return (self.rel_wilson_caption, self.rel_wilson_plot)
    else :
      return (None, None)

  def get_merging_statistics(self):
    return getattr(self, "merging_stats", None)

  def get_data_file(self):
    return getattr(self, "data_file", None)

  def original_intensities_flag(self):
    return getattr(self, "original_is_intensity_array", None)

  def get_completeness(self):
    overall = getattr(self, "completeness_overall", None)
    binned = getattr(self, "completeness_binned", None)
    return (overall, binned)

  def is_centric(self):
    return getattr(self, "centric_flag", None)
