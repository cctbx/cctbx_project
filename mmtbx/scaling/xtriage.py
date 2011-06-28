from cctbx import crystal
from cctbx.array_family import flex
from libtbx.utils import \
  Sorry, show_exception_info_if_full_testing, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx import pdb
from iotbx.option_parser import option_parser
import mmtbx.utils
import mmtbx.scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import basic_analyses, pair_analyses
from mmtbx.scaling import massage_twin_detwin_data
import libtbx.phil
from libtbx.str_utils import StringIO
from libtbx.utils import null_out
from libtbx import runtime_utils
import sys, os


master_params = iotbx.phil.parse("""\
scaling {
input {
  asu_contents
    .help = "Defines the ASU contents"
    .short_caption = ASU contents
    .style = menu_item auto_align
  {
    n_residues=None
      .type=float
      .help="Number of residues in structural unit"
      .short_caption = Number of residues in asymmetric unit
    n_bases=None
      .type=float
      .help="Number of nucleotides in structural unit"
      .short_caption = Number of nucleotides in asymmetric unit
    n_copies_per_asu=None
      .type=float
      .help="Number of copies per ASU. If not specified, Matthews analyses is performed"
      .short_caption = Number of copies in asymmetric unit
  }

  xray_data
   .help="Defines xray data"
    .style = auto_align
  {
    file_name=None
      .type=path
      .help="File name with data"
      .short_caption = Reflections
      .style = bold noauto file_type:hkl process_hkl update_d_max_min \
        child:fobs:obs_labels child:fcalc:calc_labels \
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
          .style = bold file_type:hkl process_hkl update_d_max_min \
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
         .style = file_type:pdb noauto
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
       ccp4_style_graphs=True
       .type=bool
       .help="SHall we include ccp4 style graphs?"
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
       apply_basic_filters_prior_to_twin_analysis=True
         .type=bool
         .help="Keep data cutoffs from the basic_analyses module (I/sigma,Wilson scaling,Anisotropy) when twin stats are computed."

     }
   }

   optional
    .expert_level=1
    .help="Optional data massage possibilities"
    .short_caption = Advanced data massaging options
    .style = menu_item auto_align
   {
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
  include scope libtbx.phil.interface.tracking_params
}
}
""", process_includes=True)


def print_banner(appl, out=None):
  if out is None:
    out=sys.stdout
  hashes = "#############################################################"
  def print_centered(s):
    b = max(0, len(hashes) - len(s) - 4)
    l = int(b // 2)
    r = b - l
    print >> out, "##%s%s%s##" % (" "*l, s, " "*r)
  print >> out, hashes
  print_centered(appl)
  print_centered("")
  print_centered("P.H. Zwart, R.W. Grosse-Kunstleve & P.D. Adams")
  print_centered("")
  print >> out, hashes


def print_help(appl):
  print """
----------------------------------------------------------------------------------------------------------------------------
Usage: %(appl)s file_name=myfile.sca <options>

%(appl)s performs a variety of twinning and related test on the given
x-ray data set. See CCP4 newletter number 42, July 2005 for more information.
(http://www.ccp4.ac.uk/newsletters/newsletter42/newsletter42.pdf)


The program options are summarized below

1. scope: asu_contents
   keys: * n_residues :: Number of residues per monomer/unit
         * n_bases :: Number of nucleotides per monomer/unit
         * n_copies_per_asu :: Number of copies in the ASU.

   These keywords control the determination of the absolute scale.
   If the number of residues/bases is not specified, a solvent content of 50%% is assumed.


2a.scope: xray_data
   keys: * file_name :: file name with xray data.
         * obs_labels :: labels for observed data is format is mtz of XPLOR/CNS
         * calc_labels :: optional; labels for calculated data
         * unit_cell :: overrides unit cell in reflection file (if present)
         * space_group :: overrides space group in reflection file (if prersent)
         * high_resolution :: High resolution limit of the data
         * low_resolution :: Low resolution limit of the data

   Note that the matching of specified and present labels involves a sub-string matching
   algorithm. See 'Example usage'.


2b.scope: xray_data.reference : A reference data set or structure.
   keys:  data.file_name :: file name for xray data
          structure.file_name :: file name of a PDB.
                                 Specification of a reference structure triggers RvsR cacluations.

3. scope: parameters.misc_twin_parameters.missing_symmetry
   keys: * sigma_inflation :: Sigma inflation value in scoring function.

   sigma_intensity_in_scoring_function = sigma_intensity*sigma_inflation

   Larger values of sigma inflation will tend to result in the point group selection
   algorithm to favour higher symmetry. If data is processed reasonably, the default
   should be fine.

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
   The value determined by the automatic procedure can be overruled by specification
   of the high_resolution keyword. The low resolution is set to 10A by default.


6. scope: parameters.reporting
   keys: * verbose :: verbosity level.
         * log :: log file name
         * ccp4_style_graphs :: Either True or False. Determines whether or not
                                ccp4 style logfgra plots are written to the log file


7. scope: xray_data
   keys: * file_name :: file name with xray data.
         * obs_labels :: labels for observed data is format is mtz of XPLOR/CNS
         * calc_labels :: optional; labels for calculated data
         * unit_cell :: overrides unit cell in reflection file (if present)
         * space_group :: overrides space group in reflection file (if prersent)
         * high_resolution :: High resolution limit of the data
         * low_resolution :: Low resolution limit of the data

   Note that the matching of specified and present labels involves a sub-string matching
   algorithm. See 'Example usage'.


8. scope: optional
   keys: * hklout :: output mtz or sca file
         * hklout_type :: sca, mtz or mtz_or_sca.
                          mtz_or_sca auto detects the format on the basis of the extension of hklout.
         * aniso.action :: remove_aniso or not
         * aniso.final_b :: the final b after anisotropy correction.
                            Choices are the mean or smallest eigenvalue of B-cart.
                            A user supplied B-value can be chosen as well by selection user_b_iso
                            and specifying aniso.b_iso
         * symmetry.action :: Whether to twin, detwin or leave the data alone (do nothing)
         * symmetry.twinning_parameters.twin_law :: using this twin law (h,k,l or  a,b,c or x,y,z notation)
         * symmetry.twinning_parameters.fraction :: The detwinning fraction.
         * outlier.action :: what type of outlier rejection to perform.
                             extreme: uses extreme value statistics
                             basic: uses normal wilson statistics
                             beamstop: only rejects very weak low resolution reflections
                             specific parameters can be set in the outlier.parameters scope.
                             Defaults should be fine.

   This section controls a set of processes that allows the user to perform outlier rejection,
   anisotropy correction and twinning/detwinning. It is not dependent on the results from the xtriage
   analyses.
   These options have an associated expert level of 10, and are not shown by default. Specification
   of the expert level on the command line as 'level=100' will show all available options.



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


-----------------------------------------------------------------------------------------------------------------------------

""" % vars()

class xtriage_analyses(object):
  def __init__(self,
               miller_obs,
               miller_calc = None,
               miller_ref  = None,
               parameters  = None,
               text_out    = None,
               plot_out    = None):
    self.miller_obs  = miller_obs   # array of observed data, should be intensity or amplitude
    self.miller_calc = miller_calc  # array if calculated data, need to be given
    self.miller_ref  = miller_ref   # array with 'reference' data, need not be given.

    if self.miller_obs is not None:
      self.miller_obs  = self.miller_obs.merge_equivalents().array().remove_systematic_absences()   # array of observed data, should be intensity or amplitude
    if self.miller_calc is not None:
      self.miller_calc = self.miller_calc.merge_equivalents().array().remove_systematic_absences()  # array if calculated data, need to be given
    if self.miller_ref is not None:
      self.miller_ref  = self.miller_ref.merge_equivalents().array().remove_systematic_absences()   # array with 'reference' data, need not be given.
                                    # A reference set is for instance a data set with an alternative indexing scheme
    self.text_out    = text_out     # An object with a write method, such as a multi out or StringIO.
                                    # If None, sys.stdout will be used
    if self.text_out == "silent":   # if "silent", a StringIO object will be used
     self.text_out = null_out()     # and all output is supressed

    self.plot_out    = plot_out     # as above. This will contain some ccp4 style plots. If None, no plots will be made

    self.params = parameters        # this should be a phil object like the master_params define on the top of this file
    if self.params is None:         # if nothing specified, defaults are used
      self.params = master_params.fetch(sources=[])
      self.params = self.params.extract()





    print >> self.text_out
    print >> self.text_out,"##----------------------------------------------------##"
    print >> self.text_out,"##                    Basic statistics                ##"
    print >> self.text_out,"##----------------------------------------------------##"
    # Do the basic analyses first please
    self.basic_results = basic_analyses.basic_analyses(
       self.miller_obs,
       self.params,
       out=self.text_out,
       out_plot=self.plot_out,
       miller_calc = miller_calc,
       verbose=1)
    # outliers are removed, make a new copy
    try:
      ma = self.basic_results.miller_array
      if self.params.scaling.input.parameters.misc_twin_parameters.apply_basic_filters_prior_to_twin_analysis:
        self.miller_obs = self.basic_results.miller_array.deep_copy()
        self.normalised_array = self.basic_results.normalised_miller.deep_copy()
        self.params =  self.basic_results.phil_object
    except AttributeError, e:
      print >> self.text_out, "*** ERROR ***"
      print >> self.text_out, str(e)
      show_exception_info_if_full_testing()

    print >> self.text_out
    print >> self.text_out,"##----------------------------------------------------##"
    print >> self.text_out,"##                   Twinning Analyses                ##"
    print >> self.text_out,"##----------------------------------------------------##"
    print >> self.text_out
    print >> self.text_out
    print >> self.text_out

    #Do the twinning analyses
    ## resolution check
    if (flex.min(self.miller_obs.d_spacings().data())
        > self.params.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.high_resolution):
      self.params.scaling.input.xray_data.high_resolution = flex.min(self.miller_obs.d_spacings().data())

    default_high_reso_limit_wilson_ratio = \
      self.params.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.high_resolution
    if default_high_reso_limit_wilson_ratio is None:
      default_high_reso_limit_wilson_ratio = 0.0

    d_star_sq_high_limit = default_high_reso_limit_wilson_ratio
    d_star_sq_high_limit = 1.0/((d_star_sq_high_limit+1e-6)**2.0)

    default_low_reso_limit_wilson_ratio = \
      self.params.scaling.input.parameters.misc_twin_parameters.twin_test_cuts.low_resolution


    d_star_sq_low_limit = default_low_reso_limit_wilson_ratio
    d_star_sq_low_limit = 1.0/((d_star_sq_low_limit+1e-6)**2.0)
    self.twin_results = None
    if(self.miller_obs.select_acentric().as_intensity_array().indices().size()>0):
      self.twin_results = twin_analyses.twin_analyses(
        miller_array=self.miller_obs,
        d_star_sq_low_limit=d_star_sq_low_limit,
        d_star_sq_high_limit=d_star_sq_high_limit,
        normalise=True,
        out=self.text_out,
        out_plots=self.plot_out,
        miller_calc=self.miller_calc,
        additional_parameters=self.params.scaling.input.parameters.misc_twin_parameters)

    if miller_ref is not None:
      self.reference_analyses = pair_analyses.reindexing(
        self.miller_ref,
        self.miller_obs,
        file_name=self.params.scaling.input.xray_data.file_name
      )


def run(args, command_name="phenix.xtriage", return_result=False,
    out=None):
  if (out is None) :
    out = sys.stdout
  command_line = (option_parser(
    usage=command_name+" [options] reflection_file parameters [...]",
    description="Example: %s data1.mtz" % command_name)
    .option(None, "--long_help",
      action="store_true",
      help="show more help and exit")
    .enable_show_defaults()
    .enable_symmetry_comprehensive()
    .option(None, "--weak_symmetry",
      action="store_true",
      default=False,
      help="symmetry on command line is weaker than symmetry found in files")
    .option(None, "--quiet",
      action="store_true",
      help="suppress output")
  ).process(args=args)
  co = command_line.options
  if (len(args) == 0):
    command_line.parser.show_help()
  elif (co.long_help):
    print_help(appl=command_name)
  elif (command_line.expert_level is not None):
    master_params.show(
      expert_level=command_line.expert_level,
      attributes_level=command_line.attributes_level)
  else:
    log = multi_out()
    if (not co.quiet):
      log.register(label="stdout", file_object=out)
    string_buffer = StringIO()
    string_buffer_plots = StringIO()
    log.register(label="log_buffer", file_object=string_buffer)

    print_banner(appl=command_name, out=log)
    print >> log, "#phil __OFF__"
    print >> log
    print >> log, date_and_time()
    print >> log
    print >> log

    phil_objects = []
    argument_interpreter = master_params.command_line_argument_interpreter()

    reflection_file = None

    for arg in command_line.args:
      command_line_params = None
      arg_is_processed = False
      if (os.path.isfile(arg)): ## is this a file name?
        ## Check if this is a phil file
        try:
          command_line_params = iotbx.phil.parse(file_name=arg)
        except KeyboardInterrupt: raise
        except Exception : pass
        if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True

        if (not arg_is_processed):
          ## Check if this file is a reflection file
          if command_line_params is None:
            reflection_file = reflection_file_reader.any_reflection_file(
              file_name=arg, ensure_read_access=False)
          if (reflection_file is not None):
            reflection_file = arg
            arg_is_processed = True

      ## If it is not a file, it must be a phil command
      else:
        command_line_params = argument_interpreter.process(arg=arg)
        phil_objects.append(command_line_params)
        arg_is_processed = True

      if not arg_is_processed:
        print >> log, "##--------------------------------------------------------------------##"
        print >> log, "## Unknown phil-file or phil-command:", arg
        print >> log, "##--------------------------------------------------------------------##"
        print >> log
        raise Sorry("Unknown file format or phil command: %s" % arg)

    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()
    if reflection_file is not None:
      params.scaling.input.xray_data.file_name=reflection_file
    verbose = params.scaling.input.parameters.reporting.verbose

    scope =  params.scaling.input.xray_data
    if (scope.unit_cell is None or not co.weak_symmetry):
        if command_line.symmetry.unit_cell() is not None:
          scope.unit_cell = command_line.symmetry.unit_cell()

    if (scope.space_group is None or not co.weak_symmetry):
      if command_line.symmetry.space_group_info() is not None:
        scope.space_group = command_line.symmetry.space_group_info()

    ## Check for number of residues
    reset_space_group = False
    reset_unit_cell = False

    if params.scaling.input.asu_contents.n_residues is None:
      print >> log, "##-------------------------------------------##"
      print >> log, "## WARNING:                                  ##"
      print >> log, "## Number of residues unspecified            ##"
      print >> log, "##-------------------------------------------##"

    if params.scaling.input.xray_data.file_name is None:
      print >> log,"##-------------------------------------------##"
      print >> log,"## No reflection name is defined.            ##"
      print >> log,"##-------------------------------------------##"
      raise Sorry("No reflection file defined")


    crystal_symmetry = crystal_symmetry_from_any.extract_from(
        file_name=params.scaling.input.xray_data.file_name)
    if crystal_symmetry is None:
      print >> log, "Cell and symmetry not specified in reflection file"
      if params.scaling.input.xray_data.space_group is None:
        raise Sorry("""No space group info available.
  Use keyword 'xray_data.space_group' to specify space group""" )

      if (params.scaling.input.xray_data.unit_cell is None) :
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
        else:
          reset_unit_cell=False


    if params.scaling.input.xray_data.space_group is None:
      params.scaling.input.xray_data.space_group \
        = command_line.symmetry.space_group_info()
    if params.scaling.input.xray_data.space_group is None:
      params.scaling.input.xray_data.space_group \
        = crystal_symmetry.space_group_info()
      reset_space_group = True

    if params.scaling.input.xray_data.unit_cell is None:
      params.scaling.input.xray_data.unit_cell \
        = command_line.symmetry.unit_cell()
    if params.scaling.input.xray_data.unit_cell is None:
      params.scaling.input.xray_data.unit_cell \
        = crystal_symmetry.unit_cell()
      reset_unit_cell = True

    ## Check for unit cell and spacegroup definitions
    if not reset_unit_cell:
      if params.scaling.input.xray_data.unit_cell is not None:
        print >> log,"##-------------------------------------------##"
        print >> log,"## Unit cell defined manually, will ignore"
        print >> log,"## specification in reflection file: "
        if crystal_symmetry is not None:
          print >> log,"## From file : ", \
                crystal_symmetry.unit_cell()
          print >> log,"## From input: ", \
                params.scaling.input.xray_data.unit_cell
        print >> log,"##-------------------------------------------##"

    if not reset_space_group:
      if params.scaling.input.xray_data.space_group is not None:
        print >> log,"##-------------------------------------------##"
        print >> log,"## Space group defined manually, will ignore"
        print >> log,"## specification in reflection file: "
        if crystal_symmetry is not None:
          print >> log,"## From file : ", \
                crystal_symmetry.space_group_info()
          print >> log,"## From input: ", \
                params.scaling.input.xray_data.space_group
        print >> log,"##-------------------------------------------##"

    effective_params = master_params.fetch(sources=phil_objects)
    new_params = master_params.format(python_object=params)
    if (params.scaling.input.parameters.reporting.verbose>0):
      print >> log
      print >> log
      print >> log, "Effective parameters: "
      print >> log, "#phil __ON__"
      new_params.show(out=log,expert_level=params.scaling.input.expert_level)
      print >> log, "#phil __END__"
    crystal_symmetry = crystal.symmetry(
      unit_cell = params.scaling.input.xray_data.unit_cell,
      space_group_symbol = str(params.scaling.input.xray_data.space_group) )

    ## Please check if we have a acentric space group
    if crystal_symmetry.space_group().is_centric() :
      raise Sorry("Centric space groups are not supported.")

    ## Now it time to read in reflection files somehow
    ## We do this via a reflection_file_server

    xray_data_server =  reflection_file_utils.reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      force_symmetry = True,
      reflection_files=[])

    miller_array = None
    miller_array = xray_data_server.get_xray_data(
      file_name = params.scaling.input.xray_data.file_name,
      labels = params.scaling.input.xray_data.obs_labels,
      ignore_all_zeros = True,
      parameter_scope = 'scaling.input.xray_data',
      parameter_name = 'obs_labels'
    )

    info = miller_array.info()

    if not miller_array.is_real_array():
      miller_array = abs( miller_array )
      from cctbx.xray import observation_types
      miller_array = miller_array.set_observation_type( observation_types.amplitude() )



    miller_array = miller_array.map_to_asu()

    if miller_array.observation_type() is None:
      raise Sorry("Observation type of data unkown. Please check input reflection file")

    miller_array = miller_array.select(
      miller_array.indices() != (0,0,0))

    if (miller_array.is_xray_intensity_array()):
      miller_array = miller_array.enforce_positive_amplitudes()
    elif (miller_array.is_complex_array()):
      miller_array = abs(miller_array)

    # first do a low reso cutn if applicable
    if params.scaling.input.xray_data.low_resolution is not None:
      miller_array = miller_array.resolution_filter(d_max=params.scaling.input.xray_data.low_resolution)
    if params.scaling.input.xray_data.high_resolution is not None:
      miller_array = miller_array.resolution_filter(d_min=params.scaling.input.xray_data.high_resolution)

    # make sure sigmas are okai, otherwise, cut them
    if (not miller_array.sigmas_are_sensible()):
      #clearly there is something wrong with the sigmas
      #forget about them I would say
      miller_array = miller_array.customized_copy( indices=miller_array.indices(),
                                                   data=miller_array.data(),
                                                   sigmas=None ).set_observation_type( miller_array )

    miller_array = miller_array.eliminate_sys_absent(integral_only=True, log=log)

    ## Check if Fcalc label is available
    f_calc_miller = None
    f_calc_miller_complex = None
    reference_structure = None
    reference_params = params.scaling.input.xray_data.reference
    if (params.scaling.input.xray_data.calc_labels is not None) :
      f_calc_miller = xray_data_server.get_xray_data(
        file_name = params.scaling.input.xray_data.file_name,
        labels = params.scaling.input.xray_data.calc_labels,
        ignore_all_zeros = True,
        parameter_scope = 'scaling.input.xray_data',
        parameter_name = 'calc_labels'
      )
      if not f_calc_miller.is_real_array():
        f_calc_miller = f_calc_miller.customized_copy(
          data = flex.abs( f_calc_miller.data() ) ).set_observation_type(
            f_calc_miller)
      if (f_calc_miller.is_xray_intensity_array()) :
        print >> log, "Converting %s to amplitudes" % \
          (params.scaling.input.xray_data.calc_labels)
        f_calc_miller = f_calc_miller.f_sq_as_f()
      f_calc_miller = f_calc_miller.eliminate_sys_absent(integral_only=True,
        log=log)
      if (miller_array.anomalous_flag()) :
        if (not f_calc_miller.anomalous_flag()) :
          f_calc_miller = f_calc_miller.generate_bijvoet_mates()
      elif (f_calc_miller.anomalous_flag()) :
        f_calc_miller = f_calc_miller.average_bijvoet_mates()
      f_calc_miller.set_observation_type_xray_amplitude()
    elif (reference_params.structure.file_name is not None):
      if (not os.path.isfile(reference_params.structure.file_name)) :
        raise Sorry("Can't open reference structure - not a valid file.")
      assert f_calc_miller is None
      reference_structure = iotbx.pdb.input(
        file_name=reference_params.structure.file_name).xray_structure_simple(
        crystal_symmetry = miller_array.crystal_symmetry())
      tmp_obj = mmtbx.utils.fmodel_from_xray_structure(
        xray_structure = reference_structure,
        f_obs = miller_array)
      f_calc_miller_complex = tmp_obj.f_model
      f_calc_miller = abs( tmp_obj.f_model ).eliminate_sys_absent(
        integral_only=True,
        log=log).set_observation_type_xray_amplitude()
    twin_results = None


    reference_array = None
    if params.scaling.input.xray_data.reference.data.file_name is not None:
      # do the reference analyses

      # first make a new xray data server
      user_cell=None
      user_group=None
      if params.scaling.input.xray_data.reference.data.unit_cell is not None:
        user_cell = params.scaling.input.xray_data.reference.data.unit_cell
      if params.scaling.input.xray_data.reference.data.space_group is not None:
        user_group = params.scaling.input.xray_data.reference.data.space_group
      reference_symmetry = None
      if user_cell is not None:
        if user_group is not None:
          reference_symmetry= crystal.symmetry(
            unit_cell=user_cell,
            space_group=user_group.group()
            )
      if reference_symmetry is None:
        reference_symmetry = crystal_symmetry_from_any.extract_from(
        file_name=params.scaling.input.xray_data.reference.data.file_name)

      if reference_symmetry is None:
        print >> log, "No reference unit cell and space group could be deduced"
        raise Sorry("Please provide unit cell and space group for reference data")

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
        parameter_name = 'labels'
        )

      info = reference_array.info()

      reference_array = reference_array.map_to_asu()

      reference_array = reference_array.select(
        reference_array.indices() != (0,0,0))

      reference_array = reference_array.select(
        reference_array.data() > 0 )

      if (reference_array.is_xray_intensity_array()):
        reference_array = reference_array.f_sq_as_f()
      elif (reference_array.is_complex_array()):
        reference_array = abs(reference_array)

      reference_array = reference_array.eliminate_sys_absent(integral_only=True, log=log)

    if (not miller_array.is_real_array() ):
      miller_array = abs(miller_array)

    # make sure we hold on to the 'raw' data for later usage is desired
    raw_data = miller_array.deep_copy()

    xtriage_results = None
    if (miller_array.is_real_array()):
      print >> log
      print >> log
      print >> log, "Symmetry, cell and reflection file content summary"
      print >> log
      miller_array.set_info(info=info)
      miller_array.show_comprehensive_summary(f=log)

      minimal_pass = True
      twin_pass = True
      reference_pass = False

      if params.scaling.input.xray_data.reference.data.file_name is not None:
        reference_pass = True

      xtriage_results = xtriage_analyses(miller_obs  = miller_array,
                                         miller_calc = f_calc_miller,
                                         miller_ref  = reference_array,
                                         parameters  = params,
                                         text_out    = log,
                                         plot_out    = string_buffer )

    if params.scaling.input.optional.hklout is not None:

      massage_object = massage_twin_detwin_data.massage_data(
        miller_array=raw_data,
        parameters=params.scaling.input.optional,
        out=log)
      massage_object.write_data()

    ## Append the CCP4i plots to the log StringIO object if desired
    if params.scaling.input.parameters.reporting.ccp4_style_graphs:
      print >> string_buffer, string_buffer_plots.getvalue()
    if(params.scaling.input.parameters.reporting.log is not None):
      output_file = open( params.scaling.input.parameters.reporting.log  ,'w')
      output_file.write(string_buffer.getvalue())

    if return_result :
      summary_out = StringIO()
      miller_array.show_comprehensive_summary(f=summary_out)
      return xtriage_summary(
        params=params,
        xtriage_results=xtriage_results,
        data_summary=summary_out.getvalue())
    else :
      return xtriage_results

#--- Pickle-able results object for the GUI
# TODO: regression tests
# This is *exactly* as gross as it looks.
class xtriage_summary (object) :
  def __init__ (self, params, xtriage_results, data_summary) :
    self.file_name = params.scaling.input.xray_data.file_name
    self.log_file = params.scaling.input.parameters.reporting.log
    self.file_labels = params.scaling.input.xray_data.obs_labels
    self.nresidues = params.scaling.input.asu_contents.n_residues
    self.nbases = params.scaling.input.asu_contents.n_bases
    self.data_summary = data_summary

    #-------------------------------------------------------------------
    # Part 1: basic analyses:
    #         - matthews coefficient/solvent content [ only if nres is known ]
    #         - wilson scaling (isotropic and anisotropic)
    #         - completeness/data strength [ only if I available ]
    #         - low-res completeness ???
    #         - mean intensity [ only if I available ]
    #         - ice-ring analysis
    #         - anomalous measurability [ only if Friedel pairs available ]
    basic_results = xtriage_results.basic_results
    basic_attrs = ["nresidues",
                   "nbases",
                   #"iso_scale_and_b",
                   "iso_p_scale", "iso_b_wilson", # float
                   #"aniso_scale_and_b",
                   "aniso_p_scale", # float
                   "aniso_u_star", "aniso_b_cart", # [ float ] * 6
                   "overall_b_cart"]
    # WILSON SCALING
    for attr in basic_attrs :
      setattr(self, attr, getattr(basic_results, attr, None))
    # COMPLETENESS
    data_strength = getattr(basic_results, "data_strength")
    self.completeness_table = getattr(data_strength, "table_for_gui", None)
    self.completeness_info = getattr(data_strength, "completeness_info", None)
    # WORRISOME SHELLS, MEAN INTENSITY, Z-SCORES/COMPLETENESS,
    # ANOMALOUS SIGNAL, <I/SIGI> BY SHELL
    data_stats_attrs = ["suggested_reso_for_hyss"]
    data_table_attrs = ["shell_table", "wilson_table", "zscore_table",
                        "meas_table", "i_sig_i_table"]
    for attr in (data_stats_attrs + data_table_attrs) :
      setattr(self, attr, getattr(basic_results.basic_data_stats, attr, None))
    self.low_res_info = None
    low_res_completeness = getattr(basic_results.basic_data_stats,
                                   "low_resolution_completeness", None)
    if low_res_completeness is not None :
      out = StringIO()
      low_res_completeness.show(f=out)
      self.low_res_info = out.getvalue()
    meas_anal = getattr(basic_results.basic_data_stats, "meas_anal", None)
    if meas_anal is not None :
      table = getattr(meas_anal, "meas_table", None)
      out = StringIO()
      table.show(f=out)
      self.meas_out = out.getvalue()
    else :
      self.meas_out = None
    self.meas_info = getattr(meas_anal, "message", None)
    self.low_d_cut = getattr(meas_anal, "low_d_cut", None)
    self.high_d_cut = getattr(meas_anal, "high_d_cut", None)
    rel_wilson = getattr(basic_results, "rel_wilson", None)
    if (rel_wilson is not None) :
      caption_out = StringIO()
      rel_wilson.show_summary(out=caption_out)
      self.rel_wilson_caption = caption_out.getvalue()
      self.rel_wilson_plot = rel_wilson.get_data_plot()
    outliers = getattr(basic_results.basic_data_stats, "outlier", None)
    if outliers is not None :
      self.acentric_outliers = outliers.acentric_outliers_table
      self.centric_outliers = outliers.centric_outliers_table
    # VM/%SOLV
    self.nres_known = basic_results.nres_known
    self.matthews_table = basic_results.matthews_results[4] # table
    self.matthews_info = basic_results.matthews_results
    for attr in ["defined_copies", "guessed_copies"] :
      setattr(self, attr, getattr(basic_results, attr, None))
    # ICE RINGS
    ijsco = getattr(basic_results.basic_data_stats, "ijsco")
    self.icy_shells = getattr(ijsco, "icy_shells", None)
    self.ice_warnings = getattr(ijsco, "warnings", 0)
    self.ice_comments = getattr(ijsco, "message", "")

    #-------------------------------------------------------------------
    # Part 2: twinning analyses:
    #         - translational pseudosymmetry
    #         - space group choices
    #         - nz test
    #         - l test
    #         - possible twin laws
    #         - britton plot, h test, murray-rust plot for each twin law
    twin_results = xtriage_results.twin_results
    # SYSTEMATIC ABSENCES AND SPACE GROUP
    abs_sg_anal = getattr(twin_results, "abs_sg_anal", None)
    self.sg_info = getattr(abs_sg_anal, "absence_info", None)
    self.sg_table = getattr(abs_sg_anal, "table_data", None)
    abs_table = getattr(abs_sg_anal, "absences_table", None)
    self.absence_info = getattr(abs_table, "table_text", None)
    self.absence_table = getattr(abs_table, "table_data", None)
    # TWINNING
    twin_attrs = ["nz_test_table", "l_test_table", "twin_law_names",
                  "twin_law_info"]
    for attr in twin_attrs :
      setattr(self, attr, getattr(twin_results, attr, None))
    if self.nz_test_table is not None :
      for attr in ["max_diff_ac", "max_diff_c", "sign_ac", "sign_c",
                   "mean_diff_ac", "mean_diff_c"] :
        setattr(self, "nz_test_"+attr, getattr(twin_results.nz_test,attr,None))
    if self.l_test_table is not None :
      for attr in ["parity_h", "parity_k", "parity_l", "mean_l", "mean_l2",
                   "ml_alpha"] :
        setattr(self, "l_test_"+attr, getattr(twin_results.l_test, attr, None))
    other_attrs = []
    for attr in other_attrs :
      setattr(self, attr, getattr(twin_results, attr, None))
    twin_summary = twin_results.twin_summary
    self.patterson_verdict = twin_summary.patterson_verdict.getvalue()
    self.twinning_verdict = twin_summary.twinning_verdict.getvalue()
    self.twin_law_table= getattr(twin_summary.twin_results, "table_data", None)
    self.z_score_info = getattr(twin_summary.twin_results, "z_score_info",None)
    self.intensity_stats = getattr(twin_summary.twin_results,
                                   "independent_stats", None)
    self.possible_twin_laws = getattr(twin_results, "possible_twin_laws", None)
    self.translation_pseudo_symmetry = getattr(twin_results,
      "translation_pseudo_symmetry", None)
    self.check_sg = getattr(twin_results, "check_sg", None)
    self.suggested_space_group = getattr(twin_results, "suggested_space_group",
                                         None)

  def get_relative_wilson (self) :
    if hasattr(self, "rel_wilson_caption") :
      return (self.rel_wilson_caption, self.rel_wilson_plot)
    else :
      return (None, None)

class launcher (runtime_utils.simple_target) :
  def __call__ (self) :
    return run(args=list(self.args), return_result=True)

def validate_params (params, callback=None) :
  if (params.scaling.input.xray_data.file_name is None) :
    raise Sorry("You must supply a reflection file first!")
  if (params.scaling.input.xray_data.obs_labels is None) :
    raise Sorry("Please select labels for input data.")
  d_min = params.scaling.input.xray_data.high_resolution
  d_max = params.scaling.input.xray_data.low_resolution
  if (d_min is not None) and (d_max is not None) :
    if d_min > d_max :
      raise Sorry(("The specified high-resolution cutoff (%.3f) is greater "+
        "than the low-resolution cutoff (%.3f)!") % (d_min, d_max))
  symm = crystal.symmetry(
    space_group_info=params.scaling.input.xray_data.space_group,
    unit_cell=params.scaling.input.xray_data.unit_cell,
    assert_is_compatible_unit_cell=False)
  if (not symm.is_compatible_unit_cell()):
    raise Sorry("Unit cell parameters are not consistent with the "+
        "currently set space group.  Please make sure that the symmetry "+
        "information is entered correctly.")
  return True

if (__name__ == "__main__") :
  run(sys.argv[1:])
