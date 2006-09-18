from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import basic_analyses, pair_analyses
from mmtbx.scaling import twin_detwin_data
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
import sys, os


master_params = iotbx.phil.parse("""\
scaling.input {
   parameters
   {
     asu_contents
     {
       n_residues=None
       .type=float
       n_bases=None
       .type=float
       n_copies_per_asu=None
       .type=float
     }

     misc_twin_parameters
     {
       missing_symmetry
       {
         tanh_location = 0.08
         .type=float
         tanh_slope = 50
         .type=float
       }

       twinning_with_ncs
       {
         perform_analyses = False
         .type = bool
         n_bins = 7
         .type = int
       }

       twin_test_cuts
       {
         low_resolution=10.0
         .type=float
         high_resolution=None
         .type=float
         isigi_cut=3.0
         .type=float
         completeness_cut=0.85
         .type=float
       }

     }

     reporting
     {
       verbose=1
       .type=int
       log=logfile.log
       .type=str
       ccp4_style_graphs=True
       .type=bool
     }

   }

   xray_data
   {
     file_name=None
     .type=path
     obs_labels=None
     .type=strings
     calc_labels=None
     .type=strings
     unit_cell=None
     .type=unit_cell
     space_group=None
     .type=space_group
     high_resolution=None
     .type=float
     low_resolution=None
     .type=float

     reference{
       data{
         file_name=None
         .type = path
         labels=None
         .type=strings
         unit_cell=None
         .type=unit_cell
         space_group=None
         .type=space_group
       }

     }

   }

   optional
   .expert_level=10
   {
     hklout = None
     .type=path
     twinning{
       action=detwin *None
       .type=choice
       twin_law=None
       .type=str
       fraction=None
       .type=float
     }

   }
   expert_level=1
   .type=int
   .expert_level=10
}
""")


def print_banner(out=None):
  if out is None:
    out=sys.stdout
  print >> out, "#############################################################"
  print >> out, "##                        mmtbx.xtriage                    ##"
  print >> out, "##                                                         ##"
  print >> out, "##      P.H. Zwart, R.W. Grosse-Kunstleve & P.D. Adams     ##"
  print >> out, "##                                                         ##"
  print >> out, "#############################################################"


def print_help():
  print """
----------------------------------------------------------------------------------------------------------------------------
Usage: mmtbx.xtriage file_name=myfile.sca <options>

mmtbx.xtriage performs a variety of twinning and related test on the given
x-ray data set. See CCP4 newletter number 42, July 2005 for more information.
(http://www.ccp4.ac.uk/newsletters/newsletter42/newsletter42.pdf)


The program options are summarized below

1. scope: parameters.asu_contents
   keys: * n_residues :: Number of residues per monomer/unit
         * n_bases :: Number of nucleotides per monomer/unit
         * n_copies_per_asu :: Number of copies in the ASU.

   These keywords control the determination of the absolute scale.
   If the number of residues/bases is not specified, a solvent content of 50% is assumed.


2. scope: parameters.misc_twin_parameters.missing_symmetry
   keys: * tanh_location :: tanh decision rule parameter
         * tanh_slope :: tanh decision rule parameter

   The tanh_location and tanh_slope parameter control what R-value is considered to be
   low enough to be considered a 'proper' symmetry operator. the tanh_location parameter
   corresponds to the inflection point of the approximate step function. Increasing
   tanh_location will result in large R-value thresholds.
   tanh_slope is set to 50 and should be okai.


3. scope: parameters.misc_twin_parameters.twinning_with_ncs
   keys: * perform_test :: can be set to True or False
         * n_bins :: Number of bins in determination of D_ncs

   The perform_test is by default set to False. Setting it to True triggers the
   determination of the twin fraction while taking into account NCS parallel to
   the twin axis.


4. scope: parameters.misc_twin_parameters.twin_test_cuts
   keys: * high_resolution : high resolution for twin tests
         * low_resolution: low resolution for twin tests
         * isigi_cut: I/sig(I) threshold in automatic determination
                      of high resolutiuon limit
         * completeness_cut: completeness threshold in automatic
                             determination of high resolutiuon limit

   The automatic determination of the resolution limit for the twinning test
   is determined on the basis of the completeness after removing intensities for
   which I/sigI<isigi_cut. The lowest limit obtain in this way is 3.5A.
   The value determined by the automatic procedure can be overruled by specification
   of the high_resolution keyword. The low resolution is set to 10A by default.


5. scope: parameters.reporting
   keys: * verbose :: verbosity level.
         * log :: log file name
         * ccp4_style_graphs :: Either True or False. Determines whether or not
                                ccp4 style logfgra plots are written to the log file


6. scope: xray_data
   keys: * file_name :: file name with xray data.
         * obs_labels :: labels for observed data is format is mtz of XPLOR/CNS
         * calc_labels :: optional; labels for calculated data
         * unit_cell :: overrides unit cell in reflection file (if present)
         * space_group :: overrides space group in reflection file (if prersent)
         * high_resolution :: High resolution limit of the data
         * low_resolution :: Low resolution limit of the data

   Note that the matching of specified and present labels involves a sub-string matching
   algorithm. See 'Example usage'.


7. scope: optional
   keys: * hklout :: output mtz file
         * twinning.action :: Whether to detwin the data
         * twinning.twin_law :: using this twin law (h,k,l or x,y,z notation)
         * twinning.fraction :: The detwinning fraction.
         * b_value :: the resulting Wilson B value

   The output mtz file contains an anisotropy corrected mtz file, with suspected outliers removed.
   The data is put scaled and has the specified Wilson B value.
   These options have an associated expert level of 10, and are not shown by default. Specification
   of the expert level on the command line as 'level=100' will show all available options.



Example usage:

  The commands
    mmtbx.xtriage xray_data.file_name=my_refl.mtz
    mmtbx.xtriage file_name=my_refl.mtz
    mmtbx.xtriage file=my_refl.mtz
    mmtbx.xtriage my_refl.mtz
  are equivalent.

  The commands
    mmtbx.xtriage my_refl.mtz obs=F_HG,SIGF_HG data.high=2.0 log=log.log perform=True
    mmtbx.xtriage my_refl.mtz obs=HG data.high=2.0 log=log.log perform=True
  are equivalent if the substring 'HG' is unique in the mtz file. If the labels contain character such as
  '(', use quoation marks: obs='F_HG(+)'.


-----------------------------------------------------------------------------------------------------------------------------

"""

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
                                    # A reference set is for instance a data set with an alternative indexing scheme
    self.text_out    = text_out     # An object with a write method, such as a multi out or StringIO.
                                    # If None, sys.stdout will be used
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
       out_plot=self.plot_out)
    # outliers are removed, make a new copy
    self.miller_obs = self.basic_results.miller_array.deep_copy()
    self.normalised_array = self.basic_results.normalised_miller.deep_copy()
    self.params =  self.basic_results.phil_object

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
      self.params.scaling.input.xray_data.resolution = flex.min(miller_array.d_spacings().data())

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

    self.twin_results = twin_analyses.twin_analyses(
      self.miller_obs,
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


def run(command_name, args):

  if len(args)==0:
    print_help()
  elif ( "--help" in args ):
    print_help()
  elif ( "--h" in args ):
    print_help()
  elif ("-h" in args ):
    print_help()
  else:
    log = multi_out()
    if (not "--quiet" in args):
      log.register(label="stdout", file_object=sys.stdout)
    string_buffer = StringIO()
    string_buffer_plots = StringIO()
    log.register(label="log_buffer", file_object=string_buffer)

    print_banner(log)
    print >> log, "#phil __OFF__"
    print >> log, "  This cryptic code, together with the tags __ON__ and __OFF__"
    print >> log, "  allows one to use the log file as an input file for xtriage."
    print >> log, "  Try : mmtbx.xtriage  <logfile> to give it a try!"
    print >> log
    print >> log, date_and_time()
    print >> log
    print >> log

    phil_objects = []
    argument_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_params=master_params,
      home_scope="xtriage")

    reflection_file = None

    for arg in args:
      print phil_objects
      command_line_params = None
      arg_is_processed = False
      if arg == '--quiet':
        arg_is_processed = True
        ## The associated action with this keyword is implemented above
      if (os.path.isfile(arg)): ## is this a file name?
        ## Check if this is a phil file
        try:
          command_line_params = iotbx.phil.parse(file_name=arg)
        except KeyboardInterrupt: raise
        except : pass
        if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        ## Check if this file is a reflection file
        if command_line_params is None:
          reflection_file = reflection_file_reader.any_reflection_file(
            file_name=arg, ensure_read_access=False)
        if (reflection_file is not None):
          reflection_file = arg
          arg_is_processed = True

        ## Try and see if it is a phil file
        try:
          command_line_params = iotbx.phil.parse(file_name=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except : pass
      ## If it is not a file, it must be a phil command
      else:
        try:
          command_line_params = argument_interpreter.process(arg=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except : pass



      if not arg_is_processed:
        print >> log, "##----------------------------------------------##"
        print >> log, "## Unknown phil-file or phil-command:", arg
        print >> log, "##----------------------------------------------##"
        print >> log
        raise Sorry("Unknown file format or phil command: %s" % arg)

    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()
    if reflection_file is not None:
      params.scaling.input.xray_data.file_name=reflection_file
    verbose = params.scaling.input.parameters.reporting.verbose




    ## Check for number of residues
    reset_space_group = False
    reset_unit_cell = False

    if params.scaling.input.parameters.asu_contents.n_residues is None:
      print >> log, "##-------------------------------------------##"
      print >> log, "## WARNING:                                  ##"
      print >> log, "## Number of residues unspecified            ##"
      print >> log, "##-------------------------------------------##"

    if params.scaling.input.xray_data.file_name is None:
      print >> log,"##-------------------------------------------##"
      print >> log,"## No reflection name is defined.            ##"
      print >> log,"##-------------------------------------------##"
      raise Sorry("No reflection file defined")


    crystal_symmetry = None
    crystal_symmetry = crystal_symmetry_from_any.extract_from(
        file_name=params.scaling.input.xray_data.file_name)
    if crystal_symmetry is None:
      print >> log, "Cell and symmetry not specified in reflection file"
      if params.scaling.input.xray_data.space_group is None:
        raise Sorry("""
No space group info available.
Use keyword 'xray_data.space_group' to specify space group
                    """ )

      if params.scaling.input.xray_data.unit_cell is None:
        raise Sorry("""
No unit cell info available.
Use keyword 'xray_data.unit_cell' to specify unit_cell
                    """ )


    if params.scaling.input.xray_data.space_group is None:
      params.scaling.input.xray_data.space_group \
        = crystal_symmetry.space_group_info()
      reset_space_group = True

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
    assert( not crystal_symmetry.space_group().is_centric())

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

    miller_array = miller_array.map_to_asu()

    miller_array = miller_array.select(
      miller_array.indices() != (0,0,0))

    miller_array = miller_array.select(
      miller_array.data() > 0 )

    if (miller_array.is_xray_intensity_array()):
      miller_array = miller_array.f_sq_as_f()
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

    ## Check if Fcalc label is available
    f_calc_miller = None
    if params.scaling.input.xray_data.calc_labels is not None:
      f_calc_miller = xray_data_server.get_xray_data(
        file_name = params.scaling.input.xray_data.file_name,
        labels = params.scaling.input.xray_data.calc_labels,
        ignore_all_zeros = True,
        parameter_scope = 'scaling.input.xray_data',
        parameter_name = 'calc_labels'
      )
      if not f_calc_miller.is_real_array():
        f_calc_miller = f_calc_miller.customized_copy(
          data = flex.abs( f_calc_miller.data() ) ).set_observation_type(f_calc_miller  )

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
        user_group = params.scaling.input.xray_data.reference.data.unit_group

      reference_symmetry = None
      if user_cell is not None:
        if user_group is not None:
          reference_symmetry= crystal.symmetry(
            unit_cell=user_cell,
            space_group=user_group
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

      reference_array = miller_array.map_to_asu()

      reference_array = reference_array.select(
        reference_array.indices() != (0,0,0))

      reference_array = reference_array.select(
        reference_array.data() > 0 )

      if (reference_array.is_xray_intensity_array()):
        reference_array = reference_array.f_sq_as_f()
      elif (reference_array.is_complex_array()):
        reference_array = abs(reference_array)

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
      if params.scaling.input.optional.twinning.action == "detwin":
        twin_detwin_data.detwin_data(miller_array,
                                     params,
                                     log)

    ## Append the CCP4i plots to the log StringIO object if desired
    if params.scaling.input.parameters.reporting.ccp4_style_graphs:
      print >> string_buffer, string_buffer_plots.getvalue()

    output_file = open( params.scaling.input.parameters.reporting.log  ,'w')
    output_file.write(string_buffer.getvalue())





if (__name__ == "__main__"):
  run(sys.argv[1:])
