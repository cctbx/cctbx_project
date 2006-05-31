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
from mmtbx.scaling import basic_analyses
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
     .expert_level=1
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

     }

   }
   analyses{
     verbose = 1
     .type=int
     log = logfile.log
     .type = str
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
     low_resolution_for_twin_tests=10.0
     .type=float
     high_resolution_for_twin_tests=None
     .type=float
     isigi_cut=3.0
     .type=float
     completeness_cut=0.85
     .type=float
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


The program options are defined in the following table.

  -------------------------------------------------------------------------------------
  | full parameter name                      |  meaning                               |
  -------------------------------------------------------------------------------------
  | basic.n_residues                         | Number of residues per monomer         |
  | basic.n_bases                            | Number of base pairs per monomer       |
  | basic.n_copies_per_asu                   | Number of copies in the ASU            |
  | analyses.verbose                         | Verbosity level                        |
  | analyses.log                             | Log file with CCP4i style plots        |
  | xray_data.file_name                      | File containing xray data              |
  | xray_data.obs_labels                     | Labels for experimental data           |
  | xray_data.calc_labels                    | Labels for calculated data             |
  | xray_data.unit_cell                      | Unit cell constants                    |
  | xray_data.space_group                    | Space group                            |
  | xray_data.high_resolution                | high resolution limit in wilson scaling|
  | xray_data.low_resolution                 | low resolution limit in wilson scaling |
  | xray_data.high_resolution_for_twin_tests | high resolution limit in int. stats.   |
  | xray_data.low_resolution_for_twin_tests  | low resolution limit in int. stats.    |
  | xray_data.isigi_cut                      | I/sig(I) cut in int. stats.            |
  | xray_data.completenesss_cut              | completeness cut int. stats.           |
  -------------------------------------------------------------------------------------


Detailed description:


basic.n_residues: Number of residues in monomer. If none, value is estimated
                  assuming 50% solvent content.

basic.n_bases: Number of bases in monomer. Can be used with basic.n_residues

basic.n_copies_per_asu : Number of copies in the ASU. if None, value is estimated
                         and used. If specified, value is estimated and estimate
                         is reported for your information only.

analyses.verbose : verbosity level of the output;  0 -> very brief/silent
                                                   1 -> normal

analyses.log: The name of the output log file containing CCP4I
              style loggraph plots and the same info that is dumped to the
              screen

xray_data.file_name: file containing reflection information.
                     Note that if file is in MTZ or CNS format,
                     the 'IP' or 'FP' column SHOULD be next to the
                     'SIGIP' or 'SIGFP' column.
                     For shelx files, one must specify if intensities or
                     amplitudes are in the file. See appropriate error message
                     / warning when this happens.

xray_data.obs_labels: a substring of a label indentifying the experimental data of interest.
                      For instance, if the colum labels present in the file are
                      F_CRYSTAL1,SIGF_CRYSTAL1,F_CRYSTAL2,SIGF_CRYSTAL2
                      the label variable can be set to L1 to select
                      F_CRYSTAL1,SIGF_CRYSTAL1
                      If the label name contains parenthesis, uses quotes around the (sub)string
                      to avoid confusion.

xray_data.calc_labels: A substring of a label idenitifying calculated data.
                       (Used for R vs R statistic.).
                       Calculated data needs to be in the same file as in the observed data!

xray_data.unit_cell: By default, the values in the reflection file are used.
                     These keywords must be specified for CNS or SHELX files
                     that (unfortunately) do not have this information
                     in the file.

xray_data.space_group: See xray_data.unit_cell.

xray_data.high_resolution: High Resolution limit where the data is cut.
                           The default is the highest resolution limit
                           available.

xray_data.low_resolution: Low resolution limit where the data is cut.


xray_data.high_resolution_for_twin_tests: High Resolution limit where the data is cut for the twin analyses.
                                          if None is specified, the xtriage will dtermin this limit
                                          on the basis of the parameters isigi_cut and completeness_cut

xray_data.low_resolution_for_twin_test: Low resolution limit where the data is cut for the twin analyses.
                                         The default is 10 A.

xray_data.isigi_cut and xray_data.completenesss_cut: Intensities for which I/sigI<isigi_cut are removed from the dataset.
                                                     The resolution limit for which the completeness is largern than
                                                     xray_data.completenesss_cut, is used as the high resolution limit
                                                     in the intensity statistics and twin analyses.


Example usage:

  The commands
    mmtbx.xtriage xray_data.file_name=the_ribozyme.sca
    mmtbx.xtriage file_name=the_ribozyme.sca
    mmtbx.xtriage file=the_ribozyme.sca
    mmtbx.xtriage the_ribozyme.sca
  are equivalent.

  The command
    mmtbx.xtriage the_ribozyme.mtz labels=F_HG,SIGF_HG reso=2.0 log=log.log
  will use the columns F_HG and SIGF_HG in the analyses. The resolution cut
  for twinning analyses is set at 2.0 A, the logfile with plot will be written
  to log.log.

-----------------------------------------------------------------------------------------------------------------------------

"""

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
    print >> log, date_and_time()
    print >> log
    print >> log

    phil_objects = []
    argument_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_params=master_params,
      home_scope="xtriage")

    reflection_file = None

    for arg in args:
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
    verbose = params.scaling.input.analyses.verbose




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
Use keyword 'space_group' to specify space group
                    """ )

      if params.scaling.input.xray_data.unit_cell is None:
        raise Sorry("""
No unit cell info available.
Use keyword 'unit_cell' to specify unit_cell
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
    if (params.scaling.input.analyses.verbose>0):
      print >> log
      print >> log
      print >> log, "Effective parameters: "
      new_params.show(out=log,expert_level=params.scaling.input.expert_level)

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


    if (miller_array.is_real_array()):
      print >> log
      print >> log
      print >> log, "Symmetry, cell and reflection file content summary"
      print >> log
      miller_array.set_info(info=info)
      miller_array.show_comprehensive_summary(f=log)

      minimal_pass = True
      twin_pass = True

      print >> log
      print >> log,"##----------------------------------------------------##"
      print >> log,"##                    Basic statistics                ##"
      print >> log,"##----------------------------------------------------##"
      if minimal_pass: ## this is a bit silly as it always should happen

        basic_results = basic_analyses.basic_analyses(
          miller_array,
          params,
          out=log,
          out_plot=string_buffer_plots,
          verbose=verbose)
        miller_array = basic_results.miller_array.deep_copy()
        normalised_array = basic_results.normalised_miller.deep_copy()
        params =  basic_results.phil_object

      if twin_pass:
        print >> log
        print >> log,"##----------------------------------------------------##"
        print >> log,"##                   Twinning Analyses                ##"
        print >> log,"##----------------------------------------------------##"
        print >> log

        print >> log
        print >> log
      ## resolution check
        if (flex.min(miller_array.d_spacings().data())
            > params.scaling.input.xray_data.high_resolution_for_twin_tests):
          params.scaling.input.xray_data.resolution = \
             flex.min(miller_array.d_spacings().data())

        default_high_reso_limit_wilson_ratio = \
           params.scaling.input.xray_data.high_resolution
        if default_high_reso_limit_wilson_ratio is None:
          default_high_reso_limit_wilson_ratio = 0.0

        d_star_sq_high_limit = default_high_reso_limit_wilson_ratio
        d_star_sq_high_limit = 1.0/((d_star_sq_high_limit+1e-6)**2.0)

        default_low_reso_limit_wilson_ratio = \
          params.scaling.input.xray_data.low_resolution_for_twin_tests


        d_star_sq_low_limit = default_low_reso_limit_wilson_ratio
        d_star_sq_low_limit = 1.0/((d_star_sq_low_limit+1e-6)**2.0)
        #try:
        twin_results = twin_analyses.twin_analyses(
            miller_array,
            d_star_sq_low_limit=d_star_sq_low_limit,
            d_star_sq_high_limit=d_star_sq_high_limit,
            normalise=True,
            out=log,
            out_plots=string_buffer_plots,
            verbose=verbose,
            miller_calc=f_calc_miller,
            additional_parameters=params.scaling.input.parameters.misc_twin_parameters)
        #except: pass

    if params.scaling.input.optional.hklout is not None:
      if params.scaling.input.optional.twinning.action == "detwin":
        twin_detwin_data.detwin_data(miller_array,
                                     params,
                                     log)

    ## Append the CCP4i plots to the log StringIO object.
    print >> string_buffer, string_buffer_plots.getvalue()
    output_file = open( params.scaling.input.analyses.log  ,'w')
    output_file.write(string_buffer.getvalue())





if (__name__ == "__main__"):
  run(sys.argv[1:])
