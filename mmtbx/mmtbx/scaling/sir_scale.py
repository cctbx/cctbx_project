from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling, relative_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import basic_analyses, data_statistics,pair_analyses
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
import sys, os


master_params = iotbx.phil.parse("""\
scaling.input {
   basic{
     n_residues=None
     .type=float
     n_bases=None
     .type=float
     n_copies_per_asu=None
     .type=float
     n_terms = None
     .type=int
   }

   SIR_scale
   {

     target=*ml ls loc
     .type=choice

     xray_data{

       unit_cell=None
       .type=unit_cell

       space_group=None
       .type=space_group

       high_resolution=None
       .type=float
       low_resolution=None
       .type=float

       native{
         file_name=None
         .type=path
         labels=None
         .type=strings
       }
       derivative{
         file_name=None
         .type=path
         labels=None
         .type=strings
       }

     }

   }


   output{
     log = 'logfile.log'
     .type = path

     root = SIR
     .type = strings

     format = mtz
     .type = strings

   }


}
""")


def run(args):

  if len(args)==0:
    print "no help available"
  elif ( "--help" in args ):
    print "no help available"
  elif ( "--h" in args ):
    print "no help available"
  else:
    log = multi_out()
    if (not "--quiet" in args):
      log.register(label="stdout", file_object=sys.stdout)
    string_buffer = StringIO()
    string_buffer_plots = StringIO()
    log.register(label="log_buffer", file_object=string_buffer)

    log_plots = StringIO()
    print >> log,"#phil __OFF__"
    print >> log
    print >> log, date_and_time()
    print >> log
    print >> log

    phil_objects = []
    argument_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_params=master_params,
      home_scope="scaling")

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


    ## Now please read in the reflections files

    ## get symmetry and cell data first please
    ## By default, the native cell and symmetry are used
    ## as reference
    crystal_symmetry_nat = None
    crystal_symmetry_nat = crystal_symmetry_from_any.extract_from(
      file_name=params.scaling.input.SIR_scale.xray_data.native.file_name)

    if params.scaling.input.SIR_scale.xray_data.space_group is None:
      params.scaling.input.SIR_scale.xray_data.space_group =\
        crystal_symmetry_nat.space_group_info()
      print >> log, "Using symmetry of native data"

    if params.scaling.input.SIR_scale.xray_data.unit_cell is None:
      params.scaling.input.SIR_scale.xray_data.unit_cell =\
        crystal_symmetry_nat.unit_cell()
      print >> log, "Using cell of native data"

    ## Check if a unit cell is defined
    if params.scaling.input.SIR_scale.xray_data.space_group is None:
      raise Sorry("No space group defined")
    if params.scaling.input.SIR_scale.xray_data.unit_cell is None:
      raise Sorry("No unit cell defined")


    crystal_symmetry = crystal_symmetry = crystal.symmetry(
      unit_cell =  params.scaling.input.SIR_scale.xray_data.unit_cell,
      space_group_symbol = str(
        params.scaling.input.SIR_scale.xray_data.space_group) )


    effective_params = master_params.fetch(sources=phil_objects)
    new_params = master_params.format(python_object=params)
    print >> log, "Effective parameters"
    print >> log, "#phil __ON__"
    new_params.show(out=log)
    print >> log, "#phil __END__"
    print >> log

    ## define a xray data server
    xray_data_server =  reflection_file_utils.reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      force_symmetry = True,
      reflection_files=[])

    ## Read in native data and make appropriatre selections
    miller_array_native = None
    miller_array_native = xray_data_server.get_xray_data(
      file_name = params.scaling.input.SIR_scale.xray_data.native.file_name,
      labels = params.scaling.input.SIR_scale.xray_data.native.labels,
      ignore_all_zeros = True,
      parameter_scope = 'scaling.input.SIR_scale.xray_data.native'
    )
    info_native = miller_array_native.info()
    miller_array_native=miller_array_native.map_to_asu().select(
      miller_array_native.indices()!=(0,0,0) )
    miller_array_native = miller_array_native.select(
      miller_array_native.data() > 0 )
    ## Convert to amplitudes
    if (miller_array_native.is_xray_intensity_array()):
      miller_array_native = miller_array_native.f_sq_as_f()
    elif (miller_array_native.is_complex_array()):
      miller_array_native = abs(miller_array_native)
    if not miller_array_native.is_real_array():
      raise Sorry("miller_array_native is not a real array")
    miller_array_native.set_info(info = info_native)
    ## Print info
    print >> log
    print >> log, "Native data summary"
    print >> log, "-------------------"
    miller_array_native.show_comprehensive_summary(f=log)
    print >> log
    ## Do a simple outlier analyses please
    native_outlier = data_statistics.possible_outliers(miller_array_native)
    native_outlier.show(out=log)
    miller_array_native = native_outlier.remove_outliers(
      miller_array_native )

    ## Read in derivative data and make appropriate selections
    miller_array_derivative = None
    miller_array_derivative = xray_data_server.get_xray_data(
      file_name = params.scaling.input.SIR_scale.xray_data.derivative.file_name,
      labels = params.scaling.input.SIR_scale.xray_data.derivative.labels,
      ignore_all_zeros = True,
      parameter_scope = 'scaling.input.SIR_scale.xray_data.derivative'
    )
    info_derivative = miller_array_derivative.info()
    miller_array_derivative=miller_array_derivative.map_to_asu().select(
      miller_array_derivative.indices()!=(0,0,0) )
    miller_array_derivative = miller_array_derivative.select(
      miller_array_derivative.data() > 0 )
    ## Convert to amplitudes
    if (miller_array_derivative.is_xray_intensity_array()):
      miller_array_derivative = miller_array_derivative.f_sq_as_f()
    elif (miller_array_derivative.is_complex_array()):
      miller_array_derivative = abs(miller_array_derivative)
    if not miller_array_derivative.is_real_array():
      raise Sorry("miller_array_derivative is not a real array")
    miller_array_derivative.set_info(info = info_derivative)
    ## Print data summaries
    print >> log
    print >> log
    print >> log, "Derivative data"
    print >> log, "---------------"
    miller_array_derivative.show_comprehensive_summary(f=log)
    print >> log
    ## Do a simple outlier analyses please
    derivative_outlier = data_statistics.possible_outliers(
      miller_array_derivative)
    derivative_outlier.show(out=log)
    miller_array_derivative = derivative_outlier.remove_outliers(
      miller_array_derivative )


    ## Determine the asu contents
    print >> log, "Determining ASU contents"
    print >> log, "------------------------"
    matthews_analyses =matthews.matthews_rupp(
      miller_array = miller_array_native,
      n_residues = params.scaling.input.basic.n_residues,
      n_bases = params.scaling.input.basic.n_bases,
      out=log, verbose=1)

    ## Now the asu contnets is defined, we can do the 'absolute scaling'
    ## to get the data sets on roughly the same scale
    n_residues=matthews_analyses[0]
    n_bases=matthews_analyses[1]
    n_copies_solc=matthews_analyses[2]

    if ((n_residues+n_bases)*n_copies_solc<=0):
      raise Sorry("Empty ASU, check inputs")
    print >> log
    print >> log, "Anisotropic absolute scaling of native"
    print >> log, "--------------------------------------"
    print >> log

    aniso_native = absolute_scaling.ml_aniso_absolute_scaling(
      miller_array = miller_array_native,
      n_residues = n_residues*\
        miller_array_native.space_group().order_z()*n_copies_solc,
      n_bases = n_bases*\
        miller_array_native.space_group().order_z()*n_copies_solc)
    aniso_native.show(out=log,verbose=1)
    print >> log, "removing anisotropy for native"
    u_star_correct_nat = aniso_native.u_star
    miller_array_native = absolute_scaling.anisotropic_correction(
      miller_array_native,aniso_native.p_scale, u_star_correct_nat  )



    print >> log
    print >> log, "Anisotropic absolute scaling of derivative"
    print >> log, "------------------------------------------"
    print >> log

    aniso_derivative = absolute_scaling.ml_aniso_absolute_scaling(
      miller_array = miller_array_derivative,
      n_residues = n_residues*\
        miller_array_derivative.space_group().order_z()*n_copies_solc,
      n_bases = n_bases*\
        miller_array_derivative.space_group().order_z()*n_copies_solc)
    aniso_derivative.show(out=log,verbose=1)
    print >> log, "removing anisotropy for derivative data"
    u_star_correct_der = aniso_derivative.u_star
    miller_array_derivative = absolute_scaling.anisotropic_correction(
      miller_array_derivative,aniso_native.p_scale, u_star_correct_der  )

    print >> log
    print >> log, "Checking for alternative indexing"
    print >> log, "---------------------------------"
    print >> log
    reindex_object = pair_analyses.reindexing(
      miller_array_derivative, miller_array_native)
    #miller_array_derivative = reindex_object.select_and_transform()


    ## The reflections have now been scaled and transformed properly
    ## Get the common sets and prepare to scale according to methods
    ## specified

    ls_scaling = relative_scaling.ls_rel_scale_driver(
      miller_array_native,
      miller_array_derivative)
    ls_scaling.show()











if (__name__ == "__main__"):
  run(sys.argv[1:])
