from cctbx import crystal
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import mmtbx.scaling
from mmtbx.scaling import pair_analyses
from libtbx.str_utils import StringIO
from mmtbx.scaling import pre_scale
from mmtbx.scaling import make_param
import sys, os

from mmtbx.scaling import fa_estimation
from mmtbx.scaling import random_omit



params_generator = make_param.phil_lego()
master_params = iotbx.phil.parse( params_generator.default_sad() )


def run(args):

  if len(args)==0:
    master_params.show(expert_level=0)
  elif ( "--help" in args ):
    print "no help available as yet"
  elif ( "--h" in args ):
    print "no help availableas yet"
  elif ( "--show_defaults" in args ):
    master_params.show(expert_level=0)
  elif ( "--show_defaults_all" in args ):
    master_params.show(expert_level=10)
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
    argument_interpreter = master_params.command_line_argument_interpreter(
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
        except Exception : pass
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
        except Exception : pass

      if not arg_is_processed:
        print >> log, "##----------------------------------------------##"
        print >> log, "## Unknown phil-file or phil-command:", arg
        print >> log, "##----------------------------------------------##"
        print >> log
        raise Sorry("Unknown file format or phil command: %s" % arg)


    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()

    effective_params = master_params.fetch(sources=phil_objects)
    new_params = master_params.format(python_object=params)

    ## Now please read in the reflections files

    ## get symmetry and cell data first please
    ## By default, the native cell and symmetry are used
    ## as reference
    crystal_symmetry_nat = None
    crystal_symmetry_nat = crystal_symmetry_from_any.extract_from(
      file_name=params.scaling.input.xray_data.reference.file_name)

    if params.scaling.input.xray_data.space_group is None:
      params.scaling.input.xray_data.space_group =\
        crystal_symmetry_nat.space_group_info()
      print >> log, "Using symmetry of native data"

    if params.scaling.input.xray_data.unit_cell is None:
      params.scaling.input.xray_data.unit_cell =\
        crystal_symmetry_nat.unit_cell()
      print >> log, "Using cell of native data"

    ## Check if a unit cell is defined
    if params.scaling.input.xray_data.space_group is None:
      raise Sorry("No space group defined")
    if params.scaling.input.xray_data.unit_cell is None:
      raise Sorry("No unit cell defined")


    crystal_symmetry = crystal_symmetry = crystal.symmetry(
      unit_cell =  params.scaling.input.xray_data.unit_cell,
      space_group_symbol = str(
        params.scaling.input.xray_data.space_group) )


    effective_params = master_params.fetch(sources=phil_objects)
    new_params = master_params.format(python_object=params)
    print >> log, "Effective parameters"
    print >> log, "#phil __ON__"
    new_params.show(out=log,
                    expert_level=params.scaling.input.expert_level)
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
      file_name = params.scaling.input.xray_data.reference.file_name,
      labels = params.scaling.input.xray_data.reference.labels,
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
    print >> log, "Reference data"
    print >> log, "=============="
    miller_array_native.show_comprehensive_summary(f=log)
    print >> log
    native_pre_scale = pre_scale.pre_scaler(
      miller_array_native,
      params.scaling.input.scaling_strategy.pre_scaler_protocol,
      params.scaling.input.basic)
    miller_array_native =  native_pre_scale.x1.deep_copy()
    del native_pre_scale

    scaler = fa_estimation.ano_scaling(
      miller_array_native,
      params.scaling.input.scaling_strategy.ano_protocol)

    positive_miller = scaler.x1p.deep_copy()
    negative_miller = scaler.x1n.deep_copy()


    print >> log
    print >> log, "Making delta f's"
    print >> log, "----------------"
    print >> log

    delta_gen = pair_analyses.delta_generator( positive_miller,
                                               negative_miller )
    print >> log
    print >> log, "writing mtz file"
    print >> log, "----------------"
    print >> log

    ## some assertions to make sure nothing went weerd
    assert positive_miller.observation_type() is not None
    assert negative_miller.observation_type() is not None
    assert delta_gen.abs_delta_f.observation_type() is not None

    ## Please write out the abs_delta_f array

    mtz_dataset = delta_gen.abs_delta_f.as_mtz_dataset(
      column_root_label='F'+params.scaling.input.output.outlabel)
    mtz_dataset.mtz_object().write(
      file_name=params.scaling.input.output.hklout)

    if params.scaling.input.omit.perform_omit:
      print >> log
      print >> log, "writing omit files"
      print >> log, "------------------"
      print >> log
      omit_object = random_omit.random_omit_data(
        delta_gen.abs_delta_f,
        params.scaling.input.omit )
      omit_object.write_datasets()







if (__name__ == "__main__"):
  run(sys.argv[1:])
