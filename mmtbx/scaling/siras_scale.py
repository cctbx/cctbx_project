from __future__ import absolute_import, division, print_function
from cctbx import crystal
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import mmtbx.scaling
from mmtbx.scaling import pair_analyses
from libtbx.str_utils import StringIO
from mmtbx.scaling import pre_scale, make_param
import sys, os


from mmtbx.scaling import fa_estimation


params_generator = make_param.phil_lego()
master_params = iotbx.phil.parse( params_generator.default_siras() )

def run(args):

  if len(args)==0:
    master_params.show(expert_level=0)
  elif ( "--help" in args ):
    print("no help available as yet")
  elif ( "--h" in args ):
    print("no help availableas yet")
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
    print("#phil __OFF__", file=log)
    print(file=log)
    print(date_and_time(), file=log)
    print(file=log)
    print(file=log)

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
        print("##----------------------------------------------##", file=log)
        print("## Unknown phil-file or phil-command:", arg, file=log)
        print("##----------------------------------------------##", file=log)
        print(file=log)
        raise Sorry("Unknown file format or phil command: %s" % arg)


    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()


    ## Now please read in the reflections files

    ## get symmetry and cell data first please
    ## By default, the native cell and symmetry are used
    ## as reference
    crystal_symmetry_nat = None
    print(params.scaling.input.xray_data.native.file_name)
    crystal_symmetry_nat = crystal_symmetry_from_any.extract_from(
      file_name=params.scaling.input.xray_data.native.file_name)

    if params.scaling.input.xray_data.space_group is None:
      params.scaling.input.xray_data.space_group =\
        crystal_symmetry_nat.space_group_info()
      print("Using symmetry of native data", file=log)

    if params.scaling.input.xray_data.unit_cell is None:
      params.scaling.input.xray_data.unit_cell =\
        crystal_symmetry_nat.unit_cell()
      print("Using cell of native data", file=log)

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
    print("Effective parameters", file=log)
    print("#phil __ON__", file=log)
    new_params.show(out=log,expert_level=effective_params.expert_level)
    print("#phil __END__", file=log)
    print(file=log)

    ## define a xray data server
    xray_data_server =  reflection_file_utils.reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      force_symmetry = True,
      reflection_files=[])

    ## Read in native data and make appropriatre selections
    miller_array_native = None
    miller_array_native = xray_data_server.get_xray_data(
      file_name = params.scaling.input.xray_data.native.file_name,
      labels = params.scaling.input.xray_data.native.labels,
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



    ## Read in derivative data and make appropriate selections
    miller_array_derivative = None
    miller_array_derivative = xray_data_server.get_xray_data(
      file_name = params.scaling.input.xray_data.derivative.file_name,
      labels = params.scaling.input.xray_data.derivative.labels,
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

    ## Make sure we have anomalous diffs
    assert miller_array_derivative.anomalous_flag()

    ## As this is a SIR case, we will remove any anomalous pairs from the native
    if miller_array_native.anomalous_flag():
      miller_array_native = miller_array_native.average_bijvoet_mates()\
      .set_observation_type( miller_array_native )

    ## Use this copy for anomalous diff's later
    miller_array_derivative_anom = miller_array_derivative.deep_copy()

    if miller_array_derivative.anomalous_flag():
      miller_array_derivative = miller_array_derivative.average_bijvoet_mates()\
      .set_observation_type( miller_array_derivative )

    ## Print info
    print(file=log)
    print("Native data", file=log)
    print("===========", file=log)
    miller_array_native.show_comprehensive_summary(f=log)
    print(file=log)
    native_pre_scale = pre_scale.pre_scaler(
      miller_array_native,
      params.scaling.input.scaling_strategy.pre_scaler_protocol,
      params.scaling.input.basic)
    miller_array_native =  native_pre_scale.x1.deep_copy()
    del native_pre_scale

    print(file=log)
    print("Derivative data (merged friedels)", file=log)
    print("=================================", file=log)
    miller_array_derivative.show_comprehensive_summary(f=log)
    print(file=log)
    derivative_pre_scale = pre_scale.pre_scaler(
      miller_array_derivative,
      params.scaling.input.scaling_strategy.pre_scaler_protocol,
      params.scaling.input.basic)
    miller_array_derivative =  derivative_pre_scale.x1.deep_copy()
    del derivative_pre_scale


    print(file=log)
    print("Anomalous data (non merged Friedels of derivative)", file=log)
    print("==================================================", file=log)
    miller_array_derivative_anom.show_comprehensive_summary(f=log)
    print(file=log)
    derivative_anom_pre_scale = pre_scale.pre_scaler(
      miller_array_derivative_anom,
      params.scaling.input.scaling_strategy.pre_scaler_protocol,
      params.scaling.input.basic)
    miller_array_derivative_anom =  derivative_anom_pre_scale.x1.deep_copy()


    print(file=log)
    print("Working on isomorphous differences", file=log)
    print("==================================", file=log)
    print(file=log)
    iso_scaler = fa_estimation.combined_scaling(
      miller_array_native,
      miller_array_derivative,
      params.scaling.input.scaling_strategy.iso_protocol)

    miller_array_native = iso_scaler.x1.deep_copy()
    miller_array_derivative = iso_scaler.x2.deep_copy()
    del iso_scaler

    delta_gen_iso = pair_analyses.delta_generator(
      miller_array_native,
      miller_array_derivative )

    print(file=log)
    print("Working on anomalous differences", file=log)
    print("================================", file=log)
    print(file=log)

    ano_scaler = fa_estimation.ano_scaling(
      miller_array_derivative_anom,
      params.scaling.input.scaling_strategy.ano_protocol)

    positive_miller = ano_scaler.x1p.deep_copy()
    negative_miller = ano_scaler.x1n.deep_copy()
    del ano_scaler
    delta_gen_ano = pair_analyses.delta_generator(
      positive_miller,
      negative_miller)

    print(file=log)
    print("Combining iso and ano data", file=log)
    print("==========================", file=log)
    print(file=log)

    fa = fa_estimation.naive_fa_estimation(
      delta_gen_ano.abs_delta_f,
      delta_gen_iso.abs_delta_f,
      params.scaling.input.fa_estimation
    )

    print(file=log)
    print("writing mtz file", file=log)
    print("----------------", file=log)
    print(file=log)

    ## Please write out the abs_delta_f array
    mtz_dataset = fa.fa.as_mtz_dataset(
      column_root_label='F'+params.scaling.input.output.outlabel)
    mtz_dataset.mtz_object().write(
      file_name=params.scaling.input.output.hklout)


if (__name__ == "__main__"):
  run(sys.argv[1:])
