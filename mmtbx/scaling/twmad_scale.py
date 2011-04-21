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
master_params = iotbx.phil.parse( params_generator.default_2wmad() )


def run(args):

  if len(args)==0:
    master_params.show(expert_level=100)
  elif ( "--help" in args ):
    print "no help available"
  elif ( "--h" in args ):
    print "no help available"
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
    print params.scaling.input.xray_data.wavelength1.file_name
    crystal_symmetry_nat = crystal_symmetry_from_any.extract_from(
      file_name=params.scaling.input.xray_data.wavelength1.file_name)

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
    new_params.show(out=log,expert_level=params.scaling.input.expert_level)
    print >> log, "#phil __END__"
    print >> log

    ## define a xray data server
    xray_data_server =  reflection_file_utils.reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      force_symmetry = True,
      reflection_files=[])

    ## Read in native data and make appropriate selections
    miller_array_w1 = None
    miller_array_w1 = xray_data_server.get_xray_data(
      file_name = params.scaling.input.xray_data.wavelength1.file_name,
      labels = params.scaling.input.xray_data.wavelength1.labels,
      ignore_all_zeros = True,
      parameter_scope = 'scaling.input.SIR_scale.xray_data.native'
    )
    info_native = miller_array_w1.info()
    miller_array_w1=miller_array_w1.map_to_asu().select(
      miller_array_w1.indices()!=(0,0,0) )
    miller_array_w1 = miller_array_w1.select(
      miller_array_w1.data() > 0 )
    ## Convert to amplitudes
    if (miller_array_w1.is_xray_intensity_array()):
      miller_array_w1 = miller_array_w1.f_sq_as_f()
    elif (miller_array_w1.is_complex_array()):
      miller_array_w1 = abs(miller_array_w1)
    if not miller_array_w1.is_real_array():
      raise Sorry("miller_array_native is not a real array")
    miller_array_w1.set_info(info = info_native)



    ## Read in derivative data and make appropriate selections
    miller_array_w2 = None
    miller_array_w2 = xray_data_server.get_xray_data(
      file_name = params.scaling.input.xray_data.wavelength2.file_name,
      labels = params.scaling.input.xray_data.wavelength2.labels,
      ignore_all_zeros = True,
      parameter_scope = 'scaling.input.SIR_scale.xray_data.derivative'
    )
    info_w2 = miller_array_w2.info()
    miller_array_w2=miller_array_w2.map_to_asu().select(
      miller_array_w2.indices()!=(0,0,0) )
    miller_array_w2 = miller_array_w2.select(
      miller_array_w2.data() > 0 )
    ## Convert to amplitudes
    if (miller_array_w2.is_xray_intensity_array()):
      miller_array_w2 = miller_array_w2.f_sq_as_f()
    elif (miller_array_w2.is_complex_array()):
      miller_array_w2 = abs(miller_array_w2)
    if not miller_array_w2.is_real_array():
      raise Sorry("miller_array_derivative is not a real array")
    miller_array_w2.set_info(info = info_w2)

    ## Make sure we have anomalous diffs in both files
    assert miller_array_w1.anomalous_flag()
    assert miller_array_w2.anomalous_flag()


    ## Print info
    print >> log
    print >> log, "Wavelength 1"
    print >> log, "============"
    miller_array_w1.show_comprehensive_summary(f=log)
    print >> log
    w1_pre_scale = pre_scale.pre_scaler(
      miller_array_w1,
      params.scaling.input.scaling_strategy.pre_scaler_protocol,
      params.scaling.input.basic)
    miller_array_w1 =  w1_pre_scale.x1.deep_copy()
    del w1_pre_scale

    print >> log
    print >> log, "Wavelength 2"
    print >> log, "============"
    miller_array_w2.show_comprehensive_summary(f=log)
    print >> log
    w2_pre_scale = pre_scale.pre_scaler(
      miller_array_w2,
      params.scaling.input.scaling_strategy.pre_scaler_protocol,
      params.scaling.input.basic)
    miller_array_w2 =  w2_pre_scale.x1.deep_copy()
    del w2_pre_scale


    print >> log
    print >> log, "Checking for possible reindexing schemes"
    print >> log, "----------------------------------------"
    print >> log
    print >> log, "Reindexing operator derived as described in:"
    print >> log, "Grosse-Kunstleve, Afonine, Sauter & Adams. (2005)."
    print >> log, "  IUCr Computing Commission Newsletter 5."
    print >> log

    reindex_object = pair_analyses.reindexing(
       set_a=miller_array_w1,
       set_b=miller_array_w2,
       out=log)
    miller_array_w2 = reindex_object.select_and_transform()
    miller_array_w2.map_to_asu()

    print >> log
    print >> log, "Relative scaling of 2-wavelength mad data"
    print >> log, "-----------------------------------------"
    print >> log
    scaler = fa_estimation.combined_scaling(
      miller_array_w1,
      miller_array_w2,
      params.scaling.input.scaling_strategy.iso_protocol)

    miller_array_w1 = scaler.x1.deep_copy()
    miller_array_w2 = scaler.x2.deep_copy()

    del scaler

    print >> log
    print >> log, "Estimating f\" and f' ratios"
    print >> log, "----------------------------"
    print >> log



    # now things are scaled see if we can guestimate the ratio
    fdpratio = pair_analyses.f_double_prime_ratio(
      miller_array_w1,
      miller_array_w2)

    fpfdpratio = pair_analyses.delta_f_prime_f_double_prime_ratio(
      miller_array_w1,
      miller_array_w2)

    k1 = fdpratio.ratio
    k2 = fpfdpratio.ratio

    if k1 is not None:
      print >> log
      print >> log, "  The estimate of f\"(w1)/f\"(w2) is %3.2f"\
            %(fdpratio.ratio)
    if k2 is not None:
      print >> log, "  The estimate of (f'(w1)-f'(w2))/f\"(w2) is %3.2f"\
            %(fpfdpratio.ratio)
      print >> log
      print >> log, "  The quality of these estimates depends to a large extend"
      print >> log, "  on the quality of the data. If user supplied values"
      print >> log, "  of f\" and f' are given, they will be used instead "
      print >> log, "  of the estimates."
      print >> log

    if params.scaling.input.xray_data.wavelength1.f_double_prime is not None:
      if params.scaling.input.xray_data.wavelength2.f_double_prime is not None:
        k1 = (params.scaling.input.xray_data.wavelength1.f_double_prime/
              params.scaling.input.xray_data.wavelength2.f_double_prime)
        print >> log, "    Using user specified f\" values"
        print >> log, "      user specified f\"(w1)/f\"(w2) is %3.2f"\
              %(k1)
        print >> log
    if params.scaling.input.xray_data.wavelength1.f_prime is not None:
      if params.scaling.input.xray_data.wavelength2.f_prime is not None:
        if params.scaling.input.xray_data.wavelength2.f_double_prime is not None:

          k2 = (params.scaling.input.xray_data.wavelength1.f_prime-
                params.scaling.input.xray_data.wavelength2.f_prime)\
                /params.scaling.input.xray_data.wavelength2.f_double_prime
          print >> log, "    Using user specified f\" and f' values"
          print >> log, "     user specified f\"(w1)/f\"(w2) is %3.2f"\
                %(k2)
          print >> log



    fa_gen = fa_estimation.twmad_fa_driver(miller_array_w1,
                                           miller_array_w2,
                                           k1,
                                           k2,
                                           params.scaling.input.fa_estimation)

    print >> log
    print >> log, "writing mtz file"
    print >> log, "----------------"
    print >> log

    ## Please write out the abs_delta_f array

    fa =  fa_gen.fa_values

    mtz_dataset = fa.as_mtz_dataset(
      column_root_label='F'+params.scaling.input.output.outlabel)

    mtz_dataset.mtz_object().write(
      file_name=params.scaling.input.output.hklout)




if (__name__ == "__main__"):
  run(sys.argv[1:])
