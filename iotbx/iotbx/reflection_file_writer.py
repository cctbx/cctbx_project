import iotbx.mtz
import iotbx.cns.miller_array
import iotbx.scalepack.merge
import iotbx.shelx.hklf
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx.option_parser import iotbx_option_parser

def run(args):
  command_line = (iotbx_option_parser(
    usage="iotbx.reflection_file_writer [options] reflection_file ...",
    description="Example: iotbx.reflection_file_writer w1.sca --mtz .")
    .enable_symmetry_comprehensive()
    .option(None, "--weak_symmetry",
      action="store_true",
      default=00000,
      dest="weak_symmetry",
      help="symmetry on command line is weaker than symmetry found in files")
    .option(None, "--label",
      action="store",
      type="string",
      dest="label",
      help="Substring of reflection data label or number",
      metavar="STRING")
    .option(None, "--resolution",
      action="store",
      type="float",
      dest="resolution",
      help="High resolution limit",
      metavar="FLOAT")
    .option(None, "--low_resolution",
      action="store",
      type="float",
      dest="low_resolution",
      help="Low resolution limit",
      metavar="FLOAT")
    .option(None, "--observation_type",
      choices=("amplitude", "intensity"),
      metavar="amplitude|intensity")
    .option(None, "--scale_max",
      action="store",
      type="float",
      dest="scale_max",
      help="Scales data such that the maximum is equal to the given value",
      metavar="FLOAT")
    .option(None, "--scale_factor",
      action="store",
      type="float",
      dest="scale_factor",
      help="Multiplies data with the given factor",
      metavar="FLOAT")
    .option(None, "--sca",
      action="store",
      type="string",
      dest="sca",
      help=
        "write data to Scalepack FILE ('--sca .' copies name of input file)",
      metavar="FILE")
    .option(None, "--mtz",
      action="store",
      type="string",
      dest="mtz",
      help="write data to MTZ FILE ('--mtz .' copies name of input file)",
      metavar="FILE")
    .option(None, "--cns",
      action="store",
      type="string",
      dest="cns",
      help="write data to CNS FILE ('--cns .' copies name of input file)",
      metavar="FILE")
    .option(None, "--shelx",
      action="store",
      type="string",
      dest="shelx",
      help="write data to SHELX FILE ('--shelx .' copies name of input file)",
      metavar="FILE")
  ).process()
  if (    command_line.options.scale_max is not None
      and command_line.options.scale_factor is not None):
    print
    print "--scale_max and --scale_factor are mutually exclusive."
    print
    return None
  if (len(command_line.args) == 0):
    command_line.parser.show_help()
    return None
  all_miller_arrays = reflection_file_reader.collect_arrays(
    file_names=command_line.args,
    crystal_symmetry=command_line.symmetry,
    force_symmetry=not command_line.options.weak_symmetry,
    discard_arrays=00000,
    verbose=1)
  label_table = reflection_file_utils.label_table(
    miller_arrays=all_miller_arrays)
  if (len(all_miller_arrays) == 1):
    selected_array = all_miller_arrays[0]
  elif (command_line.options.label is None):
    print
    print "Please use --label to select a reflection array."
    print "For example: --label=2"
    print "         or: --label=%s"%all_miller_arrays[1].info().split(":")[-1]
    print
    label_table.show_possible_choices()
    return None
  else:
    selected_array = label_table.match_data_label(
      label=command_line.options.label)
    if (selected_array is None):
      return None
  if (selected_array.unit_cell() is None):
    command_line.parser.show_help()
    print "Unit cell parameters unknown. Please use --symmetry or --unit_cell."
    print
    return None
  if (selected_array.space_group_info() is None):
    command_line.parser.show_help()
    print "Space group unknown. Please use --symmetry or --space_group."
    print
    return None
  if (selected_array.observation_type() is None):
    if (command_line.options.observation_type is None):
      command_line.parser.show_help()
      print "Observation type is unknown. Please use --observation_type."
      print
      return None
    if (command_line.options.observation_type == "amplitude"):
      selected_array.set_observation_type_xray_amplitude()
    else:
      selected_array.set_observation_type_xray_intensity()
  print "Selected data:"
  print " ", selected_array.info()
  print "  Observation type:", selected_array.observation_type()
  print
  d_max = command_line.options.low_resolution
  d_min = command_line.options.resolution
  if (d_max is not None or d_min is not None):
    if (d_max is not None):
      print "Applying low resolution cutoff: d_max=%.6g" % d_max
    if (d_min is not None):
      print "Applying high resolution cutoff: d_min=%.6g" % d_min
    selected_array = selected_array.resolution_filter(d_max=d_max, d_min=d_min)
    print "Number of reflections:", selected_array.indices().size()
    print
  if (command_line.options.scale_max is not None):
    print "Scaling data such that the maximum value is: %.6g" \
      % command_line.options.scale_max
    selected_array = selected_array.apply_scaling(
      target_max=command_line.options.scale_max)
    print
  if (command_line.options.scale_factor is not None):
    print "Multiplying data with the factor: %.6g" \
      % command_line.options.scale_factor
    selected_array = selected_array.apply_scaling(
      factor=command_line.options.scale_factor)
    print
  n_output_files = 0
  if (command_line.options.sca is not None):
    file_name = reflection_file_utils.construct_output_file_name(
      input_file_names=command_line.args,
      user_file_name=command_line.options.sca,
      file_type_label="Scalepack",
      file_extension="sca")
    print "Writing Scalepack file:", file_name
    iotbx.scalepack.merge.write(
      file_name=file_name,
      miller_array=selected_array)
    n_output_files += 1
    print
  if (command_line.options.mtz is not None):
    file_name = reflection_file_utils.construct_output_file_name(
      input_file_names=command_line.args,
      user_file_name=command_line.options.mtz,
      file_type_label="MTZ",
      file_extension="mtz")
    print "Writing MTZ file:", file_name
    selected_array.export_as_mtz(file_name, file_name[:-4])
    n_output_files += 1
    print
  if (command_line.options.cns is not None):
    file_name = reflection_file_utils.construct_output_file_name(
      input_file_names=command_line.args,
      user_file_name=command_line.options.cns,
      file_type_label="CNS",
      file_extension="cns")
    print "Writing CNS file:", file_name
    selected_array.export_as_cns_hkl(open(file_name, "w"), file_name)
    n_output_files += 1
    print
  if (command_line.options.shelx is not None):
    file_name = reflection_file_utils.construct_output_file_name(
      input_file_names=command_line.args,
      user_file_name=command_line.options.shelx,
      file_type_label="SHELX",
      file_extension="shelx")
    print "Writing SHELX file:", file_name
    selected_array.as_amplitude_array().export_as_shelx_hklf(
      open(file_name, "w"))
    n_output_files += 1
    print
  if (n_output_files == 0):
    command_line.parser.show_help()
    print "Please specify at least one output file format,",
    print "e.g. --mtz, --sca, etc."
    print
    return None
  return selected_array
