import iotbx.mtz
import iotbx.cns.miller_array
import iotbx.scalepack.merge
import iotbx.shelx.hklf
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx.option_parser import iotbx_option_parser
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx.python_utils.misc import plural_s
from libtbx.utils import date_and_time
import os

def run(args, simply_return_all_miller_arrays=False):
  command_line = (iotbx_option_parser(
    usage="iotbx.reflection_file_converter [options] reflection_file ...",
    description="Example: iotbx.reflection_file_converter w1.sca --mtz .")
    .enable_symmetry_comprehensive()
    .option(None, "--weak_symmetry",
      action="store_true",
      default=False,
      dest="weak_symmetry",
      help="symmetry on command line is weaker than symmetry found in files")
    .enable_resolutions()
    .option(None, "--label",
      action="store",
      type="string",
      dest="label",
      help="Substring of reflection data label or number",
      metavar="STRING")
    .option(None, "--non_anomalous",
      action="store_true",
      default=False,
      dest="non_anomalous",
      help="Averages Bijvoet mates to obtain a non-anomalous array")
    .option(None, "--change_of_basis",
      action="store",
      type="string",
      dest="change_of_basis",
      help="Change-of-basis operator: h,k,l or x,y,z"
           " or to_reference_setting, to_primitive_setting, to_niggli_cell,"
           " to_inverse_hand",
      metavar="STRING")
    .option(None, "--expand_to_p1",
      action="store_true",
      default=False,
      dest="expand_to_p1",
      help="Generates all symmetrically equivalent reflections."
           " The space group symmetry is reset to P1."
           " May be used in combination with --change_to_space_group to"
           " lower the symmetry.")
    .option(None, "--change_to_space_group",
      action="store",
      type="string",
      dest="change_to_space_group",
      help="Changes the space group and merges equivalent reflections"
           " if necessary",
      metavar="SYMBOL|NUMBER")
    .option(None, "--write_mtz_amplitudes",
      action="store_true",
      default=False,
      dest="write_mtz_amplitudes",
      help="Converts intensities to amplitudes before writing MTZ format;"
           " requires --mtz_root_label")
    .option(None, "--write_mtz_intensities",
      action="store_true",
      default=False,
      dest="write_mtz_intensities",
      help="Converts amplitudes to intensities before writing MTZ format;"
           " requires --mtz_root_label")
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
    .option(None, "--mtz_root_label",
      action="store",
      type="string",
      dest="mtz_root_label",
      help="Root label for MTZ file (e.g. Fobs)",
      metavar="STRING")
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
  ).process(args=args)
  if (    command_line.options.write_mtz_amplitudes
      and command_line.options.write_mtz_intensities):
    print
    print "--write_mtz_amplitudes and --write_mtz_intensities" \
          " are mutually exclusive."
    print
    return None
  if (   command_line.options.write_mtz_amplitudes
      or command_line.options.write_mtz_intensities):
    if (command_line.options.mtz_root_label is None):
      print
      print "--write_mtz_amplitudes and --write_mtz_intensities" \
            " require --mtz_root_label."
      print
      return None
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
    crystal_symmetry=None,
    force_symmetry=False,
    merge_equivalents=False,
    discard_arrays=False,
    verbose=1)
  if (simply_return_all_miller_arrays):
    return all_miller_arrays
  if (len(all_miller_arrays) == 0):
    print
    print "No reflection data found in input file%s." % (
      plural_s(len(command_line.args))[1])
    print
    return None
  label_table = reflection_file_utils.label_table(
    miller_arrays=all_miller_arrays)
  if (len(all_miller_arrays) == 1):
    selected_array = all_miller_arrays[0]
  elif (command_line.options.label is None):
    print
    print "Please use --label to select a reflection array."
    print "For example: --label=%s" % str(
      all_miller_arrays[1].info()).split(":")[-1]
    print
    label_table.show_possible_choices()
    return None
  else:
    selected_array = label_table.match_data_label(
      label=command_line.options.label,
      command_line_switch="--label")
    if (selected_array is None):
      return None
  print "Selected data:"
  print " ", selected_array.info()
  print "  Observation type:", selected_array.observation_type()
  print
  processed_array = selected_array.customized_copy(
    crystal_symmetry=selected_array.join_symmetry(
      other_symmetry=command_line.symmetry,
      force=not command_line.options.weak_symmetry)).set_observation_type(
        selected_array.observation_type())
  print "Input crystal symmetry:"
  crystal.symmetry.show_summary(processed_array, prefix="  ")
  print
  if (processed_array.unit_cell() is None):
    command_line.parser.show_help()
    print "Unit cell parameters unknown. Please use --symmetry or --unit_cell."
    print
    return None
  if (processed_array.space_group_info() is None):
    command_line.parser.show_help()
    print "Space group unknown. Please use --symmetry or --space_group."
    print
    return None
  if (command_line.options.change_of_basis is not None):
    print "Change of basis:"
    if   (command_line.options.change_of_basis == "to_reference_setting"):
      cb_op = processed_array.change_of_basis_op_to_reference_setting()
    elif (command_line.options.change_of_basis == "to_primitive_setting"):
      cb_op = processed_array.change_of_basis_op_to_primitive_setting()
    elif (command_line.options.change_of_basis == "to_niggli_cell"):
      cb_op = processed_array.change_of_basis_op_to_niggli_cell()
    elif (command_line.options.change_of_basis == "to_inverse_hand"):
      cb_op = processed_array.change_of_basis_op_to_inverse_hand()
    else:
      cb_op = sgtbx.change_of_basis_op(command_line.options.change_of_basis)
    if (cb_op.c_inv().t().is_zero()):
      print "  Change of basis operator in both h,k,l and x,y,z notation:"
      print "   ", cb_op.as_hkl()
    else:
      print "  Change of basis operator in x,y,z notation:"
    print "    %s [Inverse: %s]" % (cb_op.as_xyz(), cb_op.inverse().as_xyz())
    d = cb_op.c().r().determinant()
    print "  Determinant:", d
    if (d < 0 and command_line.options.change_of_basis != "to_inverse_hand"):
      print "  **************************************************************"
      print "  W A R N I N G: This change of basis operator changes the hand!"
      print "  **************************************************************"
    processed_array = processed_array.change_basis(cb_op=cb_op)
    print "  Crystal symmetry after change of basis:"
    crystal.symmetry.show_summary(processed_array, prefix="    ")
    print
  if (not processed_array.is_unique_set_under_symmetry()):
    print "Merging symmetry-equivalent values:"
    merged = processed_array.merge_equivalents()
    merged.show_summary(prefix="  ")
    print
    processed_array = merged.array()
    del merged
    processed_array.show_comprehensive_summary(prefix="  ")
    print
  if (command_line.options.expand_to_p1):
    print "Expanding symmetry and resetting space group to P1:"
    processed_array = processed_array.expand_to_p1()
    processed_array.show_comprehensive_summary(prefix="  ")
    print
  if (command_line.options.change_to_space_group is not None):
    new_space_group_info = sgtbx.space_group_info(
      symbol=command_line.options.change_to_space_group)
    print "Change to space group:", new_space_group_info
    new_crystal_symmetry = crystal.symmetry(
      unit_cell=processed_array.unit_cell(),
      space_group_info=new_space_group_info,
      assert_is_compatible_unit_cell=False)
    if (not new_crystal_symmetry.unit_cell()
              .is_similar_to(processed_array.unit_cell())):
      print "  *************"
      print "  W A R N I N G"
      print "  *************"
      print "  Unit cell parameters adapted to new space group symmetry are"
      print "  significantly different from input unit cell parameters:"
      print "      Input unit cell parameters:", \
        processed_array.unit_cell()
      print "    Adapted unit cell parameters:", \
        new_crystal_symmetry.unit_cell()
    processed_array = processed_array.customized_copy(
      crystal_symmetry=new_crystal_symmetry)
    print
    if (not processed_array.is_unique_set_under_symmetry()):
      print "  Merging values symmetry-equivalent under new symmetry:"
      merged = processed_array.merge_equivalents()
      merged.show_summary(prefix="    ")
      print
      processed_array = merged.array()
      del merged
      processed_array.show_comprehensive_summary(prefix="    ")
      print
  if (processed_array.anomalous_flag()
      and command_line.options.non_anomalous):
    if (not processed_array.is_xray_intensity_array()):
      processed_array = processed_array.average_bijvoet_mates()
    else:
      processed_array = processed_array.f_sq_as_f()
      processed_array = processed_array.average_bijvoet_mates()
      processed_array = processed_array.f_as_f_sq()
    processed_array.set_observation_type_xray_intensity()
  d_max = command_line.options.low_resolution
  d_min = command_line.options.resolution
  if (d_max is not None or d_min is not None):
    if (d_max is not None):
      print "Applying low resolution cutoff: d_max=%.6g" % d_max
    if (d_min is not None):
      print "Applying high resolution cutoff: d_min=%.6g" % d_min
    processed_array = processed_array.resolution_filter(
      d_max=d_max, d_min=d_min)
    print "Number of reflections:", processed_array.indices().size()
    print
  if (command_line.options.scale_max is not None):
    print "Scaling data such that the maximum value is: %.6g" \
      % command_line.options.scale_max
    processed_array = processed_array.apply_scaling(
      target_max=command_line.options.scale_max)
    print
  if (command_line.options.scale_factor is not None):
    print "Multiplying data with the factor: %.6g" \
      % command_line.options.scale_factor
    processed_array = processed_array.apply_scaling(
      factor=command_line.options.scale_factor)
    print
  n_output_files = 0
  if (command_line.options.sca is not None):
    file_name = reflection_file_utils.construct_output_file_name(
      input_file_names=[selected_array.info().source],
      user_file_name=command_line.options.sca,
      file_type_label="Scalepack",
      file_extension="sca")
    print "Writing Scalepack file:", file_name
    iotbx.scalepack.merge.write(
      file_name=file_name,
      miller_array=processed_array)
    n_output_files += 1
    print
  if (command_line.options.mtz is not None):
    file_name = reflection_file_utils.construct_output_file_name(
      input_file_names=[selected_array.info().source],
      user_file_name=command_line.options.mtz,
      file_type_label="MTZ",
      file_extension="mtz")
    print "Writing MTZ file:", file_name
    mtz_history_buffer = flex.std_string()
    mtz_history_buffer.append(date_and_time())
    mtz_history_buffer.append("> program: iotbx.reflection_file_converter")
    mtz_history_buffer.append("> input file name: %s" %
      os.path.basename(selected_array.info().source))
    mtz_history_buffer.append("> input directory: %s" %
      os.path.dirname(os.path.abspath(selected_array.info().source)))
    mtz_history_buffer.append("> input labels: %s" %
      selected_array.info().label_string())
    mtz_output_array = processed_array
    if (command_line.options.write_mtz_amplitudes):
      if (not mtz_output_array.is_xray_amplitude_array()):
        print "  Converting intensities to amplitudes."
        mtz_output_array = mtz_output_array.f_sq_as_f()
        mtz_history_buffer.append("> Intensities converted to amplitudes.")
    elif (command_line.options.write_mtz_intensities):
      if (not mtz_output_array.is_xray_intensity_array()):
        print "  Converting amplitudes to intensities."
        mtz_output_array = mtz_output_array.f_as_f_sq()
        mtz_history_buffer.append("> Amplitudes converted to intensities.")
    column_root_label = command_line.options.mtz_root_label
    if (column_root_label is None):
      column_root_label = file_name[:min(24,len(file_name)-4)]
    mtz_object = mtz_output_array.as_mtz_dataset(
      column_root_label=column_root_label).mtz_object()
    del mtz_output_array
    mtz_history_buffer.append("> output file name: %s" %
      os.path.basename(file_name))
    mtz_history_buffer.append("> output directory: %s" %
      os.path.dirname(os.path.abspath(file_name)))
    mtz_object.add_history(mtz_history_buffer)
    mtz_object.write(file_name=file_name)
    n_output_files += 1
    print
  if (command_line.options.cns is not None):
    file_name = reflection_file_utils.construct_output_file_name(
      input_file_names=[selected_array.info().source],
      user_file_name=command_line.options.cns,
      file_type_label="CNS",
      file_extension="cns")
    print "Writing CNS file:", file_name
    processed_array.export_as_cns_hkl(open(file_name, "w"), file_name)
    n_output_files += 1
    print
  if (command_line.options.shelx is not None):
    file_name = reflection_file_utils.construct_output_file_name(
      input_file_names=[selected_array.info().source],
      user_file_name=command_line.options.shelx,
      file_type_label="SHELX",
      file_extension="shelx")
    print "Writing SHELX file:", file_name
    processed_array.as_amplitude_array().export_as_shelx_hklf(
      open(file_name, "w"))
    n_output_files += 1
    print
  if (n_output_files == 0):
    command_line.parser.show_help()
    print "Please specify at least one output file format,",
    print "e.g. --mtz, --sca, etc."
    print
    return None
  return processed_array
