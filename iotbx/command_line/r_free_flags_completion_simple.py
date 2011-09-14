from __future__ import division
import iotbx.reflection_file_utils
import iotbx.reflection_file_reader
from libtbx.utils import Sorry, Usage
from libtbx.str_utils import show_string
import libtbx.load_env
import os.path as op
import sys

def usage():
    raise Usage("""\
%s reflection_file_name high_resolution_value

  Example: %s your.mtz 2.0
""" % (libtbx.env.dispatcher_name, libtbx.env.dispatcher_name))

def run(args):
  if (len(args) == 0 or "--help" in args or "-h" in args):
    usage()
  refl_file_name = None
  refl_file = None
  high_res = None
  for arg in args:
    if (op.isfile(arg)):
      refl_file_name = arg
      refl_file = iotbx.reflection_file_reader.any_reflection_file(
        file_name=refl_file_name)
      if (refl_file.file_type() is None):
        raise Sorry("File not recognized: %s" % show_string(arg))
    else:
      try:
        high_res = float(arg)
      except ValueError:
        raise Sorry(
          "Not a file name or high-resolution value: %s"
            % show_string(arg))
      if (high_res <= 0):
        raise Sorry(
          "High-resolution value must be greater than zero (%s given)."
            % show_string(arg))
  if (refl_file is None or high_res is None):
    usage()
  refl_file_server = iotbx.reflection_file_utils.reflection_file_server(
    reflection_files=[refl_file])
  r_free_flags, test_flag_value = refl_file_server.get_r_free_flags(
    file_name=None,
    label=None,
    test_flag_value=None,
    disable_suitability_test=False,
    parameter_scope="r_free_flags")
  print "Summary of existing R-free-flags:"
  r_free_flags.show_comprehensive_summary(prefix="  ")
  print
  print "Test flag value:", test_flag_value
  print
  assert r_free_flags.unit_cell() is not None
  assert r_free_flags.space_group_info() is not None
  assert r_free_flags.data().size() != 0
  r_free_flags = r_free_flags.array(data=r_free_flags.data()==test_flag_value)
  fraction_free = r_free_flags.data().count(True) / r_free_flags.data().size()
  print "Fraction free: %.2f %%" % (fraction_free*100)
  assert fraction_free > 0
  print
  missing_set = r_free_flags.complete_set(d_min=high_res).lone_set(
    r_free_flags.map_to_asu())
  print "Number of missing R-free-flags:", missing_set.indices().size()
  print
  missing_flags = missing_set.generate_r_free_flags(
    fraction=fraction_free,
    max_free=None,
    use_lattice_symmetry=True)
  extended_r_free_flags = r_free_flags.concatenate(other=missing_flags)
  print "Summary of extended R-free-flags:"
  extended_r_free_flags.show_comprehensive_summary(prefix="  ")
  print
  mtz_dataset = extended_r_free_flags.as_mtz_dataset(
    column_root_label="Extended-R-free-flags")
  output_file_name = op.basename(refl_file_name)
  i = output_file_name.rfind(".")
  if (i >= 0):
    output_file_name = output_file_name[:i] + "_" + output_file_name[i+1:]
  output_file_name += "_extended_r_free_flags.mtz"
  print "Writing file: %s" % show_string(output_file_name)
  mtz_dataset.mtz_object().write(file_name=output_file_name)
  print
  print """\
Miscellaneous remarks:

  - ***PLEASE INSPECT*** the extended R-free-flags with this command:

      iotbx.r_free_flags_accumulation %s

  - For use in phenix.refine, simply add the new mtz file as an
    additional command-line argument. If necessary, follow the
    phenix.refine suggestions to select the extended R-free-flag
    array.

  - You may have to remove the
      REMARK r_free_flags.md5.hexdigest
    line from your PDB file(s) to use the new flags in phenix.refine.
    (phenix.refine will tell you why and what you need to do.)

  - phenix.reflection_file_converter can be used to combine reflection
    data from another file with the extended R-free-flags.
""" % show_string(output_file_name)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
