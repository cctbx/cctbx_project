import iotbx.cns.reflection_reader
from iotbx.option_parser import iotbx_option_parser
import sys, os

def run(args):
  if (len(args) == 0): args = ["--help"]
  command_line = (iotbx_option_parser(
    usage="iotbx.cns_as_mtz [options] cns_file",
    description="Example: iotbx.cns_as_mtz scale.hkl")
    .enable_symmetry_comprehensive()
    .option("-q", "--quiet",
      action="store_true",
      default=False,
      help="suppress output")
  ).process(args=args)
  if (len(command_line.args) != 1):
    command_line.parser.show_help()
    return
  cns_file_name = command_line.args[0]
  crystal_symmetry = command_line.symmetry
  if (crystal_symmetry.unit_cell() is None):
    print
    print "*" * 79
    print "Unknown unit cell parameters."
    print "Use --symmetry or --unit_cell to define unit cell:"
    print "*" * 79
    print
    command_line.parser.show_help()
    return
  if (crystal_symmetry.space_group_info() is None):
    print
    print "*" * 79
    print "Unknown space group."
    print "Use --symmetry or --space_group to define space group:"
    print "*" * 79
    print
    command_line.parser.show_help()
    return
  if (not command_line.options.quiet):
    print "CNS file name:", cns_file_name
    print "Crystal symmetry:"
    crystal_symmetry.show_summary(prefix="  ")
  reflection_file = iotbx.cns.reflection_reader.cns_reflection_file(
    file_handle=open(cns_file_name, "r"))
  if (not command_line.options.quiet):
    reflection_file.show_summary()

if (__name__ == "__main__"):
  run(sys.argv[1:])
