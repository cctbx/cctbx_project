import libtbx.phil
from libtbx.optparse_wrapper import option_parser
import sys

def run(args, command_name="libtbx.phil", converter_registry=None):
  if (len(args) == 0): args = ["--help"]
  command_line = (option_parser(
    usage="%s [options] parameter_file ..." % command_name)
    .option(None, "--show_help",
      action="store_true",
      dest="show_help",
      help="Display help for each parameter if available.")
    .option(None, "--show_some_attributes",
      action="store_true",
      dest="show_some_attributes",
      help="Display non-default attributes for each parameter.")
    .option(None, "--show_all_attributes",
      action="store_true",
      dest="show_all_attributes",
      help="Display all attributes for each parameter.")
    .option(None, "--process_includes",
      action="store_true",
      dest="process_includes",
      help="Inline include files.")
    .option(None, "--print_width",
      action="store",
      type="int",
      dest="print_width",
      help="Width for output",
      metavar="INT")
    .option(None, "--print_prefix",
      action="store",
      type="string",
      default="",
      dest="print_prefix",
      help="Prefix string for output")
  ).process(args=args)
  attributes_level = 0
  if (command_line.options.show_all_attributes):
    attributes_level = 3
  elif (command_line.options.show_some_attributes):
    attributes_level = 2
  elif (command_line.options.show_help):
    attributes_level = 1
  prefix = command_line.options.print_prefix
  for file_name in command_line.args:
    print prefix.rstrip()
    parameters = libtbx.phil.parse(
      file_name=file_name,
      converter_registry=converter_registry,
      process_includes=command_line.options.process_includes)
    parameters.show(
      out=sys.stdout,
      prefix=prefix,
      attributes_level=attributes_level,
      print_width=command_line.options.print_width)
    print prefix.rstrip()

if (__name__ == "__main__"):
  run(sys.argv[1:])
