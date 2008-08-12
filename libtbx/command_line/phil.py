import libtbx.phil
from libtbx.utils import Sorry
from libtbx.option_parser import option_parser
import sys

def run(args, command_name="libtbx.phil", converter_registry=None):
  if (len(args) == 0): args = ["--help"]
  command_line = (option_parser(
    usage="%s [options] parameter_file ..." % command_name)
    .option(None, "--diff",
      action="store_true",
      help="Display only differences between the first file (master)"
           " and the combined definitions from all other files.")
    .option(None, "--show_help",
      action="store_true",
      help="Display help for each parameter if available.")
    .option(None, "--show_some_attributes",
      action="store_true",
      help="Display non-default attributes for each parameter.")
    .option(None, "--show_all_attributes",
      action="store_true",
      help="Display all attributes for each parameter.")
    .option(None, "--process_includes",
      action="store_true",
      help="Inline include files.")
    .option(None, "--print_width",
      action="store",
      type="int",
      help="Width for output",
      metavar="INT")
    .option(None, "--print_prefix",
      action="store",
      type="string",
      default="",
      help="Prefix string for output")
  ).process(args=args)
  co = command_line.options
  attributes_level = 0
  if (co.show_all_attributes):
    attributes_level = 3
  elif (co.show_some_attributes):
    attributes_level = 2
  elif (co.show_help):
    attributes_level = 1
  prefix = co.print_prefix
  file_names = command_line.args
  def parse(file_name):
    return libtbx.phil.parse(
      file_name=file_name,
      converter_registry=converter_registry,
      process_includes=co.process_includes)
  def show(scope):
    scope.show(
      out=sys.stdout,
      prefix=prefix,
      attributes_level=attributes_level,
      print_width=co.print_width)
  if (not co.diff):
    for file_name in file_names:
      print prefix.rstrip()
      show(scope=parse(file_name=file_name))
      print prefix.rstrip()
  else:
    if (len(file_names) < 2):
      raise Sorry("Option --diff requires at least two file names.")
    master = parse(file_name=file_names[0])
    show(scope=master.fetch_diff(
      sources=[parse(file_name=file_name)
        for file_name in file_names[1:]]))

if (__name__ == "__main__"):
  run(sys.argv[1:])
