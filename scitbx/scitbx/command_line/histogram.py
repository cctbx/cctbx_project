from scitbx.array_family import flex
from libtbx.optparse_wrapper import option_parser
import sys

def process_file(file_object, n_slots, format_cutoffs):
  data = flex.double()
  for line in file_object.read().splitlines():
    data.append(float(line))
  print "total number of data points:", data.size()
  flex.histogram(data=data, n_slots=n_slots).show(
    format_cutoffs=format_cutoffs)

def run(args, command_name="scitbx.histogram"):
  command_line = (option_parser(
    usage=command_name+" [options] [data_file ...]",
    description="Example: %s my_data --slots=20" % command_name)
    .option("-s", "--slots",
      action="store",
      type="int",
      default=10,
      help="number of histogram slots",
      metavar="INT")
    .option("-f", "--format_cutoffs",
      action="store",
      type="str",
      default="%.8g",
      help="format specifier for cutoff values",
      metavar="STR")
  ).process(args=args)
  n_slots = command_line.options.slots
  format_cutoffs = command_line.options.format_cutoffs
  if (len(command_line.args) == 0):
    process_file(
      file_object=sys.stdin,
      n_slots=n_slots,
      format_cutoffs=format_cutoffs)
  else:
    for file_name in command_line.args:
      process_file(
        file_object=open(file_name),
        n_slots=n_slots,
        format_cutoffs=format_cutoffs)

if (__name__ == "__main__"):
  run(sys.argv[1:])
