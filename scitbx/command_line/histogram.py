from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx.option_parser import option_parser
import sys

def process_file(file_object, n_slots, data_min, data_max, format_cutoffs):
  data = flex.double()
  for line in file_object.read().splitlines():
    data.append(float(line))
  print("total number of data points:", data.size())
  if (data_min is None): data_min = flex.min(data)
  if (data_max is None): data_max = flex.max(data)
  flex.histogram(
    data=data, n_slots=n_slots, data_min=data_min, data_max=data_max).show(
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
    .option(None, "--min",
      action="store",
      type="float",
      default=None,
      help="min data value in histogram",
      metavar="FLOAT")
    .option(None, "--max",
      action="store",
      type="float",
      default=None,
      help="max data value in histogram",
      metavar="FLOAT")
    .option("-f", "--format_cutoffs",
      action="store",
      type="str",
      default="%.8g",
      help="format specifier for cutoff values",
      metavar="STR")
  ).process(args=args)
  def pro(file_object):
    co = command_line.options
    process_file(
      file_object=file_object,
      n_slots=co.slots,
      data_min=co.min,
      data_max=co.max,
      format_cutoffs=co.format_cutoffs)
  if (len(command_line.args) == 0):
    pro(file_object=sys.stdin)
  else:
    for file_name in command_line.args:
      pro(file_object=open(file_name))

if (__name__ == "__main__"):
  run(sys.argv[1:])
