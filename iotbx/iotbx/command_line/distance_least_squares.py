from cctbx.geometry_restraints import distance_least_squares
from iotbx.kriber import strudat
from iotbx.option_parser import iotbx_option_parser
import sys

def run(args, distance_cutoff=3.5):
  command_line = (iotbx_option_parser(
    usage="iotbx.distance_least_squares [options] studat_file [...]",
    description="Example: iotbx.distance_least_squares strudat --tag=SOD")
    .option(None, "--tag",
      action="store",
      type="string",
      dest="tag",
      help="tag as it appears in the strudat file")
    .option(None, "--n_trials",
      action="store",
      type="int",
      default=1,
      dest="n_trials",
      help="Number of trial per structure",
      metavar="INT")
  ).process(args=args)
  if (len(command_line.args) == 0):
    command_line.parser.show_help()
    return
  for file_name in command_line.args:
    strudat_entries = strudat.read_all_entries(open(file_name))
    for entry in strudat_entries.entries:
      if (    command_line.options.tag is not None
          and command_line.options.tag != entry.tag):
        continue
      print "strudat tag:", entry.tag
      print
      distance_least_squares.distance_and_repulsion_least_squares(
        si_structure=entry.as_xray_structure(),
        distance_cutoff=distance_cutoff,
        n_trials=command_line.options.n_trials,
        connectivities=entry.connectivities(all_or_nothing=True))

if (__name__ == "__main__"):
  run(sys.argv[1:])
