from cctbx.crystal import distance_ls
from iotbx.kriber import strudat
from iotbx.option_parser import iotbx_option_parser
import sys

def run(distance_cutoff=3.5, nonbonded_distance_cutoff=5):
  command_line = (iotbx_option_parser(
    usage="iotbx.distance_least_squares [options] studat_file [...]",
    description="Example: python distance_ls.py strudat --tag=SOD")
    .option(None, "--tag",
      action="store",
      type="string",
      dest="tag",
      help="tag as it appears in the strudat file")
  ).process()
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
      distance_ls.distance_and_repulsion_least_squares(
        si_structure=entry.as_xray_structure(),
        distance_cutoff=distance_cutoff,
        nonbonded_distance_cutoff=nonbonded_distance_cutoff,
        connectivities=entry.connectivities(all_or_nothing=0001))

if (__name__ == "__main__"):
  run()
