import cctbx.crystal.pair_asu_table
from cctbx.crystal import distance_ls
from iotbx.kriber import strudat
from iotbx.option_parser import iotbx_option_parser
import sys

def run():
  command_line = (iotbx_option_parser(
    usage="iotbx.show_distances [options] studat_file [...]",
    description="Example: iotbx.show_distances strudat --tag=SOD")
    .option(None, "--tag",
      action="store",
      type="string",
      dest="tag",
      help="tag as it appears in the strudat file")
    .option(None, "--distance_cutoff",
      action="store",
      type="float",
      default=5,
      dest="distance_cutoff",
      help="Maximum distance to be considered",
      metavar="FLOAT")
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
      structure = entry.as_xray_structure()
      structure.show_summary().show_scatterers()
      print
      asu_mappings = structure.asu_mappings(
        buffer_thickness=command_line.options.distance_cutoff)
      pair_asu_table = cctbx.crystal.pair_asu_table.pair_asu_table(
        asu_mappings=asu_mappings)
      pair_asu_table.add_all_pairs(
        distance_cutoff=command_line.options.distance_cutoff)
      pairs = distance_ls.show_pairs(
        structure=structure,
        pair_asu_table=pair_asu_table)
      print

if (__name__ == "__main__"):
  run()
