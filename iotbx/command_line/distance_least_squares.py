""" Analyze non-bonded contacts in a model in strudat file
"""

from __future__ import absolute_import, division, print_function
def run(args, distance_cutoff=3.5):
  from iotbx.option_parser import option_parser
  command_line = (option_parser(
    usage="iotbx.distance_least_squares [options] studat_file [...]",
    description="Example: iotbx.distance_least_squares strudat --tag=SOD")
    .option(None, "--tag",
      action="store",
      type="string",
      help="tag as it appears in the strudat file")
    .option(None, "--repulsion_function",
      action="store",
      type="choice",
      choices=["gaussian", "cos", "prolsq"],
      default="gaussian",
      help="Nonbonded repulsion function type",
      metavar="gaussian|cos|prolsq")
    .option(None, "--bond_stretch_factor",
      action="store",
      type="float",
      default=0.1,
      help="Bond stretch factor used in max residual calculation"
           " for nonbonded cos or gaussian repulsion function",
      metavar="FLOAT")
    .option(None, "--n_trials",
      action="store",
      type="int",
      default=1,
      help="Number of trial per structure",
      metavar="INT")
    .option(None, "--n_macro_cycles",
      action="store",
      type="int",
      default=1,
      help="Number of macro cycles per trial",
      metavar="INT")
    .option(None, "--dev",
      action="store_true",
      default=False)
  ).process(args=args)
  if (len(command_line.args) == 0):
    command_line.parser.show_help()
    return
  co = command_line.options
  from cctbx.geometry_restraints import distance_least_squares as dls
  from iotbx.kriber import strudat
  for file_name in command_line.args:
    strudat_entries = strudat.read_all_entries(open(file_name))
    for entry in strudat_entries.entries:
      if (co.tag is not None and co.tag != entry.tag):
        continue
      print("strudat tag:", entry.tag)
      print()
      dls.distance_and_repulsion_least_squares(
        si_structure=entry.as_xray_structure(),
        distance_cutoff=distance_cutoff,
        nonbonded_repulsion_function_type=co.repulsion_function,
        nonbonded_max_residual_bond_stretch_factor=co.bond_stretch_factor,
        n_trials=co.n_trials,
        n_macro_cycles=co.n_macro_cycles,
        connectivities=entry.connectivities(all_or_nothing=True),
        dev=co.dev)

if (__name__ == "__main__"):
  import sys
  run(sys.argv[1:])
