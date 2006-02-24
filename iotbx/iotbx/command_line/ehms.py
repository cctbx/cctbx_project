from cctbx.sgtbx import pointgroup_tools as pgt
from cctbx import crystal
from iotbx.option_parser import iotbx_option_parser
import sys, os

def ehms( args ):
  print "Arguments: ",
  for arg in args:
    print arg,
  print
  print

  command_line = (
    iotbx_option_parser(
    usage="iotbx.ehms [options]",
    description="Explore Higher Metric Symmetry. A list of possible unit cells and spacegroups is given for the given specified unit cell and spacegroup combination.")
    .enable_symmetry_comprehensive()
    .option(None, "--max_delta",
            action = "store",
            type="float",
            default=5.0,
            dest = "max_delta",
            help = "Maximum delta for modified Le-Page algorithm",
            metavar="FLOAT")
    .option(None, "--niggli",
            action="store_true",
            dest="niggli",
            default=False,
            help="Reduce to niggli cell and assume P1")

    ).process(args=args)

  if ( command_line.symmetry.unit_cell() == None ):
    print "Sorry: Unit cell not specified."
    command_line.parser.show_help()
    return

  if ( command_line.symmetry.space_group() == None ):
    print "Sorry: space group not specified."
    command_line.parser.show_help()
    return

  if not ( command_line.symmetry.space_group().is_chiral() ):
    print "Sorry, Non chiral space groups not yet supported."
    return


  if command_line.options.niggli:
    print "*Unit cell will be niggli reduced and P1 will be assumed*"
    uc = command_line.symmetry.change_basis(
      command_line.symmetry.change_of_basis_op_to_niggli_cell() ).unit_cell()
    xs = crystal.symmetry( uc, "P 1" )
    command_line.symmetry = xs

  print
  print
  sg_explorer = pgt.space_group_graph_from_cell_and_sg(
      command_line.symmetry.unit_cell(),
      command_line.symmetry.space_group(),
      max_delta=command_line.options.max_delta)


  print "A summary of the constructed point group graph object is given below"
  print "===================================================================="
  print
  print "----------------------"
  print "Input crystal symmetry"
  print "----------------------"
  print "Unit cell: ", command_line.symmetry.unit_cell().parameters()
  print "Unit cell volume: ", command_line.symmetry.unit_cell().volume()
  print "Space group: ", command_line.symmetry.space_group_info()
  print
  print



  sg_explorer.show()


if (__name__ == "__main__"):
  ehms(sys.argv[1:])
