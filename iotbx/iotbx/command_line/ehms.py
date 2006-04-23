from cctbx import sgtbx
from cctbx.sgtbx import pointgroup_tools as pgt
from cctbx import crystal
from iotbx.option_parser import iotbx_option_parser
from libtbx.utils import Sorry, date_and_time, multi_out
from cStringIO import StringIO

import sys, os

def ehms( args ):
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
            help = "Maximum delta/obliquity used in determining the lattice symmetry, using a modified Le-Page algorithm. Default is 5.0 degrees",
            metavar="FLOAT")
    .option(None, "--start_from_p1",
            action="store_true",
            dest="niggli",
            default=False,
            help="Reduce to Niggli cell and forget the input spacegroup before higher metric symmetry is sought.")
    .option(None, "--graph",
            action="store",
            dest="graph",
            default=None,
            help="A graphical representation of the graph will be written out. Requiers Graphviz to be installed and in path.")
    .option(None, "--but_not",
            action="store",
            dest="but_not",
            default=None,
            help="Remove this particular point group from the graph and find out the consequences.")

    ).process(args=args)

  log = multi_out()
  log.register(label="stdout", file_object=sys.stdout)


  if len(args)==0:
    command_line.parser.show_help()
    return

  if ( command_line.symmetry.unit_cell() == None ):
    print >> log
    print >> log, "Sorry: Unit cell not specified."
    print >> log
    command_line.parser.show_help()
    return

  if ( command_line.symmetry.space_group_info() == None ):
    print>> log
    print>> log,  "Sorry: space group not specified."
    print>> log
    command_line.parser.show_help()
    return

  if not ( command_line.symmetry.space_group().is_chiral() ):
    print >> log, "Sorry, Non chiral space groups not yet supported."
    return


  if command_line.options.niggli:
    print >> log, "*Unit cell will be niggli reduced and P1 will be assumed*"
    uc = command_line.symmetry.change_basis(
      command_line.symmetry.change_of_basis_op_to_niggli_cell() ).unit_cell()
    xs = crystal.symmetry( uc, "P 1" )
    command_line.symmetry = xs

  sg_explorer = pgt.space_group_graph_from_cell_and_sg(
      command_line.symmetry.unit_cell(),
      command_line.symmetry.space_group(),
      max_delta=command_line.options.max_delta)


  if not (command_line.options.but_not==None):
    # remnove this point group
    sg_explorer.pg_graph.remove_point_group_and_its_super_groups_from_graph(
      command_line.options.but_not )

  print >> log, "A summary of the constructed point group graph object is given below"
  print >> log, "===================================================================="
  print >> log
  print >> log, "----------------------"
  print >> log, "Input crystal symmetry"
  print >> log, "----------------------"
  print >> log, "Unit cell: ", command_line.symmetry.unit_cell().parameters()
  print >> log, "Unit cell volume: ", command_line.symmetry.unit_cell().volume()
  print >> log, "Space group: ", command_line.symmetry.space_group_info()
  print >> log
  print >> log
  print >> log, "--------------------------"
  print >> log, "Lattice symmetry deduction"
  print >> log, "--------------------------"
  print >> log, "Niggli cell: ", sg_explorer.xs_prim_set.unit_cell().parameters()
  print >> log, "Niggli cell volume: ", sg_explorer.xs_prim_set.unit_cell().volume()
  print >> log, "Niggli transformed input symmetry: ", sg_explorer.xs_prim_set.space_group_info()
  print >> log, "Symmetry of Niggli cell: ", sgtbx.space_group_info( group = sg_explorer.pg_high )
  print >> log
  print >> log
  print >> log, "All pointgroups that are both a subgroup of the lattice symmetry and"
  print >> log, "a supergroup of the Niggli transformed input symmetry wil now be listed,"
  print >> log, "as well as their minimal supergroups/maximal subgroups and symmetry"
  print >> log, "operators that generate them."
  print >> log, "For each pointgroup, a list of compatible spacegroups will be listed."
  print >> log, "Care is takebn that there are no sysmetatic absence violation with the "
  print >> log, "provided input spacegroup."
  print >> log



  sg_explorer.show(out=log)

  if command_line.options.graph is not None:
    # check if the program 'dot' is in the path
    graphviz = os.popen3("""
dot<< EOF
EOF""", "r"
                )
    if len(graphviz[2].readlines())>0:
      raise Sorry("The program dot has not been installed or is not in the path; get it from http://www.graphviz.org")

    # it seems graphviz is there, please proceed
    buffer = StringIO()

    sg_explorer.graphviz_pg_graph(out=buffer)
    exec_command = """dot -Tpng >""" + str(command_line.options.graph) + """ << EOF""" + """
""" + buffer.getvalue()
    exec_command += """
EOF
    """
    graphviz = os.popen( exec_command, "r" )
    print >> log, "A file named" ,  command_line.options.graph, " contains a graphical representation "
    print >> log, "of the point group relations."

if (__name__ == "__main__"):
  ehms(sys.argv[1:])
