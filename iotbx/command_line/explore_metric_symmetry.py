"""Explore Metric Symmetry. A list of possible unit cells and spacegroups is
given for the given specified unit cell and spacegroup combination. If a
second unit cell is given, linear combinations of the basis vector of one
unit cell are sought that match the other."""

from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.explore_metric_symmetry

from cctbx import sgtbx
from cctbx.sgtbx import pointgroup_tools as pgt
from cctbx.sgtbx import sub_lattice_tools as slt
from cctbx import crystal
from iotbx.option_parser import option_parser
from libtbx import easy_run
from libtbx.utils import Sorry, multi_out
from libtbx.str_utils import show_string
import libtbx.path
from six.moves import cStringIO as StringIO

import sys
from six.moves import zip

def do_pointgroup_tricks(input_uc,
                         input_ls,
                         max_delta,
                         out=None):
  if out is None:
    out = sys.stdout

  sg_explorer = pgt.space_group_graph_from_cell_and_sg(
      input_uc,
      input_ls,
      max_delta)


  print("A summary of the constructed point group graph object is given below", file=out)
  print("====================================================================", file=out)
  print(file=out)
  print("----------------------", file=out)
  print("Input crystal symmetry", file=out)
  print("----------------------", file=out)
  print("Unit cell: ", input_uc.parameters(), file=out)
  print("Unit cell volume: ", input_uc.volume(), file=out)
  print("Space group: ", sgtbx.space_group_info( group=input_ls ), file=out)
  print(file=out)
  print(file=out)
  print("--------------------------", file=out)
  print("Lattice symmetry deduction", file=out)
  print("--------------------------", file=out)
  print("Niggli cell: ", sg_explorer.xs_prim_set.unit_cell().parameters(), file=out)
  print("Niggli cell volume: ", sg_explorer.xs_prim_set.unit_cell().volume(), file=out)
  print("Niggli transformed input symmetry: ", sg_explorer.xs_prim_set.space_group_info(), file=out)
  print("Symmetry of Niggli cell: ", sgtbx.space_group_info( group = sg_explorer.pg_high ), file=out)
  print(file=out)
  print(file=out)
  print("All pointgroups that are both a subgroup of the lattice symmetry and", file=out)
  print("a supergroup of the Niggli transformed input symmetry wil now be listed,", file=out)
  print("as well as their minimal supergroups/maximal subgroups and symmetry", file=out)
  print("operators that generate them.", file=out)
  print("For each pointgroup, a list of compatible spacegroups will be listed.", file=out)
  print("Care is taken that there are no systematic absence violation with the ", file=out)
  print("provided input spacegroup.", file=out)
  print(file=out)
  out.flush()

  sg_explorer.show(out=out)

  # return the object
  return sg_explorer




def make_graph_of_graph(pg_object,
                        file_name,
                        out=None):
  if out is None:
    out = sys.stdout

  dot_path = libtbx.path.full_command_path(command="dot")
  if (dot_path is None):
    raise Sorry("""\
The program "dot" is not on PATH:
  For information about "dot" visit: http://www.graphviz.org/""")

  buffer = StringIO()
  pg_object.graphviz_pg_graph(out=buffer)
  command = "%s -Tpng > %s" % (show_string(dot_path), show_string(file_name))
  # XXX warning - Fontconfig error messages cause raise_if_errors_or_output()
  # to crash even if 'dot' ran successfully.
  rc = easy_run.fully_buffered(
    command=command,
    stdin_lines=buffer.getvalue().splitlines())#.raise_if_errors_or_output()
  if (rc.return_code != 0):
    raise RuntimeError("Fatal error running %s:\n%s" % (dot_path,
      "\n".join(rc.stderr_lines)))
  print("A file named", show_string(file_name), \
    "contains a graphical representation ", file=out)
  print("of the point group relations.", file=out)


def run(args, command_name="phenix.explore_metric_symmetry"):
  command_line = (
    option_parser(
    usage=command_name+" [options]",
    description="""\
Explore Metric Symmetry. A list of possible unit cells and spacegroups is
given for the given specified unit cell and spacegroup combination. If a
second unit cell is given, linear combinations of the basis vector of one
unit cell are sought that match the other.""")

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
            default=None,
            help="A graphical representation of the graph will be written out."
                 " Requires Graphviz to be installed and on PATH.")

    .option(None, "--centring_type",
            action="store",
            type="str",
            help="Centring type, choose from P,A,B,C,I,R,F")

    .option(None, "--other_unit_cell",
            action="store",
            type="str",
            help="Other unit cell, for unit cell comparison",
            metavar="10,20,30,90,103.7,90")

    .option(None, "--other_space_group",
            action="store",
            type="str",
            help="space group for other_unit_cell, for unit cell comparison")

    .option(None, "--other_centring_type",
            action="store",
            type="str",
            help="Centring type, choose from P,A,B,C,I,R,F")

    .option(None, "--no_point_group_graph",
            action="store_true",
            dest="pg_graph",
            default=False,
            help="Do not carry out the construction of a point group graph." )

    .option(None, "--relative_length_tolerance",
            action="store",
            type="float",
            help="Tolerance for unit cell lengths to be considered equal-ish.",
            default=0.10,
            metavar="FLOAT",
            dest="rel_length_tol")

    .option(None, "--absolute_angle_tolerance",
            action="store",
            dest="abs_angle_tol",
            type="float",
            default=10.0,
            metavar="FLOAT",
            help="Angular tolerance in unit cell comparison")

     .option(None, "--max_order",
             action="store",
             type="int",
             default=1,
             metavar="INT",
             help="Maximum volume change for target cell" )
    ).process(args=args)

  log = multi_out()
  log.register(label="stdout", file_object=sys.stdout)

  allowed_centring_types={"P":"Primitive",
                          "A":"A centered",
                          "B":"B centered",
                          "C":"C centered",
                          "I":"Body centered",
                          "R":"Rombohedral",
                          "F":"Face centered"}
  if command_line.options.centring_type is not None:
    if command_line.options.centring_type not in allowed_centring_types:
      print("Sorry, the centring type %s is not known."%(command_line.options.centring_type), file=log)
      print("Choose from P,A,B,C,I,R,F ", file=log)
      return

  xs = None
  other_xs = None

  if len(args)==0:
    command_line.parser.show_help()
    return

  if ( command_line.symmetry.unit_cell() == None ):
    print(file=log)
    print("Sorry: Unit cell not specified.", file=log)
    print(file=log)
    command_line.parser.show_help()
    return

  if command_line.options.centring_type is None:
    if ( command_line.symmetry.space_group_info() == None ):
      print(file=log)
      print("Sorry: centring type or space group not specified.", file=log)
      print(file=log)
      command_line.parser.show_help()
      return
  if command_line.symmetry.space_group_info()  is not None:
    if not ( command_line.symmetry.space_group().is_chiral() ):
      print("Sorry, Non chiral space groups not yet supported.", file=log)
      return

  if command_line.options.centring_type is not None:
    xs  = crystal.symmetry(
      unit_cell=command_line.symmetry.unit_cell(),
      space_group_symbol="Hall: %s 1" %( command_line.options.centring_type )
      )
    command_line.symmetry = xs

  if command_line.options.niggli:
    print("*Unit cell will be niggli reduced and P1 will be assumed*", file=log)
    uc = command_line.symmetry.change_basis(
      command_line.symmetry.change_of_basis_op_to_niggli_cell() ).unit_cell()
    command_line.symmetry = crystal.symmetry( uc, "P 1" )

  xs = command_line.symmetry

  ############################################################################
  # ABOVE IS JUST INPUT PARSING, NOW THE ACTUAL STUFF HAPPENS
  ############################################################################


  if not command_line.options.pg_graph:
    ##############################
    #   get a point group graph  #
    ##############################

    pg_object = do_pointgroup_tricks( xs.unit_cell(),
                                      xs.space_group(),
                                      command_line.options.max_delta,
                                      log )

    ################################################
    #  make a graphical representation if desired  #
    ################################################

    if command_line.options.graph is not None:
      make_graph_of_graph(pg_object,
                          command_line.options.graph,
                          log)


  #########################################
  #  Check if other cell has been defined #
  #########################################

  if command_line.options.other_unit_cell is not None:
    print("A second unit cell has been specified. ", file=log)
    other_xs = None

    if command_line.options.other_space_group is None:
      if command_line.options.other_centring_type is None:
        raise Sorry("No space group or centring type for other cell specified.")
      else:
        other_xs = crystal.symmetry( command_line.options.other_unit_cell,
                                     space_group_symbol="Hall: %s 1" %( command_line.options.other_centring_type )
                                   )
    else:
      other_xs = crystal.symmetry( command_line.options.other_unit_cell,
                                   space_group_symbol=command_line.options.other_space_group
                                 )

    # get the graph is desired
    if not command_line.options.pg_graph:
      other_pg_object = do_pointgroup_tricks( other_xs.unit_cell(),
                                              other_xs.space_group(),
                                              command_line.options.max_delta,
                                              log )
    # do the unit cell comparison
    print(file=log)
    print(file=log)
    print("Unit cell comparison", file=log)
    print("--------------------", file=log)
    print(file=log)
    print("The unit cells will be compared. The smallest niggli cell,", file=log)
    print("will be used as a (semi-flexible) lego-block to see if it", file=log)
    print("can construct the larger Niggli cell.", file=log)
    print(file=log)
    print(file=log)

    order = command_line.options.max_order

    if order==1:
      sl_object =  slt.compare_lattice(xs_a=xs,
                                       xs_b=other_xs,
                                       max_delta=command_line.options.max_delta,
                                       out=log,
                                       relative_length_tolerance=command_line.options.rel_length_tol,
                                       absolute_angle_tolerance=command_line.options.abs_angle_tol)
    else:

      tmp_a = xs.change_basis( xs.change_of_basis_op_to_niggli_cell() )
      tmp_b = other_xs.change_basis( other_xs.change_of_basis_op_to_niggli_cell() )
      modified_xs = None
      order = command_line.options.max_order
      lego_block = None
      if ( tmp_a.unit_cell().volume() > tmp_b.unit_cell().volume() ):
        modified_xs = slt.make_list_of_target_xs_up_to_order( xs, order )
        lego_block = other_xs
      else:
        modified_xs = slt.make_list_of_target_xs_up_to_order( other_xs, order )
        lego_block = xs

      print(file=log)
      print("Volume change of largest niggli cell requested via keyword --max_order", file=log)
      print(file=log)
      print("Input crystal symmetry is tranformed to niggli setting using the operator:", file=log)
      print(modified_xs.basic_to_niggli_cb_op.as_xyz(), file=log)
      print(file=log)
      print("Comparisons for various sublattices of the target cell are listed", file=log)
      print(file=log)

      for tmp_xs,cb_op,mat in zip(modified_xs.xs_list,
                                  modified_xs.extra_cb_op,
                                  modified_xs.matrices ):
        mat=mat.as_list_of_lists()
        print("===================================================================", file=log)
        print("Niggli cell is expanded using matrix:", file=log)
        print(file=log)
        print(r"               /%4i %4i %4i  \  "%(mat[0][0],mat[0][1],mat[0][2]), file=log)
        print("          M =  |%4i %4i %4i  |  "%(mat[1][0],mat[1][1],mat[1][2]), file=log)
        print(r"               \%4i %4i %4i  /  "%(mat[2][0],mat[2][1],mat[2][2]), file=log)
        print(file=log)
        print("Change of basis operator to reference setting:", file=log)
        print("    ", cb_op.as_xyz(), file=log)
        print("resulting crystal symmetry:", file=log)
        tmp_xs.show_summary(f=log,prefix="   ")
        print(file=log)
        print(file=log)
        sl_object =  slt.compare_lattice(xs_a=tmp_xs,
                                         xs_b=lego_block,
                                         max_delta=command_line.options.max_delta,
                                         out=log,
                                         relative_length_tolerance=command_line.options.rel_length_tol,
                                         absolute_angle_tolerance=command_line.options.abs_angle_tol)


if (__name__ == "__main__"):
  run(args=sys.argv[1:])
