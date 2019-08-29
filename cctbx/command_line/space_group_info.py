from __future__ import absolute_import, division, print_function
from six.moves import range
def run(args):
  if (len(args) == 0): args = ["--help"]
  from libtbx.option_parser import libtbx_option_parser
  import libtbx.load_env
  command_line = (libtbx_option_parser(
    usage="%s [options] space_group_symbol ..." % libtbx.env.dispatcher_name)
      .option(None, "--primitive",
        action="store_true",
        help="Convert centred space groups to primitive setting.")
      .option(None, "--symops",
        action="store_true",
        help="Show list of symmetry operations.")
  ).process(args=args)
  co = command_line.options
  from cctbx import sgtbx
  args = command_line.args
  if (args == ["all"]):
    args = [str(no) for no in range(1,230+1)]
  for symbol in args:
    sgi = sgtbx.space_group_info(symbol=symbol)
    if (co.primitive):
      sgi = sgi.primitive_setting()
    sgi.show_summary()
    gr = sgi.group()
    print("  Crystal system:", gr.crystal_system())
    print("  Point group type:", gr.point_group_type())
    print("  Laue group type:", gr.laue_group_type())
    print("  Number of symmetry operations:", gr.order_z())
    print("  Lattice centering operations:", gr.n_ltr())
    if (gr.n_ltr() != 1):
      print("  Number of symmetry operations in primitive setting:", \
        gr.order_p())
    if (gr.is_centric()): s = gr(0, 1, 0)
    else:                 s = None
    print("  Center of inversion:", s)
    print("  Dimensionality of continuous allowed origin shifts:", \
      sgi.number_of_continuous_allowed_origin_shifts())
    ssi_vm = sgi.structure_seminvariants().vectors_and_moduli()
    print("  Structure-seminvariant vectors and moduli:", \
      len(ssi_vm))
    if (len(ssi_vm) != 0):
      for vm in ssi_vm:
        print("    (%2d, %2d, %2d)" % vm.v, "%2d" % vm.m)
    print("  Direct-space asymmetric unit:")
    dau = sgi.direct_space_asu()
    print("    Number of faces: %d" % len(dau.cuts))
    for cut in dau.cuts:
      print("    " + cut.as_xyz())
    print("  ADP constraint matrix:")
    _ = sgi.group().adp_constraints().gradient_sum_matrix()
    from scitbx import matrix
    print(matrix.rec(tuple(_), _.focus()).mathematica_form(
      one_row_per_line=True,
      prefix="    "))
    if (co.symops):
      print("  List of symmetry operations:")
      for s in sgi.group():
        print("    %s" % str(s))
    print()

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
