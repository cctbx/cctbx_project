def run(args):
  from cctbx import sgtbx
  if (args == ["all"]):
    args = [str(no) for no in xrange(1,230+1)]
  for symbol in args:
    sgi = sgtbx.space_group_info(symbol=symbol)
    sgi.show_summary()
    gr = sgi.group()
    print "  Crystal system:", gr.crystal_system()
    print "  Point group type:", gr.point_group_type()
    print "  Laue group type:", gr.laue_group_type()
    print "  Symmetry operations:", gr.order_z()
    print "  Lattice centering operations:", gr.n_ltr()
    if (gr.n_ltr() != 1):
      print "  Symmetry operations in primitive setting:", gr.order_p()
    if (gr.is_centric()): s = gr(0, 1, 0)
    else:                 s = None
    print "  Center of inversion:", s
    print "  Dimensionality of continuous allowed origin shifts:", \
      sgi.number_of_continuous_allowed_origin_shifts()
    ssi_vm = sgi.structure_seminvariants().vectors_and_moduli()
    print "  Structure-seminvariant vectors and moduli:", \
      len(ssi_vm)
    if (len(ssi_vm) != 0):
      for vm in ssi_vm:
        print "    (%2d, %2d, %2d)" % vm.v, "%2d" % vm.m
    print "  Direct-space asymmetric unit:"
    dau = sgi.direct_space_asu()
    print "    Number of faces: %d" % len(dau.cuts)
    for cut in dau.cuts:
      print "    " + cut.as_xyz()
    print

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
