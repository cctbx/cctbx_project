from cctbx import sgtbx

def verify_a_theorem_in_primitive_settings():
  print ("For any space group in a primitive setting, if it were to contain "
         "two elements (R|t) and (R|t+d) where t is the intrinsic part, "
         "then d is a lattice translation.")
  print "Let's try to falsify it by finding counter-examples..."
  for symbol in sgtbx.space_group_symbol_iterator():
    sg = sgtbx.space_group(symbol.hall())
    z2p_op = sg.z2p_op()
    sg_p = sg.change_basis(z2p_op)
    symbol_printed = False
    seen_rot_mx = []
    for op in sg_p:
      tr_info = sgtbx.translation_part_info(op)
      o = tr_info.origin_shift()
      if not o.is_zero() and op.r() in seen_rot_mx:
        if not symbol_printed:
          print "%s (%i) [ %s ]" % (
            symbol.hall(), symbol.number(), z2p_op.as_xyz())
          symbol_printed = True
        print "\t%s --> %s" % (op.r().as_xyz(), o)
      seen_rot_mx.append(op.r())
  print "None!"

def run():
  verify_a_theorem_in_primitive_settings()

if __name__ == '__main__':
  run()
