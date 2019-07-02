from __future__ import absolute_import, division, print_function
from cctbx import sgtbx
import scitbx.math
from libtbx.utils import format_cpu_times
import sys
from six.moves import range

special = {
  9: 8, 15: 8, 17: 4, 19: 12, 20: 4, 24: 12, 43: 4, 67: 4, 68: 4, 70: 18,
  73: 12, 74: 4, 80: 4, 85: 4, 86: 4, 88: 6, 91: 4, 92: 4, 95: 4, 96: 4,
  98: 4, 109: 4, 110: 4, 122: 4, 125: 4, 126: 4, 129: 4, 130: 4, 133: 4,
  134: 4, 137: 4, 138: 4, 141: 6, 142: 6, 151: 8, 152: 8, 153: 8, 154: 8,
  178: 8, 179: 8, 180: 8, 181: 8, 198: 21, 199: 21, 201: 18, 203: 18,
  205: 9, 206: 21, 210: 12, 212: 21, 213: 21, 214: 21, 220: 21, 222: 18,
  224: 18, 227: 18, 228: 18, 230: 21}

def run(args):
  if ("--full" in args):
    to_do = range(1,230+1)
  elif ("--special" in args):
    to_do = sorted(special.keys())
  else:
    to_do = [75, 151]
  for space_group_number in to_do:
    sgi = sgtbx.space_group_info(number=space_group_number)
    sgi.show_summary(prefix="")
    sys.stdout.flush()
    n_special = 0
    for m in scitbx.math.unimodular_generator(range=1).all():
      cb_op = sgtbx.change_of_basis_op(sgtbx.rt_mx(sgtbx.rot_mx(m,1), 1)) \
        .new_denominators(12, 144)
      cb_sgi = sgi.change_basis(cb_op=cb_op)
      cb_op_ref = cb_sgi.change_of_basis_op_to_reference_setting()
      ref_sgi = cb_sgi.change_basis(cb_op=cb_op_ref)
      assert ref_sgi.group() == sgi.group()
      c = cb_op_ref.c()
      if (c.r().is_unit_mx() and c.t().num() != (0,0,0)):
        n_special += 1
        cb_ref_sgi = sgi.change_basis(cb_op=cb_op_ref)
        print("  cb_op=%s -> %s" % (
          str(cb_op.c()), cb_ref_sgi.type().universal_hermann_mauguin_symbol()))
        sys.stdout.flush()
        # verify that c.t() is not an allowed origin shift
        assert cb_ref_sgi.group() != sgi.group()
    assert special.get(space_group_number, 0) == n_special
  print(format_cpu_times())

if (__name__ == "__main__"):
  run(sys.argv[1:])
