from __future__ import absolute_import, division, print_function
# List all axes in the unit cell.

# usage:
#   cctbx.python all_axes.py     - show axes for the 230 reference settings.
#   cctbx.python all_axes.py P2  - show axes for (e.g.) space group P2

# XXX Some further refinement is required:
# XXX   - List only the axes of highest order (e.g. only 4, not 4 and 2).
# XXX   - List only the axes with the smallest intrinsic component
# XXX     (e.g. list only 3(1), not both 3(1) and 3(2)).
# XXX See also: comment regarding shift_range below.

from cctbx import sgtbx
import sys
from six.moves import range

def str_ev(ev):
  return "[%d,%d,%d]" % ev

def rt_mx_analysis(s):
  r_info = sgtbx.rot_mx_info(s.r())
  t_info = sgtbx.translation_part_info(s)
  t_intrinsic = str(t_info.intrinsic_part().mod_positive())
  t_shift = str(t_info.origin_shift().mod_positive())
  if (r_info.type() == 1):
    return ("1", "-", "-", "-")
  if (r_info.type() == -1):
    return (str(r_info.type()), "-", "-", "(%s)" % (t_shift,))
  if (abs(r_info.type()) == 2):
    return (str(r_info.type()),
            str_ev(r_info.ev()),
            "(%s)" % (t_intrinsic,),
            "(%s)" % (t_shift,))
  return (str(r_info.type()),
          str_ev(r_info.ev()),
          "(%s)" % (t_intrinsic,),
          "(%s)" % (t_shift,))

def list_all_axes(space_group_symbol=None, space_group_info=None):
  assert space_group_symbol is None or space_group_info is None
  shift_range = 1 # XXX Works for the 230 reference settings; it is not
                  # XXX clear to me (rwgk) what value is needed in general.
  if (space_group_symbol is not None):
    space_group_info = sgtbx.space_group_info(symbol=space_group_symbol)
  space_group_info.show_summary()
  print()
  print("Rotation type, Axis direction, Intrinsic part, Origin shift")
  axes_dict = {}
  for s in space_group_info.group():
    r = s.r()
    t = s.t()
    shift = [0,0,0]
    for shift[0] in range(-shift_range,shift_range+1):
      for shift[1] in range(-shift_range,shift_range+1):
        for shift[2] in range(-shift_range,shift_range+1):
          ts = t.plus(sgtbx.tr_vec(shift, 1)).new_denominator(t.den())
          ss = sgtbx.rt_mx(r, ts)
          axes_dict[rt_mx_analysis(ss)] = 0
  axes_list = list(axes_dict.keys())
  axes_list.sort()
  for a in axes_list:
    print(a)
  print()

def run():
  if (len(sys.argv) == 1):
    for i in range(230):
      list_all_axes(i + 1)
  else:
    for symbol in sys.argv[1:]:
      list_all_axes(symbol)

if (__name__ == "__main__"):
  run()
