# List all axes in the unit cell.

# usage:
#   python all_axes.py     - show axes for the 230 reference settings.
#   python all_axes.py P2  - show axes for (e.g.) space group P2

# XXX Some further refinement is required:
# XXX   - List only the axes of highest order (e.g. only 4, not 4 and 2).
# XXX   - List only the axes with the smallest intrinsic component
# XXX     (e.g. list only 3(1), not both 3(1) and 3(2)).
# XXX See also: comment regarding shift_range below.

from cctbx_boost import sgtbx
from cctbx.misc.python_utils import list_plus

def StrEV(EV):
  return "[%d,%d,%d]" % EV

def RTMxAnalysis(M):
  RI = M.getRotMxInfo()
  TI = M.analyzeTpart()
  t_intrinsic = TI.IntrinsicPart().modPositive()
  t_shift = TI.OriginShift().modPositive()
  if (RI.Rtype() == 1):
    return ("1", "-", "-", "-")
  if (RI.Rtype() == -1):
    return (str(RI.Rtype()), "-", "-", "(%s)" % (t_shift,))
  if (abs(RI.Rtype()) == 2):
    return (str(RI.Rtype()),
            StrEV(RI.EV()),
            "(%s)" % (t_intrinsic,),
            "(%s)" % (t_shift,))
  return (str(RI.Rtype()),
          StrEV(RI.EV()),
          "(%s)" % (t_intrinsic,),
          "(%s)" % (t_shift,))

def list_all_axes(space_group_symbol):
  shift_range = 1 # XXX Works for the 230 reference settings; it is not
                  # XXX clear to me (rwgk) what value is needed in general.
  space_group = sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols(space_group_symbol))
  space_group.makeTidy()
  print space_group.Info().BuildLookupSymbol()
  axes_dict = {}
  for smx in space_group:
    r, t = smx.as_tuple(1, sgtbx.STBF)
    shift = [0,0,0]
    for shift[0] in xrange(-shift_range,shift_range+1):
      for shift[1] in xrange(-shift_range,shift_range+1):
        for shift[2] in xrange(-shift_range,shift_range+1):
          ts = list_plus(t, [s * sgtbx.STBF for s in shift])
          m = sgtbx.RTMx_from_tuple((r, ts), 1, sgtbx.STBF)
          axes_dict[RTMxAnalysis(m)] = 0
  axes_list = axes_dict.keys()
  axes_list.sort()
  for a in axes_list:
    print a
  print

if (__name__ == "__main__"):
  import sys
  if (len(sys.argv) == 1):
    for i in xrange(230):
      list_all_axes(i + 1)
  else:
    for symbol in sys.argv[1:]:
      list_all_axes(symbol)
