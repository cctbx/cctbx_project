import sys
from cctbx.sgtbx import  *
import cctbx.sgtbx.direct_space_asu.check_redundancies as asu_check

if (__name__=="__main__"):
  assert len(sys.argv) == 3
  spgr_s = sys.argv[1]
  n = int(sys.argv[2])
  print "Spacegroup = ", spgr_s,  "  Grid= ", n
  try:
    group = space_group_info(spgr_s)
  except:
    group = space_group_info(spgr_s, "Hall")
  asu = group.direct_space_asu()
  asu.show_comprehensive_summary()
  gridding = (n,n,n)
  asu_check.check_asu(
    space_group_number=group.group().type().number(),
    asu=asu,
    n=gridding,
    is_stripped_asu=False,
    soft_mode=False)
  print "Ok"

