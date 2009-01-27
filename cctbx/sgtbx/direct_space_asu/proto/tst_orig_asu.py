import cctbx
from cctbx.sgtbx import direct_space_asu
from boost import rational
from cctbx.sgtbx.direct_space_asu.cut_plane import cut
from cctbx.sgtbx import space_group
from cctbx.sgtbx import space_group_info
import sys
import boost

def loop_grid(asu, n):
  asu.show_comprehensive_summary()
  mx = asu.box_max()
  mn = asu.box_min()
  mx = (mx[0],mx[1],mx[2])
  mn = (mn[0],mn[1],mn[2])
  print "low corner= ", mn, "    high corner= ", mx
  step = ()
  box = ()
  for a,b in zip(mn,mx):
    box += tuple( [b-a] )
    step += tuple( [(b-a)/boost.rational.int(n)] )
  print "step = ", step
  
  result = 0
  i = mn[0]
  while i <= mx[0] :
    j = mn[1]
    while j <= mx[1] :
      k = mn[2]
      while k <= mx[2] :
        if asu.is_inside((i,j,k)) :
          result += 1
        k += step[2]
      j += step[1]
    i += step[0]

  volume = box[0]*box[1]*box[2]
  volume *= result
  volume /= (n*n*n)
  print "N asu points: ", result, "  volume= ", volume
  return result


def run():
  n = 20
  spgr = "P 21 21 21"
  if len(sys.argv)>1 :
    n = int(sys.argv[1])
  if len(sys.argv)>2 :
    spgr = sys.argv[2]

  print "space group= ", spgr,  "  nsteps= ", n
  grp = space_group_info(spgr)
  print grp.type().hall_symbol()
  asu = grp.direct_space_asu()

  ins = loop_grid(asu, n)
  print "N inside = ", ins

if (__name__ == "__main__"):
  run()

