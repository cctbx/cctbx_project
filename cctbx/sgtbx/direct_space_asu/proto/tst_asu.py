import cctbx
from boost.rational import int as rint
from cctbx.sgtbx import space_group_info
from cctbx.sgtbx.direct_space_asu import proto as new_asu
import sys
import boost

# Usage:
#   python tst_asu.py [ action ] [ spacegroup [ nsteps] ]
#          action =
#                    none - silently compares asus of a few selected spacegroups
#                    all  -  silently compares asus of all 230 space groups
#                    print_asu - tests one asu
#                    print_original_asu - tests one original asu
#                    print_inconsistent - prints 5 most inconsistent asus
# test groups
Groups = ('P 1', 'P 1 21 1', 'P 1 1 21', 'P 21 1 1', 'P 21 21 21')
NSteps = 11 # number of grid points in one dimenssion

def step_v(n, mn, mx):
  step = ()
  box = ()
  for a,b in zip(mn,mx):
    box += tuple( [b-a] )
    step += tuple( [(b-a)/rint(n)] )
  v = rint( n*n*n ) / (box[0]*box[1]*box[2])
  return step, v

def loop_grid(asu, n, mn, mx, asu2=None):
  step, vv = step_v(n, mn, mx)
  result = 0
  mna = list(mn)
  mxa = list(mx)
  for i in xrange(3):
    mna[i] -= step[i]
    mxa[i] += step[i]

  i = mna[0]
  while i <= mxa[0] :
    j = mna[1]
    while j <= mxa[1] :
      k = mna[2]
      while k <= mxa[2] :
        b = asu.is_inside((i,j,k))
        if b :
          result += 1
        if asu2 is not None :
          assert b == asu2.is_inside( (i,j,k) )
        k += step[2]
      j += step[1]
    i += step[0]

  volume = result / vv
  return (result,volume)

def compare(spgr, n=NSteps):
  grp = space_group_info(spgr)
  asuo = grp.direct_space_asu()
  asun = new_asu.direct_space_asu(spgr)
  mxo = tuple( asuo.box_max() )
  mno = tuple( asuo.box_min() )
  mnn = asun.box_min()
  mxn = asun.box_max()
  assert mnn == mno
  assert mxn == mxo
  v = rint(1,grp.group().order_z())
  loop_grid(asun, n, mnn, mxn, asuo)
  asun.volume_only()
  asuo = asuo.volume_only()
  mxo2 = tuple( asuo.box_max() )
  mno2 = tuple( asuo.box_min() )
  mnn2 = asun.box_min()
  mxn2 = asun.box_max()
  assert mxo2 == mxo
  assert mno2 == mno
  assert mxn2 == mxn
  assert mnn2 == mnn
  loop_grid(asun, n, mnn, mxn, asuo)

def rank_err(spgr, n=NSteps):
  grp = space_group_info(spgr)
  asun = new_asu.direct_space_asu(spgr)
  mnn = asun.box_min()
  mxn = asun.box_max()
  v = rint(1,grp.group().order_z())
  nn,vn = loop_grid(asun, n, mnn, mxn)
  step,nexp = step_v(n,mnn,mxn)
  nexp *= v
  ern = abs(nn-nexp)
  nermx = 6*n*n
  return float(ern)/float(nermx)

def sort_by_value(d):
    """ Returns the keys of dictionary d sorted by their values """
    items=d.items()
    backitems=[ [v[1],v[0]] for v in items]
    backitems.sort()
    return [ backitems[i][1] for i in range(0,len(backitems))]

def test_1(n=NSteps):
  for sg in Groups:
    compare(sg, n)

def test_2(n=NSteps):
  for i in xrange(1,231):
    compare(str(i), n)

def test_3(n=NSteps):
  errs = dict()
  for i in xrange(1,231):
    errs[i] = rank_err(str(i), n)
  sorted = sort_by_value( errs )
  ns = len(sorted)
  assert ns > 0
  ns -= 1
  print "Most disagreable space groups\nGroup     relative error"
  for i in xrange(5):
    sg = sorted[ns-i]
    print sg, "  ==  ", errs[ sg ]


def test_4(spgr, n, original=False):
  print "space group= ", spgr,  "  nsteps= ", n
  grp = space_group_info(spgr)
  if original:
    asu = grp.direct_space_asu()
  else:
    asu = new_asu.direct_space_asu(grp.type())
  asu.show_comprehensive_summary()
  mn = asu.box_min()
  mx = asu.box_max()
  print "asu box = ", mn, " : ", mx
  ins,v = loop_grid(asu, n, mn, mx)
  print "N inside = ", ins, "   volume = ", v,  "   expected volume = ", rint(1,grp.group().order_z())
  print "\n ---------- volume_only ----------------\n"
  if original:
    asu = asu.volume_only()
  else:
    asu.volume_only()
  asu.show_comprehensive_summary()
  mn = asu.box_min()
  mx = asu.box_max()
  print "asu box = ", mn, " : ", mx
  ins,v = loop_grid(asu, n, mn, mx)
  print "N inside = ", ins, "   volume = ", v,  "   expected volume = ", rint(1,grp.group().order_z())



def run():
  key = ""
  if len(sys.argv)>1 :
    key = sys.argv[1]
  if key in ("print_asu", "print_original_asu") :
    spgr = "P 21 21 21"
    n = NSteps
    if len(sys.argv)>2 :
      spgr = sys.argv[2]
    if len(sys.argv)>3 :
      n = int(sys.argv[3])
    test_4(spgr, n, key=="print_original_asu" )
  elif key == "print_inconsistent" :
    test_3()
  elif key == "all" :
    test_2()
  else :
    test_1()

if (__name__ == "__main__"):
  run()

