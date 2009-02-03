import cctbx
from boost.rational import int as rint
from cctbx.sgtbx import space_group_info
from cctbx.sgtbx.direct_space_asu import proto as new_asu
import sys
import boost

# test groups
groups = ('P 1', 'P 1 21 1', 'P 1 1 21', 'P 21 1 1', 'P 21 21 21')

def step_v(n, mn, mx):
  step = ()
  box = ()
  for a,b in zip(mn,mx):
    box += tuple( [b-a] )
    step += tuple( [(b-a)/rint(n)] )
  v = rint( n*n*n ) / (box[0]*box[1]*box[2])
  return step, v

def loop_grid(asu, n, mn, mx):
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
        if asu.is_inside((i,j,k)) :
          result += 1
        k += step[2]
      j += step[1]
    i += step[0]

  volume = result / vv
  return (result,volume)

def compare(spgr):
  n = 11 # number of grid points in one dimenssion
  grp = space_group_info(spgr)
  asuo = grp.direct_space_asu()
  asun = new_asu.direct_space_asu(spgr)
  mxo = tuple( asuo.box_max() )
  mno = tuple( asuo.box_min() )
  mnn = asun.box_min()
  mxn = asun.box_max()
  # r = rint(1)
  # tmp1 = [r, r, r]
  # tmp2 = [r, r, r]
  # asun.box_corners( tmp1, tmp2 )
  assert mnn == mno
  assert mxn == mxo
  v = rint(1,grp.group().order_z())
  no,vo = loop_grid(asuo, n, mno, mxo)
  nn,vn = loop_grid(asun, n, mnn, mxn)
  step,nexp = step_v(n,mnn,mxn)
  nexp *= v
  ern = abs(nn-nexp)
  nermx = 6*n*n
  assert no == nn
  
def rank_err(spgr):
  n = 11 # number of grid points in one dimenssion
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

def test_1():
  for sg in groups:
    compare(sg)

def test_2():
  for i in xrange(1,231):
    compare(str(i))

def test_3():
  errs = dict()
  for i in xrange(1,231):
    errs[i] = rank_err(str(i))
  sorted = sort_by_value( errs )
  n = len(sorted)
  assert n > 0
  n -= 1
  print "Most disagreable space groups"
  for i in xrange(5):
    sg = sorted[n-i]
    print sg, "  ==  ", errs[ sg ]


def test_4():
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



def run():
  test_1()

if (__name__ == "__main__"):
  run()

