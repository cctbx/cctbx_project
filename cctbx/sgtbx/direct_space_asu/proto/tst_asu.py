import sys
import StringIO
import boost
import cctbx
from boost.rational import int as rint
from cctbx.sgtbx import space_group_info
from cctbx.sgtbx.direct_space_asu import proto as new_asu
from cctbx.crystal import direct_space_asu_float_asu
from libtbx.test_utils import approx_equal, is_below_limit
from libtbx.utils import format_cpu_times

# For usage type:
#   cctbx.python tst_asu.py -h

# test groups, tests are very slow, so keep very few groups here
SpaceGroups = ('P 1 1 21',  'P 21 21 21', 'I 1 m 1')
NSteps = 11 # number of grid points in one dimenssion

cout = StringIO.StringIO()

def step_v(n, mn, mx):
  step = ()
  box = ()
  for a,b in zip(mn,mx):
    box += tuple( [b-a] )
    step += tuple( [(b-a)/rint(n)] )
  v = rint( n*n*n ) / (box[0]*box[1]*box[2])
  return step, v

def loop_grid(asu, n, mn, mx, asu2=None):
  assert (asu2 is None) or isinstance(asu, new_asu.direct_space_asu)
  first_is_new = False
  step, vv = step_v(n, mn, mx)
  grid = ()
  for s in step:
    g = rint(1) / s
    g = g.numerator()/g.denominator() + 1
    assert g > 0, "Step = "+str(step)
    grid += tuple([g])
  mna = list(mn)
  mxa = list(mx)
  for i in xrange(3):
    mna[i] -= 5*step[i]
    mna[i] *= grid[i]
    mna[i] = mna[i].numerator()/mna[i].denominator()-1
    mxa[i] += 5*step[i]
    mxa[i] *= grid[i]
    mxa[i] = mxa[i].numerator()/mxa[i].denominator()+1
  print >>cout, "grid test  step= ", step, "  min= ", mna, "  max= ", mxa, \
      "   grid=", grid
  if isinstance(asu, new_asu.direct_space_asu):
    import copy
    asu_opt = copy.copy(asu)
    asu_opt.optimize_for_grid(grid)
    assert asu_opt.is_optimized() and (not asu.is_optimized())
    first_is_new = True
    max_p = asu_opt.get_optimized_grid_limits()
    print >>cout, "grid limits: ", max_p
  result = 0

  i = mna[0]
  while i <= mxa[0] :
    ii = rint(i,grid[0])
    j = mna[1]
    while j <= mxa[1] :
      jj = rint(j,grid[1])
      k = mna[2]
      while k <= mxa[2] :
        kk = rint(k,grid[2])
        b = asu.is_inside((ii,jj,kk)) # rational test
        if b :
          result += 1
        if first_is_new:
          num = (i, j, k)
          den = grid
          where = asu.where_is(num,den) # integer test
          #TODO: implement? assert b == asu.is_inside(num,den)
          assert ( b and (where==1 or where==-1)) or ( (not b) and where==0 )
          where_opt = asu_opt.where_is( num ) # optimized test
          assert where_opt == where
          #TODO: implement? assert b == asu_opt.is_inside(num)
        if asu2 is not None :
          assert b == asu2.is_inside( (ii,jj,kk) ) # rational test
        k += 1
      j += 1
    i += 1

  volume = result / vv
  return (result,volume)

def compare(spgr, n=NSteps, verbose=False):
  if verbose:
    print cout.getvalue()
  cout.truncate(0)
  grp = space_group_info(spgr)
  print >>cout, "Comparing asus for group: ", spgr, "  n-steps= ", n
  asuo = grp.direct_space_asu()
  print >>cout, "=== Original python asu ==="
  asuo.show_comprehensive_summary(cout)
  print >>cout, ":: raw facets"
  as_raw_asu(asuo, cout)
  mxo = tuple( asuo.box_max() )
  mno = tuple( asuo.box_min() )
  print >>cout, "box  min= ", mno, "   max= ", mxo
  ### NEW C++ ASU
  asun = new_asu.direct_space_asu(grp.type())
  print >>cout, "=== New C++ asu ==="
  asun.show_comprehensive_summary(cout)
  mnn = asun.box_min()
  mxn = asun.box_max()
  print >>cout, "box  min= ", mnn, "   max= ", mxn
  assert mnn == mno
  assert mxn == mxo
  old_vertices = asuo.volume_vertices()
  new_vertices = asun.volume_vertices() # C++  sorted list
  assert len(old_vertices) == len(new_vertices)
  # TODO: the following seems to use the same ordering operation
  # as mine in C++
  old_vertices = sorted(old_vertices)
  for a,b in zip(old_vertices,new_vertices):
    print >>cout, a, " == ", b
    assert a == b, str(a)+" != "+str(b)
  ins,v = loop_grid(asun, n, mnn, mxn, asuo)
  print >>cout, "N inside = ", ins, "   volume = ", v,  \
      "   expected volume = ", rint(1,grp.group().order_z())
  ### VOLUME ONLY
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
  unit_cell = grp.any_compatible_unit_cell(volume=5000.0)
  fasun = asun.as_float_asu( unit_cell, 1.0E-6);
  fasuo = cctbx.crystal.direct_space_asu_float_asu(
    unit_cell=unit_cell,
    cuts=[cut.as_float_cut_plane() for cut in asuo.cuts],
    is_inside_epsilon=1.e-6)
  # python sucks big: still (in 2.5) there is no float.epsilon/max/tiny/etc
  assert approx_equal(fasun.is_inside_epsilon(), fasuo.is_inside_epsilon(),
      1.0E-100)
  assert len(fasun.cuts()) == len(fasuo.cuts()), \
    "%d != %d"%(len(fasuo.cuts()),len(fasun.cuts()))
  assert ((len(fasun.cuts()) < 200) & (len(fasun.cuts()) > 3)), \
    len(fasun.cuts())
  if verbose:
    print cout.getvalue()
  cout.truncate(0)

def compare_groups(groups=SpaceGroups, n=NSteps, verbose=False):
  for sg in groups:
    compare(sg, n, verbose)


def as_raw(c ):
  result = ""
  if isinstance(c, cctbx.sgtbx.direct_space_asu.cut_plane.cut):
    if not c.inclusive:
      result = "+"
    result = result + "cut(" + str(c.n)+ ", " + str(c.c) +")"
    if c.has_cuts():
      result = result + " [ " + as_raw(c.cut_expr) + " ]"
  elif isinstance(c, cctbx.sgtbx.direct_space_asu.cut_plane.cut_expression):
    result = "(" + as_raw(c.lhs) + str(c.op) + as_raw(c.rhs) + ")"
  else:
    assert False, c.__class_
  return result

def as_raw_asu(asu, f=None):
  if (f == None): f = sys.stdout
  for cut in asu.cuts:
    print >>f, "  & ", as_raw(cut)


def run():
  import libtbx.option_parser as optparse
  parser = optparse.OptionParser()
  parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
      default=False, help="be verbose")
  parser.add_option("-g", "--space_group", action="store", type="string",
      dest="space_group", help="space group symbol or number or all or all530")
  parser.add_option("-n", "--n_steps", action="store", type="int",
      dest="n_steps", default=NSteps, help="number of grid points in one"
      +" dimension of the asu box")
  parser.add_option("--groups_file", action="store", type="string",
      dest="groups_file", help="file containing space groups, one per line")

  (opts, args) = parser.parse_args()

  groups = []
  if not opts.groups_file is None:
    tmp_file = open(opts.groups_file, "r")
    for line in tmp_file.readlines(): # newlines retained
      groups.append( line.strip() ) # removes whitespace in the begining and end
    tmp_file.close()
  if (opts.space_group is None) & (len(groups)==0) :
    groups.extend(SpaceGroups)
  elif opts.space_group == "all" :
    for isg in xrange(1,231):
      groups.append(str(isg))
  elif opts.space_group == "all530":
    it = cctbx.sgtbx.space_group_symbol_iterator()
    while( True ):
      symbol = it.next()
      # TODO: the following  does not work
      #if( symbol.number()==0 ):
      #  break
      groups.append(symbol.hermann_mauguin())
      if( symbol.number()==230 ):
        break
  elif not opts.space_group is None:
    groups.append(opts.space_group)

  print >>cout, "Number of groups: ", len(groups)
  print >>cout, "Options= ", opts
  compare_groups(groups, opts.n_steps, opts.verbose)
  print format_cpu_times()


if (__name__ == "__main__"):
  try:
    run()
  except :
    log = cout.getvalue()
    if len(log) != 0:
      print "<<<<<<<< Start Log:"
      print log
      print ">>>>>>>> End Log"
    raise

