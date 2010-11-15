from cctbx.sgtbx.direct_space_asu import reference_table
from cctbx.web.asu_gallery import jv_asu
from cctbx.web.asu_gallery import jv_index
from cctbx import sgtbx
from scitbx import matrix
from scitbx.python_utils import dicts
from scitbx.python_utils import command_line
from libtbx import easy_run
from boost import rational
import sys, os

class colored_grid_point(object):

  def __init__(self, site, color):
    self.site = tuple(site)
    self.color = color

def sample_asu(asu, n=(12,12,12), volume=False, is_stripped_asu=False):
  n_redundancies = 0
  u_grid=[]
  for i in xrange(n[0]):
     b = []
     for j in xrange(n[1]):
        c = []
        for k in xrange(n[2]):
           c.append(0)
        b.append(c)
     u_grid.append(b)
  r_grid = []
  colored_grid_points = []
  for i in xrange(-n[0]//2, n[0]+1):
    b = []
    for j in xrange(-n[1]//2, n[1]+1):
      c = []
      for k in xrange(-n[2]//2, n[2]+1):
        frac = rational.vector((i,j,k), n)
        f = asu.is_inside(frac)
        fv = asu.is_inside(frac, volume_only=True)
        if (len(asu.in_which_cuts(frac)) != 0 and fv):
          colored_grid_points.append(colored_grid_point(
            frac,
            jv_asu.select_color(f)))
        if (volume):
          if (not fv): assert not f
        else:
          fv = 0
        if (f or fv):
          i_pr = i % n[0]
          j_pr = j % n[1]
          k_pr = k % n[2]
          if (u_grid[i_pr][j_pr][k_pr] != 0):
            n_redundancies += 1
            if (not is_stripped_asu):
              print "Redundancy at" , (i,j,k), (i_pr,j_pr,k_pr)
          if (f):
            u_grid[i_pr][j_pr][k_pr] = 1
            c.append(1)
          else:
            u_grid[i_pr][j_pr][k_pr] = 2
            c.append(2)
        else:
          c.append(0)
      b.append(c)
    r_grid.append(b)
  return u_grid, r_grid, colored_grid_points, n_redundancies

def check_compatibility_with_sampling_grid(asu):
  print "Volume vertices:"
  n_outside_sampling_grid = 0
  for vertex in asu.volume_vertices():
    s = ""
    for v in vertex:
      if (v < -rational.int(1,2) or v > 1):
        s = " outside sampling grid"
        n_outside_sampling_grid += 1
        break
    print "  %s%s" % (str(vertex), s)
  assert n_outside_sampling_grid == 0

def check_asu(space_group_number, asu, n, is_stripped_asu, soft_mode):
  sg_info = sgtbx.space_group_info("Hall: " + asu.hall_symbol)
  sg_info.show_summary()
  assert sg_info.type().number() == space_group_number
  print "Gridding:", n
  ops = sg_info.group()
  check_compatibility_with_sampling_grid(asu=asu)
  sys.stdout.flush()
  u_grid, r_grid, colored_grid_points, sampling_n_redundancies = sample_asu(
    asu, n, is_stripped_asu=is_stripped_asu)
  n_redundancies = 0
  redundancies = {}
  for i in xrange(n[0]):
    for j in xrange(n[1]):
      for k in xrange(n[2]):
        n_redundancies += grid_asu(
          ops=ops, n=n, u_grid=u_grid, r_grid=r_grid, i=i,j=j,k=k,
          sampling_n_redundancies=sampling_n_redundancies,
          redundancies=redundancies,
          soft_mode=soft_mode)
  print "number of redundancies: %d+%d," % (
    sampling_n_redundancies, n_redundancies),
  sg_info.show_summary()
  sys.stdout.flush()
  redundancies = sort_redundancies(redundancies)
  recolor_grid_points(
    n, colored_grid_points, redundancies, not is_stripped_asu)
  jv_asu.asu_as_jvx(
    space_group_number, asu, colored_grid_points=colored_grid_points)
  if (not is_stripped_asu):
    analyze_redundancies(asu, n, redundancies)
    if (not soft_mode):
      assert sampling_n_redundancies == 0
      assert n_redundancies == 0
  sys.stdout.flush()

class color_server(object):

  def __init__(self):
    self.color_pairs = (
      ((0,0,255), (153,204,255)),
      ((255,255,0), (255,153,0)),
      ((255,0,255), (255,102,153)),
      ((51,51,51), (178,178,178)))
    self.i = -1

  def next(self):
    if (self.i < len(self.color_pairs)-1):
      self.i += 1
    return self.color_pairs[self.i]

def recolor_grid_points(gridding, colored_grid_points, redundancies, verbose):
  color_srv = color_server()
  processed_points = {}
  for symop,pairs in redundancies:
    if (verbose):
      print "Coloring %d redundancies:" % len(pairs), symop
    sys.stdout.flush()
    colored_point_dict = {}
    for colored_point in colored_grid_points:
      colored_point_dict[colored_point.site] = colored_point
    colors = color_srv.next()
    for pair in pairs:
      for point,color in zip(pair, colors):
        frac = tuple(rational.vector(point, gridding))
        if (not frac in processed_points):
          processed_points[frac] = 1
          colored_point_dict[frac].color = color

def iround(x):
  if (x < 0):
    return int(x-0.5)
  return int(x+0.5)

def rt_plus_unit_shifts(rt, unit_shifts):
  return sgtbx.rt_mx(rt.r(), rt.t().plus(sgtbx.tr_vec(unit_shifts, 1)))

def rt_times_grid_point(rt, i_grid, n):
  grid_point = matrix.col([i_grid[i]/float(n[i]) for i in xrange(3)])
  rotat = matrix.sqr(rt.r().as_double())
  trans = matrix.col(rt.t().as_double())
  eq_pt = rotat*grid_point+trans
  eq_gpt = [0,0,0]
  unit_shifts = [0,0,0]
  for i in xrange(3):
    eg = iround(eq_pt.elems[i]*n[i])
    eq_gpt[i] = eg % n[i]
    u = (eq_gpt[i] - eg) / n[i]
    unit_shifts[i] = iround(u)
    assert abs(u - unit_shifts[i]) < 1.e-5
  return tuple(eq_gpt), rt_plus_unit_shifts(rt, unit_shifts)

def u_index_as_r_index(n, u_index, r_grid, allow_ambiguity):
  r_index = None
  for ui in (0,-1,1):
    ri = u_index[0] + ui * n[0]
    qi = ri + n[0]//2
    if (qi < 0 or qi >= n[0]//2+n[0]+1): continue
    for uj in (0,-1,1):
      rj = u_index[1] + uj * n[1]
      qj = rj + n[1]//2
      if (qj < 0 or qj >= n[1]//2+n[1]+1): continue
      for uk in (0,-1,1):
        rk = u_index[2] + uk * n[2]
        qk = rk + n[2]//2
        if (qk < 0 or qk >= n[2]//2+n[2]+1): continue
        if (r_grid[qi][qj][qk] != 0):
          if (r_index is None):
            r_index, unit_shifts = (ri,rj,rk), (ui,uj,uk)
          else:
            assert allow_ambiguity, "Double redundancy."
  assert r_index is not None
  return r_index, unit_shifts

def grid_asu(
      ops,
      n,
      u_grid,
      r_grid,
      i,j,k,
      sampling_n_redundancies,
      redundancies,
      soft_mode):
  result = 0
  marker = 0
  for rt in ops:
    eq_gpt, rtu = rt_times_grid_point(rt, (i,j,k), n)
    #assert str(rt_times_grid_point(rtu, (i,j,k), n)[1]) == str(rtu)
    if (u_grid[i][j][k] != 0):
      marker = 1
      if (eq_gpt[0] != i or eq_gpt[1] != j or eq_gpt[2] != k):
        if (u_grid[eq_gpt[0]][eq_gpt[1]][eq_gpt[2]] != 0):
          r_pivot, us_pivot = u_index_as_r_index(n, (i,j,k), r_grid,
            allow_ambiguity=sampling_n_redundancies!=0)
          u_eq, rtu = rt_times_grid_point(rt, r_pivot, n)
          r_eq, us_eq = u_index_as_r_index(n, u_eq, r_grid,
            allow_ambiguity=sampling_n_redundancies!=0)
          rtuu = rt_plus_unit_shifts(rtu, us_eq)
          s = str(rtuu)
          v = r_pivot, r_eq
          #print "Redundancy at", v, s
          result += 1
          try: redundancies[s].append(v)
          except KeyboardInterrupt: raise
          except: redundancies[s] = [v]
    else:
      if (u_grid[eq_gpt[0]][eq_gpt[1]][eq_gpt[2]] != 0):
        marker = 1
        break
  if (marker != 1):
    print "Orbit does not intersect with asymmetric unit", (i,j,k)
    assert soft_mode
  return result

def compare_redundancies(a, b):
  return cmp(len(b[1]), len(a[1]))

def sort_redundancies(redundancies):
  redundancies = redundancies.items()
  redundancies.sort(compare_redundancies)
  return redundancies

def str_ev(ev):
  return "[%d,%d,%d]" % ev

def slice(pairs, i):
  result = []
  for pair in pairs:
    result.append(pair[i])
  return result

def rt_mx_analysis(s):
  r_info = sgtbx.rot_mx_info(s.r())
  t_info = sgtbx.translation_part_info(s)
  t_intrinsic = str(t_info.intrinsic_part().mod_positive())
  t_shift = str(t_info.origin_shift().mod_positive())
  if (r_info.type() == 1):
    return ("1", "-", "-", "-")
  if (r_info.type() == -1):
    return (str(r_info.type()), "-", "-", "location=(%s)" % (t_shift,))
  if (abs(r_info.type()) == 2):
    return (str(r_info.type()),
            "axis="+str_ev(r_info.ev()),
            "(%s)" % (t_intrinsic,),
            "location=(%s)" % (t_shift,))
  sense = "+"
  if (r_info.sense() < 0):
    sense = "-"
  return (str(r_info.type())+sense,
          "axis="+str_ev(r_info.ev()),
          "(%s)" % (t_intrinsic,),
          "location=(%s)" % (t_shift,))

def analyze_redundancies(asu, n, redundancies, verbose=1):
  if (len(redundancies) == 0): return
  print "Overview:"
  for symop, pairs in redundancies:
    print symop, ": number of redundancies:", len(pairs)
    print "  ", rt_mx_analysis(sgtbx.rt_mx(symop))
  print "Details:"
  for symop, pairs in redundancies:
    print symop, ": number of redundancies:", len(pairs)
    print "  ", rt_mx_analysis(sgtbx.rt_mx(symop))
    all_cuts = dicts.with_default_factory(dict)
    not_in_cuts = {}
    for pair in pairs:
      for point in pair:
        cuts = asu.in_which_cuts(rational.vector(point, n))
        if (len(cuts) == 0):
          not_in_cuts[point] = 1
        all_cuts[tuple(cuts)][point] = 1
    print "    In cuts:"
    for cuts,points in all_cuts.items():
      print "     ",
      show_amp = False
      for cut in cuts:
        if (show_amp): print "&",
        print cut,
        show_amp = True
      print "#points: %d:" % len(points),
      print str(points.keys()[:4]).replace(" ", "")
    if (verbose):
      print "    Pairs:"
      for pair in pairs:
        print "      ", pair
    if (len(not_in_cuts) > 0):
      print "    Not in cuts:"
      for point in not_in_cuts.keys():
        print "     ", point
      raise AssertionError, "Some redundant points not in any cuts."
    print

def check_multiplicities(asu, n):
  space_group = sgtbx.space_group(asu.hall_symbol)
  all_cuts = asu.extract_all_cuts()
  print "Total number of cuts:", len(all_cuts)
  def get_code(point):
    result = 0
    bit = 1
    for cut in all_cuts:
      if (cut.evaluate(point) == 0):
        result += bit
      bit *= 2
    return result
  mults_by_code = {}
  for i in xrange(-n[0]//2, n[0]+1):
    for j in xrange(-n[1]//2, n[1]+1):
      for k in xrange(-n[2]//2, n[2]+1):
        point = rational.vector((i,j,k), n)
        if (asu.is_inside(point)):
          code = get_code(point)
          if (code != 0):
            m = space_group.multiplicity(site=point)
            mults_by_code.setdefault(code, set()).add(m)
  for code,mults in mults_by_code.items():
    if (len(mults) != 1):
      print "PROBLEM:", space_group.type().number(), mults_by_code
      break
  else:
    print "cut intersection multiplicities unique:"
    order_z = space_group.order_z()
    tab_codes = []
    for code in sorted(mults_by_code.keys()):
      m = list(mults_by_code[code])[0]
      if (m != order_z):
        print code, m
        tab_codes.append((code, m))
    print "Number of cut intersection codes:", len(tab_codes)

def test_all(n):
  for space_group_number in xrange(1, 231):
    cmd = "cctbx.python %s" % sys.argv[0] \
        + " %d,%d,%d " % n +str(space_group_number)
    print cmd
    sys.stdout.flush()
    easy_run.call(cmd)
    print
    sys.stdout.flush()

if (__name__=="__main__"):
  flags = command_line.parse_options(sys.argv[1:], [
    "show_asu",
    "strip",
    "strip_grid",
    "strip_polygons",
    "enantiomorphic",
    "soft",
    "multiplicities",
  ])
  assert len(flags.regular_args) > 0
  gridding = flags.regular_args[0].split(",")
  assert len(gridding) in (1,3)
  gridding = tuple([int(n) for n in gridding])
  if (len(gridding) == 1): gridding *= 3
  if (not os.path.isdir("asu_gallery")):
    os.mkdir("asu_gallery")
  if (len(flags.regular_args) == 1):
    if (not flags.enantiomorphic):
      test_all(gridding)
    else:
      flags.regular_args.extend([str(i) for i in
       (76, 78, 91, 95, 92, 96,
        144, 145, 151, 153, 152, 154,
        169, 170, 171, 172, 178, 179, 180, 181,
        212, 213)])
  if (len(flags.regular_args) > 1):
    for arg in flags.regular_args[1:]:
      numbers = [int(n) for n in arg.split('-')]
      assert len(numbers) in (1,2)
      if (len(numbers) == 1): numbers *= 2
      for space_group_number in xrange(numbers[0], numbers[1]+1):
        asu_original = reference_table.get_asu(space_group_number)
        assert sgtbx.space_group(asu_original.hall_symbol) \
            == sgtbx.space_group_info(number=space_group_number).group()
        asu = asu_original
        if (flags.strip or flags.strip_polygons):
          asu = asu_original.volume_only()
        print "Writing asu_gallery files"
        jv_asu.asu_as_jvx(space_group_number, asu)
        if (flags.strip_grid):
          asu = asu_original.volume_only()
        if (flags.show_asu):
          asu.show_comprehensive_summary()
        check_asu(
          space_group_number=int(space_group_number),
          asu=asu,
          n=gridding,
          is_stripped_asu=(flags.strip or flags.strip_grid),
          soft_mode=flags.soft)
        if (flags.multiplicities):
          check_multiplicities(
            asu=asu_original,
            n=gridding)
