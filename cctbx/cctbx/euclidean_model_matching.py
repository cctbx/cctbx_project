import math
from cctbx.misc import python_utils
from cctbx import matrix
from cctbx import xutils
from cctbx_boost import uctbx
from cctbx_boost import sgtbx

def filter_shift(continuous_shift_flags, shift, selector=1):
  filtered_shift = [0,0,0]
  for i in xrange(3):
    if (continuous_shift_flags[i] == selector):
      filtered_shift[i] = shift[i]
  return filtered_shift

class euclidean_match_symmetry:

  def __init__(self, crystal_space_group_info, use_K2L=1, use_L2N=0):
    python_utils.adopt_init_args(self, locals())
    self.rtmx = sgtbx.SpaceGroup()
    self.ss = sgtbx.StructureSeminvariant(crystal_space_group_info.SgOps())
    self.addl_gen = \
      self.crystal_space_group_info.getAddlGeneratorsOfEuclideanNormalizer(
        use_K2L, use_L2N)
    for g in self.addl_gen:
      # add additional generators to self.rtmx
      self.rtmx.expandSMx(g)
    self.continuous_shifts = []
    for i in xrange(self.ss.size()):
      if (self.ss.M(i) == 0):
        # collect continous allowed origin shifts
        self.continuous_shifts.append(self.ss.V(i))
      else:
        # add discrete allowed origin shifts to self.rtmx
        self.rtmx.expandLTr(self.ss.V(i), self.ss.M(i))

  def continuous_shifts_are_principal(self):
    for pa in self.continuous_shifts:
      if (not pa in ((1,0,0),(0,1,0),(0,0,1))): return False
    return True

  def filter_shift(self, shift, selector=1):
    return filter_shift(self.continuous_shift_flags, shift, selector)

  def set_continuous_shift_flags(self):
    assert self.continuous_shifts_are_principal()
    self.continuous_shift_flags = [0,0,0]
    for pa in self.continuous_shifts:
      for i in xrange(3):
        if (pa[i]): self.continuous_shift_flags[i] = 1

  def show(self, title=""):
    print "euclidean_match_symmetry:", title
    print self.rtmx.Info().BuildLookupSymbol()
    print self.continuous_shifts
    print

class labeled_position:

  def __init__(self, label, coordinates):
    python_utils.adopt_init_args(self, locals())

  def __repr__(self):
    return "%-4s %7.4f %7.4f %7.4f" % ((self.label,) + tuple(self.coordinates))

class model(xutils.crystal_symmetry):

  def __init__(self, xsym, labeled_positions):
    xutils.crystal_symmetry.__init__(self, xsym.UnitCell, xsym.SgInfo)
    self.labeled_positions = list(labeled_positions)

  def transform_to_reference_setting(self):
    CBOp = self.SgInfo.CBOp()
    cb_unit_cell = CBOp.apply(self.UnitCell)
    ref_space_group = sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols(
      self.SgInfo.SgNumber()))
    ref_space_group.CheckUnitCell(cb_unit_cell)
    cb_positions = []
    for pos in self.labeled_positions:
      cb_positions.append(labeled_position(
        pos.label,
        CBOp(pos.coordinates)))
    return model(
      xutils.crystal_symmetry(cb_unit_cell, ref_space_group.Info()),
      cb_positions)

  def __len__(self): return len(self.labeled_positions)
  def    size(self): return len(self.labeled_positions)

  def __getitem__(self, key):
    return self.labeled_positions[key]

  def show(self, title):
    print title
    print "Unit cell:", self.UnitCell
    print "Space group:", self.SgInfo.BuildLookupSymbol()
    for lp in self.labeled_positions: print lp
    print

def generate_singles(n, i):
  singles = range(n)
  del singles[i]
  return singles

def pair_sort_function(pair_a, pair_b):
  return cmp(pair_a[0], pair_b[0])

class match_refine:

  def __init__(self, tolerance,
               ref_model1, ref_model2,
               unit_cell, match_symmetry,
               equiv1,
               i_pivot1, i_pivot2,
               eucl_symop,
               unallowed_shift, initial_shift):
    python_utils.adopt_init_args(self, locals())
    self.singles1 = generate_singles(self.ref_model1.size(), self.i_pivot1)
    self.singles2 = generate_singles(self.ref_model2.size(), self.i_pivot2)
    self.pairs = [(self.i_pivot1, self.i_pivot2)]
    self.adjusted_shift = initial_shift[:]
    self.add_pairs()
    self.eliminate_weak_pairs()
    self.ref_eucl_rt = matrix.rt(
      self.eucl_symop.as_tuple()) - matrix.col(self.adjusted_shift)
    self.pairs.sort(pair_sort_function)
    self.calculate_rms()
    # clean up
    del self.unallowed_shift
    del self.initial_shift

  def add_pairs(self):
    while (len(self.singles1) and len(self.singles2)):
      shortest_dist = 2 * self.tolerance
      new_pair = 0
      for is2 in self.singles2:
        c2 = self.eucl_symop.multiply(self.ref_model2[is2].coordinates)
        c2 = python_utils.list_minus(c2, self.adjusted_shift)
        for is1 in self.singles1:
          diff = self.equiv1[is1].getShortestDifference(self.unit_cell, c2)
          dist = self.unit_cell.Length(diff)
          if (dist < shortest_dist):
            # ensure that this pair can be matched within 1 * tolerance.
            diff_allowed = self.match_symmetry.filter_shift(diff, selector=1)
            dist_allowed = self.unit_cell.Length(diff_allowed)
            if (dist_allowed < self.tolerance):
              shortest_dist = dist
              new_pair = (is1, is2)
      if (new_pair == 0):
        break
      self.pairs.append(new_pair)
      self.singles1.remove(new_pair[0])
      self.singles2.remove(new_pair[1])
      self.refine_adjusted_shift()

  def eliminate_weak_pairs(self):
    while 1:
      weak_pair = 0
      max_dist = 0
      for pair in self.pairs[1:]:
        dist = self.calculate_shortest_dist(pair)
        if (dist > max_dist):
          weak_pair = pair
          max_dist = dist
      if (weak_pair == 0): break
      if (max_dist < self.tolerance):
        dist = self.calculate_shortest_dist(self.pairs[0])
        if (dist < self.tolerance):
          break
      assert len(self.pairs) > 1
      self.pairs.remove(weak_pair)
      self.singles1.append(weak_pair[0])
      self.singles2.append(weak_pair[1])
      self.refine_adjusted_shift()

  def apply_eucl_ops(self, i_model2):
    c2 = self.eucl_symop.multiply(self.ref_model2[i_model2].coordinates)
    return python_utils.list_minus(c2, self.adjusted_shift)

  def calculate_shortest_diff(self, pair):
    c2 = self.apply_eucl_ops(pair[1])
    return self.equiv1[pair[0]].getShortestDifference(self.unit_cell, c2)

  def calculate_shortest_dist(self, pair):
    return self.unit_cell.Length(self.calculate_shortest_diff(pair))

  def calculate_shortest_diffs(self):
    shortest_diffs = []
    for pair in self.pairs:
      shortest_diffs.append(self.calculate_shortest_diff(pair))
    return shortest_diffs

  def refine_adjusted_shift(self):
    sum_diff_cart = [0,0,0]
    for diff in self.calculate_shortest_diffs():
      diff_allowed = self.match_symmetry.filter_shift(diff, selector=1)
      diff_cart = self.unit_cell.orthogonalize(diff_allowed)
      sum_diff_cart = python_utils.list_plus(sum_diff_cart, diff_cart)
    mean_diff_cart = [s / len(self.pairs) for s in sum_diff_cart]
    mean_diff_frac = self.unit_cell.fractionalize(mean_diff_cart)
    self.adjusted_shift = python_utils.list_plus(
      self.adjusted_shift, mean_diff_frac)

  def calculate_rms(self):
    sum_dist2 = 0
    for pair in self.pairs:
      sum_dist2 += self.unit_cell.Length2(self.calculate_shortest_diff(pair))
    self.rms = math.sqrt(sum_dist2 / len(self.pairs))

  def show(self):
    print "Match summary:"
    print "  Operator:"
    print "       rotation:", self.rt.r.mathematica_form()
    print "    translation:", self.rt.t.elems
    print "rms coordinate error: %.2f" % (self.rms,)
    print "  Pairs:", len(self.pairs)
    for pair in self.pairs:
      print "   ", self.ref_model1[pair[0]].label,
      print self.ref_model2[pair[1]].label,
      print "%.3f" % (self.calculate_shortest_dist(pair),)
    print "Singles model 1:", len(self.singles1)
    for s in self.singles1:
      print " ", self.ref_model1[s].label,
    print
    print "Singles model 2:", len(self.singles2)
    for s in self.singles2:
      print " ", self.ref_model2[s].label,
    print
    print

def match_sort_function(match_a, match_b):
  i = -cmp(len(match_a.pairs), len(match_b.pairs))
  if (i): return i
  return cmp(match_a.rms, match_b.rms)

def weed_refined_matches(space_group_number, refined_matches,
                         rms_penalty_per_site):
  n_matches = len(refined_matches)
  if (n_matches == 0): return
  best_rms = refined_matches[0].rms
  best_n_pairs = len(refined_matches[0].pairs)
  is_redundant = [0] * n_matches
  for i in xrange(n_matches-1):
    match_i = refined_matches[i]
    if (is_redundant[i]): continue
    if (match_i.rms < best_rms):
      best_rms = match_i.rms
      best_n_pairs = len(match_i.pairs)
    for j in xrange(i+1, n_matches):
      match_j = refined_matches[j]
      if (   match_i.pairs == match_j.pairs
          or (    rms_penalty_per_site
              and match_j.rms > best_rms * (1 - rms_penalty_per_site * (
                    best_n_pairs - len(match_j.pairs))))):
        is_redundant[j] = 1
  for i in xrange(n_matches-1, -1, -1):
    if (is_redundant[i]):
      del refined_matches[i]
  if (space_group_number == 1 and n_matches > 0):
    trivial_matches_only = True
    for match in refined_matches:
      if (len(match.pairs) > 1):
        trivial_matches_only = False
        break
    if (trivial_matches_only):
      while (refined_matches[0].pairs[0] != (0,0)): del refined_matches[0]
      while (len(refined_matches) > 1): del refined_matches[-1]
    else:
      while (len(refined_matches[-1].pairs) == 1): del refined_matches[-1]

def match_rt_from_ref_eucl_rt(model1_cbop, model2_cbop, ref_eucl_rt):
  inv_m1 = matrix.rt(model1_cbop.InvM().as_tuple())
  m2 = matrix.rt(model2_cbop.M().as_tuple())
  # X2_orig -> X2_ref -> X1_ref -> X1_orig
  #         m2   ref_eucl_rt   inv_m1
  # X1 = inv_m1 * ref_eucl_rt * m2 * X2
  return inv_m1 * ref_eucl_rt * m2

def match_models(model1, model2,
                 tolerance=1., models_are_diffraction_index_equivalent=0,
                 rms_penalty_per_site = 0.05):
  ref_model1 = model1.transform_to_reference_setting()
  ref_model2 = model2.transform_to_reference_setting()
  assert xutils.are_similar_unit_cells(ref_model1.UnitCell,ref_model2.UnitCell)
  assert ref_model1.SgOps == ref_model2.SgOps
  unit_cell = ref_model1.UnitCell
  match_symmetry = euclidean_match_symmetry(
    ref_model1.SgInfo,
    use_K2L=1, use_L2N=(not models_are_diffraction_index_equivalent))
  match_symmetry.set_continuous_shift_flags()
  equiv1 = []
  for pos in ref_model1:
    equiv1.append(sgtbx.SymEquivCoordinates(ref_model1.SgOps, pos.coordinates))
  refined_matches = []
  for i_pivot1 in xrange(ref_model1.size()):
    for i_pivot2 in xrange(ref_model2.size()):
      for eucl_symop in match_symmetry.rtmx:
        c2 = eucl_symop.multiply(ref_model2[i_pivot2].coordinates)
        diff = equiv1[i_pivot1].getShortestDifferenceUnderAllowedOriginShifts(
          unit_cell, c2, match_symmetry.continuous_shift_flags)
        unallowed_shift = match_symmetry.filter_shift(diff, selector=0)
        if (unit_cell.Length(unallowed_shift) < tolerance):
          allowed_shift = match_symmetry.filter_shift(diff, selector=1)
          match = match_refine(tolerance,
                               ref_model1, ref_model2,
                               unit_cell, match_symmetry,
                               equiv1,
                               i_pivot1, i_pivot2,
                               eucl_symop,
                               unallowed_shift, allowed_shift)
          match.rt = match_rt_from_ref_eucl_rt(
            model1.SgInfo.CBOp(),
            model2.SgInfo.CBOp(),
            match.ref_eucl_rt)
          refined_matches.append(match)
  refined_matches.sort(match_sort_function)
  weed_refined_matches(model1.SgInfo.SgNumber(), refined_matches,
                       rms_penalty_per_site)
  return refined_matches

###############################################
# The code below this line is just for testing.
###############################################

def debug_analyze_refined_matches(model1, model2, refined_matches):
  solution_counter = 0
  for match in refined_matches:
    match.show()
    debug_verify_match(match.ref_model1, match.ref_model2, match.tolerance,
                       match.ref_eucl_rt, match.pairs)
    debug_verify_match(model1, model2, match.tolerance,
                       match.rt, match.pairs)
    if (    debug_analyze_singles(model1, match.singles1)
        and debug_analyze_singles(model2, match.singles2)):
      solution_counter += 1
  print "total matches:", len(refined_matches)
  print "solutions:", solution_counter
  assert solution_counter != 0
  print

def debug_analyze_singles(model, singles):
  for i in singles:
    if (model[i].label.startswith("S")): return False
  return True

def debug_verify_match(model1, model2, tolerance, match_rt, pairs):
  adj_tolerance = tolerance * (1 + 1.e-6)
  for pair in pairs:
    c1 = model1[pair[0]].coordinates
    c2 = match_rt * model2[pair[1]].coordinates
    equiv_c2 = sgtbx.SymEquivCoordinates(model1.SgOps, c2.elems)
    diff = equiv_c2.getShortestDifference(model1.UnitCell, c1)
    assert model1.UnitCell.Length(diff) < adj_tolerance, \
      model1.SgInfo.BuildLookupSymbol()

class test_model(model):

  def __init__(self, model_id = "SBT", n_elements = 4):
    if (model_id == None): return
    self.model_id = model_id
    lp = labeled_position
    if (type(model_id) == type("")):
      if (model_id == "SBT"):
        unit_cell = uctbx.UnitCell(
          (16.8986, 16.8986, 16.8986, 61.1483, 61.1483, 61.1483))
        space_group_info = sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols(
          "R -3 m :R")).Info()
        labeled_positions = (
          lp("SI1", (-0.3584, 0.2844, 0.4622)),
          lp("SI2", (-0.2133, 0.9659, -0.6653)),
          lp("SI3", (-0.8358, 0.7, 0.3431)),
          lp("SI4", (0.4799, 1.836, 0.6598)))
      else:
        raise RuntimeError, "Unknown model_id: " + model_id
      model.__init__(self,
        xutils.crystal_symmetry(unit_cell, space_group_info),
        labeled_positions)
    else:
      from cctbx.development import debug_utils
      elements = ["S"] * n_elements
      xtal = debug_utils.random_structure(model_id, elements,
        volume_per_atom=50.,
        min_distance = 2.0,
        general_positions_only = 0)
      labeled_positions = []
      for site in xtal:
        labeled_positions.append(labeled_position(
          site.Label(), site.Coordinates()))
      model.__init__(self, xtal, labeled_positions)

  def create_new_test_model(self, new_labeled_positions):
    new_test_model = test_model(None)
    model.__init__(new_test_model,
      xutils.crystal_symmetry(self.UnitCell, self.SgInfo),
      new_labeled_positions)
    return new_test_model

  def shuffle_positions(self):
    from cctbx.development import debug_utils
    shuffled_positions = list(self.labeled_positions)
    debug_utils.random.shuffle(shuffled_positions)
    return self.create_new_test_model(shuffled_positions)

  def random_symmetry_mates(self):
    from cctbx.development import debug_utils
    new_labeled_positions = []
    for lp in self.labeled_positions:
      equiv_coor = sgtbx.SymEquivCoordinates(self.SgOps, lp.coordinates)
      i = debug_utils.random.randrange(equiv_coor.M())
      new_labeled_positions.append(labeled_position(
        lp.label, equiv_coor(i)))
    return self.create_new_test_model(new_labeled_positions)

  def apply_random_eucl_op(self, models_are_diffraction_index_equivalent = 0):
    from cctbx.development import debug_utils
    match_symmetry = euclidean_match_symmetry(
      self.SgInfo,
      use_K2L=1, use_L2N=(not models_are_diffraction_index_equivalent))
    match_symmetry.set_continuous_shift_flags()
    i = debug_utils.random.randrange(match_symmetry.rtmx.OrderZ())
    eucl_symop = match_symmetry.rtmx(i)
    shift = [debug_utils.random.random() for i in xrange(3)]
    allowed_shift = match_symmetry.filter_shift(shift, selector=1)
    new_labeled_positions = []
    for lp in self.labeled_positions:
      new_coor = eucl_symop.multiply(lp.coordinates)
      new_coor = python_utils.list_plus(new_coor, allowed_shift)
      new_labeled_positions.append(labeled_position(
        lp.label, new_coor))
    return self.create_new_test_model(new_labeled_positions)

  def add_random_positions(self, number_of_new_positions = 3, label = "R",
                           min_distance = 1.0):
    from cctbx.development import debug_utils
    existing_positions = []
    for lp in self.labeled_positions:
      existing_positions.append(lp.coordinates)
    new_positions = debug_utils.generate_positions(
      number_of_new_positions,
      self,
      sgtbx.SpecialPositionSnapParameters(
        self.UnitCell, self.SgOps, 1, min_distance),
      min_hetero_distance = min_distance,
      general_positions_only = 0,
      existing_positions = existing_positions)
    new_labeled_positions = []
    i = 0
    for coor in new_positions:
      i += 1
      new_labeled_positions.append(labeled_position(
        "%s%d" % (label, i), coor))
    return self.create_new_test_model(
      list(self.labeled_positions) + new_labeled_positions)

  def shake_positions(self, sigma = 0.2, min_distance = 1.0):
    from cctbx.development import debug_utils
    snap_parameters = sgtbx.SpecialPositionSnapParameters(
      self.UnitCell, self.SgOps, 1, min_distance)
    new_labeled_positions = []
    for lp in self.labeled_positions:
      new_coor = debug_utils.shake_position(
        self, snap_parameters, lp.coordinates,
        sigma, max_diff = min_distance * 0.99)
      new_labeled_positions.append(labeled_position(
        lp.label, new_coor))
    return self.create_new_test_model(new_labeled_positions)

def run_test(argv):
  from cctbx.development import debug_utils
  Flags = debug_utils.command_line_options(argv, (
    "RandomSeed",
    "AllSpaceGroups",
    "StaticModels",
  ))
  if (not Flags.RandomSeed): debug_utils.set_random_seed(0)
  symbols_to_stdout = 0
  if (len(argv) > Flags.n):
    symbols = argv
  else:
    symbols = debug_utils.get_test_space_group_symbols(Flags.AllSpaceGroups)
    symbols_to_stdout = 1
  if (Flags.StaticModels):
    model1 = (test_model()
      .add_random_positions(2, "A")
      .shuffle_positions()
      .random_symmetry_mates()
      .apply_random_eucl_op())
    model2 = (test_model()
      .add_random_positions(3, "B")
      .shuffle_positions()
      .random_symmetry_mates()
      .apply_random_eucl_op()
      .shake_positions())
    for i in xrange(2):
      m1 = model1
      if (i): m1 = model1.transform_to_reference_setting()
      for j in xrange(2):
        m2 = model2
        if (j): m2 = model2.transform_to_reference_setting()
        m1.show("Model1(%d)" % (i,))
        m2.show("Model2(%d)" % (j,))
        refined_matches = match_models(m1, m2, rms_penalty_per_site=0)
        debug_analyze_refined_matches(m1, m2, refined_matches)
  else:
    for RawSgSymbol in symbols:
      if (RawSgSymbol.startswith("--")): continue
      SgSymbols = sgtbx.SpaceGroupSymbols(RawSgSymbol)
      SgInfo = sgtbx.SpaceGroup(SgSymbols).Info()
      LookupSymbol = SgInfo.BuildLookupSymbol()
      sys.stdout.flush()
      print >> sys.stderr, LookupSymbol
      sys.stderr.flush()
      if (symbols_to_stdout):
        print LookupSymbol
        sys.stdout.flush()
      model_core = test_model(SgInfo)
      model1 = (model_core
        .add_random_positions(2, "A")
        )
      model2 = (model_core
        .add_random_positions(3, "B")
        .shuffle_positions()
        .random_symmetry_mates()
        .apply_random_eucl_op()
        .shake_positions()
        )
      model_core.show("Core")
      model1.show("Model1")
      model2.show("Model2")
      refined_matches = match_models(model1, model2, rms_penalty_per_site=0)
      debug_analyze_refined_matches(model1, model2, refined_matches)

if (__name__ == "__main__"):
  import sys, os
  run_test(sys.argv[1:])
  t = os.times()
  print "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])
