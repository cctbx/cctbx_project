import math
from cctbx.misc import python_utils
from cctbx import xutils
from cctbx_boost import uctbx
from cctbx_boost import sgtbx

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

  def set_continuous_shift_flags(self):
    assert self.continuous_shifts_are_principal()
    self.continuous_shift_flags = [0,0,0]
    for pa in self.continuous_shifts:
      for i in xrange(3):
        if (pa[i]): self.continuous_shift_flags[i] = 1

  def filter_shift(self, shift, selector=1):
    filtered_shift = [0,0,0]
    for i in xrange(3):
      if (self.continuous_shift_flags[i] == selector):
        filtered_shift[i] = shift[i]
    return filtered_shift

  def show(self, title=""):
    print "euclidean_match_symmetry:", title
    print self.rtmx.Info().BuildLookupSymbol()
    print self.continuous_shifts
    print

class labelled_position:

  def __init__(self, label, coordinates):
    python_utils.adopt_init_args(self, locals())

  def __repr__(self):
    return "%-4s %7.4f %7.4f %7.4f" % ((self.label,) + tuple(self.coordinates))

class model(xutils.crystal_symmetry):

  def __init__(self, xsym, labelled_positions):
    xutils.crystal_symmetry.__init__(self, xsym.UnitCell, xsym.SgInfo)
    self.labelled_positions = list(labelled_positions)

  def transform_to_reference_setting(self):
    CBOp = self.SgInfo.CBOp()
    cb_unit_cell = CBOp.apply(self.UnitCell)
    ref_space_group = sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols(
      self.SgInfo.SgNumber()))
    ref_space_group.CheckUnitCell(cb_unit_cell)
    cb_positions = []
    for pos in self.labelled_positions:
      cb_positions.append(labelled_position(
        pos.label,
        CBOp(pos.coordinates)))
    return model(
      xutils.crystal_symmetry(cb_unit_cell, ref_space_group.Info()),
      cb_positions)

  def __len__(self): return len(self.labelled_positions)
  def    size(self): return len(self.labelled_positions)

  def __getitem__(self, key):
    return self.labelled_positions[key]

  def show(self, title):
    print title
    print self.UnitCell
    print self.SgInfo.BuildLookupSymbol()
    for lp in self.labelled_positions: print lp
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
    self.add_matches()
    self.eliminate_weak_matches()
    self.pairs.sort(pair_sort_function)

  def debug_fmt_pair(self, i, j): # XXX
    return "%s, %s:" % (
      self.ref_model1[i].label,
      self.ref_model2[j].label)

  def add_matches(self):
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

  def eliminate_weak_matches(self):
    while 1:
      max_dist = 0
      weak_pair = 0
      for diff, pair in zip(self.calculate_shortest_diffs(), self.pairs):
        dist = self.unit_cell.Length(diff)
        if (dist > max_dist):
          max_dist = dist
          weak_pair = pair
      if (max_dist < self.tolerance or weak_pair == 0):
        break
      assert len(self.pairs) > 1
      self.pairs.remove(weak_pair)
      self.singles1.append(weak_pair[0])
      self.singles2.append(weak_pair[1])
      self.refine_adjusted_shift()

  def apply_eucl_ops(self, i_model2):
    c2 = self.eucl_symop.multiply(self.ref_model2[i_model2].coordinates)
    return python_utils.list_minus(c2, self.adjusted_shift)

  def calculate_shortest_diff(self, i_model1, i_model2):
    c2 = self.apply_eucl_ops(i_model2)
    return self.equiv1[i_model1].getShortestDifference(self.unit_cell, c2)

  def calculate_shortest_diffs(self):
    shortest_diffs = []
    for pair in self.pairs:
      shortest_diffs.append(self.calculate_shortest_diff(pair[0], pair[1]))
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

  def rms(self):
    sum_diff_cart = [0,0,0]
    for diff in self.calculate_shortest_diffs():
      diff_cart = self.unit_cell.orthogonalize(diff)
      sum_diff_cart = python_utils.list_plus(sum_diff_cart, diff_cart)
    mean_diff_cart = [s / len(self.pairs) for s in sum_diff_cart]
    return math.sqrt(python_utils.list_dot_product(mean_diff_cart))

  def show(self):
    print "Match summary:"
    print self.eucl_symop.as_xyz(),
    print self.adjusted_shift,
    print "%.2f" % (self.rms(),)
    for pair in self.pairs:
      print self.ref_model1[pair[0]].label,
      print self.ref_model2[pair[1]].label
    print "Singles model 1:", len(self.singles1)
    for s in self.singles1:
      print " ", self.ref_model1[s].label,
    print
    print "Singles model 2:", len(self.singles2)
    for s in self.singles2:
      print " ", self.ref_model2[s].label,
    print
    print

def match_models(model1, model2,
                 tolerance=1., models_are_diffraction_index_equivalent=0,
                 pickle_on_error = 1):
  debug_solution_counter = 0 # XXX
  debug_match_refine_counter = 0 # XXX
  ref_model1 = model1.transform_to_reference_setting()
  ref_model1.show("ref_model1")
  ref_model2 = model2.transform_to_reference_setting()
  ref_model2.show("ref_model2")
  assert ref_model1.SgOps == ref_model2.SgOps
  # XXX unit cells must also be compatible
  unit_cell = ref_model1.UnitCell
  SgInfo = ref_model1.SgInfo
  match_symmetry = euclidean_match_symmetry(
    SgInfo, use_K2L=1, use_L2N=(not models_are_diffraction_index_equivalent))
  match_symmetry.show()
  match_symmetry.set_continuous_shift_flags()
  equiv1 = []
  for pos in ref_model1:
    equiv1.append(sgtbx.SymEquivCoordinates(ref_model1.SgOps, pos.coordinates))
  for i_pivot1 in xrange(ref_model1.size()):
    for i_pivot2 in xrange(ref_model2.size()):
      for eucl_symop in match_symmetry.rtmx:
        c2 = eucl_symop.multiply(ref_model2[i_pivot2].coordinates)
        diff = equiv1[i_pivot1].getShortestDifferenceUnderAllowedOriginShifts(
          unit_cell, c2, match_symmetry.continuous_shift_flags)
        unallowed_shift = match_symmetry.filter_shift(diff, selector=0)
        if (unit_cell.Length(unallowed_shift) < tolerance):
          allowed_shift = match_symmetry.filter_shift(diff, selector=1)
          matches = match_refine(tolerance,
                                 ref_model1, ref_model2,
                                 unit_cell, match_symmetry,
                                 equiv1,
                                 i_pivot1, i_pivot2,
                                 eucl_symop,
                                 unallowed_shift, allowed_shift)
          debug_match_refine_counter += 1
          matches.show() # XXX
          debug_verify_matches(ref_model1, ref_model2, tolerance,
                               eucl_symop,
                               matches.adjusted_shift, matches.pairs)
          if (    debug_analyze_singles(model1, matches.singles1)
              and debug_analyze_singles(model2, matches.singles2)):
            debug_solution_counter += 1
  print "match refine:", debug_match_refine_counter
  print "solutions:", debug_solution_counter
  assert debug_solution_counter != 0
  if (not (debug_solution_counter != 0)):
    print "ERROR: no solution:", SgInfo.BuildLookupSymbol(), SgInfo.SgNumber()
    if (pickle_on_error):
      import cPickle
      file_name = "no_sol_%03d.pickle" % (SgInfo.SgNumber(),)
      f = open(file_name, "w")
      del model1.SgInfo
      del model2.SgInfo
      cPickle.dump((model1, model2), f)
      f.close()

###############################################
# The code below this line is only for testing.
###############################################

def debug_analyze_singles(model, singles):
  for i in singles:
    if (model[i].label.startswith("S")): return False
  return True

def debug_verify_matches(model1, model2, tolerance,
                         eucl_symop, shift, pairs):
  adj_tolerance = tolerance * (1 + 1.e-6)
  for pair in pairs:
    c1 = model1[pair[0]].coordinates
    c2 = model2[pair[1]].coordinates
    c2 = eucl_symop.multiply(c2)
    c2 = python_utils.list_minus(c2, shift)
    equiv_c2 = sgtbx.SymEquivCoordinates(model1.SgOps, c2)
    diff = equiv_c2.getShortestDifference(model1.UnitCell, c1)
    #print diff, model1.UnitCell.Length(diff)
    assert model1.UnitCell.Length(diff) < adj_tolerance, \
      model1.SgInfo.BuildLookupSymbol()

class test_model(model):

  def __init__(self, model_id = "SBT", n_elements = 4):
    if (model_id == None): return
    self.model_id = model_id
    lp = labelled_position
    if (type(model_id) == type("")):
      if (model_id == "SBT"):
        unit_cell = uctbx.UnitCell(
          (16.8986, 16.8986, 16.8986, 61.1483, 61.1483, 61.1483))
        space_group_info = sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols(
          "R -3 m :R")).Info()
        labelled_positions = (
          lp("SI1", (-0.3584, 0.2844, 0.4622)),
          lp("SI2", (-0.2133, 0.9659, -0.6653)),
          lp("SI3", (-0.8358, 0.7, 0.3431)),
          lp("SI4", (0.4799, 1.836, 0.6598)))
      else:
        raise RuntimeError, "Unknown model_id: " + model_id
      model.__init__(self,
        xutils.crystal_symmetry(unit_cell, space_group_info),
        labelled_positions)
    else:
      from cctbx.development import debug_utils
      elements = ["S"] * n_elements
      xtal = debug_utils.random_structure(model_id, elements,
        volume_per_atom=50.,
        min_distance = 2.0,
        general_positions_only = 0)
      labelled_positions = []
      for site in xtal:
        labelled_positions.append(labelled_position(
          site.Label(), site.Coordinates()))
      model.__init__(self, xtal, labelled_positions)

  def create_new_test_model(self, new_labelled_positions):
    new_test_model = test_model(None)
    model.__init__(new_test_model,
      xutils.crystal_symmetry(self.UnitCell, self.SgInfo),
      new_labelled_positions)
    return new_test_model

  def shuffle_positions(self):
    from cctbx.development import debug_utils
    shuffled_positions = list(self.labelled_positions)
    debug_utils.random.shuffle(shuffled_positions)
    return self.create_new_test_model(shuffled_positions)

  def random_symmetry_mates(self):
    from cctbx.development import debug_utils
    new_labelled_positions = []
    for lp in self.labelled_positions:
      equiv_coor = sgtbx.SymEquivCoordinates(self.SgOps, lp.coordinates)
      i = debug_utils.random.randrange(equiv_coor.M())
      new_labelled_positions.append(labelled_position(
        lp.label, equiv_coor(i)))
    return self.create_new_test_model(new_labelled_positions)

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
    new_labelled_positions = []
    for lp in self.labelled_positions:
      new_coor = eucl_symop.multiply(lp.coordinates)
      new_coor = python_utils.list_plus(new_coor, allowed_shift)
      new_labelled_positions.append(labelled_position(
        lp.label, new_coor))
    return self.create_new_test_model(new_labelled_positions)

  def add_random_positions(self, number_of_new_positions = 3, label = "R",
                           min_distance = 1.0):
    from cctbx.development import debug_utils
    existing_positions = []
    for lp in self.labelled_positions:
      existing_positions.append(lp.coordinates)
    new_positions = debug_utils.generate_positions(
      number_of_new_positions,
      self,
      sgtbx.SpecialPositionSnapParameters(
        self.UnitCell, self.SgOps, 1, min_distance),
      min_hetero_distance = min_distance,
      general_positions_only = 0,
      existing_positions = existing_positions)
    new_labelled_positions = []
    i = 0
    for coor in new_positions:
      i += 1
      new_labelled_positions.append(labelled_position(
        "%s%d" % (label, i), coor))
    return self.create_new_test_model(
      list(self.labelled_positions) + new_labelled_positions)

  def shake_positions(self, sigma = 0.2, min_distance = 1.0):
    from cctbx.development import debug_utils
    snap_parameters = sgtbx.SpecialPositionSnapParameters(
      self.UnitCell, self.SgOps, 1, min_distance)
    new_labelled_positions = []
    for lp in self.labelled_positions:
      new_coor = debug_utils.shake_position(
        self, snap_parameters, lp.coordinates,
        sigma, max_diff = min_distance * 0.99)
      new_labelled_positions.append(labelled_position(
        lp.label, new_coor))
    return self.create_new_test_model(new_labelled_positions)

def run_test(argv):
  from cctbx.development import debug_utils
  Flags = debug_utils.command_line_options(argv, (
    "RandomSeed",
    "AllSpaceGroups",
    "StaticModels",
    "load",
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
    model1.show("Model1")
    model2.show("Model2")
    match_models(model1, model2)
  elif (Flags.load):
    import cPickle
    for file_name in argv:
      if (file_name.startswith("--")): continue
      f = open(file_name, "r")
      model1, model2 = cPickle.load(f)
      f.close()
      model1.SgInfo = model1.SgOps.Info()
      model2.SgInfo = model1.SgOps.Info()
      model1.show("Model1")
      model2.show("Model2")
      match_models(model1, model2, pickle_on_error=0)
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
      match_models(model1, model2)

if (__name__ == "__main__"):
  import sys, os
  run_test(sys.argv[1:])
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
