import math
from cctbx.misc import python_utils
from cctbx import xutils
from cctbx_boost import uctbx
from cctbx_boost import sgtbx

class euclidean_match_symmetry:

  def __init__(self, crystal_space_group_info, use_K2L=1, use_L2N=0):
    python_utils.adopt_init_args(self, locals())
    self.rtmx = crystal_space_group_info.SgOps()
    self.ss = sgtbx.StructureSeminvariant(self.rtmx)
    self.addl_gen = \
      self.crystal_space_group_info.getAddlGeneratorsOfEuclideanNormalizer(
        use_K2L, use_L2N)
    for g in self.addl_gen:
      # add additional generators to self.rtmx
      self.rtmx.expandSMx(g)
    self.polar_axes = []
    for i in xrange(self.ss.size()):
      if (self.ss.M(i) == 0):
        # collect continous allowed origin shifts
        self.polar_axes.append(self.ss.V(i))
      else:
        # add discrete allowed origin shifts to self.rtmx
        self.rtmx.expandLTr(self.ss.V(i), self.ss.M(i))

  def polar_axes_are_principal(self):
    for pa in self.polar_axes:
      if (not pa in ((1,0,0),(0,1,0),(0,0,1))): return False
    return True

  def set_sum_polar_axes(self):
    assert self.polar_axes_are_principal()
    self.sum_polar_axes = [0,0,0]
    for pa in self.polar_axes:
      for i in xrange(3): self.sum_polar_axes[i] += pa[i]

  def continuous_shift(self, reference, other):
    # definition: reference = other + shift
    shift = [0,0,0]
    for i in xrange(3):
      shift[i] = self.sum_polar_axes[i] * (reference[i] - other[i])
    return shift

  def show(self, title=""):
    print "euclidean_match_symmetry:", title
    print self.rtmx.Info().BuildLookupSymbol()
    print self.polar_axes
    print

class labelled_position:

  def __init__(self, label, coordinates):
    python_utils.adopt_init_args(self, locals())

class model(xutils.crystal_symmetry):

  def __init__(self, xsym, labelled_positions):
    xutils.crystal_symmetry.__init__(self, xsym.UnitCell, xsym.SgInfo)
    self.labelled_positions = labelled_positions

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
    for lp in self.labelled_positions:
      print lp.label, lp.coordinates
    print

def generate_singles(n, i):
  singles = range(n)
  del singles[i]
  return singles

class match_refine:

  def __init__(self, tolerance,
               ref_model1, ref_model2,
               unit_cell, match_symmetry,
               i_pivot1, i_pivot2,
               symop, continuous_shift):
    python_utils.adopt_init_args(self, locals())
    self.singles1 = generate_singles(self.ref_model1.size(), self.i_pivot1)
    self.singles2 = generate_singles(self.ref_model2.size(), self.i_pivot2)
    self.pairs = [(self.i_pivot1, self.i_pivot2)]
    self.adjusted_shift = continuous_shift[:]
    self.add_matches()
    self.eliminate_weak_matches()

  def add_matches(self):
    while (len(self.singles1) and len(self.singles2)):
      shortest_distance = 2 * self.tolerance
      new_pair = 0
      for is2 in self.singles2:
        c2 = self.symop.multiply(self.ref_model2[is2].coordinates)
        c2 = python_utils.list_plus(c2, self.adjusted_shift)
        for is1 in self.singles1:
          c1 = self.ref_model1[is1].coordinates
          distance = self.unit_cell.Distance(c1, c2)
          if (distance < shortest_distance):
            shortest_distance = distance
            new_pair = (is1, is2)
            self.singles1.remove(is1)
            self.singles2.remove(is2)
      if (new_pair == 0):
        break
      self.pairs.append(new_pair)
      self.refine_adjusted_shift()

  def refine_adjusted_shift(self):
    # XXX It is not clear that this is the correct least squares fit
    # We might need special orthogonalization matrices that ensure
    # that the polar axes are parallel or perpendicular to the
    # cartesian basis system.
    sum_diff_cart = [0,0,0]
    for pair in self.pairs:
      cf1 = self.ref_model1[pair[0]].coordinates
      cf2 = self.ref_model2[pair[1]].coordinates
      diff_allowed = self.match_symmetry.continuous_shift(cf1, cf2)
      diff_cart = self.unit_cell.orthogonalize(diff_allowed)
      sum_diff_cart = python_utils.list_plus(sum_diff_cart, diff_cart)
    mean_diff_cart = [s / len(self.pairs) for s in sum_diff_cart]
    mean_diff_frac = self.unit_cell.fractionalize(mean_diff_cart)
    self.adjusted_shift = python_utils.list_plus(
      self.adjusted_shift, mean_diff_frac)

  def eliminate_weak_matches(self):
    while 1:
      max_diff = 0
      weak_pair = 0
      for pair in self.pairs:
        cf1 = self.ref_model1[pair[0]].coordinates
        cf2 = self.ref_model2[pair[1]].coordinates
        diff = self.unit_cell.Distance(cf1, cf2)
        if (diff > max_diff):
          max_diff = diff
          weak_pair = pair
      if (max_diff < self.tolerance or weak_pair == 0):
        break
      self.pairs.remove(weak_pair)
      self.refine_adjusted_shift()

  def rms(self):
    sum_diff_cart = [0,0,0]
    for pair in self.pairs:
      cf1 = self.ref_model1[pair[0]].coordinates
      cf2 = self.ref_model2[pair[1]].coordinates
      diff_frac = python_utils.list_minus(cf1, cf2)
      diff_cart = self.unit_cell.orthogonalize(diff_frac)
      sum_diff_cart = python_utils.list_plus(sum_diff_cart, diff_cart)
    mean_diff_cart = [s / len(self.pairs) for s in sum_diff_cart]
    return math.sqrt(python_utils.list_dot_product(mean_diff_cart))

def match_models(model1, model2,
                 tolerance=1., models_are_diffraction_index_equivalent=0):
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
  match_symmetry.set_sum_polar_axes()
  i_pivot1 = -1
  for pivot1 in ref_model1:
    i_pivot1 += 1
    i_pivot2 = -1
    for pivot2 in ref_model2:
      i_pivot2 += 1
      for symop in match_symmetry.rtmx:
        eq_pivot2 = symop.multiply(pivot2.coordinates)
        continuous_shift = match_symmetry.continuous_shift(
          pivot1.coordinates, eq_pivot2)
        sh_eq_pivot2 = python_utils.list_plus(eq_pivot2, continuous_shift)
        if (  unit_cell.modShortDistance(pivot1.coordinates, sh_eq_pivot2)
            < tolerance):
          matches = match_refine(tolerance,
                                 ref_model1, ref_model2,
                                 unit_cell, match_symmetry,
                                 i_pivot1, i_pivot2,
                                 symop, continuous_shift)
          print matches.symop.as_xyz(),
          print matches.adjusted_shift,
          print matches.rms()
          for pair in matches.pairs:
            print ref_model1[pair[0]].label,
            print ref_model1[pair[1]].label
          print

class test_model(model):

  def __init__(self, model_id = "SBT"):
    self.model_id = model_id
    lp = labelled_position
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

  # XXX P1 case (3 polar directions)
  # XXX 2 polar directions
  # XXX 1 poloar direction

  # XXX method to permute sites
  # XXX method to slightly move sites (but maintain site symmetry)
  # XXX add a few random sites

def run_test(argv):
  from cctbx.development import debug_utils
  Flags = debug_utils.command_line_options(argv, (
    "RandomSeed",
    "AllSpaceGroups",
    "MatchSymmetry",
    "MatchModels",
  ))
  if (not Flags.RandomSeed): debug_utils.set_random_seed(0)
  symbols_to_stdout = 0
  if (len(argv) > Flags.n):
    symbols = argv
  else:
    symbols = debug_utils.get_test_space_group_symbols(Flags.AllSpaceGroups)
    symbols_to_stdout = 1
  if (Flags.MatchModels):
    model = test_model()
    model.show("Original")
    match_models(model, model)
    return
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
    if (Flags.MatchSymmetry):
      match_symmetry = euclidean_match_symmetry(SgInfo)
      match_symmetry.show()
      assert match_symmetry.polar_axes_are_principal()

if (__name__ == "__main__"):
  import sys
  run_test(sys.argv[1:])
