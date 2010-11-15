from cctbx import euclidean_model_matching as emma
from cctbx import crystal
from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from scitbx import matrix
from libtbx.test_utils import approx_equal, show_diff
from cStringIO import StringIO
import random
import sys

def verify_match(model1, model2, tolerance, match_rt, pairs):
  adj_tolerance = tolerance * (1 + 1.e-6)
  for pair in pairs:
    c1 = model1[pair[0]].site
    c2 = match_rt * model2[pair[1]].site
    equiv_c2 = sgtbx.sym_equiv_sites(model1.site_symmetry(c2.elems))
    dist_info = sgtbx.min_sym_equiv_distance_info(equiv_c2, c1)
    assert dist_info.dist() < adj_tolerance, str(model1.space_group_info())

def analyze_singles(model, singles):
  for i in singles:
    if (model[i].label.startswith("S")): return False
  return True

def analyze_refined_matches(model1, model2, refined_matches, verbose):
  solution_counter = 0
  for match in refined_matches:
    if (0 or verbose):
      match.show()
    verify_match(match.ref_model1, match.ref_model2, match.tolerance,
                 match.ref_eucl_rt, match.pairs)
    verify_match(model1, model2, match.tolerance,
                 match.rt, match.pairs)
    if (    analyze_singles(model1, match.singles1)
        and analyze_singles(model2, match.singles2)):
      solution_counter += 1
  if (0 or verbose):
    print "total matches:", len(refined_matches)
    print "solutions:", solution_counter
    print
  assert solution_counter != 0

class test_model(emma.model):

  def __init__(self, model_id="SBT", n_elements=4):
    if (model_id is None): return
    self.model_id = model_id
    pos = emma.position
    if (type(model_id) == type("")):
      if (model_id == "SBT"):
        emma.model.__init__(self,
          crystal.special_position_settings(crystal.symmetry(
            (16.8986, 16.8986, 16.8986, 61.1483, 61.1483, 61.1483),
            "R -3 m :R")),
          (pos("SI1", (-0.3584, 0.2844, 0.4622)),
           pos("SI2", (-0.2133, 0.9659, -0.6653)),
           pos("SI3", (-0.8358, 0.7, 0.3431)),
           pos("SI4", (0.4799, 1.836, 0.6598))))
      else:
        raise RuntimeError, "Unknown model_id: " + model_id
    else:
      structure = random_structure.xray_structure(
        model_id,
        elements=["S"]*n_elements,
        volume_per_atom=50.,
        min_distance=2.0)
      positions = []
      for scatterer in structure.scatterers():
        positions.append(emma.position(scatterer.label, scatterer.site))
      emma.model.__init__(self, structure, positions)

  def create_new_test_model(self, new_positions):
    new_test_model = test_model(None)
    emma.model.__init__(new_test_model, self, new_positions)
    return new_test_model

  def shuffle_positions(self):
    shuffled_positions = list(self.positions())
    random.shuffle(shuffled_positions)
    return self.create_new_test_model(shuffled_positions)

  def random_symmetry_mates(self):
    new_positions = []
    for pos in self.positions():
      equiv = sgtbx.sym_equiv_sites(self.site_symmetry(pos.site))
      i = random.randrange(equiv.coordinates().size())
      new_positions.append(emma.position(pos.label, equiv.coordinates()[i]))
    return self.create_new_test_model(new_positions)

  def apply_random_eucl_op(self, models_are_diffraction_index_equivalent=0):
    match_symmetry = emma.euclidean_match_symmetry(
      self.space_group_info(),
      use_k2l=True, use_l2n=(not models_are_diffraction_index_equivalent))
    i = random.randrange(match_symmetry.rt_mx.order_z())
    eucl_symop = match_symmetry.rt_mx(i)
    shift = [0.5 - random.random() for i in xrange(3)]
    allowed_shift = matrix.col(match_symmetry.filter_shift(shift, selector=1))
    new_positions = []
    for pos in self.positions():
      new_site = matrix.col(eucl_symop * pos.site) + allowed_shift
      new_positions.append(emma.position(pos.label, new_site.elems))
    return self.create_new_test_model(new_positions)

  def add_random_positions(self, number_of_new_positions=3, label="R",
                           min_distance=1.0):
    existing_sites = []
    for pos in self.positions():
      existing_sites.append(pos.site)
    new_sites = random_structure.random_sites(
      special_position_settings=self,
      existing_sites=existing_sites,
      n_new=number_of_new_positions,
      min_hetero_distance=min_distance,
      general_positions_only=False)
    new_positions = []
    i = 0
    for site in new_sites:
      i += 1
      new_positions.append(emma.position("%s%d" % (label, i), site))
    return self.create_new_test_model(self.positions() + new_positions)

  def shake_positions(self, gauss_sigma=0.2, min_distance=1.0):
    new_positions = []
    for pos in self.positions():
      new_coor = random_structure.random_modify_site(
        special_position_settings=self,
        site=pos.site,
        gauss_sigma=gauss_sigma,
        max_distance=min_distance*0.99)
      new_positions.append(emma.position(pos.label, new_coor))
    return self.create_new_test_model(new_positions)

  def random_hand(self):
    if (random.random() < 0.5):
      cb_op = sgtbx.change_of_basis_op()
    else:
      cb_op = self.space_group_info().type().change_of_hand_op()
    return self.change_basis(cb_op).reset_cb_op()

def run_call_back(flags, space_group_info):
  verbose = flags.Verbose
  if (flags.StaticModels):
    model1 = (test_model()
      .add_random_positions(2, "A")
      .shuffle_positions()
      .random_symmetry_mates()
      .apply_random_eucl_op()
      )
    model2 = (test_model()
      .add_random_positions(3, "B")
      .shuffle_positions()
      .random_symmetry_mates()
      .apply_random_eucl_op()
      .shake_positions()
      .random_hand()
      )
    for i in (0,1):
      m1 = model1
      if (i): m1 = model1.transform_to_reference_setting().reset_cb_op()
      for j in (0,1):
        m2 = model2
        if (j): m2 = model2.transform_to_reference_setting().reset_cb_op()
        if (0 or verbose):
          m1.show("Model1(%d)" % (i,))
          m2.show("Model2(%d)" % (j,))
        model_matches = emma.model_matches(m1, m2, rms_penalty_per_site=0)
        analyze_refined_matches(m1, m2, model_matches.refined_matches, verbose)
    return False
  model_core = test_model(space_group_info)
  model1 = (model_core
    .add_random_positions(2, "A")
    )
  model2 = (model_core
    .add_random_positions(3, "B")
    .shuffle_positions()
    .random_symmetry_mates()
    .apply_random_eucl_op()
    .shake_positions()
    .random_hand()
    )
  if (0 or verbose):
    model_core.show("Core")
    model1.show("Model1")
    model2.show("Model2")
  model_matches = emma.model_matches(model1, model2, rms_penalty_per_site=0)
  analyze_refined_matches(
    model1, model2, model_matches.refined_matches, verbose)
  assert model_matches.consensus_model().size() >= model1.size()-2
  assert model_matches.consensus_model(i_model=2).size() >= model2.size()-3
  model1.expand_to_p1()
  model2.as_xray_structure()
  for i1,i2,m1,m2 in [(1,2,model1,model2),(2,1,model2,model1)]:
    m2_t = model_matches.transform_model(i_model=i2)
    assert m1.unit_cell().is_similar_to(m2_t.unit_cell())
    assert m1.space_group() == m2_t.space_group()
    for pair in model_matches.refined_matches[0].pairs:
      site1 = m1.positions()[pair[i1-1]].site
      site2 = m2_t.positions()[pair[i2-1]].site
      equiv_sites1 = sgtbx.sym_equiv_sites(m1.site_symmetry(site1))
      dist_info = sgtbx.min_sym_equiv_distance_info(equiv_sites1, site2)
      assert dist_info.dist() < model_matches.tolerance + 1.e-6
      site2_closest = dist_info.sym_op() * site2
      assert approx_equal(
        m1.unit_cell().distance(site1, site2_closest),
        dist_info.dist())
    if (i1 == 1):
      singles = model_matches.refined_matches[0].singles2
    else:
      singles = model_matches.refined_matches[0].singles1
    for i_site2 in singles:
      site2 = m2_t.positions()[i_site2].site
      for i_site1 in xrange(len(model1.positions())):
        site1 = m1.positions()[i_site1].site
        equiv_sites1 = sgtbx.sym_equiv_sites(m1.site_symmetry(site1))
        dist_info = sgtbx.min_sym_equiv_distance_info(equiv_sites1, site2)
        if (dist_info.dist() < model_matches.tolerance - 1.e-6):
          ok = False
          for pair in model_matches.refined_matches[0].pairs:
            if (pair[i1-1] == i_site1):
              ok = True
          assert ok

def run():
  match_symmetry = emma.euclidean_match_symmetry(
    space_group_info=sgtbx.space_group_info(symbol="P1"),
    use_k2l=True,
    use_l2n=False)
  out = StringIO()
  match_symmetry.show(title="test", f=out)
  assert not show_diff(out.getvalue(), """\
euclidean_match_symmetry: test
P -1
((1, 0, 0), (0, 1, 0), (0, 0, 1))
""")
  #
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back, (
    "StaticModels",))

if (__name__ == "__main__"):
  run()
