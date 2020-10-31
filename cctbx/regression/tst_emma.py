from __future__ import absolute_import, division, print_function
from cctbx import euclidean_model_matching as emma
from iotbx.command_line.emma import get_emma_model_from_pdb
from cctbx import xray
from cctbx import crystal
from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
import iotbx.pdb
from scitbx import matrix
from libtbx.test_utils import approx_equal, show_diff
from six.moves import cStringIO as StringIO
import random
import sys
from six.moves import range
from six.moves import zip

target_p1="""
CRYST1  144.039  144.039  178.924  90.00  90.00  90.00 P 1
HETATM 8695 ZN    ZN A   1      17.869  52.603 -22.252  1.00 71.42          ZN
HETATM 8696 ZN    ZN A   2      13.880  35.387 -29.691  1.00 52.39          ZN
HETATM 8697 ZN    ZN B   1     -18.309  55.887 -22.399  1.00 69.35          ZN
HETATM 8698 ZN    ZN B   2     -16.873  38.206 -14.862  1.00 52.32          ZN
HETATM 8699 ZN    ZN C   1      17.122  50.509  37.216  1.00 71.99          ZN
HETATM 8700 ZN    ZN C   2      13.610  33.302  29.476  1.00 64.25          ZN
HETATM 8701 ZN    ZN D   1     -18.170  53.456  36.514  1.00 69.79          ZN
HETATM 8702 ZN    ZN D   2     -17.078  36.031  44.336  1.00 54.56          ZN
HETATM 8703 ZN    ZN E   1      19.385  86.977  52.280  1.00 71.29          ZN
HETATM 8704 ZN    ZN E   2      36.676  83.705  44.841  1.00 64.49          ZN
HETATM 8705 ZN    ZN F   1      17.852  51.246  51.495  1.00 69.29          ZN
HETATM 8706 ZN    ZN F   2      35.249  52.846  59.559  1.00 57.36          ZN
"""

target_p1_inverse="""
CRYST1  144.039  144.039  178.924  90.00  90.00  90.00 P 1
ATOM      1 ZN    ZN A   1     -17.869 -52.603  22.252  1.00 71.42          ZN
ATOM      2 ZN    ZN A   2     -13.880 -35.387  29.691  1.00 52.39          ZN
ATOM      3 ZN    ZN B   1      18.309 -55.887  22.399  1.00 69.35          ZN
ATOM      4 ZN    ZN B   2      16.873 -38.206  14.862  1.00 52.32          ZN
ATOM      5 ZN    ZN C   1     -17.122 -50.509 -37.216  1.00 71.99          ZN
ATOM      6 ZN    ZN C   2     -13.610 -33.302 -29.476  1.00 64.25          ZN
ATOM      7 ZN    ZN D   1      18.170 -53.456 -36.514  1.00 69.79          ZN
ATOM      8 ZN    ZN D   2      17.078 -36.031 -44.336  1.00 54.56          ZN
ATOM      9 ZN    ZN E   1     -19.385 -86.977 -52.280  1.00 71.29          ZN
ATOM     10 ZN    ZN E   2     -36.676 -83.705 -44.841  1.00 64.49          ZN
ATOM     11 ZN    ZN F   1     -17.852 -51.246 -51.495  1.00 69.29          ZN
ATOM     12 ZN    ZN F   2     -35.249 -52.846 -59.559  1.00 57.36          ZN
"""

target_p1_partial="""
CRYST1  144.039  144.039  178.924  90.00  90.00  90.00 P 1
ATOM      1 ZN    ZN A   1      38.446  73.180   3.309  1.00 71.42          ZN
ATOM      2 ZN    ZN A   2      34.457  55.964  -4.130  1.00 52.39          ZN
ATOM      3 ZN    ZN B   1       2.268  76.464   3.162  1.00 69.35          ZN
ATOM      4 ZN    ZN B   2       3.704  58.783  10.699  1.00 52.32          ZN
ATOM      5 ZN    ZN C   1      37.699  71.086  62.777  1.00 71.99          ZN
ATOM      6 ZN    ZN C   2      34.187  53.879  55.037  1.00 64.25          ZN
ATOM      7 ZN    ZN D   1       2.407  74.033  62.075  1.00 69.79          ZN
ATOM      8 ZN    ZN D   2       3.499  56.608  69.897  1.00 54.56          ZN
ATOM      9 ZN    ZN E   1      39.962 107.554  77.841  1.00 71.29          ZN
ATOM     10 ZN    ZN E   2      57.253 104.282  70.402  1.00 64.49          ZN
ATOM     11 ZN    ZN F   1      38.429  71.823  77.056  1.00 69.29          ZN
ATOM     12 ZN    ZN F   2      55.826  73.423  85.120  1.00 57.36          ZN
"""

target_p1_inverse_partial="""
CRYST1  144.039  144.039  178.924  90.00  90.00  90.00 P 1
ATOM      1 ZN    ZN A   1       2.708 -32.026  47.813  1.00 71.42          ZN
ATOM      2 ZN    ZN A   2       6.697 -14.810  55.252  1.00 52.39          ZN
ATOM      3 ZN    ZN B   1      38.886 -35.310  47.960  1.00 69.35          ZN
ATOM      4 ZN    ZN B   2      37.450 -17.629  40.423  1.00 52.32          ZN
ATOM      5 ZN    ZN C   1       3.455 -29.932 -11.655  1.00 71.99          ZN
ATOM      6 ZN    ZN C   2       6.967 -12.725  -3.915  1.00 64.25          ZN
ATOM      7 ZN    ZN D   1      38.747 -32.879 -10.953  1.00 69.79          ZN
ATOM      8 ZN    ZN D   2      37.655 -15.454 -18.775  1.00 54.56          ZN
ATOM      9 ZN    ZN E   1       1.192 -66.400 -26.719  1.00 71.29          ZN
ATOM     10 ZN    ZN E   2     -16.099 -63.128 -19.280  1.00 64.49          ZN
ATOM     11 ZN    ZN F   1       2.725 -30.669 -25.934  1.00 69.29          ZN
ATOM     12 ZN    ZN F   2     -14.672 -32.269 -33.998  1.00 57.36          ZN
"""

target_p43212="""
CRYST1  144.039  144.039  178.924  90.00  90.00  90.00 P 43 21 2
HETATM 8695 ZN    ZN A   1      17.869  52.603 -22.252  1.00 71.42          ZN
HETATM 8696 ZN    ZN A   2      13.880  35.387 -29.691  1.00 52.39          ZN
HETATM 8697 ZN    ZN B   1     -18.309  55.887 -22.399  1.00 69.35          ZN
HETATM 8698 ZN    ZN B   2     -16.873  38.206 -14.862  1.00 52.32          ZN
HETATM 8699 ZN    ZN C   1      17.122  50.509  37.216  1.00 71.99          ZN
HETATM 8700 ZN    ZN C   2      13.610  33.302  29.476  1.00 64.25          ZN
HETATM 8701 ZN    ZN D   1     -18.170  53.456  36.514  1.00 69.79          ZN
HETATM 8702 ZN    ZN D   2     -17.078  36.031  44.336  1.00 54.56          ZN
HETATM 8703 ZN    ZN E   1      19.385  86.977  52.280  1.00 71.29          ZN
HETATM 8704 ZN    ZN E   2      36.676  83.705  44.841  1.00 64.49          ZN
HETATM 8705 ZN    ZN F   1      17.852  51.246  51.495  1.00 69.29          ZN
HETATM 8706 ZN    ZN F   2      35.249  52.846  59.559  1.00 57.36          ZN
"""

target_p43212_inverse="""
CRYST1  144.039  144.039  178.924  90.00  90.00  90.00 P 41 21 2
ATOM      1 ZN    ZN A   1     -17.869 -52.603  22.252  1.00 71.42          ZN
ATOM      2 ZN    ZN A   2     -13.880 -35.387  29.691  1.00 52.39          ZN
ATOM      3 ZN    ZN B   1      18.309 -55.887  22.399  1.00 69.35          ZN
ATOM      4 ZN    ZN B   2      16.873 -38.206  14.862  1.00 52.32          ZN
ATOM      5 ZN    ZN C   1     -17.122 -50.509 -37.216  1.00 71.99          ZN
ATOM      6 ZN    ZN C   2     -13.610 -33.302 -29.476  1.00 64.25          ZN
ATOM      7 ZN    ZN D   1      18.170 -53.456 -36.514  1.00 69.79          ZN
ATOM      8 ZN    ZN D   2      17.078 -36.031 -44.336  1.00 54.56          ZN
ATOM      9 ZN    ZN E   1     -19.385 -86.977 -52.280  1.00 71.29          ZN
ATOM     10 ZN    ZN E   2     -36.676 -83.705 -44.841  1.00 64.49          ZN
ATOM     11 ZN    ZN F   1     -17.852 -51.246 -51.495  1.00 69.29          ZN
ATOM     12 ZN    ZN F   2     -35.249 -52.846 -59.559  1.00 57.36          ZN
"""

target_p43212_half="""
CRYST1  144.039  144.039  178.924  90.00  90.00  90.00 P 43 21 2
ATOM      1 ZN    ZN A   1      89.888 124.623  67.210  1.00 71.42          ZN
ATOM      2 ZN    ZN A   2      85.899 107.406  59.771  1.00 52.39          ZN
ATOM      3 ZN    ZN B   1      53.710 127.906  67.063  1.00 69.35          ZN
ATOM      4 ZN    ZN B   2      55.146 110.225  74.600  1.00 52.32          ZN
ATOM      5 ZN    ZN C   1      89.141 122.528 126.678  1.00 71.99          ZN
ATOM      6 ZN    ZN C   2      85.629 105.321 118.938  1.00 64.25          ZN
ATOM      7 ZN    ZN D   1      53.849 125.475 125.976  1.00 69.79          ZN
ATOM      8 ZN    ZN D   2      54.941 108.050 133.798  1.00 54.56          ZN
ATOM      9 ZN    ZN E   1      91.404 158.996 141.742  1.00 71.29          ZN
ATOM     10 ZN    ZN E   2     108.695 155.724 134.303  1.00 64.49          ZN
ATOM     11 ZN    ZN F   1      89.871 123.266 140.957  1.00 69.29          ZN
ATOM     12 ZN    ZN F   2     107.268 124.865 149.021  1.00 57.36          ZN
"""

target_p43212_inverse_half="""
CRYST1  144.039  144.039  178.924  90.00  90.00  90.00 P 41 21 2
ATOM      1 ZN    ZN A   1      54.150  19.416 111.714  1.00 71.42          ZN
ATOM      2 ZN    ZN A   2      58.139  36.632 119.153  1.00 52.39          ZN
ATOM      3 ZN    ZN B   1      90.328  16.132 111.861  1.00 69.35          ZN
ATOM      4 ZN    ZN B   2      88.892  33.813 104.324  1.00 52.32          ZN
ATOM      5 ZN    ZN C   1      54.897  21.510  52.246  1.00 71.99          ZN
ATOM      6 ZN    ZN C   2      58.409  38.717  59.986  1.00 64.25          ZN
ATOM      7 ZN    ZN D   1      90.189  18.563  52.948  1.00 69.79          ZN
ATOM      8 ZN    ZN D   2      89.097  35.988  45.126  1.00 54.56          ZN
ATOM      9 ZN    ZN E   1      52.634 -14.958  37.182  1.00 71.29          ZN
ATOM     10 ZN    ZN E   2      35.343 -11.686  44.621  1.00 64.49          ZN
ATOM     11 ZN    ZN F   1      54.167  20.773  37.967  1.00 69.29          ZN
ATOM     12 ZN    ZN F   2      36.770  19.173  29.903  1.00 57.36          ZN
"""

pdb6="""
CRYST1   25.000   25.000   25.000  90.00 100.00  90.00 P 1 21 1
ATOM      1  N   GLN R  92      56.308  21.768  24.816  1.00  0.00           N
ATOM      2  CA  GLN R  92      56.603  22.162  26.183  1.00  2.50           C
ATOM      3  C   GLN R  92      56.704  23.694  26.279  1.00  0.00           C
ATOM      4  O   GLN R  92      57.301  24.360  25.424  1.00  0.00           O
ATOM      5  CB  GLN R  92      57.895  21.483  26.649  1.00  0.00           C
ATOM      6  CG  GLN R  92      58.220  21.612  28.124  1.00  0.00           C
"""
pdb5="""
CRYST1   25.000   25.000   25.000  90.00 100.00  90.00 P 1 21 1
ATOM    664  N   GLN R  92      10.579  -3.235   0.175  1.00 16.04      P9   N
ATOM    665  CA  GLN R  92      10.924  -2.975   1.574  1.00 16.32      P9   C
ATOM    666  CB  GLN R  92      12.268  -3.629   1.888  0.65 16.46      P9   C
ATOM    667  CG  GLN R  92      12.787  -3.408   3.288  0.55 16.88      P9   C
ATOM    668  CD  GLN R  92      14.105  -4.123   3.509  0.51 17.19      P9   C
"""

pdb5_out=\
"""REMARK Number of scatterers: 5
REMARK At special positions: 0
REMARK Cartesian coordinates
CRYST1   25.000   25.000   25.000  90.00 100.00  90.00 P 1 21 1
SCALE1      0.040000  0.000000  0.007053        0.00000
SCALE2      0.000000  0.040000  0.000000        0.00000
SCALE3      0.000000  0.000000  0.040617        0.00000
ATOM      1  N    N      1      10.579  21.833   0.175  1.00 16.04           N
ATOM      2  CA   CA     2      10.924  22.093   1.574  1.00 16.32           C
ATOM      3  CB   CB     3      12.268  21.439   1.888  0.65 16.46           C
ATOM      4  CG   CG     4      12.787  21.660   3.288  0.55 16.88           C
ATOM      5  CD   CD     5      14.105  20.945   3.509  0.51 17.19           C
END"""


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
    print("total matches:", len(refined_matches))
    print("solutions:", solution_counter)
    print()
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
        raise RuntimeError("Unknown model_id: " + model_id)
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
    shift = [0.5 - random.random() for i in range(3)]
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
      for i_site1 in range(len(model1.positions())):
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

def tst_pdb_output():
  print("Testing pdb-output option")
  xray_scatterer = xray.scatterer( scattering_type = 'SE')
  for sg,target_list in zip(
     ['p1','p43212'],
     [
        [target_p1,target_p1_inverse,target_p1_partial,
           target_p1_inverse_partial],
        [target_p43212,target_p43212_inverse,target_p43212_half,
           target_p43212_inverse_half]
     ]
     ):
    print("Testing group of targets in %s" %(sg))
    for t1 in target_list:
      e1=get_emma_model_from_pdb(pdb_records=t1)
      for t2 in target_list:
        e2=get_emma_model_from_pdb(pdb_records=t2)
        match_list=e1.best_superpositions_on_other(
          e2)
        match=match_list[0]
        assert match
        offset_e2=match.get_transformed_model2()

        # make sure that offset_i2 is pretty much the same as e1 now.
        new_match_list=e1.best_superpositions_on_other(
          offset_e2)
        new_match=new_match_list[0]
        assert new_match
        assert approx_equal(new_match.rms,0.,eps=0.01)
        assert len(new_match.pairs)==12
        assert approx_equal(
           new_match.rt.r,matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1)))
        assert approx_equal(new_match.rt.t.transpose(),matrix.col((0, 0, 0)))

  print("Testing pdb-output option with different-sized entries")
  e1=get_emma_model_from_pdb(pdb_records=pdb6)
  e2=get_emma_model_from_pdb(pdb_records=pdb5)
  pdb_inp_e2=iotbx.pdb.input(source_info=None, lines=pdb5)

  match_list=e1.best_superpositions_on_other(e2)
  match=match_list[0]
  assert match
  output_pdb="output.pdb"
  offset_e2=match.get_transformed_model2(output_pdb=output_pdb,
    template_pdb_inp=pdb_inp_e2)
  with open(output_pdb) as f:
    offset_e2_text_lines=f.readlines()
  for o,e in zip(offset_e2_text_lines, pdb5_out.splitlines()):
     o=o.strip()
     e=e.strip()
     if o != e:
       print(o)
       print(e)
       assert o==e

if (__name__ == "__main__"):
  run()
  tst_pdb_output()
