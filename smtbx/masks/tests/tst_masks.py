from __future__ import division

from cctbx import miller
from cctbx import euclidean_model_matching as emma
from cctbx.array_family import flex
from cctbx.masks import flood_fill
from cctbx.xray import structure
from libtbx.test_utils import approx_equal, show_diff
import libtbx.utils
from smtbx import masks
from smtbx.refinement import constraints
import smtbx.utils

import cStringIO

def exercise_masks():
  mt = flex.mersenne_twister(seed=0)
  xs_ref = structure.from_shelx(
    file=cStringIO.StringIO(YAKRUY_ins))
  mi = xs_ref.crystal_symmetry().build_miller_set(
    d_min=0.5, anomalous_flag=False)
  fo = mi.structure_factors_from_scatterers(
    xs_ref, algorithm="direct").f_calc().as_amplitude_array()
  k = 0.05 + 10 * mt.random_double()
  fo = fo.customized_copy(data=fo.data()*k)
  fo2 = fo.f_as_f_sq()
  acetonitrile_sel = xs_ref.label_selection(
    'N4', 'C20', 'C21', 'H211', 'H212', 'H213')
  xs_no_sol = xs_ref.deep_copy_scatterers().select(
    acetonitrile_sel, negate=True)
  # check what happens when no voids are found
  mask = masks.mask(xs_ref, fo2)
  mask.compute(solvent_radius=1.2,
               shrink_truncation_radius=1.2,
               resolution_factor=1/2,
               atom_radii_table={'C':1.70, 'B':1.63, 'N':1.55, 'O':1.52})
  assert mask.structure_factors() is None
  assert mask.n_voids() == 0
  assert mask.n_solvent_grid_points() == 0
  assert mask.f_mask() is None
  assert mask.f_model() is None
  assert mask.modified_intensities() is None
  assert mask.f_000 is None
  s = cStringIO.StringIO()
  mask.show_summary(log=s)
  assert not show_diff(s.getvalue(), """\
solvent_radius: 1.20
shrink_truncation_radius: 1.20
van der Waals radii:
    B     C     H     N     O
 1.63  1.70  1.20  1.55  1.52

Total solvent accessible volume / cell = 0.0 Ang^3 [0.0%]

gridding: (30,45,54)
""")
  # and now with some voids
  fo2_complete = fo2.sort()
  fo2_missing_1 = fo2.select_indices(flex.miller_index([(0,0,1),
                                                        ]), negate=True)
  mt = flex.mersenne_twister(seed=0)
  fo2_incomplete = fo2.select(mt.random_bool(fo2.size(), 0.95))

  for fo2, use_space_group_symmetry in zip(
    (fo2_complete, fo2_complete, fo2_missing_1, fo2_incomplete),
    (True, False, True, True)):
    if fo2 is fo2_complete: use_complete_set=False
    else: use_complete_set=True
    mask = masks.mask(xs_no_sol, fo2, use_complete_set=use_complete_set)
    mask.compute(solvent_radius=1.2,
                 shrink_truncation_radius=1.2,
                 resolution_factor=1/3,
                 #atom_radii_table={'C':1.70, 'B':1.63, 'N':1.55, 'O':1.52},
                 use_space_group_symmetry=use_space_group_symmetry)
    n_voids = flex.max(mask.mask.data) - 1
    f_mask = mask.structure_factors()
    f_model = mask.f_model()
    modified_fo = mask.modified_intensities().as_amplitude_array()
    f_obs_minus_f_model = fo.common_set(f_model).f_obs_minus_f_calc(f_obs_factor=1/k, f_calc=f_model)
    diff_map = miller.fft_map(mask.crystal_gridding, f_obs_minus_f_model)
    diff_map.apply_volume_scaling()
    stats = diff_map.statistics()
    assert n_voids == 2
    assert approx_equal(n_voids, mask.n_voids())
    assert mask.n_solvent_grid_points() == 42148
    if fo2 is fo2_complete:
      # check the difference map has no large peaks/holes
      assert max(stats.max(), abs(stats.min())) < 0.11
    # expected electron count: 44
    assert approx_equal(mask.f_000_s, 44, eps=1)
    assert modified_fo.r1_factor(mask.f_calc.common_set(modified_fo), k) < 0.006
    assert fo.common_set(fo2).r1_factor(f_model.common_set(fo2), k) < 0.006

  s = cStringIO.StringIO()
  mask.show_summary(log=s)
  assert not show_diff(s.getvalue(), """\
solvent_radius: 1.20
shrink_truncation_radius: 1.20
van der Waals radii:
    C     H     N     O
 1.77  1.20  1.50  1.45

Total solvent accessible volume / cell = 146.5 Ang^3 [16.3%]
Total electron count / cell = 43.0

gridding: (45,72,80)
Void #Grid points Vol/A^3 Vol/%  Centre of mass (frac)   Eigenvectors (frac)
   1        21074    73.3   8.1  ( 0.267, 0.461, 0.672)  1  ( 0.982, 0.126, 0.142)
                                                         2  (-0.166, 0.206, 0.964)
                                                         3  (-0.092, 0.970,-0.223)
   2        21074    73.3   8.1  (-0.267, 0.539, 0.328)  1  ( 0.982, 0.126, 0.142)
                                                         2  (-0.166, 0.206, 0.964)
                                                         3  (-0.092, 0.970,-0.223)

Void  Vol/Ang^3  #Electrons
   1       73.3         21.5
   2       73.3         21.5
""")
  cif_block = mask.as_cif_block()

  fo2 = fo.f_as_f_sq()
  # this bit is necessary until we have constraints, as
  # otherwise the hydrogens just disappear into the ether.
  xs = xs_no_sol.deep_copy_scatterers()
  h_selection = xs.element_selection('H')
  orig_flags = xs.scatterer_flags()
  flags = orig_flags.deep_copy()
  for flag, is_h in zip(flags, h_selection):
    if is_h:
      flag.set_grads(False)
  xs.set_scatterer_flags(flags)

  # first refine with no mask
  xs = exercise_least_squares(xs, fo2, mask=None)
  xs.set_scatterer_flags(orig_flags)
  for i in range(1):
    # compute improved mask/f_mask
    mask = masks.mask(xs, fo2)
    mask.compute(solvent_radius=1.2,
                 shrink_truncation_radius=1.2,
                 atom_radii_table={'C':1.70, 'B':1.63, 'N':1.55, 'O':1.52},
                 resolution_factor=1/3)
    mask.structure_factors()
    xs = exercise_least_squares(xs, fo2, mask)
  # again exclude hydrogens from tests because of lack of constraints
  emma_ref = xs_no_sol.select(h_selection, negate=True).as_emma_model()
  match = emma.model_matches(emma_ref, xs.select(
    h_selection, negate=True).as_emma_model()).refined_matches[0]
  assert approx_equal(match.rms, 0, eps=1e-3)

def exercise_least_squares(xray_structure, fo_sq, mask=None):
  from smtbx.refinement import least_squares
  fo_sq = fo_sq.customized_copy(sigmas=flex.double(fo_sq.size(),1.))
  xs = xray_structure.deep_copy_scatterers()
  if mask is not None:
    f_mask = mask.f_mask()
  else:
    f_mask = None
  connectivity_table = smtbx.utils.connectivity_table(xs)
  reparametrisation = constraints.reparametrisation(
    structure=xs,
    geometrical_constraints=[],
    connectivity_table=connectivity_table)
  normal_eqns = least_squares.normal_equations(
    fo_sq, reparametrisation,
    f_mask=f_mask,
    weighting_scheme="default")
  objectives = []
  scales = []
  fo_sq_max = flex.max(fo_sq.data())
  for i in xrange(3):
    normal_eqns.build_up()
    objectives.append(normal_eqns.objective)
    scales.append(normal_eqns.scale_factor)
    gradient_relative_norm = normal_eqns.gradient.norm()/fo_sq_max
    normal_eqns.solve_and_apply_shifts()
    shifts = normal_eqns.shifts
  return xs

YAKRUY_ins = """
CELL 0.71073   7.086  10.791  12.850 104.16 105.87  95.86
ZERR    2.00   0.001   0.003   0.002   0.02   0.02   0.01
LATT  1
SFAC C H N O
UNIT 42 40 8 4
O1    4    0.263832    0.309851    0.916785    11.00000    0.03486    0.02047 =
         0.02294    0.00910    0.01221    0.00848
O2    4    0.259857    0.121826    0.750368    11.00000    0.03493    0.02585 =
         0.01911    0.00907    0.01212    0.00849
N1    3    0.244867   -0.219643    1.065706    11.00000    0.02111    0.01983 =
         0.02052    0.00529    0.00767    0.00528
N2    3    0.238512   -0.124207    1.256893    11.00000    0.02403    0.01605 =
         0.01797    0.00354    0.00780    0.00428
N3    3    0.236863   -0.003326    1.323510    11.00000    0.02340    0.01743 =
         0.02104    0.00455    0.00803    0.00515
N4    3    0.313652    0.621206    0.720884    11.00000    0.05212    0.03488 =
         0.03702    0.00851    0.01917    0.00007
C1    1    0.258790    0.192265    0.938230    11.00000    0.01801    0.02150 =
         0.02273    0.00950    0.00637    0.00518
C2    1    0.258684    0.087261    0.845526    11.00000    0.01829    0.02676 =
         0.01813    0.00919    0.00712    0.00512
C3    1    0.258107   -0.035539    0.857494    11.00000    0.01763    0.02320 =
         0.01795    0.00465    0.00660    0.00526
H3    2    0.265419   -0.106192    0.797272    11.00000    0.02235
C4    1    0.252521   -0.058937    0.960660    11.00000    0.01442    0.02188 =
         0.01933    0.00569    0.00576    0.00362
C5    1    0.249654    0.044792    1.051770    11.00000    0.01428    0.02074 =
         0.01927    0.00546    0.00489    0.00333
C6    1    0.254907    0.171444    1.038876    11.00000    0.01823    0.01960 =
         0.01938    0.00441    0.00607    0.00410
H6    2    0.260661    0.244361    1.101313    11.00000    0.02028
C7    1    0.244340    0.011338    1.151454    11.00000    0.01528    0.01760 =
         0.02039    0.00560    0.00559    0.00358
C8    1    0.242078   -0.118214    1.151436    11.00000    0.01632    0.02066 =
         0.01923    0.00609    0.00594    0.00402
C9    1    0.250274   -0.186707    0.973843    11.00000    0.01810    0.02164 =
         0.01942    0.00361    0.00684    0.00426
H9    2    0.252942   -0.255996    0.908935    11.00000    0.01373
C10   1    0.238915    0.077502    1.260637    11.00000    0.01646    0.02046 =
         0.02005    0.00660    0.00460    0.00293
C11   1    0.241179   -0.230214    1.303905    11.00000    0.01923    0.01988 =
         0.02379    0.01016    0.00832    0.00634
C12   1    0.301101   -0.204248    1.421081    11.00000    0.02960    0.02212 =
         0.02187    0.00694    0.00858    0.00529
H12   2    0.347647   -0.113158    1.469110    11.00000    0.03067
C13   1    0.297066   -0.306624    1.468541    11.00000    0.03346    0.03026 =
         0.02258    0.01182    0.00920    0.00828
H13   2    0.338439   -0.286559    1.549032    11.00000    0.03493
C14   1    0.238170   -0.433815    1.401003    11.00000    0.03483    0.02587 =
         0.03168    0.01598    0.01317    0.01089
H14   2    0.239203   -0.504609    1.435015    11.00000    0.03255
C15   1    0.182288   -0.458453    1.284537    11.00000    0.03573    0.01878 =
         0.02946    0.00736    0.01196    0.00646
H15   2    0.141861   -0.547777    1.233813    11.00000    0.03479
C16   1    0.182000   -0.357827    1.235343    11.00000    0.02803    0.02211 =
         0.02063    0.00598    0.00799    0.00558
H16   2    0.140443   -0.375280    1.154110    11.00000    0.02660
C17   1    0.228172    0.216669    1.308393    11.00000    0.02843    0.01824 =
         0.02238    0.00458    0.00869    0.00462
H171  2    0.343696    0.277860    1.305773    11.00000    0.03122
H172  2    0.233564    0.235853    1.390463    11.00000    0.03214
H173  2    0.103775    0.237732    1.266193    11.00000    0.03633
C18   1    0.262699    0.418567    1.006264    11.00000    0.03735    0.01968 =
         0.02610    0.00669    0.01069    0.00556
H181  2    0.384186    0.436044    1.074408    11.00000    0.02881
H182  2    0.136744    0.404648    1.027622    11.00000    0.02787
H183  2    0.269188    0.493767    0.976527    11.00000    0.02735
C19   1    0.246327    0.019380    0.652743    11.00000    0.03131    0.03014 =
         0.01841    0.00608    0.01024    0.00667
H191  2    0.369862   -0.025013    0.666318    11.00000    0.03203
H192  2    0.122807   -0.042860    0.628185    11.00000    0.03308
H193  2    0.249416    0.060472    0.594291    11.00000    0.03051
C20   1    0.295027    0.510549    0.696308    11.00000    0.03015    0.03545 =
         0.02324    0.01128    0.01118    0.00361
C21   1    0.271551    0.369674    0.665707    11.00000    0.03896    0.03186 =
         0.04785    0.01623    0.01900    0.00949
H211  2    0.372427    0.342541    0.630638    11.00000    0.06494
H212  2    0.140789    0.331654    0.610485    11.00000    0.06386
H213  2    0.292381    0.339113    0.730897    11.00000    0.08983
HKLF 4
END
"""

def run():
  libtbx.utils.show_times_at_exit()
  exercise_masks()

if __name__ == '__main__':
  run()
