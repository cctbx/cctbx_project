from __future__ import absolute_import, division, print_function

from cctbx import miller
from cctbx.development import debug_utils, random_structure
from cctbx import translation_search
from scitbx.array_family import flex
from scitbx import matrix as mat
from libtbx.test_utils import approx_equal
from libtbx.utils import progress_displayed_as_fraction
import random
import sys

from cctbx.symmetry_search import symmetrised_shifted_structure_factors
from six.moves import range


def exercise_value(f_c_in_p1, f_o, flags, n_sampled):
  """ where we check against the trusted correlation map implemented
  in cctbx.translation_search """
  crystal_gridding = f_o.crystal_gridding(
    symmetry_flags=translation_search.symmetry_flags(
      is_isotropic_search_model=False,
      have_f_part=False),
    resolution_factor=1/3)
  correlation_map = translation_search.fast_nv1995(
    gridding=crystal_gridding.n_real(),
    space_group=f_o.space_group(),
    anomalous_flag=f_o.anomalous_flag(),
    miller_indices_f_obs=f_o.indices(),
    f_obs=f_o.data(),
    f_part=flex.complex_double(), ## no sub-structure is already fixed
    miller_indices_p1_f_calc=f_c_in_p1.indices(),
    p1_f_calc=f_c_in_p1.data()).target_map()
  nx, ny, nz = n = crystal_gridding.n_real()
  def sampled_indices():
    from random import randrange
    yield (0,0,0)
    for i in range(n_sampled - 1):
      yield randrange(0, n[0]), randrange(0, n[1]), randrange(0, n[2])
  f_o_sq = f_o.as_intensity_array()
  for i,j,k in sampled_indices():
    x = (i/nx, j/ny, k/nz)
    if random.random() > 0.5: f_o_or_f_o_sq = f_o
    else: f_o_or_f_o_sq = f_o_sq
    gos = symmetrised_shifted_structure_factors(
      f_o_or_f_o_sq, f_c_in_p1, x).misfit(f_o_or_f_o_sq)
    assert approx_equal(gos.correlation, correlation_map[i,j,k]), (i,j,k)

def exercise_gradient(f_c_in_p1, f_o, flags, n_sampled):
  """ where we use finite differences to check our derivatives """
  from random import random
  from scitbx.math import approx_equal_relatively
  f_o_sq = f_o.as_intensity_array()
  for i in range(n_sampled):
    x = mat.col((random(), random(), random()))
    gos = symmetrised_shifted_structure_factors(
      f_o_sq, f_c_in_p1, x, compute_gradient=True).misfit(f_o_sq)
    for j,e in enumerate([ (1,0,0), (0,1,0), (0,0,1) ]):
      h = mat.col(e)*1e-3
      fm3, fm2, fm1, f1, f2, f3 = [
        symmetrised_shifted_structure_factors(
          f_o_sq, f_c_in_p1, y).misfit(f_o_sq).value
        for y in (x-3*h, x-2*h, x-h, x+h, x+2*h, x+3*h) ]
      finite_diff = ((f3 - fm3)/60 - 3/20*(f2 - fm2) + 3/4*(f1 - fm1))/abs(h)
      assert approx_equal_relatively(gos.gradient[j], finite_diff,
                                     relative_error=1e-2,
                                     near_zero_threshold=1e-6), \
             (j, i, tuple(x), gos.gradient[j], finite_diff)

def exercise(flags, space_group_info, n_sampled):
  symbol = space_group_info.type().hall_symbol()
  print(symbol, end=' ')
  if flags.fix_seed:
    random.seed(0)
  if not flags.include_high_symmetry:
    if space_group_info.group().order_p() > 8:
      if len(symbol) > 15: print()
      print("  [ Omitted, rerun with --include_high_symmetry to override ]")
      return
  print()
  n = int(flags.repeats)
  if n == 0: n = 1
  progress = progress_displayed_as_fraction(n)
  for i in range(n):
    xs = random_structure.xray_structure(
      space_group_info=space_group_info,
      elements=['C']*5 + ['O']*2 + ['N'],
      use_u_iso=True,
      random_u_iso=True,
      random_u_iso_scale=0.04,
      use_u_aniso=False)
    f_c = miller.build_set(
      xs, anomalous_flag=False, d_min=0.8
      ).structure_factors_from_scatterers(
        xs, algorithm='direct').f_calc()
    f_o = f_c.as_amplitude_array()
    f_c_in_p1 = f_c.expand_to_p1()
    exercise_value(f_c_in_p1, f_o, flags, n_sampled)
    exercise_gradient(f_c_in_p1, f_o, flags, n_sampled)
    progress.advance()
  progress.done()

def run():
  debug_utils.parse_options_loop_space_groups(
    sys.argv[1:],
    call_back=exercise,
    keywords=('include_high_symmetry', 'fix_seed', 'repeats'),
    symbols_to_stderr=False,
    n_sampled=50,
  )

if __name__ == '__main__':
  run()
