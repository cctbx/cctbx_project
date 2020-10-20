"""\
Essence of bulk-solvent mask calculation, testing interaction of
symmetry summation in maptbx.structure_factors.from_map with handling
of grid points. This script demonstrates that grid points redundant
under symmetry must be kept at zero and that all grid points inside
the asu have to be weighted by site-multiplicity/order_z.
As an aside, this script also exercises space_group.multiplicity(site).
"""

from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from cctbx import miller
from cctbx import sgtbx
from cctbx.development import debug_utils
from scitbx import fftpack
from scitbx.array_family import flex
from libtbx.math_utils import iround
from libtbx.test_utils import approx_equal
from libtbx.utils import n_dim_index_from_one_dim
import boost_adaptbx.boost.rational
import random
from six.moves import zip

def exercise(space_group_info, redundancy_counter=0):
  n_real = (12,12,12)
  miller_max = (2,2,2)
  gt = maptbx.grid_tags(n_real)
  uc = space_group_info.any_compatible_unit_cell(volume=1000)
  fl = sgtbx.search_symmetry_flags(use_space_group_symmetry=True)
  gt.build(space_group_info.type(), fl)
  fft = fftpack.real_to_complex_3d(n_real)
  map0 = flex.double(flex.grid(fft.m_real()).set_focus(fft.n_real()), 0)
  weight_map = map0.deep_copy()
  map = map0.deep_copy()
  ta = gt.tag_array()
  order_z = space_group_info.group().order_z()
  problems_expected = (redundancy_counter != 0)
  for ijk in flex.nested_loop(n_real):
    t = ta[ijk]
    if (t < 0):
      xyz = [i/n for i,n in zip(ijk, n_real)]
      ss = sgtbx.site_symmetry(
        unit_cell=uc,
        space_group=space_group_info.group(),
        original_site=xyz,
        min_distance_sym_equiv=1e-5)
      m = space_group_info.group().multiplicity(
        site=boost_adaptbx.boost.rational.vector(ijk, n_real))
      assert m == ss.multiplicity()
      w = m / order_z
      weight_map[ijk] = w
      map[ijk] = w
    elif (redundancy_counter != 0):
      redundancy_counter -= 1
      ijk_asu = n_dim_index_from_one_dim(i1d=t, sizes=n_real)
      assert ta.accessor()(ijk_asu) == t
      map[ijk] = map[ijk_asu]
  sf_map = fft.forward(map)
  del map
  mi = miller.index_generator(
    space_group_info.type(), False, miller_max).to_array()
  assert mi.size() != 0
  from_map = maptbx.structure_factors.from_map(
    space_group=space_group_info.group(),
    anomalous_flag=False,
    miller_indices=mi,
    complex_map=sf_map,
    conjugate_flag=True)
  sf = [iround(abs(f)) for f in from_map.data()]
  if (sf != [0]*len(sf)):
    assert problems_expected
    return
  else:
    not problems_expected
  #
  map_p1 = map0.deep_copy()
  map_sw = map0.deep_copy()
  for ijk in flex.nested_loop(n_real):
    t = ta[ijk]
    if (t < 0):
      v = random.random()*2-1
      map_p1[ijk] = v
      map_sw[ijk] = v * weight_map[ijk]
    else:
      ijk_asu = n_dim_index_from_one_dim(i1d=t, sizes=n_real)
      assert ta.accessor()(ijk_asu) == t
      assert map_p1[ijk_asu] != 0
      map_p1[ijk] = map_p1[ijk_asu]
  #
  # fft followed by symmetry summation in reciprocal space
  sf_map_sw = fft.forward(map_sw)
  del map_sw
  sf_sw = maptbx.structure_factors.from_map(
    space_group=space_group_info.group(),
    anomalous_flag=False,
    miller_indices=mi,
    complex_map=sf_map_sw,
    conjugate_flag=True).data()
  del sf_map_sw
  #
  # symmetry expansion in real space (done above already) followed fft
  sf_map_p1 = fft.forward(map_p1)
  del map_p1
  sf_p1 = maptbx.structure_factors.from_map(
    space_group=sgtbx.space_group(),
    anomalous_flag=False,
    miller_indices=mi,
    complex_map=sf_map_p1,
    conjugate_flag=True).data()
  del sf_map_p1
  #
  corr = flex.linear_correlation(x=flex.abs(sf_sw), y=flex.abs(sf_p1))
  assert corr.is_well_defined
  assert approx_equal(corr.coefficient(), 1)

def run_call_back(flags, space_group_info):
  exercise(space_group_info=space_group_info)

def run(args):
  exercise(
    space_group_info=sgtbx.space_group_info(symbol="P222"),
    redundancy_counter=1)
  debug_utils.parse_options_loop_space_groups(args, run_call_back)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
