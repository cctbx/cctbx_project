from __future__ import absolute_import, division, print_function

from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
import boost_adaptbx.boost.python as bp
from six.moves import range
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")
from cctbx_asymmetric_map_ext import *
from cctbx.array_family import flex
from cctbx import maptbx
import mmtbx.masks
import random

if (1):
  random.seed(0)
  flex.set_random_seed(0)

def run_group(symbol, preprocess_against_shallow):
  """
  Regions of a symmetric whole-cell mask, of the same mask with symmetry
  related regions merged, and of the mask restored from its asymmetric unit
  must be the same. After merging, the map of region numbers must obey the
  space group symmetry. The space group is given to connectivity so that removing
  shallow voxels (preprocess_against_shallow) keeps the mask symmetric: the
  three grid axes alone are not a symmetric set of directions in trigonal and
  hexagonal space groups.
  """
  group = space_group_info(symbol)
  print("\n== space group %d"%symbol)
  xrs = random_structure.xray_structure(
    space_group_info       = group,
    volume_per_atom        = 15.,
    general_positions_only = False,
    elements               = ('C', 'N', 'O', 'H')*31,
    min_distance           = 1.0)
  sgt = xrs.space_group().type()
  #
  cg = maptbx.crystal_gridding(
    unit_cell        = xrs.unit_cell(),
    space_group_info = xrs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = 0.4)
  n_real = cg.n_real()
  mask_p1 = mmtbx.masks.mask_from_xray_structure(
    xray_structure        = xrs,
    p1                    = True,
    for_structure_factors = True,
    n_real                = n_real,
    in_asu                = False).mask_data
  maptbx.unpad_in_place(map=mask_p1)
  assert flex.min(mask_p1)==0
  assert flex.max(mask_p1)==1
  #
  co = maptbx.connectivity(
    map_data                   = mask_p1,
    threshold                  = 0.01,
    preprocess_against_shallow = preprocess_against_shallow,
    wrapping                   = True,
    space_group                = xrs.space_group())
  #
  print("Regions in P1")
  regions_p1 = list(co.regions())
  s1 = flex.sum(flex.int(regions_p1))
  print(regions_p1, s1)
  conn_map_p1 = co.result().as_double()
  print(flex.min(conn_map_p1), flex.max(conn_map_p1))
  #
  print("Merge symmetry related")
  co.merge_symmetry_related_regions(space_group = xrs.space_group())
  conn_map_p1_merged = co.result().as_double()
  regions_p1_merged = list(co.regions())
  s2 = flex.sum(flex.int(regions_p1_merged))
  print(list(regions_p1_merged), s2)
  amap = asu_map_ext.asymmetric_map(sgt, conn_map_p1_merged)
  conn_map_asu = amap.data()
  conn_map_p1_restored = amap.symmetry_expanded_map()
  print(flex.min(conn_map_asu), flex.max(conn_map_asu))
  # All symmetry copies of a region have the same number after merging: the
  # map of region numbers obeys the space group symmetry.
  restored_unpadded = conn_map_p1_restored.deep_copy()
  maptbx.unpad_in_place(map=restored_unpadded)
  assert restored_unpadded.all_eq(conn_map_p1_merged)
  # Merging changes numbers only: the solvent points are the same, and the
  # number of regions cannot increase.
  assert ((conn_map_p1 > 0) == (conn_map_p1_merged > 0)).all_eq(True)
  assert len(regions_p1_merged) <= len(regions_p1)
  assert regions_p1_merged[0] == regions_p1[0]
  #
  mask_p1_1 = conn_map_p1_restored.set_selected(conn_map_p1_restored>0.01, 1)
  maptbx.unpad_in_place(map=mask_p1_1)
  co = maptbx.connectivity(
    map_data                   = mask_p1_1,
    threshold                  = 0.01,
    preprocess_against_shallow = preprocess_against_shallow,
    wrapping                   = True,
    space_group                = xrs.space_group())
  print("Restored")
  regions_p1_restored = list(co.regions())
  s3 = flex.sum(flex.int(regions_p1_restored))
  print(regions_p1_restored, s3)
  conn_map_p1_restored = co.result().as_double()
  print(flex.min(conn_map_p1_restored), flex.max(conn_map_p1_restored))
  assert regions_p1 == regions_p1_restored
  #
  assert s1 == s2
  assert s2 == s3

def run(preprocess_against_shallow):
  for i in range(1,231):
    run_group(i,preprocess_against_shallow);

if (__name__ == "__main__"):
  for ppas in [True, False]:
    run(preprocess_against_shallow = ppas)
