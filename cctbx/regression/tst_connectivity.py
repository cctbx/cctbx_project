from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import time
import iotbx.pdb
from cctbx import maptbx
from cctbx import miller
from six.moves import range
from six.moves import zip
from libtbx.test_utils import approx_equal
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_asymmetric_map_ext")
from cctbx_asymmetric_map_ext import *


def getvs(cmap, threshold, wrap=True):
  co = maptbx.connectivity(map_data=cmap, threshold=threshold, wrapping=wrap)
  map_result = co.result()
  regs = co.regions()
  coors = co.maximum_coors()
  vals = co.maximum_values()
  assert len(list(regs)) == len(list(coors)) == len(list(vals))
  # check dimensions
  assert cmap.all() == map_result.all()
  v=[0,0,0]
  for i in range(3):
    v[i] = (map_result==i).count(True)
  return v, list(co.regions())

def exercise1():
  pdb_str="""
CRYST1   10.000  10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       2.000   2.000   2.000  1.00 20.00           C
END
"""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  cg = maptbx.crystal_gridding(unit_cell=xrs.unit_cell(),
      pre_determined_n_real=(100,100,100),
      space_group_info=xrs.space_group_info())
  fc = xrs.structure_factors(d_min = 1., algorithm = "direct").f_calc()
  fft_map = miller.fft_map(crystal_gridding=cg, fourier_coefficients=fc)
  map_data = fft_map.real_map_unpadded()

  # pass map and threshold value
  co = maptbx.connectivity(map_data=map_data, threshold=100.)
  # get 'map' of the same size with integers: 0 where below threshold,
  # 1,2,3... - for connected regions
  map_result = co.result()
  # to find out the number of connected region for particular point:
  assert map_result[0,0,0] == 0    # means under threshold
  assert map_result[20,20,20] == 1 # blob 1

  # get 1d array of integer volumes and transform it to list.
  volumes = list(co.regions())
  # find max volume (except volume of 0-region which will be probably max)
  max_volume = max(volumes[1:])
  # find number of the region with max volume
  max_index = volumes.index(max_volume)
  v=[0,0,0]
  for i in range(3):
    # !!! Do not do this because it's extremely slow! Used for test purposes.
    v[i] = (map_result==i).count(True)

  assert v[2] == 0
  assert v[1] < 15000
  assert v[0]+v[1]+v[2] == 1000000
  assert volumes == v[:2]

def exercise3():
  pdb_str="""
CRYST1   10.000  10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       2.000   2.000   2.000  1.00  2.00           C
HETATM    1  C    C      1       3.500   2.000   2.000  1.00  2.00           C
END
"""

  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  cg = maptbx.crystal_gridding(unit_cell=xrs.unit_cell(),
      pre_determined_n_real=(100,100,100),
      space_group_info=xrs.space_group_info())
  fc = xrs.structure_factors(d_min = 1., algorithm = "direct").f_calc()
  fft_map = miller.fft_map(crystal_gridding=cg, fourier_coefficients=fc)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  #all filled
  v, volumes = getvs(map_data, -100)
  v2, volumes2 = getvs(map_data, -100, False)
  assert v == v2 == [0, 1000000, 0]
  assert v[:2] == v2[:2] == volumes == volumes2
  # can see one blob
  v, volumes = getvs(map_data, 5)
  assert v[0]+v[1]+v[2] == 1000000
  assert v[2] == 0
  assert v[:2] == volumes
  # can see separate, approx equal volume bloobs
  v, volumes = getvs(map_data, 10)
  assert v[0]+v[1]+v[2] == 1000000
  assert abs(v[1] - v[2]) < 5
  assert v == volumes
  # nothing to see
  v, volumes = getvs(map_data, 1000)
  assert v == [1000000, 0, 0]
  assert v[:1] == volumes

def exercise4():
  cmap = flex.double(flex.grid(100,100,100))
  cmap.fill(1)
  for i in range(10,20):
    for j in range(10,20):
      for k in range(10,20):
        cmap[i,j,k] = 10
  v, volumes = getvs(cmap, 5)
  v2, volumes2 = getvs(cmap, 5, False)
  assert v == v2 == [999000, 1000, 0]
  assert v[:2] == volumes == volumes2
  #print "all filled"
  v, volumes = getvs(cmap, -5)
  v2, volumes2 = getvs(cmap, -5, False)
  assert v == v2
  assert volumes == volumes2
  assert v == [0,1000000,0]
  assert v[:2] == volumes
  #print "none filled"
  v, volumes = getvs(cmap, 20)
  v2, volumes2 = getvs(cmap, 20, False)
  assert v == v2
  assert volumes == volumes2
  assert v == [1000000,0,0]
  assert v[:1] == volumes

def exercise41():
  cmap = flex.int(flex.grid(100,100,100))
  cmap.fill(1)
  for i in range(10,20):
    for j in range(10,20):
      for k in range(10,20):
        cmap[i,j,k] = 10
  v, volumes = getvs(cmap, 5)
  assert v == [999000, 1000, 0]
  assert v[:2] == volumes
  #print "all filled"
  v, volumes = getvs(cmap, -5)
  assert v == [0,1000000,0]
  assert v[:2] == volumes
  #print "none filled"
  v, volumes = getvs(cmap, 20)
  assert v == [1000000,0,0]
  assert v[:1] == volumes


def exercise5():
  #print "corner blob"
  cmap = flex.double(flex.grid(100,100,100))
  cmap.fill(0)
  for i in range(100):
    for j in range(100):
      for k in range(100):
        if (i<10 or i>=90) and (j<10 or j>=90) and (k<10 or k>=90):
          cmap[i,j,k] = 10
          #print i,j,k
  v, volumes = getvs(cmap, 5)
  assert v == [992000, 8000, 0]
  assert v[:2] == volumes

  # test wrapping = false - borders are not transparent
  v, volumes = getvs(cmap, 5, False)
  assert volumes == [992000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]

  #print "2 blobs"
  cmap.fill(0)
  for i in range(100):
    for j in range(100):
      for k in range(100):
        if (5<i<10) and (5<j<10) and (5<k<10):
          cmap[i,j,k] = 10
        if (15<i<20) and (15<j<20) and (15<k<20):
          cmap[i,j,k] = 20
  v, volumes = getvs(cmap, 5)
  assert v == [999872,64,64]
  assert v == volumes
  v, volumes = getvs(cmap, 15)
  assert v == [999936, 64,0]
  assert v[:2] == volumes

  #print "endless blob"
  cmap.fill(0)
  for j in range(100):
    for k in range(100):
      cmap[5,j,k] = 10
  v, volumes = getvs(cmap, 5)
  assert v == [990000, 10000, 0]
  assert v[:2] == volumes

def exercise6():
  cmap = flex.double(flex.grid(100,100,100))
  #print "corner touch"
  cmap.fill(0)
  cmap[1,1,1] = cmap[2,2,2] = 10
  v, volumes = getvs(cmap, 5)
  assert v == [999998, 1, 1]
  assert v == volumes
  #print "edges touch"
  cmap.fill(0)
  cmap[1,1,1] = cmap[2,2,1] = 10
  v, volumes = getvs(cmap, 5)
  assert v == [999998, 1, 1]
  assert v == volumes
  #print "face touch"
  cmap.fill(0)
  cmap[1,1,1] = cmap[2,1,1] = 10
  v, volumes = getvs(cmap, 5)
  assert v == [999998, 2, 0]
  assert v[:2] == volumes

def exercise_volume_cutoff():
  cmap = flex.double(flex.grid(100,100,100))
  cmap.fill(0)
  for i in range(100):
    for j in range(100):
      for k in range(100):
        if (5<i<10) and (5<j<10) and (5<k<10):
          cmap[i,j,k] = 10
        if (15<i<25) and (15<j<25) and (15<k<25):
          cmap[i,j,k] = 20

  co = maptbx.connectivity(map_data=cmap, threshold=5)
  map_result = co.result()
  volumes = list(co.regions())
  #print volumes
  #[999207, 64, 729]
  vol_mask = co.volume_cutoff_mask(volume_cutoff=10)
  assert (vol_mask==1).count(True) == 793
  assert (vol_mask==0).count(True) == 999207
  vol_mask = co.volume_cutoff_mask(volume_cutoff=100)
  assert (vol_mask==1).count(True) == 729
  assert (vol_mask==0).count(True) == 999271
  vol_mask = co.volume_cutoff_mask(volume_cutoff=1000)
  assert (vol_mask==1).count(True) == 0
  assert (vol_mask==0).count(True) == 1000000

def exercise_max_values():
  cmap = flex.double(flex.grid(100,100,100))
  cmap.fill(0)
  for i in range(100):
    for j in range(100):
      for k in range(100):
        if (5<i<10) and (5<j<10) and (5<k<10):
          cmap[i,j,k] = 10
        if (15<i<25) and (15<j<25) and (15<k<25):
          cmap[i,j,k] = 20

  cmap[7,7,7] = 15
  cmap[20,20,20] = 25
  co = maptbx.connectivity(map_data=cmap, threshold=5)
  m_coors = list(co.maximum_coors())
  m_vals = list(co.maximum_values())
  vols = list(co.regions())
  assert len(m_coors) == len(m_vals) == len(vols)
  assert m_coors == [(0, 0, 0), (7, 7, 7), (20, 20, 20)]
  assert m_vals == [0.0, 15.0, 25.0]

def debug_printing(co):
  print("volumes    :",  list(co.regions()))
  print("values     :",  list(co.maximum_values()))
  print("coordinates:",  list(co.maximum_coors()))
  print("============")


def exercise_noise_elimination_two_cutoffs():
  # Purpose: eliminate noise.
  # We want to delete small blobs from the map. On the particular contouring
  # (cutoff) level we can set a threshold for volume and say: all blobs that
  # have volume less than threshold value should be deleted.
  # One more point is that we want to delete them with their 'root', meaning
  # that we are lowering threshold level and put zeros on that bigger regions.
  # But we are zeroing only those which are not merged with big good blobs.
  # Everything under second contouring level also will be zero.
  # ======================
  # From another point of view.
  # We know some threshold value for volume of good blobs on t1 contouring
  # level. We want to keep only them and clear out everything else. But the
  # keeping and clearing should be done at lower t2 contouring level.
  #
  # The result (res_mask) is 3d integer array sized as original map.
  # res_mask contain 0 for noise, 1 for valuable information.
  # Mask corresponding to t2 contouring level.
  #
  # The option "zero_all_interblob_region" by default is True, and this means
  # that everything below threshold on t2 level will be 0. If
  # zero_all_interblob_region=False then everything below threshold on t2
  # level will be 1.
  #
  #map preparation for test
  cmap = flex.double(flex.grid(100,2,2))
  cmap.fill(10)
  for i in range(10,40):
    cmap[i,1,1] = i
  for i,v in zip(range(40,60), range(40,20,-1)):
    cmap[i,1,1] = v
  for i,v in zip(range(60,70), range(20,30)):
    cmap[i,1,1] = v
  for i,v in zip(range(70,90), range(30,10,-1)):
    cmap[i,1,1] = v
  #for i in range(100):
  #  print "%d   : %d" % (i,  cmap[i,1,1])

  co1 = maptbx.connectivity(map_data=cmap, threshold=25)
  co2 = maptbx.connectivity(map_data=cmap, threshold=22)
  co3 = maptbx.connectivity(map_data=cmap, threshold=18)

  # Example 1. We have one good blob (volume>12) and one bad (volume < 12).
  # After lowering contour level they are still separate, so we want to keep
  # only big first blob, which has volume=35 on t2 contour level.
  # Here is actual call to get a mask.
  res_mask = co2.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=12,
      zero_all_interblob_region=True)
  assert (res_mask!=0).count(True) == 35
  # 2 good ===> 2 separate
  res_mask = co2.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=8)
  assert (res_mask!=0).count(True) == 50
  # 1 good, 1 bad ===> 1 big
  res_mask = co3.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=12)
  assert (res_mask!=0).count(True) == 63
  # 2 good ===> 1 big
  res_mask = co3.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=8)
  assert (res_mask!=0).count(True) == 63
  # 2 bad ===> 1 big
  res_mask = co3.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=30)
  assert (res_mask!=0).count(True) == 0

  # extreme case: nothing above t1 ==> result: everything is 0 on the mask
  co1 = maptbx.connectivity(map_data=cmap, threshold=40)
  co2 = maptbx.connectivity(map_data=cmap, threshold=22)
  res_mask = co2.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=10)
  assert (res_mask!=0).count(True) == 0

  # extreme case: everything above t1 ==> result is undefined.

  # =================================================================
  # same as above, but zero_all_interblob_region = False
  # In the first test we have 1 good blob and one bad blob. Bad one
  # will have volume=15 on t2 contouring level so we want to have 385 non-zeros
  # on resulting mask
  co1 = maptbx.connectivity(map_data=cmap, threshold=25)
  co2 = maptbx.connectivity(map_data=cmap, threshold=22)
  co3 = maptbx.connectivity(map_data=cmap, threshold=18)
  res_mask = co2.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=12,
      zero_all_interblob_region=False)
  #for i in range(100):
  #  print "%d   : %d | %d" % (i,  cmap[i,1,1], res_mask[i,1,1])
  assert (res_mask!=0).count(True) == 385

  # 2 good ===> 2 separate
  res_mask = co2.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=8,
      zero_all_interblob_region=False)
  assert (res_mask==1).count(True) == 400
  # 1 good, 1 bad ===> 1 big
  res_mask = co3.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=12,
      zero_all_interblob_region=False)
  assert (res_mask!=0).count(True) == 400
  # 2 good ===> 1 big
  res_mask = co3.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=8,
      zero_all_interblob_region=False)
  assert (res_mask!=0).count(True) == 400
  # 2 bad ===> 1 big
  res_mask = co3.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=30,
      zero_all_interblob_region=False)
  assert (res_mask!=0).count(True) == 337

  # extreme case: nothing above t1, something above t2 ==> result:
  # everything between blobs on t2 will be 1.
  co1 = maptbx.connectivity(map_data=cmap, threshold=40)
  co2 = maptbx.connectivity(map_data=cmap, threshold=22)
  res_mask = co2.noise_elimination_two_cutoffs(
      connectivity_object_at_t1=co1,
      elimination_volume_threshold_at_t1=10,
      zero_all_interblob_region=False)
  assert (res_mask!=0).count(True) == 350

def exercise_get_blobs_boundaries():
  cmap = flex.double(flex.grid(100,100,100))
  cmap.fill(1)
  for i in range(10,20):
    for j in range(10,20):
      for k in range(10,20):
        cmap[i,j,k] = 10
  co = maptbx.connectivity(map_data=cmap, threshold=5)
  # raw function:
  boundaries = co.get_blobs_boundaries()
  # how to use this:
  # boundaries[min/max, n_blob, x/y/z]
  blob_0_min_boundaries = \
      (boundaries[0,0,0], boundaries[0,0,1], boundaries[0,0,1])
  blob_0_max_boundaries = \
      (boundaries[1,0,0], boundaries[1,0,1], boundaries[1,0,1])
  # 0th blob - under the limit, covering almost whole cell
  assert blob_0_min_boundaries == (0,0,0)
  assert blob_0_max_boundaries == (99,99,99)
  # 1st blob - covers coordinates from 10 to 19 by construction
  blob_1_min_boundaries = \
      (boundaries[0,1,0], boundaries[0,1,1], boundaries[0,1,1])
  blob_1_max_boundaries = \
      (boundaries[1,1,0], boundaries[1,1,1], boundaries[1,1,1])
  assert blob_1_min_boundaries == (10,10,10)
  assert blob_1_max_boundaries == (19,19,19)
  # convinient get_blobs_boundaries_tuples
  minb, maxb = co.get_blobs_boundaries_tuples()
  assert minb == [(0,0,0), (10,10,10)]
  assert maxb == [(99,99,99), (19,19,19)]

  # ==============================
  # two blobs test
  # just add a blob to the previous cmap
  for i in range(50,70):
    for j in range(50,80):
      for k in range(50,90):
        cmap[i,j,k] = 10
  co = maptbx.connectivity(map_data=cmap, threshold=5)
  minb, maxb = co.get_blobs_boundaries_tuples()
  assert minb == [(0,0,0), (10,10,10), (50,50,50)]
  assert maxb == [(99,99,99), (19,19,19), (69,79,89)]

def exercise_expand_mask():
  # case 1: standard
  cmap = flex.double(flex.grid(30,30,30))
  cmap.fill(1)
  for i in range(10,20):
    for j in range(10,20):
      for k in range(10,20):
        cmap[i,j,k] = 10
  co = maptbx.connectivity(map_data=cmap, threshold=5)
  new_mask = co.expand_mask(id_to_expand=1, expand_size=1)
  for i in range(30):
    for j in range(30):
      for k in range(30):
        assert new_mask[i,j,k] == (i in range(9,21) and
            j in range(9,21) and k in range(9,21))

  # case 2: over boundaries
  cmap = flex.double(flex.grid(30,30,30))
  cmap.fill(1)
  cmap[1,1,1] = 10
  co = maptbx.connectivity(map_data=cmap, threshold=5)
  new_mask = co.expand_mask(id_to_expand=1, expand_size=2)
  for i in range(30):
    for j in range(30):
      for k in range(30):
        assert new_mask[i,j,k] == (i in [29,0,1,2,3] and
            j in [29,0,1,2,3] and k in [29,0,1,2,3])

def exercise_wrapping():
  cmap = flex.double(flex.grid(30,30,30))
  cmap.fill(1)
  for i in range(0,5):
    for j in range(0,5):
      for k in range(0,5):
        cmap[i,j,k] = 10
  for i in range(0,5):
    for j in range(25,30):
      for k in range(0,5):
        cmap[i,j,k] = 10
  for i in range(0,5):
    for j in range(0,5):
      for k in range(25,30):
        cmap[i,j,k] = 10

  for i in range(25,30):
    for j in range(0,5):
      for k in range(0,5):
        cmap[i,j,k] = 10
  for i in range(25,30):
    for j in range(25,30):
      for k in range(0,5):
        cmap[i,j,k] = 10
  for i in range(25,30):
    for j in range(0,5):
      for k in range(25,30):
        cmap[i,j,k] = 10

  n_in_blob = cmap.count(10)
  co = maptbx.connectivity(map_data=cmap, threshold=5, wrapping=True)
  dres = co.result().as_double()
  regs = list(co.regions())
  minb, maxb = co.get_blobs_boundaries_tuples()
  assert n_in_blob == 750
  assert regs == [26250, 750]

def exercise_preprocess_against_shallow():
  # case 1: simple
  cmap = flex.double(flex.grid(30,30,30))
  cmap.fill(1)
  for i in range(10,20):
    for j in range(10,20):
      for k in range(10,20):
        cmap[i,j,k] = 10
  for i in range(10,20):
    cmap[i,5,5] = 10
  co = maptbx.connectivity(map_data=cmap, threshold=5)
  minb, maxb = co.get_blobs_boundaries_tuples()
  assert minb == [(0, 0, 0), (10, 5, 5), (10, 10, 10)]
  assert maxb == [(29, 29, 29), (19, 5, 5), (19, 19, 19)]
  co = maptbx.connectivity(map_data=cmap, threshold=5, preprocess_against_shallow=True)
  minb, maxb = co.get_blobs_boundaries_tuples()
  assert minb == [(0, 0, 0), (10, 10, 10)]
  assert maxb == [(29, 29, 29), (19, 19, 19)] # note dissapearance of (10,5,5)(19,5,5)
  # check new map values
  for i in range(10,20):
    assert approx_equal(cmap[i,5,5], 4)


  # case 2: wrapping
  cmap = flex.double(flex.grid(30,30,30))
  cmap.fill(1)
  for i in range(10,20):
    for j in range(10,20):
      for k in range(10,20):
        cmap[i,j,k] = 10
  for i in range(10,20):
    for j in range(10,20):
      cmap[i,j,0] = 10
      cmap[i,j,29] = 10
  # standard, no wrap, 4 regions
  co = maptbx.connectivity(map_data=cmap, threshold=5, wrapping=False, preprocess_against_shallow=False)
  minb, maxb = co.get_blobs_boundaries_tuples()
  assert minb == [(0, 0, 0), (10, 10, 0), (10, 10, 10), (10, 10, 29)]
  assert maxb == [(29, 29, 29), (19, 19, 0), (19, 19, 19), (19, 19, 29)]
  # 2 small regions merged
  co = maptbx.connectivity(map_data=cmap, threshold=5, wrapping=True, preprocess_against_shallow=False)
  minb, maxb = co.get_blobs_boundaries_tuples()
  assert minb == [(0, 0, 0), (10, 10, 0), (10, 10, 10)]
  assert maxb == [(29, 29, 29), (19, 19, 29), (19, 19, 19)]
  # with wrapping the region preserved
  co = maptbx.connectivity(map_data=cmap, threshold=5, wrapping=True, preprocess_against_shallow=True)
  minb, maxb = co.get_blobs_boundaries_tuples()
  assert minb == [(0, 0, 0), (10, 10, 0), (10, 10, 10)]
  assert maxb == [(29, 29, 29), (19, 19, 29), (19, 19, 19)]

  # without wrapping - no
  co = maptbx.connectivity(map_data=cmap, threshold=5, wrapping=False, preprocess_against_shallow=True)
  minb, maxb = co.get_blobs_boundaries_tuples()
  assert minb == [(0, 0, 0), (10, 10, 10)]
  assert maxb ==[(29, 29, 29), (19, 19, 19)]

  # case 3: blob has a spike that needs to be 'shaved off'
  cmap = flex.double(flex.grid(30,30,30))
  cmap.fill(1)
  for i in range(10,20):
    for j in range(10,20):
      for k in range(10,20):
        cmap[i,j,k] = 10
  for i in range(0,10):
    cmap[i,15,15] = 10

  co = maptbx.connectivity(map_data=cmap, threshold=5, preprocess_against_shallow=False)
  volumes = list(co.regions())
  assert volumes == [25990, 1010]
  co = maptbx.connectivity(map_data=cmap, threshold=5, preprocess_against_shallow=True)
  volumes = list(co.regions())
  assert volumes == [26000, 1000]
  for i in range(0,10):
    assert approx_equal(cmap[i,15,15], 4)

  # case 4: need at least two passes
  cmap = flex.double(flex.grid(30,30,30))
  cmap.fill(1)
  cmap[10,10,10] = 10
  cmap[9,10,10] = 10
  cmap[11,10,10] = 10
  cmap[10,9,10] = 10
  cmap[10,11,10] = 10
  cmap[10,10,9] = 10
  cmap[10,10,11] = 10
  co = maptbx.connectivity(map_data=cmap, threshold=5, preprocess_against_shallow=False)
  volumes = list(co.regions())
  assert volumes == [26993, 7]
  co = maptbx.connectivity(map_data=cmap, threshold=5, preprocess_against_shallow=True)
  volumes = list(co.regions())
  assert volumes == [27000]
  assert co.preprocessing_changed_voxels == 7
  assert co.preprocessing_n_passes == 3

def write_ccp4_map(fname, unit_cell, space_group, map_data):
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
      file_name=fname,
      unit_cell=unit_cell,
      space_group=space_group,
      #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
      #gridding_last=n_real,  # This causes a bug (map gets shifted)
      map_data=map_data,
      labels=flex.std_string([""]))


def exercise_symmetry_related_regions():
  pdb_str="""
CRYST1   10.000  10.000   10.000  90.00  90.00  90.00 P 4
HETATM    1  C    C      1       2.000   2.000   2.000  1.00 20.00           C
HETATM    2  C    C      2       4.000   4.000   4.000  1.00 20.00           C
END
"""

  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  # xrs.show_summary()
  d_min = 1.
  fc = xrs.structure_factors(d_min=d_min).f_calc()
  symmetry_flags = maptbx.use_space_group_symmetry
  fftmap = fc.fft_map(symmetry_flags = symmetry_flags)
  rmup = fftmap.real_map_unpadded()
  # print ('rmup size', rmup.accessor().focus())
  # This produces 4 separate blobs
  co = maptbx.connectivity(
      map_data=rmup,
      threshold=400,
      wrapping=False,
      preprocess_against_shallow=False)
  original_regions = list(co.regions())
  # print ('regions', original_regions)
  assert len(original_regions) == 5
  beg_mask = co.result()
  # dv_mask = co.volume_cutoff_mask(0).as_double() ???
  # write_ccp4_map('volume_mask_1000.ccp4', fc.unit_cell(), fc.space_group(), dv_mask)

  co.merge_symmetry_related_regions(space_group=xrs.space_group())
  new_mask = co.result()
  assert beg_mask.count(0) == new_mask.count(0)
  assert beg_mask.count(1) + beg_mask.count(3) == new_mask.count(1)
  assert beg_mask.count(2) + beg_mask.count(4) == new_mask.count(2)
  assert sum(original_regions[1:]) == sum(original_regions[1:])


  new_regions = list(co.regions())
  assert len(new_regions) == 3
  assert list(co.maximum_values()) == []
  assert list(co.maximum_coors()) == []

  # ======================================================================
  # At this threshold 2 carbons merge. But one of the blob is cutted,
  # therefore producing 3 separate regions in unit cell
  co = maptbx.connectivity(
      map_data=rmup,
      # threshold=1000,
      threshold=1.1,
      wrapping=False,
      preprocess_against_shallow=True)
  original_regions = list(co.regions())
  assert len(original_regions) == 4
  # print ('regions', original_regions)
  beg_mask = co.result()
  # Particular numbers here seem to be platform-dependent
  # These should work on Mac
  # assert beg_mask.count(0) == 29019
  # assert beg_mask.count(1) == 1885
  # assert beg_mask.count(2) == 1714
  # assert beg_mask.count(3) == 150

  # assert original_regions == [29019, 1885, 1714, 150]
  co.merge_symmetry_related_regions(space_group=xrs.space_group())
  new_mask = co.result()
  # assert new_mask.count(0) == 29019
  # assert new_mask.count(1) == 3749
  assert beg_mask.count(0) == new_mask.count(0)
  assert beg_mask.count(1) + beg_mask.count(2) + beg_mask.count(3) == new_mask.count(1)

  new_regions = list(co.regions())
  assert len(new_regions) == 2
  # print('new regs', new_regions)
  # assert new_regions == [29019, 3749]
  assert list(co.maximum_values()) == []
  assert list(co.maximum_coors()) == []


def exercise_work_in_asu():
  pdb_str="""
CRYST1   10.000  10.000   10.000  90.00  90.00  90.00  P 4
HETATM    1  C    C      1       2.000   2.000   2.000  1.00 20.00           C
HETATM    2  C    C      2       4.000   4.000   4.000  1.00 20.00           C
END
"""

  from time import time
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  # xrs.show_summary()
  d_min = 1
  fc = xrs.structure_factors(d_min=d_min).f_calc()
  symmetry_flags = maptbx.use_space_group_symmetry
  fftmap = fc.fft_map(symmetry_flags = symmetry_flags)
  # rmup = fftmap.real_map_unpadded()
  rm = fftmap.real_map().deep_copy()
  maptbx.unpad_in_place(rm)
  mmm = rm.as_1d().min_max_mean()
  print (mmm.min, mmm.max, mmm.mean)
  # rmup = fftmap.real_map_unpadded()
  # print (dir(rm))
  print ("full size:", fftmap.real_map().accessor().focus())
  print(rm[0,0,0])
  # print (type(rm))
  # print (dir(rm))
  # STOP()
  # print(rmup[0,0,0])
  amap0  = asymmetric_map(xrs.space_group().type(), rm)
  # print(dir(amap0))
  mmm = amap0.data().as_1d().min_max_mean()
  print (mmm.min, mmm.max, mmm.mean)
  amap_data = amap0.data()
  write_ccp4_map('amap.ccp4', xrs.unit_cell(), xrs.space_group(), amap_data)
  write_ccp4_map('rm.ccp4', xrs.unit_cell(), xrs.space_group(), rm)
  # for i in range(50):
  #   print(i, amap_data[i,0,0])
  exp_map = amap0.symmetry_expanded_map()
  print(exp_map[0,0,0])
  # for i in range(32):
  #   for j in range(32):
  #     for k in range(32):
  #       assert approx_equal(rm[i,j,k], exp_map[i,j,k])

  # print(dir(amap0))
  # STOP()
  # This produces 2 separate blobs
  sg = xrs.space_group()
  print (dir(sg))
  print (sg.all_ops())
  print (sg.info())
  print ("amap0 size:", amap0.data().accessor().focus())
  # STOP()
  print (type(amap0.data()))
  threshold = 0.
  preprocess_against_shallow = True
  print ('threshold:', threshold)
  print ('preprocess_against_shallow', preprocess_against_shallow)
  t0 = time()
  co_amap = maptbx.connectivity(
      map_data=amap0.data(),
      # threshold=threshold,
      # space_group=xrs.space_group(),
      # uc_dimensions=exp_map.accessor().focus(),
      # wrapping=False,
      preprocess_against_shallow=preprocess_against_shallow)
  t1 = time()
  print ('amap time:', t1-t0)
  original_regions = list(co_amap.regions())
  print ('start regions:', original_regions)
  print ('max coords', list(co_amap.maximum_coors()))
  print ('max vals', list(co_amap.maximum_values()))

  # print(dir(exp_map))
  print(type(exp_map))
  print ("exp_map size:", exp_map.accessor().focus())
  t0 = time()
  co_full = maptbx.connectivity(
      map_data=rm,
      threshold=threshold,
      wrapping=False,
      preprocess_against_shallow=preprocess_against_shallow)
  t1 = time()
  print ('full time:', t1-t0)

  original_regions = list(co_full.regions())
  print ('start regions:', original_regions)
  print ('max coords', list(co_full.maximum_coors()))
  print ('max vals', list(co_full.maximum_values()))



  # STOP()
  # co.experiment_with_symmetry(
  #     space_group=xrs.space_group(),
  #     uc_dims=exp_map.accessor().focus())


  co_full.merge_symmetry_related_regions(
      space_group=xrs.space_group(),
      uc_dims=exp_map.accessor().focus())
  new_regions = list(co_full.regions())
  print ('new regions:', new_regions)
  print ('max coords', list(co_full.maximum_coors()))
  print ('max vals', list(co_full.maximum_values()))

if __name__ == "__main__":
  t0 = time.time()
  exercise1()  # examples of usage are here!
  exercise3()
  exercise4()
  exercise41()
  exercise5()
  exercise6()
  exercise_volume_cutoff()
  exercise_max_values()
  exercise_noise_elimination_two_cutoffs() # example and comment
  exercise_get_blobs_boundaries()
  exercise_expand_mask()
  exercise_wrapping()
  exercise_preprocess_against_shallow()
  exercise_symmetry_related_regions()
  # exercise_work_in_asu()
  print("OK time =%8.3f"%(time.time() - t0))
