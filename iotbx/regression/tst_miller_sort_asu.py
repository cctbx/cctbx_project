
from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx import miller

def exercise_sort():
  expected_unsorted_data="""    4 I 4
  1  0  0  0  1  0  0  0  1
  0  0  0
  0 -1  0  1  0  0  0  0  1
  0  0  0
  0  1  0 -1  0  0  0  0  1
  0  0  0
 -1  0  0  0 -1  0  0  0  1
  0  0  0
   0   0   4   0   0   4     0 1 0  1 80331.8  8648.0
   0   0  12   0   0  12     0 1 0  1104526.7 11623.3
   0   0  -4   0   0   4     0 2 0  1 44883.6  4527.6
   0   0  -8   0   0   8     0 2 0  1 41134.1  4431.9
   0   0   8   0   0   8     0 1 0  1 50401.8  5464.6
   0   0   8   0   0   8     0 1 0  1 53386.0  5564.1
   0   0  -8   0   0   8     0 2 0  1119801.4 12231.2
   0   0  12   0   0  12     0 1 0  1105312.9 12602.2
   0   0  16   0   0  16     0 1 0  1 14877.6  2161.5
"""
  expected_sorted_data="""    4 I 4
  1  0  0  0  1  0  0  0  1
  0  0  0
  0 -1  0  1  0  0  0  0  1
  0  0  0
  0  1  0 -1  0  0  0  0  1
  0  0  0
 -1  0  0  0 -1  0  0  0  1
  0  0  0
   0   0   4   0   0   4     0 1 0  1 80331.8  8648.0
   0   0  -4   0   0   4     0 2 0  1 44883.6  4527.6
   0   0   8   0   0   8     0 1 0  1 50401.8  5464.6
   0   0   8   0   0   8     0 1 0  1 53386.0  5564.1
   0   0  -8   0   0   8     0 2 0  1 41134.1  4431.9
   0   0  -8   0   0   8     0 2 0  1119801.4 12231.2
   0   0  12   0   0  12     0 1 0  1104526.7 11623.3
   0   0  12   0   0  12     0 1 0  1105312.9 12602.2
   0   0  16   0   0  16     0 1 0  1 14877.6  2161.5
"""

  from cctbx.xray import observation_types
  from cctbx.array_family import flex
  xs = crystal.symmetry((113.949,113.949,32.474,90.000,90.000,90.000), "I4")
  mi = flex.miller_index((
  (0,0,4),
  (0,0,12),
  (0,0,-4),
  (0,0,-8),
  (0,0,8),
  (0,0,8),
  (0,0,-8),
  (0,0,12),
  (0,0,16),
   ))
  data = flex.double(( 80331.8, 104526.7, 44883.6, 41134.1, 50401.8, 53386.0, 119801.4, 105312.9, 14877.6,))
  sigmas = flex.double((8648.0, 11623.3, 4527.6, 4431.9, 5464.6, 5564.1, 12231.2, 12602.2, 2161.5,))
  ms = miller.set(xs, mi,anomalous_flag=True)
  i_obs = miller.array(ms, data=data, sigmas=sigmas).set_observation_type(
    observation_types.intensity() )

  i_obs.export_as_scalepack_unmerged(file_name='unsorted.sca',
     scale_intensities_for_scalepack_merge=True)
  i_obs=i_obs.sort(by_value="asu_indices")
  i_obs.export_as_scalepack_unmerged(file_name='sorted.sca',
     scale_intensities_for_scalepack_merge=True)
  with open('unsorted.sca') as f:
    unsorted = f.read()
  with open('sorted.sca') as f:
    sorted = f.read()
  assert unsorted==expected_unsorted_data
  assert sorted==expected_sorted_data

def exercise():
  exercise_sort()
  print("OK")

if __name__ == "__main__" :
  exercise()
#---end
