from __future__ import division
import os
import sys
import libtbx.load_env
from libtbx.utils import Sorry
from iotbx import reflection_file_reader
from libtbx import easy_run
from libtbx.test_utils import Exception_expected, open_tmp_file
from iotbx.file_reader import any_file, sort_by_file_type, group_files

def exercise_sort() :
  unmerged_data = """    4 i4            
  1  0  0  0  1  0  0  0  1
  0  0  0
 -1  0  0  0 -1  0  0  0  1
  0  0  0
  0 -1  0  1  0  0  0  0  1
  0  0  0
  0  1  0 -1  0  0  0  0  1
  0  0  0
   0   0   4   0   0   4   188 1 1  1 80331.8  8648.0
   0   0  12   0   0  12    22 1 0  1119801.4 12231.2
   0   0  -4   0   0   4     8 2 1  1104526.7 11623.3
   0   0  -8   0   0   8   198 2 0  1 50401.8  5464.6
   0   0   8   0   0   8    19 1 0  1 44883.6  4527.6
   0   0   8   0   0   8   184 1 1  1 41134.1  4431.9
   0   0  -8   0   0   8     4 2 1  1 53386.0  5564.1
   0   0  12   0   0  12   180 1 1  1105312.9 12602.2
   0   0  16   0   0  16    26 1 0  1 14877.6  2161.5
"""

  merged_data = """    1
 -987
   113.949   113.949    32.474    90.000    90.000    90.000 i4
  65   2   1  6044.5   555.3  5871.7   534.6
  64  11   1   438.1   251.3   602.3   247.0
  64  10   2  7917.3   815.7  7919.8   696.0
"""

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
   0   0  12   0   0  12     0 1 0  1119801.4 12231.2
   0   0  -4   0   0   4     0 2 0  1104526.7 11623.3
   0   0  -8   0   0   8     0 2 0  1 50401.8  5464.6
   0   0   8   0   0   8     0 1 0  1 44883.6  4527.6
   0   0   8   0   0   8     0 1 0  1 41134.1  4431.9
   0   0  -8   0   0   8     0 2 0  1 53386.0  5564.1
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
   0   0  -4   0   0   4     0 2 0  1104526.7 11623.3
   0   0   8   0   0   8     0 1 0  1 44883.6  4527.6
   0   0   8   0   0   8     0 1 0  1 41134.1  4431.9
   0   0  -8   0   0   8     0 2 0  1 50401.8  5464.6
   0   0  -8   0   0   8     0 2 0  1 53386.0  5564.1
   0   0  12   0   0  12     0 1 0  1119801.4 12231.2
   0   0  12   0   0  12     0 1 0  1105312.9 12602.2
   0   0  16   0   0  16     0 1 0  1 14877.6  2161.5
"""
  f1 = open("merged.sca","w")
  f1.write(merged_data)
  f1.close()
  merged=reflection_file_reader.any_reflection_file("merged.sca")
  crystal_symmetry=merged.as_miller_arrays()[0].crystal_symmetry()

  f = open("unmerged.sca", "w")
  f.write(unmerged_data)
  f.close()
  unmerged=reflection_file_reader.any_reflection_file("unmerged.sca")

  arrays=unmerged.as_miller_arrays()
  i_obs=unmerged.as_miller_arrays(merge_equivalents=False)[0].customized_copy(
      crystal_symmetry=crystal_symmetry)
  
  i_obs.export_as_scalepack_unmerged(file_name='unsorted.sca',
     scale_intensities_for_scalepack_merge=True)
  i_obs=i_obs.sort(by_value="asu_indices")
  i_obs.export_as_scalepack_unmerged(file_name='sorted.sca',
     scale_intensities_for_scalepack_merge=True)
  unsorted=open('unsorted.sca').read()
  sorted=open('sorted.sca').read()
  assert unsorted==expected_unsorted_data
  assert sorted==expected_sorted_data

def exercise () :
  exercise_sort()
  print "OK"

if __name__ == "__main__" :
  exercise()
#---end
