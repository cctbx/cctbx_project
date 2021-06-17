
from __future__ import absolute_import, division, print_function
from iotbx.reflection_file_reader import any_reflection_file
from iotbx.command_line import export_scalepack_unmerged
from libtbx.utils import null_out

prefix = "tst_export_scalepack_unmerged"

def exercise_xds() : pass # TODO
def exercise_mtz() : pass # TODO

def exercise_scalepack():
  sca_in = """\
    4 i41
  1  0  0  0  1  0  0  0  1
  0  0  0
 -1  0  0  0 -1  0  0  0  1
  0  0  0
  0 -1  0  1  0  0  0  0  1
  0  6  3
  0  1  0 -1  0  0  0  0  1
  0  6  3
   0   0   4   0   0   4   188 1 1  1 80331.8  8648.0
   0   0  -4   0   0   4     8 2 1  1104526.7 11623.3
   0   0   8   0   0   8    19 1 0  1 44883.6  4527.6
   0   0   8   0   0   8   184 1 1  1 41134.1  4431.9
   0   0  -8   0   0   8     4 2 1  1 53386.0  5564.1
   0   0  -8   0   0   8   198 2 0  1 50401.8  5464.6
   0   0  12   0   0  12    22 1 0  1119801.4 12231.2
   0   0  12   0   0  12   180 1 1  1105312.9 12602.2
   0   0  16   0   0  16    26 1 0  1 14877.6  2161.5
   1   0   3   1   0   3    10 1 0  1 19407.0  2301.4
   1   0   3   1   0   3   184 1 1  1 19337.0  2204.0
  -1   0   3   1   0   3    19 1 0  2 19513.6  2220.3
  -1   0   3   1   0   3   193 1 1  2 20702.3  2227.8
   0  -1   3   1   0   3    11 1 0  3 20041.6  2360.2
   0  -1   3   1   0   3   185 1 1  3 20043.9  2296.8
   0   1   3   1   0   3    18 1 0  4 19677.4  2125.2
  -1   0  -3   1   0   3     4 2 1  1 17056.4  1829.7
  -1   0  -3   1   0   3   189 2 0  1 16545.2  2065.9
   1   0  -3   1   0   3    13 2 1  2 17169.5  1768.9
   1   0  -3   1   0   3   198 2 0  2 16195.8  1927.5
"""
  sca_file = prefix + "_in.sca"
  with open(sca_file, "w") as f:
    f.write(sca_in)
  ofn = export_scalepack_unmerged.run(args=[sca_file], out=null_out())[0]
  assert (ofn == "tst_export_scalepack_unmerged_in_unmerged.sca")
  hkl_in = any_reflection_file(ofn).file_content()
  i_obs = hkl_in.as_miller_arrays(merge_equivalents=False)[0]
  assert (not i_obs.is_unique_set_under_symmetry())
  assert (list(hkl_in.batch_numbers) == [188, 8, 19, 184, 4, 198, 22, 180, 26,
          10, 184, 19, 193, 11, 185, 18, 4, 189, 13, 198])

def exercise_cif():
  # .cif input
  cif_in = """\
data_r2etdsf
#
_cell.entry_id      2etd
_cell.length_a      80.4540
_cell.length_b      85.2590
_cell.length_c      53.3970
_cell.angle_alpha   90.0000
_cell.angle_beta    90.0000
_cell.angle_gamma   90.0000
#
_symmetry.entry_id               2etd
_symmetry.space_group_name_H-M   'C 2 2 2'
#
loop_
_refln.crystal_id
_refln.wavelength_id
_refln.scale_group_code
_refln.index_h
_refln.index_k
_refln.index_l
_refln.status
_refln.intensity_meas
_refln.intensity_sigma
1 1 1    0    2   -8  o    1544.6    130.4
1 1 1    0    2    9  o    3450.1    264.9
1 1 1    0   -2    9  o    3243.5    268.5
1 1 1    0   -2   -9  o    3475.4    265.0
1 1 1    0    2   -9  o    3420.1    260.3
1 1 1    0    2   10  o      58.8     31.8
1 1 1    0   -2  -10  o     131.7     50.7
1 1 1    0    2  -10  o      97.9     48.7
1 1 1    0   -2   11  o    9808.3    953.5
1 1 1    0   -2  -11  o   11486.1    970.3
1 1 1    0    2  -11  o   11278.1    967.8
1 1 1    0    2   12  o    1368.9    150.7
1 1 1    0   -2   12  o    1620.9    148.7
1 1 1    0   -2  -12  o    1293.5    147.6
1 1 1    0    2  -12  o    1619.5    155.7
1 1 1    0   -2   13  o   47438.3   4104.1
1 1 1    0    2   14  o    8577.5   1188.0
1 1 1    0   -2   14  o    7996.7   1179.9
1 1 1    0   -2  -14  o    9333.7   1178.1
1 1 1    0    2  -14  o    9642.3   1197.6
1 1 1    0    2   15  o   13577.7   1852.4
1 1 1    0   -2   15  o   14100.9   1852.8
1 1 1    0   -2  -15  o   14184.2   1871.6
1 1 1    0    2   16  o     135.6     76.3
1 1 1    0   -2   16  o     117.0     60.4
"""
  cif_file = prefix + ".cif"
  with open(cif_file, "w") as f:
    f.write(cif_in)
  hkl_orig = any_reflection_file(cif_file).file_content()
  i_obs_orig = hkl_orig.as_miller_arrays(merge_equivalents=False)[0]
  ofn = export_scalepack_unmerged.run(args=[cif_file], out=null_out())[0]
  assert (ofn == prefix + "_unmerged.sca")
  hkl_new = any_reflection_file(ofn).file_content()
  i_obs_new = hkl_new.as_miller_arrays(merge_equivalents=False)[0]
  assert (not i_obs_new.is_unique_set_under_symmetry())
  assert i_obs_new.indices().all_eq(i_obs_orig.indices())

if (__name__ == "__main__"):
  exercise_scalepack()
  exercise_cif()
  print("OK")
