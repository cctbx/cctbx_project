from __future__ import division, print_function
import time, os
from iotbx import reflection_file_reader
#from iotbx import cif
from libtbx import easy_run

# Expected to fail
# Source of error might be in cif.reader
# The miller_arrays in cif.reader(input_string=cif_string_5pti).as_miller_arrays()
# don't group _refln.F_calc,_refln.phase_calc together in the second data set


def exercise_cif_as_mtz_neutron_joint():
  #
  # applies to pdb codes 5PTI and 5RSA
  #
  cif_file_name = 'tst_cif_as_mtz_neutron.cif'
  mtz_file_name = 'tst_cif_as_mtz_neutron.mtz'
  with open(cif_file_name, 'w') as f:
    f.write(cif_string_5pti)

  cmds = [
      "phenix.cif_as_mtz",
      "%s" % cif_file_name,
      '--output_file_name=%s' % mtz_file_name]

  cmd = " ".join(cmds)
  print(cmd)
  easy_run.call(cmd)

  assert(os.path.isfile(mtz_file_name))

  miller_arrays = reflection_file_reader.any_reflection_file(
    file_name = mtz_file_name).as_miller_arrays()

#  expected_ma_labels = [
#  'FC,PHIFC',
#  'FOBS',
#  'R-free-flags',
#  'FOBS-1',
#  'FC-1,PHIC',
#  'R-free-flags-1'
#  ]

  labels = []
  for ma in miller_arrays:
    label = ma.info().label_string()
    labels.append(label)

  assert('FOBS-1' in labels)
  # label name is prob irrelevant, but main point is that there shoulb be no PHIC
  # associated to the observed array

#  miller_arrays2 = cif.reader(input_string=cif_string_5pti).as_miller_arrays()
#  for ma in miller_arrays2:
#    print(ma.info().label_string())


cif_string_5pti = """
data_r5ptisf
#
_cell.entry_id      5pti
_cell.length_a      74.100
_cell.length_b      23.400
_cell.length_c      28.900
_cell.angle_alpha   90.000
_cell.angle_beta    90.000
_cell.angle_gamma   90.000
#
_entry.id   5pti
#
_exptl_crystal.id   1
#
_reflns_scale.group_code   1
#
_symmetry.entry_id               5pti
_symmetry.space_group_name_H-M   'P 21 21 21'
#
#
loop_
_refln.crystal_id
_refln.wavelength_id
_refln.scale_group_code
_refln.index_h
_refln.index_k
_refln.index_l
_refln.status
_refln.F_meas_au
_refln.F_calc
_refln.phase_calc
1 1 1   12    0    0  o   187.00     139.97    0.00
1 1 1   14    0    0  o    84.50     127.35    0.00
1 1 1   16    0    0  o   717.50     685.54  180.00
1 1 1   18    0    0  o   253.00     229.50  180.00
1 1 1   20    0    0  o   349.50     294.53  180.00
1 1 1   22    0    0  o   193.50     222.86    0.00
1 1 1   24    0    0  o    95.50      74.90  180.00
1 1 1   26    0    0  o   309.50     333.77  180.00
1 1 1   28    0    0  o   131.00     123.72  180.00
1 1 1   30    0    0  o    96.50      63.88  180.00
1 1 1   32    0    0  o   127.50     110.84  180.00
1 1 1   34    0    0  o   198.00     230.39    0.00
1 1 1   36    0    0  o    25.00      30.19  180.00
1 1 1   38    0    0  o    97.00     115.32    0.00
1 1 1   40    0    0  o    97.00      56.27    0.00
1 1 1   42    0    0  o    14.00       2.70    0.00
1 1 1   44    0    0  o    63.00      38.24  180.00
1 1 1   46    0    0  o    83.00      83.22  180.00
1 1 1   48    0    0  o    56.50      64.21    0.00
1 1 1   50    0    0  o    53.50      59.85    0.00
1 1 1   52    0    0  o    26.50      23.11  180.00
1 1 1   54    0    0  o    22.00       7.41  180.00
1 1 1   56    0    0  o    12.50      12.48  180.00
1 1 1   58    0    0  o    16.00      13.94    0.00
1 1 1   60    0    0  o    18.50      14.35    0.00
1 1 1   64    0    0  o    16.50       7.87    0.00
1 1 1   66    0    0  o    14.50      20.83    0.00
1 1 1   10    1    0  o   209.00     268.30  180.00
1 1 1   11    1    0  o     8.00      20.85   90.00
1 1 1   12    1    0  o    34.00      26.35    0.00
1 1 1   13    1    0  o    83.50     109.88   90.00
1 1 1   14    1    0  o   531.50     542.91    0.00
#END
data_r5ptiAsf
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details             '  The following is NEUTRON DIFFRACTION DATA'
#
_entry.id   5pti
#
_exptl_crystal.id   1
#
_reflns_scale.group_code   1
#
#
loop_
_refln.crystal_id
_refln.wavelength_id
_refln.scale_group_code
_refln.index_h
_refln.index_k
_refln.index_l
_refln.status
_refln.F_meas_au
_refln.F_calc
_refln.phase_calc
1 1 1   12    0    0  o    22.33      22.25    0.00
1 1 1   14    0    0  o    13.64      10.69    0.00
1 1 1   16    0    0  o   121.84     119.95  180.00
1 1 1   18    0    0  o    14.43      14.25  180.00
1 1 1   20    0    0  o    41.56      40.00  180.00
1 1 1   22    0    0  o    10.74      13.40  180.00
1 1 1   24    0    0  o    22.12      12.54    0.00
1 1 1   26    0    0  o    58.41      62.15  180.00
1 1 1   28    0    0  o    20.79      13.98  180.00
1 1 1   30    0    0  o    25.86      24.56  180.00
1 1 1   32    0    0  o    21.82      19.66  180.00
1 1 1   34    0    0  o    38.21      36.26    0.00
1 1 1   36    0    0  o    15.54      11.88  180.00
1 1 1   38    0    0  o    13.19      13.77    0.00
1 1 1   40    0    0  o    17.90      18.77    0.00
1 1 1   10    1    0  o    59.72      62.73  180.00
1 1 1   11    1    0  o    26.68      26.27  270.00
1 1 1   12    1    0  o    16.32      29.94  180.00
1 1 1   14    1    0  o    74.44      74.78    0.00
1 1 1   16    1    0  o    54.65      43.53    0.00
1 1 1   17    1    0  o    65.02      70.74   90.00
1 1 1   20    1    0  o    19.83      15.54    0.00
1 1 1   21    1    0  o    37.57      41.73  270.00
1 1 1   22    1    0  o    17.50      21.70    0.00
1 1 1   23    1    0  o    23.32      33.85  270.00
1 1 1   24    1    0  o    12.68      10.56    0.00
1 1 1   26    1    0  o    15.99      11.97    0.00
1 1 1   29    1    0  o    14.00      11.96   90.00
1 1 1   30    1    0  o    13.34       6.28  180.00
1 1 1   31    1    0  o    27.42      25.47   90.00
1 1 1   32    1    0  o    14.05      13.24  180.00
1 1 1   33    1    0  o    14.77      15.13   90.00
1 1 1   34    1    0  o    16.72      12.70  180.00
1 1 1    9    2    0  o    20.73      16.32   90.00
1 1 1   10    2    0  o     7.28      18.44    0.00
1 1 1   11    2    0  o    34.37      45.06   90.00
1 1 1   12    2    0  o    80.47      66.92  180.00
1 1 1   13    2    0  o    64.54      61.32   90.00
1 1 1   14    2    0  o    13.21       8.68  180.00
1 1 1   15    2    0  o    21.26      19.30  270.00
1 1 1   16    2    0  o    51.45      40.70  180.00
1 1 1   17    2    0  o    43.97      39.43  270.00
1 1 1   18    2    0  o    36.28      34.01  180.00
1 1 1   19    2    0  o    31.65      35.55   90.00
1 1 1   20    2    0  o    43.51      36.46    0.00
1 1 1   22    2    0  o    27.88      34.95  360.00
#END OF REFLECTIONS
"""

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_cif_as_mtz_neutron_joint()
  print("OK. Time: %-8.4f"%(time.time()-t0))
