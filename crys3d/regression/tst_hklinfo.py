# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.test_utils import contains_substring
import os, subprocess

os.environ['PYTHONIOENCODING'] = 'UTF-8'

test1_cifstr = """
data_r4jglsf
#
_audit.revision_id    1_0
_audit.creation_date  2013-04-03
_audit.update_record
;The data loops contain the following:
 structure factors used in refinement, from crystal 185601 that is a Met-Inhibition crystal.
 unmerged original index intensities of data set used for refinement and phasing, this data set comes from crystal 185601 that is a Met-Inhibition crystal, wavelength 1.
 unmerged original index intensities of data set used for phasing, this data set comes from crystal 185601 that is a Met-Inhibition crystal, wavelength 2.
 unmerged original index intensities of data set used for phasing, this data set comes from crystal 185601 that is a Met-Inhibition crystal, wavelength 3.
 experimental phases from crystal 185601 that is a Met-Inhibition crystal.
 density modified experimental phases from crystal 185601 that is a Met-Inhibition crystal.
;
#
_cell.entry_id      4jgl
_cell.length_a      79.322
_cell.length_b      79.322
_cell.length_c      50.113
_cell.angle_alpha   90.000
_cell.angle_beta    90.000
_cell.angle_gamma   120.000
#
loop_
_diffrn_radiation_wavelength.id
_diffrn_radiation_wavelength.wavelength
1 0.918401
2 0.979338
3 0.979108
#
_entry.id   4jgl
#
_exptl_crystal.id   1
#
_reflns_scale.group_code   1
#
_symmetry.entry_id               4jgl
_symmetry.space_group_name_H-M   'P 61'
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
_refln.F_meas_sigma_au
_refln.F_calc
_refln.phase_calc
1 1 1    0    0    3 o   1055.8   38.2     1482.3     7.0
1 1 1    0    0   12 x        ?      ?     1790.1   149.3
1 1 1    3   32    3 o    204.3    3.6      199.8   328.7
1 1 1    3   32    4 o     37.8    5.6       30.8    26.8
1 1 1    3   32    5 o    216.6    3.9      184.1   304.7
1 1 1    3   32    6 o     54.2    3.7       65.5   335.2
1 1 1    3   32    7 f     85.2    3.0      103.4    80.7
1 1 1    3   32    8 o    138.0    3.5      118.6   305.6
1 1 1    3   32    9 o    139.2    3.3      122.3   302.1
1 1 1    3   32   10 o     33.0    5.5       21.9   357.2
1 1 1    3   32   11 o    102.9    3.5      105.7   129.7
1 1 1    3   32   12 o     57.9    5.2       52.7   182.7
1 1 1   54    1    0 f     20.6    9.6       28.9   180.0
1 1 1   54    1    4 o     15.6    5.4       19.4   331.4
1 1 1   54    1    5 o     25.2    6.6       16.9    69.0
#END
data_r4jglAsf
#
_cell.entry_id      4jgl
_cell.length_a      79.3216
_cell.length_b      79.3216
_cell.length_c      50.1126
_cell.angle_alpha   90.0000
_cell.angle_beta    90.0000
_cell.angle_gamma   120.0000
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details
;   unmerged original index intensities of data set used for refinement  and phasing.
;
#
_entry.id   4jgl
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
_refln.intensity_meas
_refln.intensity_sigma
1 1 1    0    0    2  o       5.3     11.4
1 1 1    0    0   -2  o      13.4     13.3
1 1 1    0    0    3  o      -5.9      9.6
1 1 1    0    0   -3  o      14.8     14.3
1 1 1    2   -2  -36  o    1056.4    235.7
1 1 1    2    0  -36  o    1176.5    218.3
1 1 1    0   -2   36  o     788.1    195.6
1 1 1    2   -2   36  o     984.7    196.5
1 1 1   -1   55    5  o     285.9    179.5
1 1 1    1  -55   -5  o    -105.2    188.9
1 1 1  -55   54    5  o     558.6    268.4
1 1 1   55  -54   -5  o     340.6    311.3
1 1 1  -54   -1    5  o     153.0    129.0
1 1 1  -54   -1    5  o      -8.7    101.3
1 1 1    1  -55    5  o      52.0    249.1
1 1 1   -1   55   -5  o     -80.1    130.8
1 1 1   55  -54    5  o     576.9    272.9
1 1 1  -55   54   -5  o     219.8    205.6
#END
data_r4jglBsf
#
_cell.entry_id      4jgl
_cell.length_a      79.3132
_cell.length_b      79.3132
_cell.length_c      50.1025
_cell.angle_alpha   90.0000
_cell.angle_beta    90.0000
_cell.angle_gamma   120.0000
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details
'   unmerged original index intensities of data set used for phasing.'
#
_entry.id   4jgl
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
_refln.intensity_meas
_refln.intensity_sigma
1 2 1    0    0    2  o       2.5      7.6
1 2 1    0    0   -2  o       0.5     11.0
1 2 1    0    0    3  o      20.5     17.2
1 2 1   -1    0  -12  o   22476.8   2327.6
1 2 1    0   -1  -12  o   22730.3   1962.0
1 2 1    1   -1  -12  o   23595.6   1782.0
1 2 1   -1    0   12  o   24913.4   1960.4
1 2 1    1    0  -12  o   23227.8   1777.6
1 2 1   -1    1  -12  o   22906.0   2033.0
1 2 1    0   -1  -14  o   10141.0    877.2
1 2 1   -1    1   14  o    8471.0    838.7
1 2 1    1   -1  -14  o   10248.7    748.4
1 2 1   -1    0   14  o    7835.0    786.6
1 2 1    1  -55    4  o     -43.1     74.1
1 2 1    1  -55    5  o     105.6     60.3
1 2 1   55  -54    5  o     303.0    276.2
1 2 1  -55   54   -5  o     -44.3    127.2
#END
data_r4jglCsf
#
_cell.entry_id      4jgl
_cell.length_a      79.4003
_cell.length_b      79.4003
_cell.length_c      50.1616
_cell.angle_alpha   90.0000
_cell.angle_beta    90.0000
_cell.angle_gamma   120.0000
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details
'   unmerged original index intensities of data set used for phasing.'
#
_entry.id   4jgl
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
_refln.intensity_meas
_refln.intensity_sigma
1 3 1    0    0    2  o       2.2     14.2
1 3 1    0    0   -2  o      12.7     14.3
1 3 1    0    0    3  o      -9.5      3.8
1 3 1    0    0   -3  o      -8.2      3.9
1 3 1    0    0    4  o      20.2     22.1
1 3 1   11  -11   20  o    7381.3    861.9
1 3 1  -11   11  -20  o    8143.7    797.9
1 3 1   11    0   20  o    7939.2    822.9
1 3 1   11    0   20  o    8401.4    889.9
1 3 1  -49   -2    5  o    -682.5    697.3
1 3 1   50    1    1  o     164.3    164.9
1 3 1  -50   -1   -1  o    1390.9    464.9
1 3 1   -1   51    1  o    -589.3    274.5
1 3 1    1  -51   -1  o     903.4    414.8
1 3 1  -51   50    1  o      63.5     67.2
1 3 1  -50   -1    1  o   -1052.0    570.9
1 3 1   50    1   -1  o     174.0    158.8
1 3 1    1  -51    1  o     198.8    479.5
1 3 1   -1   51   -1  o    -259.0    239.2
#END
data_r4jglDsf
#
_cell.entry_id      4jgl
_cell.length_a      79.3130
_cell.length_b      79.3130
_cell.length_c      50.1030
_cell.angle_alpha   90.0000
_cell.angle_beta    90.0000
_cell.angle_gamma   120.0000
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details             '  experimental phases.'
#
_entry.id   4jgl
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
_refln.pdbx_HL_A_iso
_refln.pdbx_HL_B_iso
_refln.pdbx_HL_C_iso
_refln.pdbx_HL_D_iso
1 1 1    0    0    3  o    0.28   -1.43    0.00    0.00
1 1 1    3    0   12  o    0.02    0.24    0.00    0.00
1 1 1    0    0   18  o    2.89    4.23    0.73   -1.26
1 1 1    0    0   24  o   -0.08   -4.39    0.86    1.74
1 1 1    0    0   30  o   23.29  -37.82   13.58    1.94
1 1 1    0    0   36  o    3.35    3.12    0.29   -0.96
1 1 1    0    1    0  o       ?       ?       ?       ?
1 1 1    0    1    1  o       ?       ?       ?       ?
1 1 1   12   11    7  o   -2.50    7.76    1.03    0.31
1 1 1   12   11    8  o  -50.91  153.61   31.01   32.58
1 1 1   12   11    9  o    5.73   -8.35    0.09    2.76
1 1 1   12   11   10  o  -27.65   30.85   -1.02    9.22
1 1 1   53    3    3  o   -0.31    0.46    0.00    0.00
1 1 1   53    3    4  o    0.32   -0.06    0.00    0.00
1 1 1   54    1    0  o   -0.27    0.00    0.00    0.00
1 1 1   54    1    1  o    0.47   -0.01    0.00    0.00
1 1 1   54    1    2  o   -0.18   -0.34   -0.01   -0.02
1 1 1   54    1    3  o   -0.02   -0.28    0.00    0.00
1 1 1   54    1    4  o    0.35    0.49   -0.01   -0.01
1 1 1   54    1    5  o    0.15   -0.48    0.00    0.00
#END
data_r4jglEsf
#
_cell.entry_id      4jgl
_cell.length_a      79.3130
_cell.length_b      79.3130
_cell.length_c      50.1030
_cell.angle_alpha   90.0000
_cell.angle_beta    90.0000
_cell.angle_gamma   120.0000
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details             '  density modified experimental phases.'
#
_entry.id   4jgl
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
_refln.phase_meas
_refln.pdbx_FWT
_refln.fom
1 1 1    0    0    3  o   41.29 1238.431  0.99
1 1 1    0    0   12  o  125.51 1995.770  1.00
1 1 1    0    0   18  o   99.07 807.187  0.98
1 1 1    0    0   24  o -122.68 491.014  0.97
1 1 1    0    0   30  o  -25.97 161.558  0.99
1 1 1    0    0   36  o   31.40  56.578  0.81
1 1 1   53    2    1  o -126.70  22.452  0.60
1 1 1   53    2    2  o   92.10  17.818  0.61
1 1 1   53    2    3  o  -93.89   7.406  0.27
1 1 1   53    2    4  o -168.56  10.492  0.43
1 1 1   53    3    3  o  127.04  10.290  0.35
1 1 1   53    3    4  o   -0.28   4.982  0.18
1 1 1   54    1    0  o  180.00   7.752  0.33
1 1 1   54    1    1  o   13.94   5.584  0.21
1 1 1   54    1    2  o -111.70   4.379  0.16
1 1 1   54    1    3  o  -79.85   4.164  0.15
1 1 1   54    1    4  o   56.34   8.051  0.32
1 1 1   54    1    5  o  -69.69   7.166  0.29
#END OF REFLECTIONS

"""


phil1str = """
selected_info {
  labels = False
  description = True
  wavelength = True
  n_reflections = True
  span = False
  minmax_data = True
  minmax_sigmas = False
  d_minmax = True
  data_sigdata_max=True
  unit_cell = True
  space_group = True
  n_centrics = False
  is_anomalous = False
  is_symmetry_unique = False
  n_sys_abs = False
  data_completeness = True
  data_compl_infty = False
  ano_completeness = False
  ano_mean_diff = False
  n_bijvoet = False
  n_singletons = False
}
wrap_labels = -20
"""


expected1str = """
Starting job
===============================================================================
33 Miller arrays in this dataset:
       Type      |   λ/Å   |  #HKLs  |     min,max data       |MaxDatSigDat|  d_min,d_max/Å   |     unit cell (a/Å, b/Å, c/Å, α%s, β%s, γ%s)      |   space group      |Data compl.|
         Integer |0.918401 |      15 |        1.0,         1.0|        nan |    1.25,   16.704| 79.322, 79.322, 50.113,     90,     90,    120 |     P 61 (No. 169) |0.00029017 |
         Integer |0.918401 |      15 |        1.0,         1.0|        nan |    1.25,   16.704| 79.322, 79.322, 50.113,     90,     90,    120 |     P 61 (No. 169) |0.00029017 |
         Integer |0.918401 |      15 |        1.0,         1.0|        nan |    1.25,   16.704| 79.322, 79.322, 50.113,     90,     90,    120 |     P 61 (No. 169) |0.00029017 |
          String |0.918401 |      15 |        nan,         nan|        nan |    1.25,   16.704| 79.322, 79.322, 50.113,     90,     90,    120 |     P 61 (No. 169) |0.00029017 |
       Amplitude |0.918401 |      14 |       15.6,      1055.8|      56.75 |    1.25,   16.704| 79.322, 79.322, 50.113,     90,     90,    120 |     P 61 (No. 169) |0.00034303 |
      Map coeffs |0.918401 |      15 |       16.9,      1790.1|        nan |    1.25,   16.704| 79.322, 79.322, 50.113,     90,     90,    120 |     P 61 (No. 169) |0.00029017 |
         Integer |0.918401 |      18 |        1.0,         1.0|        nan |    1.25,   25.056| 79.322, 79.322, 50.113,     90,     90,    120 |                  ? |       nan |
         Integer |0.918401 |      18 |        1.0,         1.0|        nan |    1.25,   25.056| 79.322, 79.322, 50.113,     90,     90,    120 |                  ? |       nan |
         Integer |0.918401 |      18 |        1.0,         1.0|        nan |    1.25,   25.056| 79.322, 79.322, 50.113,     90,     90,    120 |                  ? |       nan |
          String |0.918401 |      18 |        nan,         nan|        nan |    1.25,   25.056| 79.322, 79.322, 50.113,     90,     90,    120 |                  ? |       nan |
       Intensity |0.918401 |      18 |     -105.2,      1176.5|      5.389 |    1.25,   25.056| 79.322, 79.322, 50.113,     90,     90,    120 |                  ? |       nan |
         Integer |0.979338 |      17 |        1.0,         1.0|        nan |    1.25,   25.051| 79.313, 79.313, 50.102,     90,     90,    120 |                  ? |       nan |
         Integer |0.979338 |      17 |        2.0,         2.0|        nan |    1.25,   25.051| 79.313, 79.313, 50.102,     90,     90,    120 |                  ? |       nan |
         Integer |0.979338 |      17 |        1.0,         1.0|        nan |    1.25,   25.051| 79.313, 79.313, 50.102,     90,     90,    120 |                  ? |       nan |
          String |0.979338 |      17 |        nan,         nan|        nan |    1.25,   25.051| 79.313, 79.313, 50.102,     90,     90,    120 |                  ? |       nan |
       Intensity |0.979338 |      17 |      -44.3,  2.4913e+04|      13.69 |    1.25,   25.051| 79.313, 79.313, 50.102,     90,     90,    120 |                  ? |       nan |
         Integer |0.979108 |      19 |        1.0,         1.0|        nan |   1.361,   25.081|   79.4,   79.4, 50.162,     90,     90,    120 |                  ? |       nan |
         Integer |0.979108 |      19 |        3.0,         3.0|        nan |   1.361,   25.081|   79.4,   79.4, 50.162,     90,     90,    120 |                  ? |       nan |
         Integer |0.979108 |      19 |        1.0,         1.0|        nan |   1.361,   25.081|   79.4,   79.4, 50.162,     90,     90,    120 |                  ? |       nan |
          String |0.979108 |      19 |        nan,         nan|        nan |   1.361,   25.081|   79.4,   79.4, 50.162,     90,     90,    120 |                  ? |       nan |
       Intensity |0.979108 |      19 |    -1052.0,      8401.4|      10.21 |   1.361,   25.081|   79.4,   79.4, 50.162,     90,     90,    120 |                  ? |       nan |
         Integer |0.918401 |      20 |        1.0,         1.0|        nan |    1.25,   68.687| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |
         Integer |0.918401 |      20 |        1.0,         1.0|        nan |    1.25,   68.687| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |
         Integer |0.918401 |      20 |        1.0,         1.0|        nan |    1.25,   68.687| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |
          String |0.918401 |      20 |        nan,         nan|        nan |    1.25,   68.687| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |
       HL coeffs |0.918401 |      18 |      -1.87,      166.29|        nan |    1.25,   16.701| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |
         Integer |0.918401 |      18 |        1.0,         1.0|        nan |    1.25,   16.701| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |
         Integer |0.918401 |      18 |        1.0,         1.0|        nan |    1.25,   16.701| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |
         Integer |0.918401 |      18 |        1.0,         1.0|        nan |    1.25,   16.701| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |
          String |0.918401 |      18 |        nan,         nan|        nan |    1.25,   16.701| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |
  Floating-point |0.918401 |      18 |    -168.56,       180.0|        nan |    1.25,   16.701| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |
  Floating-point |0.918401 |      18 |      4.164,      1995.8|        nan |    1.25,   16.701| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |
  Floating-point |0.918401 |      18 |       0.15,         1.0|        nan |    1.25,   16.701| 79.313, 79.313, 50.103,     90,     90,    120 |                  ? |       nan |

===============================================================================
Job complete
""" %("\u00b0", "\u00b0", "\u00b0")



test2_cifstr = """
data_mytest
#
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
  1  x,y,z
  2  -y,x-y,z+1/3
  3  -x+y,-x,z+2/3
#
_diffrn.id                  1
_diffrn.crystal_id          1
#
_diffrn_radiation_wavelength.id           1
_diffrn_radiation_wavelength.wavelength   0.997432
#
_reflns_scale.group_code   1
#
_space_group.crystal_system       trigonal
_space_group.IT_number            144
_space_group.name_H-M_alt         'P 31'
_space_group.name_Hall            ' P 31'
_symmetry.space_group_name_H-M    'P 31'
_symmetry.space_group_name_Hall   ' P 31'
_symmetry.Int_Tables_number       144
_cell.length_a                    50.000
_cell.length_b                    50.000
_cell.length_c                    40.000
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 120.000
_cell.volume                      86602.540
loop_
  _refln.crystal_id
  _refln.wavelength_id
  _refln.scale_group_code
  _refln.index_h
  _refln.index_k
  _refln.index_l
  _refln.I
  _refln.SIGI
  _refln.MyMap
  _refln.PHIMyMap
  _refln.Oink
  _refln.bleep
  _refln.FOM
  _refln.A_calc
  _refln.B_calc
1 1 1   1  -2   3   11.205  13.695    1   180    24.9  12.429     0    2.3   -3.2
1 1 1   0   0  -4    6.353   6.353    3   120  14.521       ?   0.1   12.4    4.7
1 1 1   1   2   3   26.167  24.921  0.8    90       ?  -3.328  0.25   -1.5   -4.8
1 1 1   0   1   2    14.94   6.225    2    60   3.738   3.738  0.35    7.9   10.2
1 1 1   1   0   2     2.42  11.193    4     0  22.429     4.9   0.4    8.3   22.5
1 1 1  -1   1  -2   24.921  26.167    1   -60  28.635       ?   0.5      ?      ?
1 1 1   2  -2  -2   16.185   8.715    1   -90       ?  -5.521   0.6   -9.4    8.6
1 1 1  -2   1   0   11.798   4.538    5  -120   3.328  10.738  0.75   -1.4    6.3
1 1 1   1   0  -2   21.183  27.413    ?     ?   3.738       ?     ?    9.1    7.7
1 1 1   0   0   2     4.98  21.165    6   135   19.92   19.92     1    3.9   -2.1
1 1 1  -1  -2   3    456.5    36.3    ?     ?       ?       ?     ?   -3.1    5.9
1 1 1   0   0   4   654.36   9.123    ?     ?       ?       ?     ?    6.9   13.8
1 1 1   1   2  -3  -78.234   76.37    ?     ?       ?       ?     ?      ?      ?
1 1 1   0  -1   2   369.78   11.29    ?     ?       ?       ?     ?    0.7   -3.2
1 1 1  -1   0   2  672.899   56.14    ?     ?       ?       ?     ?    7.8    0.2
1 1 1  -1   1   2   90.316    78.9    ?     ?       ?       ?     ?    2.8   11.5
1 1 1   0   0  -3        ?       ?    ?     ?       ?  38.635     ?    7.1    3.8

"""

phil2str = """
selected_info {
  labels = True
  description = True
  wavelength = False
  n_reflections = True
  span = True
  minmax_data = True
  minmax_sigmas = False
  data_sigdata_max = False
  d_minmax = True
  unit_cell = False
  space_group = True
  n_centrics = False
  is_anomalous = True
  is_symmetry_unique = False
  n_sys_abs = False
  data_completeness = True
  data_compl_infty = False
  ano_completeness = False
  ano_mean_diff = False
  n_bijvoet = False
  n_singletons = False
}
wrap_labels = -20
"""

expected2str ="""
Starting job
===============================================================================
9 Miller arrays in this dataset:
 Labels          |       Type      |  #HKLs  |               Span              |     min,max data       |  d_min,d_max/Å   |   space group      |Anomalous|Data compl.|
mytest,_refln.crystal_id
                 |         Integer |      17 |         (-2, -2, -4), (2, 2, 4) |        1.0,         1.0|    10.0,     25.0|     P 31 (No. 144) |    True |       nan |
mytest,_refln.wavelength_id
                 |         Integer |      17 |         (-2, -2, -4), (2, 2, 4) |        1.0,         1.0|    10.0,     25.0|     P 31 (No. 144) |    True |       nan |
mytest,_refln.scale_group_code
                 |         Integer |      17 |         (-2, -2, -4), (2, 2, 4) |        1.0,         1.0|    10.0,     25.0|     P 31 (No. 144) |    True |       nan |
mytest,_refln.I,_refln.SIGI
                 |  Floating-point |      16 |         (-2, -2, -4), (2, 2, 4) |    -78.234,       672.9|    10.0,     25.0|     P 31 (No. 144) |    True |       nan |
mytest,_refln.MyMap,_refln.PHIMyMap
                 |      Map coeffs |       9 |         (-2, -2, -4), (2, 2, 3) |        0.8,         6.0|    10.0,     25.0|     P 31 (No. 144) |    True |  0.067308 |
mytest,_refln.Oink
                 |  Floating-point |       8 |         (-2, -2, -4), (1, 1, 3) |      3.328,      28.635|    10.0,     25.0|     P 31 (No. 144) |    True |       nan |
mytest,_refln.bleep
                 |  Floating-point |       8 |         (-2, -2, -3), (2, 2, 3) |     -5.521,      38.635|  10.337,     25.0|     P 31 (No. 144) |   False |   0.13462 |
mytest,_refln.FOM
                 |  Floating-point |       9 |         (-2, -2, -4), (2, 2, 3) |        0.0,         1.0|    10.0,     25.0|     P 31 (No. 144) |    True |  0.067308 |
mytest,_refln.A_calc,_refln.B_calc
                 |      Map coeffs |      15 |         (-2, -2, -4), (2, 2, 4) |     3.2757,      23.982|    10.0,     25.0|     P 31 (No. 144) |    True |       nan |

===============================================================================
Job complete
"""



testhklstr = """!FORMAT=XDS_ASCII    MERGE=FALSE    FRIEDEL'S_LAW=TRUE
!OUTPUT_FILE=XDS_ASCII.HKL        DATE=10-Aug-2012
!Generated by CORRECT   (VERSION  March 15, 2012)
!PROFILE_FITTING= TRUE
!NAME_TEMPLATE_OF_DATA_FRAMES=../../img/x4_1_?????.cbf CBF
!DATA_RANGE=       1     360
!ROTATION_AXIS=  0.999993 -0.000046 -0.003827
!OSCILLATION_RANGE=  1.000000
!STARTING_ANGLE=     0.000
!STARTING_FRAME=       1
!INCLUDE_RESOLUTION_RANGE=    50.000     0.798
!SPACE_GROUP_NUMBER=    1
!UNIT_CELL_CONSTANTS=     7.484    13.661    21.131  94.941  96.598  99.633
!UNIT_CELL_A-AXIS=     4.886    -4.416     3.556
!UNIT_CELL_B-AXIS=    -8.390    -9.505    -5.088
!UNIT_CELL_C-AXIS=    11.220     2.222   -17.768
!REFLECTING_RANGE_E.S.D.=     0.837
!BEAM_DIVERGENCE_E.S.D.=     0.136
!X-RAY_WAVELENGTH=  0.700010
!INCIDENT_BEAM_DIRECTION=  0.003855  0.003610  1.428541
!FRACTION_OF_POLARIZATION=   0.990
!POLARIZATION_PLANE_NORMAL=  0.000000  1.000000  0.000000
!AIR=  0.000107
!SILICON=  1.347949
!SENSOR_THICKNESS=  0.320000
!DETECTOR=PILATUS
!OVERLOAD=   1048500
!DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
!DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
!DETECTOR_DISTANCE=   170.465
!ORGX=   1248.63  ORGY=   1252.18
!NX=  2463  NY=  2527    QX=  0.172000  QY=  0.172000
!VARIANCE_MODEL=  1.494E+00  1.180E-03
!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=12
!ITEM_H=1
!ITEM_K=2
!ITEM_L=3
!ITEM_IOBS=4
!ITEM_SIGMA(IOBS)=5
!ITEM_XD=6
!ITEM_YD=7
!ITEM_ZD=8
!ITEM_RLP=9
!ITEM_PEAK=10
!ITEM_CORR=11
!ITEM_PSI=12
!END_OF_HEADER
     0     0    -1  4.691E+02  2.124E+01  1232.6  1282.1    271.8 0.02777 100  30   16.89
     0     0    -1  4.685E+02  2.084E+01  1232.7  1226.9     94.1 0.02777 100  27 -164.84
     0     0     1  4.884E+02  2.164E+01  1270.1  1281.9     91.8 0.02777  99  22  163.53
     0     0     1  5.193E+02  2.330E+01  1270.0  1227.1    274.1 0.02777 100  25  -15.58
     0     0    -2  4.858E+04  2.040E+03  1213.9  1199.3     95.3 0.05557 100  20 -165.49
     0     0    -2  4.852E+04  2.038E+03  1213.8  1309.6    270.6 0.05553 100  16   17.54
     0     0     2  4.915E+04  2.064E+03  1288.9  1309.5     90.7 0.05557 100  16  162.88
     0     0     2  4.800E+04  2.017E+03  1288.7  1199.6    275.3 0.05558 100  17  -14.92
     0     0    -3  2.103E+04  8.846E+02  1195.0  1171.6     96.5 0.08343 100  24 -166.14
     0     0    -3  2.041E+04  8.589E+02  1194.8  1337.3    269.4 0.08328 100  27   18.20
"""


phil3str = """
selected_info {
  labels = True
  description = True
  wavelength = False
  n_reflections = True
  span = False
  minmax_data = True
  minmax_sigmas = False
  d_minmax = True
  unit_cell = False
  space_group = True
  n_centrics = False
  is_anomalous = True
  is_symmetry_unique = False
  n_sys_abs = False
  data_completeness = True
  data_compl_infty = True
  ano_completeness = True
  ano_mean_diff = False
  n_bijvoet = True
  n_singletons = True
}
merge_equivalents = True
wrap_labels = -20
"""


expected3str = """
Starting job
===============================================================================
1 Miller arrays in this dataset:
 Labels          |       Type      |  #HKLs  |     min,max data       |  d_min,d_max/Å   |   space group      |Anomalous|Data compl.|Compl.inf.|Ano.complete|#Bijvoets|#Singletons|
 iobs,sigma_iobs |       Intensity |       5 |     468.79,  4.8562e+04|   6.956,   20.869|        P 1 (No. 1) |    True |   0.22727 |  0.22727 |    0.18182 |       2 |         1 |

===============================================================================
Job complete
"""


testscastr = """    1
 -987
    44.151    44.151    44.151    80.000    80.000    80.000 h3
   1   0   0    -0.2     0.4
   2  -1   0     0.2     0.6
   2  -1   2  2471.3   119.7
   2   0   1  3166.3   237.6  2602.5   113.7
   3   0   1 10052.2   496.7 10107.3   504.2
   3   0   2 15486.2   469.7 15548.0   590.2
   3   0   3  4730.2   167.5  4454.7   188.8
   3   0   4  9840.6   543.6  8136.2   375.4
   4  -3   1 23731.2   553.5 20271.4   584.9
   4  -3   2  1311.5    57.1  1342.2    49.4
   4  -3   3   556.5    22.2  1231.4    56.0
   4  -3   4 17415.2   642.7 21023.1  1031.7
   4  -3   5 27394.1  1103.3 29405.3  2926.8
   4  -3   6  1462.1    75.0  1507.1   170.9
   4  -2   0  3763.1   166.4
   4  -2   1 15626.8   900.8 15210.5   706.0
   4  -2   2  5633.4   164.3  6092.0   186.0
   4  -2   3   236.7    12.3   428.2    17.7
   4  -2   4   837.0    40.5   762.2    45.4
   5  -3   2  6569.5   193.7  7542.3   254.6
   5  -3   3  6283.6   190.9  7280.4   241.6
   5  -3   4   567.3    36.4   949.3    46.8
   5  -3   5    77.0    43.2    26.9    58.7
   5  -3   6  3170.0   185.7  3134.1   175.3
   5  -3   7  1859.2   124.2  1936.3   158.1
   5  -3   8  3124.1   253.3  3047.0   216.0
   5  -3   9  2296.7   138.1
   5  -3  18   -11.8    53.3
   5  -2   0 10279.5   312.4
   6  -5   1  7179.4   200.2  7433.3   194.9
   6  -5   2  8662.2   388.1 10568.4   365.1
   6  -5   3 15441.5   598.0 14331.0   551.3
"""

phil4str = """
selected_info {
  labels = True
  description = True
  wavelength = False
  n_reflections = True
  span = False
  minmax_data = True
  minmax_sigmas = True
  d_minmax = True
  unit_cell = False
  space_group = True
  data_sigdata=True
  n_centrics = False
  is_anomalous = True
  is_symmetry_unique = False
  n_sys_abs = False
  data_completeness = True
  data_compl_infty = False
  ano_completeness = True
  ano_mean_diff = False
  n_bijvoet = True
  n_singletons = True
}
wrap_labels = -10
"""

expected4str = """
Starting job
===============================================================================
1 Miller arrays in this dataset:
 Labels          |       Type      |  #HKLs  |     min,max data       |     min,max sigmas     | DatSigDat|  d_min,d_max/Å   |   space group      |Anomalous|Data compl.|Ano.complete|#Bijvoets|#Singletons|
I(+),SIGI(+),I(-),SIGI(-)
                 |       Intensity |      57 |      -11.8,  2.9405e+04|        0.4,      2926.8|    21.45 |   2.293,   43.002|   R 3 :R (No. 146) |    True | 0.0059511 |  0.0052203 |      25 |         7 |

===============================================================================
Job complete
"""


def test_hklinfo_run(cifstr, filename, philstr, expectedstr):
  with open(filename, "w") as fh:
    fh.write(cifstr)
  with open("philinput.txt", "w") as fh:
    fh.write(philstr)

  obj = subprocess.Popen("cctbx.HKLinfo %s philinput.txt" %filename,
                         shell=True,
                         stdin = subprocess.PIPE,
                         stdout = subprocess.PIPE,
                         stderr = subprocess.STDOUT)
  out,err = obj.communicate()
  fullstdout = out.decode().replace("\r\n", "\n") # omit \r\n line endings on Windows
  #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
  assert (contains_substring( fullstdout, expectedstr ) )


# =============================================================================
if (__name__ == '__main__'):
  test_hklinfo_run(test1_cifstr, "test-sf.cif", phil1str, expected1str)
  test_hklinfo_run(test2_cifstr, "test-sf.cif", phil2str, expected2str)
  test_hklinfo_run(testhklstr, "test.hkl", phil3str, expected3str)
  test_hklinfo_run(testscastr, "test.sca", phil4str, expected4str)


  print('OK')
