# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx import easy_run
import iotbx.cif
from libtbx.test_utils import approx_equal, show_diff
from six.moves import cStringIO as StringIO
import iotbx.mtz
import os

sf_5r82 = """\
data_r5r82sf
#
_audit.revision_id     1_0
_audit.creation_date   2020-03-11
_audit.update_record   "Initial release"
#
_cell.entry_id      5r82
_cell.length_a      112.665
_cell.length_b      52.848
_cell.length_c      44.468
_cell.angle_alpha   90.000
_cell.angle_beta    102.966
_cell.angle_gamma   90.000
#
_diffrn.id        1
_diffrn.details   "data from final refinement with ligand, final.mtz"
_diffrn.crystal_id 1
#
_diffrn_radiation_wavelength.id           1
_diffrn_radiation_wavelength.wavelength   0.9126
#
_entry.id   5r82
#
_exptl_crystal.id   1
#
_reflns_scale.group_code   1
#
_symmetry.entry_id               5r82
_symmetry.space_group_name_H-M   "C 1 2 1"
_symmetry.Int_Tables_number      5
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
_refln.F_calc_au
_refln.phase_calc
_refln.pdbx_FWT
_refln.pdbx_PHWT
_refln.pdbx_DELFWT
_refln.pdbx_DELPHWT
_refln.fom
1 1 1 -86 0  8  o ?       ?     14.57   0.00   8.89    0.00   0.00    0.00   0.00
1 1 1 -85 1  3  o ?       ?     74.69   159.82 45.89   159.75 0.00    0.00   0.00
1 1 1 -85 1  4  o ?       ?     18.89   136.18 11.62   136.53 0.00    0.00   0.00
1 1 1 -85 1  5  o ?       ?     40.54   44.81  24.89   44.97  0.00    0.00   0.00
1 1 1 -85 1  6  o ?       ?     14.55   265.49 8.91    267.15 0.00    0.00   0.00
1 1 1 -85 1  7  o ?       ?     38.21   15.80  23.59   15.42  0.00    0.00   0.00
1 1 1 -85 1  8  o ?       ?     47.35   181.29 29.24   181.19 0.00    0.00   0.00
1 1 1 -85 1  9  o 56.21   26.85 38.50   89.75  41.26   89.77  17.50   89.77  0.58
1 1 1 -85 1  10 o 65.09   31.91 27.40   104.11 36.22   104.55 19.20   104.55 0.41
1 1 1 -85 1  11 o 54.28   26.34 24.33   27.77  29.23   27.86  14.13   27.86  0.41
1 1 1 -85 1  12 o 58.89   28.38 21.26   158.09 27.96   158.29 14.93   158.29 0.35
1 1 1 -85 3  4  o ?       ?     16.70   93.70  10.28   93.61  0.00    0.00   0.00
1 1 1 -85 3  5  o ?       ?     18.82   319.51 11.62   319.61 0.00    0.00   0.00
1 1 1 -85 3  6  o ?       ?     9.64    143.57 6.02    143.97 0.00    0.00   0.00
1 1 1 -85 3  7  o ?       ?     34.35   260.46 21.18   260.07 0.00    0.00   0.00
1 1 1 -85 3  8  o ?       ?     30.57   8.43   18.77   8.70   0.00    0.00   0.00
1 1 1 -85 3  9  o 66.81   31.32 64.85   81.01  57.84   80.86  17.89   80.86  0.73
1 1 1 -85 3  10 o 59.74   28.33 64.04   132.52 49.17   132.36 9.69    132.36 0.74
1 1 1 -85 3  11 o 50.02   24.64 5.34    93.51  6.70    95.13  3.34    95.13  0.10
1 1 1 -85 5  5  o ?       ?     29.69   339.80 18.22   340.06 0.00    0.00   0.00
1 1 1 -85 5  6  o ?       ?     52.33   264.97 32.14   265.05 0.00    0.00   0.00
1 1 1 -85 5  7  o ?       ?     63.54   232.41 39.02   232.26 0.00    0.00   0.00
1 1 1 -85 5  8  o ?       ?     25.64   55.65  15.80   55.84  0.00    0.00   0.00
1 1 1 -85 5  9  o 44.64   22.72 27.71   201.58 23.51   202.12 6.57    202.12 0.45
1 1 1 -85 5  10 o ?       ?     48.38   118.56 29.58   118.56 0.00    0.00   0.00
1 1 1 -84 0  1  o ?       ?     3.94    180.00 2.34    180.00 0.00    0.00   0.00
1 1 1 -84 0  2  o ?       ?     9.64    180.00 5.88    180.00 0.00    0.00   0.00
1 1 1 -84 0  3  o ?       ?     57.42   180.00 35.50   180.00 0.00    0.00   0.00
1 1 1 -84 0  4  o ?       ?     5.23    0.00   3.15    0.00   0.00    0.00   0.00
1 1 1 -84 0  5  o ?       ?     3.86    180.00 2.39    180.00 0.00    0.00   0.00
1 1 1 -84 0  6  o ?       ?     21.93   0.00   13.65   0.00   0.00    0.00   0.00
1 1 1 -84 0  7  o ?       ?     7.45    0.00   4.61    0.00   0.00    0.00   0.00
1 1 1 -84 0  8  o 73.12   45.70 7.69    0.00   15.82   0.00   11.03   0.00   0.14
#
data_r5r82Asf
#
_cell.entry_id      5r82
_cell.length_a      112.665
_cell.length_b      52.848
_cell.length_c      44.468
_cell.angle_alpha   90.000
_cell.angle_beta    102.966
_cell.angle_gamma   90.000
#
_diffrn.id        1
_diffrn.details   "data from original reflections, data.mtz"
_diffrn.crystal_id 1
#
_diffrn_radiation_wavelength.id           1
_diffrn_radiation_wavelength.wavelength   0.9126
#
_entry.id   5r82
#
_exptl_crystal.id   1
#
_reflns_scale.group_code   1
#
_symmetry.entry_id               5r82
_symmetry.space_group_name_H-M   "C 1 2 1"
_symmetry.Int_Tables_number      5
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
1 1 1 -86 0  7  o ?      ?
1 1 1 -86 0  8  o ?      ?
1 1 1 -85 1  3  o ?      ?
1 1 1 -85 1  4  o ?      ?
1 1 1 -85 1  5  o ?      ?
1 1 1 -85 1  6  o ?      ?
1 1 1 -85 1  7  o ?      ?
1 1 1 -85 1  8  o ?      ?
1 1 1 -85 1  9  o 1.28   0.61
1 1 1 -85 1  10 o 1.47   0.72
1 1 1 -85 1  11 o 1.21   0.59
1 1 1 -85 1  12 o 1.30   0.63
1 1 1 -85 3  4  o ?      ?
1 1 1 -85 3  5  o ?      ?
1 1 1 -85 3  6  o ?      ?
1 1 1 -85 3  7  o ?      ?
1 1 1 -85 3  8  o ?      ?
1 1 1 -85 3  9  o 1.53   0.72
1 1 1 -85 3  10 o 1.35   0.64
1 1 1 -85 3  11 o 1.12   0.55
1 1 1 -85 3  12 o ?      ?
1 1 1 -85 5  5  o ?      ?
1 1 1 -85 5  6  o ?      ?
1 1 1 -85 5  7  o ?      ?
1 1 1 -85 5  8  o ?      ?
1 1 1 -85 5  9  o 1.02   0.52
1 1 1 -85 5  10 o ?      ?
1 1 1 -84 0  1  o ?      ?
1 1 1 -84 0  2  o ?      ?
1 1 1 -84 0  3  o ?      ?
1 1 1 -84 0  4  o ?      ?
1 1 1 -84 0  5  o ?      ?
1 1 1 -84 0  6  o ?      ?
1 1 1 -84 0  7  o ?      ?
1 1 1 -84 0  8  o 1.71   1.07
#
data_r5r82Bsf
#
_cell.entry_id      5r82
_cell.length_a      112.665
_cell.length_b      52.848
_cell.length_c      44.468
_cell.angle_alpha   90.000
_cell.angle_beta    102.970
_cell.angle_gamma   90.000
#
_diffrn.id        1
_diffrn.details   "data for ligand evidence map (PanDDA event map), event_map_1.mtz"
_diffrn.crystal_id 1
#
_diffrn_radiation_wavelength.id           1
_diffrn_radiation_wavelength.wavelength   0.9126
#
_entry.id   5r82
#
_exptl_crystal.id   1
#
_reflns_scale.group_code   1
#
_symmetry.entry_id               5r82
_symmetry.space_group_name_H-M   "P 1"
_symmetry.Int_Tables_number      1
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
_refln.pdbx_PHWT
1 1 1 -85 -6  7  o 183.35   1.00 42.48
1 1 1 -85 -6  8  o 204.28   1.00 118.73
1 1 1 -85 -5  5  o 408.51   1.00 26.47
1 1 1 -85 -5  6  o 368.34   1.00 -17.28
1 1 1 -85 -5  7  o 105.28   1.00 -108.00
1 1 1 -85 -5  8  o 447.04   1.00 49.71
1 1 1 -85 -5  9  o 151.27   1.00 76.05
1 1 1 -85 -5  10 o 177.73   1.00 -11.88
1 1 1 -85 -4  4  o 176.60   1.00 137.96
1 1 1 -85 -4  5  o 94.48    1.00 45.34
1 1 1 -85 -4  6  o 269.47   1.00 3.56
1 1 1 -85 -4  7  o 146.26   1.00 103.92
1 1 1 -85 -4  8  o 251.33   1.00 152.35
1 1 1 -85 -4  9  o 189.51   1.00 -158.44
1 1 1 -85 -4  10 o 340.08   1.00 50.25
1 1 1 -85 -4  11 o 85.03    1.00 42.71
1 1 1 -85 -3  4  o 458.44   1.00 -148.81
1 1 1 -85 -3  5  o 928.82   1.00 25.99
1 1 1 -85 -3  6  o 266.68   1.00 -145.70
1 1 1 -85 -3  7  o 524.04   1.00 20.59
1 1 1 -85 -3  8  o 824.38   1.00 126.20
1 1 1 -85 -3  9  o 334.69   1.00 -9.45
1 1 1 -85 -3  10 o 367.92   1.00 -34.78
1 1 1 -85 -3  11 o 622.58   1.00 74.86
1 1 1 -85 -2  3  o 197.98   1.00 -127.10
1 1 1 -85 -2  4  o 140.13   1.00 -111.56
#
#END OF REFLECTIONS
"""

sf_5r82_mtz_results = {'r5r82sf': """\
Title: phenix.cif_as_mtz
Space group symbol from file: C2
Space group number from file: 5
Space group from matrices: C 1 2 1 (No. 5)
Point group symbol from file: 2
Number of crystals: 2
Number of Miller indices: 33
Resolution range: 1.34113 1.30997
History:
Crystal 1:
  Name: HKL_base
  Project: HKL_base
  Id: 0
  Unit cell: (112.665, 52.848, 44.468, 90, 102.966, 90)
  Number of datasets: 1
  Dataset 1:
    Name: HKL_base
    Id: 0
    Wavelength: 0
    Number of columns: 0
Crystal 2:
  Name: crystal_0
  Project: project_0
  Id: 2
  Unit cell: (112.665, 52.848, 44.468, 90, 102.966, 90)
  Number of datasets: 1
  Dataset 1:
    Name: dataset
    Id: 1
    Wavelength: 0.9126
    Number of columns: 13
    label        #valid  %valid     min    max type
    H                33 100.00%  -86.00 -84.00 H: index h,k,l
    K                33 100.00%    0.00   5.00 H: index h,k,l
    L                33 100.00%    1.00  12.00 H: index h,k,l
    R-free-flags     33 100.00%    1.00   1.00 I: integer
    FOBS              9  27.27%   44.64  73.12 F: amplitude
    SIGFOBS           9  27.27%   22.72  45.70 Q: standard deviation
    FC               33 100.00%    3.86  74.69 F: amplitude
    PHIFC            33 100.00% -178.71 180.00 P: phase angle in degrees
    FWT              33 100.00%    2.34  57.84 F: amplitude
    PHWT             33 100.00% -178.81 180.00 P: phase angle in degrees
    DELFWT           33 100.00%    0.00  19.20 F: amplitude
    PHDELWT          33 100.00% -157.88 158.29 P: phase angle in degrees
    FOM              33 100.00%    0.00   0.74 F: amplitude
""",
'r5r82Asf':
"""Title: phenix.cif_as_mtz
Space group symbol from file: C2
Space group number from file: 5
Space group from matrices: C 1 2 1 (No. 5)
Point group symbol from file: 2
Number of crystals: 2
Number of Miller indices: 35
Resolution range: 1.34113 1.3096
History:
Crystal 1:
  Name: HKL_base
  Project: HKL_base
  Id: 0
  Unit cell: (112.665, 52.848, 44.468, 90, 102.966, 90)
  Number of datasets: 1
  Dataset 1:
    Name: HKL_base
    Id: 0
    Wavelength: 0
    Number of columns: 0
Crystal 2:
  Name: crystal_0
  Project: project_0
  Id: 2
  Unit cell: (112.665, 52.848, 44.468, 90, 102.966, 90)
  Number of datasets: 1
  Dataset 1:
    Name: dataset
    Id: 1
    Wavelength: 0.9126
    Number of columns: 6
    label        #valid  %valid    min    max type
    H                35 100.00% -86.00 -84.00 H: index h,k,l
    K                35 100.00%   0.00   5.00 H: index h,k,l
    L                35 100.00%   1.00  12.00 H: index h,k,l
    R-free-flags     35 100.00%   1.00   1.00 I: integer
    FOBS              9  25.71%   1.02   1.71 F: amplitude
    SIGFOBS           9  25.71%   0.52   1.07 Q: standard deviation
""", 'r5r82Bsf': """\
Title: phenix.cif_as_mtz
Space group symbol from file: P1
Space group number from file: 1
Space group from matrices: P 1 (No. 1)
Point group symbol from file: 1
Number of crystals: 2
Number of Miller indices: 26
Resolution range: 1.3216 1.31054
History:
Crystal 1:
  Name: HKL_base
  Project: HKL_base
  Id: 0
  Unit cell: (112.665, 52.848, 44.468, 90, 102.97, 90)
  Number of datasets: 1
  Dataset 1:
    Name: HKL_base
    Id: 0
    Wavelength: 0
    Number of columns: 0
Crystal 2:
  Name: crystal_0
  Project: project_0
  Id: 2
  Unit cell: (112.665, 52.848, 44.468, 90, 102.97, 90)
  Number of datasets: 1
  Dataset 1:
    Name: dataset
    Id: 1
    Wavelength: 0.9126
    Number of columns: 7
    label        #valid  %valid     min    max type
    H                26 100.00%  -85.00 -85.00 H: index h,k,l
    K                26 100.00%   -6.00  -2.00 H: index h,k,l
    L                26 100.00%    3.00  11.00 H: index h,k,l
    R-free-flags     26 100.00%    1.00   1.00 I: integer
    FOBS             26 100.00%   85.03 928.82 F: amplitude
    SIGFOBS          26 100.00%    1.00   1.00 Q: standard deviation
    PHWT             26 100.00% -158.44 152.35 P: phase angle in degrees
"""}

def tst_1(prefix="tst_split_data_cif_1"):
  with open(prefix+'.cif', 'w') as f:
    f.write(sf_5r82)
  easy_run.fully_buffered('iotbx.split_data_cif %s.cif' % prefix)
  original_model = iotbx.cif.reader(input_string=sf_5r82).model()
  block_names = original_model.keys()
  for bn in block_names:
    assert os.path.isfile('%s.cif_%s_000.mtz' % (prefix, bn))
    assert os.path.isfile('%s.cif_%s_000.cif' % (prefix, bn))

  # testing output cif
  for bn in block_names:
    m = iotbx.cif.reader(file_path='%s.cif_%s_000.cif' % (prefix, bn)).model()
    for k in original_model[bn].keys():
      assert approx_equal(m[bn][k], original_model[bn][k])

  # testing output mtz
  for bn in block_names:
    mtz_obj = iotbx.mtz.object('%s.cif_%s_000.mtz' % (prefix, bn))
    # print('working with', '%s.cif_%s_000.mtz' % (prefix, bn))
    strio = StringIO()
    mtz_obj.show_summary(out=strio)
    # mtz_obj.show_column_data(out=strio)
    val = strio.getvalue()
    assert not show_diff(val, sf_5r82_mtz_results[bn])






if __name__ == '__main__':
  tst_1()
  print("OK")
