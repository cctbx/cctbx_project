
from __future__ import division
from __future__ import print_function
from libtbx import easy_run

def exercise():
  from mmtbx.command_line import cif_as_mtz
  from iotbx import file_reader
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
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
#
loop_
_diffrn_radiation_wavelength.id
_diffrn_radiation_wavelength.wavelength
1 1.0163
2 0.9797
3 0.9796
#_entry.id   2etd
#
_exptl_crystal.id   1
#
_reflns_scale.group_code   1
#
_symmetry.entry_id               2etd
_symmetry.space_group_name_H-M   'C 2 2 2'
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
1 1 1    0    2    8 o    149.0    3.1      558.1     0.0
1 1 1    0    2    9 o    217.9    4.3      119.3   180.0
1 1 1    0    2   10 o     33.3    5.3      340.6   180.0
1 1 1    0    2   11 o    389.0   10.0      692.8     0.0
1 1 1    0    2   12 o    143.5    3.7       54.5   180.0
1 1 1    0    2   13 f    801.4   36.0     1185.6     0.0
1 1 1    0    2   14 o    351.5   11.8      403.6   360.0
1 1 1    0    2   15 o    439.4   17.2      644.6   180.0
1 1 1    0    2   16 o     39.5    6.8      161.6   180.0
#END
data_r2etdAsf
#
_cell.entry_id      2etd
_cell.length_a      80.454
_cell.length_b      85.259
_cell.length_c      53.397
_cell.angle_alpha   90.000
_cell.angle_beta    90.000
_cell.angle_gamma   90.000
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details
;   unmerged original index intensities of data set used for refinement and phasing.
;
#
_entry.id   2etd
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
#END
data_r2etdBsf
#
_cell.entry_id      2etd
_cell.length_a      80.419
_cell.length_b      85.227
_cell.length_c      53.389
_cell.angle_alpha   90.000
_cell.angle_beta    90.000
_cell.angle_gamma   90.000
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details
'   unmerged original index intensities of data set used for phasing.'
#
_entry.id   2etd
#
_exptl_crystal.id   1
#
_reflns_scale.group_code   1
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
1 2 1    0    2    8  o    1095.4    110.5
1 2 1    0   -2    8  o    1105.1    117.9
1 2 1    0   -2   -8  o    1055.4    108.6
1 2 1    0    2   -8  o    1105.1    118.0
1 2 1    0    2    9  o    3795.1    366.5
1 2 1    0   -2    9  o    3925.3    371.7
1 2 1    0   -2   -9  o    3726.2    366.9
1 2 1    0    2   -9  o    3899.3    371.8
1 2 1    0    2   10  o      27.0     40.7
1 2 1    0   -2   10  o      80.7     34.8
1 2 1    0   -2  -10  o      18.3     40.3
1 2 1    0    2  -10  o      16.5     45.4
1 2 1    0    2   11  o   10830.4    977.6
1 2 1    0   -2   11  o    9647.1    972.1
1 2 1    0   -2  -11  o   10854.6    983.4
1 2 1    0    2  -11  o   10369.4    976.0
1 2 1    0    2   12  o    1143.7    131.7
1 2 1    0   -2   12  o    1378.0    140.1
1 2 1    0   -2  -12  o    1090.9    134.2
1 2 1    0    2  -12  o    1290.6    135.6
1 2 1    0   -2   13  o   44492.5   4083.0
1 2 1    0    2  -13  o   44931.1   4102.4
1 2 1    0    2   14  o    6646.0    952.7
1 2 1    0   -2   14  o    6516.6    949.1
1 2 1    0   -2  -14  o    7589.5    952.0
1 2 1    0    2  -14  o    7559.7    954.2
1 2 1    0    2   16  o     256.6     68.6
1 2 1    0   -2   16  o     137.1     67.8
1 2 1    0   -2  -16  o     307.6     69.9
1 2 1    0    2  -16  o     265.4     71.5
"""
  open("r2etd-sf.cif", "w").write(cif_in)
  args = ["phenix.cif_as_mtz", "r2etd-sf.cif", "--merge"]
  rc = easy_run.fully_buffered(" ".join(args)).raise_if_errors().return_code
  assert (rc == 0)
  hkl_in = file_reader.any_file("r2etd-sf.mtz")
  l_w = []
  for array in hkl_in.file_server.miller_arrays :
    l_w.append((array.info().label_string(), "%.4f" % array.info().wavelength))
  assert (l_w == [
    ('FC,PHIFC', '1.0163'),
    ('FOBS,SIGFOBS', '1.0163'),
    ('R-free-flags', '1.0163'),
    ('I(+),SIGI(+),I(-),SIGI(-)', '1.0163'),
    ('I2(+),SIGI2(+),I2(-),SIGI2(-)', '0.9797')
  ])

if (__name__ == "__main__"):
  exercise()
  print("OK")
