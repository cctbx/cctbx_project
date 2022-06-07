
from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from mmtbx.command_line import cif_as_mtz
from iotbx import mtz
from iotbx import file_reader
from libtbx.test_utils import approx_equal

def exercise1():
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
  with open("r2etd-sf.cif", "w") as f:
    f.write(cif_in)
  args = ["phenix.cif_as_mtz", "r2etd-sf.cif", "--merge"]
  rc = easy_run.fully_buffered(" ".join(args)).raise_if_errors().return_code
  assert (rc == 0)
  hkl_in = file_reader.any_file("r2etd-sf.mtz")
  l_w = []
  for array in hkl_in.file_server.miller_arrays :
    l_w.append((array.info().label_string(), "%.4f" % array.info().wavelength))
  assert (sorted(l_w) == sorted([
    ('FOBS,SIGFOBS', '1.0163'),
    ('FC,PHIFC', '1.0163'),
    ('R-free-flags', '1.0163'),
    ('I(+),SIGI(+),I(-),SIGI(-)', '1.0163'),
    ('I2(+),SIGI2(+),I2(-),SIGI2(-)', '0.9797')
  ]) )


r1vjzsf_cif_snip = """
data_r1vjzsf
#
loop_
_audit.revision_id
_audit.creation_date
_audit.update_record
1   2004-04-13 'initial processing'
1_1 2009-01-27 'correct token name  _refln.pdbx_HLA (_HLB, _HLC, _HLD)'
1_2 2013-12-04 'Format standardization'
#
_cell.entry_id      1vjz
_cell.length_a      64.6840
_cell.length_b      64.6840
_cell.length_c      202.1890
_cell.angle_alpha   90.0000
_cell.angle_beta    90.0000
_cell.angle_gamma   90.0000
#
loop_
_diffrn_radiation_wavelength.id
_diffrn_radiation_wavelength.wavelength
1 0.98050
2 0.964063
3 0.979292
#
_entry.id   1vjz
#
_exptl_crystal.id   1
#
_reflns_scale.group_code   1
#
_symmetry.entry_id               1vjz
_symmetry.space_group_name_H-M   'P 43 21 2'
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
1 1 1    0    0    4 x        ?      ?          ?       ?
1 1 1    0    0    8 x        ?      ?     7109.5   180.0
1 1 1    0    0   12 x        ?      ?     8288.9     0.0
1 1 1    0    0   16 x        ?      ?       61.6   180.0
1 1 1    0    0   88 x        ?      ?      152.2   180.0
1 1 1    0    0   92 x        ?      ?      283.3   180.0
1 1 1    0    0   96 x        ?      ?        4.0     0.0
1 1 1    1    0    4 x        ?      ?          ?       ?
1 1 1    1    0    5 x        ?      ?          ?       ?
1 1 1    1    0    6 o      2.9    2.2     2024.1   180.0
1 1 1    1    0   20 x        ?      ?      156.2    90.0
1 1 1    1    0   21 x        ?      ?       55.4   225.0
1 1 1    1    0   22 x        ?      ?      372.1     0.0
1 1 1    1    0   23 x        ?      ?      787.1   315.0
1 1 1    1    0   24 x        ?      ?     1019.9    90.0
1 1 1   29    1   32 o     65.1   31.3       51.5   162.3
1 1 1   29    1   33 o    114.8   41.1      158.8   333.1
1 1 1   29    1   34 o     56.3   26.4       91.5   288.3
1 1 1   29    1   35 o     65.4   30.2       25.1   310.6
1 1 1   29    1   36 o    153.4   30.3      186.2   324.9
1 1 1   29    1   37 o     93.5   34.8       98.6    95.9
1 1 1   29    1   38 o    169.8   28.2      162.2   118.8
1 1 1   29    2    0 o     34.4   24.5       34.0     0.0
1 1 1   31    3    0 x        ?      ?       60.1   180.0
1 1 1   31    3    1 x        ?      ?       67.5   280.1
1 1 1   31    5    4 x        ?      ?       34.2   340.2
1 1 1   31    5    5 x        ?      ?       87.2   143.7
1 1 1   31    5    6 x        ?      ?       98.8   296.4
1 1 1   31    5    7 x        ?      ?       82.1   258.8
1 1 1   31    5    8 x        ?      ?       26.0    52.6
#END
data_r1vjzAsf
#
_cell.entry_id      1vjz
_cell.length_a      64.6840
_cell.length_b      64.6840
_cell.length_c      202.1890
_cell.angle_alpha   90.0000
_cell.angle_beta    90.0000
_cell.angle_gamma   90.0000
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details
'   unmerged original index intensities of refinement data set from crystal 9787'
#
_entry.id   1vjz
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
1 1 1   -1   -1   10  o     951.9     69.6
1 1 1   -1   -1   11  o    9456.0    609.5
1 1 1    1    1  -11  o    9620.4    605.4
1 1 1   -1   -1   12  o   16133.9   1012.7
1 1 1    1    1  -12  o   15721.5   1009.6
1 1 1    1    1  -13  o    1927.0    136.7
1 1 1    1    1  -14  o    3807.9    256.6
1 1 1    1    1  -15  o   15116.1    974.6
1 1 1    1    1  -16  o     234.3     31.4
1 1 1   29  -10   17  o    -252.6    243.5
1 1 1   29  -10   18  o     156.6    249.9
1 1 1   29  -10   19  o      15.8    240.5
1 1 1   29  -10   20  o     -44.9    252.3
1 1 1  -29   10  -21  o     237.6    205.8
1 1 1   29  -10   21  o     355.8    255.2
1 1 1  -29   10  -22  o      39.6    184.6
1 1 1   29  -10   22  o    -420.4    253.3
1 1 1  -29   10  -23  o     574.3    226.0
1 1 1   29  -10   23  o     147.9    293.5
1 1 1   30   -6   21  o    -767.8    264.6
1 1 1   30   -6   22  o     371.8    262.9
#END
data_r1vjzBsf
#
_cell.entry_id      1vjz
_cell.length_a      64.299
_cell.length_b      64.299
_cell.length_c      201.930
_cell.angle_alpha   90.0000
_cell.angle_beta    90.0000
_cell.angle_gamma   90.0000
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details
;   unmerged original index intensities of data set from crystal 9774 for phasing
;
#
_entry.id   1vjz
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
2 2 1    0    0   -8  o    1242.5     76.9
2 2 1    0    0   16  o   16497.6    696.8
2 2 1    0    0  -16  o   16886.5    627.7
2 2 1    0    0   20  o    1152.8     84.5
2 2 1    0    0  -20  o    1245.8     91.2
2 2 1    0    0   24  o   56057.3   2118.1
2 2 1    0    0   24  o   59770.5   2186.0
2 2 1   -6  -25    3  o      22.2     92.9
2 2 1   -6   25    4  o     277.8    206.5
2 2 1   -6  -25   -4  o      68.3     87.2
2 2 1   -6   25   -4  o     160.7    147.9
2 2 1   -6  -25    4  o      91.5     93.8
2 2 1   -6   25    5  o     -91.5    130.8
2 2 1   -6  -25   -5  o      97.6    104.4
2 2 1   -6   25   -5  o     130.0    106.8
2 2 1   -6  -25    5  o      -9.6    109.3
2 2 1    6   25    5  o      67.3     80.6
2 2 1   -6   25   -6  o     -41.7    118.2
2 2 1   -6  -25    6  o      30.7    107.1
#END
data_r1vjzCsf
#
_cell.entry_id      1vjz
_cell.length_a      64.635
_cell.length_b      64.635
_cell.length_c      202.270
_cell.angle_alpha   90.0000
_cell.angle_beta    90.0000
_cell.angle_gamma   90.0000
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details
;   unmerged original index intensities of data set from crystal 9774 for phasing
;
#
_entry.id   1vjz
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
2 3 1    0    0    8  o     933.7     60.3
2 3 1    0    0   -8  o     984.1     57.7
2 3 1    0    0   16  o   15152.9    586.0
2 3 1    0    0  -16  o   14897.3    566.8
2 3 1    0    0   20  o     772.2     57.1
2 3 1    0    0  -20  o     816.3     77.3
2 3 1    5  -25   -8  o      57.2    193.8
2 3 1    5  -25    9  o     -43.3    123.1
2 3 1    5   25   -9  o      59.0    148.2
2 3 1   -5   25   -9  o      46.7    207.2
2 3 1    5  -25   -9  o     181.2    215.8
2 3 1    6   25    4  o     215.7    297.9
2 3 1    6   25   -5  o      17.4    108.4
2 3 1   -6   25   -5  o      24.5    124.0
2 3 1    6   25   -6  o      46.3    101.5
2 3 1   -6   25   -6  o    2189.4    287.3
2 3 1    6   25   -7  o      71.2    113.9
2 3 1   -6   25   -7  o     -34.1    124.7
2 3 1   -6   25   -8  o     240.5    148.9
#END
data_r1vjzDsf
#
_cell.entry_id      1vjz
_cell.length_a      64.299
_cell.length_b      64.299
_cell.length_c      201.930
_cell.angle_alpha   90.0000
_cell.angle_beta    90.0000
_cell.angle_gamma   90.0000
#
_diffrn.id                  1
_diffrn.crystal_id          1
_diffrn.ambient_temp        ?
_diffrn.crystal_treatment   ?
_diffrn.details             '   phases from crystal 9774'
#
_entry.id   1vjz
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
2 2 1    0    0    4  o       ?       ?       ?       ?
2 2 1    0    0    8  o  -28.66    0.00    0.00    0.00
2 2 1    0    0   12  o       ?       ?       ?       ?
2 2 1    0    0   16  o  -39.51    0.00    0.00    0.00
2 2 1    0    0   20  o  -59.59    0.00    0.00    0.00
2 2 1    0    0   24  o    0.91    0.00    0.00    0.00
2 2 1    0    0   28  o  -19.74    0.00    0.00    0.00
2 2 1    0    0   32  o  -56.97    0.00    0.00    0.00
2 2 1   25    4   15  o    0.02    0.27    0.00    0.01
2 2 1   25    4   16  o       ?       ?       ?       ?
2 2 1   25    5    0  o    0.13    0.00    0.00    0.00
2 2 1   25    5    1  o   -0.20    0.12    0.00    0.00
2 2 1   25    5   13  o       ?       ?       ?       ?
2 2 1   25    6    0  o   -0.24    0.00    0.00    0.00
2 2 1   25    6    1  o   -0.42    0.45    0.01    0.00
2 2 1   25    6    2  o    0.06   -0.49    0.00    0.00
2 2 1   25    6    3  o    0.48    0.32    0.00   -0.01
2 2 1   25    6    4  o   -0.36   -0.21    0.00    0.00
2 2 1   25    6    5  o   -0.06   -0.05    0.00    0.00
2 2 1   25    6    6  o   -0.51   -0.49    0.00    0.00
2 2 1   25    6    7  o       ?       ?       ?       ?
2 2 1   25    6    8  o       ?       ?       ?       ?
#END OF REFLECTIONS
"""

def exercise2():
  """
  Check individual unit cell dimensions and wavelengths of cif datasets are retained when
  reading an mtz file originating from a cif file
  """
  with open("r1vjzsf_snip.cif", "w") as f:
    f.write(r1vjzsf_cif_snip)
  args = ["phenix.cif_as_mtz", "r1vjzsf_snip.cif", "--merge"]
  rc = easy_run.fully_buffered(" ".join(args)).raise_if_errors().return_code
  assert (rc == 0)
  hkl_in = file_reader.any_file("r1vjzsf_snip.mtz")
  l_w = []
  for array in hkl_in.file_server.miller_arrays:
    uc = array.unit_cell().parameters()
    l_w.append((array.info().label_string(),
                "%.5f"%array.info().wavelength, "%.4f"%uc[0], "%.4f"%uc[1], "%.4f"%uc[2], "%.3f"%uc[3], "%.3f"%uc[4], "%.3f"%uc[5]))
  assert (sorted(l_w) == sorted([
    ('R-free-flags',                  '0.98050', '64.6840', '64.6840', '202.1890', '90.000', '90.000', '90.000'),
    ('FOBS,SIGFOBS',                  '0.98050', '64.6840', '64.6840', '202.1890', '90.000', '90.000', '90.000'),
    ('FC,PHIFC',                      '0.98050', '64.6840', '64.6840', '202.1890', '90.000', '90.000', '90.000'),
    ('I(+),SIGI(+),I(-),SIGI(-)',     '0.98050', '64.6840', '64.6840', '202.1890', '90.000', '90.000', '90.000'),
    ('I2(+),SIGI2(+),I2(-),SIGI2(-)', '0.96406', '64.2990', '64.2990', '201.9300', '90.000', '90.000', '90.000'),
    ('I3(+),SIGI3(+),I3(-),SIGI3(-)', '0.97929', '64.6350', '64.6350', '202.2700', '90.000', '90.000', '90.000'),
    ('R-free-flags2',                 '0.96406', '64.2990', '64.2990', '201.9300', '90.000', '90.000', '90.000'),
    ('HLA,HLB,HLC,HLD',               '0.96406', '64.2990', '64.2990', '201.9300', '90.000', '90.000', '90.000')
  ]) )


def exercise3():
  """
  Check individual unit cell dimensions and wavelengths of cif datasets are retained in
  mtz crystals within the same mtz file
  """
  with open("r1vjzsf_snip.cif", "w") as f:
    f.write(r1vjzsf_cif_snip)
  args = ["phenix.cif_as_mtz", "r1vjzsf_snip.cif", "--merge"]
  rc = easy_run.fully_buffered(" ".join(args)).raise_if_errors().return_code
  assert (rc == 0)
  m = mtz.object(file_name="r1vjzsf_snip.mtz")
  assert (len(m.crystals()) == 5)
  eps = 1.e-5
  assert approx_equal(m.crystals()[1].datasets()[0].wavelength(), 0.98049998, eps)
  assert approx_equal(m.crystals()[2].datasets()[0].wavelength(), 0.96406001, eps)
  assert approx_equal(m.crystals()[3].datasets()[0].wavelength(), 0.97929000, eps)
  assert approx_equal(m.crystals()[4].datasets()[0].wavelength(), 0.96406000, eps)

  uc = m.crystals()[1].unit_cell_parameters()
  fuc = "%.4f"%uc[0], "%.4f"%uc[1], "%.4f"%uc[2], "%.3f"%uc[3], "%.3f"%uc[4], "%.3f"%uc[5]
  assert (('64.6840', '64.6840', '202.1890', '90.000', '90.000', '90.000') == fuc)

  uc = m.crystals()[2].unit_cell_parameters()
  fuc = "%.4f"%uc[0], "%.4f"%uc[1], "%.4f"%uc[2], "%.3f"%uc[3], "%.3f"%uc[4], "%.3f"%uc[5]
  assert (('64.2990', '64.2990', '201.9300', '90.000', '90.000', '90.000') == fuc)

  uc = m.crystals()[3].unit_cell_parameters()
  fuc = "%.4f"%uc[0], "%.4f"%uc[1], "%.4f"%uc[2], "%.3f"%uc[3], "%.3f"%uc[4], "%.3f"%uc[5]
  assert (('64.6350', '64.6350', '202.2700', '90.000', '90.000', '90.000') == fuc)

  uc = m.crystals()[4].unit_cell_parameters()
  fuc = "%.4f"%uc[0], "%.4f"%uc[1], "%.4f"%uc[2], "%.3f"%uc[3], "%.3f"%uc[4], "%.3f"%uc[5]
  assert (('64.2990', '64.2990', '201.9300', '90.000', '90.000', '90.000') == fuc)



if (__name__ == "__main__"):
  exercise1()
  exercise2()
  exercise3()
  print("OK")
