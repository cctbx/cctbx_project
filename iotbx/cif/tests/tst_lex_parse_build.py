from cctbx.array_family import flex
from cctbx import miller
from iotbx import cif
from iotbx.cif import CifParserError
from iotbx.cif.builders import CifBuilderError
from iotbx.reflection_file_reader import any_reflection_file
from libtbx.test_utils import \
     approx_equal, show_diff, Exception_expected, open_tmp_file
from cStringIO import StringIO
import sys


def exercise_miller_arrays_as_cif_block():
  from iotbx.cif import reader
  cif_model = reader(input_string=cif_miller_array,
                     builder=cif.builders.cif_model_builder()).model()
  ma_builder = cif.builders.miller_array_builder(cif_model['global'])
  ma1 = ma_builder.arrays()['_refln_F_squared_meas']
  mas_as_cif_block = cif.miller_arrays_as_cif_block(
    ma1, array_type='meas')
  mas_as_cif_block.add_miller_array(
    ma1.array(data=flex.complex_double([1-1j]*ma1.size())), array_type='calc')
  mas_as_cif_block.add_miller_array(
    ma1.array(data=flex.complex_double([1-2j]*ma1.size())), column_names=[
      '_refln_A_calc', '_refln_B_calc'])
  for key in ('_refln_F_squared_meas', '_refln_F_squared_sigma',
              '_refln_F_calc', '_refln_phase_calc',
              '_refln_A_calc', '_refln_A_calc'):
    assert key in mas_as_cif_block.cif_block.keys()
  #
  mas_as_cif_block = cif.miller_arrays_as_cif_block(
    ma1, column_names=['_diffrn_refln_intensity_net',
                       '_diffrn_refln_intensity_sigma'],
         miller_index_prefix='_diffrn_refln_')
  mas_as_cif_block.add_miller_array(
    ma1.array(data=flex.std_string(ma1.size(), 'om')),
    column_name='_diffrn_refln_intensity_u')
  for key in ('_diffrn_refln_intensity_net', '_diffrn_refln_intensity_sigma',
              '_diffrn_refln_intensity_u'):
    assert key in mas_as_cif_block.cif_block.keys()
  #
  try: reader(input_string=cif_global)
  except CifParserError, e: pass
  else: raise Exception_expected
  cif_model = reader(input_string=cif_global, strict=False).model()
  assert not show_diff(str(cif_model), """\
data_1
_c                                3
_d                                4
""")


def exercise_lex_parse_build():
  exercise_parser(cif.reader, cif.builders.cif_model_builder)
  cm = cif.reader(input_string=cif_quoted_string).model()
  assert cm['global']['_a'] == 'a"b'
  assert cm['global']['_b'] == "a dog's life"
  stdout = sys.stdout
  s = StringIO()
  sys.stdout = s
  try: cif.reader(input_string=cif_invalid_missing_value)
  except CifParserError: pass
  else: raise Exception_expected
  r = cif.reader(
    input_string=cif_invalid_missing_value, raise_if_errors=False)
  assert r.error_count() == 1
  try: cif.reader(input_string=cif_invalid_string)
  except CifParserError: pass
  else: raise Exception_expected
  a = cif.reader(input_string=cif_cod)
  assert a.error_count() == 0
  try: cif.reader(input_string=cif_invalid_semicolon_text_field)
  except CifParserError: pass
  else: raise Exception_expected
  d = cif.reader(input_string=cif_valid_semicolon_text_field)
  assert d.error_count() == 0
  assert d.model()['1']['_a'] == '\n1\n'
  e = cif.reader(input_string=cif_unquoted_string_semicolon)
  assert not show_diff(str(e.model()), """\
data_1
_a                                ;1
_b                                ;
_c                                2
""")
  cif_str_1 = """\
data_1
_a 1
"""
  cif_str_2 = """\
data_2
_b 2
"""
  cm = cif.reader(input_string=cif_str_1).model()
  assert cm.keys() == ['1']
  cif.reader(input_string=cif_str_2, cif_object=cm).model()
  assert cm.keys() == ['1', '2']

  sys.stdout = stdout

  arrays = miller.array.from_cif(file_object=StringIO(
    cif_miller_array_template %(
      '_refln_F_calc', '_refln_F_meas', '_refln_F_sigma')),
                                 data_block_name='global')
  assert sorted(arrays.keys()) == ['_refln_F_calc', '_refln_F_meas']
  assert arrays['_refln_F_calc'].sigmas() is None
  assert isinstance(arrays['_refln_F_meas'].sigmas(), flex.double)
  arrays = miller.array.from_cif(file_object=StringIO(
    cif_miller_array_template %(
      '_refln_A_calc', '_refln_B_calc', '_refln_F_meas')),
                                 data_block_name='global')
  assert sorted(arrays.keys()) == ['_refln_A_calc', '_refln_F_meas']
  assert arrays['_refln_A_calc'].is_complex_array()
  arrays = miller.array.from_cif(file_object=StringIO(
    cif_miller_array_template %(
      '_refln_A_meas', '_refln_B_meas', '_refln_F_meas')),
                                 data_block_name='global')
  assert sorted(arrays.keys()) == ['_refln_A_meas', '_refln_F_meas']
  assert arrays['_refln_A_meas'].is_complex_array()
  arrays = miller.array.from_cif(file_object=StringIO(
    cif_miller_array_template %(
      '_refln_intensity_calc', '_refln_intensity_meas',
      '_refln_intensity_sigma')),
                                 data_block_name='global')
  assert sorted(arrays.keys()) == [
    '_refln_intensity_calc', '_refln_intensity_meas']
  arrays = miller.array.from_cif(file_object=StringIO(
    cif_miller_array_template %(
      '_refln_F_calc', '_refln_phase_calc', '_refln_F_sigma')),
                                 data_block_name='global')
  assert sorted(arrays.keys()) == ['_refln_F_calc']
  assert arrays['_refln_F_calc'].is_complex_array()

  for data_block_name in (None, "global"):
    miller_arrays = cif.reader(file_object=StringIO(
      cif_miller_array_template %(
        '_refln_F_calc',
        '_refln_F_meas',
        '_refln_F_sigma'))).as_miller_arrays(data_block_name=data_block_name)
    assert " ".join(sorted([str(ma.info()) for ma in miller_arrays])) \
      == "cif:global,_refln_F_calc cif:global,_refln_F_meas,_refln_F_sigma"
  f = open_tmp_file(suffix="cif")
  f.write(cif_miller_array_template %(
        '_refln_F_calc',
        '_refln_F_meas',
        '_refln_F_sigma'))
  f.close()
  miller_arrays = any_reflection_file(file_name=f.name).as_miller_arrays()
  assert len(miller_arrays) == 2

def exercise_parser(reader, builder):
  cif_model = reader(
    input_string=cif_xray_structure, builder=builder()).model()
  xs_builder = cif.builders.crystal_structure_builder(cif_model['global'])
  xs1 = xs_builder.structure
  # also test construction of cif model from xray structure
  xs_cif_block = xs1.as_cif_block()
  xs2 = cif.builders.crystal_structure_builder(xs_cif_block).structure
  sio = StringIO()
  xs1.as_cif_simple(out=sio)
  xs3 = reader(input_string=sio.getvalue()).build_crystal_structures(
    data_block_name='global')
  xs_iso = xs1.deep_copy_scatterers()
  xs_iso.convert_to_isotropic()
  xs4 = cif.builders.crystal_structure_builder(xs_iso.as_cif_block()).structure

  for xs in (xs1, xs2, xs3):
    sc = xs.scatterers()
    assert list(sc.extract_labels()) == ['o','c']
    assert list(sc.extract_scattering_types()) == ['O','C']
    assert approx_equal(sc.extract_occupancies(), (0.8, 1))
    assert approx_equal(sc.extract_sites(), ((0.5,0,0),(0,0,0)))
    assert approx_equal(sc.extract_u_star(),
      [(-1, -1, -1, -1, -1, -1), (1e-3, 5e-4, (1e-3)/3, 0, 0, 0)])
    assert approx_equal(sc.extract_u_iso(), (0.1, -1))
    assert approx_equal(xs.unit_cell().parameters(),
                        (10,20,30,90,90,90))
    assert str(xs.space_group_info()) == 'C 1 2/m 1'
  #
  cif_model = reader(
    input_string=cif_miller_array, builder=builder()).model()
  ma_builder = cif.builders.miller_array_builder(cif_model['global'])
  ma1 = ma_builder.arrays()['_refln_F_squared_meas']
  if isinstance(cif_model, cif.model.cif):
    assert (ma_builder.arrays()['_refln_observed_status'].data() ==
            flex.std_string(['o'] * ma1.size()))
  # also test construction of cif model from miller array
  ma2 = cif.builders.miller_array_builder(
    ma1.as_cif_block(array_type='meas')).arrays()['_refln_F_squared_meas']
  for ma in (ma1, ma2):
    sio = StringIO()
    ma.show_array(sio)
    assert not show_diff(sio.getvalue(), """\
(1, 0, 0) 748.71 13.87
(2, 0, 0) 1318.51 24.29
(3, 0, 0) 1333.51 33.75
(4, 0, 0) 196.58 10.85
(5, 0, 0) 3019.71 55.29
(6, 0, 0) 1134.38 23.94
(7, 0, 0) 124.01 15.16
(8, 0, 0) -1.22 10.49
(9, 0, 0) 189.09 20.3
(10, 0, 0) 564.68 35.61
(-10, 1, 0) 170.23 22.26
""")
    sio = StringIO()
    ma.show_summary(sio)
    assert not show_diff(sio.getvalue(), """\
Miller array info: cif:_refln_F_squared_meas,_refln_F_squared_sigma
Observation type: xray.intensity
Type of data: double, size=11
Type of sigmas: double, size=11
Number of Miller indices: 11
Anomalous flag: False
Unit cell: (7.9999, 9.3718, 14.7362, 82.625, 81.527, 81.726)
Space group: P -1 (No. 2)
""")


cif_xray_structure = """\
data_global
loop_
    _symmetry_equiv_pos_as_xyz
    'x, y, z'
    '-x, y, -z'
    '-x, -y, -z'
    'x, -y, z'
    'x+1/2, y+1/2, z'
    '-x+1/2, y+1/2, -z'
    '-x+1/2, -y+1/2, -z'
    'x+1/2, -y+1/2, z'
_cell_length_a 10
_cell_length_b 20
_cell_length_c 30
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_ # comment
    _atom_site_label # another comment
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_U_iso_or_equiv
    _atom_site_occupancy
    _atom_site_type_symbol
    'o' 0.5 0 0 0.1 0.8 'O'
    'c' 0 0 0 0.2 1 'C'
loop_
    _atom_site_aniso_label
    _atom_site_aniso_U_11
    _atom_site_aniso_U_22
    _atom_site_aniso_U_33
    _atom_site_aniso_U_12
    _atom_site_aniso_U_13
    _atom_site_aniso_U_23
    'c' 0.1 0.2 0.3 0 0 0
"""

cif_miller_array = """\
data_global
loop_
 _symmetry_equiv_pos_as_xyz
 'x, y, z'
 '-x, -y, -z'

_cell_length_a     7.9999
_cell_length_b     9.3718
_cell_length_c    14.7362
_cell_angle_alpha  82.625
_cell_angle_beta   81.527
_cell_angle_gamma  81.726

loop_
 _refln_index_h
 _refln_index_k
 _refln_index_l
 _refln_F_squared_calc
 _refln_F_squared_meas
 _refln_F_squared_sigma
 _refln_observed_status
   1   0   0      756.07      748.71     13.87 o
   2   0   0     1266.94     1318.51     24.29 o
   3   0   0     1381.53     1333.51     33.75 o
   4   0   0      194.77      196.58     10.85 o
   5   0   0     3102.74     3019.71     55.29 o
   6   0   0     1145.05     1134.38     23.94 o
   7   0   0      103.98      124.01     15.16 o
   8   0   0       16.94       -1.22     10.49 o
   9   0   0      194.74      189.09     20.30 o
  10   0   0      581.30      564.68     35.61 o
 -10   1   0      148.83      170.23     22.26 o
"""


cif_miller_array_template = """\
data_global

loop_
 _symmetry_equiv_pos_as_xyz
 'x, y, z'
 '-x, -y, -z'

_cell_length_a     7.9999
_cell_length_b     9.3718
_cell_length_c    14.7362
_cell_angle_alpha  82.625
_cell_angle_beta   81.527
_cell_angle_gamma  81.726

loop_
 _refln_index_h
 _refln_index_k
 _refln_index_l
 %s
 %s
 %s
   1 0 0 1.2 1.3 0.1
   2 0 0 2.3 2.4 0.2
   3 0 0 3.4 3.5 0.3
   4 0 0 4.5 6.7 0.4

# comment with WS before and after

"""

cif_cod = """\
data_global
_a 1
_[b] 2
_c e43
_d # comment
   # another comment
'1 2'
_e )
"""

cif_quoted_string = """\
data_global
# both of these are legal apparently
_a "a"b"
_b 'a dog's life'
"""

cif_invalid_loop = """\
data_global
loop_
  _a
  _b
_c 1
"""

cif_invalid_missing_value = """\
data_global
_a 1
_b
_c 3
"""

cif_invalid_string = """\
data_global
_a 'no closing quote
_b 1
"""

cif_invalid_semicolon_text_field = """\
data_1
_a ;
1
;
"""

cif_valid_semicolon_text_field = """\
data_1
_a
;
1
;
"""

cif_unquoted_string_semicolon = """\
data_1
_a ;1
_b ;
_c 2
"""

cif_global = """\
global_
_a 1
_b 2

data_1
_c 3
_d 4
"""

def exercise_partial_crystal_symmetry():
  def get_inp(u, s):
    result = ["data_test"]
    if (u):
      result.append("_cell_length_a 12.605(3)")
    if (s):
      result.append("_symmetry_space_group_name_Hall  '-P 2yn'")
    return "\n".join(result)
  def get_cs(input_string):
    cif_model = cif.reader(input_string=input_string).model()
    return cif.builders.crystal_symmetry_builder(
      cif_block=cif_model["test"]).crystal_symmetry
  cs = get_cs(get_inp(False, False))
  assert cs.unit_cell() is None
  assert cs.space_group_info() is None
  cs = get_cs(get_inp(False, True))
  assert cs.unit_cell() is None
  assert str(cs.space_group_info()) == "P 1 21/n 1"
  try:
    get_cs(get_inp(True, False))
  except CifBuilderError, e:
    assert str(e) == "Not all unit cell parameters are given in the cif file"
  else: raise Exception_expected

def exercise_mmcif_structure_factors():
  miller_arrays = cif.reader(input_string=r3adsrf).as_miller_arrays()
  assert len(miller_arrays) == 11
  hl_coeffs = find_miller_array_from_labels(
    miller_arrays, ','.join(['_refln.pdbx_HL_A_iso',
                             '_refln.pdbx_HL_B_iso',
                             '_refln.pdbx_HL_C_iso',
                             '_refln.pdbx_HL_D_iso']))
  assert hl_coeffs.is_hendrickson_lattman_array()
  assert hl_coeffs.size() == 10
  f_meas_au = find_miller_array_from_labels(
    miller_arrays, '_refln.F_meas_au,_refln.F_meas_sigma_au')
  assert f_meas_au.is_xray_amplitude_array()
  assert f_meas_au.size() == 10
  assert f_meas_au.sigmas() is not None
  assert f_meas_au.space_group_info().symbol_and_number() == 'C 1 2 1 (No. 5)'
  assert approx_equal(f_meas_au.unit_cell().parameters(),
                      (163.97, 45.23, 110.89, 90.0, 131.64, 90.0))
  pdbx_I_plus = find_miller_array_from_labels(
    miller_arrays, '_refln.pdbx_I_plus,_refln.pdbx_I_plus_sigma')
  assert pdbx_I_plus.is_xray_intensity_array()
  assert pdbx_I_plus.size() == 12
  pdbx_I_minus = find_miller_array_from_labels(
    miller_arrays, '_refln.pdbx_I_minus,_refln.pdbx_I_minus_sigma')
  assert pdbx_I_minus.is_xray_intensity_array()
  assert pdbx_I_minus.size() == 9
  assert pdbx_I_minus.unit_cell() is None     # no symmetry information in
  assert pdbx_I_minus.space_group() is None   # this CIF block
  #
  miller_arrays = cif.reader(input_string=r3ad7sf).as_miller_arrays()
  assert len(miller_arrays) == 7

def find_miller_array_from_labels(miller_arrays, labels):
  for ma in miller_arrays:
    if labels in str(ma.info()):
      return ma
  raise RuntimeError("Could not find miller array with labels %s" %labels)

r3adsrf = """
data_r3adrsf
_cell.length_a  163.970
_cell.length_b  45.230
_cell.length_c  110.890
_cell.angle_alpha  90.000
_cell.angle_beta  131.640
_cell.angle_gamma  90.000

#
_symmetry.entry_id    3adr
_symmetry.space_group_name_H-M 'C 1 2 1'

loop_
_refln.crystal_id
_refln.wavelength_id
_refln.scale_group_code
_refln.status
_refln.index_h
_refln.index_k
_refln.index_l
_refln.F_meas_au
_refln.F_meas_sigma_au
_refln.pdbx_HL_A_iso
_refln.pdbx_HL_B_iso
_refln.pdbx_HL_C_iso
_refln.pdbx_HL_D_iso
1 1 1 o    0    2    2   484.70  11.74  -3.168  -3.662   3.247   3.261
1 1 1 o    0    2    3   337.40   4.60   2.287  -3.629   0.053  -0.597
1 1 1 o    0    2    4   735.40   8.43   0.476  -8.605  -4.861   0.921
1 1 1 o    0    2    5   433.70   4.76   2.453  -1.663   0.000  -9.540
1 1 1 o    0    2    6   435.20   4.96  -9.081   3.661   1.511  -2.736
1 1 1 o    0    2    7   427.50   4.56   8.677  -0.013   2.307  -4.220
1 1 1 o    0    2    8   661.00   7.26  -6.468   5.112  -5.411  -6.731
1 1 1 o    0    2    9   211.90   2.64  -2.279   0.914   0.455  -6.038
1 1 1 f    0    2   10   466.80   5.01  -0.005   6.330  -1.916  -0.598
1 1 1 o    0    2   11   424.20   4.66  -0.938   7.011  -1.603   1.724
#END
data_r3adrAsf
#
#
#
loop_
_refln.crystal_id
_refln.wavelength_id
_refln.scale_group_code
_refln.index_h
_refln.index_k
_refln.index_l
_refln.pdbx_I_plus
_refln.pdbx_I_plus_sigma
_refln.pdbx_I_minus
_refln.pdbx_I_minus_sigma
1 3 1   87    5  -46       40.2     40.4        6.7     63.9
1 3 1   87    5  -45       47.8     29.7       35.1     30.5
1 3 1   87    5  -44       18.1     33.2        0.5     34.6
1 3 1   87    5  -43        6.1     45.4       12.9     51.6
1 3 1   87    5  -42       -6.6     45.6      -15.5     55.8
1 3 1   87    7  -37        6.3     43.4          ?        ?
1 3 1   87    7  -36      -67.2     55.4          ?        ?
1 3 1   88    2  -44          0       -1       35.0     38.5
1 3 1   88    2  -43          0       -1       57.4     41.5
1 3 1   88    4  -45       -1.0     46.1       -9.1     45.6
1 3 1   88    4  -44      -19.8     49.2        0.3     34.7
1 3 1   88    6  -44       -1.8     34.8          ?        ?
#END OF REFLECTIONS
"""

r3ad7sf = """
data_r3ad7sf
_cell.length_a             198.9488
_cell.length_b             198.9488
_cell.length_c             196.7646
_cell.angle_alpha           90.0000
_cell.angle_beta            90.0000
_cell.angle_gamma          120.0000

loop_
_symmetry_equiv.id
_symmetry_equiv.pos_as_xyz
1 'X,  Y,  Z'
2 'X-Y,  X,  Z+5/6'
3 '-Y,  X-Y,  Z+2/3'
4 '-X,  -Y,  Z+1/2'
5 '-X+Y,  -X,  Z+1/3'
6 'Y,  -X+Y,  Z+1/6'
7 '-Y,  -X,  -Z+1/6'
8 'X-Y,  -Y,  -Z'
9 'X,  X-Y,  -Z+5/6'
10 'Y,  X,  -Z+2/3'
11 '-X+Y,  Y,  -Z+1/2'
12 '-X,  -X+Y,  -Z+1/3'

loop_
_refln.wavelength_id
_refln.crystal_id
_refln.scale_group_code
_refln.index_h
_refln.index_k
_refln.index_l
_refln.status
_refln.F_meas_au
_refln.F_meas_sigma_au
_refln.F_calc
_refln.phase_calc
_refln.fom
1 1 1       0    0    6 o     267.2   12.1    353.2   180.0  0.04
1 1 1       0    0   12 o    4700.0  113.0   3962.4   360.0  1.00
1 1 1       0    0   18 o   10214.0  222.9   8775.9   360.0  1.00
1 1 1       0    0   24 o    8268.9  192.9  11557.1   180.0  1.00
1 1 1       0    0   30 o    3274.6   77.5   2214.5     0.0  1.00
1 1 1       0    0   36 o     317.8   30.4   1993.3     0.0  0.05
1 1 1       0    0   42 o   12026.4  286.2   6514.5   180.0  1.00
1 1 1       0    0   48 o    1972.6   51.4   1357.9   180.0  0.91
"""

def exercise():
  if not cif.has_antlr3:
    print "Skipping tst_lex_parse_build.py (antlr3 is not available)"
    return
  exercise_miller_arrays_as_cif_block()
  exercise_lex_parse_build()
  exercise_partial_crystal_symmetry()
  exercise_mmcif_structure_factors()

if __name__ == '__main__':
  exercise()
  print "OK"
