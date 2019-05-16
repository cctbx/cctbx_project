from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import crystal, miller, sgtbx, uctbx
from iotbx import cif
from iotbx.cif import CifParserError
from iotbx.cif.builders import CifBuilderError
from iotbx import crystal_symmetry_from_any
from iotbx.reflection_file_reader import any_reflection_file
from libtbx.test_utils import \
     approx_equal, show_diff, Exception_expected, open_tmp_file
from cStringIO import StringIO
import sys
from six.moves import zip


def exercise_miller_arrays_as_cif_block():
  from iotbx.cif import reader
  cif_model = reader(input_string=cif_miller_array,
                     builder=cif.builders.cif_model_builder()).model()
  ma_builder = cif.builders.miller_array_builder(cif_model['global'])
  ma1 = ma_builder.arrays()['_refln_F_squared_meas']
  mas_as_cif_block = cif.miller_arrays_as_cif_block(
    ma1, array_type='meas', format="corecif")
  mas_as_cif_block.add_miller_array(
    ma1.array(data=flex.complex_double([1-1j]*ma1.size())), array_type='calc')
  mas_as_cif_block.add_miller_array(
    ma1.array(data=flex.complex_double([1-2j]*ma1.size())), column_names=[
      '_refln_A_calc', '_refln_B_calc'])
  for key in ('_refln_F_squared_meas', '_refln_F_squared_sigma',
              '_refln_F_calc', '_refln_phase_calc',
              '_refln_A_calc', '_refln_A_calc'):
    assert (key in mas_as_cif_block.cif_block.keys()), key
  #
  mas_as_cif_block = cif.miller_arrays_as_cif_block(
    ma1, array_type='meas', format="mmcif")
  mas_as_cif_block.add_miller_array(
    ma1.array(data=flex.complex_double([1-1j]*ma1.size())), array_type='calc')
  for key in ('_refln.F_squared_meas', '_refln.F_squared_sigma',
              '_refln.F_calc', '_refln.phase_calc',
              '_space_group_symop.operation_xyz',
              '_cell.length_a', '_refln.index_h'):
    assert key in mas_as_cif_block.cif_block.keys()
  #
  mas_as_cif_block = cif.miller_arrays_as_cif_block(
    ma1, column_names=['_diffrn_refln_intensity_net',
                       '_diffrn_refln_intensity_sigma'],
         miller_index_prefix='_diffrn_refln')
  mas_as_cif_block.add_miller_array(
    ma1.array(data=flex.std_string(ma1.size(), 'om')),
    column_name='_diffrn_refln_intensity_u')
  for key in ('_diffrn_refln_intensity_net', '_diffrn_refln_intensity_sigma',
              '_diffrn_refln_intensity_u'):
    assert key in list(mas_as_cif_block.cif_block.keys())
  #
  try: reader(input_string=cif_global)
  except CifParserError as e: pass
  else: raise Exception_expected
  cif_model = reader(input_string=cif_global, strict=False).model()
  assert not show_diff(str(cif_model), """\
data_1
_c                                3
_d                                4
""")
  # exercise adding miller arrays with non-matching indices
  cs = crystal.symmetry(unit_cell=uctbx.unit_cell((10, 10, 10, 90, 90, 90)),
                        space_group_info=sgtbx.space_group_info(symbol="P1"))
  mi = flex.miller_index(((1,0,0), (1,2,3), (2,3,4)))
  ms1 = miller.set(cs, mi)
  ma1 = miller.array(ms1, data=flex.double((1,2,3)))
  mas_as_cif_block = cif.miller_arrays_as_cif_block(
    ma1, column_name="_refln.F_meas_au")
  ms2 = miller.set(cs, mi[:2])
  ma2 = miller.array(ms2, data=flex.complex_double([1-2j]*ms2.size()))
  mas_as_cif_block.add_miller_array(
    ma2, column_names=("_refln.F_calc_au", "_refln.phase_calc")),
  ms3 = miller.set(cs, flex.miller_index(((1,0,0), (5,6,7), (2,3,4))))
  ma3 = miller.array(ms3, data=flex.double((4,5,6)))
  mas_as_cif_block.add_miller_array(ma3, column_name="_refln.F_squared_meas")
  ms4 = miller.set(cs, flex.miller_index(((1,2,3), (5,6,7), (1,1,1), (1,0,0), (2,3,4))))
  ma4 = ms4.d_spacings()
  mas_as_cif_block.add_miller_array(ma4, column_name="_refln.d_spacing")
  # extract arrays from cif block and make sure we get back what we started with
  arrays = cif.builders.miller_array_builder(mas_as_cif_block.cif_block).arrays()
  recycled_arrays = (arrays['_refln.F_meas_au'],
                     arrays['_refln.F_calc_au'],
                     arrays['_refln.F_squared_meas'],
                     arrays['_refln.d_spacing'])
  for orig, recycled in zip((ma1, ma2, ma3, ma4), recycled_arrays):
    assert orig.size() == recycled.size()
    recycled = recycled.customized_copy(anomalous_flag=orig.anomalous_flag())
    orig, recycled = orig.common_sets(recycled)
    assert orig.indices().all_eq(recycled.indices())
    assert approx_equal(orig.data(), recycled.data(), eps=1e-5)
  #
  cif_model = reader(input_string=r3adrsf,
                     builder=cif.builders.cif_model_builder()).model()
  cs = cif.builders.crystal_symmetry_builder(cif_model["r3adrsf"]).crystal_symmetry

  ma_builder = cif.builders.miller_array_builder(
    cif_model['r3adrAsf'],
    base_array_info=miller.array_info(crystal_symmetry_from_file=cs))
  miller_arrays = list(ma_builder.arrays().values())
  assert len(miller_arrays) == 4
  mas_as_cif_block = cif.miller_arrays_as_cif_block(
      miller_arrays[0].map_to_asu(),
      column_names=miller_arrays[0].info().labels,
      format="corecif")
  for array in miller_arrays[1:]:
    labels = array.info().labels
    if len(labels) > 1 :
      for label in labels :
        if label.startswith("wavelength_id"):
          labels.remove(label)
    mas_as_cif_block.add_miller_array(
      array=array.map_to_asu(), column_names=array.info().labels)
  s = StringIO()
  print(mas_as_cif_block.refln_loop, file=s)
  assert not show_diff(s.getvalue(), """\
loop_
  _refln_index_h
  _refln_index_k
  _refln_index_l
  _refln.crystal_id
  _refln.scale_group_code
  _refln.wavelength_id
  _refln.pdbx_I_plus
  _refln.pdbx_I_plus_sigma
  _refln.pdbx_I_minus
  _refln.pdbx_I_minus_sigma
  -87  5  46  1  1  3   40.2  40.4    6.7  63.9
  -87  5  45  1  1  3   47.8  29.7   35.1  30.5
  -87  5  44  1  1  3   18.1  33.2    0.5  34.6
  -87  5  43  1  1  3    6.1  45.4   12.9  51.6
  -87  5  42  1  1  3   -6.6  45.6  -15.5  55.8
  -87  7  37  1  1  3    6.3  43.4      ?     ?
  -87  7  36  1  1  3  -67.2  55.4      ?     ?
  -88  2  44  1  1  3      0    -1     35  38.5
  -88  2  43  1  1  3      0    -1   57.4  41.5
  -88  4  45  1  1  3     -1  46.1   -9.1  45.6
  -88  4  44  1  1  3  -19.8  49.2    0.3  34.7
  -88  6  44  1  1  3   -1.8  34.8      ?     ?

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
  assert list(cm.keys()) == ['1']
  cif.reader(input_string=cif_str_2, cif_object=cm).model()
  assert list(cm.keys()) == ['1', '2']
  try: cm = cif.reader(input_string=cif_invalid_loop).model()
  except CifParserError: pass
  else: raise Exception_expected
  try: cm = cif.reader(input_string=cif_invalid_loop_2).model()
  except CifParserError: pass
  else: raise Exception_expected

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
    ma1.as_cif_block(array_type='meas')).arrays()['_refln.F_squared_meas']
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
  ma1.show_summary(sio)
  # Miller array info: cif:_refln.F_squared_meas,_refln.F_squared_sigma
  assert not show_diff(sio.getvalue(), """\
Miller array info: cif:_refln_F_squared_meas,_refln_F_squared_sigma
Observation type: xray.intensity
Type of data: double, size=11
Type of sigmas: double, size=11
Number of Miller indices: 11
Anomalous flag: False
Unit cell: (7.999, 9.372, 14.736, 82.625, 81.527, 81.726)
Space group: P -1 (No. 2)
""")
  sio = StringIO()
  ma2.show_summary(sio)
  # Miller array info: cif:_refln.F_squared_meas,_refln.F_squared_sigma
  assert not show_diff(sio.getvalue(), """\
Miller array info: cif:_refln.F_squared_meas,_refln.F_squared_sigma
Observation type: xray.intensity
Type of data: double, size=11
Type of sigmas: double, size=11
Number of Miller indices: 11
Anomalous flag: False
Unit cell: (7.999, 9.372, 14.736, 82.625, 81.527, 81.726)
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

_cell_length_a      7.999
_cell_length_b      9.372
_cell_length_c     14.736
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

cif_invalid_loop_2 = """\
data_1
loop_
_a _b _c
1 2 3 4 5
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

def exercise_atom_type_loop():
  from cctbx import xray
  cif_model = cif.reader(input_string=cif_xray_structure).model()
  xs = cif.builders.crystal_structure_builder(cif_model['global']).structure
  xs.set_inelastic_form_factors(photon=0.71073, table="henke")
  loop = cif.atom_type_cif_loop(xray_structure=xs, format="mmcif")
  s = StringIO()
  print(loop, file=s)
  assert not show_diff(
    "\n".join([li.rstrip() for li in s.getvalue().splitlines()]), """\
loop_
  _atom_type.symbol
  _atom_type.scat_dispersion_real
  _atom_type.scat_dispersion_imag
  _atom_type.scat_Cromer_Mann_a1
  _atom_type.scat_Cromer_Mann_a2
  _atom_type.scat_Cromer_Mann_a3
  _atom_type.scat_Cromer_Mann_a4
  _atom_type.scat_Cromer_Mann_a5
  _atom_type.scat_Cromer_Mann_a6
  _atom_type.scat_Cromer_Mann_b1
  _atom_type.scat_Cromer_Mann_b2
  _atom_type.scat_Cromer_Mann_b3
  _atom_type.scat_Cromer_Mann_b4
  _atom_type.scat_Cromer_Mann_b5
  _atom_type.scat_Cromer_Mann_b6
  _atom_type.scat_Cromer_Mann_c
  _atom_type.scat_source
  _atom_type.scat_dispersion_source
  C  0.00347  0.00161  2.18189  1.77612  1.08772  0.64146  0.20789  0.10522  13.45337  32.57901  0.74729  0.25125  80.97993  0.05873  0.0
;
6-Gaussian fit: Grosse-Kunstleve RW, Sauter NK, Adams PD:
Newsletter of the IUCr Commission on Crystallographic Computing 2004, 3, 22-31.
;
  'Henke, Gullikson and Davis, At. Data and Nucl. Data Tables, 1993, 54, 2'
  O  0.01158  0.00611  2.91262  2.58808  0.98057  0.69663  0.68508  0.13677  14.48462   6.03818  0.42255  0.15446  35.53892  0.03841  0.0
;
6-Gaussian fit: Grosse-Kunstleve RW, Sauter NK, Adams PD:
Newsletter of the IUCr Commission on Crystallographic Computing 2004, 3, 22-31.
;
  'Henke, Gullikson and Davis, At. Data and Nucl. Data Tables, 1993, 54, 2'
""")

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
  except CifBuilderError as e:
    assert str(e) == "Not all unit cell parameters are given in the cif file"
  else: raise Exception_expected

def exercise_crystal_symmetry():
  cm = cif.reader(input_string=p1_sym_ops).model()
  cs_builder = cif.builders.crystal_symmetry_builder(cm["r1e5xsf"])
  assert cs_builder.crystal_symmetry.space_group_info().symbol_and_number() \
         == 'P 1 (No. 1)'
  file_object = open_tmp_file(suffix=".cif")
  file_object.write(p1_sym_ops)
  file_object.close()
  cs = crystal_symmetry_from_any.extract_from(file_name=file_object.name)
  assert cs.space_group_info().symbol_and_number() == 'P 1 (No. 1)'


def exercise_mmcif_structure_factors():
  miller_arrays = cif.reader(input_string=r3adrsf).as_miller_arrays()
  assert len(miller_arrays) == 16
  hl_coeffs = find_miller_array_from_labels(
    miller_arrays, ','.join([
      'scale_group_code=1', 'crystal_id=2', 'wavelength_id=3',
      '_refln.pdbx_HL_A_iso', '_refln.pdbx_HL_B_iso',
      '_refln.pdbx_HL_C_iso', '_refln.pdbx_HL_D_iso']))
  assert hl_coeffs.is_hendrickson_lattman_array()
  assert hl_coeffs.size() == 2
  mas_as_cif_block = cif.miller_arrays_as_cif_block(
    hl_coeffs, column_names=('_refln.pdbx_HL_A_iso', '_refln.pdbx_HL_B_iso',
                             '_refln.pdbx_HL_C_iso', '_refln.pdbx_HL_D_iso'))
  abcd = []
  for key in ('_refln.pdbx_HL_A_iso', '_refln.pdbx_HL_B_iso',
              '_refln.pdbx_HL_C_iso', '_refln.pdbx_HL_D_iso'):
    assert key in list(mas_as_cif_block.cif_block.keys())
    abcd.append(flex.double(mas_as_cif_block.cif_block[key]))
  hl_coeffs_from_cif_block = flex.hendrickson_lattman(*abcd)
  assert approx_equal(hl_coeffs.data(), hl_coeffs_from_cif_block)
  f_meas_au = find_miller_array_from_labels(
    miller_arrays, ','.join([
      'scale_group_code=1', 'crystal_id=1', 'wavelength_id=1',
      '_refln.F_meas_au', '_refln.F_meas_sigma_au']))
  assert f_meas_au.is_xray_amplitude_array()
  assert f_meas_au.size() == 5
  assert f_meas_au.sigmas() is not None
  assert f_meas_au.space_group_info().symbol_and_number() == 'C 1 2 1 (No. 5)'
  assert approx_equal(f_meas_au.unit_cell().parameters(),
                      (163.97, 45.23, 110.89, 90.0, 131.64, 90.0))
  pdbx_I_plus_minus = find_miller_array_from_labels(
    miller_arrays, ','.join(
      ['_refln.pdbx_I_plus', '_refln.pdbx_I_plus_sigma',
       '_refln.pdbx_I_minus', '_refln.pdbx_I_minus_sigma']))
  assert pdbx_I_plus_minus.is_xray_intensity_array()
  assert pdbx_I_plus_minus.anomalous_flag()
  assert pdbx_I_plus_minus.size() == 21
  assert pdbx_I_plus_minus.unit_cell() is None     # no symmetry information in
  assert pdbx_I_plus_minus.space_group() is None   # this CIF block
  #
  miller_arrays = cif.reader(input_string=r3ad7sf).as_miller_arrays()
  assert len(miller_arrays) == 11
  f_calc = find_miller_array_from_labels(
    miller_arrays, ','.join([
      'crystal_id=2', 'wavelength_id=1', '_refln.F_calc', '_refln.phase_calc']))
  assert f_calc.is_complex_array()
  assert f_calc.size() == 4
  #
  miller_arrays = cif.reader(input_string=integer_observations).as_miller_arrays()
  assert len(miller_arrays) == 2
  assert isinstance(miller_arrays[0].data(), flex.double)
  assert isinstance(miller_arrays[0].sigmas(), flex.double)
  #
  miller_arrays = cif.reader(input_string=r3v56sf).as_miller_arrays()
  assert len(miller_arrays) == 2
  for ma in miller_arrays: assert ma.is_complex_array()
  assert miller_arrays[0].info().labels == [
    'r3v56sf', '_refln.pdbx_DELFWT', '_refln.pdbx_DELPHWT']
  assert miller_arrays[1].info().labels == [
    'r3v56sf', '_refln.pdbx_FWT', '_refln.pdbx_PHWT']


def find_miller_array_from_labels(miller_arrays, labels):
  for ma in miller_arrays:
    if labels in str(ma.info()):
      return ma
  raise RuntimeError("Could not find miller array with labels %s" %labels)

r3adrsf = """
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
2 3 1 o    0    2    7   427.50   4.56   8.677  -0.013   2.307  -4.220
2 3 1 o    0    2    8   661.00   7.26  -6.468   5.112  -5.411  -6.731
2 3 4 o    0    2    9   211.90   2.64  -2.279   0.914   0.455  -6.038
2 3 4 f    0    2   10   466.80   5.01  -0.005   6.330  -1.916  -0.598
2 3 4 o    0    2   11   424.20   4.66  -0.938   7.011  -1.603   1.724
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
1 2 1       0    0   30 o    3274.6   77.5   2214.5     0.0  1.00
1 2 1       0    0   36 o     317.8   30.4   1993.3     0.0  0.05
1 2 1       0    0   42 o   12026.4  286.2   6514.5   180.0  1.00
1 2 1       0    0   48 o    1972.6   51.4   1357.9   180.0  0.91
"""

r3v56sf = """
data_r3v56sf

_cell.length_a  121.6330
_cell.length_b  121.6330
_cell.length_c  157.2140
_cell.angle_alpha  90.0000
_cell.angle_beta  90.0000
_cell.angle_gamma  120.0000

loop_
_symmetry_equiv.id
_symmetry_equiv.pos_as_xyz
1  'X,  Y,  Z'
2  'X-Y,  X,  Z+5/6'
3  '-Y,  X-Y,  Z+2/3'
4  '-X,  -Y,  Z+1/2'
5  '-X+Y,  -X,  Z+1/3'
6  'Y,  -X+Y,  Z+1/6'

loop_
_refln.index_h
_refln.index_k
_refln.index_l
_refln.pdbx_FWT
_refln.pdbx_PHWT
_refln.pdbx_DELFWT
_refln.pdbx_DELPHWT
0    0    6 1793.541 297.538 758.340 297.538
0    0   12 926.361   8.616  79.780   8.616
0    0   18 2789.494   3.915 158.391 183.915
0    0   24  18.715 196.581 133.229 196.581
0    0   30 1040.906 174.397 370.046 174.397
0    0   36 2065.720 180.659 163.982   0.659
0    0   42 1850.858 199.047 1496.942 199.047
0    0   48 882.296 190.106 289.614  10.106
0    1    6 161.907 235.025 332.891 235.025
0    1    7 168.822  39.648 511.322  39.648
0    1    8  42.963  55.666 347.014  55.666
0    1    9  81.428  68.620  89.690 248.620
"""

integer_observations = """
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
_refln.index_h
_refln.index_k
_refln.index_l
_refln.status
_refln.F_meas_au
_refln.F_meas_sigma_au
0    0    6 o     267   12
0    0   12 o    4700  113
0    0   18 o   10214  222
0    0   24 o    8268  192
0    0   30 o    3274   77
0    0   36 o     317   30
0    0   42 o   12026  286
0    0   48 o    1972   51
"""

p1_sym_ops = """\
data_r1e5xsf
_cell.entry_id          1E5X
_cell.length_a          57.760
_cell.length_b          62.140
_cell.length_c          76.590
_cell.angle_alpha       109.48
_cell.angle_beta        97.61
_cell.angle_gamma       112.74
_cell.formula_units_Z   2

_symmetry_equiv.id           1
_symmetry_equiv.pos_as_xyz   X,Y,Z
"""

def exercise_detect_binary():
  binary_string = '\xff\xf8\x00\x00\x00\x00\x00\x00'
  from iotbx.cif import reader
  try: reader(input_string=binary_string)
  except CifParserError as e: pass
  else: raise Exception_expected

def exercise_syntax_errors():
  empty_loop_str = """\
data_fred
loop_
loop_
_a
_b
1 2
3 4
"""
  try: cif.reader(input_string=empty_loop_str)
  except CifParserError as e: pass
  else: raise Exception_expected
  bad_semicolon_text_field = """\
data_sucrose
_a 1
_exptl_absorpt_process_details
;
Final HKLF 4 output contains 64446 reflections, Rint = 0.0650
 (47528 with I > 3sig(I), Rint = 0.0624);
"""
  try: cif.reader(input_string=bad_semicolon_text_field)
  except CifParserError as e: pass
  else: raise Exception_expected

def exercise_missing_atom_site_type_symbol():
  m = cif.reader(input_string=cif_with_atom_site_type_symbol).model()
  xs = cif.builders.crystal_structure_builder(m['9000000']).structure
  sc = xs.scatterers()
  assert len(sc) == 18
  assert sc[0].label == 'Fe1'
  assert sc[0].scattering_type == 'Fe'
  assert sc[-1].label == 'O7'
  assert sc[-1].scattering_type == 'O'

cif_with_atom_site_type_symbol = """\
#------------------------------------------------------------------------------
#$Date$
#$Revision$
#$URL: svn://www.crystallography.net/cod/cif/9/00/00/9000000.cif $
#------------------------------------------------------------------------------
#
# This file is available in the Crystallography Open Database (COD),
# http://www.crystallography.net/. The original data for this entry
# were provided the American Mineralogist Crystal Structure Database,
# http://rruff.geo.arizona.edu/AMS/amcsd.php
#
# The file may be used within the scientific community so long as
# proper attribution is given to the journal article from which the
# data were obtained.
#
data_9000000
loop_
_publ_author_name
'Finger, L. W.'
_publ_section_title
;
 The crystal structure and cation distribution of a grunerite
 Locality: Wabush iron formation, Labrador, Canada
;
_journal_name_full
'Mineralogical Society of America Special Paper'
_journal_page_first              95
_journal_page_last               100
_journal_volume                  2
_journal_year                    1969
_chemical_formula_sum            'F0.5 Fe6.1 H1.5 Mg0.9 O23.5 Si8'
_chemical_name_mineral           Grunerite
_space_group_IT_number           12
_symmetry_space_group_name_Hall  '-C 2y'
_symmetry_space_group_name_H-M   'C 1 2/m 1'
_cell_angle_alpha                90
_cell_angle_beta                 101.892
_cell_angle_gamma                90
_cell_length_a                   9.5642
_cell_length_b                   18.393
_cell_length_c                   5.3388
_cell_volume                     919.015
_exptl_crystal_density_diffrn    3.521
_[local]_cod_chemical_formula_sum_orig '(Fe6.1 Mg.9) Si8 O23.5 (F.5 H1.5)'
_cod_database_code               9000000
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
1/2+x,1/2+y,z
x,-y,z
1/2+x,1/2-y,z
-x,y,-z
1/2-x,1/2+y,-z
-x,-y,-z
1/2-x,1/2-y,-z
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
Fe1 0.00532 0.00943 0.00498 0.00000 0.00198 0.00000
Mg1 0.00532 0.00943 0.00498 0.00000 0.00198 0.00000
Fe2 0.00532 0.00823 0.00622 0.00000 0.00198 0.00000
Mg2 0.00532 0.00823 0.00622 0.00000 0.00198 0.00000
Fe3 0.00621 0.00857 0.00622 0.00000 0.00074 0.00000
Mg3 0.00621 0.00857 0.00622 0.00000 0.00074 0.00000
Fe4 0.00843 0.01645 0.01106 0.00000 0.00421 0.00000
Mg4 0.00843 0.01645 0.01106 0.00000 0.00421 0.00000
Si1 0.00399 0.00703 0.00553 -0.00026 0.00050 0.00000
Si2 0.00444 0.00651 0.00747 -0.00166 0.00074 0.00000
O1 0.00577 0.01028 0.00899 0.00000 0.00173 0.00146
O2 0.00399 0.01028 0.00940 -0.00087 0.00272 -0.00146
OH3 0.01864 0.01200 0.01272 0.00000 0.00570 0.00000
F3 0.01864 0.01200 0.01272 0.00000 0.00570 0.00000
O4 0.00887 0.00686 0.00567 -0.00087 -0.00050 0.00097
O5 0.00355 0.01371 0.01244 -0.00262 0.00223 0.00730
O6 0.00399 0.02057 0.00760 0.00349 -0.00025 -0.00438
O7 0.00843 0.00000 0.02019 0.00000 0.00520 0.00000
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
Fe1 0.00000 0.08781 0.50000 0.84800 0.00646
Mg1 0.00000 0.08781 0.50000 0.15200 0.00646
Fe2 0.00000 0.17936 0.00000 0.77300 0.00646
Mg2 0.00000 0.17936 0.00000 0.22700 0.00646
Fe3 0.00000 0.00000 0.00000 0.88800 0.00709
Mg3 0.00000 0.00000 0.00000 0.11200 0.00709
Fe4 0.00000 0.25741 0.50000 0.98500 0.01165
Mg4 0.00000 0.25741 0.50000 0.01500 0.01165
Si1 0.28670 0.08360 0.27070 1.00000 0.00557
Si2 0.29930 0.16670 0.77800 1.00000 0.00621
O1 0.11200 0.08820 0.20440 1.00000 0.00849
O2 0.12530 0.17350 0.71420 1.00000 0.00747
O-H3 0.11470 0.00000 0.70350 0.75000 0.01381
F3 0.11470 0.00000 0.70350 0.25000 0.01381
O4 0.38390 0.24160 0.76890 1.00000 0.00735
O5 0.34830 0.12750 0.05190 1.00000 0.00975
O6 0.34780 0.11820 0.55300 1.00000 0.01089
O7 0.33760 0.00000 0.27000 1.00000 0.00937
"""

def exercise():
  exercise_missing_atom_site_type_symbol()
  exercise_syntax_errors()
  exercise_detect_binary()
  exercise_crystal_symmetry()
  exercise_miller_arrays_as_cif_block()
  exercise_lex_parse_build()
  exercise_partial_crystal_symmetry()
  exercise_mmcif_structure_factors()
  exercise_atom_type_loop()

if __name__ == '__main__':
  exercise()
  print("OK")
