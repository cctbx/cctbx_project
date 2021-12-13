from __future__ import absolute_import, division, print_function
import iotbx.cif
from libtbx.utils import Sorry
from libtbx.test_utils import show_diff


def test_repeated_loop_1():
  cif_txt_1 = """\
data_test
loop_
  _geom_bond_atom_site_label_1
  _geom_bond_atom_site_label_2
  _geom_bond_distance
  _geom_bond_site_symmetry_2
  Si  O  1.6160  4_554
  Si  O  1.6160  2_554
  """
  cif_txt_2 = """\
loop_
  _geom_bond_atom_site_label_1
  _geom_bond_atom_site_label_2
  _geom_bond_distance
  _geom_bond_site_symmetry_2
  Si  O  1.6160  3_664
  Si  O  1.6160  5_664
  """
  r = iotbx.cif.reader(input_string=cif_txt_1)
  try:
    r = iotbx.cif.reader(input_string=cif_txt_1+cif_txt_2)
  except Sorry as e:
    assert not show_diff(str(e), "Loop containing tags _geom_bond_atom_site_label_1, _geom_bond_atom_site_label_2, _geom_bond_distance, _geom_bond_site_symmetry_2 appears repeated")
  else:
    assert 0, "Sorry must be raised above."

def test_repeated_loop_2():
  cif_txt_1 = """\
data_test
loop_
  _atom_site_label
  _atom_site_type_symbol
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  Rb1  Rb           0.5           0         0.25
  Rb2  Rb    0.31563(3)         0.5   0.14226(3)

loop_
  _atom_site_label
  _atom_site_U_iso_or_equiv
  Rb1 0.5(2)
  Rb2 0.5(3)
  """
  try:
    r = iotbx.cif.reader(input_string=cif_txt_1)
  except Sorry as e:
    assert not show_diff(str(e), "Loop containing tags _atom_site_label, _atom_site_U_iso_or_equiv appears repeated")
  else:
    assert 0, "Sorry must be raised above."

def test_repeated_loop_3_with_typo():
  cif_txt_1 = """\
  data_test
loop_
  _atom_site_displace_Fourier_id
  _atom_site_displace_Fourier_atom_site_label
  _atom_site_displace_Fourier_axis
  _atom_site_displace_Fourier_wave_vector_seq_id
  Rb1y1 Rb1  y  1
  Rb2x1 Rb2  x  1
  Rb2y1 Rb2  y  1
  Rb2z1 Rb2  z  1

loop_
  _atom_site_displace_Fourier_parm_id
  _atom_site_displace_Fourier_param_cos
  _atom_site_displace_Fourier_param_sin
  Rb1y1   -0.0043(3)           0
  Rb2x1            0  -0.0001(1)
  Rb2y1   -0.0027(2)           0
  Rb2z1            0   0.0001(1)
  """
  try:
    r = iotbx.cif.reader(input_string=cif_txt_1)
  except Sorry as e:
    assert not show_diff(str(e), "Loop containing tags _atom_site_displace_Fourier_parm_id, _atom_site_displace_Fourier_param_cos, _atom_site_displace_Fourier_param_sin appears repeated")
  else:
    assert 0, "Sorry must be raised above."

def test_repeated_tags():
  cif_txt_1 = """\
data_RbMo
_audit_block_code                 RbMo

_diffrn_ambient_temperature 200
_diffrn_ambient_temperature 300
  """
  try:
    r = iotbx.cif.reader(input_string=cif_txt_1)
  except Sorry as e:
    assert not show_diff(str(e), "Data item _diffrn_ambient_temperature received multiple values")
  else:
    assert 0, "Sorry must be raised above."

if __name__ == '__main__':
  # from: https://github.com/cctbx/cctbx_project/issues/662
  test_repeated_loop_1()
  test_repeated_loop_2()
  test_repeated_loop_3_with_typo()
  test_repeated_tags()
