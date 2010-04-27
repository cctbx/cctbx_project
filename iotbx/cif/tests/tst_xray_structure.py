import iotbx.cif.xray_structure
from cctbx import xray
from cctbx import crystal
from cctbx import adptbx
from cctbx.array_family import flex
from libtbx.test_utils import show_diff
from cStringIO import StringIO
import sys

def exercise_format_float():
  ff = iotbx.cif.xray_structure.format_float
  assert ff(value=None) == "."
  assert ff(value=0) == "0"
  assert ff(value=1) == "1"
  assert ff(value=-1) == "-1"
  assert ff(value=100000) == "100000"
  assert ff(value=-100000) == "-100000"
  assert ff(value=1000000) == "1+6"
  assert ff(value=-1000000) == "-1+6"
  assert ff(value=1e-6) == "1-6"
  assert ff(value=-1e-6) == "-1-6"
  assert ff(value=1e12) == "1+12"
  assert ff(value=-1e12) == "-1+12"
  assert ff(value=1e-12) == "1-12"
  assert ff(value=-1e-12) == "-1-12"
  assert ff(value=1.2345678) == "1.23457"
  expected = """\
1.23457-12
1.23457-11
1.23457-10
1.23457-9
1.23457-8
1.23457-7
1.23457-6
1.23457-5
0.000123457
0.00123457
0.0123457
0.123457
1.23457
12.3457
123.457
1234.57
12345.7
123457
1.23457+6
1.23457+7
1.23457+8
1.23457+9
1.23457+10
1.23457+11
1.23457+12
""".splitlines()
  for e in xrange(-12,13):
    assert not show_diff(ff(value=1.2345678*10**e), expected[e+12])
    assert not show_diff(ff(value=-1.2345678*10**e), "-"+expected[e+12])

def exercise_as_cif_simple():
  cs = crystal.symmetry((10, 20, 30, 90, 90, 90), "C 2/m")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer((
    xray.scatterer(
      label="o",
      site=(0.5, 0, 0),
      u=0.1,
      occupancy=0.8),
    xray.scatterer(
      label="c",
      site=(0, 0, 0),
      u=adptbx.u_cart_as_u_star(cs.unit_cell(), (0.1,0.2,0.3,0,0,0)))))
  xs = xray.structure(sp, scatterers)
  sio = StringIO()
  xs.as_cif_simple(out=sio)
  assert not show_diff(sio.getvalue(), """\
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
loop_
    _atom_site_label
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
""")

def run(args):
  assert len(args) == 0
  exercise_format_float()
  exercise_as_cif_simple()
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
