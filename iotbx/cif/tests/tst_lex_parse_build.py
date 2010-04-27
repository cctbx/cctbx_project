from cctbx.array_family import flex
from iotbx import cif
import libtbx.load_env
from libtbx.test_utils import show_diff
from libtbx.utils import time_log
from cStringIO import StringIO

def exercise():
  readers = [cif.python_reader]
  #if libtbx.env.has_module('antlr'):
    #readers.append(cif.fast_reader)
  #else:
    #print "Skipping compiled CIF reader tests"
  builders = [cif.builders.cif_model_builder]
  if libtbx.env.has_module('PyCifRW'):
    builders.append(cif.builders.PyCifRW_model_builder)
  else:
    print "Skipping PyCifRW builder tests"
  for reader in readers:
    for builder in builders:
      cif_model = reader(
        input_string=cif_xray_structure, builder=builder()).model()
      xs_builder = cif.builders.crystal_structure_builder(cif_model['global'])
      xs1 = xs_builder.structure
      # also test construction of cif model from xray structure
      xs_cif_block = xs1.as_cif_block()
      xs2 = cif.builders.crystal_structure_builder(xs_cif_block).structure
      for xs in (xs1, xs2):
        sio = StringIO()
        xs.show_scatterers(sio)
        assert not show_diff(sio.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso, Ustar as Uiso
o    O      2 ( 0.5000  0.0000  0.0000) 0.80 0.1000 [ - ]
c    C      2 ( 0.0000  0.0000  0.0000) 1.00 [ - ] 0.2000
     u_cart =  0.100  0.200  0.300  0.000  0.000  0.000
""")
        sio = StringIO()
        xs.show_summary(sio)
        assert not show_diff(sio.getvalue(), """\
Number of scatterers: 2
At special positions: 2
Unit cell: (10, 20, 30, 90, 90, 90)
Space group: C 1 2/m 1 (No. 12)
""")
      #
      cif_model = reader(
        input_string=cif_miller_array, builder=builder()).model()
      ma_builder = cif.builders.miller_array_builder(cif_model['global'])
      ma1 = ma_builder.array
      # also test construction of cif model from miller array
      #ma_cif_block = cif.miller_array_as_cif_block(ma1).cif_block
      #ma2 = cif.builders.miller_array_builder(ma_cif_block).array
      for ma in (ma1,):
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
Miller array info: None
Observation type: xray.intensity
Type of data: double, size=11
Type of sigmas: double, size=11
Number of Miller indices: 11
Anomalous flag: None
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
"""

cif_miller_array = """\
#
# h,k,l, Fc-squared, Fo-squared, sigma(Fo-squared) and status flag
#
data_global
_exptl_crystal_F_000       570.00
_reflns_d_resolution_high  0.7701

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

_shelx_F_squared_multiplier     1.000

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

if __name__ == '__main__':
  exercise()
  print "OK"
