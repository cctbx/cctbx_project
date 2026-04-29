from __future__ import absolute_import, division, print_function
import time
from libtbx.test_utils import assert_lines_in_text
import iotbx.pdb
import iotbx.cif
import mmtbx.model
from mmtbx import monomer_library
from libtbx.utils import Sorry
from six.moves import cStringIO as StringIO

# ------------------------------------------------------------------------------

# from https://github.com/wwPDB/extended-wwPDB-identifier-examples
# https://github.com/wwPDB/extended-wwPDB-identifier-examples/blob/main/Models/7fgz-extended_CCD_code-model.cif
mm_cif = '''
data_XXXX
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   2140 N  N   . LYS A 1 261 ? 0.399   -10.171 39.802 1.00 40.11 ? 279 LYS A N   1
ATOM   2141 C  CA  . LYS A 1 261 ? -0.169  -9.988  41.173 1.00 43.86 ? 279 LYS A CA  1
ATOM   2142 C  C   . LYS A 1 261 ? 0.687   -9.011  41.991 1.00 41.94 ? 279 LYS A C   1
ATOM   2143 O  O   . LYS A 1 261 ? 1.044   -7.920  41.556 1.00 39.32 ? 279 LYS A O   1
ATOM   2144 C  CB  . LYS A 1 261 ? -0.260  -11.336 41.902 1.00 46.47 ? 279 LYS A CB  1
ATOM   2145 C  CG  . LYS A 1 261 ? -1.583  -12.074 41.713 1.00 49.13 ? 279 LYS A CG  1
ATOM   2146 C  CD  . LYS A 1 261 ? -1.611  -13.468 42.315 1.00 51.03 ? 279 LYS A CD  1
ATOM   2147 C  CE  . LYS A 1 261 ? -2.923  -13.799 42.993 1.00 52.86 ? 279 LYS A CE  1
ATOM   2148 N  NZ  . LYS A 1 261 ? -3.209  -12.856 44.100 1.00 54.19 ? 279 LYS A NZ  1
HETATM 2150 N  N1  . 7ZTVU C 3 .   ? -7.743  -6.355  8.243  1.00 18.72 ? 302 7ZTVU A N1  1
HETATM 2151 C  C2  . 7ZTVU C 3 .   ? -8.462  -5.534  9.265  1.00 16.68 ? 302 7ZTVU A C2  1
HETATM 2152 C  C3  . 7ZTVU C 3 .   ? -8.092  -5.865  10.711 1.00 17.35 ? 302 7ZTVU A C3  1
HETATM 2153 N  N4  . 7ZTVU C 3 .   ? -7.975  -7.334  10.767 1.00 17.04 ? 302 7ZTVU A N4  1
HETATM 2154 C  C5  . 7ZTVU C 3 .   ? -6.781  -7.689  10.027 1.00 17.11 ? 302 7ZTVU A C5  1
HETATM 2155 C  C6  . 7ZTVU C 3 .   ? -7.363  -7.730  8.633  1.00 16.97 ? 302 7ZTVU A C6  1
HETATM 2156 C  C7  . 7ZTVU C 3 .   ? -8.071  -7.951  12.064 1.00 16.77 ? 302 7ZTVU A C7  1
HETATM 2157 C  C8  . 7ZTVU C 3 .   ? -8.490  -9.440  12.025 1.00 16.60 ? 302 7ZTVU A C8  1
HETATM 2158 O  O8  . 7ZTVU C 3 .   ? -8.391  -10.140 10.781 1.00 14.14 ? 302 7ZTVU A O8  1
HETATM 2159 C  C9  . 7ZTVU C 3 .   ? -8.438  -6.311  6.935  1.00 20.13 ? 302 7ZTVU A C9  1
HETATM 2160 C  C10 . 7ZTVU C 3 .   ? -7.646  -6.965  5.796  1.00 22.62 ? 302 7ZTVU A C10 1
HETATM 2161 S  S   . 7ZTVU C 3 .   ? -8.372  -6.708  4.279  1.00 24.90 ? 302 7ZTVU A S   1
HETATM 2162 O  O1S . 7ZTVU C 3 .   ? -7.579  -7.371  3.224  1.00 25.47 ? 302 7ZTVU A O1S 1
HETATM 2163 O  O2S . 7ZTVU C 3 .   ? -8.533  -5.211  4.078  1.00 25.01 ? 302 7ZTVU A O2S 1
HETATM 2164 O  O3S . 7ZTVU C 3 .   ? -9.691  -7.326  4.282  1.00 26.42 ? 302 7ZTVU A O3S 1
HETATM 2165 N  N1  . 7ZTVU D 3 .   ? -11.233 1.732   10.446 1.00 74.61 ? 303 7ZTVU A N1  1
HETATM 2166 C  C2  . 7ZTVU D 3 .   ? -11.682 2.495   9.264  1.00 77.03 ? 303 7ZTVU A C2  1
HETATM 2167 C  C3  . 7ZTVU D 3 .   ? -11.907 3.979   9.617  1.00 76.82 ? 303 7ZTVU A C3  1
HETATM 2168 N  N4  . 7ZTVU D 3 .   ? -12.335 4.282   11.022 1.00 73.77 ? 303 7ZTVU A N4  1
HETATM 2169 C  C5  . 7ZTVU D 3 .   ? -12.654 3.123   11.908 1.00 69.90 ? 303 7ZTVU A C5  1
HETATM 2170 C  C6  . 7ZTVU D 3 .   ? -12.333 1.719   11.409 1.00 69.07 ? 303 7ZTVU A C6  1
HETATM 2171 C  C7  . 7ZTVU D 3 .   ? -11.415 5.232   11.725 1.00 71.46 ? 303 7ZTVU A C7  1
HETATM 2172 C  C8  . 7ZTVU D 3 .   ? -12.126 5.976   12.871 1.00 70.74 ? 303 7ZTVU A C8  1
HETATM 2173 O  O8  . 7ZTVU D 3 .   ? -11.360 5.921   14.096 1.00 61.89 ? 303 7ZTVU A O8  1
HETATM 2174 C  C9  . 7ZTVU D 3 .   ? -10.756 0.357   10.142 1.00 74.40 ? 303 7ZTVU A C9  1
HETATM 2175 C  C10 . 7ZTVU D 3 .   ? -10.121 -0.289  11.396 1.00 75.92 ? 303 7ZTVU A C10 1
HETATM 2176 S  S   . 7ZTVU D 3 .   ? -8.757  -1.214  11.091 1.00 77.53 ? 303 7ZTVU A S   1
HETATM 2177 O  O1S . 7ZTVU D 3 .   ? -7.735  -0.305  10.525 1.00 74.01 ? 303 7ZTVU A O1S 1
HETATM 2178 O  O2S . 7ZTVU D 3 .   ? -9.042  -2.313  10.141 1.00 83.78 ? 303 7ZTVU A O2S 1
HETATM 2179 O  O3S . 7ZTVU D 3 .   ? -8.264  -1.864  12.321 1.00 70.44 ? 303 7ZTVU A O3S 1
'''

ligand_cif = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 7ZTVU         7ZTVU '2-[4-(2-hydroxyethyl)piperazin-1-yl]ethanesulfonic acid' ligand 32 15 .
#
data_comp_7ZTVU
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 7ZTVU         N1     N   NT     0    .      -0.5253    0.0932   -1.5549
 7ZTVU         C2     C   CH2    0    .       0.9262    0.1022   -1.5401
 7ZTVU         C3     C   CH2    0    .       1.4857    0.1109   -0.1192
 7ZTVU         N4     N   NT     0    .       0.9247    1.1836    0.6816
 7ZTVU         C5     C   CH2    0    .      -0.5267    1.1746    0.6669
 7ZTVU         C6     C   CH2    0    .      -1.0863    1.1660   -0.7540
 7ZTVU         C7     C   CH2    0    .       1.3954    1.0634    2.0495
 7ZTVU         C8     C   CH2    0    .       2.1337    2.3395    2.4478
 7ZTVU         O8     O   OH1    0    .       2.5312    2.2499    3.7877
 7ZTVU         C9     C   CH2    0    .      -1.0063   -1.1833   -1.0588
 7ZTVU         C10    C   CH2    0    .      -1.8024   -1.8881   -2.1549
 7ZTVU         S      S   S      0    .      -2.4815   -3.4426   -1.5057
 7ZTVU         O1S    O   OS     0    .      -3.5911   -3.1842   -0.5146
 7ZTVU         O2S    O   OS     0    .      -3.2186   -4.2088   -2.5779
 7ZTVU         O3S    O   OS    -1    .      -1.3827   -4.2775   -0.8928
 7ZTVU         H21    H   HCH2   0    .       1.2877   -0.7819   -2.0541
 7ZTVU         H22    H   HCH2   0    .       1.2768    0.9870   -2.0606
 7ZTVU         H31    H   HCH2   0    .       1.2556   -0.8376    0.3543
 7ZTVU         H32    H   HCH2   0    .       2.5623    0.2348   -0.1669
 7ZTVU         H51    H   HCH2   0    .      -0.8883    2.0588    1.1809
 7ZTVU         H52    H   HCH2   0    .      -0.8774    0.2898    1.1874
 7ZTVU         H61    H   HCH2   0    .      -0.8561    2.1144   -1.2275
 7ZTVU         H62    H   HCH2   0    .      -2.1629    1.0420   -0.7063
 7ZTVU         H71    H   HCH2   0    .       2.0691    0.2167    2.1269
 7ZTVU         H72    H   HCH2   0    .       0.5490    0.9138    2.7113
 7ZTVU         H81    H   HCH2   0    .       1.4752    3.1923    2.3222
 7ZTVU         H82    H   HCH2   0    .       3.0087    2.4624    1.8185
 7ZTVU         HO8    H   HOH1   0    .       3.0956    2.9773    3.9980
 7ZTVU         H91    H   HCH2   0    .      -0.1620   -1.8008   -0.7712
 7ZTVU         H92    H   HCH2   0    .      -1.6444   -1.0192   -0.1971
 7ZTVU        H101    H   HCH2   0    .      -1.1503   -2.1015   -2.9952
 7ZTVU        H102    H   HCH2   0    .      -2.6144   -1.2471   -2.4811
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 7ZTVU   N1      C2    single        1.452 0.020
 7ZTVU   N1      C6    single        1.452 0.020
 7ZTVU   N1      C9    single        1.452 0.020
 7ZTVU   C2      C3    single        1.527 0.020
 7ZTVU   C2      H21   single        1.085 0.020
 7ZTVU   C2      H22   single        1.085 0.020
 7ZTVU   C3      N4    single        1.452 0.020
 7ZTVU   C3      H31   single        1.085 0.020
 7ZTVU   C3      H32   single        1.085 0.020
 7ZTVU   N4      C5    single        1.452 0.020
 7ZTVU   N4      C7    single        1.452 0.020
 7ZTVU   C5      C6    single        1.527 0.020
 7ZTVU   C5      H51   single        1.085 0.020
 7ZTVU   C5      H52   single        1.085 0.020
 7ZTVU   C6      H61   single        1.085 0.020
 7ZTVU   C6      H62   single        1.085 0.020
 7ZTVU   C7      C8    single        1.527 0.020
 7ZTVU   C7      H71   single        1.085 0.020
 7ZTVU   C7      H72   single        1.085 0.020
 7ZTVU   C8      O8    single        1.401 0.020
 7ZTVU   C8      H81   single        1.085 0.020
 7ZTVU   C8      H82   single        1.085 0.020
 7ZTVU   O8      HO8   single        0.944 0.020
 7ZTVU   C9      C10   single        1.527 0.020
 7ZTVU   C9      H91   single        1.085 0.020
 7ZTVU   C9      H92   single        1.085 0.020
 7ZTVU   C10     S     single        1.816 0.020
 7ZTVU   C10    H101   single        1.085 0.020
 7ZTVU   C10    H102   single        1.085 0.020
 7ZTVU   S       O1S   deloc         1.510 0.020
 7ZTVU   S       O2S   deloc         1.510 0.020
 7ZTVU   S       O3S   deloc         1.510 0.020
"""

def test1():
  """
  Creating restraints for long residue name
  """
  #dumping to disk if needed:
  # for name, s in [('model', mm_cif), ('restr', ligand_cif)]:
  #   with open('%s.cif' % name, 'w') as f:
  #     f.write(s)
  inp = iotbx.pdb.input(lines=mm_cif.split("\n"), source_info=None)
  cif_object = iotbx.cif.reader(input_string = ligand_cif).model()
  cif_objects = [('bla.cif', cif_object)]
  model = mmtbx.model.manager(
      model_input = inp,
      restraint_objects = cif_objects)
  try:
    model.process(make_restraints=True)
    geo_str = model.restraints_as_geo()
  except Sorry as e:
    pass
  # print(geo_str)
  for l in [
    'bond pdb=" C10 7ZTVU A 303 "',
    '     pdb=" S   7ZTVU A 303 "',
    'nonbonded pdb=" C2  7ZTVU A 302 "',
    '          pdb=" O8  7ZTVU A 302 "',
    ]:
    # print(geo_str.find(l),l)
    assert_lines_in_text(geo_str, l)
  model_cif = model.model_as_mmcif()
   # print(model_cif)
  for l in [
    'HETATM 10 C10 . 7ZTVU A 302 ? -7.64600 -6.96500 5.79600 1.000 22.62000 C ? B ? . C10 1',
    '7ZTVU N1 C2 single 1.452 0.020',
    'data_comp_7ZTVU',
    ]:
    assert_lines_in_text(model_cif, l)
  assert not model.can_be_output_as_pdb()

def test2():
  """
  Try creating restraints for long residue name without cif.
  Check error message formatting (atom.quote() function)
  """
  #dumping to disk if needed:
  # for name, s in [('model', mm_cif), ('restr', ligand_cif)]:
  #   with open('%s.cif' % name, 'w') as f:
  #     f.write(s)
  inp = iotbx.pdb.input(lines=mm_cif.split("\n"), source_info=None)
  cif_object = iotbx.cif.reader(input_string = ligand_cif).model()
  cif_objects = [('bla.cif', cif_object)]
  mlog = StringIO()
  model = mmtbx.model.manager(
      model_input = inp,
      log = mlog)
  try:
    model.process(make_restraints=True)
  except Sorry as e:
    mlog_txt = mlog.getvalue()
    # print(str(e))
    # print(mlog_txt)
    assert_lines_in_text(mlog_txt, """ "HETATM   15  C7  7ZTVU A 302 .*.     C  " """)


if (__name__ == "__main__"):
  t0 = time.time()
  mon_lib_srv = None
  try:
    mon_lib_srv = monomer_library.server.server()
  except: # intentional
    print("Can not initialize monomer_library, skipping test.")
  if mon_lib_srv is not None:
    test1()
    test2()
  print("OK. Time: %8.3f"%(time.time()-t0))
