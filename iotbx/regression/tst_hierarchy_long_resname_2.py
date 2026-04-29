from __future__ import absolute_import, division, print_function
import time
from libtbx.test_utils import assert_lines_in_text
import iotbx.pdb
import iotbx.cif
import mmtbx.model
from mmtbx import monomer_library

# ------------------------------------------------------------------------------

# from https://github.com/wwPDB/extended-wwPDB-identifier-examples
# https://github.com/wwPDB/extended-wwPDB-identifier-examples/blob/main/Models/7fgz-extended_CCD_code-model.cif
mmcif_str = '''
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
HETATM 2180 O  O   . longHOH E 4 .   ? -3.209  -12.856 46.100 1.00 30.11 ? 401 longHOH A O   1
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
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 longHOH  longHOH  'water                    '  ligand  3  1  .  2012-10-17  2022-03-18
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as GRAND
;

data_comp_longHOH
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
 longHOH  O   O  OH2   0  -0.881  -23.1198  18.3926  -21.5996
 longHOH  H1  H  HOH2  0   0.440  -22.1601  18.3725  -21.6036
 longHOH  H2  H  HOH2  0   0.440  -23.3469  18.3903  -20.6669

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 longHOH  O  H1  single  0.850  0.020  0.980
 longHOH  O  H2  single  0.850  0.020  0.980

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 longHOH  H2  O  H1  103.91  3.000
"""

def test1():
  """
  Creating restraints for long residue name:
    HOH changed to longHOH, cif copied from modules/chem_data/geostd/h/data_HOH.cif
  """

  inp = iotbx.pdb.input(lines=mmcif_str.split("\n"), source_info=None)
  cif_object = iotbx.cif.reader(input_string = ligand_cif).model()
  cif_objects = [('bla.cif', cif_object)]
  model = mmtbx.model.manager(model_input = inp, restraint_objects = cif_objects)
  model.process(make_restraints=True)
  geo_str = model.restraints_as_geo()
  # print(geo_str)
  for l in [
      'bond pdb=" NZ  LYS A 279 "',
      '     pdb=" O   longHOH A 401 "',
      'nonbonded pdb=" CD  LYS A 279 "',
      '          pdb=" O   longHOH A 401 "']:
    assert_lines_in_text(geo_str, l)
  model_cif = model.model_as_mmcif()
  # print(model_cif)
  for l in [
      'longHOH H2 O H1 103.91 3.000',
      'data_comp_longHOH',
      '   HETATM 10 O . longHOH A 401 ? -3.20900 -12.85600 46.10000 1.000 30.11000 O ? B ? . O 1']:
    assert_lines_in_text(model_cif, l)
  model_pdb = model.model_as_pdb()
  assert not model.can_be_output_as_pdb()

if (__name__ == "__main__"):
  t0 = time.time()
  mon_lib_srv = None
  try:
    mon_lib_srv = monomer_library.server.server()
  except: # intentional
    print("Can not initialize monomer_library, skipping test.")
  if mon_lib_srv is not None:
    test1()
  print("OK. Time: %8.3f"%(time.time()-t0))
