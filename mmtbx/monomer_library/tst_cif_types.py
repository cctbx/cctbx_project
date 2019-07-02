from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import cif_types
from six.moves import cStringIO as StringIO

def exercise():
  chem_comp = cif_types.chem_comp(
    id="Id",
    three_letter_code="TLC",
    name="Name",
    group="Group",
    number_atoms_all=22,
    number_atoms_nh=11,
    desc_level="")
  comp_comp_id = cif_types.comp_comp_id(source_info=None, chem_comp=chem_comp)
  for i,a in enumerate("ABC"):
    comp_comp_id.atom_list.append(cif_types.chem_comp_atom(
      atom_id="I"+a,
      type_symbol="T"+a,
      type_energy="E"+a,
      partial_charge=i))
  comp_comp_id.bond_list.append(cif_types.chem_comp_bond(
   atom_id_1="IA",
   atom_id_2="IC",
   type="single",
   value_dist="1",
   value_dist_esd="2"))
  comp_comp_id.bond_list.append(cif_types.chem_comp_bond(
   atom_id_1="IB",
   atom_id_2="IC",
   type="double",
   value_dist="3",
   value_dist_esd="4"))
  s = StringIO()
  comp_comp_id.show(s)
  assert s.getvalue() == """\
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
IA TA EA 0.0
IB TB EB 1.0
IC TC EC 2.0

loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
IA IC single 1.0 2.0 .
IB IC double 3.0 4.0 .

"""
  chem_mod = cif_types.chem_mod(
    id="MI",
    name="Name",
    comp_id="Comp Id",
    group_id="Group Id")
  mod_mod_id = cif_types.mod_mod_id(source_info=None, chem_mod=chem_mod)
  mod_mod_id.atom_list.append(cif_types.chem_mod_atom(
    function="add",
    atom_id="",
    new_atom_id="ID",
    new_type_symbol="TD",
    new_type_energy="TD",
    new_partial_charge=5))
  c = comp_comp_id.apply_mod(mod_mod_id)
  assert len(c.atom_list) == 4
  assert len(c.bond_list) == 2
  s = StringIO()
  c.show(s)
  assert s.getvalue().splitlines()[8] == "ID TD TD 5.0"
  mod_mod_id.atom_list[0] = cif_types.chem_mod_atom(
    function="change",
    atom_id="IA",
    new_atom_id="ID",
    new_type_symbol="TD",
    new_type_energy="ED",
    new_partial_charge=5)
  c = comp_comp_id.apply_mod(mod_mod_id)
  assert len(c.atom_list) == 3
  assert len(c.bond_list) == 2
  s = StringIO()
  c.show(s)
  assert s.getvalue().splitlines()[5] == "ID TD ED 5.0"
  mod_mod_id.atom_list[0] = cif_types.chem_mod_atom(
    function="change",
    atom_id="IA",
    new_atom_id="IA",
    new_type_symbol="TD",
    new_type_energy="ED",
    new_partial_charge=5)
  c = comp_comp_id.apply_mod(mod_mod_id)
  assert len(c.atom_list) == 3
  assert len(c.bond_list) == 2
  s = StringIO()
  c.show(s)
  assert s.getvalue().splitlines()[5] == "IA TD ED 5.0"
  mod_mod_id.atom_list[0] = cif_types.chem_mod_atom(
    function="delete",
    atom_id="IC",
    new_atom_id="",
    new_type_symbol="",
    new_type_energy="",
    new_partial_charge="")
  c = comp_comp_id.apply_mod(mod_mod_id)
  assert len(c.atom_list) == 2
  assert len(c.bond_list) == 0
  mod_mod_id.atom_list = []
  mod_mod_id.bond_list.append(cif_types.chem_mod_bond(
    function="add",
    atom_id_1="IA",
    atom_id_2="IB",
    new_type="triple",
    new_value_dist=5,
    new_value_dist_esd=6))
  c = comp_comp_id.apply_mod(mod_mod_id)
  assert len(c.atom_list) == 3
  assert len(c.bond_list) == 3
  s = StringIO()
  c.show(s)
  assert s.getvalue().splitlines()[-2] == "IA IB triple 5.0 6.0 ."
  mod_mod_id.bond_list[0] = cif_types.chem_mod_bond(
    function="change",
    atom_id_1="IA",
    atom_id_2="IC",
    new_type="quadruple",
    new_value_dist=7,
    new_value_dist_esd=8)
  c = comp_comp_id.apply_mod(mod_mod_id)
  s = StringIO()
  c.show(s)
  assert s.getvalue().splitlines()[-3] == "IA IC quadruple 7.0 8.0 ."
  print("OK")

if (__name__ == "__main__"):
  exercise()
