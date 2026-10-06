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
_chem_comp_atom.charge
IA TA EA 0.0 .
IB TB EB 1.0 .
IC TC EC 2.0 .

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
  assert s.getvalue().splitlines()[9] == "ID TD TD 5.0 ."
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
  assert s.getvalue().splitlines()[6] == "ID TD ED 5.0 ."
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
  assert s.getvalue().splitlines()[6] == "IA TD ED 5.0 ."
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

charge_cif = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 CHG  CHG  'charges' ligand 6 6 .
 NOC  NOC  'no charge column' ligand 2 2 .

data_comp_CHG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
 CHG  A1  O  O   -1  -0.5
 CHG  A2  C  C    0   0.1
 CHG  A3  N  N   +1   0.4
 CHG  A4  C  C    ?   0.0
 CHG  A5  C  C    .   0.0
 CHG  A6  O  O   1.0  0.0

data_comp_NOC
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
 NOC  B1  C  C   0.0
 NOC  B2  O  O   0.0
"""

def exercise_charge():
  """_chem_comp_atom.charge: parsed as written (untyped), copied, written back."""
  import copy
  import iotbx.cif
  from mmtbx.monomer_library import server
  comps = dict([(c.chem_comp.id, c) for c in server.convert_comp_list(
    source_info=None, cif_object=iotbx.cif.reader(input_string=charge_cif).model())])
  chg, noc = comps["CHG"], comps["NOC"]
  assert [a.charge for a in chg.atom_list] == ["-1", "0", "+1", "?", "", "1.0"]
  assert [a.formal_charge() for a in chg.atom_list] == [-1, 0, 1, None, None, 1]
  assert [a.charge for a in noc.atom_list] == [None, None]
  assert [a.formal_charge() for a in noc.atom_list] == [None, None]
  assert [cif_types.formal_charge_and_problem(v) for v in ("0.5", "1-", "x")] == [
    (None, "charge '0.5' not integral"), (-1, None), (None, "charge 'x' not a number")]
  # copy keeps the value
  assert [copy.copy(a).charge for a in chg.atom_list] == [a.charge for a in chg.atom_list]
  # written and read back: '.', '' and a missing column are written as '.'
  s = StringIO()
  chg.show(s)
  lines = s.getvalue().splitlines()
  assert "_chem_comp_atom.charge" in lines
  assert [l.split()[-1] for l in lines[6:12]] == ["-1", "0", "+1", "?", ".", "1.0"], lines
  block = "data_comp_list\nloop_\n_chem_comp.id\n_chem_comp.three_letter_code\n" \
    "_chem_comp.name\n_chem_comp.group\n_chem_comp.number_atoms_all\n" \
    "_chem_comp.number_atoms_nh\n_chem_comp.desc_level\n CHG CHG x ligand 6 6 .\n\n" \
    "data_comp_CHG\n" + s.getvalue()
  (again,) = list(server.convert_comp_list(source_info=None,
    cif_object=iotbx.cif.reader(input_string=block).model()))
  assert [a.charge for a in again.atom_list] == ["-1", "0", "+1", "?", "", "1.0"]
  s = StringIO()
  noc.show(s)
  assert [l.split()[-1] for l in s.getvalue().splitlines()[6:8]] == [".", "."]
  # a chem_mod "add" atom has no charge; "change" keeps the original's
  mod = cif_types.mod_mod_id(source_info=None, chem_mod=cif_types.chem_mod(id="M",
    name="", comp_id="CHG", group_id=""))
  mod.atom_list.append(cif_types.chem_mod_atom(function="change", atom_id="A1",
    new_atom_id="", new_type_symbol="", new_type_energy="OC", new_partial_charge=""))
  mod.atom_list.append(cif_types.chem_mod_atom(function="add", atom_id="",
    new_atom_id="A7", new_type_symbol="H", new_type_energy="H", new_partial_charge=0))
  c = chg.apply_mod(mod)
  d = dict([(a.atom_id, a) for a in c.atom_list])
  assert d["A1"].charge == "-1" and d["A1"].type_energy == "OC"
  assert d["A7"].charge is None

if (__name__ == "__main__"):
  exercise()
  exercise_charge()
  print("OK")
