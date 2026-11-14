"""
iotbx.pdb.mmcif.xray_structure_from_atom_site against the hierarchy route,
hierarchy_from_atom_site(block).extract_xray_structure(): the same scatterers
(label, scattering type, site, U, flags, occupancy, fp, fdp, multiplicity),
compared in label order since the direct route keeps the row order, over the
fixtures and model files of tst_mmcif_hierarchy_parity.py.
"""
import os
import libtbx.load_env
from libtbx.test_utils import Exception_expected
import iotbx.cif
from iotbx.cif.builders import CifBuilderError
from iotbx.pdb import mmcif
from iotbx.pdb import tst_mmcif_hierarchy as fixtures
from iotbx.regression import tst_mmcif_hierarchy_parity as parity


def scatterer_rows(xray_structure):
  result = []
  for sc in xray_structure.scatterers():
    result.append((sc.label, sc.scattering_type, tuple(sc.site), sc.u_iso,
      tuple(sc.u_star), sc.flags.use_u_iso(), sc.flags.use_u_aniso(),
      sc.occupancy, sc.fp, sc.fdp, sc.multiplicity()))
  return sorted(result)


def compare_block(block, label):
  hierarchy = mmcif.hierarchy_from_atom_site(block)
  # a block without a cell gets the box extract_xray_structure would make
  crystal_symmetry = hierarchy.generate_crystal_symmetry(
    mmcif.crystal_symmetry_builder(block).crystal_symmetry)
  ref = hierarchy.extract_xray_structure(crystal_symmetry=crystal_symmetry)
  new = mmcif.xray_structure_from_atom_site(block,
    crystal_symmetry=crystal_symmetry)
  assert new.crystal_symmetry().is_similar_symmetry(ref.crystal_symmetry())
  if scatterer_rows(new) != scatterer_rows(ref):
    raise AssertionError("%s: scatterers differ" % label)
  return new


def exercise_fixtures():
  for name, text in parity.parity_fixtures():
    compare_block(fixtures.first_block(text), "fixture " + name)


def exercise_model_files():
  n_aniso = 0
  for module, rel in parity.PARITY_FILES:
    fn = os.path.join(libtbx.env.dist_path(module), rel)
    assert os.path.isfile(fn), fn
    model = iotbx.cif.reader(file_path=fn).model()
    new = compare_block(list(model.values())[0], os.path.basename(fn))
    n_aniso += new.use_u_aniso().count(True)
  assert n_aniso > 0


scattering_type_header = """\
data_test
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
"""

# a blank element takes its type from the name (two names, so the type
# lookup must not be shared between them), IS with a charge, an atom name
# with the "const" scattering type
scattering_type_rows = """\
HETATM 1 ' ' FE   . HEM B 2 1 ? 1.0 2.0 3.0 1.0 20.0 ? 1 HEM A FE   1
HETATM 2 ' ' CA   . HEM B 2 1 ? 2.0 2.0 3.0 1.0 20.0 ? 1 HEM A CA   1
HETATM 3 IS  IS1  . UNL C 3 2 ? 3.0 2.0 3.0 1.0 20.0 2 2 UNL A IS1  1
HETATM 4 C   PEAK . UNL C 3 2 ? 4.0 2.0 3.0 1.0 20.0 ? 2 UNL A PEAK 1
"""

unknown_element_row = """\
HETATM 5 Xx  X1   . UNL C 3 2 ? 5.0 2.0 3.0 1.0 20.0 ? 2 UNL A X1   1
"""


def exercise_scattering_types():
  block = fixtures.first_block(scattering_type_header + scattering_type_rows)
  new = compare_block(block, "scattering types")
  assert sorted(new.scatterers().extract_scattering_types()) == [
    "Ca", "Fe", "IS2", "const"]
  block = fixtures.first_block(
    scattering_type_header + scattering_type_rows + unknown_element_row)
  hierarchy = mmcif.hierarchy_from_atom_site(block)
  crystal_symmetry = hierarchy.generate_crystal_symmetry(None)
  ref = hierarchy.extract_xray_structure(crystal_symmetry=crystal_symmetry,
    enable_scattering_type_unknown=True)
  new = mmcif.xray_structure_from_atom_site(block,
    crystal_symmetry=crystal_symmetry, enable_scattering_type_unknown=True)
  assert scatterer_rows(new) == scatterer_rows(ref)
  assert "unknown" in new.scatterers().extract_scattering_types()
  try:
    mmcif.xray_structure_from_atom_site(block,
      crystal_symmetry=crystal_symmetry)
  except CifBuilderError as e:
    assert str(e).startswith("Unknown chemical element type")
    assert '"Xx"' in str(e)
  else:
    raise Exception_expected


def exercise_errors():
  block = fixtures.first_block(
    fixtures.atom_site_header + fixtures.input_two_atoms)
  try:
    mmcif.xray_structure_from_atom_site(block)
  except CifBuilderError as e:
    assert str(e) == "crystal symmetry is missing"
  else:
    raise Exception_expected


def run():
  exercise_fixtures()
  exercise_model_files()
  exercise_scattering_types()
  exercise_errors()


if __name__ == "__main__":
  run()
  print("OK")
