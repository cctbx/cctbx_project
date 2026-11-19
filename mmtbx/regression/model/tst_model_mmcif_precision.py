"""mmtbx.model.manager must keep mmCIF atom values at deposited precision.

Regression for GitHub issue #1269: the manager obtains its hierarchy as
deepcopy(model_input).construct_hierarchy(). For an mmCIF input the hierarchy
is already cached on the input by the time of that copy (crystal_symmetry()
builds it), so the copy went through root.__getstate__, i.e. through PDB
text, and every mmCIF model loaded through the manager or the DataManager had
coordinates rounded to 1e-3 A and occupancies / B to 1e-2 (the DataManager
builds its models through the same call). Pickling a manager
went the same way and left the hierarchy rounded while the xray_structure was
exact.
"""
import pickle
from copy import deepcopy

import iotbx.pdb
import mmtbx.model
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times

_CIF = """\
data_precise
_cell.length_a      10.000
_cell.length_b      10.000
_cell.length_c      10.000
_cell.angle_alpha   90.000
_cell.angle_beta    90.000
_cell.angle_gamma   90.000
_symmetry.space_group_name_H-M   'P 1'
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
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N  N  . GLY A 1 1 1.12345 2.23456 3.34567 0.987 20.1234 1 A 1
ATOM 2 C  CA . GLY A 1 1 2.00001 3.00002 4.00003 0.5   21.5555 1 A 1
"""

_EXPECTED = [
  ((1.12345, 2.23456, 3.34567), 0.987, 20.1234),
  ((2.00001, 3.00002, 4.00003), 0.5,   21.5555),
]


def _values(hierarchy):
  return [(a.xyz, a.occ, a.b) for a in hierarchy.atoms()]


def _make_input():
  return iotbx.pdb.input(source_info=None, lines=_CIF.splitlines())


def exercise_manager_from_cif_input():
  """manager(model_input=cif_input).get_hierarchy() must carry the deposited
  values, not PDB-rounded ones."""
  m = mmtbx.model.manager(model_input=_make_input())
  assert _values(m.get_hierarchy()) == _EXPECTED, _values(m.get_hierarchy())


def exercise_manager_copies_are_exact_and_consistent():
  """deepcopy / pickle of a manager must keep the hierarchy exact and in
  agreement with the xray_structure that travels with it."""
  m = mmtbx.model.manager(model_input=_make_input())
  m.get_xray_structure()
  for label, clone in (
      ("deep_copy", m.deep_copy()),
      ("deepcopy", deepcopy(m)),
      ("pickle", pickle.loads(pickle.dumps(m)))):
    h = clone.get_hierarchy()
    assert _values(h) == _EXPECTED, (label, _values(h))
    assert approx_equal(
      clone.get_xray_structure().sites_cart(), h.atoms().extract_xyz(),
      eps=1e-10), label


def run():
  exercise_manager_from_cif_input()
  exercise_manager_copies_are_exact_and_consistent()


if __name__ == "__main__":
  run()
  print(format_cpu_times())
  print("OK")
