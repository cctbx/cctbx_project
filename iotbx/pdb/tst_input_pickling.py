"""Pickling / deepcopy regression tests for the PDB and mmCIF input classes.

Regression for: ``RuntimeError: Incomplete pickle support
(__getstate_manages_dict__ not set)``.

``iotbx_pdb_ext.input`` is pickled by value (Boost.Python ``enable_pickling``,
reconstructed from ``__getinitargs__``). As soon as ``scale_matrix()`` /
``xray_structure_simple()`` caches ``_scale_matrix`` on the instance, the
instance ``__dict__`` becomes non-empty, and Boost.Python refused to
pickle/deepcopy the object unless the class declares
``__getstate_manages_dict__``. Because ``mmtbx.model.manager`` deepcopies its
``model_input`` in the constructor, ``model_manager(pdb_inp)`` crashed whenever
``xray_structure_simple()`` had been called on the input first.

The ``cif_input`` test is a forward-looking guard. ``iotbx.pdb.mmcif.cif_input``
is pure Python today, so it pickles its ``__dict__`` natively and the test
passes trivially. It is planned to move to C++; when it does it will hit the
same Boost.Python limitation, and whoever makes that move must add equivalent
dict-managing pickle support so this test keeps passing.
"""

import pickle
from copy import deepcopy

import iotbx.pdb
import iotbx.pdb.mmcif
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times

# P1, a=b=c=10 (fractionalization matrix = diag(0.1)); two atoms. The SCALE
# records are what make scale_matrix() cache a real _scale_matrix.
_PDB = "\n".join([
  "CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1",
  "SCALE1      0.100000  0.000000  0.000000        0.00000",
  "SCALE2      0.000000  0.100000  0.000000        0.00000",
  "SCALE3      0.000000  0.000000  0.100000        0.00000",
  "ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00 20.00           N",
  "ATOM      2  CA  GLY A   1       2.000   3.000   4.000  1.00 20.00           C",
  "END",
])

# Same model in mmCIF. The _atom_sites.fract_transf_matrix/_vector records are
# what make cif_input.scale_matrix() cache a real _scale_matrix.
_CIF = """\
data_test
_cell.length_a      10.000
_cell.length_b      10.000
_cell.length_c      10.000
_cell.angle_alpha   90.000
_cell.angle_beta    90.000
_cell.angle_gamma   90.000
_symmetry.space_group_name_H-M   'P 1'
_symmetry.Int_Tables_number      1
_atom_sites.fract_transf_matrix[1][1]   0.100000
_atom_sites.fract_transf_matrix[1][2]   0.000000
_atom_sites.fract_transf_matrix[1][3]   0.000000
_atom_sites.fract_transf_matrix[2][1]   0.000000
_atom_sites.fract_transf_matrix[2][2]   0.100000
_atom_sites.fract_transf_matrix[2][3]   0.000000
_atom_sites.fract_transf_matrix[3][1]   0.000000
_atom_sites.fract_transf_matrix[3][2]   0.000000
_atom_sites.fract_transf_matrix[3][3]   0.100000
_atom_sites.fract_transf_vector[1]      0.000000
_atom_sites.fract_transf_vector[2]      0.000000
_atom_sites.fract_transf_vector[3]      0.000000
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
ATOM 1 N  N  . GLY A 1 1 1.000 2.000 3.000 1.00 20.00 1 A 1
ATOM 2 C  CA . GLY A 1 1 2.000 3.000 4.000 1.00 20.00 1 A 1
"""


def _make_pdb():
  return iotbx.pdb.input(source_info="test", lines=_PDB.splitlines())


def _make_cif():
  return iotbx.pdb.mmcif.cif_input(lines=_CIF)


def _fresh_clones(make_input):
  """(label, clone) for a deepcopy and a pickle of a freshly-built input whose
  _scale_matrix cache has been populated.

  A fresh input per clone is deliberate: cif_input builds its builder/hierarchy
  lazily, so calling atoms() on one shared original between the two legs would
  leave them exercising different states (and make the pickle leg incidentally
  depend on the builder being picklable). Each clone here carries only the
  _scale_matrix cache this test targets.
  """
  def cached():
    inp = make_input()
    inp.scale_matrix()                             # populate the cache under test
    return inp
  return [
    ("deepcopy", deepcopy(cached())),
    ("pickle", pickle.loads(pickle.dumps(cached(), pickle.HIGHEST_PROTOCOL))),
  ]


def _assert_roundtrips(make_input):
  """deepcopy and pickle of an input with a populated _scale_matrix cache must
  both succeed (they used to raise) and round-trip the model faithfully."""
  ref = make_input()
  assert "_scale_matrix" not in ref.__dict__, list(ref.__dict__)   # clean before
  scale_matrix = ref.scale_matrix()                                # populates cache
  assert scale_matrix is not None, \
    "fixture must populate the scale-matrix cache"
  assert "_scale_matrix" in ref.__dict__, list(ref.__dict__)       # cached after
  n_atoms = ref.atoms().size()
  crystal_symmetry = ref.crystal_symmetry()
  for label, clone in _fresh_clones(make_input):
    assert clone.atoms().size() == n_atoms, label
    assert clone.crystal_symmetry().is_similar_symmetry(crystal_symmetry), label
    recomputed = clone.scale_matrix()
    assert approx_equal(recomputed[0], scale_matrix[0]), (label, recomputed)
    assert approx_equal(recomputed[1], scale_matrix[1]), (label, recomputed)


def exercise_pdb_input_pickling_after_scale_matrix_cache():
  """A pdb input must stay pickle/deepcopy-able after scale_matrix() caches
  _scale_matrix on it. Without __getstate_manages_dict__ this raised
  "Incomplete pickle support"."""
  _assert_roundtrips(_make_pdb)


def exercise_cif_input_pickling_after_scale_matrix_cache():
  """Forward guard for the planned C++ move of cif_input: the same pickle/
  deepcopy contract must hold once scale_matrix() caches _scale_matrix."""
  _assert_roundtrips(_make_cif)


def exercise_pdb_input_pickling_preserves_suppressed_scale_matrix():
  """xray_structure_simple(crystal_symmetry=...) caches _scale_matrix = None to
  ignore the file's SCALE records -- a decision NOT recoverable from the source
  lines. deepcopy/pickle must preserve that None; a __getstate__ that dropped
  the __dict__ would let the clone silently recompute a non-None matrix (and a
  clone of a malformed-SCALE input would raise where the original,
  short-circuited on the cached None, never did)."""
  # Sanity: the file on its own yields a real matrix; the override suppresses it.
  assert _make_pdb().scale_matrix() is not None
  ref = _make_pdb()
  ref.xray_structure_simple(crystal_symmetry=ref.crystal_symmetry())
  assert ref.scale_matrix() is None, ref.scale_matrix()
  for label, clone in (
      ("deepcopy", deepcopy(ref)),
      ("pickle", pickle.loads(pickle.dumps(ref, pickle.HIGHEST_PROTOCOL)))):
    assert clone.scale_matrix() is None, (label, clone.scale_matrix())


def exercise():
  exercise_pdb_input_pickling_after_scale_matrix_cache()
  exercise_cif_input_pickling_after_scale_matrix_cache()
  exercise_pdb_input_pickling_preserves_suppressed_scale_matrix()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
