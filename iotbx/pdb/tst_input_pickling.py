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

The ``*_lossless`` tests are a regression for mmCIF coordinates being rounded
to PDB precision by ``copy.deepcopy`` / ``pickle`` of a hierarchy (GitHub
issue #1269). ``root.__getstate__`` serialised through PDB text whenever the
model fit the PDB format, which rounds xyz to 1e-3 A, occupancy and B to 1e-2,
U to 1e-4 and drops fp/fdp entirely. ``mmtbx.model.manager`` deepcopies its
``model_input`` to obtain its hierarchy, so every mmCIF model loaded through
it carried PDB-precision atoms.
"""

import pickle
from copy import deepcopy

import iotbx.pdb
import iotbx.pdb.mmcif
from iotbx.pdb import hierarchy
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

# Two atoms whose values PDB text cannot represent: 5-decimal coordinates,
# 3-decimal occupancy, 4-decimal B, 5-decimal anisotropic U. Every one of
# them is deposited like this in the PDB archive today.
_CIF_PRECISE = """\
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
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.type_symbol
_atom_site_anisotrop.pdbx_label_atom_id
_atom_site_anisotrop.pdbx_label_alt_id
_atom_site_anisotrop.pdbx_label_comp_id
_atom_site_anisotrop.pdbx_label_asym_id
_atom_site_anisotrop.pdbx_label_seq_id
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
1 N N . GLY A 1 0.12345 0.23456 0.34567 0.01234 -0.02345 0.03456
"""


def _make_precise_cif():
  return iotbx.pdb.mmcif.cif_input(lines=_CIF_PRECISE)


def _precise_hierarchy():
  """Hierarchy of _CIF_PRECISE with the fields that have no mmCIF _atom_site
  home (fp, fdp, resolution, the sigmas) set in code on the first atom, so
  that every numeric field holds a value the text would round or drop, after
  checking that the fixture really carries values PDB text cannot hold."""
  h = _make_precise_cif().construct_hierarchy()
  atoms = h.atoms()
  atoms[0].fp = 0.12345
  atoms[0].fdp = 0.6789
  atoms[0].resolution = 2.567
  atoms[0].sigxyz = (0.011, 0.022, 0.033)
  atoms[0].sigocc = 0.044
  atoms[0].sigb = 0.055
  assert atoms[0].xyz == (1.12345, 2.23456, 3.34567), atoms[0].xyz
  assert atoms[0].occ == 0.987, atoms[0].occ
  assert atoms[0].b == 20.1234, atoms[0].b
  assert atoms[0].uij_is_defined() and atoms[0].uij[0] == 0.12345, atoms[0].uij
  assert not atoms[1].uij_is_defined()
  return h


def _atom_state(h):
  """Everything per atom that a copy must reproduce: labels and the exact
  numeric values (compared with ==, no tolerance)."""
  return [
    (a.id_str(), a.xyz, a.sigxyz, a.occ, a.sigocc, a.b, a.sigb, a.uij,
     a.fp, a.fdp, a.resolution)
    for a in h.atoms()]


def _assert_same_hierarchy(ref, clone, label):
  assert clone.is_similar_hierarchy(ref), label
  assert list(clone.info) == list(ref.info), label
  assert _atom_state(clone) == _atom_state(ref), (label, _atom_state(clone))


def exercise_hierarchy_deepcopy_is_lossless():
  """copy.deepcopy(root) must reproduce every numeric value exactly instead of
  rounding to PDB text precision."""
  ref = _precise_hierarchy()
  _assert_same_hierarchy(ref, deepcopy(ref), "deepcopy")


def exercise_hierarchy_pickle_is_lossless():
  """pickle of root must reproduce every numeric value exactly, at every
  protocol in use (tst_hierarchy pickles with protocol 1)."""
  ref = _precise_hierarchy()
  for protocol in (1, pickle.HIGHEST_PROTOCOL):
    clone = pickle.loads(pickle.dumps(ref, protocol))
    _assert_same_hierarchy(ref, clone, "pickle protocol %d" % protocol)


def exercise_hierarchy_pickle_is_lossless_outside_pdb_format():
  """The mmCIF text branch of __getstate__ (model does not fit PDB format)
  must be exact too, not merely 5-decimal."""
  ref = _precise_hierarchy()
  ref.only_model().only_chain().id = "LONGID"
  assert not ref.fits_in_pdb_format()
  _assert_same_hierarchy(ref, pickle.loads(pickle.dumps(ref)), "pickle")


def exercise_cif_input_deepcopy_keeps_hierarchy_precision():
  """deepcopy(cif_input).construct_hierarchy() is how mmtbx.model.manager gets
  its hierarchy. With the hierarchy already cached on the input, the copy went
  through root.__getstate__ and lost precision."""
  inp = _make_precise_cif()
  ref = inp.construct_hierarchy()                 # populate the cache
  clone = deepcopy(inp).construct_hierarchy()
  _assert_same_hierarchy(ref, clone, "deepcopy(cif_input)")
  # The copy must be independent of the original's cached hierarchy.
  clone.atoms()[0].set_xyz((9, 9, 9))
  assert ref.atoms()[0].xyz == (1.12345, 2.23456, 3.34567)


def exercise_hierarchy_copy_refuses_stray_attributes():
  """Existing contract (used to make forward-compatible hierarchies
  unpicklable): a root carrying a Python attribute other than the known
  caches refuses to pickle. The same must hold for copy.deepcopy, which is now
  a separate code path."""
  h = _precise_hierarchy()
  h._not_a_known_cache = 1
  for label, copier in (
      ("pickle", lambda x: pickle.dumps(x)),
      ("deepcopy", lambda x: deepcopy(x))):
    try:
      copier(h)
    except AssertionError:
      pass
    else:
      raise AssertionError(
        "%s of a root with a stray attribute succeeded" % label)


def exercise_hierarchy_setstate_accepts_older_pickle_versions():
  """Pickles written before the exact numeric state was added (versions 1 and
  2) must still load; they carry text only."""
  ref = _make_pdb().construct_hierarchy()
  pdb_text = ref.as_pdb_string()
  for state in (
      (1, ref.info, pdb_text),
      (2, hierarchy.pickle_import_trigger(), ref.info, pdb_text)):
    clone = hierarchy.root()
    clone.__setstate__(state)
    assert clone.is_similar_hierarchy(ref), state[0]
    assert clone.as_pdb_string() == pdb_text, state[0]


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
  exercise_hierarchy_deepcopy_is_lossless()
  exercise_hierarchy_pickle_is_lossless()
  exercise_hierarchy_pickle_is_lossless_outside_pdb_format()
  exercise_cif_input_deepcopy_keeps_hierarchy_precision()
  exercise_hierarchy_copy_refuses_stray_attributes()
  exercise_hierarchy_setstate_accepts_older_pickle_versions()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
