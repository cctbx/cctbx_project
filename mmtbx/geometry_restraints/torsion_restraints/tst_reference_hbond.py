from __future__ import absolute_import, division, print_function

import sys

from cctbx import geometry_restraints
from cctbx.array_family import flex
from cctbx.geometry_restraints.linking_class import linking_class
from libtbx.utils import null_out
from scitbx.matrix import col
from six.moves import cStringIO as StringIO
import iotbx.pdb
import mmtbx.model

from mmtbx.geometry_restraints.torsion_restraints import reference_model as rm_mod
from mmtbx.geometry_restraints.torsion_restraints.reference_model import (
  _detect_reference_hbonds,
  _ensure_hydrogens_on_reference,
  reference_model,
  reference_model_params,
)


def exercise_origin_id_registered():
  """'reference hydrogen bonds' must be registered as a distinct origin_id
  (different from the existing 'hydrogen bonds' origin used by SS restraints)
  so that bond/angle proxies derived from reference-model H-bonds can be
  tagged and counted independently."""
  lc = linking_class()
  oid = lc.get_origin_id('reference hydrogen bonds')
  assert isinstance(oid, int), \
    "expected int origin_id, got %r" % (oid,)
  assert oid > 0, "origin_id must be > 0 (0 is reserved for 'covalent geometry')"
  # distinct from the existing SS 'hydrogen bonds' origin_id
  ss_oid = lc.get_origin_id('hydrogen bonds')
  assert oid != ss_oid, \
    "new origin_id must differ from SS 'hydrogen bonds' (%d)" % ss_oid


def exercise_phil_defaults():
  """The new hydrogen_bonds { ... } sub-block inside reference_model parses
  with the documented defaults."""
  p = reference_model_params.extract().reference_model
  assert hasattr(p, 'hydrogen_bonds'), \
    "reference_model PHIL must contain hydrogen_bonds sub-block"
  hb = p.hydrogen_bonds
  assert hb.enabled is False, hb.enabled
  assert hb.target == 'ideal', hb.target
  assert hb.restrain_angles is False, hb.restrain_angles
  assert abs(hb.sigma_bond - 0.05) < 1.e-9, hb.sigma_bond
  assert abs(hb.sigma_angle - 5.0) < 1.e-9, hb.sigma_angle
  assert abs(hb.slack_bond - 0.0) < 1.e-9, hb.slack_bond
  assert abs(hb.ideal_distance - 2.9) < 1.e-9, hb.ideal_distance
  assert abs(hb.ideal_angle - 180.0) < 1.e-9, hb.ideal_angle
  assert hb.add_hydrogens_if_missing is True, hb.add_hydrogens_if_missing
  assert abs(hb.partner_distance_cutoff - 5.0) < 1.e-9, hb.partner_distance_cutoff


# A tiny single-residue PDB used to instantiate a real GRM cheaply.
_minimal_pdb_str = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.458  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      11.984  10.717  11.232  1.00 20.00           C
ATOM      4  O   ALA A   1      11.247  10.969  12.184  1.00 20.00           O
ATOM      5  CB  ALA A   1      11.954  10.690   8.737  1.00 20.00           C
END
"""


def _build_minimal_model_and_geometry():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=_minimal_pdb_str)
  model = mmtbx.model.manager(model_input=pdb_inp)
  model.process(make_restraints=True)
  return model, model.get_restraints_manager().geometry


def exercise_get_hbond_proxies_returns_empty_when_disabled():
  """reference_model.get_hbond_proxies returns ([], []) when the
  hydrogen_bonds.enabled flag is False, without doing any work."""
  model, geometry = _build_minimal_model_and_geometry()
  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.hydrogen_bonds.enabled = False
  rm = reference_model(
    model=model,
    reference_hierarchy_list=[model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  bp, ap = rm.get_hbond_proxies(
    geometry=geometry, sites_cart=model.get_sites_cart())
  assert bp == [] and ap == [], (bp, ap)
  # GRM count method is still reachable and reports zero
  assert geometry.get_n_reference_hbond_proxies() == 0


# Two-residue fragment with an explicit H atom on N of the second residue.
_pdb_with_h = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.458  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      11.984  10.717  11.232  1.00 20.00           C
ATOM      4  O   ALA A   1      11.247  10.969  12.184  1.00 20.00           O
ATOM      5  CB  ALA A   1      11.954  10.690   8.737  1.00 20.00           C
ATOM      6  N   ALA A   2      13.241  11.014  11.218  1.00 20.00           N
ATOM      7  H   ALA A   2      13.700  10.823  10.500  1.00 20.00           H
ATOM      8  CA  ALA A   2      13.876  11.703  12.330  1.00 20.00           C
ATOM      9  C   ALA A   2      13.347  13.122  12.483  1.00 20.00           C
ATOM     10  O   ALA A   2      12.158  13.375  12.371  1.00 20.00           O
ATOM     11  CB  ALA A   2      15.385  11.748  12.140  1.00 20.00           C
END
"""


def exercise_ensure_hydrogens_idempotent_on_h_bearing():
  """_ensure_hydrogens_on_reference returns the input hierarchy unchanged
  when it already passes mmtbx.model.manager.has_hd()."""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=_pdb_with_h)
  ref_h = pdb_inp.construct_hierarchy()
  n_atoms_before = ref_h.atoms_size()
  result = _ensure_hydrogens_on_reference(ref_h)
  assert result.atoms_size() == n_atoms_before, \
    "expected unchanged hierarchy, got atoms_size %d -> %d" % (
      n_atoms_before, result.atoms_size())
  # Idempotent on the object (return the same hierarchy without copying)
  assert result is ref_h, "expected same hierarchy instance"


# Two-residue fragment WITHOUT explicit hydrogens. _ensure_hydrogens_on_reference
# should add them.
_pdb_without_h = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.458  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      11.984  10.717  11.232  1.00 20.00           C
ATOM      4  O   ALA A   1      11.247  10.969  12.184  1.00 20.00           O
ATOM      5  CB  ALA A   1      11.954  10.690   8.737  1.00 20.00           C
ATOM      6  N   ALA A   2      13.241  11.014  11.218  1.00 20.00           N
ATOM      7  CA  ALA A   2      13.876  11.703  12.330  1.00 20.00           C
ATOM      8  C   ALA A   2      13.347  13.122  12.483  1.00 20.00           C
ATOM      9  O   ALA A   2      12.158  13.375  12.371  1.00 20.00           O
ATOM     10  CB  ALA A   2      15.385  11.748  12.140  1.00 20.00           C
END
"""


# 9-residue alpha-helical poly-Gly with explicit H atoms. Borrowed from
# mmtbx/nci/tst_hbond.py (pdb_str_00) where it is the standard alpha-helix
# fixture for H-bond detection tests.
_pdb_helix_with_h = """
ATOM      1  N   GLY A   1      -5.606  -2.251 -12.878  1.00  0.00           N
ATOM      2  CA  GLY A   1      -5.850  -1.194 -13.852  1.00  0.00           C
ATOM      3  C   GLY A   1      -5.186  -1.524 -15.184  1.00  0.00           C
ATOM      4  O   GLY A   1      -5.744  -1.260 -16.249  1.00  0.00           O
ATOM      0  H1  GLY A   1      -6.104  -2.109 -12.154  1.00  0.00           H   new
ATOM      0  H2  GLY A   1      -5.819  -3.038 -13.234  1.00  0.00           H   new
ATOM      0  H3  GLY A   1      -4.746  -2.252 -12.650  1.00  0.00           H   new
ATOM      0  HA2 GLY A   1      -6.805  -1.081 -13.981  1.00  0.00           H   new
ATOM      0  HA3 GLY A   1      -5.507  -0.352 -13.515  1.00  0.00           H   new
ATOM      1  N   GLY A   2      -3.992  -2.102 -15.115  1.00  0.00           N
ATOM      2  CA  GLY A   2      -3.261  -2.499 -16.313  1.00  0.00           C
ATOM      3  C   GLY A   2      -3.961  -3.660 -17.011  1.00  0.00           C
ATOM      4  O   GLY A   2      -4.016  -3.716 -18.240  1.00  0.00           O
ATOM      0  H   GLY A   2      -3.585  -2.274 -14.378  1.00  0.00           H   new
ATOM      0  HA2 GLY A   2      -3.191  -1.745 -16.920  1.00  0.00           H   new
ATOM      0  HA3 GLY A   2      -2.356  -2.756 -16.075  1.00  0.00           H   new
ATOM      1  N   GLY A   3      -4.492  -4.585 -16.219  1.00  0.00           N
ATOM      2  CA  GLY A   3      -5.216  -5.731 -16.755  1.00  0.00           C
ATOM      3  C   GLY A   3      -6.531  -5.289 -17.389  1.00  0.00           C
ATOM      4  O   GLY A   3      -6.939  -5.814 -18.425  1.00  0.00           O
ATOM      0  H   GLY A   3      -4.443  -4.566 -15.361  1.00  0.00           H   new
ATOM      0  HA2 GLY A   3      -4.669  -6.185 -17.416  1.00  0.00           H   new
ATOM      0  HA3 GLY A   3      -5.392  -6.369 -16.046  1.00  0.00           H   new
ATOM      1  N   GLY A   4      -7.189  -4.323 -16.758  1.00  0.00           N
ATOM      2  CA  GLY A   4      -8.442  -3.785 -17.273  1.00  0.00           C
ATOM      3  C   GLY A   4      -8.205  -3.003 -18.561  1.00  0.00           C
ATOM      4  O   GLY A   4      -9.007  -3.065 -19.492  1.00  0.00           O
ATOM      0  H   GLY A   4      -6.924  -3.963 -16.024  1.00  0.00           H   new
ATOM      0  HA2 GLY A   4      -9.066  -4.509 -17.439  1.00  0.00           H   new
ATOM      0  HA3 GLY A   4      -8.848  -3.207 -16.608  1.00  0.00           H   new
ATOM      1  N   GLY A   5      -7.099  -2.269 -18.604  1.00  0.00           N
ATOM      2  CA  GLY A   5      -6.735  -1.498 -19.787  1.00  0.00           C
ATOM      3  C   GLY A   5      -6.358  -2.423 -20.939  1.00  0.00           C
ATOM      4  O   GLY A   5      -6.687  -2.157 -22.094  1.00  0.00           O
ATOM      0  H   GLY A   5      -6.541  -2.204 -17.953  1.00  0.00           H   new
ATOM      0  HA2 GLY A   5      -7.477  -0.932 -20.051  1.00  0.00           H   new
ATOM      0  HA3 GLY A   5      -5.990  -0.912 -19.580  1.00  0.00           H   new
ATOM      1  N   GLY A   6      -5.665  -3.509 -20.614  1.00  0.00           N
ATOM      2  CA  GLY A   6      -5.268  -4.493 -21.614  1.00  0.00           C
ATOM      3  C   GLY A   6      -6.485  -5.236 -22.153  1.00  0.00           C
ATOM      4  O   GLY A   6      -6.565  -5.533 -23.345  1.00  0.00           O
ATOM      0  H   GLY A   6      -5.413  -3.695 -19.813  1.00  0.00           H   new
ATOM      0  HA2 GLY A   6      -4.804  -4.051 -22.343  1.00  0.00           H   new
ATOM      0  HA3 GLY A   6      -4.645  -5.125 -21.223  1.00  0.00           H   new
ATOM      1  N   GLY A   7      -7.430  -5.532 -21.267  1.00  0.00           N
ATOM      2  CA  GLY A   7      -8.660  -6.212 -21.655  1.00  0.00           C
ATOM      3  C   GLY A   7      -9.529  -5.303 -22.518  1.00  0.00           C
ATOM      4  O   GLY A   7     -10.158  -5.756 -23.474  1.00  0.00           O
ATOM      0  H   GLY A   7      -7.377  -5.346 -20.429  1.00  0.00           H   new
ATOM      0  HA2 GLY A   7      -8.447  -7.022 -22.144  1.00  0.00           H   new
ATOM      0  HA3 GLY A   7      -9.151  -6.479 -20.862  1.00  0.00           H   new
ATOM      1  N   GLY A   8      -9.559  -4.021 -22.172  1.00  0.00           N
ATOM      2  CA  GLY A   8     -10.324  -3.039 -22.930  1.00  0.00           C
ATOM      3  C   GLY A   8      -9.706  -2.819 -24.306  1.00  0.00           C
ATOM      4  O   GLY A   8     -10.416  -2.660 -25.299  1.00  0.00           O
ATOM      0  H   GLY A   8      -9.139  -3.697 -21.495  1.00  0.00           H   new
ATOM      0  HA2 GLY A   8     -11.241  -3.341 -23.027  1.00  0.00           H   new
ATOM      0  HA3 GLY A   8     -10.352  -2.199 -22.445  1.00  0.00           H   new
ATOM      1  N   GLY A   9      -8.378  -2.810 -24.356  1.00  0.00           N
ATOM      2  CA  GLY A   9      -7.658  -2.641 -25.613  1.00  0.00           C
ATOM      3  C   GLY A   9      -7.843  -3.861 -26.508  1.00  0.00           C
ATOM      4  O   GLY A   9      -7.980  -3.734 -27.725  1.00  0.00           O
ATOM      0  H   GLY A   9      -7.872  -2.901 -23.667  1.00  0.00           H   new
ATOM      0  HA2 GLY A   9      -7.978  -1.847 -26.070  1.00  0.00           H   new
ATOM      0  HA3 GLY A   9      -6.714  -2.505 -25.435  1.00  0.00           H   new
"""


def exercise_detect_reference_hbonds():
  """_detect_reference_hbonds returns donor/H/acceptor i_seqs and geometry
  values for each H-bond found in the reference hierarchy."""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=_pdb_helix_with_h)
  ref_h = pdb_inp.construct_hierarchy()
  hbonds = _detect_reference_hbonds(ref_h, log=null_out())
  assert len(hbonds) > 0, "expected at least one H-bond in poly-Gly helix"
  atoms = ref_h.atoms()
  for d_iseq, h_iseq, a_iseq, d_DA, d_HA, a_DHA in hbonds:
    # donor heavy atom is N, H atom is H/H1/H2/H3, acceptor is O (alpha helix)
    d_elt = atoms[d_iseq].element.strip()
    h_elt = atoms[h_iseq].element.strip()
    a_elt = atoms[a_iseq].element.strip()
    assert d_elt == 'N', (d_iseq, d_elt)
    assert h_elt == 'H' or h_elt == 'D', (h_iseq, h_elt)
    assert a_elt == 'O', (a_iseq, a_elt)
    assert 2.5 < d_DA < 3.5, d_DA
    assert 1.5 < d_HA < 2.8, d_HA
    assert 120.0 < a_DHA <= 180.0, a_DHA


def _build_helix_model():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=_pdb_helix_with_h)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  return model


def exercise_end_to_end_as_found_no_angles_no_ss():
  """End-to-end with self-reference, target=as_found, no SS overlap, no
  angle restraints. Expects bond proxies tagged with the new origin_id and
  distance_ideal matching the measured D-A distance in the reference."""
  model = _build_helix_model()
  geometry = model.get_restraints_manager().geometry
  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.hydrogen_bonds.enabled = True
  p.reference_model.hydrogen_bonds.target = 'as_found'
  p.reference_model.hydrogen_bonds.restrain_angles = False
  rm = reference_model(
    model=model,
    reference_hierarchy_list=[model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)
  bp, ap = rm.get_hbond_proxies(
    geometry=geometry, sites_cart=model.get_sites_cart())
  assert len(bp) > 0, "expected at least one bond proxy from 9-res helix"
  assert ap == [], "restrain_angles=False should yield no angle proxies"
  ref_hb_oid = linking_class().get_origin_id('reference hydrogen bonds')
  for proxy in bp:
    assert proxy.origin_id == ref_hb_oid, proxy.origin_id
    assert 2.5 < proxy.distance_ideal < 3.5, proxy.distance_ideal


def exercise_end_to_end_target_ideal():
  """With target='ideal' all emitted proxies use ideal_distance regardless
  of the geometry measured in the reference."""
  model = _build_helix_model()
  geometry = model.get_restraints_manager().geometry
  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.hydrogen_bonds.enabled = True
  p.reference_model.hydrogen_bonds.target = 'ideal'
  p.reference_model.hydrogen_bonds.ideal_distance = 2.85
  rm = reference_model(
    model=model,
    reference_hierarchy_list=[model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)
  bp, _ = rm.get_hbond_proxies(
    geometry=geometry, sites_cart=model.get_sites_cart())
  assert len(bp) > 0
  ref_hb_oid = linking_class().get_origin_id('reference hydrogen bonds')
  for proxy in bp:
    assert proxy.origin_id == ref_hb_oid
    assert abs(proxy.distance_ideal - 2.85) < 1.e-6, proxy.distance_ideal


def exercise_function_rename_and_dual_dispatch():
  """add_reference_model_restraints_if_requested (formerly
  add_reference_dihedral_restraints_if_requested) dispatches independently
  on params.enabled (torsion master) and params.hydrogen_bonds.enabled.
  Either, both, or neither may be active."""
  # Old name must be gone, new name must exist
  assert hasattr(rm_mod, 'add_reference_model_restraints_if_requested')
  assert not hasattr(rm_mod, 'add_reference_dihedral_restraints_if_requested')
  fn = rm_mod.add_reference_model_restraints_if_requested
  # Case A: both flags enabled
  model = _build_helix_model()
  geometry = model.get_restraints_manager().geometry
  p = rm_mod.reference_model_params.extract()
  p.reference_model.enabled = True
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.hydrogen_bonds.enabled = True
  fn(model=model, geometry=geometry,
     params=p.reference_model, selection=None, log=null_out())
  assert geometry.reference_dihedral_manager is not None
  assert geometry.get_n_reference_hbond_proxies() > 0, \
    "expected reference H-bond proxies attached"
  # Case B: only hydrogen_bonds enabled, master torsion flag off
  model = _build_helix_model()
  geometry = model.get_restraints_manager().geometry
  p = rm_mod.reference_model_params.extract()
  p.reference_model.enabled = False
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.hydrogen_bonds.enabled = True
  fn(model=model, geometry=geometry,
     params=p.reference_model, selection=None, log=null_out())
  # reference_dihedral_manager IS attached (needed by GRM hbond method) but
  # no dihedral proxies are applied (master torsion flag off would gate the
  # actual refinement-time application elsewhere; here we just check that
  # h-bond proxies were added).
  assert geometry.reference_dihedral_manager is not None
  assert geometry.get_n_reference_hbond_proxies() > 0
  # Case C: neither flag -> no-op
  model = _build_helix_model()
  geometry = model.get_restraints_manager().geometry
  p = rm_mod.reference_model_params.extract()
  p.reference_model.enabled = False
  p.reference_model.hydrogen_bonds.enabled = False
  fn(model=model, geometry=geometry,
     params=p.reference_model, selection=None, log=null_out())
  assert geometry.reference_dihedral_manager is None
  assert geometry.get_n_reference_hbond_proxies() == 0


def exercise_geo_output_labels_new_origin_id():
  """The .geo output contains the new origin_id headers ('Reference H-bond'
  for bonds, 'Reference H-bond angles' for angles) registered in
  auto_linking_types.py."""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=_pdb_helix_with_h)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.reference_model.use_starting_model_as_reference = True
  params.reference_model.hydrogen_bonds.enabled = True
  params.reference_model.hydrogen_bonds.restrain_angles = True
  model.process(pdb_interpretation_params=params, make_restraints=True)
  geometry = model.get_restraints_manager().geometry
  buf = StringIO()
  geometry.write_geo_file(
    sites_cart=model.get_sites_cart(),
    site_labels=[a.id_str() for a in model.get_hierarchy().atoms()],
    file_descriptor=buf)
  geo_text = buf.getvalue()
  assert 'Reference H-bond' in geo_text, \
    "expected 'Reference H-bond' header in .geo output"


def exercise_integration_via_model_process():
  """Full path through mmtbx.model.manager.process(). Setting
  reference_model.hydrogen_bonds.enabled in pdb_interpretation_params should
  produce reference H-bond proxies in the final geometry."""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=_pdb_helix_with_h)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.reference_model.use_starting_model_as_reference = True
  params.reference_model.hydrogen_bonds.enabled = True
  model.process(pdb_interpretation_params=params, make_restraints=True)
  geometry = model.get_restraints_manager().geometry
  n_ref_hb = geometry.get_n_reference_hbond_proxies()
  assert n_ref_hb > 0, \
    "expected reference H-bond proxies via model.process(), got %d" % n_ref_hb


def _strip_h_from_pdb_lines(pdb_str):
  """Helper: return pdb_str with H/D ATOM lines removed."""
  keep = []
  for line in pdb_str.split("\n"):
    if line.startswith("ATOM") and line[76:78].strip() in ('H', 'D'):
      continue
    keep.append(line)
  return "\n".join(keep)


def _build_two_chain_working_and_one_chain_ref():
  """Build a working model with two helical chains (A and B, NCS-related by
  a 30-A translation) and a reference hierarchy containing just chain A.

  Returns (working_model, ref_hierarchy).
  """
  # Parse the single-chain helix
  pdb_inp_a = iotbx.pdb.input(source_info=None, lines=_pdb_helix_with_h)
  hier_a = pdb_inp_a.construct_hierarchy()
  # Make chain B = chain A translated by (30, 0, 0)
  hier_b = hier_a.deep_copy()
  hier_b.only_model().chains()[0].id = 'B'
  xyz_b = hier_b.atoms().extract_xyz()
  shift = col((30.0, 0.0, 0.0))
  xyz_b_shifted = flex.vec3_double([
    tuple(col(c) + shift) for c in xyz_b])
  hier_b.atoms().set_xyz(xyz_b_shifted)
  # Concatenate A and B into one model object
  hier_ab = hier_a.deep_copy()
  for ch in hier_b.only_model().chains():
    hier_ab.only_model().append_chain(ch.detached_copy())
  hier_ab.atoms().reset_i_seq()
  pdb_str_ab = hier_ab.as_pdb_string()
  pdb_inp_ab = iotbx.pdb.input(source_info=None, lines=pdb_str_ab)
  work_model = mmtbx.model.manager(model_input=pdb_inp_ab, log=null_out())
  work_model.process(make_restraints=True)
  return work_model, hier_a


def exercise_ncs_intra_chain_expansion():
  """Working has two NCS-related chains, reference has one chain; one
  intra-chain reference H-bond produces two bond proxies in the working
  model (one per NCS copy)."""
  work_model, ref_hier = _build_two_chain_working_and_one_chain_ref()
  geometry = work_model.get_restraints_manager().geometry
  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = False
  p.reference_model.hydrogen_bonds.enabled = True
  rm = reference_model(
    model=work_model,
    reference_hierarchy_list=[ref_hier],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)
  bp, _ = rm.get_hbond_proxies(
    geometry=geometry, sites_cart=work_model.get_sites_cart())
  # The reference produces some number of H-bonds; each should appear in
  # both working chains -> total = 2 * (# unique ref H-bonds).
  assert len(bp) > 0
  assert len(bp) % 2 == 0, \
    "expected even number of proxies (2 per ref H-bond), got %d" % len(bp)
  # Each working pair should be within one chain (intra-chain restraint).
  atoms = work_model.get_hierarchy().atoms()
  chains_used = set()
  for proxy in bp:
    i, j = proxy.i_seqs
    chain_i = atoms[i].parent().parent().parent().id
    chain_j = atoms[j].parent().parent().parent().id
    assert chain_i == chain_j, \
      "expected intra-chain pair, got %s-%s" % (chain_i, chain_j)
    chains_used.add(chain_i)
  assert chains_used == {'A', 'B'}, chains_used
  ref_hb_oid = linking_class().get_origin_id('reference hydrogen bonds')
  for proxy in bp:
    assert proxy.origin_id == ref_hb_oid


def exercise_partner_distance_cutoff_rejects():
  """A tight partner_distance_cutoff rejects all candidates (zero proxies).
  Under NCS, the same filter rejects the K^2 - K cross-copy candidates and
  keeps only the K within-copy proxies."""
  # Case A: extremely tight cutoff -> no proxies
  model = _build_helix_model()
  geometry = model.get_restraints_manager().geometry
  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.hydrogen_bonds.enabled = True
  p.reference_model.hydrogen_bonds.partner_distance_cutoff = 0.1
  rm = reference_model(
    model=model,
    reference_hierarchy_list=[model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)
  bp, _ = rm.get_hbond_proxies(
    geometry=geometry, sites_cart=model.get_sites_cart())
  assert len(bp) == 0, \
    "expected zero proxies with cutoff=0.1, got %d" % len(bp)

  # Case B: NCS - cutoff just above H-bond distance (4.0 A) keeps within-copy
  # pairs (~2.9 A). Cross-copy pairs (separated by NCS translation 30 A) are
  # well above 4 A and must be rejected.
  work_model, ref_hier = _build_two_chain_working_and_one_chain_ref()
  geometry = work_model.get_restraints_manager().geometry
  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = False
  p.reference_model.hydrogen_bonds.enabled = True
  p.reference_model.hydrogen_bonds.partner_distance_cutoff = 4.0
  rm = reference_model(
    model=work_model,
    reference_hierarchy_list=[ref_hier],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)
  bp, _ = rm.get_hbond_proxies(
    geometry=geometry, sites_cart=work_model.get_sites_cart())
  atoms = work_model.get_hierarchy().atoms()
  for proxy in bp:
    i, j = proxy.i_seqs
    assert atoms[i].parent().parent().parent().id == \
           atoms[j].parent().parent().parent().id, \
      "cross-NCS pair leaked through cutoff filter"


def exercise_warning_when_angles_requested_but_no_working_h():
  """When restrain_angles=True is requested but the working model has no
  H atoms, add_reference_model_restraints_if_requested prints a clear
  warning next to the bond/angle count line."""
  import os, tempfile
  # Working: H-less helix
  working_pdb = _strip_h_from_pdb_lines(_pdb_helix_with_h)
  work_inp = iotbx.pdb.input(source_info=None, lines=working_pdb)
  work_model = mmtbx.model.manager(model_input=work_inp, log=null_out())
  work_model.process(make_restraints=True)
  geometry = work_model.get_restraints_manager().geometry
  # Reference: H-bearing helix written to a temp file (so the non-self path
  # in add_reference_model_restraints_if_requested is exercised; self-ref
  # with H-augmentation does not currently support an H-less working model)
  tmpdir = tempfile.mkdtemp(prefix='tst_ref_hbond_warning_')
  ref_path = os.path.join(tmpdir, 'ref.pdb')
  with open(ref_path, 'w') as f:
    f.write(_pdb_helix_with_h)
  try:
    p = reference_model_params.extract()
    p.reference_model.use_starting_model_as_reference = False
    p.reference_model.file = [ref_path]
    p.reference_model.hydrogen_bonds.enabled = True
    p.reference_model.hydrogen_bonds.restrain_angles = True
    buf = StringIO()
    rm_mod.add_reference_model_restraints_if_requested(
      model=work_model, geometry=geometry,
      params=p.reference_model, selection=None, log=buf)
    msg = buf.getvalue()
    # Collapse whitespace so wrapped multi-line print is matched intact.
    flat = ' '.join(msg.split())
    expected = (
      'WARNING: restrain_angles=True was requested but the working '
      'model has no hydrogen atoms; all D-H-A angle restraints '
      'were skipped.')
    assert expected in flat, msg
  finally:
    os.remove(ref_path)
    os.rmdir(tmpdir)


def exercise_h_less_working_skips_angles_keeps_bonds():
  """When restrain_angles=True but the working model has no H, the angle
  proxy is silently skipped for each H-bond while the bond proxy is still
  emitted. No exception."""
  # Working: H-less version of the 9-residue helix
  working_pdb = _strip_h_from_pdb_lines(_pdb_helix_with_h)
  work_inp = iotbx.pdb.input(source_info=None, lines=working_pdb)
  work_model = mmtbx.model.manager(model_input=work_inp, log=null_out())
  work_model.process(make_restraints=True)
  geometry = work_model.get_restraints_manager().geometry
  # Reference: H-bearing version of the same helix
  ref_inp = iotbx.pdb.input(source_info=None, lines=_pdb_helix_with_h)
  ref_h = ref_inp.construct_hierarchy()
  # Configure: H-bonds enabled, angles requested, non-self-reference
  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = False
  p.reference_model.hydrogen_bonds.enabled = True
  p.reference_model.hydrogen_bonds.restrain_angles = True
  rm = reference_model(
    model=work_model,
    reference_hierarchy_list=[ref_h],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)
  bp, ap = rm.get_hbond_proxies(
    geometry=geometry, sites_cart=work_model.get_sites_cart())
  assert len(bp) > 0, "expected bond proxies even when working has no H"
  assert len(ap) == 0, \
    "expected no angle proxies when working model lacks H, got %d" % len(ap)
  ref_hb_oid = linking_class().get_origin_id('reference hydrogen bonds')
  for proxy in bp:
    assert proxy.origin_id == ref_hb_oid


def exercise_angle_proxies_from_restrain_angles():
  """With restrain_angles=True, D-H-A angle proxies are emitted alongside
  bond proxies, tagged with the new origin_id and with angle_ideal close
  to the measured a_DHA when target='as_found'."""
  model = _build_helix_model()
  geometry = model.get_restraints_manager().geometry
  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.hydrogen_bonds.enabled = True
  p.reference_model.hydrogen_bonds.target = 'as_found'
  p.reference_model.hydrogen_bonds.restrain_angles = True
  rm = reference_model(
    model=model,
    reference_hierarchy_list=[model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)
  bp, ap = rm.get_hbond_proxies(
    geometry=geometry, sites_cart=model.get_sites_cart())
  assert len(bp) > 0
  assert len(ap) > 0
  ref_hb_oid = linking_class().get_origin_id('reference hydrogen bonds')
  atoms = model.get_hierarchy().atoms()
  for proxy in ap:
    assert proxy.origin_id == ref_hb_oid
    # i_seqs is [donor, h, acceptor]
    d_iseq, h_iseq, a_iseq = proxy.i_seqs
    assert atoms[d_iseq].element.strip() == 'N', atoms[d_iseq].name
    assert atoms[h_iseq].element.strip() in ('H', 'D'), atoms[h_iseq].name
    assert atoms[a_iseq].element.strip() == 'O', atoms[a_iseq].name
    assert 120.0 < proxy.angle_ideal <= 180.0, proxy.angle_ideal


def exercise_ss_filter_excludes_overlap():
  """Pairs already restrained as SS hydrogen bonds are skipped by
  reference_model.get_hbond_proxies (no double-restraining)."""
  # 1st pass: discover the i_seqs of one reference H-bond
  model = _build_helix_model()
  geometry = model.get_restraints_manager().geometry
  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.hydrogen_bonds.enabled = True
  rm = reference_model(
    model=model,
    reference_hierarchy_list=[model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)
  bp, _ = rm.get_hbond_proxies(
    geometry=geometry, sites_cart=model.get_sites_cart())
  assert len(bp) > 0
  baseline_pairs = set(frozenset(proxy.i_seqs) for proxy in bp)
  victim_pair = next(iter(baseline_pairs))
  i, j = tuple(victim_pair)

  # 2nd pass: rebuild and inject SS H-bond proxy on the victim pair
  model = _build_helix_model()
  geometry = model.get_restraints_manager().geometry
  ss_oid = linking_class().get_origin_id('hydrogen bonds')
  ss_proxy = geometry_restraints.bond_simple_proxy(
    i_seqs=[i, j],
    distance_ideal=2.9,
    weight=400.0,
    origin_id=ss_oid)
  geometry.add_new_bond_restraints_in_place([ss_proxy],
    sites_cart=model.get_sites_cart())
  rm = reference_model(
    model=model,
    reference_hierarchy_list=[model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)
  bp2, _ = rm.get_hbond_proxies(
    geometry=geometry, sites_cart=model.get_sites_cart())
  filtered_pairs = set(frozenset(proxy.i_seqs) for proxy in bp2)
  assert victim_pair not in filtered_pairs, \
    "victim pair %s should have been filtered out" % (victim_pair,)
  # the other reference H-bonds should still be there
  assert len(filtered_pairs) == len(baseline_pairs) - 1, \
    (len(filtered_pairs), len(baseline_pairs))


def exercise_collect_ss_hbond_pairs():
  """GRM.get_hbond_proxies_iseqs returns pairs for proxies with origin_id
  'hydrogen bonds' (the SS origin). This is the basis for the SS overlap
  filter used by the new reference-H-bond builder."""
  model, geometry = _build_minimal_model_and_geometry()
  ss_oid = linking_class().get_origin_id('hydrogen bonds')
  proxy = geometry_restraints.bond_simple_proxy(
    i_seqs=[0, 3],
    distance_ideal=3.0,
    weight=400.0,
    origin_id=ss_oid)
  sites_cart = model.get_sites_cart()
  geometry.add_new_bond_restraints_in_place([proxy], sites_cart=sites_cart)
  pairs = geometry.get_hbond_proxies_iseqs()
  pair_set = set(frozenset(p) for p in pairs)
  assert frozenset({0, 3}) in pair_set, pair_set


def exercise_ensure_hydrogens_adds_h_to_h_less():
  """_ensure_hydrogens_on_reference places hydrogens on an H-less reference
  hierarchy and returns the augmented one."""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=_pdb_without_h)
  ref_h = pdb_inp.construct_hierarchy()
  # Confirm precondition: this hierarchy has no H
  m0 = mmtbx.model.manager(model_input=None,
                           pdb_hierarchy=ref_h.deep_copy(),
                           log=null_out())
  assert not m0.has_hd(), "test precondition: ref_h must lack H"
  augmented = _ensure_hydrogens_on_reference(ref_h, log=null_out())
  assert augmented is not None
  # The returned hierarchy passes has_hd()
  m1 = mmtbx.model.manager(model_input=None,
                           pdb_hierarchy=augmented.deep_copy(),
                           log=null_out())
  assert m1.has_hd(), "returned hierarchy should have H atoms"
  assert augmented.atoms_size() > ref_h.atoms_size(), \
    "atoms_size should grow: %d -> %d" % (
      ref_h.atoms_size(), augmented.atoms_size())


def run(args):
  assert not args, args
  exercise_origin_id_registered()
  exercise_phil_defaults()
  exercise_get_hbond_proxies_returns_empty_when_disabled()
  exercise_ensure_hydrogens_idempotent_on_h_bearing()
  exercise_ensure_hydrogens_adds_h_to_h_less()
  exercise_detect_reference_hbonds()
  exercise_collect_ss_hbond_pairs()
  exercise_end_to_end_as_found_no_angles_no_ss()
  exercise_end_to_end_target_ideal()
  exercise_ss_filter_excludes_overlap()
  exercise_angle_proxies_from_restrain_angles()
  exercise_h_less_working_skips_angles_keeps_bonds()
  exercise_warning_when_angles_requested_but_no_working_h()
  exercise_ncs_intra_chain_expansion()
  exercise_partner_distance_cutoff_rejects()
  exercise_function_rename_and_dual_dispatch()
  exercise_integration_via_model_process()
  exercise_geo_output_labels_new_origin_id()
  print("OK")


if __name__ == "__main__":
  run(sys.argv[1:])
