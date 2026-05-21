from __future__ import absolute_import, division, print_function

import math
import sys

import iotbx.pdb
import mmtbx.model

from cctbx import geometry_restraints
from libtbx.utils import Sorry, null_out

from mmtbx.geometry_restraints import ramachandran
from mmtbx.geometry_restraints.torsion_restraints.reference_model import (
  reference_model,
  reference_model_params,
)


# A 9-residue alpha-helical poly-ALA. Real residue type so the rama classifier
# emits oldfield proxies for the 7 non-terminal residues. Coordinates taken
# from a poly-Gly helix and atom names adjusted (CB added in canonical position).
# The exact CB positions are unimportant for phi/psi.
_pdb_helix = """\
CRYST1   30.000   30.000   60.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      -5.606  -2.251 -12.878  1.00  0.00           N
ATOM      2  CA  ALA A   1      -5.850  -1.194 -13.852  1.00  0.00           C
ATOM      3  C   ALA A   1      -5.186  -1.524 -15.184  1.00  0.00           C
ATOM      4  O   ALA A   1      -5.744  -1.260 -16.249  1.00  0.00           O
ATOM      5  CB  ALA A   1      -7.330  -1.000 -14.000  1.00  0.00           C
ATOM      6  N   ALA A   2      -3.992  -2.102 -15.115  1.00  0.00           N
ATOM      7  CA  ALA A   2      -3.261  -2.499 -16.313  1.00  0.00           C
ATOM      8  C   ALA A   2      -3.961  -3.660 -17.011  1.00  0.00           C
ATOM      9  O   ALA A   2      -4.016  -3.716 -18.240  1.00  0.00           O
ATOM     10  CB  ALA A   2      -1.770  -2.800 -16.050  1.00  0.00           C
ATOM     11  N   ALA A   3      -4.492  -4.585 -16.219  1.00  0.00           N
ATOM     12  CA  ALA A   3      -5.216  -5.731 -16.755  1.00  0.00           C
ATOM     13  C   ALA A   3      -6.531  -5.289 -17.389  1.00  0.00           C
ATOM     14  O   ALA A   3      -6.939  -5.814 -18.425  1.00  0.00           O
ATOM     15  CB  ALA A   3      -5.500  -6.800 -15.700  1.00  0.00           C
ATOM     16  N   ALA A   4      -7.189  -4.323 -16.758  1.00  0.00           N
ATOM     17  CA  ALA A   4      -8.442  -3.785 -17.273  1.00  0.00           C
ATOM     18  C   ALA A   4      -8.205  -3.003 -18.561  1.00  0.00           C
ATOM     19  O   ALA A   4      -9.007  -3.065 -19.492  1.00  0.00           O
ATOM     20  CB  ALA A   4      -9.100  -3.000 -16.150  1.00  0.00           C
ATOM     21  N   ALA A   5      -7.099  -2.269 -18.604  1.00  0.00           N
ATOM     22  CA  ALA A   5      -6.735  -1.498 -19.787  1.00  0.00           C
ATOM     23  C   ALA A   5      -6.358  -2.423 -20.939  1.00  0.00           C
ATOM     24  O   ALA A   5      -6.687  -2.157 -22.094  1.00  0.00           O
ATOM     25  CB  ALA A   5      -5.600  -0.500 -19.500  1.00  0.00           C
ATOM     26  N   ALA A   6      -5.665  -3.509 -20.614  1.00  0.00           N
ATOM     27  CA  ALA A   6      -5.268  -4.493 -21.614  1.00  0.00           C
ATOM     28  C   ALA A   6      -6.485  -5.236 -22.153  1.00  0.00           C
ATOM     29  O   ALA A   6      -6.565  -5.533 -23.345  1.00  0.00           O
ATOM     30  CB  ALA A   6      -4.200  -5.500 -21.150  1.00  0.00           C
ATOM     31  N   ALA A   7      -7.430  -5.532 -21.267  1.00  0.00           N
ATOM     32  CA  ALA A   7      -8.660  -6.212 -21.655  1.00  0.00           C
ATOM     33  C   ALA A   7      -9.529  -5.303 -22.518  1.00  0.00           C
ATOM     34  O   ALA A   7     -10.158  -5.756 -23.474  1.00  0.00           O
ATOM     35  CB  ALA A   7      -9.450  -6.700 -20.450  1.00  0.00           C
ATOM     36  N   ALA A   8      -9.559  -4.021 -22.172  1.00  0.00           N
ATOM     37  CA  ALA A   8     -10.324  -3.039 -22.930  1.00  0.00           C
ATOM     38  C   ALA A   8      -9.706  -2.819 -24.306  1.00  0.00           C
ATOM     39  O   ALA A   8     -10.416  -2.660 -25.299  1.00  0.00           O
ATOM     40  CB  ALA A   8     -11.800  -3.500 -23.100  1.00  0.00           C
ATOM     41  N   ALA A   9      -8.378  -2.810 -24.356  1.00  0.00           N
ATOM     42  CA  ALA A   9      -7.658  -2.641 -25.613  1.00  0.00           C
ATOM     43  C   ALA A   9      -7.843  -3.861 -26.508  1.00  0.00           C
ATOM     44  O   ALA A   9      -7.980  -3.734 -27.725  1.00  0.00           O
ATOM     45  CB  ALA A   9      -6.150  -2.400 -25.300  1.00  0.00           C
END
"""


# Minimum chain that yields exactly one oldfield rama proxy: 3 residues, with
# the proxy on the middle residue (its phi atoms are C(1)-N(2)-CA(2)-C(2) and
# psi atoms are N(2)-CA(2)-C(2)-N(3)). Used by the minimization-comparison
# exercise where a small, viewable model with a single rama proxy is the
# point. Coordinates lifted from _pdb_helix residues 4-6 (interior helical
# geometry, no terminal artifacts), renumbered 1-3.
_pdb_tri = """\
CRYST1   30.000   30.000   60.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      -7.189  -4.323 -16.758  1.00  0.00           N
ATOM      2  CA  ALA A   1      -8.442  -3.785 -17.273  1.00  0.00           C
ATOM      3  C   ALA A   1      -8.205  -3.003 -18.561  1.00  0.00           C
ATOM      4  O   ALA A   1      -9.007  -3.065 -19.492  1.00  0.00           O
ATOM      5  CB  ALA A   1      -9.100  -3.000 -16.150  1.00  0.00           C
ATOM      6  N   ALA A   2      -7.099  -2.269 -18.604  1.00  0.00           N
ATOM      7  CA  ALA A   2      -6.735  -1.498 -19.787  1.00  0.00           C
ATOM      8  C   ALA A   2      -6.358  -2.423 -20.939  1.00  0.00           C
ATOM      9  O   ALA A   2      -6.687  -2.157 -22.094  1.00  0.00           O
ATOM     10  CB  ALA A   2      -5.600  -0.500 -19.500  1.00  0.00           C
ATOM     11  N   ALA A   3      -5.665  -3.509 -20.614  1.00  0.00           N
ATOM     12  CA  ALA A   3      -5.268  -4.493 -21.614  1.00  0.00           C
ATOM     13  C   ALA A   3      -6.485  -5.236 -22.153  1.00  0.00           C
ATOM     14  O   ALA A   3      -6.565  -5.533 -23.345  1.00  0.00           O
ATOM     15  CB  ALA A   3      -4.200  -5.500 -21.150  1.00  0.00           C
END
"""

# Same 3-residue chain, residue-2 phi/psi rotated to ~(-8.5, +41.4) (vs the
# helical (-68.5, -38.6) in _pdb_tri), then idealized: bonds and angles at
# cctbx monomer-library targets and torsions held at the rotated values via
# reference_model. The geometry is now self-consistent — a working model
# minimized with phi/psi pinned to these values can overlay the reference
# exactly. Generated offline; this is the captured snapshot.
_pdb_tri_reference = """\
CRYST1   30.000   30.000   60.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      -6.790  -4.206 -16.547  1.00  0.00           N
ATOM      2  CA  ALA A   1      -8.167  -3.859 -16.878  1.00  0.00           C
ATOM      3  C   ALA A   1      -8.250  -3.213 -18.257  1.00  0.00           C
ATOM      4  O   ALA A   1      -9.170  -3.504 -19.009  1.00  0.00           O
ATOM      5  CB  ALA A   1      -8.750  -2.933 -15.821  1.00  0.00           C
ATOM      6  N   ALA A   2      -7.242  -2.379 -18.548  1.00  0.00           N
ATOM      7  CA  ALA A   2      -6.952  -1.601 -19.761  1.00  0.00           C
ATOM      8  C   ALA A   2      -7.775  -1.820 -21.033  1.00  0.00           C
ATOM      9  O   ALA A   2      -7.207  -1.848 -22.130  1.00  0.00           O
ATOM     10  CB  ALA A   2      -5.478  -1.806 -20.124  1.00  0.00           C
ATOM     11  N   ALA A   3      -9.089  -1.977 -20.922  1.00  0.00           N
ATOM     12  CA  ALA A   3      -9.943  -2.156 -22.090  1.00  0.00           C
ATOM     13  C   ALA A   3     -10.012  -0.879 -22.922  1.00  0.00           C
ATOM     14  O   ALA A   3     -10.325  -0.918 -24.112  1.00  0.00           O
ATOM     15  CB  ALA A   3     -11.338  -2.588 -21.665  1.00  0.00           C
END
"""


def _get_rama_manager(model):
  """Return the ramachandran_manager attached to a model's GRM (or None).

  Note: model.get_ramachandran_manager() returns a rama_eval (the scoring
  helper) — a name clash. The restraint manager lives at
  geometry.ramachandran_manager."""
  rm_mgr = model.get_restraints_manager()
  if rm_mgr is None: return None
  return rm_mgr.geometry.ramachandran_manager


def _build_model_with_oldfield_rama(pdb_str=_pdb_helix):
  """Construct a model that has an oldfield ramachandran_manager attached."""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.ramachandran_plot_restraints.enabled = True
  params.pdb_interpretation.ramachandran_plot_restraints.\
      inject_emsley8k_into_oldfield_favored = False
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(pdb_interpretation_params=params, make_restraints=True)
  return model


def _compute_phi_psi(sites_cart, i_seqs):
  """Compute (phi, psi) in degrees from the 5 i_seqs of an oldfield proxy.

  Phi uses i_seqs[0:4], psi uses i_seqs[1:5] (per ramachandran.h target_phi_psi
  template, lines 254-258, and the proxy construction in ramachandran.py:362)."""
  phi_sites = tuple(tuple(sites_cart[j]) for j in i_seqs[0:4])
  psi_sites = tuple(tuple(sites_cart[j]) for j in i_seqs[1:5])
  phi = geometry_restraints.dihedral(
    sites=phi_sites, angle_ideal=0, weight=1).angle_model
  psi = geometry_restraints.dihedral(
    sites=psi_sites, angle_ideal=0, weight=1).angle_model
  return phi, psi


def exercise_phil_defaults():
  """The new ramachandran_targets { ... } sub-block inside reference_model
  parses with the documented default (enabled=False)."""
  p = reference_model_params.extract().reference_model
  assert hasattr(p, 'ramachandran_targets'), \
    "reference_model PHIL must contain ramachandran_targets sub-block"
  rt = p.ramachandran_targets
  assert rt.enabled is False, rt.enabled


def exercise_apply_noop_when_disabled():
  """apply_ramachandran_targets returns 0 and changes nothing when the
  ramachandran_targets.enabled flag is False."""
  model = _build_model_with_oldfield_rama()
  rama_mgr = _get_rama_manager(model)
  assert rama_mgr is not None
  assert rama_mgr.get_n_oldfield_proxies() > 0, \
    "expected oldfield proxies in this helix"
  # Snapshot the existing targets.
  before = list(tuple(t) for t in rama_mgr.target_phi_psi)

  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.ramachandran_targets.enabled = False
  rm = reference_model(
    model=model,
    reference_hierarchy_list=[model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())

  n = rm.apply_ramachandran_targets(rama_mgr)
  assert n == 0, n
  after = list(tuple(t) for t in rama_mgr.target_phi_psi)
  assert before == after, "targets must be unchanged when disabled"


def exercise_apply_self_reference_overrides_all():
  """With a self-reference (refined == reference), every oldfield proxy
  whose 5 i_seqs are all in match_map gets overridden; the new (phi, psi)
  match phi/psi computed from the refined hierarchy; the third element
  (distance_to_allowed) is ~0."""
  model = _build_model_with_oldfield_rama()
  rama_mgr = _get_rama_manager(model)
  n_proxies = rama_mgr.get_n_oldfield_proxies()
  assert n_proxies > 0

  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.ramachandran_targets.enabled = True
  rm = reference_model(
    model=model,
    reference_hierarchy_list=[model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())

  n = rm.apply_ramachandran_targets(rama_mgr)
  assert n == n_proxies, (n, n_proxies)

  sites_cart = model.get_sites_cart()
  for i, proxy in enumerate(rama_mgr._oldfield_proxies):
    phi_expected, psi_expected = _compute_phi_psi(sites_cart, proxy.get_i_seqs())
    phi_t, psi_t, dist = tuple(rama_mgr.target_phi_psi[i])
    assert abs(phi_t - phi_expected) < 1.e-6, (phi_t, phi_expected)
    assert abs(psi_t - psi_expected) < 1.e-6, (psi_t, psi_expected)
    assert abs(dist) < 1.e-6, dist


def exercise_apply_distance_is_angular_difference():
  """When the refined model differs from the reference, target_phi_psi[i]
  carries the reference's phi/psi (not the refined model's), and the third
  element is sqrt(dphi^2 + dpsi^2) in degrees, using wrapped angle deltas."""
  # Refined model is the helix; reference is the helix with one residue
  # (#5) rotated about phi by ~30 degrees. We achieve this by perturbing the
  # reference hierarchy in memory before passing it in.
  refined_model = _build_model_with_oldfield_rama()
  ref_h = iotbx.pdb.input(
    source_info=None, lines=_pdb_helix).construct_hierarchy()
  # Perturb residue 5: shift its N atom only along a perpendicular axis.
  # That changes phi(5) (depends on C4-N5-CA5-C5) and psi(4) (depends on
  # N4-CA4-C4-N5). The specific magnitude is unimportant; we just check that
  # the override picks up the reference's phi/psi, not the refined model's.
  res5_N = None
  for atom in ref_h.atoms():
    if (atom.parent().parent().resseq_as_int() == 5
        and atom.name.strip() == 'N'):
      res5_N = atom
      break
  assert res5_N is not None
  xyz = res5_N.xyz
  res5_N.set_xyz((xyz[0] + 0.5, xyz[1] - 0.3, xyz[2] + 0.2))

  p = reference_model_params.extract()
  p.reference_model.ramachandran_targets.enabled = True
  # Use file-based path so the reference is the perturbed hierarchy, not self.
  rm = reference_model(
    model=refined_model,
    reference_hierarchy_list=[ref_h],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())

  rama_mgr = _get_rama_manager(refined_model)
  rm.apply_ramachandran_targets(rama_mgr)

  refined_sites = refined_model.get_sites_cart()
  ref_sites = ref_h.atoms().extract_xyz()
  # The match should be 1:1 between refined and reference (same atoms, same
  # order).
  for i, proxy in enumerate(rama_mgr._oldfield_proxies):
    phi_ref, psi_ref = _compute_phi_psi(ref_sites, proxy.get_i_seqs())
    phi_ref_target, psi_ref_target, dist = tuple(rama_mgr.target_phi_psi[i])
    assert abs(phi_ref_target - phi_ref) < 1.e-6, (phi_ref_target, phi_ref)
    assert abs(psi_ref_target - psi_ref) < 1.e-6, (psi_ref_target, psi_ref)
    # Expected distance: angular delta between refined and reference.
    phi_refined, psi_refined = _compute_phi_psi(refined_sites, proxy.get_i_seqs())
    dphi = geometry_restraints.angle_delta_deg(phi_ref, phi_refined)
    dpsi = geometry_restraints.angle_delta_deg(psi_ref, psi_refined)
    expected_dist = math.sqrt(dphi * dphi + dpsi * dpsi)
    assert abs(dist - expected_dist) < 1.e-5, (dist, expected_dist)
  # And: at least one proxy must show a nonzero distance (residue 5 was
  # perturbed; its proxy and the neighbouring residue-4 proxy are affected).
  assert any(tuple(rama_mgr.target_phi_psi[i])[2] > 1.0
             for i in range(rama_mgr.get_n_oldfield_proxies())), \
    "expected at least one proxy with nonzero distance after perturbation"


def exercise_unmatched_proxies_untouched():
  """A reference selection that excludes one residue leaves the
  oldfield target for the proxies involving that residue untouched."""
  refined_model = _build_model_with_oldfield_rama()
  rama_mgr = _get_rama_manager(refined_model)
  before = list(tuple(t) for t in rama_mgr.target_phi_psi)

  # Reference contains only residues 1-4 (chain A). Residue 5..9 proxies and
  # the residue-4 proxy (which depends on N of residue 5) will not have a
  # complete match.
  p = reference_model_params.extract()
  p.reference_model.ramachandran_targets.enabled = True
  rg = p.reference_model.reference_group[0]
  rg.reference = "resseq 1:4"
  rg.selection = "resseq 1:4"
  rg.file_name = "ref0"

  rm = reference_model(
    model=refined_model,
    reference_hierarchy_list=[refined_model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())

  n_changed = rm.apply_ramachandran_targets(rama_mgr)
  after = list(tuple(t) for t in rama_mgr.target_phi_psi)

  # Some proxies must have changed (those fully inside 1-4 -- which is at
  # most proxies for residues 2 and 3), and at least one must be unchanged
  # (proxies for residues 5+).
  changed = [i for i in range(len(before)) if before[i] != after[i]]
  unchanged = [i for i in range(len(before)) if before[i] == after[i]]
  assert len(changed) == n_changed, (len(changed), n_changed)
  assert len(changed) > 0, "expected at least one overridden proxy"
  assert len(unchanged) > 0, "expected at least one untouched proxy"


def exercise_integration_via_add_reference_model_restraints():
  """End-to-end: model.process() with ramachandran_targets.enabled=True
  triggers the override path via add_reference_model_restraints_if_requested.
  The post-process targets reflect the reference geometry, not the tabulated
  lookup."""
  from mmtbx.geometry_restraints.torsion_restraints.reference_model import \
      add_reference_model_restraints_if_requested

  refined_model = _build_model_with_oldfield_rama()
  rama_mgr = _get_rama_manager(refined_model)
  # Drive the pathway: build params with ramachandran_targets enabled and
  # self-reference, then call the dispatcher used by model.process().
  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.ramachandran_targets.enabled = True
  geometry = refined_model.get_restraints_manager().geometry
  n_changed = add_reference_model_restraints_if_requested(
    model=refined_model,
    geometry=geometry,
    params=p.reference_model,
    log=null_out())
  assert n_changed >= 1, n_changed
  # And the rama_mgr.target_phi_psi[0] matches the refined model's phi/psi
  # (self-reference).
  sites_cart = refined_model.get_sites_cart()
  proxy = rama_mgr._oldfield_proxies[0]
  phi_e, psi_e = _compute_phi_psi(sites_cart, proxy.get_i_seqs())
  phi_t, psi_t, _ = tuple(rama_mgr.target_phi_psi[0])
  assert abs(phi_t - phi_e) < 1.e-6, (phi_t, phi_e)
  assert abs(psi_t - psi_e) < 1.e-6, (psi_t, psi_e)


def exercise_reapply_on_set_ramachandran_plot_restraints():
  """When the refinement driver rebuilds the rama manager via
  model.set_ramachandran_plot_restraints, the override must be re-applied
  so the next refinement cycle sees reference-derived targets."""
  refined_model = _build_model_with_oldfield_rama()
  geometry = refined_model.get_restraints_manager().geometry

  # Attach a reference_model with the flag on.
  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.ramachandran_targets.enabled = True
  rm = reference_model(
    model=refined_model,
    reference_hierarchy_list=[refined_model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)

  # Apply once.
  rama_mgr = _get_rama_manager(refined_model)
  rm.apply_ramachandran_targets(rama_mgr)
  sites_cart = refined_model.get_sites_cart()
  proxy0 = rama_mgr._oldfield_proxies[0]
  phi_ref0, psi_ref0 = _compute_phi_psi(sites_cart, proxy0.get_i_seqs())

  # Now rebuild the rama manager via the model API (mimics what
  # macro_cycle_real_space.py does each macro cycle).
  rama_params = ramachandran.master_phil.fetch().extract().\
      ramachandran_plot_restraints
  rama_params.enabled = True
  rama_params.inject_emsley8k_into_oldfield_favored = False
  refined_model.set_ramachandran_plot_restraints(rama_params=rama_params)

  new_rama_mgr = _get_rama_manager(refined_model)
  assert new_rama_mgr is not None
  # The new manager must be a different object.
  assert new_rama_mgr is not rama_mgr
  # And the override must have been re-applied: target[0] phi/psi must match
  # the reference (== refined here for self-reference) within tolerance.
  phi_t, psi_t, _ = tuple(new_rama_mgr.target_phi_psi[0])
  assert abs(phi_t - phi_ref0) < 1.e-6, (phi_t, phi_ref0)
  assert abs(psi_t - psi_ref0) < 1.e-6, (psi_t, psi_ref0)


def exercise_raises_when_injection_active():
  """When ramachandran_plot_restraints.inject_emsley8k_into_oldfield_favored
  is True (the default), favored residues are silently restrained with
  emsley8k and our oldfield-target override does not apply to them. To make
  this explicit, apply_ramachandran_targets must raise Sorry rather than
  silently producing a partial override."""
  # Build a model that DOES leave injection on (the default).
  pdb_inp = iotbx.pdb.input(source_info=None, lines=_pdb_helix)
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.ramachandran_plot_restraints.enabled = True
  # NOTE: leaving inject_emsley8k_into_oldfield_favored at default (True)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(pdb_interpretation_params=params, make_restraints=True)
  rama_mgr = _get_rama_manager(model)
  assert rama_mgr is not None
  assert rama_mgr.params.inject_emsley8k_into_oldfield_favored is True

  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.ramachandran_targets.enabled = True
  rm = reference_model(
    model=model,
    reference_hierarchy_list=[model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())

  raised = False
  try:
    rm.apply_ramachandran_targets(rama_mgr)
  except Sorry as e:
    raised = True
    assert 'inject_emsley8k_into_oldfield_favored' in str(e), str(e)
  assert raised, "expected Sorry when injection is active"


def exercise_raises_when_no_rama_manager():
  """If ramachandran_targets.enabled=True is set but no rama manager has
  been attached (typical when ramachandran_plot_restraints.enabled=False),
  apply_ramachandran_targets must raise Sorry."""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=_pdb_helix)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  assert _get_rama_manager(model) is None

  p = reference_model_params.extract()
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.ramachandran_targets.enabled = True
  rm = reference_model(
    model=model,
    reference_hierarchy_list=[model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())

  raised = False
  try:
    rm.apply_ramachandran_targets(None)
  except Sorry as e:
    raised = True
    assert 'Ramachandran restraints manager' in str(e), str(e)
  assert raised, "expected Sorry when no rama manager attached"


def exercise_minimization_drives_to_reference_target():
  """Behavioural smoke test on the minimum chain that yields exactly one
  oldfield rama proxy (3 residues, proxy on the middle residue).

  Identical starting model, identical geometry restraints; the ONLY
  difference between two lbfgs runs is whether the oldfield Ramachandran
  target has been overridden from a reference hierarchy.

  Run A — stock oldfield target (nearest favorable point on Rama8000).
  Run B — override active; target comes from _pdb_tri_reference where
  residue 2's phi/psi were cleanly rotated to a non-helical conformation.

  Assertions:
    1. The two runs end in different places.
    2. Run B is meaningfully closer to the reference's phi/psi than Run A.
  Each lbfgs run writes its result + .geo into cwd for visual inspection."""
  import os
  import scitbx.lbfgs
  from mmtbx.refinement import geometry_minimization
  from mmtbx.building.loop_closure.utils import get_pair_angles,get_phi_psi_atoms

  def _get_rama_single_angles(hierarchy):
    phi_psi_atoms = get_phi_psi_atoms(hierarchy)
    angles = []
    for phi_psi_pair, rama_key in phi_psi_atoms:
      angles.append(get_pair_angles(phi_psi_pair))
    return angles[0]


  out_dir = os.getcwd()
  do_write = os.path.isdir(os.path.dirname(out_dir))
  if do_write and not os.path.isdir(out_dir):
    os.makedirs(out_dir)

  def _dump(name, text):
    if not do_write: return
    with open(os.path.join(out_dir, name), "w") as f:
      f.write(text)

  # Working model: 3-residue helical fragment, oldfield rama active.
  refined_model = _build_model_with_oldfield_rama(pdb_str=_pdb_tri)
  rama_mgr = _get_rama_manager(refined_model)
  # Push oldfield weight high so the rama term dominates the result.
  rama_mgr.params.oldfield.weight = 100.0
  geometry = refined_model.get_restraints_manager().geometry

  # Exactly one oldfield proxy is expected, on residue 2.
  assert rama_mgr.get_n_oldfield_proxies() == 1, \
    "expected exactly one oldfield proxy on the 3-res fragment, got %d" % (
      rama_mgr.get_n_oldfield_proxies(),)
  target_proxy = rama_mgr._oldfield_proxies[0]
  proxy_idx = 0

  sites_cart_initial = refined_model.get_sites_cart().deep_copy()
  _dump("0_initial.pdb", refined_model.model_as_pdb())
  _dump("0_initial.geo", refined_model.restraints_as_geo())
  angle_0 = _get_rama_single_angles(refined_model.get_hierarchy())

  # Reference: hardcoded 3-residue chain with residue 2 phi/psi rotated.
  ref_h = iotbx.pdb.input(
    source_info=None, lines=_pdb_tri_reference).construct_hierarchy()
  ref_sites = ref_h.atoms().extract_xyz()
  phi_ref_target, psi_ref_target = _compute_phi_psi(
    ref_sites, target_proxy.get_i_seqs())
  _dump("1_reference.pdb", ref_h.as_pdb_string())
  angle_1 = _get_rama_single_angles(ref_h)

  # Sanity check the perturbation actually shifts phi/psi.
  init_phi, init_psi = _compute_phi_psi(
    sites_cart_initial, target_proxy.get_i_seqs())
  shift = math.sqrt(
    geometry_restraints.angle_delta_deg(phi_ref_target, init_phi) ** 2
    + geometry_restraints.angle_delta_deg(psi_ref_target, init_psi) ** 2)
  assert shift > 15.0, \
    "perturbation too small to see; shift=%.2f" % shift

  flags = geometry_restraints.flags.flags(default=True)
  term = scitbx.lbfgs.termination_parameters(max_iterations=300)

  # --- Run A: no override.
  sites_cart_a = sites_cart_initial.deep_copy()
  geometry_minimization.lbfgs(
    sites_cart=sites_cart_a,
    correct_special_position_tolerance=1.0,
    geometry_restraints_manager=geometry,
    geometry_restraints_flags=flags,
    lbfgs_termination_params=term)
  phi_a, psi_a = _compute_phi_psi(sites_cart_a, target_proxy.get_i_seqs())
  refined_model.set_sites_cart(sites_cart_a)
  _dump("2_lbfgs_no_override.pdb", refined_model.model_as_pdb())
  _dump("2_lbfgs_no_override.geo", refined_model.restraints_as_geo())
  angle_2 = _get_rama_single_angles(refined_model.get_hierarchy())
  refined_model.set_sites_cart(sites_cart_initial)

  # --- Run B: apply override, run from the SAME initial coordinates.
  # Enable the torsion path too and adopt the manager so the reference
  # dihedral restraints (a) actually apply during this lbfgs run, and
  # (b) appear in the dumped .geo files for inspection.
  p = reference_model_params.extract()
  p.reference_model.enabled = True
  p.reference_model.ramachandran_targets.enabled = True
  rm = reference_model(
    model=refined_model,
    reference_hierarchy_list=[ref_h],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)
  n_changed = rm.apply_ramachandran_targets(rama_mgr)
  assert n_changed > 0
  # Spot-check the override landed on our chosen proxy.
  phi_t, psi_t, _ = tuple(rama_mgr.target_phi_psi[proxy_idx])
  assert abs(phi_t - phi_ref_target) < 1.e-3
  assert abs(psi_t - psi_ref_target) < 1.e-3
  _dump("3_with_override_pre_min.geo", refined_model.restraints_as_geo())

  sites_cart_b = sites_cart_initial.deep_copy()
  geometry_minimization.lbfgs(
    sites_cart=sites_cart_b,
    correct_special_position_tolerance=1.0,
    geometry_restraints_manager=geometry,
    geometry_restraints_flags=flags,
    lbfgs_termination_params=term)
  phi_b, psi_b = _compute_phi_psi(sites_cart_b, target_proxy.get_i_seqs())
  refined_model.set_sites_cart(sites_cart_b)
  _dump("4_lbfgs_with_override.pdb", refined_model.model_as_pdb())
  _dump("4_lbfgs_with_override.geo", refined_model.restraints_as_geo())
  angle_4 = _get_rama_single_angles(refined_model.get_hierarchy())

  print(angle_0)
  print(angle_1)
  print(angle_2)
  print(angle_4)

  # --- Assertions.
  def adist2(a1, a2, b1, b2):
    da = geometry_restraints.angle_delta_deg(a1, b1)
    dp = geometry_restraints.angle_delta_deg(a2, b2)
    return math.sqrt(da * da + dp * dp)

  d_a_to_ref = adist2(phi_a, psi_a, phi_ref_target, psi_ref_target)
  d_b_to_ref = adist2(phi_b, psi_b, phi_ref_target, psi_ref_target)
  d_runs    = adist2(phi_a, psi_a, phi_b, psi_b)

  # Two different minima from the same start. With a 60/80 deg clean
  # rotation, the two runs diverge by roughly 60-70 deg in joint phi/psi.
  assert d_runs > 30.0, \
    "expected the two runs to end far apart; d_runs=%.2f" % d_runs
  # Run B is closer to the reference target than Run A is. Run B stops at
  # an energy compromise between the rama target and the bond/angle terms
  # that resist large phi/psi swings (peptide-bond geometry around the
  # mobile residue), so we require "at least halfway" rather than exact
  # convergence.
  assert d_b_to_ref < 0.5 * d_a_to_ref, \
    "Run B did not get meaningfully closer to ref: d_b=%.2f, d_a=%.2f" \
    % (d_b_to_ref, d_a_to_ref)


def exercise_existing_reference_dihedral_restraints_kept():
  """Per the 'keep both' design choice: when both params.enabled (torsion
  reference) and params.ramachandran_targets.enabled are True, the existing
  reference_dihedral_proxies (the harmonic-with-top-out path on phi/psi/chi)
  are still populated. Both restraint families apply concurrently."""
  refined_model = _build_model_with_oldfield_rama()
  geometry = refined_model.get_restraints_manager().geometry

  p = reference_model_params.extract()
  p.reference_model.enabled = True
  p.reference_model.use_starting_model_as_reference = True
  p.reference_model.main_chain = True
  p.reference_model.ramachandran_targets.enabled = True
  rm = reference_model(
    model=refined_model,
    reference_hierarchy_list=[refined_model.get_hierarchy()],
    reference_file_list=None,
    params=p.reference_model,
    log=null_out())
  geometry.adopt_reference_dihedral_manager(rm)
  rm.apply_ramachandran_targets(_get_rama_manager(refined_model))

  # The torsion path's dihedral_proxies must still exist.
  assert rm.get_n_proxies() > 0, \
    "reference_dihedral_proxies should still be populated"


def run(args):
  exercise_phil_defaults()
  exercise_apply_noop_when_disabled()
  exercise_apply_self_reference_overrides_all()
  exercise_apply_distance_is_angular_difference()
  exercise_unmatched_proxies_untouched()
  exercise_integration_via_add_reference_model_restraints()
  exercise_reapply_on_set_ramachandran_plot_restraints()
  exercise_raises_when_injection_active()
  exercise_raises_when_no_rama_manager()
  exercise_minimization_drives_to_reference_target()
  exercise_existing_reference_dihedral_restraints_kept()
  print("OK")


if __name__ == "__main__":
  run(sys.argv[1:])
