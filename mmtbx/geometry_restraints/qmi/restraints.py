"""Assemble ``geometry_restraints.edits`` from a QM-optimised metal cluster.

Reads the GFN2-xTB-optimised geometry of an endo_exo sub-model and emits metal
bond + angle restraints in the Phenix ``geometry_restraints.edits`` PHIL form
(plan Appendix A): one bond per metal-donor contact, one angle per
donor-metal-donor pair, with the target taken from the optimised geometry.

Donor identification uses the *input* (crystal) coordinates (coordination is
a property of the deposited structure), while the restraint targets are
measured in the *optimised* geometry.

Each restraint sigma is the modified-Seminario curvature (plan Sec 5F); when the
curvature is unavailable (Hessian failure, or still an imaginary mode after the
one perturb-and-reoptimise) it falls back to a conservative M-L default
(:data:`DEFAULT_BOND_SIGMA` / :data:`DEFAULT_ANGLE_SIGMA`) flagged
``conservative_fallback``.  Atom selections are built from atom identity (chain /
resseq / resname / name / altloc) so the ``.edits`` file applies to the full
input model at refinement time; donors on symmetry images (renamed chains in the
endo_exo sub-model) cannot be expressed this way and are flagged, not emitted.
"""

from __future__ import absolute_import, division, print_function

from libtbx.utils import null_out
from scitbx import matrix

from mmtbx.geometry_restraints.qmi import util


# Biological donor elements (carbon excluded by default: at ~3 A the His ring
# carbons sit at the cutoff and would be picked up as spurious donors).  Same
# whitelist as the qmi coordination-shell convention.
DEFAULT_DONOR_ELEMENTS = frozenset({
  "O", "N", "S", "SE", "F", "CL", "BR", "I", "P",
})

# Distance defining a first-shell contact (Harding / CheckMyMetal convention);
# the single value lives in qmi.util so the pipeline cannot silently diverge.
DEFAULT_COORDINATION_CUTOFF_A = util.DEFAULT_COORDINATION_CUTOFF_A

# Conservative M-L fallback sigmas, used when the modified-Seminario curvature is
# unavailable (plan Sec 5F/9.2). Deliberately looser than the covalent-calibrated
# Phenix defaults (0.02 A / 3 deg): metal-ligand bonds are softer (lower force
# constants), so an M-L bond we could not characterise should not be restrained as
# tightly as a covalent one. An engineering default, NOT a database value; it sits
# inside the computed Seminario clamp range so a failed Hessian is never tighter
# than a soft computed bond.
DEFAULT_BOND_SIGMA = 0.05      # A
DEFAULT_ANGLE_SIGMA = 5.0      # deg


def detect_first_shell(submodel, metal_iseq, sites_cart=None,
                       cutoff_A=DEFAULT_COORDINATION_CUTOFF_A,
                       donor_elements=DEFAULT_DONOR_ELEMENTS):
  """Return the donor i_seqs within *cutoff_A* of the metal.

  Coordination is read from *sites_cart* (default: the sub-model's input
  coordinates, the crystal geometry), not the optimised geometry.  Hydrogens
  are excluded, and only elements in *donor_elements* are admitted.

  Parameters
  ----------
  submodel : mmtbx.model.manager
  metal_iseq : int
  sites_cart : scitbx.array_family.flex.vec3_double or None
      Coordinates to measure against; default ``submodel.get_sites_cart()``.
  cutoff_A : float
  donor_elements : iterable of str

  Returns
  -------
  list of int
      Donor i_seqs, in ascending order.
  """
  if sites_cart is None:
    sites_cart = submodel.get_sites_cart()
  accept = frozenset(e.upper() for e in donor_elements)
  atoms = submodel.get_hierarchy().atoms()
  # All metal-atom distances in one vectorised flex pass; the element/H test
  # stays a per-atom mask.
  dists = (sites_cart - matrix.col(sites_cart[metal_iseq])).norms()
  donors = []
  for i, atom in enumerate(atoms):
    if i == metal_iseq or dists[i] > cutoff_A:
      continue
    elem = atom.element.strip().upper()
    if atom.element_is_hydrogen() or elem not in accept:
      continue
    donors.append(i)
  return donors


def _selection_string(atom):
  """A CCTBX atom-selection string identifying *atom* in the full model.

  Built from chain / resseq / resname / name (+ altloc when present): the
  identity the endo_exo sub-model preserves for identity-image atoms, so the
  selection resolves in the original input model.
  """
  ag = atom.parent()
  rg = ag.parent()
  chain = rg.parent()
  terms = [
    f"chain {chain.id.strip()}",
    f"resseq {rg.resseq.strip()}",
    f"resname {ag.resname.strip()}",
    f"name {atom.name.strip()}",
  ]
  altloc = ag.altloc.strip()
  if altloc:
    terms.append(f'altloc "{altloc}"')
  return " and ".join(terms)


def _is_symmetry_image(atom, original_chain_ids):
  """True if *atom* sits on a chain that endo_exo created for a symmetry image
  (i.e. a chain id not present in the original model), so its selection string
  would not resolve in the input model."""
  return atom.parent().parent().parent().id.strip() not in original_chain_ids


def metal_ligand_restraints(submodel, opt_sites_cart, metal_iseq, donor_iseqs,
                            bond_sigmas=None, angle_sigmas=None,
                            bond_sigma=DEFAULT_BOND_SIGMA,
                            angle_sigma=DEFAULT_ANGLE_SIGMA,
                            original_chain_ids=None, log=None):
  """Build bond + angle restraint records from the optimised geometry.

  One bond per metal-donor contact (``distance_ideal`` = optimised M-donor
  distance) and one angle per donor-metal-donor pair (``angle_ideal`` =
  optimised donor-metal-donor angle).  Donors whose selection cannot be
  expressed against the input model (symmetry images on renamed chains) are
  skipped and counted in the returned ``n_skipped_symmetry``.

  Parameters
  ----------
  submodel : mmtbx.model.manager
      Provides atom identity for the selection strings.
  opt_sites_cart : scitbx.array_family.flex.vec3_double
      Optimised geometry (sub-model atom order); the restraint targets.
  metal_iseq : int
  donor_iseqs : iterable of int
  bond_sigmas : dict {donor_iseq: sigma_A} or None
      Per-bond modified-Seminario sigmas (plan Sec 5F); falls back to
      *bond_sigma* for any donor absent from the dict (or when None).
  angle_sigmas : dict {frozenset({i, j}): sigma_deg} or None
      Per-angle Seminario sigmas; falls back to *angle_sigma*.
  bond_sigma, angle_sigma : float
      Conservative M-L fallback sigmas used when no Seminario value is supplied
      (flagged ``conservative_fallback``).
  original_chain_ids : iterable of str or None
      Chain ids present in the original input model; used to skip
      symmetry-image donors.  ``None`` disables the check (all donors emitted).
  log : file-like or None

  Returns
  -------
  libtbx.group_args
      Fields ``bonds`` (list of dict), ``angles`` (list of dict),
      ``n_skipped_symmetry`` (int).
  """
  from libtbx import group_args
  if log is None:
    log = null_out()
  atoms = submodel.get_hierarchy().atoms()
  metal_atom = atoms[metal_iseq]
  metal_xyz = matrix.col(opt_sites_cart[metal_iseq])
  metal_sel = _selection_string(metal_atom)

  if original_chain_ids is not None:
    original_chain_ids = frozenset(c.strip() for c in original_chain_ids)

  usable = []
  n_skipped = 0
  for d in donor_iseqs:
    if (original_chain_ids is not None
        and _is_symmetry_image(atoms[d], original_chain_ids)):
      donor_id = atoms[d].id_str().replace("pdb=", "")
      print(f"    skipping symmetry-image donor {donor_id} "
            f"(no input-model selection)", file=log)
      n_skipped += 1
      continue
    usable.append(d)

  bonds = []
  for d in usable:
    vec = matrix.col(opt_sites_cart[d]) - metal_xyz
    s = bond_sigmas.get(d) if bond_sigmas else None
    bonds.append(dict(
      atom_selection_1=metal_sel,
      atom_selection_2=_selection_string(atoms[d]),
      donor_element=atoms[d].element.strip(),
      distance_ideal=vec.length(),
      sigma=s if s is not None else bond_sigma,
      sigma_origin="hessian" if s is not None else "conservative_fallback",
    ))

  angles = []
  for a in range(len(usable)):
    va = matrix.col(opt_sites_cart[usable[a]]) - metal_xyz
    for b in range(a + 1, len(usable)):
      vb = matrix.col(opt_sites_cart[usable[b]]) - metal_xyz
      s = (angle_sigmas.get(frozenset((usable[a], usable[b])))
           if angle_sigmas else None)
      angles.append(dict(
        atom_selection_1=_selection_string(atoms[usable[a]]),
        atom_selection_2=metal_sel,
        atom_selection_3=_selection_string(atoms[usable[b]]),
        angle_ideal=va.angle(vb, deg=True),
        sigma=s if s is not None else angle_sigma,
        sigma_origin="hessian" if s is not None else "conservative_fallback",
      ))

  return group_args(bonds=bonds, angles=angles, n_skipped_symmetry=n_skipped)


def as_edits_phil(restraints, header=None):
  """Render restraint records as a ``geometry_restraints.edits`` PHIL string.

  Parameters
  ----------
  restraints : libtbx.group_args
      As returned by :func:`metal_ligand_restraints` (``bonds`` + ``angles``).
  header : str or None
      Optional comment lines prepended to the file (each line is ``# ...``).

  Returns
  -------
  str
  """
  lines = []
  if header:
    for line in header.splitlines():
      lines.append(f"# {line}")
  lines += ["refinement {", "  geometry_restraints {", "    edits {"]
  for b in restraints.bonds:
    lines += [
      "      bond {",
      "        action = *add",
      f"        atom_selection_1 = {b['atom_selection_1']}",
      f"        atom_selection_2 = {b['atom_selection_2']}",
      f"        distance_ideal = {b['distance_ideal']:.4f}",
      f"        sigma = {b['sigma']:.4f}",
      "      }",
    ]
  for a in restraints.angles:
    lines += [
      "      angle {",
      "        action = *add",
      f"        atom_selection_1 = {a['atom_selection_1']}",
      f"        atom_selection_2 = {a['atom_selection_2']}",
      f"        atom_selection_3 = {a['atom_selection_3']}",
      f"        angle_ideal = {a['angle_ideal']:.2f}",
      f"        sigma = {a['sigma']:.2f}",
      "      }",
    ]
  lines += ["    }", "  }", "}", ""]
  return "\n".join(lines)
