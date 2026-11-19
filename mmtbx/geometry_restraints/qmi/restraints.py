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
resid / resname / name / altloc) so the ``.edits`` file applies to the full
input model at refinement time; donors on symmetry images (renamed chains in the
endo_exo sub-model) cannot be expressed this way and are flagged, not emitted.
"""

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
                       coord_cutoff=None,
                       donor_elements=DEFAULT_DONOR_ELEMENTS):
  """Return the donor i_seqs coordinating the metal.

  Coordination is read from *sites_cart* (default: the sub-model's input
  coordinates, the crystal geometry), not the optimised geometry.  A donor
  qualifies when its metal distance is within
  ``coord_cutoff(metal_element, donor_element)``.  Hydrogens are excluded, and
  only elements in *donor_elements* are admitted.

  Parameters
  ----------
  submodel : mmtbx.model.manager
  metal_iseq : int
  sites_cart : scitbx.array_family.flex.vec3_double or None
      Coordinates to measure against; default ``submodel.get_sites_cart()``.
  coord_cutoff : callable or None
      ``(metal_element, donor_element) -> max distance (A)``; see
      :func:`util.make_coordination_cutoff`. Default: the flat fixed cutoff.
  donor_elements : iterable of str

  Returns
  -------
  list of int
      Donor i_seqs, in ascending order.
  """
  if sites_cart is None:
    sites_cart = submodel.get_sites_cart()
  if coord_cutoff is None:
    coord_cutoff = util.make_coordination_cutoff()
  accept = frozenset(e.upper() for e in donor_elements)
  atoms = submodel.get_hierarchy().atoms()
  metal_element = atoms[metal_iseq].element
  # All metal-atom distances in one vectorised flex pass; the element/H test
  # stays a per-atom mask.
  dists = (sites_cart - matrix.col(sites_cart[metal_iseq])).norms()
  donors = []
  for i, atom in enumerate(atoms):
    if i == metal_iseq:
      continue
    elem = atom.element.strip().upper()
    if atom.element_is_hydrogen() or elem not in accept:
      continue
    if dists[i] > coord_cutoff(metal_element, atom.element):
      continue
    donors.append(i)
  return donors


def _selection_string(labels):
  """A CCTBX atom-selection string identifying an atom in the full model, from
  its *labels* (``atom.fetch_labels()``).

  Built from chain / resid / resname / name (+ altloc when present): the
  identity the endo_exo sub-model preserves for identity-image atoms, so the
  selection resolves in the original input model.  resid carries the insertion
  code, so residues 301 and 301A are told apart.  The single formatter for
  both a live atom and a symmetry-image donor's stored ASU-parent labels, so
  the two paths cannot diverge.
  """
  terms = [
    f"chain {labels.chain_id.strip()}",
    f"resid {labels.resid().strip()}",
    f"resname {labels.resname.strip()}",
    f"name {labels.name.strip()}",
  ]
  altloc = labels.altloc.strip()
  if altloc:
    terms.append(f'altloc "{altloc}"')
  return " and ".join(terms)


def is_symmetry_image(atom, original_chain_ids):
  """True if *atom* sits on a chain that endo_exo created for a symmetry image
  (i.e. a chain id not present in the original model), so its selection string
  would not resolve in the input model."""
  return atom.parent().parent().parent().id.strip() not in original_chain_ids


def metal_ligand_restraints(submodel, opt_sites_cart, metal_iseq, donor_iseqs,
                            bond_sigmas=None, angle_sigmas=None,
                            angle_origins=None, emit_angles=True,
                            bond_sigma=DEFAULT_BOND_SIGMA,
                            angle_sigma=DEFAULT_ANGLE_SIGMA,
                            original_chain_ids=None, sym_image_provenance=None,
                            log=None):
  """Build bond + angle restraint records from the optimised geometry.

  One bond per metal-donor contact (``distance_ideal`` = optimised M-donor
  distance) and one angle per donor-metal-donor pair (``angle_ideal`` =
  optimised donor-metal-donor angle).  A symmetry-image donor (on a renamed
  chain) is emitted as a bond to its ASU parent via a ``symmetry_operation`` when
  *sym_image_provenance* supplies it (angles cannot cross symmetry, so it gets no
  angle); a symmetry-image donor without provenance is skipped and counted in the
  returned ``n_skipped_symmetry``.

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
      Per-angle sigmas (plan Sec 5F): trusted Seminario curvatures plus the
      site's widen-only ``seminario_fallback`` for contaminated angles; a key
      absent from the dict (contaminated, with no trusted site angle to widen to)
      takes the conservative *angle_sigma* default, so ``emit_angles`` always
      yields every angle.  ``None`` means the whole-Hessian conservative fallback.
  angle_origins : dict {frozenset({i, j}): str} or None
      Per-angle ``sigma_origin`` (``hessian`` / ``seminario_fallback``) parallel
      to *angle_sigmas*; defaults to ``hessian`` for present keys.
  emit_angles : bool
      Whether to emit donor-metal-donor angle restraints at all (the CShM gate;
      plan Sec 7).
  bond_sigma, angle_sigma : float
      Conservative M-L fallback sigmas used when no Seminario value is supplied
      (flagged ``conservative_fallback``).
  original_chain_ids : iterable of str or None
      Chain ids present in the original input model; used to identify
      symmetry-image donors.  ``None`` disables the check (all donors emitted as
      plain bonds).
  sym_image_provenance : dict {donor_iseq: (asu_parent_labels, symmetry_op_xyz)} or None
      Provenance for symmetry-image donors (from endo_exo); asu_parent_labels is
      the ASU parent's ``atom.fetch_labels()``. A
      symmetry-image donor present here is emitted as a bond to its ASU parent
      with ``symmetry_operation`` set; one absent here is skipped.
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
  metal_sel = _selection_string(metal_atom.fetch_labels())

  if original_chain_ids is not None:
    original_chain_ids = frozenset(c.strip() for c in original_chain_ids)

  prov = sym_image_provenance or {}
  usable = []          # non-symmetry donors: plain bonds + angles
  sym_bonds = []       # (donor_iseq, asu_parent_labels, op_xyz): symmetry bonds only
  n_skipped = 0
  for d in donor_iseqs:
    if (original_chain_ids is not None
        and is_symmetry_image(atoms[d], original_chain_ids)):
      entry = prov.get(d)
      if entry is None:
        donor_id = atoms[d].id_str().replace("pdb=", "")
        print(f"    skipping symmetry-image donor {donor_id} "
              f"(no provenance for a symmetry_operation bond)", file=log)
        n_skipped += 1
        continue
      sym_bonds.append((d, entry[0], entry[1]))
      continue
    usable.append(d)

  # Per-donor selection string and metal-relative vector, computed once; the
  # angle loop would otherwise rebuild each donor's selection string O(N) times.
  sels = [_selection_string(atoms[d].fetch_labels()) for d in usable]
  vecs = [matrix.col(opt_sites_cart[d]) - metal_xyz for d in usable]

  def _bond_record(donor, selection, distance, sym_op):
    s = bond_sigmas.get(donor) if bond_sigmas else None
    return dict(
      atom_selection_1=metal_sel,
      atom_selection_2=selection,
      symmetry_operation=sym_op,
      donor_element=atoms[donor].element.strip(),
      distance_ideal=distance,
      sigma=s if s is not None else bond_sigma,
      sigma_origin="hessian" if s is not None else "conservative_fallback",
    )

  bonds = [_bond_record(d, sels[idx], vecs[idx].length(), None)
           for idx, d in enumerate(usable)]
  # Symmetry-image donors: restrain the metal to symmetry_operation * ASU-parent
  # (angles cannot cross symmetry, so these get bonds only).
  for d, labels, op_xyz in sym_bonds:
    distance = (matrix.col(opt_sites_cart[d]) - metal_xyz).length()
    bonds.append(
      _bond_record(d, _selection_string(labels), distance, op_xyz))

  angles = []
  if emit_angles:
    for a in range(len(usable)):
      for b in range(a + 1, len(usable)):
        key = frozenset((usable[a], usable[b]))
        if angle_sigmas is None:            # whole-Hessian conservative fallback
          s, origin = angle_sigma, "conservative_fallback"
        elif key in angle_sigmas:           # trusted Seminario or seminario_fallback
          s = angle_sigmas[key]
          origin = angle_origins.get(key, "hessian") if angle_origins else "hessian"
        else:                               # contaminated, no site fallback:
          s, origin = angle_sigma, "conservative_fallback"   # conservative default
        angles.append(dict(
          atom_selection_1=sels[a],
          atom_selection_2=metal_sel,
          atom_selection_3=sels[b],
          angle_ideal=vecs[a].angle(vecs[b], deg=True),
          sigma=s,
          sigma_origin=origin,
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
    ]
    if b.get('symmetry_operation'):
      lines.append(f"        symmetry_operation = {b['symmetry_operation']}")
    lines += [
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
