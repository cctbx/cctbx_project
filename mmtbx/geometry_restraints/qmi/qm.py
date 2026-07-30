"""xTB engine integration for the qmi metal-restraint pipeline.

Drives the in-tree ``mmtbx.geometry_restraints.xtb_manager`` *directly* as a
constrained-optimisation runner (plan Sec 5E.1), not through the QMR
``quantum_restraints_manager`` three-restraints-manager pipeline, whose
in-macro-cycle quantum-refinement design (ligand/buffer/model GRMs +
``super_cell`` buffer + ``atom.tmp`` back-mapping) is the wrong fit for the
endo_exo cluster + static ``.edits`` emission this module wants.

Responsibilities:

- ``ensure_restraints``: build the endo_exo sub-model's geometry restraints
  manager (it arrives without one), restoring cap elements across the
  ``process(make_restraints=True)`` call so the monomer library matches.
- ``determine_frozen_atoms``: the frozen/free partition for the QM
  optimisation, ported from the calibration ``prepare_qm_calcs.py``: relax the
  first shell, freeze the protein scaffold and the truncation boundary.
- ``run_constrained_opt``: construct an ``xtb_manager`` from the sub-model's
  atoms, inject the frozen partition via ``set_ligand_atoms`` (the boolean
  optimise-mask), run the constrained GFN2-xTB optimisation, and return the
  optimised geometry.

The freeze partition is injected as a per-atom optimise mask (True = relax,
False = ``$fix``).  ``run_constrained_opt`` drives ``get_opt`` with
``optimise_h=False`` so the mask is authoritative; otherwise
``base_qm_manager._is_atom_for_opt(optimise_h=True)`` force-frees *every*
hydrogen, which would release the cap H atoms that must stay frozen.
"""

from __future__ import absolute_import, division, print_function

import contextlib
import os

from libtbx import group_args
from libtbx.utils import null_out, Sorry
from scitbx import matrix
from iotbx.pdb import common_residue_names_get_class

from mmtbx.geometry_restraints import xtb_manager
from mmtbx.geometry_restraints import orca_manager
from mmtbx.geometry_restraints.qmi import util


# Cofactor (macrocycle) residue names: a cofactor residue touching the free
# zone is relaxed *whole* so its doming/ruffling degrees of freedom survive
# the constrained optimisation.  (heme variants / chlorophylls / cobalamins /
# F430 / siroheme; same set used by the calibration pipeline.)
COFACTOR_RESNAMES = frozenset({
  "HEM", "HEC", "HEA", "HEB", "HEO", "HES", "CYC", "DHE", "FCO",
  "HDM", "HEV", "HE5", "HE6", "HEG",
  "BCL", "BCB", "BPB", "CLA", "CL0", "CL1", "CL2", "CHL",
  "B12", "COB", "COY", "BTQ", "MM5", "CBY",
  "F43", "F430",
  "SRM", "SH0",
})

# Defaults shared with the calibration pipeline.
DEFAULT_FREE_RADIUS = 5.0      # A; atoms within this of the seed relax
# ALPB named solvent nearest the protein-interior dielectric (eps ~ 4); xTB has
# no arbitrary-epsilon option, so the nearest named solvent stands in (plan
# Sec 5E.1).  Passed to the xtb_manager constructor so the opt is solvated.
DEFAULT_XTB_SOLVENT = "ether"
DEFAULT_XTB_METHOD = "GFN2-xTB"
# ORCA backend: r2SCAN-3c is a composite method (def2-mTZVP + geometric
# counterpoise) that reproduces transition-metal geometries at a fraction of a
# hybrid's cost; CPCM(diethylether) matches the xTB ether solvent's dielectric.
DEFAULT_ORCA_METHOD = "r2SCAN-3c"
DEFAULT_ORCA_SOLVENT = "CPCM(diethylether)"
DEFAULT_ORCA_OPT = "LooseOpt"


def engine_options(engine="xtb", method=None, basis_set="",
                   solvent_model=None, optimisation=DEFAULT_ORCA_OPT,
                   hessian="xtb", optimise_hydrogens=False, exe=None):
  """Bundle the QM-backend configuration shared by :func:`run_constrained_opt`
  and :func:`run_hessian`.

  *engine* is ``'xtb'`` or ``'orca'``.  Unset *method* / *solvent_model* fall
  back to the engine's convention (GFN2-xTB / ether for xTB, r2SCAN-3c /
  CPCM(diethylether) for ORCA).  *optimisation* is the ORCA convergence preset
  (``'LooseOpt'`` or ``'Opt'``); *hessian* selects the force-constant source
  (``'xtb'`` computes the Hessian with GFN2-xTB even under the ORCA optimiser,
  ``'partial'`` / ``'full'`` use an ORCA ``NumFreq``); *optimise_hydrogens*
  ``True`` keeps the non-cap hydrogens free in the ORCA optimisation as xTB
  does, ``False`` freezes them too except the coordinating-water hydrogens,
  cutting the DFT cost (the capping hydrogens are frozen either way); *exe*
  overrides the ORCA executable (else the ``PHENIX_ORCA`` environment
  variable)."""
  if method is None:
    method = DEFAULT_XTB_METHOD if engine == "xtb" else DEFAULT_ORCA_METHOD
  if solvent_model is None:
    solvent_model = DEFAULT_XTB_SOLVENT if engine == "xtb" else DEFAULT_ORCA_SOLVENT
  return group_args(
    engine=engine, method=method, basis_set=basis_set or "",
    solvent_model=solvent_model, optimisation=optimisation,
    hessian=hessian, optimise_hydrogens=optimise_hydrogens, exe=exe)


@contextlib.contextmanager
def _orca_exe_env(exe):
  """Point ``PHENIX_ORCA`` at *exe* for the duration of the block (restoring the
  prior value), so orca_manager's ``get_cmd`` finds the executable.  A falsy
  *exe* leaves the environment untouched (the caller relies on ``PHENIX_ORCA``
  already being set)."""
  if not exe:
    yield
    return
  saved = os.environ.get("PHENIX_ORCA")
  os.environ["PHENIX_ORCA"] = exe
  try:
    yield
  finally:
    if saved is None:
      os.environ.pop("PHENIX_ORCA", None)
    else:
      os.environ["PHENIX_ORCA"] = saved


class _orca_manager(orca_manager.orca_manager):
  """orca_manager wired for QMI, kept as a QMI-local subclass so the shared
  manager stays untouched: run in parallel when ``nproc > 1`` and honour
  ``self.optimisation`` (the base always emits a serial LooseOpt job)."""

  def _input_header(self, gradients_only=False):
    header = super()._input_header(gradients_only=gradients_only)
    optimisation = getattr(self, 'optimisation', None)
    if not gradients_only and optimisation and optimisation != 'LooseOpt':
      header = header.replace('LooseOpt', optimisation, 1)
    if self.nproc > 1:
      header = f'%pal nprocs {self.nproc} end\n' + header
    return header

_BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT"}
# Protein residue classes qualifying for the coordinating-backbone rule.
_AA_CLASSES = util.AA_RESIDUE_CLASSES
# Carboxylate O pairs, for the optional monodentate far-O freeze.
_CARBOXYLATE_O = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}


def ensure_restraints(submodel, cap_iseqs=(), cap_original_elements=()):
  """Build *submodel*'s geometry restraints manager if it has none.

  Idempotent: a no-op when restraints already exist.  The materialised
  sub-model from endo_exo arrives without restraints (the parent's bond proxies
  do not apply atom-for-atom to the symmetry-expanded hierarchy), and the bond
  graph is needed downstream for charge accounting.

  endo_exo cap atoms arrive as ``element='H'`` carrying the anchor's name (e.g.
  CB, C2).  Before ``process(make_restraints=True)`` the original element is
  restored, else the monomer library raises "Conflicting scattering type
  symbols" for non-standard residues (it forgives the element/name conflict
  only for the 20 standard amino acids + MSE); afterwards the element is reset
  to ``'H'`` for the QM run.

  Parameters
  ----------
  submodel : mmtbx.model.manager
      The endo_exo sub-model; processed in place.
  cap_iseqs : sequence of int
      Positional indices of cap atoms in *submodel*.
  cap_original_elements : sequence of str
      Element each cap carried before endo_exo replaced it with H, parallel to
      *cap_iseqs*.
  """
  if submodel.get_restraints_manager() is not None:
    return
  atoms = submodel.get_hierarchy().atoms()
  for iseq, orig_el in zip(cap_iseqs, cap_original_elements):
    atoms[iseq].set_element(orig_el)
  try:
    interpretation_params = submodel.get_current_pdb_interpretation_params()
    interpretation_params.pdb_interpretation.restraints_library.mcl = False
    interpretation_params.pdb_interpretation.automatic_linking.link_metals = False
    submodel.process(
      pdb_interpretation_params=interpretation_params,
      make_restraints=True,
    )
  finally:
    for iseq in cap_iseqs:
      atoms[iseq].set_element('H')


def determine_frozen_atoms(submodel, seed_iseq, cap_iseqs, cap_anchor_iseqs=(),
                           free_radius=DEFAULT_FREE_RADIUS,
                           coord_cutoff=None,
                           freeze_waters=False, freeze_far_carboxylate_o=False,
                           log=None):
  """Return ``(frozen_iseqs, free_iseqs)`` for the QM optimisation.

  The first shell must relax while the protein scaffold and the truncation
  boundary stay put, so the *optimised* first-shell geometry is the restraint
  signal.  Ported from the calibration ``prepare_qm_calcs.py``.

  Free atoms:

  - the seed metal and everything within *free_radius* of it;
  - every atom of a cofactor residue (:data:`COFACTOR_RESNAMES`) that has any
    atom in the free zone, freed whole so doming/ruffling survives;
  - a standard-amino-acid backbone atom (N/CA/C/O/OXT) only when it sits within
    ``coord_cutoff(metal_element, atom_element)`` of the metal (a coordinating
    backbone carbonyl O), plus that O's partner carbonyl C and the residue CA
    (so the C=O reorients via the soft CA-C rotation, hinged at the still-frozen
    amide N);
  - every non-cap H/D, regardless of distance (their positions are placed, not
    measured, so the optimiser should settle them).

  Frozen atoms (the returned set):

  - cap atoms (always; they encode the QM-region boundary);
  - each cap's heavy-atom anchor (*cap_anchor_iseqs*, as recorded by endo_exo),
    which pins the severed boundary bond vector at both ends;
  - non-coordinating standard-AA backbone, and anything beyond *free_radius*
    not in a freed cofactor;
  - with *freeze_waters*: the O of every water (a coordinating water has no
    H-bond network in the truncated cluster and would drift; the 2 H stay
    free).
  - with *freeze_far_carboxylate_o*: the non-coordinating O of a MONODENTATE
    Asp/Glu carboxylate (one O inside the coordination cutoff, its sibling
    outside), so the optimiser cannot pull the sibling into a spurious
    bidentate.

  Parameters
  ----------
  submodel : mmtbx.model.manager
      The endo_exo sub-model (restraints manager not required).
  seed_iseq : int
      Positional index of the seed metal atom.
  cap_iseqs : sequence of int
      Positional indices of cap atoms.
  cap_anchor_iseqs : sequence of int
      Positional indices of the heavy atoms the caps are bonded to (their
      anchors), as recorded by endo_exo (``seed['cap_anchor_iseqs']``).
  free_radius : float
  coord_cutoff : callable or None
      ``(metal_element, donor_element) -> max distance (A)``; see
      :func:`util.make_coordination_cutoff`. Default: the flat fixed cutoff.
  freeze_waters : bool
  freeze_far_carboxylate_o : bool
      Freeze the non-coordinating O of a monodentate Asp/Glu carboxylate.
  log : file-like or None

  Returns
  -------
  (set of int, set of int)
      ``(frozen_iseqs, free_iseqs)``, a partition of ``range(n_atoms)``.
  """
  if log is None:
    log = null_out()
  if coord_cutoff is None:
    coord_cutoff = util.make_coordination_cutoff()
  hier = submodel.get_hierarchy()
  atoms = hier.atoms()
  n = atoms.size()

  seed_xyz = matrix.col(atoms[seed_iseq].xyz)
  metal_element = atoms[seed_iseq].element
  near_seed = {
    i for i, atom in enumerate(atoms)
    if (matrix.col(atom.xyz) - seed_xyz).length() <= free_radius
  }
  free_iseqs = set(near_seed)

  cofactor_residues_freed = set()
  for chain in hier.chains():
    for rg in chain.residue_groups():
      for ag in rg.atom_groups():
        resname = ag.resname.strip().upper()
        if resname not in COFACTOR_RESNAMES:
          continue
        ag_iseqs = [a.i_seq for a in ag.atoms()]
        if any(i in near_seed for i in ag_iseqs):
          for i in ag_iseqs:
            free_iseqs.add(i)
          cofactor_residues_freed.add(
            (chain.id.strip(), rg.resseq.strip(), resname))

  coord_unfreezes = []
  carbonyl_unit_iseqs = set()
  carbonyl_unit_log = []
  for i, atom in enumerate(atoms):
    resname = atom.parent().resname
    name = atom.name.strip().upper()
    if (common_residue_names_get_class(resname) in _AA_CLASSES
        and name in _BACKBONE_NAMES):
      d = (matrix.col(atom.xyz) - seed_xyz).length()
      if d > coord_cutoff(metal_element, atom.element):
        free_iseqs.discard(i)
      else:
        ag = atom.parent()
        rg = ag.parent()
        ch = rg.parent()
        coord_unfreezes.append((ch.id.strip(), rg.resseq.strip(), name, d))
        # Free the partner carbonyl C and the residue CA along with the O,
        # so the C=O reorients via the soft CA-C (psi) rotation, hinged at
        # the still-frozen amide N.  The next residue's peptide N stays
        # frozen, keeping the rigid planar peptide bond from distorting in
        # the truncated cluster.
        if name in ("O", "OXT"):
          for sib in ag.atoms():
            if sib.name.strip().upper() in ("C", "CA"):
              carbonyl_unit_iseqs.add(sib.i_seq)
          carbonyl_unit_log.append((ch.id.strip(), rg.resseq.strip()))

  # Applied after the discard pass so the result is independent of atom order
  # (the carbonyl C / CA are backbone atoms usually beyond coord_cutoff, so
  # their own iterations discard them).
  free_iseqs |= carbonyl_unit_iseqs

  # Free every non-cap H/D: positions are placed, not measured.  Caps (also H)
  # are subtracted just below, so this must precede the cap removal.
  h_iseqs = {
    i for i, atom in enumerate(atoms)
    if atom.element_is_hydrogen()
  }
  free_iseqs |= h_iseqs

  cap_set = set(cap_iseqs)
  free_iseqs -= cap_set

  # Freeze the heavy atom each cap is bonded to: a frozen cap alone pins only
  # one end of the severed boundary bond.  endo_exo records each cap's anchor
  # when it cuts the region (cap_anchor_iseqs), so no bond graph is needed here.
  # Applied last so the boundary anchor wins over near-seed / coordinating-
  # carbonyl freeing.
  cap_anchor_set = set(cap_anchor_iseqs)
  anchor_over_carbonyl = cap_anchor_set & carbonyl_unit_iseqs
  free_iseqs -= cap_anchor_set

  water_o_iseqs = set()
  if freeze_waters:
    for i, atom in enumerate(atoms):
      if (common_residue_names_get_class(atom.parent().resname) == "common_water"
          and not atom.element_is_hydrogen()):
        water_o_iseqs.add(i)
    free_iseqs -= water_o_iseqs

  # Optional: freeze the non-coordinating O of a MONODENTATE carboxylate so the
  # optimiser cannot swing it in to a spurious bidentate (a semiempirical over-
  # coordination tendency, worse in a truncated cluster). Only fires when one
  # carboxylate O is inside the coordination cutoff and its sibling is outside;
  # a genuine bidentate (both O inside) has no far O, so it is untouched.
  carboxylate_far_o = []
  if freeze_far_carboxylate_o:
    for ag in hier.atom_groups():
      pair = _CARBOXYLATE_O.get(ag.resname.strip().upper())
      if pair is None:
        continue
      abn = {a.name.strip().upper(): a for a in ag.atoms()}
      o1, o2 = abn.get(pair[0]), abn.get(pair[1])
      if o1 is None or o2 is None:
        continue
      cut = coord_cutoff(metal_element, o1.element)
      d1 = (matrix.col(o1.xyz) - seed_xyz).length()
      d2 = (matrix.col(o2.xyz) - seed_xyz).length()
      far = o2 if d1 <= d2 else o1
      if min(d1, d2) <= cut < max(d1, d2) and far.i_seq in free_iseqs:
        free_iseqs.discard(far.i_seq)
        rg = ag.parent()
        ch = rg.parent()
        carboxylate_far_o.append((ch.id.strip(), rg.resseq.strip(),
                                  ag.resname.strip(), far.name.strip(),
                                  max(d1, d2)))

  frozen_iseqs = set(range(n)) - free_iseqs

  print(f"    free atoms: {len(free_iseqs)}/{n} "
        f"(frozen: {len(frozen_iseqs)})", file=log)
  print(f"    cap atoms frozen: {len(cap_set)}", file=log)
  print(f"    cap-anchor heavy atoms frozen: {len(cap_anchor_set)}", file=log)
  if anchor_over_carbonyl:
    print(f"    note: {len(anchor_over_carbonyl)} carbonyl-unit atom(s) "
          f"re-frozen as cap anchor", file=log)
  print(f"    H/D atoms freed: {len(h_iseqs - cap_set)} "
        f"(of {len(h_iseqs)}; {len(h_iseqs & cap_set)} are caps)", file=log)
  if freeze_waters:
    print(f"    water O frozen (freeze_waters): {len(water_o_iseqs)}",
          file=log)
  if carboxylate_far_o:
    cx = ", ".join(f"{rn} {ch}{rs}/{nm}(d={d:.2f})"
                   for ch, rs, rn, nm, d in carboxylate_far_o)
    print(f"    monodentate carboxylate far O frozen: {cx}", file=log)
  if cofactor_residues_freed:
    cof_str = ", ".join(
      f"{rn} {ch}{rs}"
      for ch, rs, rn in sorted(cofactor_residues_freed))
    print(f"    cofactor residues fully free: {cof_str}", file=log)
  if coord_unfreezes:
    coord_str = ", ".join(
      f"{ch}{rs}/{nm}(d={d:.2f})"
      for ch, rs, nm, d in coord_unfreezes)
    print(f"    backbone donors kept free (coordinating): {coord_str}",
          file=log)
  if carbonyl_unit_log:
    c_str = ", ".join(f"{ch}{rs}"
                      for ch, rs in sorted(set(carbonyl_unit_log)))
    print(f"    carbonyl C+CA freed (partner of coordinating O): {c_str}",
          file=log)

  return frozen_iseqs, free_iseqs


def write_constraints_summary(submodel, frozen_iseqs, free_iseqs, output_path):
  """Write *output_path* (a ``constraints.txt``): the frozen vs free inventory
  for the QM optimisation, one line per atom (index, FROZEN/free, PDB record),
  with a count header. Mirrors the calibration ``prepare_qm_calcs.py``."""
  atoms = submodel.get_hierarchy().atoms()
  frozen_iseqs = set(frozen_iseqs)
  with open(output_path, "w") as fh:
    fh.write("# QM optimization constraint list\n")
    fh.write(f"# Total atoms: {atoms.size()}\n")
    fh.write(f"# Free:   {len(free_iseqs)}\n")
    fh.write(f"# Frozen: {len(frozen_iseqs)}\n\n")
    fh.write(f"{'idx':>4s}  {'state':6s}  PDB record\n")
    fh.write(f"{'-'*4}  {'-'*6}  {'-'*60}\n")
    for i, atom in enumerate(atoms):
      state = "FROZEN" if i in frozen_iseqs else "free"
      fh.write(f"{i:>4d}  {state:6s}  {atom.format_atom_record().rstrip()}\n")


def optimise_mask_from_frozen(n_atoms, frozen_iseqs):
  """Return the per-atom optimise mask ``xtb_manager.set_ligand_atoms`` wants.

  ``mask[i]`` is ``True`` where atom *i* should relax (``i not in
  frozen_iseqs``) and ``False`` where it should be held by ``$fix``.  This
  inverts the frozen set from :func:`determine_frozen_atoms` into the boolean
  selection the manager consumes; the indices are positional into the same
  atoms array (``i_seq == positional index`` for the materialised sub-model).
  """
  frozen = set(frozen_iseqs)
  return [i not in frozen for i in range(n_atoms)]


def run_constrained_opt(submodel, frozen_iseqs, charge, multiplicity,
                        engine_opts=None, nproc=1, preamble="qmi", work_dir=None,
                        robust=False, cleanup=False, file_read=True, log=None):
  """Run a constrained QM optimisation of *submodel* and return the optimised
  geometry.

  Builds the QM manager for ``engine_opts.engine`` (``'xtb'`` -> xtb_manager,
  ``'orca'`` -> orca_manager) from the sub-model's atoms, injects the
  frozen/free partition as the optimise mask (``set_ligand_atoms``), and runs
  ``get_opt`` with ``optimise_h=False`` so the mask, not the manager's
  blanket "free all H" default, decides every atom (the cap H atoms must stay
  frozen).  Both backends freeze the same scaffold: the xTB ``$fix`` block and
  the ORCA ``%geom Constraints {C i C}`` block are driven by the one mask.
  *charge* and *multiplicity* fix the electronic state.

  Parameters
  ----------
  submodel : mmtbx.model.manager
      Capped endo_exo sub-model.
  frozen_iseqs : iterable of int
      Atoms to hold via ``$fix`` (from :func:`determine_frozen_atoms`).
  charge : int
      Total cluster charge passed to xTB ``--chrg`` (the whole system, not
      just the metal oxidation state).
  multiplicity : int
      Spin multiplicity 2S+1; the engine receives ``multiplicity - 1`` unpaired
      electrons (xTB ``--uhf``, ORCA ``* xyz charge mult``).
  engine_opts : libtbx.group_args or None
      Backend configuration from :func:`engine_options` (engine / method /
      basis_set / solvent_model / optimisation / exe).  ``None`` uses the xTB
      defaults.
  nproc : int
  preamble : str
      Namespaces the QM scratch/output files.
  work_dir : str or None
      Directory to run in (the engine writes files into the cwd); ``None`` runs
      in the current directory.
  robust : bool
      Use the xtb_manager robustness preset (looser opt + Fermi smearing) for
      clusters that otherwise stall/crash the SCF (xTB engine only).
  cleanup : bool or str
      Forwarded to ``xtb_manager.get_opt`` (``'all'`` / ``'most'`` / False).
  file_read : bool
      Reuse a matching cached optimisation in this namespace if present
      (``True``). Pass ``False`` to always run fresh, e.g. for a re-optimisation
      whose start geometry is not bit-reproducible across runs, where a stale
      cache would otherwise trip ``check_file_read_safe``.
  log : file-like or None

  Returns
  -------
  (scitbx.array_family.flex.vec3_double, base_qm_manager.base_qm_manager)
      The optimised geometry of **all** atoms (sub-model atom order), and the
      manager instance (for energy / charge / timing access).
  """
  if log is None:
    log = null_out()
  if engine_opts is None:
    engine_opts = engine_options()
  atoms = submodel.get_atoms()
  mask = optimise_mask_from_frozen(len(atoms), frozen_iseqs)

  if engine_opts.engine == 'orca':
    qmm = _orca_manager(
      atoms, engine_opts.method, engine_opts.basis_set,
      engine_opts.solvent_model, charge, multiplicity,
      nproc=nproc, preamble=preamble)
    qmm.optimisation = engine_opts.optimisation
  else:
    qmm = xtb_manager.xtb_manager(
      atoms, engine_opts.method, '',   # basis_set: unused by xTB
      engine_opts.solvent_model, charge, multiplicity,
      nproc=nproc, preamble=preamble)
    qmm.robust = robust
  qmm.program = engine_opts.engine
  qmm.program_goal = 'opt'   # makes get_opt set optimise_ligand=True
  qmm.set_ligand_atoms(mask)

  # ORCA finds its executable through PHENIX_ORCA; the xTB branch passes
  # exe=None so the context manager is a no-op.
  with util.in_directory(work_dir), _orca_exe_env(
      engine_opts.exe if engine_opts.engine == 'orca' else None):
    # optimise_h=False: the injected mask is authoritative; without it the
    # manager force-frees every H, releasing the frozen cap H atoms.
    _opt_subset, full_coords = qmm.get_opt(
      optimise_h=False, cleanup=cleanup, file_read=file_read, log=log)

  return full_coords, qmm


def run_hessian(submodel, sites_cart, charge, multiplicity, engine_opts=None,
                hessian_atoms=None, nproc=1, preamble="qmi_hess", work_dir=None,
                log=None):
  """Compute the Cartesian Hessian at *sites_cart* and return it parsed.

  Evaluated at the constrained-optimisation geometry (plan Sec 5F).  Only the
  metal + first-shell blocks are consumed by :mod:`sigma`, so the frozen
  scaffold's non-stationarity does not enter and no ``$fix`` / constraint is
  applied here.

  With ``engine_opts.hessian == 'xtb'`` (the default, including under the ORCA
  optimiser) the Hessian is a GFN2-xTB ``--hess`` at *sites_cart*: fast, and its
  curvature is a restraint weight, not a geometry.  With ``'partial'`` /
  ``'full'`` it is an ORCA ``NumFreq`` (``'partial'`` displaces only
  *hessian_atoms*, the metal and its donors, keeping the cost tractable; ORCA
  still writes the full 3N x 3N Hessian).

  Parameters
  ----------
  sites_cart : scitbx.array_family.flex.vec3_double
      The geometry the Hessian is computed at (Angstrom, sub-model atom order).
  engine_opts : libtbx.group_args or None
      Backend configuration from :func:`engine_options`; ``None`` uses xTB.
  hessian_atoms : iterable of int or None
      Atoms displaced by an ORCA partial ``NumFreq`` (metal + donor iseqs);
      ignored by the xTB and full-ORCA paths.

  Returns
  -------
  scitbx.array_family.flex.double
      The (3N, 3N) Hessian in Hartree/Bohr^2, sub-model atom order.
  """
  if log is None:
    log = null_out()
  if engine_opts is None:
    engine_opts = engine_options()
  from mmtbx.geometry_restraints.qmi import sigma
  n_atoms = len(submodel.get_atoms())
  if engine_opts.engine == 'orca' and engine_opts.hessian in ('partial', 'full'):
    hess_path = _run_orca_hessian(
      submodel, sites_cart, charge, multiplicity, engine_opts, hessian_atoms,
      nproc=nproc, preamble=preamble, work_dir=work_dir, log=log)
    return sigma.read_orca_hessian(hess_path, n_atoms)
  # xTB Hessian: always GFN2-xTB in ether, even when ORCA drove the opt (the
  # r2SCAN-3c label engine_opts carries then is meaningless to xtb).
  if engine_opts.engine == 'xtb':
    method, solvent_model = engine_opts.method, engine_opts.solvent_model
  else:
    method, solvent_model = DEFAULT_XTB_METHOD, DEFAULT_XTB_SOLVENT
  hess_path = _run_xtb_hessian(
    submodel, sites_cart, charge, multiplicity, method, solvent_model,
    nproc=nproc, preamble=preamble, work_dir=work_dir, log=log)
  return sigma.read_xtb_hessian(hess_path, n_atoms)


def _run_xtb_hessian(submodel, sites_cart, charge, multiplicity, method,
                     solvent_model, nproc=1, preamble="qmi_hess", work_dir=None,
                     log=None):
  """Run ``xtb --hess`` at *sites_cart*; return the ``.hessian`` file path.

  Reuses ``xtb_manager`` for exe discovery and the GFN / solvent / OpenMP
  conventions, but writes its own ``.xyz`` (at *sites_cart*, not the sub-model's
  input coordinates) and a minimal xcontrol (charge/spin/GFN, no ``$opt``)."""
  if log is None:
    log = null_out()
  from libtbx import easy_run
  exe = xtb_manager.get_exe()
  if not exe:
    raise Sorry("xtb executable not found. Set the PHENIX_XTB environment "
                "variable or install xtb in the build.")
  atoms = submodel.get_atoms()
  # A throwaway manager only to reuse its GFN/solvent/namespace/OMP helpers.
  qmm = xtb_manager.xtb_manager(
    atoms, method, '', solvent_model, charge, multiplicity,
    nproc=nproc, preamble=preamble)
  gfn = qmm._get_gfn()
  solvent = qmm._get_solvent()
  threads = qmm._get_threads()
  ns = qmm.get_namespace()
  # ns already ends in '_hess' (the preamble names this the Hessian job), which
  # namespaces these and the .hessian/.charges/... xTB writes; no '.hess' infix.
  xyz_name = f"{ns}.xyz"
  inp_name = f"{ns}.inp"
  log_name = f"{ns}.log"
  hess_name = f"{ns}.hessian"

  inp_lines = [f"# qmi hessian ({preamble})",
               f"$chrg {int(charge)}",
               f"$spin {int(multiplicity) - 1}"]
  if gfn != 'ff':
    inp_lines += ["$gfn", f"   method={gfn}", "$end"]

  xyz_lines = [str(len(atoms)), f"qmi hessian ({preamble})"]
  for atom, xyz in zip(atoms, sites_cart):
    elem = atom.element.strip()
    if not elem:
      raise Sorry(f"atom {atom.quote()} has blank element column; "
                  "cannot write xTB xyz block")
    xyz_lines.append(f"{elem:<2s} {xyz[0]:14.6f} {xyz[1]:14.6f} {xyz[2]:14.6f}")

  gfn_flag = '--gfnff' if gfn == 'ff' else f'--gfn {gfn}'
  solvent_flag = f' --alpb {solvent}' if solvent else ''
  cmd = (f"{exe} {xyz_name} --input {inp_name} --namespace {ns} "
         f"{gfn_flag}{solvent_flag} -P {threads} --hess >& {log_name}")

  with util.in_directory(work_dir):
    with open(inp_name, 'w') as f:
      f.write("\n".join(inp_lines) + "\n")
    with open(xyz_name, 'w') as f:
      f.write("\n".join(xyz_lines) + "\n")
    # A prior run's Hessian in this namespace must not survive an aborted --hess
    # and be read back as this run's output.
    if os.path.exists(hess_name):
      os.remove(hess_name)
    saved = {k: os.environ.get(k) for k in qmm._omp_env()}
    os.environ.update(qmm._omp_env())
    try:
      print(f"  Starting (hessian): {cmd}", file=log)
      easy_run.go(cmd)
    finally:
      for k, v in saved.items():
        if v is None:
          os.environ.pop(k, None)
        else:
          os.environ[k] = v
    # Trust the Hessian only if xtb terminated normally (mirrors the opt path's
    # xtb_manager check); an abort otherwise leaves a stale/partial file.
    qmm._check_log(log_name, log=log)
    if not os.path.exists(hess_name):
      raise Sorry(f"xtb --hess did not produce {hess_name}; check {log_name}")
    return os.path.abspath(hess_name)


def _run_orca_hessian(submodel, sites_cart, charge, multiplicity, engine_opts,
                      hessian_atoms, nproc=1, preamble="qmi_hess", work_dir=None,
                      log=None):
  """Run an ORCA ``NumFreq`` at *sites_cart*; return the ``.hess`` file path.

  ``engine_opts.hessian == 'partial'`` displaces only *hessian_atoms* (metal +
  donors) via ``%freq Partial_Hess {...} end``; ORCA still writes the full
  3N x 3N Hessian, so :func:`sigma.read_orca_hessian` reads it unchanged.  A
  ``'full'`` request displaces every atom (much costlier)."""
  if log is None:
    log = null_out()
  from libtbx import easy_run
  atoms = submodel.get_atoms()
  base = f"orca_{preamble}"
  inp_name = f"{base}.in"
  log_name = f"{base}.log"
  hess_name = f"{base}.hess"

  header = ["%maxcore 2048"]
  if nproc > 1:
    header.append(f"%pal nprocs {int(nproc)} end")
  keywords = " ".join(w for w in (engine_opts.method, engine_opts.basis_set,
                                  engine_opts.solvent_model, "NumFreq") if w)
  lines = header + [f"! {keywords}"]
  if engine_opts.hessian == 'partial' and hessian_atoms:
    sel = " ".join(str(int(i)) for i in sorted(set(hessian_atoms)))
    lines += ["%freq", f"  Partial_Hess {{{sel}}} end", "end"]
  lines.append(f"* xyz {int(charge)} {int(multiplicity)}")
  for atom, xyz in zip(atoms, sites_cart):
    elem = atom.element.strip()
    if not elem:
      raise Sorry(f"atom {atom.quote()} has blank element column; "
                  "cannot write ORCA coordinate block")
    lines.append(f" {elem} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}")
  lines.append("*")

  with util.in_directory(work_dir), _orca_exe_env(engine_opts.exe):
    exe = os.environ.get("PHENIX_ORCA")
    if not exe:
      raise Sorry("ORCA executable not found. Set qm.orca.exe or the "
                  "PHENIX_ORCA environment variable.")
    with open(inp_name, 'w') as f:
      f.write("\n".join(lines) + "\n")
    # A prior run's Hessian must not survive an aborted NumFreq and be read back.
    if os.path.exists(hess_name):
      os.remove(hess_name)
    cmd = f"{exe} {inp_name} >& {log_name}"
    print(f"  Starting (hessian): {cmd}", file=log)
    easy_run.go(cmd)
    if not os.path.exists(hess_name):
      raise Sorry(f"ORCA NumFreq did not produce {hess_name}; check {log_name}")
    return os.path.abspath(hess_name)
