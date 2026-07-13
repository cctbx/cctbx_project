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

import os
from contextlib import contextmanager

from libtbx.utils import null_out, Sorry
from scitbx import matrix

from mmtbx.geometry_restraints import xtb_manager


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
DEFAULT_COORD_CUTOFF = 3.0     # A; a backbone donor inside this coordinates
# ALPB named solvent nearest the protein-interior dielectric (eps ~ 4); xTB has
# no arbitrary-epsilon option, so the nearest named solvent stands in (plan
# Sec 5E.1).  Passed to the xtb_manager constructor so the opt is solvated.
DEFAULT_SOLVENT = "ether"
DEFAULT_METHOD = "GFN2-xTB"

_BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT"}
_STD_AA = {
  "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
  "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


@contextmanager
def _in_directory(path):
  """chdir into *path* for the block, restoring the prior cwd after.  xTB (via
  ``xtb_manager``) writes its scratch / output files into the current working
  directory, so a per-run scratch dir keeps concurrent sites from colliding and
  keeps the source tree clean.  ``path=None`` is a no-op (stay in cwd)."""
  if path is None:
    yield
    return
  saved = os.getcwd()
  os.chdir(path)
  try:
    yield
  finally:
    os.chdir(saved)


def ensure_restraints(submodel, cap_iseqs=(), cap_original_elements=()):
  """Build *submodel*'s geometry restraints manager if it has none.

  Idempotent: a no-op when restraints already exist.  The materialised
  sub-model from endo_exo arrives without restraints (the parent's bond proxies
  do not apply atom-for-atom to the symmetry-expanded hierarchy), and the bond
  graph is needed downstream: by :func:`determine_frozen_atoms` (cap-anchor
  detection) and by charge accounting.

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
    submodel.process(make_restraints=True)
  finally:
    for iseq in cap_iseqs:
      atoms[iseq].set_element('H')


def determine_frozen_atoms(submodel, seed_iseq, cap_iseqs,
                           free_radius=DEFAULT_FREE_RADIUS,
                           coord_cutoff_A=DEFAULT_COORD_CUTOFF,
                           freeze_water_o=False, log=None):
  """Return ``(frozen_iseqs, free_iseqs)`` for the QM optimisation.

  The first shell must relax while the protein scaffold and the truncation
  boundary stay put, so the *optimised* first-shell geometry is the restraint
  signal.  Ported from the calibration ``prepare_qm_calcs.py``.

  Free atoms:

  - the seed metal and everything within *free_radius* of it;
  - every atom of a cofactor residue (:data:`COFACTOR_RESNAMES`) that has any
    atom in the free zone, freed whole so doming/ruffling survives;
  - a standard-amino-acid backbone atom (N/CA/C/O/OXT) only when it sits within
    *coord_cutoff_A* of the metal (a coordinating backbone carbonyl O), plus
    that O's partner carbonyl C and the residue CA (so the C=O reorients via
    the soft CA-C rotation, hinged at the still-frozen amide N);
  - every non-cap H/D, regardless of distance (their positions are placed, not
    measured, so the optimiser should settle them).

  Frozen atoms (the returned set):

  - cap atoms (always; they encode the QM-region boundary);
  - the heavy atom each cap is bonded to (the cap anchor, which pins the boundary
    bond vector at both ends; needs the bond graph, so call
    :func:`ensure_restraints` first);
  - non-coordinating standard-AA backbone, and anything beyond *free_radius*
    not in a freed cofactor;
  - with *freeze_water_o*: the O of every water (a coordinating water has no
    H-bond network in the truncated cluster and would drift; the 2 H stay
    free).

  Parameters
  ----------
  submodel : mmtbx.model.manager
      Sub-model with restraints attached (see :func:`ensure_restraints`).
  seed_iseq : int
      Positional index of the seed metal atom.
  cap_iseqs : sequence of int
      Positional indices of cap atoms.
  free_radius, coord_cutoff_A : float
  freeze_water_o : bool
  log : file-like or None

  Returns
  -------
  (set of int, set of int)
      ``(frozen_iseqs, free_iseqs)``, a partition of ``range(n_atoms)``.
  """
  if log is None:
    log = null_out()
  hier = submodel.get_hierarchy()
  atoms = hier.atoms()
  n = atoms.size()

  seed_xyz = matrix.col(atoms[seed_iseq].xyz)
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

  # Positional index per atom, so a coordinating backbone O can free its
  # partner carbonyl C (located by atom-group membership).
  pos_by_id = {atom.memory_id(): i for i, atom in enumerate(atoms)}
  coord_unfreezes = []
  carbonyl_unit_iseqs = set()
  carbonyl_unit_log = []
  for i, atom in enumerate(atoms):
    resname = atom.parent().resname.strip().upper()
    name = atom.name.strip().upper()
    if resname in _STD_AA and name in _BACKBONE_NAMES:
      d = (matrix.col(atom.xyz) - seed_xyz).length()
      if d > coord_cutoff_A:
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
              carbonyl_unit_iseqs.add(pos_by_id[sib.memory_id()])
          carbonyl_unit_log.append((ch.id.strip(), rg.resseq.strip()))

  # Applied after the discard pass so the result is independent of atom order
  # (the carbonyl C / CA are backbone atoms usually beyond coord_cutoff, so
  # their own iterations discard them).
  free_iseqs |= carbonyl_unit_iseqs

  # Free every non-cap H/D: positions are placed, not measured.  Caps (also H)
  # are subtracted just below, so this must precede the cap removal.
  h_iseqs = {
    i for i, atom in enumerate(atoms)
    if atom.element.strip().upper() in ("H", "D")
  }
  free_iseqs |= h_iseqs

  cap_set = set(cap_iseqs)
  free_iseqs -= cap_set

  # Freeze the heavy atom each cap is bonded to: a frozen cap alone pins only
  # one end of the severed boundary bond.  Each cap is a terminal H, so its
  # single bonded heavy neighbour is the anchor.  Needs the bond graph; if
  # restraints were not built, the anchors are skipped (logged).  Applied last
  # so the boundary anchor wins over near-seed / coordinating-carbonyl freeing.
  cap_anchor_iseqs = set()
  rm = submodel.get_restraints_manager()
  if cap_set and rm is not None:
    bond_proxies, _asu = rm.geometry.get_all_bond_proxies(
      sites_cart=submodel.get_sites_cart())
    for bp in bond_proxies:
      a, b = bp.i_seqs
      if a in cap_set and b not in cap_set:
        cap_anchor_iseqs.add(b)
      elif b in cap_set and a not in cap_set:
        cap_anchor_iseqs.add(a)
  anchor_over_carbonyl = cap_anchor_iseqs & carbonyl_unit_iseqs
  free_iseqs -= cap_anchor_iseqs

  water_o_iseqs = set()
  if freeze_water_o:
    for i, atom in enumerate(atoms):
      if (atom.parent().resname.strip().upper() in ("HOH", "DOD")
          and atom.element.strip().upper() not in ("H", "D")):
        water_o_iseqs.add(i)
    free_iseqs -= water_o_iseqs

  frozen_iseqs = set(range(n)) - free_iseqs

  print(f"    free atoms: {len(free_iseqs)}/{n} "
        f"(frozen: {len(frozen_iseqs)})", file=log)
  print(f"    cap atoms frozen: {len(cap_set)}", file=log)
  if cap_set and rm is None:
    print("    WARNING: no restraints manager; cap-anchor atoms not frozen "
          "(boundary bonds pinned at one end only)", file=log)
  else:
    print(f"    cap-anchor heavy atoms frozen: {len(cap_anchor_iseqs)}",
          file=log)
  if anchor_over_carbonyl:
    print(f"    note: {len(anchor_over_carbonyl)} carbonyl-unit atom(s) "
          f"re-frozen as cap anchor", file=log)
  print(f"    H/D atoms freed: {len(h_iseqs - cap_set)} "
        f"(of {len(h_iseqs)}; {len(h_iseqs & cap_set)} are caps)", file=log)
  if freeze_water_o:
    print(f"    water O frozen (freeze_water_o): {len(water_o_iseqs)}",
          file=log)
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
                        solvent_model=DEFAULT_SOLVENT, method=DEFAULT_METHOD,
                        nproc=1, preamble="qmi", work_dir=None, robust=False,
                        cleanup=False, log=None):
  """Run a constrained GFN2-xTB optimisation of *submodel* and return the
  optimised geometry.

  Builds an ``xtb_manager`` from the sub-model's atoms, injects the
  frozen/free partition as the optimise mask (``set_ligand_atoms``), and runs
  ``get_opt`` with ``optimise_h=False`` so the mask, not the manager's
  blanket "free all H" default, decides every atom (the cap H atoms must stay
  frozen).  *charge* and *multiplicity* fix the electronic state; the engine
  is run spin-restricted with ``$spin = multiplicity - 1``.

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
      Spin multiplicity 2S+1; xTB receives ``--uhf (multiplicity - 1)``.
  solvent_model : str
      ALPB named solvent (default ``'ether'``); empty/None disables solvation.
  method : str
      xTB parametrisation (default ``'GFN2-xTB'``).
  nproc : int
  preamble : str
      Namespaces the xTB scratch/output files.
  work_dir : str or None
      Directory to run in (xTB writes files into the cwd); ``None`` runs in the
      current directory.
  robust : bool
      Use the xtb_manager robustness preset (looser opt + Fermi smearing) for
      clusters that otherwise stall/crash the SCF.
  cleanup : bool or str
      Forwarded to ``xtb_manager.get_opt`` (``'all'`` / ``'most'`` / False).
  log : file-like or None

  Returns
  -------
  (scitbx.array_family.flex.vec3_double, xtb_manager.xtb_manager)
      The optimised geometry of **all** atoms (sub-model atom order), and the
      manager instance (for energy / charge / timing access).
  """
  if log is None:
    log = null_out()
  atoms = submodel.get_atoms()
  mask = optimise_mask_from_frozen(len(atoms), frozen_iseqs)

  qmm = xtb_manager.xtb_manager(
    atoms,
    method,
    '',              # basis_set: unused by xTB
    solvent_model,
    charge,
    multiplicity,
    nproc=nproc,
    preamble=preamble,
  )
  qmm.program = 'xtb'
  qmm.program_goal = 'opt'   # makes get_opt set optimise_ligand=True
  qmm.robust = robust
  qmm.set_ligand_atoms(mask)

  with _in_directory(work_dir):
    # optimise_h=False: the injected mask is authoritative; without it the
    # manager force-frees every H, releasing the frozen cap H atoms.
    _opt_subset, full_coords = qmm.get_opt(
      optimise_h=False, cleanup=cleanup, log=log)

  return full_coords, qmm


def run_hessian(submodel, sites_cart, charge, multiplicity,
                solvent_model=DEFAULT_SOLVENT, method=DEFAULT_METHOD,
                nproc=1, preamble="qmi_hess", work_dir=None, log=None):
  """Run ``xtb --hess`` at *sites_cart* and return the Hessian file path.

  Computes the full Cartesian Hessian (Hartree/Bohr^2) at the given geometry
  (the constrained-optimisation result: plan Sec 5F evaluates the Hessian at
  the optimised geometry).  No ``$fix`` is applied: the metal and donors were
  free during the opt, so their gradient is ~0 and their Hessian blocks are at
  a true local minimum; only those inner-sphere blocks are consumed by
  :mod:`sigma`, so the frozen scaffold's non-stationarity does not enter.

  Reuses ``xtb_manager`` for exe discovery and the GFN / solvent / OpenMP
  conventions, but writes its own ``.xyz`` (at *sites_cart*, not the
  sub-model's input coordinates) and a minimal xcontrol (charge/spin/GFN, no
  ``$opt``).

  Parameters mirror :func:`run_constrained_opt`.  *sites_cart* is the geometry
  the Hessian is computed at (Angstrom, sub-model atom order).

  Returns
  -------
  str
      Absolute path to the ``<namespace>.hessian`` file (read it with
      :func:`mmtbx.geometry_restraints.qmi.sigma.read_xtb_hessian`).
  """
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
  xyz_name = f"{ns}.xyz"
  inp_name = f"{ns}.hess.inp"
  log_name = f"{ns}.hess.log"
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

  with _in_directory(work_dir):
    with open(inp_name, 'w') as f:
      f.write("\n".join(inp_lines) + "\n")
    with open(xyz_name, 'w') as f:
      f.write("\n".join(xyz_lines) + "\n")
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
    if not os.path.exists(hess_name):
      raise Sorry(f"xtb --hess did not produce {hess_name}; check {log_name}")
    return os.path.abspath(hess_name)
