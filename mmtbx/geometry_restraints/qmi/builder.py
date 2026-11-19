"""QMRestraintsBuilder: the per-site QM metal-restraint pipeline (Phase 1).

Ties the pieces together for one model: extract a capped cluster per metal
site with the ``endo_exo`` builder (in-memory), fix the electronic state, run a
constrained GFN2-xTB optimisation, and emit ``geometry_restraints.edits``.

The electronic state comes from the educated-guess engine (:mod:`electronic`,
plan Sec 6) using the coordination geometry from the CShM classifier
(:mod:`shape`, Sec 5C), or from a user ``metal_state`` override (authoritative).
The total cluster charge is derived from the metal oxidation state with
``quantum_interface.electrons`` (which walks the bond graph + per-element
valences), so deprotonated/charged ligands are accounted for.  Sigmas are the
modified-Seminario curvature sigmas (Sec 5F), with a conservative M-L fallback.
"""

from __future__ import absolute_import, division, print_function

import os

from libtbx import group_args
from libtbx.utils import null_out, Sorry

from mmtbx.geometry_restraints.endo_exo.builder import (
  QMRegionBuilder, default_params)
from scitbx import matrix
from scitbx.array_family import flex
from scitbx.linalg import eigensystem
from cctbx.eltbx import tiny_pse
from iotbx.pdb import common_residue_names_get_class

from mmtbx.geometry_restraints import quantum_interface
from mmtbx.geometry_restraints.qmi import metals as qmi_metals
from mmtbx.geometry_restraints.qmi import donor_perception
from mmtbx.geometry_restraints.qmi import donors, electronic, shape, util
from mmtbx.geometry_restraints.qmi import (
  bvs, provenance, qm, restraints, sanity, sigma)

# Angle-gate coplanarity cutoff (plan Sec 7): lambda_min(sum u u^T)/n below this
# means the bonded-donor directions are (near-)coplanar. 0 for an ideal planar
# shell, 0.2 for a square pyramid (the tightest 3D-spanning shell), ~0.33 for
# tetrahedral/octahedral, so 0.1 separates planar from 3D with margin.
_MIN_OUT_OF_PLANE_COVERAGE = 0.1


def _directions_span_3d(unit_dirs):
  """Whether the metal->donor unit vectors span 3D (so the M-L bonds triangulate
  the metal). The metal-position bond stiffness is ``sum_i u_i u_i^T``; its
  smallest eigenvalue is the softest-mode stiffness, ~0 for a coplanar/collinear
  shell (a soft out-of-plane mode the bonds cannot hold). Normalised by the count
  so the cutoff is coordination-number independent."""
  acc = [0.0] * 9
  for u in unit_dirs:
    for i in range(3):
      for j in range(3):
        acc[3 * i + j] += u[i] * u[j]
  m = flex.double(acc)
  m.reshape(flex.grid(3, 3))
  return (min(eigensystem.real_symmetric(m).values()) / len(unit_dirs)
          >= _MIN_OUT_OF_PLANE_COVERAGE)


def _endo_exo_params(metal_element=None, selection=None, radius=4.0):
  """Build an endo_exo parameter object."""
  params = default_params()
  params.buffer.radius = radius
  if selection:
    params.selection = list(selection)
  if metal_element:
    params.element_filter = list(metal_element)
  return params


class QMRestraintsBuilder(object):
  """Generate QM metal restraints for one model (Phase-1 pipeline).

  Parameters
  ----------
  metal_element : iterable of str or None
      Restrict the metal scan to these element(s); default (None) scans every
      metal endo_exo finds, and each is handled if the electronic engine can
      guess its state (``electronic.OXIDATION_PRIORS``) or the user overrides it.
  selection : iterable of str or None
      Explicit endo_exo seed selection string(s); overrides ``metal_element``.
  radius : float
      Seed-inclusion radius (A) for the endo_exo region (its ``buffer.radius``);
      the final capped cluster can extend beyond it.
  max_depth : int
      endo_exo BFS depth.
  solvent_model : str
      ALPB named solvent for xTB (default ``'ether'``); empty disables.
  temperature : float
      Temperature (K) for the modified-Seminario force-constant to sigma
      conversion (plan Sec 5F); the global sigma tightness scalar.
  metal_states : iterable of group_args or None
      Per-site electronic-state overrides. Each entry has a required
      ``.selection`` (an atom-selection string identifying the metal site(s)),
      plus ``.oxidation_state`` and ``.multiplicity``. Authoritative (plan
      Sec 6): a matched override replaces the engine guess, and a full override
      lets a metal the engine cannot guess be run. Each field falls back
      independently to the engine when left unset; selections must be disjoint
      (an unmatched site just uses the engine guess).
  work_dir : str or None
      Directory the xTB jobs run in (xTB writes many scratch/output files into
      the cwd, so they are collected here rather than the current directory).
      Created if absent; ``None`` keeps them in the current directory.
  predict_only : bool
      Run only the electronic-state prediction per site and stop (no QM, no
      restraints).
  states : str
      Which states to run per site: ``'predicted'`` (the predictor commits one),
      ``'reasonable'`` (scan the accessible oxidations x the geometry-allowed
      spins), or ``'all'`` (scan the accessible oxidations x both spin extremes,
      geometry-blind). Both scans run each state in its own subdirectory.
  nproc : int
      Threads for each xTB optimisation and Hessian.
  perturb_attempts : int
      Perturb-and-reoptimise attempts to escape an imaginary-mode saddle before
      the conservative fallback (plan Sec 5F); each attempt uses a larger kick.
  use_bvs : bool
      Commit the BVS most-likely oxidation instead of the predictor's (spin is
      re-derived for it); a user metal_state oxidation still wins. Applies to the
      committing (states='predicted') path only.
  freeze_waters : str
      Which relax-zone water O to freeze during the optimisation (the H stay
      free): ``'all'`` (every water), ``'peripheral'`` (only non-coordinating
      waters; a metal-coordinating water is held by the metal bond and relaxes,
      so its M-L distance is an optimised restraint), or ``'none'``.
  freeze_far_carboxylate_o : bool
      Freeze the non-coordinating O of a monodentate Asp/Glu carboxylate during
      the optimisation, so it cannot swing in to a spurious bidentate.
  angle_restraints : str
      Donor-metal-donor angle-restraint policy (plan Sec 7): ``'auto'`` emits
      them only for geometries the M-L distances cannot hold (angularly-ambiguous
      isomer at the same distances: everything but a regular high-CN polyhedron);
      ``'always'`` / ``'never'`` are unconditional.
  coord_cutoff : callable or None
      Shared metal-donor coordination cutoff (see
      :func:`util.make_coordination_cutoff`), threaded to donor detection and
      qm-region freezing. None uses the flat fixed cutoff.
  logger : file-like or None
  """

  def __init__(self, metal_element=None, selection=None, radius=4.0,
               max_depth=3, solvent_model=qm.DEFAULT_XTB_SOLVENT,
               temperature=sigma.DEFAULT_TEMPERATURE_K,
               metal_states=None, work_dir="qm_work_dir",
               predict_only=False, states="predicted", nproc=1,
               perturb_attempts=1, use_bvs=False, freeze_waters="peripheral",
               freeze_far_carboxylate_o=False, optimise_all=False,
               angle_restraints="auto", coord_cutoff=None,
               engine="xtb", orca_method=None, orca_basis_set=None,
               orca_solvent_model=None, orca_optimisation="looser",
               orca_hessian="xtb",
               orca_optimise_hydrogens="free_atoms", orca_exe=None, logger=None):
    self.metal_element = list(metal_element) if metal_element else None
    self.selection = list(selection) if selection else None
    self.radius = radius
    self.max_depth = max_depth
    self.solvent_model = solvent_model
    self.temperature = temperature
    self.metal_states = list(metal_states) if metal_states else []
    self.work_dir = work_dir
    self.predict_only = predict_only
    # predict_only runs no QM, so a state scan is meaningless: force the single
    # committing path (the program's validate() rejects the combo up front too).
    self.states = "predicted" if predict_only else states   # 'predicted' | 'reasonable' | 'all'
    self.nproc = nproc
    self.perturb_attempts = perturb_attempts
    self.use_bvs = use_bvs
    self.freeze_waters = freeze_waters
    self.freeze_far_carboxylate_o = freeze_far_carboxylate_o
    self.optimise_all = optimise_all   # True => no $fix, whole cluster relaxes
    self.angle_restraints = angle_restraints   # 'auto' | 'always' | 'never'
    self.coord_cutoff = coord_cutoff   # shared metal-donor cutoff callable, or None
    self.engine = engine               # 'xtb' | 'orca'
    # One QM-backend bundle threaded to every run_constrained_opt / run_hessian
    # call. For xTB the orca_* arguments are ignored; for ORCA the unset ones
    # fall back to the r2SCAN-3c / CPCM(diethylether) defaults.
    if engine == "orca":
      # One user knob (orca_optimisation) drives ORCA's two mechanisms: 'looser'
      # requests a plain Opt and relaxes its tolerances via an explicit %geom
      # block; every other value is a built-in ORCA convergence keyword (no
      # %geom override).
      if orca_optimisation == "looser":
        opt_keyword, geom_conv = "Opt", "looser"
      else:
        opt_keyword, geom_conv = orca_optimisation, None
      self._engine_opts = qm.engine_options(
        engine="orca", method=orca_method, basis_set=orca_basis_set,
        solvent_model=orca_solvent_model, optimisation=opt_keyword,
        convergence=geom_conv, hessian=orca_hessian,
        optimise_hydrogens=orca_optimise_hydrogens, exe=orca_exe)
    else:
      self._engine_opts = qm.engine_options(
        engine="xtb", solvent_model=solvent_model)
    self.logger = logger if logger is not None else null_out()
    self._model_stem = "qmi_model"
    self._results = []
    self._provenance = []   # global provenance comment lines (set by run())

  def run(self, model, model_name="qmi_model", provenance=None):
    """Run the pipeline on *model*; return a list of per-site result objects.

    Each result is a ``libtbx.group_args`` with fields ``site_id``,
    ``metal_element``, ``oxidation_state``, ``multiplicity``, ``charge``,
    ``n_donors``, ``restraints`` (the :func:`restraints.metal_ligand_restraints`
    record), ``edits`` (the PHIL string), and ``skipped`` (None, or a reason).

    *provenance* is an optional list of global comment lines (plan Sec 9.3,
    from :func:`provenance.global_block`) prepended to every emitted site's
    ``.edits`` header.
    """
    self._provenance = list(provenance) if provenance else []
    self._model_stem = os.path.splitext(os.path.basename(model_name))[0]
    write_regions = self.work_dir is not None and not self.predict_only
    if write_regions:
      # qm_work_dir holds the input (protonated) model plus one subdirectory per
      # site (qm_work_dir/<stem>_<chain>_<resseq>) holding that site's QM region
      # (region.pdb), xTB files, and optimised geometry. The .edits files still
      # go to the current directory.
      self.work_dir = os.path.abspath(self.work_dir)
      if not os.path.isdir(self.work_dir):
        os.makedirs(self.work_dir)
      model.get_hierarchy().write_pdb_file(
        crystal_symmetry=model.crystal_symmetry(),
        file_name=os.path.join(self.work_dir,
                               f"{self._model_stem}_protonated.pdb"))
    else:
      self.work_dir = None   # predict_only (or no work dir) writes no QM files

    params = _endo_exo_params(
      metal_element=self.metal_element, selection=self.selection,
      radius=self.radius)
    region_builder = QMRegionBuilder(params, logger=self.logger)
    seed_results = region_builder.run(model)

    original_chain_ids = [c.id.strip()
                          for c in model.get_hierarchy().chains()]

    self._results = []
    for i, seed in enumerate(seed_results, start=1):
      self._results.extend(
        self._process_site(i, seed, original_chain_ids))
    return self._results

  def get_results(self):
    return self._results

  def combined_edits(self, header=None):
    """One ``geometry_restraints.edits`` PHIL string covering every emitted
    site's bonds + angles (sites that were skipped contribute nothing).

    The header is the global provenance block, each emitted site's provenance
    block, then *header* (a summary line)."""
    bonds, angles, site_lines = [], [], []
    for r in self._results:
      if r.restraints is None:
        continue
      bonds += r.restraints.bonds
      angles += r.restraints.angles
      site_lines += r.provenance_lines
    lines = self._provenance + site_lines + ([header] if header else [])
    return restraints.as_edits_phil(
      group_args(bonds=bonds, angles=angles),
      header="\n".join(lines) if lines else None)

  # ------------------------------------------------------------------

  def _lookup_override(self, submodel, metal_iseq):
    """The user ``metal_state`` override whose ``selection`` picks this site's
    metal, or None. Each override's selection is tested against the sub-model
    (which keeps the original chain/resseq labelling); a block that sets
    neither value is ignored. Raises if two selections both match the metal.
    """
    match = None
    for scope in self.metal_states:
      if not scope.selection:
        continue
      if scope.oxidation_state is None and scope.multiplicity is None:
        continue  # nothing forced
      try:
        picked = submodel.selection(scope.selection)[metal_iseq]
      except Exception as exc:
        raise Sorry(f"Invalid metal_state.selection "
                    f"'{scope.selection}': {exc}")
      if picked:
        if match is not None:
          raise Sorry(f"Two metal_state blocks both match the metal picked by "
                      f"'{match.selection}' and '{scope.selection}'; "
                      f"selections must be disjoint.")
        match = scope
    return match

  def _resolve_state(self, metal_element, override, geometry, cn, summary,
                     bvs_est=None):
    """Combine a user ``metal_state`` *override* (or None) with the committing
    electronic predictor into a single (oxidation_state, multiplicity) state,
    given the CShM *geometry* class, *cn*, and the donor *summary* (from
    :func:`donors.summarize`, or None).

    A user override is authoritative for whichever field it sets; the rest come
    from ``electronic.predict_oxidation`` / ``predict_spin`` on refinement-robust
    features only (donor summary + geometry + CN, never bond lengths). qmi ALWAYS
    commits: the only skip is a metal outside the vetted priors with no override
    (``predict_oxidation`` returns None), which then needs a metal_state block.

    With ``use_bvs`` and a *bvs_est*, the BVS most-likely oxidation replaces the
    predictor's (a user oxidation still wins); the spin is re-derived for it.

    Returns a group_args; ``skipped`` is a reason string when the state cannot be
    committed (fields then None), else None with the committed state,
    ``state_source``, ``low_spin``, ``confidence``, and the per-part reasons.
    """
    def skip(reason):
      return group_args(oxidation_state=None, multiplicity=None, low_spin=None,
                        state_source=None, confidence=None,
                        oxidation_reason=None, spin_reason=None, skipped=reason)

    user_ox = override.oxidation_state if override is not None else None
    user_mult = override.multiplicity if override is not None else None

    # Oxidation: a user value wins; else BVS (when use_bvs) or the predictor.
    if user_ox is not None:
      ox, ox_conf, ox_reason = user_ox, electronic.CONF_HIGH, "user-supplied"
    elif self.use_bvs and bvs_est is not None and bvs_est.best is not None:
      ox, ox_conf, ox_reason = (bvs_est.best, electronic.CONF_MEDIUM,
                                "BVS most-likely")
    else:
      pred = electronic.predict_oxidation(metal_element, geometry, cn, summary)
      if pred is None:
        return skip(f"{metal_element} has no electronic prior; add a "
                    f"metal_state block with oxidation_state and multiplicity")
      ox, ox_conf, ox_reason = pred.oxidation, pred.confidence, pred.reason

    # Multiplicity: a user value wins, else the ligand-field spin for that
    # oxidation (heme-aware; high-spin off-heme).
    if user_mult is not None:
      mult, low_spin = user_mult, None
      spin_conf, spin_reason = electronic.CONF_HIGH, "user-supplied"
    else:
      try:
        sp = electronic.predict_spin(metal_element, ox, geometry, summary)
        d = electronic.d_electron_count(metal_element, ox)
        unpaired = electronic.unpaired_electrons(geometry, d,
                                                 low_spin=sp.low_spin)
      except (KeyError, ValueError):
        return skip(f"{metal_element}{ox:+d} gives an impossible d-count; check "
                    f"the oxidation_state override")
      mult, low_spin = unpaired + 1, sp.low_spin
      spin_conf, spin_reason = sp.confidence, sp.reason

    if (user_ox is None and self.use_bvs and bvs_est is not None
        and bvs_est.best is not None):
      source = "bvs"           # BVS set the oxidation (matches the branch above)
    elif user_ox is not None or user_mult is not None:
      source = "user-supplied"
    else:
      source = "guessed"

    # An empty donor summary is not the same as an inconclusive one. When
    # classify() placed no donor into any motif (summary.motif_evidence False),
    # the predictor decided with NO donor information and its reason string
    # ("no ferric signature") reads as a finding rather than as an absence. Cap
    # the confidence and say so. A user-supplied state is unaffected: it did not
    # come from the donors. See docs/donor_classification_generalisation.md.
    confidence = electronic.min_confidence(ox_conf, spin_conf)
    # Only when the committed oxidation actually came FROM the donor summary. A
    # user override, a BVS-sourced state, a single-accessible-state metal (Zn,
    # Mg) and Cu (decided on CN/geometry) are all unaffected by empty donor
    # evidence, so capping them would contradict their own reason strings.
    from_donors = (source == "guessed"
                   and electronic.oxidation_uses_donors(metal_element))
    # Evidence can come from either pass: a motif the residue names placed, or an
    # anionic group only the bond graph could see. A site with neither decided on
    # nothing. The bond-graph term counts only where the element's rules actually
    # read q_hard: Mn keys on motifs alone, so a recovered charge tells it
    # nothing and must not lift its cap.
    informed = (getattr(summary, "motif_evidence", True)
                or (getattr(summary, "n_anionic_unplaced", 0) > 0
                    and electronic.oxidation_uses_donor_charge(metal_element)))
    blind = from_donors and summary is not None and not informed
    if blind:
      confidence = electronic.CONF_LOW
      ox_reason = f"{ox_reason}; NO donor was classified, so this is a default"

    return group_args(
      oxidation_state=ox, multiplicity=mult, low_spin=low_spin,
      state_source=source, confidence=confidence,
      oxidation_reason=ox_reason, spin_reason=spin_reason, skipped=None)

  def _donor_perception(self, submodel, metal_iseq, metal_element, cap_iseqs,
                        cap_original_elements):
    """Bond-graph perception of the cluster, or None when it is unavailable.

    Built on a COPY of the sub-model, deliberately.  Interpreting a model to get
    its bond graph also normalises it (symmetric side chains are renamed), and
    the real sub-model must reach ``region.pdb`` and the frozen partition exactly
    as endo_exo produced it, so the perception is not allowed to perturb
    anything the QM path sees.  The extra interpretation costs a few tenths of a
    second per site.

    The metal needs a charge before ``electrons`` will run (Fe/Mn/Co/Ni are
    absent from its default table), which would be circular here since the
    oxidation state is what we are trying to predict.  It is not: with
    ``link_metals`` off the metal is an isolated node in the graph, so no donor
    charge depends on the value.  A placeholder is passed and only the donors are
    read back.

    Never raises: a cluster the monomer library cannot interpret (an unknown
    ligand) leaves the predictor on the residue-name path it used before, rather
    than turning a prediction into a crash.
    """
    try:
      probe = submodel.deep_copy()
      qm.ensure_restraints(probe, cap_iseqs=cap_iseqs,
                           cap_original_elements=cap_original_elements)
      dist = quantum_interface.electrons(
        probe, return_atom_valences=True, log=null_out(),
        specific_atom_charges=[group_args(
          atom_selection=f"element {metal_element}", charge=2)])
      return donor_perception.from_electron_distribution(
        probe.get_hierarchy().atoms(), dist, metal_iseqs=[metal_iseq])
    except Exception as exc:                  # Sorry, AssertionError, KeyError
      # electrons.py raises bare asserts, whose message is empty, so take the
      # type when there is no text: indexing splitlines() on "" raises inside
      # the handler and destroys the failure this exists to absorb.
      detail = (str(exc) or type(exc).__name__).splitlines()[0][:70]
      print(f"  donor charge unavailable ({type(exc).__name__}: {detail}); "
            f"residue names only", file=self.logger)
      return None

  def _donor_summary(self, atoms, donor_iseqs, perception=None):
    """Classify the first-shell donors into the coordination summary the
    committing predictor reads (:func:`donors.summarize`): heme (FeN4) topology
    plus donor hardness, from residue/element identity only (no bond lengths).
    None when there are no donors.

    *perception* adds each donor's bond-graph resonance-group charge, which
    :func:`donors.summarize` uses only where the residue names ran out. A group
    on a ring the Lewis walk never solved carries no charge here: its value is
    where that walk stopped, not what it concluded."""
    descs = []
    for i in donor_iseqs:
      a = atoms[i]
      c = donors.classify(a.element, a.parent().resname)
      gid = gq = None
      if perception is not None and perception.group_is_resolved(i):
        gid = perception.resonance_group(i)
        gq = perception.group_charge(i)
      descs.append(group_args(element=c.element, kind=c.kind,
                              charge=c.charge, strong_axial=c.strong_axial,
                              res_id=util.residue_key(a),
                              group_id=gid, group_charge=gq))
    return donors.summarize(descs) if descs else None

  def _bvs_estimate(self, atoms, metal_iseq, metal_element, donor_iseqs):
    """BVS estimate of the most self-consistent oxidation over the metal's
    accessible states, from the observed first-shell distances (plan Sec 6.4).
    Returns None when the metal has no accessible states or no donor matched a
    vendored parameter. Advisory cross-check; never changes the committed state.
    """
    ox_states = electronic.accessible_oxidation_states(metal_element)
    if not ox_states:
      return None
    m = matrix.col(atoms[metal_iseq].xyz)
    contacts = [(qmi_metals.normalize_element(atoms[i].element),
                 (matrix.col(atoms[i].xyz) - m).length())
                for i in donor_iseqs]
    return bvs.estimate_oxidation(metal_element, tuple(ox_states), contacts)

  def _process_site(self, index, seed, original_chain_ids):
    """Run the QM-restraint pipeline for one endo_exo seed; return a LIST of
    per-site results (one normally, one per state in a reasonable/all scan)."""
    submodel = seed.model
    seed_iseqs = seed.seed_iseqs
    cap_iseqs = seed.cap_iseqs
    cap_original_elements = seed.cap_original_elements
    cap_anchor_iseqs = seed.cap_anchor_iseqs
    sym_image_provenance = seed.sym_image_provenance
    atoms = submodel.get_hierarchy().atoms()

    if not seed_iseqs:
      return [self._skip_result(index, seed,
                                "no seed atom returned by endo_exo")]
    # Phase-1 mononuclear gate: exactly one metal in the cluster.
    n_metals = qmi_metals.count_metals(submodel)
    if n_metals != 1:
      return [self._skip_result(
        index, seed, f"polynuclear_out_of_scope ({n_metals} metals in cluster)")]

    metal_iseq = seed_iseqs[0]
    metal_atom = atoms[metal_iseq]
    metal_element = qmi_metals.normalize_element(metal_atom.element)
    _rg = metal_atom.parent().parent()
    chain_id = _rg.parent().id.strip()
    resseq = _rg.resseq.strip()
    site_tag = f"{metal_element}_{chain_id}{resseq}"
    # First shell from the crystal geometry, then the CShM geometry class; both
    # feed the electronic guess (the class sets the ligand-field spin).
    donor_iseqs = restraints.detect_first_shell(
      submodel, metal_iseq, coord_cutoff=self.coord_cutoff)
    geom = shape.classify_geometry(atoms, metal_iseq, donor_iseqs)
    # Angle-restraint gate (plan Sec 7): the metal->donor directions of the donors
    # that will actually get a bond. A symmetry-image donor is bonded to its ASU
    # parent (via symmetry_operation) when endo_exo supplied provenance; one
    # without it is skipped and cannot triangulate. _emit_angles then tests
    # whether those bonds pin the metal in 3D.
    metal_xyz = matrix.col(atoms[metal_iseq].xyz)
    bond_dirs = []
    for d in donor_iseqs:
      if (restraints.is_symmetry_image(atoms[d], original_chain_ids)
          and d not in sym_image_provenance):
        continue                       # skipped symmetry donor: no bond emitted
      v = matrix.col(atoms[d].xyz) - metal_xyz
      if v.length() > 0:
        bond_dirs.append(v.normalize())
    emit_angles = self._emit_angles(bond_dirs)
    # BVS cross-check (plan Sec 6.4, advisory): estimate the most self-consistent
    # oxidation over the accessible states, from the OBSERVED first-shell
    # distances. Reported before the skip so an ambiguous metal (e.g. Fe) shows
    # its BVS estimate to guide an override; it drives the committed oxidation
    # only under use_bvs, and is trustworthy only at high resolution / undamaged
    # metal.
    bvs_est = self._bvs_estimate(atoms, metal_iseq, metal_element, donor_iseqs)
    if bvs_est is not None:
      vstr = ", ".join(f"{q:+d}:{v:.2f}"
                       for q, v in sorted(bvs_est.valences.items()))
      if bvs_est.best is None:
        print(f"Site {index}: BVS no verdict, parameter coverage differs "
              f"between candidate states (V {vstr})", file=self.logger)
      else:
        print(f"Site {index}: BVS most-likely ox={bvs_est.best:+d} "
              f"(V {vstr}; observed geometry, high-resolution only)",
              file=self.logger)

    # Donor perception: the residue-name classification, plus the resonance-group
    # charge the bond graph already knows and the name whitelist cannot see (a
    # citrate, a phosphate). Runs on every path including predict_only, since a
    # prediction is exactly what it improves.
    perception = self._donor_perception(submodel, metal_iseq, metal_element,
                                        cap_iseqs, cap_original_elements)
    summary = self._donor_summary(atoms, donor_iseqs, perception=perception)

    # The frozen/free partition is state-independent and, with endo_exo's cap
    # anchors, needs no restraints manager, so compute it ONCE per site and share
    # it between the predicted and scan QM paths.  predict_only runs no QM (and
    # forces states='predicted', so the scan branch is never taken there).
    frozen_iseqs = site_dir = None
    if not self.predict_only:
      frozen_iseqs, free_iseqs = qm.determine_frozen_atoms(
        submodel, metal_iseq, cap_iseqs, cap_anchor_iseqs=cap_anchor_iseqs,
        coord_cutoff=self.coord_cutoff,
        freeze_waters=self.freeze_waters,
        freeze_far_carboxylate_o=self.freeze_far_carboxylate_o, log=self.logger)
      if self.optimise_all:
        frozen_iseqs = []   # freeze nothing: the whole cluster relaxes
      elif (self.engine == 'orca'
            and self._engine_opts.optimise_hydrogens == 'free_atoms'):
        # 'all' leaves the base partition (every non-cap H free) as xTB runs it;
        # 'free_atoms' freezes each H whose heavy atom is frozen, sparing the DFT
        # the floppy scaffold-H modes. Recompute the free set so the constraints
        # summary stays accurate.
        frozen_iseqs = self._freeze_fixed_atom_hydrogens(
          submodel, frozen_iseqs, donor_iseqs, cap_iseqs)
        frozen_set = set(frozen_iseqs)
        free_iseqs = [i for i in range(submodel.get_hierarchy().atoms().size())
                      if i not in frozen_set]
      site_dir = self._site_dir(chain_id, resseq)
      submodel.get_hierarchy().write_pdb_file(
        crystal_symmetry=submodel.crystal_symmetry(),
        file_name=os.path.join(site_dir, "region.pdb"))
      self._write_constraints(submodel, frozen_iseqs, free_iseqs, site_dir)

    if self.states in ("reasonable", "all"):
      return self._process_state_scan(
        index, seed, submodel, metal_iseq, donor_iseqs, cap_iseqs,
        cap_original_elements, geom, metal_element, site_tag,
        original_chain_ids, bvs_est, summary, frozen_iseqs, site_dir,
        emit_angles, sym_image_provenance)

    override = self._lookup_override(submodel, metal_iseq)
    coord_n = self._coordination_number(geom)
    state = self._resolve_state(metal_element, override, geom.shape, coord_n,
                                summary, bvs_est=bvs_est)
    if state.skipped:
      return [self._skip_result(index, seed, state.skipped,
                                metal_element=metal_element, site_tag=site_tag)]

    conf_note = (f" [confidence {state.confidence}]"
                 if state.state_source in ("guessed", "bvs") else "")
    print(f"Site {index}: {metal_element}, {state.state_source} "
          f"ox={state.oxidation_state:+d} mult={state.multiplicity} "
          f"(shape={geom.shape or 'unclassified'}, CN={coord_n}){conf_note}",
          file=self.logger)
    if state.state_source in ("guessed", "bvs"):
      prefix = "predicted: " if state.state_source == "guessed" else ""
      print(f"  {prefix}oxidation ({state.oxidation_reason}); "
            f"spin ({state.spin_reason})", file=self.logger)
    # Advisory BVS against the committed state. Reported here rather than at the
    # BVS line above, which runs before the state is resolved. Also covers
    # predict_only, which returns just below.
    if (bvs_est is not None and bvs_est.best is not None
        and bvs_est.best != state.oxidation_state
        and state.state_source != "bvs"):
      why = " (close call)" if bvs_est.close_call else ""
      print(f"  BVS disagrees: observed distances favour "
            f"ox={bvs_est.best:+d}{why}; committed "
            f"ox={state.oxidation_state:+d} stands (BVS is advisory)",
            file=self.logger)

    if self.predict_only:   # stop before the QM optimisation; prediction only
      return [group_args(
        site_id=seed.label, site_tag=site_tag,
        metal_element=metal_element, oxidation_state=state.oxidation_state,
        multiplicity=state.multiplicity, low_spin=state.low_spin, charge=None,
        state_source=state.state_source, confidence=state.confidence,
        oxidation_reason=state.oxidation_reason, spin_reason=state.spin_reason,
        bvs=bvs_est, n_donors=len(donor_iseqs), restraints=None, edits=None,
        sigma_origin=None, skipped=None)]

    # Restraints manager: needed only for charge accounting (the frozen
    # partition above needs no bond graph).
    qm.ensure_restraints(submodel, cap_iseqs=cap_iseqs,
                         cap_original_elements=cap_original_elements)
    charge = self._total_charge(submodel, metal_element,
                                state.oxidation_state, state.multiplicity)
    return [self._run_state(
      index, seed, submodel, metal_iseq, donor_iseqs, frozen_iseqs,
      metal_element, site_tag, original_chain_ids, bvs_est,
      state.oxidation_state, state.multiplicity, charge,
      site_dir, state.state_source, emit_angles=emit_angles,
      sym_image_provenance=sym_image_provenance,
      low_spin=state.low_spin, confidence=state.confidence,
      oxidation_reason=state.oxidation_reason, spin_reason=state.spin_reason,
      donor_summary=summary)]

  def _coordination_number(self, geom):
    """Coordinative coordination number for the electronic predictor: the
    effective vertex count with each symmetric bidentate carboxylate counted as
    its two donor O (its formal denticity), but an anisobidentate far O left out
    as a long semicoordinate contact. The bond-based angle gate uses raw_cn
    instead (every emitted M-L bond)."""
    return geom.cn + geom.n_bidentate_carboxylate

  def _emit_angles(self, bond_dirs):
    """Whether to emit donor-metal-donor angle restraints for this site (plan
    Sec 7). ``always`` / ``never`` are unconditional; ``auto`` emits when the
    emitted M-L bonds leave the metal's position under-determined: <= 2 bonds
    (linear/bent, metal free on a circle), or >= 3 bonds whose directions are
    (near-)coplanar (a trigonal-planar / square-planar shell, where the bonds
    give zero out-of-plane stiffness and the metal can pucker freely). A bond set
    that spans 3D triangulates the metal and needs no angle (a rubredoxin
    Fe(SCys)4 tetrahedron is held by its four bonds alone). *bond_dirs* are the
    metal->donor unit vectors of the donors that actually get a bond (a chelate's
    second O and a symmetry-image donor with provenance both count)."""
    if self.angle_restraints == "always":
      return True
    if self.angle_restraints == "never":
      return False
    if len(bond_dirs) <= 2:
      return True                      # metal on a circle/line: angle is free
    return not _directions_span_3d(bond_dirs)   # coplanar shell: soft out-of-plane

  def _process_state_scan(self, index, seed, submodel, metal_iseq, donor_iseqs,
                          cap_iseqs, cap_original_elements, geom, metal_element,
                          site_tag, original_chain_ids, bvs_est, summary,
                          frozen_iseqs, site_dir, emit_angles,
                          sym_image_provenance):
    """states=reasonable/all scan: report the predictor's guess (informational;
    not used to pick), then run every scanned (oxidation, multiplicity) state,
    each in its own qm_work_dir/<stem>_<chain>_<resseq>/ox<n>_mult<m>/ subdir.
    'all' is geometry-blind (both spin extremes); 'reasonable' keeps only the
    spin(s) the site's ligand field permits. Returns one result per state (the
    predicted state, if it is among those scanned, is tagged)."""
    coord_n = self._coordination_number(geom)
    pred = electronic.predict_state(metal_element, geom.shape, coord_n, summary)
    if pred is not None:
      print(f"Site {index}: {metal_element}, predictor -> "
            f"ox={pred.oxidation_state:+d} mult={pred.multiplicity} "
            f"[confidence {pred.confidence}] "
            f"(shape={geom.shape or 'unclassified'}, CN={coord_n})",
            file=self.logger)
      print(f"  predicted: oxidation ({pred.oxidation_reason}); "
            f"spin ({pred.spin_reason})", file=self.logger)

    # 'reasonable' prunes to the geometry-allowed spins; 'all' is geometry-blind.
    geometry = geom.shape if self.states == "reasonable" else None
    states = electronic.enumerate_states(metal_element, geometry=geometry)
    if not states:
      return [self._skip_result(
        index, seed, f"{metal_element} has no electronic prior; cannot enumerate "
        f"accessible states", metal_element=metal_element, site_tag=site_tag)]

    print(f"Site {index}: {self.states}-states scan of {len(states)} state(s): "
          + ", ".join(f"ox{ox:+d}/mult{m}" for ox, m in states),
          file=self.logger)
    # ensure_restraints for the per-state charge accounting below; the frozen
    # partition and site_dir were computed once in _process_site.
    qm.ensure_restraints(submodel, cap_iseqs=cap_iseqs,
                         cap_original_elements=cap_original_elements)

    results = []
    for ox, mult in states:
      is_pred = (pred is not None and ox == pred.oxidation_state
                 and mult == pred.multiplicity)
      state_dir = None
      if site_dir is not None:
        state_dir = os.path.join(site_dir, f"ox{ox:+d}_mult{mult}")
        if not os.path.isdir(state_dir):
          os.makedirs(state_dir)
      charge = self._total_charge(submodel, metal_element, ox, mult)
      try:
        results.append(self._run_state(
          index, seed, submodel, metal_iseq, donor_iseqs, frozen_iseqs,
          metal_element, site_tag, original_chain_ids, bvs_est, ox, mult, charge,
          state_dir,
          f"{self.states}-states (predicted)" if is_pred else f"{self.states}-states",
          emit_angles=emit_angles, sym_image_provenance=sym_image_provenance,
          confidence=(pred.confidence if is_pred else None),
          oxidation_reason=(pred.oxidation_reason if is_pred else None),
          spin_reason=(pred.spin_reason if is_pred else None),
          donor_summary=summary))
      except Exception as exc:
        # One state's QM failing must not abort the scan: states=all deliberately
        # enumerates SCF-fragile states. Record the skip and run the rest.
        results.append(self._skip_result(
          index, seed, f"ox={ox:+d} mult={mult} QM failed: {exc}",
          metal_element=metal_element, site_tag=site_tag))
    return results

  def _total_charge(self, submodel, metal_element, oxidation, multiplicity):
    """Total cluster charge for xTB ``--chrg``: the whole system given the
    metal's oxidation state, via ``quantum_interface.electrons`` (bond graph +
    per-element valences, so charged/deprotonated ligands count)."""
    sel = f"element {metal_element}"
    sac = [group_args(atom_selection=sel, charge=oxidation)]
    sam = [group_args(atom_selection=sel, multiplicity=multiplicity)]
    return quantum_interface.electrons(
      submodel, specific_atom_charges=sac,
      specific_atom_multiplicities=sam, log=self.logger)

  def _skip_result(self, index, seed, reason, metal_element=None, site_tag=None):
    """A per-site result marking the site skipped (no restraints)."""
    print(f"Site {index} skipped: {reason}", file=self.logger)
    return group_args(
      site_id=seed.label, site_tag=site_tag,
      metal_element=metal_element, oxidation_state=None, multiplicity=None,
      charge=None, n_donors=None, restraints=None, edits=None,
      sigma_origin=None, skipped=reason)

  def _site_dir(self, chain_id, resseq):
    """The per-site subdirectory qm_work_dir/<stem>_<chain>_<resseq> (created),
    or None when there is no work dir."""
    if self.work_dir is None:
      return None
    d = os.path.join(self.work_dir, f"{self._model_stem}_{chain_id}_{resseq}")
    if not os.path.isdir(d):
      os.makedirs(d)
    return d

  def _write_constraints(self, submodel, frozen_iseqs, free_iseqs, site_dir):
    """Write site_dir/constraints.txt (the frozen vs free atom inventory) beside
    region.pdb; a no-op when there is no work dir."""
    if site_dir is None:
      return
    qm.write_constraints_summary(
      submodel, frozen_iseqs, free_iseqs,
      os.path.join(site_dir, "constraints.txt"))

  def _freeze_fixed_atom_hydrogens(self, submodel, frozen_iseqs, donor_iseqs,
                                   cap_iseqs):
    """Return *frozen_iseqs* with the 'free_atoms' rule: a hydrogen is frozen iff
    its heavy atom is frozen.

    From the base partition (every non-cap hydrogen free) this freezes only the
    hydrogens riding a frozen heavy atom, so no free heavy atom drags a frozen H
    behind it (which would distort its X-H bond). Capping hydrogens are left as the
    base partition froze them (together with their anchors, in determine_frozen_atoms);
    coordinating-water hydrogens are freed even though their O is frozen (the O-H
    orientation shapes the M-O interaction)."""
    atoms = submodel.get_atoms()
    def res_key(atom):
      # A value key: atom.parent() hands back a fresh proxy each call, so object
      # identity is not stable; (chain, resseq, icode, altloc, resname) is.
      ag = atom.parent()
      if ag is None:
        return None
      rg = ag.parent()
      chain = rg.parent() if rg is not None else None
      return (None if chain is None else chain.id,
              None if rg is None else rg.resseq,
              None if rg is None else rg.icode,
              ag.altloc, ag.resname)
    water_donor_keys = set()
    for i in donor_iseqs:
      ag = atoms[i].parent()
      if ag is not None and \
         common_residue_names_get_class(ag.resname) == 'common_water':
        water_donor_keys.add(res_key(atoms[i]))
    caps = set(cap_iseqs)
    frozen_heavy = set(i for i in frozen_iseqs
                       if not atoms[i].element_is_hydrogen())
    frozen = set(frozen_iseqs)
    for j, atom in enumerate(atoms):
      if not atom.element_is_hydrogen():
        continue
      if j in caps:
        continue                    # capping H: keep the base partition's frozen state
      if res_key(atom) in water_donor_keys:
        frozen.discard(j)           # coordinating-water H: always free
        continue
      heavy = self._bonded_heavy_iseq(atoms, j)
      if heavy is not None and heavy in frozen_heavy:
        frozen.add(j)               # rides a frozen atom: freeze it
      else:
        frozen.discard(j)           # rides a free atom: free it
    return sorted(frozen)

  def _bonded_heavy_iseq(self, atoms, h_iseq):
    """iseq of the nearest non-hydrogen atom to *h_iseq* (the heavy atom it is
    bonded to), or None if none lies within a covalent X-H distance."""
    hx, hy, hz = atoms[h_iseq].xyz
    best, best_d2 = None, 1.6 * 1.6   # generous X-H bound (covers S-H)
    for k, a in enumerate(atoms):
      if k == h_iseq or a.element_is_hydrogen():
        continue
      ax, ay, az = a.xyz
      d2 = (ax - hx) ** 2 + (ay - hy) ** 2 + (az - hz) ** 2
      if d2 < best_d2:
        best, best_d2 = k, d2
    return best

  def _run_state(self, index, seed, submodel, metal_iseq, donor_iseqs,
                 frozen_iseqs, metal_element, site_tag, original_chain_ids,
                 bvs_est, oxidation, multiplicity, charge, work_dir,
                 state_source, emit_angles=True, sym_image_provenance=None,
                 low_spin=None, confidence=None,
                 oxidation_reason=None, spin_reason=None, donor_summary=None):
    """Run the constrained opt + Hessian + restraint assembly for a single
    (oxidation, multiplicity, charge) state in *work_dir* (writing optimized.pdb
    and qmi.edits there), returning the per-state result group_args."""
    print(f"Site {index}: calculating {metal_element} "
          f"ox={oxidation:+d} mult={multiplicity}", file=self.logger)
    try:
      opt_sites, _qmm = qm.run_constrained_opt(
        submodel, frozen_iseqs, charge=charge, multiplicity=multiplicity,
        engine_opts=self._engine_opts, nproc=self.nproc, preamble="qmi",
        work_dir=work_dir, log=self.logger)
    except Sorry:
      # The default optimiser can stall on a floppy cluster: the inertial engine
      # never drives a flat PES below the gradient threshold, so no converged
      # geometry is written. Retry once with the robustness preset (looser
      # optlevel + Fermi smearing), which reshapes the PES enough to converge
      # while still honouring $fix. A distinct namespace keeps the stalled
      # attempt's files from being read back.
      print(f"Site {index}: constrained opt did not converge; retrying with "
            f"the robust preset", file=self.logger)
      opt_sites, _qmm = qm.run_constrained_opt(
        submodel, frozen_iseqs, charge=charge, multiplicity=multiplicity,
        engine_opts=self._engine_opts, nproc=self.nproc,
        preamble="qmi_robust", robust=True, work_dir=work_dir, log=self.logger)

    # Modified-Seminario sigmas from the Hessian at the optimised geometry, with
    # ONE perturb-and-reoptimise if the inner sphere is not a true minimum (plan
    # Sec 5F). Returns the final (possibly re-optimised) geometry.
    (opt_sites, bond_sigmas, angle_sigmas, angle_origins,
     sigma_origin, n_perturb) = self._sigmas_with_reopt(
      index, submodel, opt_sites, metal_iseq, donor_iseqs, frozen_iseqs, charge,
      multiplicity, work_dir)

    if work_dir is not None:
      # The (final) optimised QM region, beside region.pdb (the input geometry).
      opt_hier = submodel.get_hierarchy().deep_copy()
      opt_hier.atoms().set_xyz(opt_sites)
      opt_hier.write_pdb_file(
        crystal_symmetry=submodel.crystal_symmetry(),
        file_name=os.path.join(work_dir, "optimized.pdb"))

    site_restraints = restraints.metal_ligand_restraints(
      submodel, opt_sites, metal_iseq, donor_iseqs,
      bond_sigmas=bond_sigmas, angle_sigmas=angle_sigmas,
      angle_origins=angle_origins, emit_angles=emit_angles,
      original_chain_ids=original_chain_ids,
      sym_image_provenance=sym_image_provenance, log=self.logger)
    # Sanity bounds (plan Sec 9.4): cross-check each emitted M-donor distance vs
    # the empirical target and bond-valence windows; flag, never alter.
    mode_changes = sanity.coordination_mode_changes(
      submodel, opt_sites, metal_iseq, donor_iseqs, self.coord_cutoff)
    site_sanity = sanity.check_site(
      metal_element, [(b["donor_element"], b["distance_ideal"])
                      for b in site_restraints.bonds],
      mode_changes=mode_changes)
    if (site_sanity.n_atypical or site_sanity.n_implausible
        or site_sanity.n_mode_changed):
      print(f"Site {index}: sanity-bound flags: "
            + "; ".join(site_sanity.details), file=self.logger)
    # Provenance (plan Sec 9.3): the global block (set by run()) plus this site's
    # committed state / source / reasons / origin / sanity, as the .edits header.
    site_lines = provenance.site_block(
      metal_element, site_tag, oxidation, multiplicity, charge, state_source,
      confidence, oxidation_reason, spin_reason, bvs_est, sigma_origin,
      site_restraints.n_skipped_symmetry, n_perturb, sanity=site_sanity,
      donor_summary=donor_summary)
    edits = restraints.as_edits_phil(
      site_restraints, header="\n".join(self._provenance + site_lines))
    if work_dir is not None:
      with open(os.path.join(work_dir, "qmi.edits"), "w") as fh:
        fh.write(edits)

    skip_note = (f" ({site_restraints.n_skipped_symmetry} symmetry donor(s) "
                 f"skipped)" if site_restraints.n_skipped_symmetry else "")
    print(f"Site {index} ox={oxidation:+d} mult={multiplicity}: "
          f"{len(site_restraints.bonds)} bond + "
          f"{len(site_restraints.angles)} angle restraint(s){skip_note}",
          file=self.logger)

    return group_args(
      site_id=seed.label, site_tag=site_tag,
      metal_element=metal_element, oxidation_state=oxidation,
      multiplicity=multiplicity, charge=charge, low_spin=low_spin,
      state_source=state_source, confidence=confidence,
      oxidation_reason=oxidation_reason, spin_reason=spin_reason, bvs=bvs_est,
      n_donors=len(donor_iseqs), restraints=site_restraints, edits=edits,
      sigma_origin=sigma_origin, sanity=site_sanity,
      provenance_lines=site_lines, skipped=None)

  def _sigmas_with_reopt(self, index, submodel, opt_sites, metal_iseq,
                         donor_iseqs, frozen_iseqs, charge, multiplicity,
                         work_dir):
    """Hessian -> modified-Seminario sigmas at *opt_sites*, with up to
    perturb_attempts perturb-and-reoptimise cycles if the inner sphere carries an
    imaginary mode (plan Sec 5F).

    Returns ``(opt_sites, bond_sigmas, angle_sigmas, angle_origins,
    sigma_origin, n_perturb)``: *opt_sites* is the final (possibly re-optimised)
    geometry, *angle_origins* is the per-angle ``sigma_origin`` map parallel to
    *angle_sigmas*, *sigma_origin* is ``"hessian"`` at a true minimum, else
    ``"conservative_fallback"`` (bond/angle sigmas None, so restraints applies the
    conservative M-L defaults), and *n_perturb* is the number of perturbation
    cycles run. A QM or Hessian failure also degrades to the conservative
    fallback."""
    masses = [tiny_pse.table(a.element.strip()).weight()
              for a in submodel.get_atoms()]

    def hessian_sigmas(sites, preamble):
      # A distinct namespace per call: run_hessian guards only on output-file
      # existence, so a shared namespace would let a failed second --hess read
      # the first geometry's stale Hessian (missing file now raises instead).
      hessian = qm.run_hessian(
        submodel, sites, charge=charge, multiplicity=multiplicity,
        engine_opts=self._engine_opts,
        hessian_atoms=[metal_iseq] + list(donor_iseqs),
        nproc=self.nproc, preamble=preamble, work_dir=work_dir,
        log=self.logger)
      sig = sigma.seminario_sigmas(
        hessian, metal_iseq, donor_iseqs, sites,
        temperature=self.temperature, log=self.logger)
      chk = sigma.imaginary_mode_check(
        hessian, metal_iseq, donor_iseqs, masses, log=self.logger)
      return hessian, sig, chk

    n_perturb = 0
    try:
      hessian, sig, chk = hessian_sigmas(opt_sites, "qmi_hess")
      # Perturb along the imaginary mode and re-optimise, up to perturb_attempts
      # times with a LARGER kick each time (a repeated identical kick is a fixed-
      # point no-op). A distinct namespace per attempt keeps the opt/Hessian cache
      # from replaying an earlier attempt.
      for attempt in range(1, self.perturb_attempts + 1):
        if chk.stable:
          break
        n_perturb = attempt
        step = attempt * sigma.DEFAULT_PERTURB_STEP_A
        print(f"Site {index}: inner-sphere imaginary mode; perturb "
              f"({step:.2f} A) and re-optimise "
              f"(attempt {attempt}/{self.perturb_attempts}, plan Sec 5F)",
              file=self.logger)
        perturbed = sigma.imaginary_mode_perturbation(
          hessian, metal_iseq, donor_iseqs, masses, opt_sites, step=step,
          eigen=chk.eigen)
        opt_sites = self._reoptimise(submodel, perturbed, frozen_iseqs, charge,
                                     multiplicity, work_dir, f"qmi_reopt{attempt}")
        hessian, sig, chk = hessian_sigmas(opt_sites, f"qmi_reopt{attempt}_hess")
      if not chk.stable:
        print(f"Site {index}: still imaginary after {self.perturb_attempts} "
              f"perturb-reoptimise attempt(s); conservative fallback sigmas",
              file=self.logger)
        return opt_sites, None, None, None, "conservative_fallback", n_perturb
      return (opt_sites, sig.bond_sigmas, sig.angle_sigmas, sig.angle_origins,
              "hessian", n_perturb)
    except Exception as exc:
      print(f"Site {index}: Hessian/Seminario failed ({exc}); "
            f"conservative fallback sigmas", file=self.logger)
      return opt_sites, None, None, None, "conservative_fallback", n_perturb

  def _reoptimise(self, submodel, start_sites, frozen_iseqs, charge,
                  multiplicity, work_dir, preamble):
    """Constrained re-optimisation from *start_sites* (the perturbed geometry).

    Uses a per-attempt *preamble* so xtb_manager's coordinate cache does not
    replay an earlier optimisation, and restores the sub-model's original
    coordinates afterwards so a states scan re-uses the same submodel
    independently. Runs with ``file_read=False``: the perturbed start geometry
    is not bit-reproducible across runs, so a stale re-optimisation left in a
    persistent work_dir would differ in the low decimals and trip
    ``check_file_read_safe`` instead of being reused."""
    atoms = submodel.get_hierarchy().atoms()
    saved = atoms.extract_xyz()
    atoms.set_xyz(start_sites)
    try:
      opt_sites, _qmm = qm.run_constrained_opt(
        submodel, frozen_iseqs, charge=charge, multiplicity=multiplicity,
        engine_opts=self._engine_opts, nproc=self.nproc,
        preamble=preamble, work_dir=work_dir, file_read=False, log=self.logger)
    finally:
      atoms.set_xyz(saved)
    return opt_sites
