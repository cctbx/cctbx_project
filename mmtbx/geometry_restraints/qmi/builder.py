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
modified-Seminario curvature sigmas (Sec 5F), with a placeholder fallback.
"""

from __future__ import absolute_import, division, print_function

import os

import libtbx.phil
from libtbx import group_args
from libtbx.utils import null_out, Sorry

from mmtbx.programs.endo_exo import master_phil_str as endo_exo_master_phil_str
from mmtbx.geometry_restraints.endo_exo.builder import QMRegionBuilder
from scitbx import matrix

from mmtbx.geometry_restraints import quantum_interface
from mmtbx.geometry_restraints.qmi import metals as qmi_metals
from mmtbx.geometry_restraints.qmi import donors, electronic, shape
from mmtbx.geometry_restraints.qmi import bvs, qm, restraints, sigma


def _default_output_filename(prefix=None, suffix=None):
  """Filename-stem callable for QMRegionBuilder when running in-memory (we set
  ``write_files=False``, but endo_exo still derives a stem per seed)."""
  return f"{prefix or ''}{suffix or ''}"


def _endo_exo_params(metal_element=None, selection=None, radius=4.0,
                    max_depth=3):
  """Build an endo_exo parameter object for an in-memory run (no file output)."""
  params = libtbx.phil.parse(endo_exo_master_phil_str).extract()
  params.write_files = False
  params.buffer.radius = radius
  params.max_search_depth = max_depth
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
  logger : file-like or None
  """

  def __init__(self, metal_element=None, selection=None, radius=4.0,
               max_depth=3, solvent_model=qm.DEFAULT_SOLVENT,
               temperature=sigma.DEFAULT_TEMPERATURE_K,
               metal_states=None, work_dir="qm_work_dir",
               logger=None):
    self.metal_element = list(metal_element) if metal_element else None
    self.selection = list(selection) if selection else None
    self.radius = radius
    self.max_depth = max_depth
    self.solvent_model = solvent_model
    self.temperature = temperature
    self.metal_states = list(metal_states) if metal_states else []
    self.work_dir = work_dir
    self.logger = logger if logger is not None else null_out()
    self._model_stem = "qmi_model"
    self._results = []

  def run(self, model, model_name="qmi_model"):
    """Run the pipeline on *model*; return a list of per-site result objects.

    Each result is a ``libtbx.group_args`` with fields ``site_id``,
    ``metal_element``, ``oxidation_state``, ``multiplicity``, ``charge``,
    ``n_donors``, ``restraints`` (the :func:`restraints.metal_ligand_restraints`
    record), ``edits`` (the PHIL string), and ``skipped`` (None, or a reason).
    """
    self._model_stem = os.path.splitext(os.path.basename(model_name))[0]
    if self.work_dir is not None:
      # Collect the many per-site xTB scratch/output files in one directory
      # (created if absent); the .edits files still go to the current directory.
      self.work_dir = os.path.abspath(self.work_dir)
      if not os.path.isdir(self.work_dir):
        os.makedirs(self.work_dir)

    params = _endo_exo_params(
      metal_element=self.metal_element, selection=self.selection,
      radius=self.radius, max_depth=self.max_depth)
    region_builder = QMRegionBuilder(params, logger=self.logger)
    seed_results = region_builder.run(
      model, model_name=model_name,
      default_output_filename=_default_output_filename)

    original_chain_ids = [c.id.strip()
                          for c in model.get_hierarchy().chains()]

    self._results = []
    for i, seed in enumerate(seed_results, start=1):
      self._results.append(
        self._process_site(i, seed, original_chain_ids))
    return self._results

  def get_results(self):
    return self._results

  def combined_edits(self, header=None):
    """One ``geometry_restraints.edits`` PHIL string covering every emitted
    site's bonds + angles (sites that were skipped contribute nothing)."""
    bonds, angles = [], []
    for r in self._results:
      if r.restraints is None:
        continue
      bonds += r.restraints.bonds
      angles += r.restraints.angles
    return restraints.as_edits_phil(
      group_args(bonds=bonds, angles=angles), header=header)

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

  def _resolve_state(self, metal_element, override, geometry, cn, summary):
    """Combine a user ``metal_state`` *override* (or None) with the committing
    electronic predictor into a single (oxidation_state, multiplicity) state,
    given the CShM *geometry* class, *cn*, and the donor *summary* (from
    :func:`donors.summarize`, or None).

    A user override is authoritative for whichever field it sets; the rest come
    from ``electronic.predict_oxidation`` / ``predict_spin`` on refinement-robust
    features only (donor summary + geometry + CN, never bond lengths). qmi ALWAYS
    commits: the only skip is a metal outside the vetted priors with no override
    (``predict_oxidation`` returns None), which then needs a metal_state block.

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

    # Oxidation: a user value wins, else the committing predictor.
    if user_ox is not None:
      ox, ox_conf, ox_reason = user_ox, electronic.CONF_HIGH, "user-supplied"
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

    user = user_ox is not None or user_mult is not None
    return group_args(
      oxidation_state=ox, multiplicity=mult, low_spin=low_spin,
      state_source="user-supplied" if user else "guessed",
      confidence=electronic.min_confidence(ox_conf, spin_conf),
      oxidation_reason=ox_reason, spin_reason=spin_reason, skipped=None)

  def _donor_summary(self, atoms, donor_iseqs):
    """Classify the first-shell donors into the coordination summary the
    committing predictor reads (:func:`donors.summarize`): heme (FeN4) topology
    plus donor hardness, from residue/element identity only (no bond lengths).
    None when there are no donors."""
    descs = []
    for i in donor_iseqs:
      a = atoms[i]
      ag = a.parent()
      rg = ag.parent()
      res_id = (rg.parent().id.strip(), rg.resseq.strip(), rg.icode.strip(),
                ag.resname.strip())
      c = donors.classify(a.element, ag.resname, a.name)
      descs.append(group_args(element=c.element, kind=c.kind, field=c.field,
                              charge=c.charge, strong_axial=c.strong_axial,
                              res_id=res_id))
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
    contacts = [(atoms[i].element.strip().capitalize(),
                 (matrix.col(atoms[i].xyz) - m).length())
                for i in donor_iseqs]
    return bvs.estimate_oxidation(metal_element, tuple(ox_states), contacts)

  def _process_site(self, index, seed, original_chain_ids):
    """Run the QM-restraint pipeline for one endo_exo seed result."""
    submodel = seed["model"]
    seed_iseqs = seed["seed_iseqs"]
    cap_iseqs = seed["cap_iseqs"]
    cap_original_elements = seed.get("cap_original_elements", ())
    atoms = submodel.get_hierarchy().atoms()

    def _skip(reason, **extra):
      print(f"Site {index} skipped: {reason}", file=self.logger)
      return group_args(
        site_id=seed.get("file_name", f"site_{index}"),
        site_tag=extra.get("site_tag"),
        metal_element=extra.get("metal_element"),
        oxidation_state=None, multiplicity=None, charge=None,
        n_donors=None, restraints=None, edits=None, sigma_origin=None,
        skipped=reason)

    if not seed_iseqs:
      return _skip("no seed atom returned by endo_exo")
    # Phase-1 mononuclear gate: exactly one metal in the cluster.
    n_metals = qmi_metals.count_metals(submodel)
    if n_metals != 1:
      return _skip(f"polynuclear_out_of_scope ({n_metals} metals in cluster)")

    metal_iseq = seed_iseqs[0]
    metal_atom = atoms[metal_iseq]
    metal_element = metal_atom.element.strip().capitalize()
    _rg = metal_atom.parent().parent()
    site_tag = f"{metal_element}_{_rg.parent().id.strip()}{_rg.resseq.strip()}"
    override = self._lookup_override(submodel, metal_iseq)
    # First shell from the crystal geometry, then the CShM geometry class; both
    # feed the electronic guess (the class sets the ligand-field spin).
    donor_iseqs = restraints.detect_first_shell(submodel, metal_iseq)
    geom = shape.classify_geometry(atoms, metal_iseq, donor_iseqs)
    # BVS cross-check (plan Sec 6.4, advisory): estimate the most self-consistent
    # oxidation over the accessible states, from the OBSERVED first-shell
    # distances. Reported before the skip so an ambiguous metal (e.g. Fe) shows
    # its BVS estimate to guide an override; it never overrides the electronic
    # guess, and is trustworthy only at high resolution / undamaged metal.
    bvs_est = self._bvs_estimate(atoms, metal_iseq, metal_element, donor_iseqs)
    if bvs_est is not None:
      vstr = ", ".join(f"{q:+d}:{v:.2f}"
                       for q, v in sorted(bvs_est.valences.items()))
      print(f"Site {index}: BVS most-likely ox={bvs_est.best:+d} "
            f"(V {vstr}; observed geometry, high-resolution only)",
            file=self.logger)

    summary = self._donor_summary(atoms, donor_iseqs)
    state = self._resolve_state(metal_element, override, geom.shape, geom.cn,
                                summary)
    if state.skipped:
      return _skip(state.skipped, metal_element=metal_element)

    conf_note = (f" [confidence {state.confidence}]"
                 if state.state_source == "guessed" else "")
    print(f"Site {index}: {metal_element}, {state.state_source} "
          f"ox={state.oxidation_state:+d} mult={state.multiplicity} "
          f"(shape={geom.shape or 'unclassified'}, CN={geom.cn}){conf_note}",
          file=self.logger)
    if state.state_source == "guessed":
      print(f"  predicted: oxidation ({state.oxidation_reason}); "
            f"spin ({state.spin_reason})", file=self.logger)

    # Restraints manager (needed for charge accounting + cap-anchor freezing).
    qm.ensure_restraints(submodel, cap_iseqs=cap_iseqs,
                         cap_original_elements=cap_original_elements)

    charge = self._total_charge(submodel, metal_element, state)
    frozen_iseqs, _free = qm.determine_frozen_atoms(
      submodel, metal_iseq, cap_iseqs, log=self.logger)

    # Intermediary xTB files share the .edits naming scheme
    # (<model>_<Metal>_<Chain><Resseq>); xtb_manager prepends 'xtb_'.
    preamble = f"qmi_{self._model_stem}_{site_tag}"
    opt_sites, qmm = qm.run_constrained_opt(
      submodel, frozen_iseqs, charge=charge,
      multiplicity=state.multiplicity, solvent_model=self.solvent_model,
      preamble=preamble, work_dir=self.work_dir,
      log=self.logger)

    # Donors were detected above (crystal geometry); targets come from the
    # optimised geometry.
    # Modified-Seminario sigmas from the Hessian at the optimised geometry
    # (plan Sec 5F); degrade gracefully to placeholder sigmas on any failure.
    bond_sigmas = angle_sigmas = None
    sigma_origin = "placeholder"
    try:
      hess_path = qm.run_hessian(
        submodel, opt_sites, charge=charge,
        multiplicity=state.multiplicity, solvent_model=self.solvent_model,
        preamble=f"{preamble}_hess", work_dir=self.work_dir,
        log=self.logger)
      n_atoms = submodel.get_hierarchy().atoms().size()
      hessian = sigma.read_xtb_hessian(hess_path, n_atoms)
      sigmas = sigma.seminario_sigmas(
        hessian, metal_iseq, donor_iseqs, opt_sites,
        temperature=self.temperature, log=self.logger)
      bond_sigmas, angle_sigmas = sigmas.bond_sigmas, sigmas.angle_sigmas
      sigma_origin = "hessian"
    except Exception as exc:
      print(f"Site {index}: Hessian/Seminario failed ({exc}); "
            f"using placeholder sigmas", file=self.logger)

    site_restraints = restraints.metal_ligand_restraints(
      submodel, opt_sites, metal_iseq, donor_iseqs,
      bond_sigmas=bond_sigmas, angle_sigmas=angle_sigmas,
      original_chain_ids=original_chain_ids, log=self.logger)
    sigma_label = ("modified-Seminario sigmas" if sigma_origin == "hessian"
                   else "placeholder sigmas")
    header = (f"qmi Phase-1 restraints: {metal_element} "
              f"ox={state.oxidation_state:+d} mult={state.multiplicity} "
              f"charge={charge:+d} "
              f"(GFN2-xTB/{self.solvent_model or 'gas'}; {sigma_label})")
    edits = restraints.as_edits_phil(site_restraints, header=header)

    skip_note = (f" ({site_restraints.n_skipped_symmetry} symmetry "
                 f"donor(s) skipped)"
                 if site_restraints.n_skipped_symmetry else "")
    print(f"Site {index}: {len(site_restraints.bonds)} bond + "
          f"{len(site_restraints.angles)} angle restraint(s){skip_note}",
          file=self.logger)

    return group_args(
      site_id=seed.get("file_name", f"site_{index}"),
      site_tag=site_tag,
      metal_element=metal_element,
      oxidation_state=state.oxidation_state,
      multiplicity=state.multiplicity, charge=charge, low_spin=state.low_spin,
      state_source=state.state_source, confidence=state.confidence,
      oxidation_reason=state.oxidation_reason, spin_reason=state.spin_reason,
      bvs=bvs_est,
      n_donors=len(donor_iseqs), restraints=site_restraints,
      edits=edits, sigma_origin=sigma_origin, skipped=None)

  def _total_charge(self, submodel, metal_element, state):
    """Total cluster charge for xTB ``--chrg``: the whole system given the
    metal's resolved oxidation state (engine guess or ``metal_state`` override),
    via ``quantum_interface.electrons`` (bond graph + per-element valences, so
    charged/deprotonated ligands count)."""
    sel = f"element {metal_element}"
    sac = [group_args(atom_selection=sel, charge=state.oxidation_state)]
    sam = [group_args(atom_selection=sel, multiplicity=state.multiplicity)]
    return quantum_interface.electrons(
      submodel, specific_atom_charges=sac,
      specific_atom_multiplicities=sam, log=self.logger)
