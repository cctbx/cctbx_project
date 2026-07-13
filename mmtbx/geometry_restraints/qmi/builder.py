"""QMRestraintsBuilder: the per-site QM metal-restraint pipeline (Phase 1).

Ties the pieces together for one model: extract a capped cluster per metal
site with the ``endo_exo`` builder (in-memory), fix the electronic state, run a
constrained GFN2-xTB optimisation, and emit ``geometry_restraints.edits``.

Phase-1 scope (plan Sec 11): the electronic-guess engine (Sec 6) is **out of
the path**: the (oxidation, multiplicity) state is *hardcoded* per metal via
:data:`PHASE1_STATES` (Zn(II) d10 singlet; Mn(II) d5 high-spin sextet).  The
total cluster charge is derived from that hardcoded metal oxidation state with
``quantum_interface.electrons`` (which walks the bond graph + per-element
valences), so deprotonated/charged ligands are accounted for.  Sigmas are the
placeholders in :mod:`restraints`; modified-Seminario curvature sigmas (Sec 5F)
come later.
"""

from __future__ import absolute_import, division, print_function

import libtbx.phil
from libtbx import group_args
from libtbx.utils import null_out, Sorry

from mmtbx.programs.endo_exo import master_phil_str as endo_exo_master_phil_str
from mmtbx.geometry_restraints.endo_exo.builder import QMRegionBuilder
from mmtbx.geometry_restraints import quantum_interface
from mmtbx.geometry_restraints.qmi import metals as qmi_metals
from mmtbx.geometry_restraints.qmi import qm, restraints, sigma


# Hardcoded Phase-1 electronic states: metal element -> (oxidation_state,
# multiplicity).  The electronic-guess engine (plan Sec 6) replaces this table
# in Phase 4; until then only these metals are handled (others are skipped).
#   Zn(II): d10, closed shell    -> singlet (mult 1, 0 unpaired) [Phase 1a]
#   Mn(II): d5, high spin        -> sextet  (mult 6, 5 unpaired) [Phase 1b]
PHASE1_STATES = {
  "Zn": group_args(oxidation_state=2, multiplicity=1),
  "Mn": group_args(oxidation_state=2, multiplicity=6),
}


def _default_output_filename(prefix=None, suffix=None):
  """Filename-stem callable for QMRegionBuilder when running in-memory (we set
  ``write_files=False``, but endo_exo still derives a stem per seed)."""
  return f"{prefix or ''}{suffix or ''}"


def _endo_exo_params(metal_element=None, selection=None, radius=5.0,
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
      metal endo_exo finds, and each is handled if it is in :data:`PHASE1_STATES`.
  selection : iterable of str or None
      Explicit endo_exo seed selection string(s); overrides ``metal_element``.
  radius : float
      endo_exo buffer radius (A).
  max_depth : int
      endo_exo BFS depth.
  solvent_model : str
      ALPB named solvent for xTB (default ``'ether'``); empty disables.
  temperature : float
      Temperature (K) for the modified-Seminario force-constant to sigma
      conversion (plan Sec 5F); the global sigma tightness scalar.
  work_dir : str or None
      Directory the xTB jobs run in (xTB writes scratch files into the cwd).
  logger : file-like or None
  """

  def __init__(self, metal_element=None, selection=None, radius=5.0,
               max_depth=3, solvent_model=qm.DEFAULT_SOLVENT,
               temperature=sigma.DEFAULT_TEMPERATURE_K, work_dir=None,
               logger=None):
    self.metal_element = list(metal_element) if metal_element else None
    self.selection = list(selection) if selection else None
    self.radius = radius
    self.max_depth = max_depth
    self.solvent_model = solvent_model
    self.temperature = temperature
    self.work_dir = work_dir
    self.logger = logger if logger is not None else null_out()
    self._results = []

  def run(self, model, model_name="qmi_model"):
    """Run the pipeline on *model*; return a list of per-site result objects.

    Each result is a ``libtbx.group_args`` with fields ``site_id``,
    ``metal_element``, ``oxidation_state``, ``multiplicity``, ``charge``,
    ``n_donors``, ``restraints`` (the :func:`restraints.metal_ligand_restraints`
    record), ``edits`` (the PHIL string), and ``skipped`` (None, or a reason).
    """
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
    state = PHASE1_STATES.get(metal_element)
    if state is None:
      return _skip(f"metal {metal_element} not in Phase-1 scope "
                   f"(no hardcoded state)", metal_element=metal_element)

    print(f"Site {index}: {metal_element}, hardcoded "
          f"ox={state.oxidation_state:+d} mult={state.multiplicity}",
          file=self.logger)

    # Restraints manager (needed for charge accounting + cap-anchor freezing).
    qm.ensure_restraints(submodel, cap_iseqs=cap_iseqs,
                         cap_original_elements=cap_original_elements)

    charge = self._total_charge(submodel, metal_element, state)
    frozen_iseqs, _free = qm.determine_frozen_atoms(
      submodel, metal_iseq, cap_iseqs, log=self.logger)

    opt_sites, qmm = qm.run_constrained_opt(
      submodel, frozen_iseqs, charge=charge,
      multiplicity=state.multiplicity, solvent_model=self.solvent_model,
      preamble=f"qmi_site{index:03d}", work_dir=self.work_dir,
      log=self.logger)

    # Donors from the crystal geometry; targets from the optimised geometry.
    donor_iseqs = restraints.detect_first_shell(submodel, metal_iseq)

    # Modified-Seminario sigmas from the Hessian at the optimised geometry
    # (plan Sec 5F); degrade gracefully to placeholder sigmas on any failure.
    bond_sigmas = angle_sigmas = None
    sigma_origin = "placeholder"
    try:
      hess_path = qm.run_hessian(
        submodel, opt_sites, charge=charge,
        multiplicity=state.multiplicity, solvent_model=self.solvent_model,
        preamble=f"qmi_site{index:03d}_hess", work_dir=self.work_dir,
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
      multiplicity=state.multiplicity, charge=charge,
      n_donors=len(donor_iseqs), restraints=site_restraints,
      edits=edits, sigma_origin=sigma_origin, skipped=None)

  def _total_charge(self, submodel, metal_element, state):
    """Total cluster charge for xTB ``--chrg``: the whole system given the
    metal's hardcoded oxidation state, via ``quantum_interface.electrons``
    (bond graph + per-element valences, so charged/deprotonated ligands
    count)."""
    sel = f"element {metal_element}"
    sac = [group_args(atom_selection=sel, charge=state.oxidation_state)]
    sam = [group_args(atom_selection=sel, multiplicity=state.multiplicity)]
    return quantum_interface.electrons(
      submodel, specific_atom_charges=sac,
      specific_atom_multiplicities=sam, log=self.logger)
