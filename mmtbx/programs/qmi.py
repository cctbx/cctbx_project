"""mmtbx.qmi: QM-based geometry restraints for mononuclear metal sites.

For each mononuclear metal site in the model, extracts a capped QM cluster
(via endo_exo), runs one constrained GFN2-xTB optimisation at a fixed
electronic state, derives per-restraint sigmas from the modified-Seminario
force constants, and writes one Phenix ``geometry_restraints.edits`` file per
site (``<model_stem>_<site>_qmi.edits``), or a single combined file with
``combine=True``.

Phase-1 scope: only mononuclear sites are handled. The electronic state is
chosen by the educated-guess engine (``electronic.OXIDATION_PRIORS``, currently
Zn/Mn/Mg/Fe/Cu) from the coordination geometry, or by a user ``metal_state``
override; a metal with no prior, or an ambiguous one with no dominant default
and no override (e.g. Fe), is skipped. The model must already be protonated
(standard residue names), or pass ``protonation.enable=True`` to add hydrogens
(reduce2 + water placement + donor deprotonation) during the run.

Usage (via dispatcher):

    mmtbx.development.qmi model.cif
    mmtbx.development.qmi model.cif sites.element_filter=Zn qm_region.radius=5
    mmtbx.development.qmi model.cif sites.selection="chain A and resseq 262"

Programmatic use::

    from mmtbx.programs.qmi import Program
    from iotbx.cli_parser import run_program
    run_program(program_class=Program)
"""

from __future__ import absolute_import, division, print_function

import os

try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate

from libtbx.utils import Sorry

from mmtbx.geometry_restraints.qmi.builder import QMRestraintsBuilder


master_phil_str = """
sites
  .short_caption = Metal sites
  .help = "Which metal site(s) to build restraints for, and optional per-site \
electronic-state overrides."
{
  selection = None
    .type = str
    .multiple = True
    .short_caption = Site selection
    .help = "Atom selection string(s) for the metal site(s); each entry seeds \
one site. If omitted, every metal in the structure is used."
  element_filter = None
    .type = str
    .multiple = True
    .short_caption = Element filter
    .help = "Restrict the metal scan to these element(s) (e.g. \
element_filter=Zn element_filter=Mn). Case-insensitive; only consulted when no \
selection is given."
  metal_state
    .multiple = True
    .short_caption = Metal state override
    .help = "Set the oxidation state and/or spin multiplicity for specific \
metal site(s), overriding the built-in guess. Repeat the block for different \
sites; a site matched by no block uses the built-in guess."
  {
    selection = None
      .type = str
      .help = "Atom selection for the site(s) this state applies to (e.g. \
'element Zn' covers every Zn). Required."
    oxidation_state = None
      .type = int
      .help = "Metal oxidation state for the matched site(s). Falls back to the \
built-in value when unset."
    multiplicity = None
      .type = int(value_min=1)
      .help = "Spin multiplicity (2S+1): number of unpaired electrons + 1 \
(e.g. 1 = singlet, 2 = doublet, 6 = high-spin d5). Falls back to the built-in \
value when unset."
  }
}
qm_region
  .short_caption = QM region
  .help = "Size of the QM cluster grown around each metal site."
{
  radius = 4.0
    .type = float(value_min=0)
    .short_caption = Region radius (A)
    .help = "Seed-inclusion radius (Angstrom): atoms within this distance of a \
metal are pulled into the QM region. The final capped cluster can extend beyond \
this radius."
  max_search_depth = 3
    .type = int(value_min=0)
    .expert_level = 2
    .short_caption = Max search depth
    .help = "Maximum number of bonds to expand outward when growing the QM \
region around each site."
}
protonation
  .short_caption = Protonation
  .help = "Add hydrogens to the model before building restraints."
{
  enable = False
    .type = bool
    .short_caption = Protonate model
    .help = "Add hydrogens (reduce2), place water hydrogens, and deprotonate \
metal-coordinating donors (Cys/Tyr/Asp/Glu; His tautomer-aware) before \
building restraints. When False (default), the input model must already be \
protonated or the run aborts."
  water
    .short_caption = Water hydrogens
    .help = "Water hydrogen placement options; only used when enable=True."
  {
    refine {
      max_sweeps = 5
        .type = int(value_min=0)
        .help = "Maximum passes that re-place each water against its final \
environment to relax water-water clashes (0 disables). Stops early once a pass \
clears fewer than 'tolerance' close contacts."
      tolerance = 1
        .type = int(value_min=0)
        .help = "Stop refining once a pass removes fewer than this many close \
(<2.0 A) H-H contacts (1 = stop at a plateau, 0 = run all max_sweeps). The best \
pass is always kept."
    }
    basin {
      rounds = 0
        .type = int(value_min=0)
        .expert_level = 2
        .help = "Extra rounds after refinement that randomly re-orient \
still-clashing waters and keep the best result (0 = off). Deterministic; helps \
only where a better orientation exists."
    }
    lone_pair = False
      .type = bool
      .expert_level = 2
      .help = "Aim each O-H at an acceptor's lone pair rather than at the \
acceptor atom, for better hydrogen-bond angles."
    joint = False
      .type = bool
      .expert_level = 2
      .help = "Optimize both O-H of a water together instead of one after the \
other, for a better fit to two acceptors at once. Slower; falls back to the \
simpler method in crowded pockets."
    stats = False
      .type = bool
      .expert_level = 2
      .help = "After placement, print the per-pass water-H clash summary and \
list the residual inter-water H-H contacts (grouped by 1.5/1.8/2.0 A)."
  }
}
qm
  .short_caption = QM / restraints
  .help = "QM optimisation and restraint-tightness settings."
{
  solvent = ether
    .type = str
    .short_caption = Solvent model
    .help = "Named solvent model for the xTB optimisation (default 'ether', \
closest to the protein-interior dielectric). Empty string disables solvation."
  temperature = 298.15
    .type = float(value_min=0)
    .short_caption = Temperature (K)
    .help = "Temperature (K) used to convert QM force constants into restraint \
sigmas: higher values give looser restraints (sigma scales with \
sqrt(temperature))."
  nproc = 1
    .type = int(value_min=1)
    .short_caption = Threads per xTB job
    .help = "Threads for each xTB optimisation and Hessian."
  perturb_attempts = 1
    .type = int(value_min=1)
    .short_caption = Perturb-reoptimise attempts
    .help = "Perturb-and-reoptimise attempts when the inner sphere has an \
imaginary mode (each attempt uses a larger kick)."
}
output
  .style = menu_item auto_align
  .short_caption = Output
  .help = "Where and whether to write the restraint files."
{
  write_files = True
    .type = bool
    .short_caption = Write .edits files
    .help = "If True (default), write the .edits restraint file(s) to disk. \
Set False to consume results in memory via get_results()."
  combine = False
    .type = bool
    .short_caption = Combine into one file
    .help = "If True, write all sites' restraints to a single combined \
<model_stem>_qmi.edits file instead of one file per site."
  qm_work_dir = "qm_work_dir"
    .type = path
    .short_caption = QM working directory
    .help = "Directory to hold the intermediate per-site xTB files (created if \
absent). Set to '.' to write them in the current directory instead."
}
predict_only = False
  .type = bool
  .short_caption = Predict only
  .help = "Run only the electronic-state prediction (oxidation/spin + \
confidence) for each site and stop; no restraints are computed or written."
states = *predicted reasonable all
  .type = choice
  .short_caption = States to run
  .help = "Which metal states to run per site. predicted: the predictor commits \
one (oxidation, spin). reasonable: scan every accessible oxidation crossed with \
the spin(s) the site's ligand-field geometry permits (tetrahedral stays high-\
spin). all: scan every accessible oxidation crossed with both spin extremes, \
ignoring geometry. reasonable and all run each state in its own qm_work_dir \
subdirectory with its own restraints, and still report what the predictor thinks."
use_bvs = False
  .type = bool
  .short_caption = Oxidation from BVS
  .help = "Commit the BVS most-likely oxidation instead of the predictor's."
"""


class Program(ProgramTemplate):
  """Generate QM metal restraints and write per-site geometry_restraints.edits.

  The per-site pipeline is delegated to
  :class:`mmtbx.geometry_restraints.qmi.builder.QMRestraintsBuilder`.
  """

  description = '''
Generate QM-based geometry restraints (bond + angle targets with per-restraint
sigmas) for mononuclear metal sites, as per-site Phenix
geometry_restraints.edits files. Targets come from a constrained GFN2-xTB
optimisation of each site; sigmas from the modified-Seminario force constants.
Requires xtb.
'''

  datatypes = ['model', 'phil']
  master_phil_str = master_phil_str

  def validate(self):
    """Validate user inputs before :meth:`run`."""
    if not self.data_manager.has_models():
      raise Sorry('No model provided. Please supply a PDB or mmCIF file.')
    if self.params.predict_only and self.params.states != 'predicted':
      raise Sorry(f'predict_only cannot be combined with states='
                  f'{self.params.states} (predict_only runs the predictor only).')
    model = self.data_manager.get_model()
    for sel_str in [s for s in (self.params.sites.selection or []) if s]:
      try:
        model.selection(sel_str)
      except Exception as e:
        raise Sorry(f"Invalid selection string '{sel_str}': {e}")
    for ms in (self.params.sites.metal_state or []):
      if (ms.selection is None and ms.oxidation_state is None
          and ms.multiplicity is None):
        continue  # empty default block from a PHIL .multiple scope
      if not ms.selection:
        raise Sorry("Each metal_state block requires a `selection` "
                    "(use e.g. 'element Zn' to cover every Zn).")
      try:
        model.selection(ms.selection)
      except Exception as e:
        raise Sorry(f"Invalid metal_state.selection '{ms.selection}': {e}")
      if ms.oxidation_state is None and ms.multiplicity is None:
        raise Sorry(f"metal_state block for '{ms.selection}' sets neither "
                    f"oxidation_state nor multiplicity.")
    if not self.params.predict_only:
      # predict_only is QM-free (state from donor geometry, not H), so it needs
      # neither protonation nor xtb.
      self._require_hydrogens(model)
      from mmtbx.geometry_restraints import xtb_manager
      if not xtb_manager.get_exe():
        raise Sorry("xtb executable not found. qmi requires xtb; set the "
                    "PHENIX_XTB environment variable or install xtb in the "
                    "build.")
    if self.params.protonation.enable:
      try:
        import mmtbx.programs.reduce2  # noqa: F401
      except ImportError:
        raise Sorry("protonation.enable=True requires mmtbx.reduce2, which is "
                    "not available in this build.")

  def _require_hydrogens(self, model):
    """Protonation guard: with protonation off the model must already be
    protonated. has_hd() is global (protein H alone satisfies it), so also
    require every water to carry a hydrogen: a bare-O water would enter the QM
    region unprotonated (two H for water, one for a hydroxo ligand; only a
    zero-H water is rejected)."""
    if self.params.protonation.enable:
      return
    if not model.has_hd():
      raise Sorry("The model does not appear to be protonated and qmi needs "
                  "hydrogens. Protonate it first (e.g. mmtbx.reduce2) or set "
                  "protonation.enable=True to add hydrogens during this run.")
    if self._has_unprotonated_water(model):
      raise Sorry("The model has a water with no hydrogens; qmi needs every "
                  "water protonated (two H for water, one for a hydroxo ligand). "
                  "Protonate the waters or set protonation.enable=True to place "
                  "water hydrogens.")

  def _has_unprotonated_water(self, model):
    """True if any water residue carries no H/D (a bare O). has_hd() is global
    (protein H alone passes it), so it misses a reduce-only model whose waters
    are bare; such a water would enter the QM region unprotonated."""
    water = model.get_hierarchy().select(model.selection("water"))
    for rg in water.residue_groups():
      if not any(a.element_is_hydrogen() for a in rg.atoms()):
        return True
    return False

  def run(self):
    """Build the restraints and (optionally) write the per-site .edits."""
    model = self.data_manager.get_model()
    element_filter = [e for e in (self.params.sites.element_filter or [])
                      if e and e.strip()]
    selection = [s for s in (self.params.sites.selection or []) if s]
    if self.params.protonation.enable:
      from mmtbx.geometry_restraints.qmi import protonation
      from mmtbx.geometry_restraints.qmi import metals as qmi_metals
      # Deprotonate metal-coordinating donors for every metal element
      # (qmi.metals.METALS), optionally narrowed by element_filter. This covers
      # engine-guessed metals AND any a user metal_state override forces. It is
      # chemically safe as a broad scope because deprotonate_metal_donors gates
      # the metal-induced Cys/Tyr rules to anion-inducing metals internally,
      # while carboxylate/His (metal-independent) apply to all.
      scope = set(qmi_metals.METALS)
      if element_filter:
        scope &= {e.strip().capitalize() for e in element_filter}
      wp = self.params.protonation.water
      model = protonation.prepare_model(
        model, metals=scope, log=self.logger,
        water_kwargs=dict(
          n_refine=wp.refine.max_sweeps, refine_tol=wp.refine.tolerance,
          n_basin=wp.basin.rounds, lone_pair_directed=wp.lone_pair,
          joint=wp.joint,
          report_stats=wp.stats))

    self._builder = QMRestraintsBuilder(
      metal_element=element_filter or None,
      selection=selection or None,
      radius=self.params.qm_region.radius,
      max_depth=self.params.qm_region.max_search_depth,
      solvent_model=self.params.qm.solvent,
      temperature=self.params.qm.temperature,
      metal_states=self.params.sites.metal_state,
      work_dir=self.params.output.qm_work_dir,
      predict_only=self.params.predict_only,
      states=self.params.states,
      nproc=self.params.qm.nproc,
      perturb_attempts=self.params.qm.perturb_attempts,
      use_bvs=self.params.use_bvs,
      logger=self.logger,
    )
    self._model_name = self.data_manager.get_default_model_name()
    self._results = self._builder.run(
      model, model_name=self._model_name,
      provenance=None if self.params.predict_only
      else self._global_provenance())
    self._show_summary()
    if (self.params.output.write_files and not self.params.predict_only
        and self.params.states == 'predicted'):
      self._write_edits()

  def get_results(self):
    """Return the list of per-site result objects (see QMRestraintsBuilder)."""
    return self._results

  # ------------------------------------------------------------------

  def _global_provenance(self):
    """Global provenance comment lines for the .edits header (plan Sec 9.3):
    tool/method, xtb + cctbx versions, model, protonation source/opt-out."""
    from mmtbx.geometry_restraints.qmi import provenance
    proto = ("reduce2 + naiad water-H + metal-donor deprotonation"
             if self.params.protonation.enable
             else "user-supplied (protonation off)")
    return provenance.global_block(
      method="GFN2-xTB", solvent=self.params.qm.solvent,
      temperature=self.params.qm.temperature, model_name=self._model_name,
      protonation=proto)

  def _show_summary(self):
    print('\nqmi metal-restraint summary:', file=self.logger)
    for r in self._results:
      if r.skipped:
        print(f'  {r.site_id}: skipped ({r.skipped})', file=self.logger)
        continue
      state = r.state_source
      if r.state_source in ("guessed", "bvs"):
        state = f"{r.state_source}/{r.confidence}-confidence"
      if r.restraints is None:   # predict_only: no QM, no restraints
        print(f'  {r.site_tag}: {r.metal_element} ox={r.oxidation_state:+d} '
              f'mult={r.multiplicity} state={state} (prediction only)',
              file=self.logger)
        continue
      print(f'  {r.site_tag}: {r.metal_element} ox={r.oxidation_state:+d} '
            f'mult={r.multiplicity} charge={r.charge:+d} '
            f'state={state} '
            f'({len(r.restraints.bonds)} bond, {len(r.restraints.angles)} '
            f'angle; sigma={r.sigma_origin})', file=self.logger)

  def _write_edits(self):
    emitted = [r for r in self._results if r.restraints is not None]
    if not emitted:
      print('No restraints emitted (all sites skipped); no .edits written.',
            file=self.logger)
      return
    stem = os.path.splitext(os.path.basename(self._model_name))[0]
    if self.params.output.combine:
      out = f'{stem}_qmi.edits'
      n_bonds = sum(len(r.restraints.bonds) for r in emitted)
      n_angles = sum(len(r.restraints.angles) for r in emitted)
      header = (f'qmi restraints for {self._model_name}: {len(emitted)} '
                f'site(s), {n_bonds} bond + {n_angles} angle')
      with open(out, 'w') as fh:
        fh.write(self._builder.combined_edits(header=header))
      print(f'  wrote {out} ({len(emitted)} site(s), {n_bonds} bond, '
            f'{n_angles} angle)', file=self.logger)
      return
    for r in emitted:
      out = f'{stem}_{r.site_tag}_qmi.edits'
      with open(out, 'w') as fh:
        fh.write(r.edits)
      print(f'  wrote {out} ({len(r.restraints.bonds)} bond, '
            f'{len(r.restraints.angles)} angle)', file=self.logger)
