"""mmtbx.qmi: QM-based geometry restraints for mononuclear metal sites.

For each mononuclear metal site in the model, extracts a capped QM cluster
(via endo_exo), runs one constrained GFN2-xTB optimisation at a fixed
electronic state, derives per-restraint sigmas from the modified-Seminario
force constants, and writes one Phenix ``geometry_restraints.edits`` file per
site (``<model_stem>_<site>_qmi.edits``), or a single combined file with
``combine=True``.

Phase-1 scope: the electronic state is hardcoded per metal (Zn(II) singlet,
Mn(II) high-spin sextet); other metals are skipped. The input model must be
protonated (e.g. with mmtbx.reduce2) and use standard residue names; in-model
protonation and the educated-guess electronic engine come later.

Usage (via dispatcher):

    mmtbx.development.qmi model.cif
    mmtbx.development.qmi model.cif element_filter=Zn buffer_radius=3.0
    mmtbx.development.qmi model.cif selection="chain A and resseq 262"

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
selection = None
  .type = str
  .multiple = True
  .help = "Atom selection string(s) for the metal site(s). Each entry seeds \
one site. If no selection is given, every metal in the structure is used \
(one per metal)."
element_filter = None
  .type = str
  .multiple = True
  .help = "Restrict the metal scan to these element(s) (e.g. element_filter=Zn \
element_filter=Mn). Case-insensitive. Only consulted when no `selection` is \
given."
metal_state
  .multiple = True
  .help = "Force the electronic state of the metal site(s) matched by \
`selection`. Authoritative: a matched block overrides the built-in guess, and \
supplying both values lets qmi run a metal that has no built-in state. Repeat \
the block to set different states for different sites; a site matched by no \
block falls through to the built-in guess."
{
  selection = None
    .type = str
    .help = "Atom-selection string identifying the metal site(s) this state \
applies to; the block is used for a site when its metal atom falls in the \
selection (use e.g. 'element Zn' to cover every Zn). Required."
  oxidation_state = None
    .type = int
    .help = "Metal oxidation state for the matched site(s). Falls back to the \
built-in value when left unset."
  multiplicity = None
    .type = int(value_min=1)
    .help = "Spin multiplicity 2S+1 for the matched site(s) (1 = singlet, \
2 = doublet, 6 = high-spin d5 sextet; = number of unpaired electrons + 1). \
Falls back to the built-in value when left unset."
}
buffer_radius = 3.0
  .type = float(value_min=0)
  .help = "Radius (Angstrom) of the endo_exo buffer region grown around each \
metal seed to build the QM cluster."
max_search_depth = 3
  .type = int(value_min=0)
  .help = "endo_exo BFS depth from any atom within the QM region."
protonate = False
  .type = bool
  .help = "Add hydrogens with mmtbx.reduce2, place water hydrogens (naiad), and \
deprotonate metal-coordinating donors (Cys/Tyr/Asp/Glu stripped to anionic; \
His tautomer-aware) before building restraints. When False (default), the \
model must already be protonated or the run aborts."
water_protonation
  .help = "Water-H placement quality knobs (naiad); only used when \
protonate=True."
{
  refine {
    max_sweeps = 5
      .type = int(value_min=0)
      .help = "Maximum relaxation sweeps re-placing each water against the \
final environment to relax water-water clashes. 0 disables refinement; \
sweeping stops early once a sweep clears fewer than 'tolerance' close contacts."
    tolerance = 1
      .type = int(value_min=0)
      .help = "Stop refining once a sweep removes fewer than this many close \
(<2.0 A) H-H contacts (1 = stop at a plateau, 0 = run all max_sweeps). The best \
sweep is always kept."
  }
  basin {
    rounds = 0
      .type = int(value_min=0)
      .help = "Basin-hopping rounds after refinement (0 = off): each round \
randomly re-orients still-clashing waters and relaxes, keeping the best. \
Deterministic (seeded); helps only where a better orientation exists."
  }
  lone_pair = False
    .type = bool
    .help = "Aim each O-H at an acceptor lone-pair lobe (from its \
bonded-neighbour geometry) rather than its nucleus, for better \
donor-H-acceptor angles."
  joint = False
    .type = bool
    .help = "Optimize both O-H of a water together (best orientation \
satisfying two acceptors at once) instead of placing H1 greedily then H2 on \
its cone. Slower; falls back to greedy in crowded pockets."
  stats = False
    .type = bool
    .help = "After placement, print the per-sweep water-H clash summary (count \
of H-H contacts below 2.0/1.8/1.5 A between different waters, and the closest)."
  stats_worst = None
    .type = int(value_min=0)
    .help = "List the N closest residual water-H contacts with their residue \
IDs (implies stats=True)."
}
solvent = ether
  .type = str
  .help = "ALPB named solvent for the xTB optimisation (nearest to the \
protein-interior dielectric; empty string disables solvation)."
temperature = 298.15
  .type = float(value_min=0)
  .help = "Temperature (K) for the modified-Seminario force-constant to sigma \
conversion; the global restraint-tightness scalar (sigma proportional to \
sqrt(temperature))."
write_files = True
  .type = bool
  .help = "If True (default), write the .edits file(s). Set False for \
in-memory use (consume Program(...).get_results() directly)."
combine = False
  .type = bool
  .help = "If True, write all sites' restraints to a single combined \
<model_stem>_qmi.edits file instead of one file per site."
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
    model = self.data_manager.get_model()
    for sel_str in [s for s in (self.params.selection or []) if s]:
      try:
        model.selection(sel_str)
      except Exception as e:
        raise Sorry(f"Invalid selection string '{sel_str}': {e}")
    for ms in (self.params.metal_state or []):
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
    self._require_hydrogens(model)
    from mmtbx.geometry_restraints import xtb_manager
    if not xtb_manager.get_exe():
      raise Sorry("xtb executable not found. qmi requires xtb; set the "
                  "PHENIX_XTB environment variable or install xtb in the "
                  "build.")
    if self.params.protonate:
      try:
        import mmtbx.programs.reduce2  # noqa: F401
      except ImportError:
        raise Sorry("protonate=True requires mmtbx.reduce2, which is not "
                    "available in this build.")

  def _require_hydrogens(self, model):
    """Protonation guard shared by validate() and run(): abort on an
    unprotonated model unless protonate=True. Called from run() too, so a
    programmatic caller that skips run_program's validate() is still guarded."""
    if not model.has_hd() and not self.params.protonate:
      raise Sorry("The model does not appear to be protonated and qmi needs "
                  "hydrogens. Protonate it first (e.g. mmtbx.reduce2) or set "
                  "protonate=True to add hydrogens during this run.")

  def run(self):
    """Build the restraints and (optionally) write the per-site .edits."""
    model = self.data_manager.get_model()
    self._require_hydrogens(model)
    p = self.params
    element_filter = [e for e in (p.element_filter or []) if e and e.strip()]
    selection = [s for s in (p.selection or []) if s]
    if p.protonate:
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
      wp = p.water_protonation
      model = protonation.prepare_model(
        model, metals=scope, log=self.logger,
        water_kwargs=dict(
          n_refine=wp.refine.max_sweeps, refine_tol=wp.refine.tolerance,
          n_basin=wp.basin.rounds, lone_pair_directed=wp.lone_pair,
          joint=wp.joint,
          report_stats=wp.stats or (wp.stats_worst is not None),
          n_worst=wp.stats_worst or 0))

    self._builder = QMRestraintsBuilder(
      metal_element=element_filter or None,
      selection=selection or None,
      radius=p.buffer_radius,
      max_depth=p.max_search_depth,
      solvent_model=p.solvent,
      temperature=p.temperature,
      metal_states=p.metal_state,
      logger=self.logger,
    )
    self._model_name = self.data_manager.get_default_model_name()
    self._results = self._builder.run(model, model_name=self._model_name)
    self._show_summary()
    if p.write_files:
      self._write_edits()

  def get_results(self):
    """Return the list of per-site result objects (see QMRestraintsBuilder)."""
    return self._results

  # ------------------------------------------------------------------

  def _show_summary(self):
    print('\nqmi metal-restraint summary:', file=self.logger)
    for r in self._results:
      if r.skipped:
        print(f'  {r.site_id}: skipped ({r.skipped})', file=self.logger)
        continue
      print(f'  {r.site_tag}: {r.metal_element} ox={r.oxidation_state:+d} '
            f'mult={r.multiplicity} charge={r.charge:+d} '
            f'state={r.state_source} '
            f'({len(r.restraints.bonds)} bond, {len(r.restraints.angles)} '
            f'angle; sigma={r.sigma_origin})', file=self.logger)

  def _write_edits(self):
    emitted = [r for r in self._results if r.restraints is not None]
    if not emitted:
      print('No restraints emitted (all sites skipped); no .edits written.',
            file=self.logger)
      return
    stem = os.path.splitext(os.path.basename(self._model_name))[0]
    if self.params.combine:
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
