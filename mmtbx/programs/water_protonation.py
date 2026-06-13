"""H-bond-aware placement of hydrogens on water (HOH/DOD) residues.

Thin DataManager/PHIL wrapper around
:func:`mmtbx.hydrogens.water_protonation.place_water_hydrogens`. The
command-line dispatcher is ``mmtbx.naiad``.
"""

from __future__ import absolute_import, division, print_function

import os

from libtbx import group_args
from libtbx.program_template import ProgramTemplate
from libtbx.str_utils import make_sub_header

from mmtbx.hydrogens import water_protonation

master_phil_str = '''
oh_distance = *auto neutron xray
  .type = choice(multi=False)
  .short_caption = O-H bond length
  .help = "O-H bond length to place: neutron (0.984 A), xray (0.957 A), or auto to pick from the experiment type (neutron diffraction or D atoms present gives neutron, else X-ray)."
element = *auto H D
  .type = choice(multi=False)
  .short_caption = Hydrogen element
  .help = "Element for the placed water hydrogens: H, D, or auto (D for DOD residues, H for HOH)."
reorient_existing = False
  .type = bool
  .short_caption = Reorient existing water H
  .help = "Strip any H already on waters and re-place them, instead of leaving already-protonated waters untouched (the default)."
lone_pair = False
  .type = bool
  .short_caption = Lone-pair-directed placement
  .help = "Aim each O-H at an acceptor lone-pair lobe (derived from its bonded-neighbour geometry) rather than its nucleus, for better D-H...A angles."
joint = False
  .type = bool
  .short_caption = Joint H1/H2 optimization
  .help = "Optimize both O-H together: pick the orientation that best satisfies two acceptors at once, instead of placing H1 greedily then H2 on its cone. Slower; falls back to greedy in crowded pockets."
refine
  .short_caption = Refinement sweeps
{
  max_sweeps = 5
    .type = int(value_min=0)
    .short_caption = Maximum relaxation sweeps
    .help = "Maximum relaxation sweeps re-placing each water against the final environment to relax water-water clashes. Each sweep costs about one placement pass; 0 disables refinement. Sweeping stops early once a sweep removes fewer than tolerance close contacts, so a generous cap is safe."
  tolerance = 1
    .type = int(value_min=0)
    .short_caption = Early-stop tolerance
    .help = "Stop refining once a sweep removes fewer than this many close (<2.0 A) H-H contacts (1 = stop at a true plateau, larger = stop sooner on diminishing returns, 0 = run all max_sweeps). The best sweep is always kept."
}
basin
  .short_caption = Basin-hopping
{
  rounds = 0
    .type = int(value_min=0)
    .short_caption = Basin-hopping rounds
    .help = "Basin-hopping rounds after refinement (0 = off): each round randomly re-orients the still-clashing waters and relaxes, keeping the best. Deterministic (seeded); helps only where a better orientation exists."
}
stats = False
  .type = bool
  .short_caption = Report water-H clashes
  .help = "After placement, print the per-sweep water-H clash summary (count of H-H contacts below 2.0/1.8/1.5 A between different waters, and the closest contact)."
stats_worst = None
  .type = int(value_min=0)
  .short_caption = List worst contacts
  .help = "List the N closest residual water-H contacts with their residue IDs (implies stats=True)."
output {
  format = *mmcif pdb
    .type = choice(multi=False)
    .short_caption = Output format
    .help = "Output format. mmCIF by default (PDB format breaks on large structures); pdb falls back to mmCIF when the structure does not fit the standard PDB format."
}
'''


class Program(ProgramTemplate):
  description = '''
mmtbx.naiad: H-bond-aware placement of hydrogens on water residues.

Adds the two H atoms to every bare HOH/DOD oxygen in a model, orienting each
proton toward a nearby H-bond acceptor while staying clash-free against the
whole structure (including H placed on other waters) and keeping off metal
cations. Map-free and library-free -- placement is purely from geometry.

Inputs:
  PDB or mmCIF file containing an atomic model.
Output:
  Model with water hydrogens added (mmCIF by default), written to
  <model-stem>_waters_protonated.<ext> unless output.file_name is given.

By default it is idempotent: waters that already carry H are left untouched
(reorient_existing=True strips and re-places them).
'''
  datatypes = ['model', 'phil']
  master_phil_str = master_phil_str
  data_manager_options = ['model_skip_expand_with_mtrix',
                          'model_skip_ss_annotations']

  # ----------------------------------------------------------------------------

  def validate(self):
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)
    # Non-negativity of the int params is enforced by value_min=0 in the PHIL.
    if self.params.stats_worst is not None:
      self.params.stats = True   # listing offenders only makes sense with stats

  # ----------------------------------------------------------------------------

  def run(self):
    model = self.data_manager.get_model()
    hier = model.get_hierarchy()
    pdb_in = model.get_model_input()

    # O-H distance: honour an explicit choice, else infer from the structure.
    if self.params.oh_distance == "auto":
      neutron, source = water_protonation._detect_neutron(pdb_in, hier)
      print("O-H distance: %s (auto: %s)"
            % ("neutron" if neutron else "X-ray", source), file=self.logger)
    else:
      neutron = self.params.oh_distance == "neutron"
      print("O-H distance: %s (forced)" % self.params.oh_distance,
            file=self.logger)
    oh = water_protonation._WATER_OH_NEUTRON if neutron \
        else water_protonation._WATER_OH_XRAY

    # Element: "auto" leaves the per-residue choice (DOD->D, HOH->H) to the
    # placer (element=None); "H"/"D" force it.
    placer_element = None if self.params.element == "auto" else self.params.element
    print("water hydrogen element: %s"
          % ("auto (H for HOH, D for DOD)" if placer_element is None
             else placer_element), file=self.logger)

    # Consistency warning: deuterium is modelled only from neutron (or joint)
    # data, whose O-D distances are the longer neutron value. Placing D at the
    # shorter X-ray length is geometrically inconsistent, so warn. The reverse
    # (H at the neutron length, e.g. a hydrogenous neutron structure) is
    # legitimate, so it is not flagged.
    places_d = (placer_element == "D"
                or (placer_element is None
                    and any(ag.resname.strip().upper() == "DOD"
                            for ag in hier.atom_groups())))
    if places_d and not neutron:
      print("warning: placing D (deuterium) at the X-ray O-H length (%.3f A); "
            "deuterium is normally neutron-derived and uses the longer neutron "
            "distance (%.3f A). Pass oh_distance=neutron for consistent "
            "geometry." % (water_protonation._WATER_OH_XRAY,
                           water_protonation._WATER_OH_NEUTRON),
            file=self.logger)

    n_before = self._count_water_h(hier)
    report_stats = self.params.stats
    n_worst = self.params.stats_worst if self.params.stats_worst is not None else 0

    # Stream the clash table as the sweeps complete. The summary line and
    # header print lazily on the first state (after the greedy pass, when the
    # H count is final), then each row is printed as its sweep finishes.
    header = {"printed": False}
    def on_state(label, stats):
      if not header["printed"]:
        print(self._summary(n_before, self._count_water_h(hier)),
              file=self.logger)
        print("water H clash report (%d placed; H-H between different "
              "waters):" % stats[0], file=self.logger)
        header["printed"] = True
      water_protonation._clash_row(label, stats, self.logger)

    make_sub_header('Placing water hydrogens', out=self.logger)
    kept_label = water_protonation.place_water_hydrogens(
      hier,
      oh_length          = oh,
      element            = placer_element,
      n_refine           = self.params.refine.max_sweeps,
      refine_tol         = self.params.refine.tolerance,
      n_basin            = self.params.basin.rounds,
      reorient_existing  = self.params.reorient_existing,
      lone_pair_directed = self.params.lone_pair,
      joint              = self.params.joint,
      on_state           = on_state if report_stats else None)

    n_after = self._count_water_h(hier)
    self.n_added = n_after - n_before
    self.kept_label = kept_label
    if not report_stats:
      print(self._summary(n_before, n_after), file=self.logger)
    elif header["printed"]:
      if kept_label is not None:
        print("  kept: %s" % kept_label, file=self.logger)
      if n_worst:
        worst = water_protonation._worst_water_clashes(hier, n_worst)
        if worst:
          print("  worst %d contacts:" % len(worst), file=self.logger)
          for d, id_a, id_b in worst:
            print("    %.2f A  %s  <->  %s" % (d, id_a, id_b),
                  file=self.logger)
        else:
          print("  no contacts < 2.0 A", file=self.logger)

    self.model = model
    self._write_output(model, hier)

  # ----------------------------------------------------------------------------

  def _write_output(self, model, hier):
    """Serialize the (in-place mutated) hierarchy and write via DataManager.

    mmCIF by default; 'pdb' is honoured only if the structure fits the
    standard PDB format, else it falls back to mmCIF with a note. The output
    extension is forced to match the format actually written.
    """
    # New H atoms come back with blank serial numbers; the mmCIF writer
    # rejects blanks as "invalid number literal". Re-number first.
    hier.atoms().reset_serial()
    cs = model.crystal_symmetry()

    as_pdb = (self.params.output.format == "pdb"
              and hier.fits_in_pdb_format(use_hybrid36=False))
    if self.params.output.format == "pdb" and not as_pdb:
      print("note: structure does not fit standard PDB format; "
            "writing mmCIF instead", file=self.logger)

    ext = "pdb" if as_pdb else "cif"
    self.output_file_name = self._output_file_name(ext)
    if as_pdb:
      txt = hier.as_pdb_string(crystal_symmetry=cs)
    else:
      txt = hier.as_mmcif_string(crystal_symmetry=cs)
    self.data_manager.write_model_file(
      model_str = txt,
      filename  = self.output_file_name,
      format    = ext)
    print("Wrote file: %s" % self.output_file_name, file=self.logger)

  def _output_file_name(self, ext):
    """``output.file_name`` (extension forced to ``ext``), else the default
    ``<model-stem>_waters_protonated.<ext>``."""
    # output.file_name and output.filename are PHIL aliases; honour either.
    given = self.params.output.file_name or self.params.output.filename
    if given is not None:
      return os.path.splitext(given)[0] + "." + ext
    stem = os.path.splitext(os.path.basename(
      self.data_manager.get_default_model_name()))[0]
    return "%s_waters_protonated.%s" % (stem, ext)

  @staticmethod
  def _count_water_h(hier):
    return sum(
      1 for ag in hier.atom_groups()
      if ag.resname.strip().upper() in ("HOH", "DOD")
      for a in ag.atoms()
      if a.element.strip().upper() in ("H", "D"))

  def _summary(self, n_before, n_now):
    added = n_now - n_before
    # With reorient_existing the existing H were stripped and re-placed, so
    # report the reoriented count (the net change alone reads as "+0").
    if self.params.reorient_existing and n_before:
      return ("water H/D atoms: %d -> %d (reoriented %d, +%d)"
              % (n_before, n_now, n_before, added))
    return "water H/D atoms: %d -> %d (+%d)" % (n_before, n_now, added)

  # ----------------------------------------------------------------------------

  def get_results(self):
    return group_args(
      model            = self.model,
      n_added          = self.n_added,
      kept_label       = self.kept_label,
      output_file_name = self.output_file_name)
