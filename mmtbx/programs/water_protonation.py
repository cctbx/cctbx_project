"""H-bond-aware placement of hydrogens on water residues.

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
  .help = "After placement, print the per-sweep water-H clash summary and list the residual inter-water H-H contacts, grouped by the 1.5/1.8/2.0 A thresholds."
'''


class Program(ProgramTemplate):
  description = '''
mmtbx.naiad: H-bond-aware placement of hydrogens on water residues.

Adds the two H atoms to every bare water oxygen in a model, orienting each
proton toward a nearby H-bond acceptor while staying clash-free against the
whole structure (including H placed on other waters) and keeping off metal
cations. Map-free and library-free -- placement is purely from geometry.

Inputs:
  PDB or mmCIF file containing an atomic model.
Output:
  Model with water hydrogens added, written to
  <model-stem>_waters_protonated.<ext> unless output.file_name is given. The
  format follows the input (output.target_output_format overrides it).

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

  # ----------------------------------------------------------------------------

  def run(self):
    model = self.data_manager.get_model()
    hier = model.get_hierarchy()
    pdb_in = model.get_model_input()

    # O-H distance: honour an explicit choice, else infer from the structure.
    if self.params.oh_distance == "auto":
      neutron, source = water_protonation._detect_neutron(pdb_in, hier)
      kind = "neutron" if neutron else "X-ray"
      print(f"O-H distance: {kind} (auto: {source})", file=self.logger)
    else:
      neutron = self.params.oh_distance == "neutron"
      print(f"O-H distance: {self.params.oh_distance} (forced)",
            file=self.logger)
    oh = water_protonation._WATER_OH_NEUTRON if neutron \
        else water_protonation._WATER_OH_XRAY

    # Element: "auto" leaves the per-residue choice (DOD->D, HOH->H) to the
    # placer (element=None); "H"/"D" force it.
    placer_element = None if self.params.element == "auto" else self.params.element
    element_desc = ("auto (H for HOH, D for DOD)" if placer_element is None
                    else placer_element)
    print(f"water hydrogen element: {element_desc}", file=self.logger)

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
      xray = water_protonation._WATER_OH_XRAY
      neut = water_protonation._WATER_OH_NEUTRON
      print(f"warning: placing D (deuterium) at the X-ray O-H length "
            f"({xray:.3f} A); deuterium is normally neutron-derived and uses "
            f"the longer neutron distance ({neut:.3f} A). Pass "
            f"oh_distance=neutron for consistent geometry.", file=self.logger)

    n_before = self._count_water_h(hier)
    report_stats = self.params.stats

    # Stream the clash table as the sweeps complete. The summary line and
    # header print lazily on the first state (after the greedy pass, when the
    # H count is final), then each row is printed as its sweep finishes.
    header = {"printed": False}
    def on_state(label, stats):
      if not header["printed"]:
        print(self._summary(n_before, self._count_water_h(hier)),
              file=self.logger)
        print(f"water H clash report ({stats[0]} placed; H-H between "
              f"different waters):", file=self.logger)
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
        print(f"  kept: {kept_label}", file=self.logger)
      self._print_residual_contacts(hier)

    self.model = model
    self._write_output(model)

  # ----------------------------------------------------------------------------

  def _write_output(self, model):
    """Write the (in-place mutated) model via the DataManager.

    Format follows the input, or ``output.target_output_format`` when set,
    dropping to mmCIF for structures too large for the standard PDB format.
    The writer fixes the extension to match the format actually written.
    """
    # New H atoms come back with blank serial numbers; the mmCIF writer
    # rejects blanks as "invalid number literal". Re-number first.
    model.get_hierarchy().atoms().reset_serial()
    self.output_file_name = self.data_manager.write_model_file(
      model, filename=self._output_file_name())
    print(f"Wrote file: {self.output_file_name}", file=self.logger)

  def _output_file_name(self):
    """``output.file_name`` if given, else the default
    ``<model-stem>_waters_protonated``. The writer appends the extension."""
    # output.file_name and output.filename are PHIL aliases; honour either.
    given = self.params.output.file_name or self.params.output.filename
    if given is not None:
      return given
    stem = os.path.splitext(os.path.basename(
      self.data_manager.get_default_model_name()))[0]
    return f"{stem}_waters_protonated"

  @staticmethod
  def _count_water_h(hier):
    return sum(
      1 for ag in hier.atom_groups()
      if water_protonation._is_water(ag.resname)
      for a in ag.atoms()
      if a.element.strip().upper() in ("H", "D"))

  def _summary(self, n_before, n_now):
    added = n_now - n_before
    # With reorient_existing the existing H were stripped and re-placed, so
    # report the reoriented count (the net change alone reads as "+0").
    if self.params.reorient_existing and n_before:
      return (f"water H/D atoms: {n_before} -> {n_now} "
              f"(reoriented {n_before}, +{added})")
    return f"water H/D atoms: {n_before} -> {n_now} (+{added})"

  def _print_residual_contacts(self, hier):
    """List the residual inter-water H-H contacts (< 2.0 A) with residue IDs,
    grouped into the 1.5/1.8/2.0 A bands and closest-first within each band."""
    contacts = water_protonation._worst_water_clashes(hier)
    if not contacts:
      print("  no contacts < 2.0 A", file=self.logger)
      return
    print("  residual contacts:", file=self.logger)
    for label, lo, hi in (("< 1.5 A", 0.0, 1.5),
                          ("1.5-1.8 A", 1.5, 1.8),
                          ("1.8-2.0 A", 1.8, 2.0)):
      band = [c for c in contacts if lo <= c[0] < hi]
      print(f"    {label} ({len(band)}):", file=self.logger)
      for d, id_a, id_b in band:
        print(f"      {d:.2f} A  {id_a}  <->  {id_b}", file=self.logger)

  # ----------------------------------------------------------------------------

  def get_results(self):
    return group_args(
      model            = self.model,
      n_added          = self.n_added,
      kept_label       = self.kept_label,
      output_file_name = self.output_file_name)
