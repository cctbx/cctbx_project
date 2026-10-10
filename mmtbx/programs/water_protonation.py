"""H-bond-aware placement of hydrogens on water residues.

Thin DataManager/PHIL wrapper around
:func:`mmtbx.hydrogens.water_protonation.place_water_hydrogens`. The
command-line dispatcher is ``mmtbx.development.naiad``.
"""

from __future__ import absolute_import, division, print_function

import os

from iotbx.pdb.utils import check_for_missing_elements
from libtbx import group_args
from libtbx.program_template import ProgramTemplate
from libtbx.str_utils import make_sub_header
from libtbx.utils import Sorry, null_out

import mmtbx.model
from mmtbx.hydrogens import water_protonation

_MAX_LISTED_WATERS = 20

master_phil_str = '''
water_geometry = *auto xray neutron gas_phase
  .type = choice(multi=False)
  .short_caption = Water geometry
  .help = "O-H length and H-O-H angle to place: xray (0.850 A, 103.91 deg) or neutron (0.980 A, 103.91 deg), the cctbx restraint targets; gas_phase (0.957 A, 104.5 deg); or auto to pick xray or neutron from the experiment type (neutron diffraction or D atoms present gives neutron, else X-ray)."
element = *auto H D
  .type = choice(multi=False)
  .short_caption = Hydrogen element
  .help = "Element for the placed water hydrogens: H, D, or auto (the element the water already carries, else D for DOD residues and H for HOH)."
existing_h = *keep complete reorient
  .type = choice(multi=False)
  .short_caption = Waters that already carry H
  .help = "What to do with a water that already carries H: keep leaves it untouched, complete builds the missing partner on the cone of the existing proton for a water carrying a single H, reorient strips all water H and re-places both."
symmetry = True
  .type = bool
  .short_caption = Honour crystal symmetry
  .help = "Add the atoms that crystal symmetry places near a water to its environment, so H at a lattice contact avoid the neighbouring asymmetric units instead of pointing into them. Ignored for a model without a unit cell and space group, and for an electron microscopy model."
lone_pair = False
  .type = bool
  .short_caption = Lone-pair-directed placement
  .help = "Aim each O-H at an acceptor lone-pair lobe (derived from its bonded-neighbour geometry) rather than its nucleus, for better D-H...A angles."
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
output {
  suffix = _waters_protonated
    .type = str
    .help = "Suffix string added to automatically generated output filenames"
  serial = None
    .type = int
    .help = "Serial number added to automatically generated output filenames"
}
'''


class Program(ProgramTemplate):
  description = '''
mmtbx.development.naiad: H-bond-aware placement of hydrogens on water residues.

Adds the two H atoms to every bare water oxygen in a model, orienting each
proton toward a nearby H-bond acceptor while staying clash-free against the
whole structure (including H placed on other waters) and keeping off metal
cations. Map-free and library-free: placement is purely from geometry.

Inputs:
  PDB or mmCIF file containing an atomic model.
Output:
  Model with water hydrogens added, written to
  <model-stem>_waters_protonated.<ext> unless output.file_name is given. The
  format follows the input (output.target_output_format overrides it).

By default it is idempotent: waters that already carry H are left untouched.
That includes a water carrying a single H, since omitting the second proton
is a common way of writing hydroxide; those are reported rather than
completed. existing_h=complete builds the missing partner on the existing
proton's cone, and existing_h=reorient strips all water H and re-places both.
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
    model = self.data_manager.get_model()
    n_models = model.get_number_of_models()
    if n_models > 1:
      raise Sorry(f"Multi-model files are not supported ({n_models} models).")
    # Water, H and acceptor detection read the element column.
    try:
      check_for_missing_elements(
        model.get_hierarchy(),
        file_name=self.data_manager.get_default_model_name())
    except AssertionError as e:
      raise Sorry(str(e))
    self._warn_if_environment_unprotonated(model)

  # ----------------------------------------------------------------------------

  def run(self):
    model = self.data_manager.get_model()
    hier = model.get_hierarchy()
    pdb_in = model.get_model_input()

    # Water geometry: honour an explicit choice, else infer from the structure.
    geometry = self.params.water_geometry
    if geometry == "auto":
      neutron, source = water_protonation._detect_neutron(pdb_in, hier)
      geometry = "neutron" if neutron else "xray"
      print(f"Water geometry: {geometry} (auto: {source})", file=self.logger)
    else:
      print(f"Water geometry: {geometry} (forced)", file=self.logger)

    # Element: "auto" leaves the per-water choice (the element it carries,
    # else DOD->D, HOH->H) to the placer (element=None); "H"/"D" force it.
    placer_element = None if self.params.element == "auto" else self.params.element
    element_desc = ("auto (as carried, else H for HOH, D for DOD)"
                    if placer_element is None else placer_element)
    print(f"water hydrogen element: {element_desc}", file=self.logger)

    # Consistency warning: deuterium is modelled only from neutron (or joint)
    # data, whose O-D distances are the longer neutron value. The reverse (H at
    # the neutron length) is legitimate and is not flagged.
    places_d = (placer_element == "D"
                or (placer_element is None
                    and any(ag.resname.strip().upper() == "DOD"
                            or (ag.atoms().extract_element(strip=True)
                                == "D").count(True)
                            for ag in hier.atom_groups()
                            if water_protonation._is_water(ag.resname))))
    if places_d and geometry == "xray":
      xray = water_protonation._WATER_GEOMETRY["xray"][0]
      neut = water_protonation._WATER_GEOMETRY["neutron"][0]
      print(f"warning: placing D (deuterium) at the X-ray O-H length "
            f"({xray:.3f} A); deuterium is normally neutron-derived and uses "
            f"the longer neutron distance ({neut:.3f} A). Pass "
            f"water_geometry=neutron for consistent geometry.",
            file=self.logger)

    n_before = self._count_water_h(hier)
    report_stats = self.params.stats

    # Stream the clash table as the sweeps complete: the summary and header
    # print lazily on the first state, when the H count is final, then each
    # row as its sweep finishes.
    header = {"printed": False}
    def on_state(label, stats):
      if not header["printed"]:
        print(self._summary(n_before, self._count_water_h(hier)),
              file=self.logger)
        print(f"water H clash report ({stats[0]} placed; H-H between "
              f"different waters):", file=self.logger)
        header["printed"] = True
      water_protonation._clash_row(label, stats, self.logger)

    cs = None
    if self.params.symmetry:
      cs = model.crystal_symmetry()
      if cs is None or cs.unit_cell() is None or cs.space_group_info() is None:
        print("Crystal symmetry: none in the model, treating it as isolated",
              file=self.logger)
        cs = None
      elif pdb_in.get_experiment_type().is_electron_microscopy():
        # The cell of an electron microscopy model is the map box, not a
        # lattice.
        print("Crystal symmetry: ignored for electron microscopy, treating "
              "the model as isolated", file=self.logger)
        cs = None
      else:
        print(f"Crystal symmetry: {cs.space_group_info()}", file=self.logger)

    make_sub_header('Placing water hydrogens', out=self.logger)
    result = water_protonation.place_water_hydrogens(
      hier,
      geometry           = geometry,
      element            = placer_element,
      n_refine           = self.params.refine.max_sweeps,
      refine_tol         = self.params.refine.tolerance,
      n_basin            = self.params.basin.rounds,
      existing_h         = self.params.existing_h,
      lone_pair_directed = self.params.lone_pair,
      crystal_symmetry   = cs,
      on_state           = on_state if report_stats else None)

    n_after = self._count_water_h(hier)
    self.n_added = n_after - n_before
    self.kept_label = result.kept_label
    self.partial_waters = result.partial_waters
    if not report_stats:
      print(self._summary(n_before, n_after), file=self.logger)
    elif header["printed"]:
      if result.kept_label is not None:
        print(f"  kept: {result.kept_label}", file=self.logger)
      self._print_residual_contacts(hier, cs)
    self._print_partial_waters(result.partial_waters)

    self._write_output(model)
    # Placement added atoms to the hierarchy in place, past the loaded
    # model's caches; the result is a model built on it.
    hier.atoms().reset_i_seq()
    self.model = mmtbx.model.manager(
      model_input=None, pdb_hierarchy=hier,
      crystal_symmetry=model.crystal_symmetry(), log=null_out())

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
    """``output.file_name`` if given, else built from ``output.prefix``
    (default: the model file's stem), ``output.suffix`` and
    ``output.serial``. The writer appends the extension."""
    prefix = self.params.output.prefix
    if prefix is None:
      prefix = os.path.splitext(os.path.basename(
        self.data_manager.get_default_model_name()))[0]
    return self.get_default_output_filename(prefix=prefix)

  @staticmethod
  def _count_water_h(hier):
    return sum(water_protonation._n_hd(ag)
               for ag in hier.atom_groups()
               if water_protonation._is_water(ag.resname))

  def _summary(self, n_before, n_now):
    added = n_now - n_before
    # The net change alone reads as "+0" after a reorient.
    if self.params.existing_h == "reorient" and n_before:
      return (f"water H/D atoms: {n_before} -> {n_now} "
              f"(reoriented {n_before}, +{added})")
    return f"water H/D atoms: {n_before} -> {n_now} (+{added})"

  def _warn_if_environment_unprotonated(self, model):
    """Warn when nothing outside the waters carries a hydrogen."""
    hier = model.get_hierarchy()
    if not any(not water_protonation._is_water(ag.resname)
               for ag in hier.atom_groups()):
      return  # solvent-only model: there is nothing else to protonate
    if water_protonation.count_environment_hydrogens(hier):
      return
    print("warning: the model has no hydrogens outside its waters. Placement "
          "tests candidate positions against the surrounding atoms and reads "
          "a bonded H to tell a donor N from an acceptor, so without them the "
          "orientations degrade towards random. Protonating the rest of the "
          "model first is advisable.", file=self.logger)

  def _print_partial_waters(self, partial):
    """Report the waters that carried exactly one H on input.

    Metal-coordinating ones are annotated with the cation and its distance.
    """
    if not partial:
      return
    tail = {
      "kept": ("carry a single H and were left untouched (possible "
               "hydroxides); existing_h=complete adds the missing partner"),
      "completed": ("carried a single H and were completed; any intended as "
                    "hydroxide are now water"),
      "stripped": "carried a single H and were stripped and re-protonated",
    }
    # One mode can produce more than one outcome.
    for action in ("completed", "stripped", "kept"):
      group = [p for p in partial if p[2] == action]
      if not group:
        continue
      n_metal = sum(1 for _, metal, _ in group if metal is not None)
      head = f"{len(group)} water(s) {tail[action]}"
      if n_metal:
        head += f" ({n_metal} metal-coordinated)"
      print(f"  {head}:", file=self.logger)
      for rid, metal, _ in group[:_MAX_LISTED_WATERS]:
        note = f"  [{metal[0]} {metal[1]:.2f} A]" if metal is not None else ""
        print(f"    {rid}{note}", file=self.logger)
      if len(group) > _MAX_LISTED_WATERS:
        print(f"    ... and {len(group) - _MAX_LISTED_WATERS} more",
              file=self.logger)

  def _print_residual_contacts(self, hier, crystal_symmetry):
    """List the residual inter-water H-H contacts (< 2.0 A) with residue IDs,
    grouped into the 1.5/1.8/2.0 A bands and closest-first within each band.
    Contacts with symmetry equivalents carry the operator."""
    contacts = water_protonation._worst_water_clashes(
      hier, crystal_symmetry=crystal_symmetry)
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
      partial_waters   = self.partial_waters,
      output_file_name = self.output_file_name)
