"""Validate ligands in a model"""
from __future__ import absolute_import, division, print_function
import json
import os
import traceback
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import mmtbx.validation.ligands
from mmtbx.validation import validate_ligands
from mmtbx.hydrogens import place_and_optimize_hydrogens
from mmtbx.hydrogens import reduce_hydrogen
import iotbx.pdb
from libtbx.utils import null_out, Sorry
from libtbx.str_utils import make_sub_header
from libtbx import group_args


master_phil_str = """
include scope mmtbx.validation.validate_ligands.master_params
include scope mmtbx.probe.Helpers.probe_phil_parameters
scattering_table = *n_gaussian wk1995 it1992 neutron electron
  .type = choice
  .short_caption = Scattering table
  .help = Scattering table for structure factors calculations
run_reduce2 = True
  .type = bool
add_flip_movers = True
  .type = bool
  .short_caption = Flip Asn/Gln/His
  .help = "When reduce2 places H, let it flip Asn and Gln amides and His rings
          (heavy atoms swapped, as in MolProbity) and choose the His
          protonation. The validated model, and the model written with
          save_reduce2_model, carry the flipped heavy atoms. The flips are
          listed in the log. False: no heavy atom moves, no Asn/Gln amide or
          His ring is flipped, also not to bind a metal ion. The His
          protonation is chosen either way."
include_symmetric_copies = True
  .type = bool
  .short_caption = Crystal contacts in H placement
  .help = "When reduce2 places H, let H orientations and flips account for
          atoms of symmetry-related molecules (crystal contacts), using the
          model's own crystal symmetry. Not used with map input, whose box is
          not a lattice, or without a unit cell. With map input the metal
          contacts of H placement are not searched across the box either."
metal_bound_amines_as_nh2 = True
  .type = bool
  .short_caption = Metal-bound amines as NH2
  .help = "When reduce2 places H, place a primary amine (NH3+) that coordinates
          a metal as NH2, as reduce2 does. Listed in the log."
save_reduce2_model = False
  .type = bool
save_map_coeffs = False
  .type = bool
  .short_caption = Save map coefficients (MTZ)
  .help = "Write 2mFo-DFc and mFo-DFc map coefficients (from fmodel) to <basename>_map_coeffs.mtz. Requires reflection data."
verbose = False
  .type = bool
gui
  .help = GUI-specific parameters
{
  output_dir = None
    .type = path
    .style = output_dir

  data_column_label = None
    .type = str
    .style = noauto renderer:draw_any_label_widget
    .input_size = 200

  free_column_label = None
    .type = str
    .style = noauto renderer:draw_any_label_widget
    .input_size = 200
}
"""

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.validate_ligands model.pdb data.mtz
phenix.validate_ligands model.pdb map.mrc [resolution=3.0]
phenix.validate_ligands model.pdb
phenix.validate_ligands model.pdb data.mtz interactions=True

Print out basic statistics for residue(s) with the given code(s), including
RSCC.

interactions=True adds the ligand interaction profile (experimental): H-bonds,
salt bridges, clashes and vdW contacts per ligand, in the log and in
<basename>_ligand_interactions.json. Its settings are
validate_ligands.ligand_interactions (e.g.
ligand_interactions.salt_bridge.criterion=charge_centre); its probe2 settings are
validate_ligands.ligand_interactions.probe, separate from the probe scope that
reduce2 uses.

To validate against X-ray data, pass the reflection file, so that omit maps can
be computed. To validate against a cryo-EM map, supply a map file. The ligand
is then compared against an Fcalc map.
'''

  datatypes = ['model', 'phil', 'restraint', 'miller_array', 'real_map']

  master_phil_str = master_phil_str
  data_manager_options = ['model_skip_expand_with_mtrix',
                          'model_skip_ss_annotations']

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs...\n', file=self.logger)

    # allow only one model
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)

    if (self.data_manager.get_default_miller_array_name() is not None and
        self.data_manager.get_default_real_map_name() is not None):
      raise Sorry(
        'Please supply either reflection data (for X-ray) \n'
        'or a map (for cryo-EM), not both.')

  # ---------------------------------------------------------------------------

  def add_hydrogens(self, model, crystal_symmetry=None,
                    input_crystal_symmetry='from_model'):
    '''
    Place H atoms with reduce2

    crystal_symmetry: the input's own symmetry for the symmetry copies, or
    None (no copies); input_crystal_symmetry: the symmetry of placement's
    metal contacts (None for map input); the flips, amines placed as NH2 and
    copies go to self.reduce2_report (show_reduce2_changes()).
    '''
    make_sub_header('Placing hydrogen atoms with reduce2', out=self.logger)
    self.reduce2_report = group_args(flips=[], metal_amines_nh2=[],
      symmetry_copies=[])
    try:
      self.working_model = place_and_optimize_hydrogens(
        model           = model,
        do_flips        = self.params.add_flip_movers,
        keep_existing_H = False,
        optimize_his_protonation = True,
        probe_phil      = self.params.probe,
        stop_for_unknowns = False,
        raise_on_missing = False,
        log             = self.logger,
        metal_bound_amines_as_nh2 = self.params.metal_bound_amines_as_nh2,
        crystal_symmetry = crystal_symmetry,
        input_crystal_symmetry = input_crystal_symmetry,
        report          = self.reduce2_report)
      self.working_model.unset_riding_h_manager()
    except Exception:
      msg = traceback.format_exc()
      print('Reduce2 failed.\n' + msg, file=self.logger)
      raise Sorry('Hydrogen placement (reduce2) failed; '
                  'see log for details.')

  # ---------------------------------------------------------------------------

  @staticmethod
  def flip_label(key):
    '''
    'GLN B 14A conformer A (amide)', with 'model N ' in front for a MODEL id,
    from a flip key of place_and_optimize_hydrogens()'s report.
    '''
    model_id, chain_id, resseq, icode, resname, altloc, kind = key
    label = '%s %s %s%s' % (resname.strip(), chain_id.strip(), resseq.strip(),
      icode.strip())
    if altloc.strip(): label += ' conformer %s' % altloc.strip()
    if model_id.strip(): label = 'model %s %s' % (model_id.strip(), label)
    return '%s (%s)' % (label, kind)

  def show_reduce2_changes(self, ligand_manager, near=5.0):
    '''
    Log what reduce2 changed beyond adding H: whether symmetry copies were
    used, the Asn/Gln/His whose final model is flipped against the input
    (none with add_flip_movers=False), each
    marked when one of its flipped atoms (OD1/ND2, OE1/NE2, the His ring) is
    within near A of a heavy atom of a validated ligand, in the model without
    symmetry; and the amines placed as NH2.
    '''
    report = getattr(self, 'reduce2_report', None)
    if report is None: return
    from mmtbx.hydrogens import FLIPPED_ATOMS
    make_sub_header('Changes by reduce2 besides adding H', out=self.logger)
    lines = report.symmetry_copies if self._copies_reason == 'used' else []
    if lines == ['used, none within 3.6 A']:
      print('Symmetry copies (crystal contacts): used, none within 3.6 A',
        file=self.logger)
    else:
      print('Symmetry copies (crystal contacts): %s' % self._copies_reason,
        file=self.logger)
      # the Optimizer's summary; without copies the reason above says why
      for line in lines:
        print('  ' + line, file=self.logger)
    if not self.params.add_flip_movers:
      print('Flips: off (add_flip_movers=False), none made',
        file=self.logger)
    else:
      print('Flipped Asn/Gln/His: %d' % len(report.flips), file=self.logger)
    ph = self.working_model.get_hierarchy()
    atoms = ph.atoms()
    ligands = [(lr.id_str, [a for a in atoms.select(lr.ligand_isel)
      if a.element.strip().upper() not in ('H', 'D')])
      for lr in ligand_manager]
    for key in report.flips:
      model_id, chain_id, resseq, icode, resname, altloc, kind = key
      names = FLIPPED_ATOMS[resname.strip().upper()]
      flipped = [a for m in ph.models() if m.id == model_id
        for c in m.chains() if c.id == chain_id
        for rg in c.residue_groups()
        if rg.resseq == resseq and rg.icode == icode
        for ag in rg.atom_groups()
        if ag.resname == resname and ag.altloc == altloc
        for a in ag.atoms() if a.name.strip() in names]
      close = []
      for id_str, lig in ligands:
        d = min([a.distance(b) for a in flipped for b in lig] or [None])
        if d is not None and d <= near:
          close.append('%s %.1f A' % (id_str, d))
      text = '  ' + self.flip_label(key)
      if close:
        text += '  * within %.0f A of ligand: %s' % (near, '; '.join(close))
      print(text, file=self.logger)
    if report.metal_amines_nh2:
      print('Metal-bound amines placed as NH2: %d' %
        len(report.metal_amines_nh2), file=self.logger)
      for text, contact in report.metal_amines_nh2:
        print('  %s  (%s)' % (text, contact), file=self.logger)
    elif not self.params.metal_bound_amines_as_nh2:
      print('Metal-bound amines as NH2: off (metal_bound_amines_as_nh2=False)',
        file=self.logger)
    else:
      print('Metal-bound amines placed as NH2: 0', file=self.logger)

  # ---------------------------------------------------------------------------

  @staticmethod
  def _remove_element_x(model):
    '''Drop pseudo-atoms with element "X"; they choke pdb_interpretation.'''
    return model.select(~model.selection('element X'))

  # ---------------------------------------------------------------------------

  def set_map_resolution(self, mmm):
    '''
    Settle the resolution used to compute the cryo-EM Fcalc map.
    '''
    d_min = self.params.validate_ligands.resolution
    if d_min is not None:
      print('Using supplied map resolution: %.2f A' % d_min, file=self.logger)
      return
    try:
      d_min = mmm.resolution()
    except Exception:
      d_min = None
    if d_min is None:
      raise Sorry('Could not determine the resolution of the input map. '
                  'Please supply it with validate_ligands.resolution=<d_min>.')
    print('Map resolution estimated from the map: %.2f A '
          '(override with validate_ligands.resolution=)' % d_min,
          file=self.logger)
    self.params.validate_ligands.resolution = d_min

  # ---------------------------------------------------------------------------

  def check_ligands(self, model):
    make_sub_header('Check if input model has ligands', out=self.logger)
    get_class = iotbx.pdb.common_residue_names_get_class
    exclude = ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
               "modified_rna_dna", "ccp4_mon_lib_rna_dna", "common_water",
                "common_element"]
    self.has_ligands = False
    model = self._remove_element_x(model)
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          if (get_class(name=ag.resname) in exclude): continue
          print('Found ligand: ', ag.resname, file=self.logger)
          self.has_ligands = True

  # ---------------------------------------------------------------------------

  def run(self):
    has_miller = False
    has_map = False
    fmodel = None
    map_manager = None
    self.ligand_manager = None
    self.model_fn_reduce2 = None
    self.interactions_fn = None
    model_fn = self.data_manager.get_default_model_name()
    self._original_model_fn = model_fn
    data_fn = self.data_manager.get_default_miller_array_name()
    map_fn = self.data_manager.get_default_real_map_name()

    print('Using model file:', model_fn, file=self.logger)
    if data_fn is not None:
      print('Using reflection file:', data_fn, file=self.logger)
      has_miller = True
    if map_fn is not None:
      print('Using map file', map_fn, file=self.logger)
      has_map = True

    # get model object from input file
    m = self.data_manager.get_model()
    m.set_log(log = null_out())
    # before process() and reduce2, which box a model without symmetry
    output_cs = reduce_hydrogen.get_output_crystal_symmetry(m)
    if self.data_manager.has_restraints():
      m.set_stop_for_unknowns(False)
      #m.set_log(log = null_out())
      m.process(make_restraints=False)

    # stop if multi-model file
    if(len(m.get_hierarchy().models())>1):
      raise Sorry('Multi-model files currently not supported.')

    self.check_ligands(model = m)
    if not self.has_ligands:
      print('No ligands found. Exiting.', file=self.logger)
      return

    if self.params.validate_ligands.ligand_code:
      print('\nFocusing on the following ligands only:', file=self.logger)
      for lc in self.params.validate_ligands.ligand_code:
        print('\t', lc, file=self.logger)

    # get rid of element X as it will choke pdb_interpretation
    if ' X' in m.get_hierarchy().atoms().extract_element():
      print('\nFound atoms with element "X" in model. Removing...',
        file=self.logger)
      m = self._remove_element_x(m)

    if has_map:
      # Adopt the map's crystal symmetry before anything else touches the model.
      # A cryo-EM model usually carries only a dummy CRYST1 (1 1 1 P 1), and
      # with no real symmetry reduce2 puts a P1 box around the model and moves
      # the coordinates into it. map_model_manager later relabels the model with
      # the map's symmetry but does not move it back.
      m.set_crystal_symmetry(
        self.data_manager.get_real_map(map_fn).crystal_symmetry())

    self.working_model = None

    # copies from the input's own symmetry, before the box; none for a map
    copies_cs = None
    if has_map:
      self._copies_reason = 'none: map input (the map box is not a lattice)'
    elif not self.params.include_symmetric_copies:
      self._copies_reason = 'off (include_symmetric_copies=False)'
    elif reduce_hydrogen.contact_symmetry(output_cs) is None:
      self._copies_reason = ('none: no unit cell, or the CRYST1 1 1 1 '
        'placeholder')
    else:
      copies_cs = output_cs
      self._copies_reason = 'used'
    # the map box is not a lattice: no metal contacts across it either
    if self.params.run_reduce2:
      self.add_hydrogens(model = m, crystal_symmetry = copies_cs,
        input_crystal_symmetry = None if has_map else output_cs)
    else:
      self.working_model = m

    if self.working_model is None:
      raise Sorry('Could not create model object.')

    # Register the working model (with reduce2 H placed and element-X atoms
    # removed) so downstream managers use it instead of the original input
    # model. This does not change the DataManager default model.
    # In-memory DataManager key for the working model; never written to disk.
    _model_fn = '%s_validate_ligands_working.pdb' % os.path.splitext(
      os.path.basename(model_fn))[0]
    self.data_manager.add_model(_model_fn, self.working_model)

    if has_map:
      mmm = self.data_manager.get_map_model_manager(
        model_file = _model_fn)
      map_manager = mmm.map_manager()
      self.working_model = mmm.model()
      # Map inputs are cryo-EM; electron scattering is physically correct here
      # regardless of the phil default (scattering_table only applies to the
      # fmodel/reflection path).
      print('Map input: treating as a cryo-EM map, using electron scattering '
            'factors.', file=self.logger)
      self.working_model.setup_scattering_dictionaries(scattering_table='electron')
      # keep the registered model in sync with the boxed map model
      self.data_manager.add_model(_model_fn, self.working_model)
      self.set_map_resolution(mmm = mmm)
      # the map's box is not a lattice
      if self.params.validate_ligands.interactions:
        self.params.validate_ligands.ligand_interactions.symmetry = False

    ro = self.working_model.get_restraint_objects()
    if ro is None: ro=[]
    ro_no_duplicates = []
    if ro:
      seen = set()
      for name, _ro in ro:
        if name not in seen:
          ro_no_duplicates.append((name, _ro))
          seen.add(name)

    #  print(_ro[1]['comp_list']['_chem_comp.three_letter_code'][0])

    # get fmodel object if reflection data were provided
    if has_miller:
      make_sub_header(' Creating fmodel object ', out=self.logger)
      fmodel_params = self.data_manager.get_fmodel_params()
      fmodel_params.xray_data.r_free_flags.required = False
      fmodel_params.xray_data.r_free_flags.ignore_r_free_flags = True
      self.data_manager.set_fmodel_params(fmodel_params)
      fmodel = self.data_manager.get_fmodel(
        scattering_table = self.params.scattering_table,
        model_filename   = _model_fn)
      print('\n', file = self.logger)
      fmodel.update_all_scales()
      fmodel.show(log=self.logger, show_header=False)
      print ("r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free()),
        file=self.logger)

    #self.data_manager.write_real_map_file(map_manager,filename="my_map.map")
    #self.data_manager.write_model_file(self.working_model,filename="my_model.pdb")

    self.working_model.set_restraint_objects(ro_no_duplicates)
    self.working_model.set_stop_for_unknowns(False)
    pi = self.working_model.get_current_pdb_interpretation_params()
    pi.pdb_interpretation.allow_polymer_cross_special_position = True
    try:
      self.working_model.process(
        make_restraints=True,
        pdb_interpretation_params = pi)
    except Exception as e:
      print(e, file=self.logger)
      print('Could not process model to create restraints.', file=self.logger)
      return


    # split(".")[0] strips any remaining extension after splitext, e.g. the
    # ".ent" left over from a compound extension like "pdb1avd.ent.gz".
    basename = os.path.splitext(os.path.basename(model_fn))[0].split(".")[0]
    self.model_fn_reduce2 = "%s_newH.cif" % basename
    if self.params.save_reduce2_model:
      # the reduce2 P1 box is not written; the map path keeps the map's cs
      write_cs = True
      if not has_map:
        if output_cs is None:
          write_cs = False
        else:
          self.working_model.set_unit_cell_crystal_symmetry(output_cs)
      self.data_manager.set_overwrite(True)
      self.data_manager.write_model_file(self.working_model,
        filename=self.model_fn_reduce2, format='cif', output_cs=write_cs)

    if self.params.save_map_coeffs:
      if fmodel is not None:
        map_coeffs_fn = "%s_map_coeffs.mtz" % basename
        mtz_object = validate_ligands.map_coefficients_as_mtz_object(fmodel)
        self.data_manager.set_overwrite(True)
        self.data_manager.write_miller_array_file(
          mtz_object, filename=map_coeffs_fn)
        print('Wrote map coefficients:', map_coeffs_fn, file=self.logger)
      else:
        print('save_map_coeffs requested but no reflection data were provided; '
              'skipping map output.', file=self.logger)

    ligand_manager = validate_ligands.manager(
      model = self.working_model,
      fmodel = fmodel,
      map_manager = map_manager,
      params = self.params.validate_ligands,
      log   = self.logger)
    ligand_manager.run()
    if self.params.run_reduce2:
      self.show_reduce2_changes(ligand_manager)
    ligand_manager.show_ligand_counts()
    ligand_manager.show_fragmentation()
    ligand_manager.show_sites_within()
    ligand_manager.show_nonbonded_overlaps()
    ligand_manager.show_table(out=self.logger)
    if self.params.validate_ligands.interactions:
      ligand_manager.show_interactions(out=self.logger)
      self.interactions_fn = "%s_ligand_interactions.json" % basename
      with open(self.interactions_fn, "w") as fh:
        json.dump(ligand_manager.interactions_as_dict(), fh, indent=1, sort_keys=True)
      print('\nWrote ligand interactions: %s' % self.interactions_fn, file=self.logger)

    self.ligand_manager = ligand_manager

  # ---------------------------------------------------------------------------

  def get_results(self):
    if (self.params.run_reduce2
        and self.params.save_reduce2_model
        and self.model_fn_reduce2 is not None):
      model_to_open = os.path.abspath(self.model_fn_reduce2)
    else:
      model_to_open = getattr(self, '_original_model_fn', None)
    ligand_results = None
    if self.ligand_manager is not None:
      ligand_results = [lr.as_picklable_snapshot()
                        for lr in self.ligand_manager]
    extra = {}
    if self.params.validate_ligands.interactions:
      extra['interactions_fn'] = getattr(self, 'interactions_fn', None)
    return group_args(
      working_model_fn = model_to_open,
      ligand_manager   = self.ligand_manager,
      ligand_results   = ligand_results,
      **extra)
