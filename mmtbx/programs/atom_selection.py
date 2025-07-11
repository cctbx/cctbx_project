"""Tool for selecting some atoms and optionally
  write out them as a PDB file"""
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate
from libtbx.utils import plural_s
from cctbx.array_family import flex
from cctbx import crystal, uctbx, xray
from libtbx.utils import Sorry
from cctbx.maptbx.box import shift_and_box_model

class Program(ProgramTemplate):

  description = '''
phenix.pdb_atom_selection: tool for selecting some atoms and optionally
  write out them as a PDB file

Usage examples:
  phenix.pdb_atom_selection model.pdb "chain A"
  phenix.pdb_atom_selection model.pdb "chain A" --write-pdb-file=sel.pdb
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """
  atom_selection_program {
    inselection = None
      .type = atom_selection
      .help = what to select
      .multiple = True
    cryst1_replacement_buffer_layer = None
      .type = float
      .help = replace original symmetry with P1 and size that is sufficient to \
        place molecule and surround it with this amount of empty space
    write_pdb_file = None
      .type = path
  }
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(raise_sorry=True)
    if (self.params.atom_selection_program.inselection is None or
        len(self.params.atom_selection_program.inselection) == 0):
      raise Sorry("Need selections")

  # ---------------------------------------------------------------------------
  def run(self):
    # I'm guessing self.data_manager, self.params and self.logger
    # are already defined here...

    # this must be mmtbx.model.manager?
    model = self.data_manager.get_model()
    atoms = model.get_atoms()
    all_bsel = flex.bool(atoms.size(), False)
    for selection_string in self.params.atom_selection_program.inselection:
      print("Selecting '%s'" % selection_string, file=self.logger)
      isel = model.iselection(string=selection_string)
      all_bsel.set_selected(isel, True)
      if self.params.atom_selection_program.write_pdb_file is None:
        print("  %d atom%s selected" % plural_s(isel.size()), file=self.logger)
        for atom in atoms.select(isel):
          print ("    %s" % atom.format_atom_record(), file=self.logger)
    print("", file=self.logger)
    if self.params.atom_selection_program.write_pdb_file is not None:
      ss_ann = model.get_ss_annotation()
      if not model.crystal_symmetry() or \
        (not model.crystal_symmetry().unit_cell()):
        model = shift_and_box_model(model, shift_model=False)
      selected_model = model.select(all_bsel)
      if(ss_ann is not None):
        selected_model.set_ss_annotation(ss_ann.\
            filter_annotation(
                hierarchy=selected_model.get_hierarchy(),
                asc=selected_model.get_atom_selection_cache(),
                remove_short_annotations=False,
                remove_3_10_helices=False,
                remove_empty_annotations=True,
                concatenate_consecutive_helices=False,
                split_helices_with_prolines=False,
                filter_sheets_with_long_hbonds=False))
      if self.params.atom_selection_program.cryst1_replacement_buffer_layer is not None:
        box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
            sites_cart=selected_model.get_atoms().extract_xyz(),
            buffer_layer=self.params.atom_selection_program.cryst1_replacement_buffer_layer)
        sp = crystal.special_position_settings(box.crystal_symmetry())
        sites_frac = box.sites_frac()
        xrs_box = selected_model.get_xray_structure().replace_sites_frac(box.sites_frac())
        xray_structure_box = xray.structure(sp, xrs_box.scatterers())
        selected_model.set_xray_structure(xray_structure_box)
      written_fname = self.data_manager.write_model_file(
          model_str=selected_model,
          filename=self.params.atom_selection_program.write_pdb_file)
      print("Wrote file: %s" % written_fname, file=self.logger)
      print("", file=self.logger)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return None
