"""Tools for PDB model manipulations"""
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from libtbx.program_template import ProgramTemplate
from mmtbx import pdbtools
from libtbx import Auto
import os
import mmtbx.pdbtools
from cctbx import uctbx

class Program(ProgramTemplate):

  description = '''
phenix.pdbtools tools for PDB model manipulations.

Usage examples:
  phenix.pdbtools model.pdb sites.shake=0.4
  phenix.pdbtools model.cif remove="element H"
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """\
include scope mmtbx.pdbtools.master_params

output {
  prefix = None
    .type = str
  suffix = _modified
    .type = str
  serial = None
    .type = int
  overwrite = True
    .type = bool
}
# temporary GUI PHIL
include scope libtbx.phil.interface.tracking_params
gui
  .help = "GUI-specific parameter required for output directory"
{
  output_dir = None
  .type = path
  .style = output_dir
}
"""

  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)

  def run(self):
    self.model = self.data_manager.get_model()
    cs = self.model.crystal_symmetry()
    if(cs is None or cs.is_empty() or cs.is_nonsense()):
      print("Crystal symmetry undefined, creating fake P1 box.")
      box_crystal_symmetry = \
        uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
          sites_cart   = self.model.get_sites_cart(),
          buffer_layer = 5).crystal_symmetry()
      self.model.set_crystal_symmetry(crystal_symmetry = box_crystal_symmetry)
    print('Performing manipulations', file=self.logger)
    self.model = mmtbx.pdbtools.modify(
      model  = self.model,
      params = self.params.modify,
      log    = self.logger).get_results().model
    # Write output model file
    input_file_name_base = os.path.basename(
      self.data_manager.get_default_model_name())[:-4]
    if(  self.model.input_model_format_cif()) or (
       not self.model.can_be_output_as_pdb()):
          extension = ".cif"
    elif(self.model.input_model_format_pdb()): extension = ".pdb"
    else:
      assert self.model.input_model_format_pdb() or \
        self.model.input_model_format_cif()
    if(self.params.output.prefix is not None):
      output_file_name = self.params.output.prefix
      if(self.params.output.suffix is not None):
        output_file_name = output_file_name + self.params.output.suffix
    else:
      output_file_name = input_file_name_base + self.params.output.suffix
    output_file_name = output_file_name + extension
    ofn = self.get_default_output_filename(
      prefix=output_file_name,
      suffix=None,
      serial=Auto)
    print('Writing output model to %s' %ofn, file=self.logger)
    output_cs=True
    if(cs is None): output_cs = False

    self.data_manager.write_model_file(self.model, ofn, output_cs=output_cs)
    self.result = ofn

  def get_results(self):
    return self.result

# So master_phil_str can be called
master_phil_str = Program.master_phil_str
