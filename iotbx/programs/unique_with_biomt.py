# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate

from libtbx import Auto
import os

# =============================================================================

class Program(ProgramTemplate):

  description = '''
iotbx.unique_with_biomt: Tool to make an mmCIF file with only part of model
  and _pdbx_struct_assembly* records defining symmetry operations to reconstruct
  the rest of the model.

Usage examples:
  iotbx.unique_with_biomt model.pdb
  iotbx.unique_with_biomt model.pdb chain_id_to_leave='A'
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """\
  chain_id_to_leave = ''
    .type = str
    .help = If chain with this id is present, it will be kept as a unique part \
      of the model.
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(expected_n=1, exact_count=False, raise_sorry=True)
    m = self.data_manager.get_model()
    print ('Inputs OK', file=self.logger)

  # ---------------------------------------------------------------------------

  def run(self):
    # print(dir(self.data_manager))
    m = self.data_manager.get_model()
    m.search_for_ncs()
    m.setup_ncs_constraints_groups(filter_groups=False)
    if not m.can_be_unique_with_biomt():
      print("Model cannot be reduced.", file=self.logger)
      return

    if len(self.params.chain_id_to_leave) > 0:
      # make sure master has the correct chain
      nrgl = m.get_ncs_groups()
      found = False
      master_atom = m.get_hierarchy().atoms()[nrgl[0].master_iselection[0]]
      if master_atom.parent().parent().parent().id !=self.params.chain_id_to_leave:
        for i_c, c in enumerate(nrgl[0].copies):
          c_atom = m.get_hierarchy().atoms()[c.iselection[0]]
          if c_atom.parent().parent().parent().id == self.params.chain_id_to_leave:
            nrgl[0].make_nth_copy_master(i_c)
            found = True
      if not found:
        print("Chain id specified in chain_id_to_leave not found.", file=self.logger)
        print("Proceeding with a default value.", file=self.logger)
    cif_txt = m.model_as_mmcif(try_unique_with_biomt=True)
    inp_fn = os.path.basename(self.data_manager.get_default_model_name())[:-4]
    fn = "%s.cif" % self.get_default_output_filename(
        prefix='%s_' % inp_fn,
        suffix='unique_biomt',
        serial=Auto)
    print("Saving:", fn, file=self.logger)
    self.data_manager.write_model_file(model_str=cif_txt, filename=fn, format='cif')






  # ---------------------------------------------------------------------------
  def get_results(self):
    return None
