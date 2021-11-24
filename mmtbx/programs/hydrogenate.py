from __future__ import absolute_import, division, print_function
import os
from libtbx.program_template import ProgramTemplate
#from libtbx.utils import null_out
from libtbx import group_args
from libtbx.str_utils import make_sub_header
from mmtbx.hydrogens import reduce_hydrogen

master_phil_str = '''
use_neutron_distances = False
  .type = bool
  .help = Use neutron distances

keep_existing_H = False
  .type = bool
  .help = Keep existing H atoms in the model

n_terminal_charge = *residue_one first_in_chain no_charge
  .type = choice(multi=False)
  .help = Mode for placing H3 at terminal nitrogen.

output
  .style = menu_item auto_align
{
  file_name_prefix = None
    .type = path
    .short_caption = Prefix for file name
    .help = Prefix for file name
    .input_size = 400
}
'''

# ------------------------------------------------------------------------------

class Program(ProgramTemplate):
  description = '''
Add hydrogens.

Inputs:
  PDB or mmCIF file containing atomic model
  Ligand CIF file, if needed
'''
  datatypes = ['model', 'restraint', 'phil']
  master_phil_str = master_phil_str
  data_manager_options = ['model_skip_expand_with_mtrix',
                          'model_skip_ss_annotations']

# ------------------------------------------------------------------------------

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

# ------------------------------------------------------------------------------

  def run(self):
    self.model = self.data_manager.get_model()
    #
    make_sub_header('Add H atoms', out=self.logger)
    reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
      model = self.model,
      use_neutron_distances = self.params.use_neutron_distances,
      n_terminal_charge = self.params.n_terminal_charge)
    #import line_profiler
    #lp = line_profiler.LineProfiler(reduce_add_h_obj.run)
    #lp.enable()
    reduce_add_h_obj.run()
    #lp.disable()
    #lp.print_stats()
    self.model = reduce_add_h_obj.get_model()
    reduce_add_h_obj.show(log = self.logger)
    #
    make_sub_header('Optimize H atoms', out=self.logger)
    self.model = reduce_hydrogen.optimize(model=self.model)
    #
    if(self.params.output.file_name_prefix is not None):
      base = self.params.output.file_name_prefix
    else:
      fp = self.data_manager.get_default_model_name()
      base = os.path.splitext(os.path.basename(fp))[0]
    of = open("%s_hydrogenate.pdb"%base,"w")
    of.write(self.model.model_as_pdb())
    of.close()

# ------------------------------------------------------------------------------

  def get_results(self):
    return group_args(model = self.model)

