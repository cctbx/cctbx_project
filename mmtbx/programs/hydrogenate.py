"""Add hydrogens to a model"""

from __future__ import absolute_import, division, print_function
import os, time
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

add_h_to_water = False
  .type = bool

add_d_to_water = False
  .type = bool

print_time = False
  .type = bool
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

  # ----------------------------------------------------------------------------

  def validate(self):
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)

  # ----------------------------------------------------------------------------

  def run(self):
    t0 = time.time()
    self.model = self.data_manager.get_model()
    if self.data_manager.has_restraints():
      self.model.set_stop_for_unknowns(False)
      self.model.process(make_restraints=False)
    time_get_model = round(time.time()-t0, 2)

    make_sub_header('Add H atoms', out=self.logger)
    hydrogenate_obj = reduce_hydrogen.place_hydrogens(
      model                 = self.model,
      use_neutron_distances = self.params.use_neutron_distances,
      n_terminal_charge     = self.params.n_terminal_charge,
      print_time            = self.params.print_time)
    #import line_profiler
    #lp = line_profiler.LineProfiler(reduce_add_h_obj.run)
    #lp.enable()
    hydrogenate_obj.run()
    #lp.disable()
    #lp.print_stats()
    self.model = hydrogenate_obj.get_model()
    hydrogenate_obj.show(log = self.logger)

    if self.params.print_time:
      print("Get model obj in program template:", time_get_model)
      hydrogenate_obj.print_times()

    # Why is this done here and in the class?
    if self.params.add_h_to_water:
      self.model.add_hydrogens(1., occupancy=1.)
    elif self.params.add_d_to_water:
      self.model.add_hydrogens(1., element="D", occupancy=1.)

    if self.params.output.prefix is None:
      self.params.output.prefix = os.path.split(os.path.splitext(
        self.data_manager.get_default_model_name())[0])[1]

    if self.data_manager.get_model().input_model_format_cif():
      self.output_file_name = self.params.output.prefix+"_hydrogenate.cif"
      self.data_manager.write_model_file(
        model_str = self.model.model_as_mmcif(),
        filename  = self.output_file_name)
    else:
      self.output_file_name = self.params.output.prefix+"_hydrogenate.pdb"
      self.data_manager.write_model_file(
        model_str = self.model.model_as_pdb(),
        filename  = self.output_file_name)

    print("Wrote file: %s" % self.output_file_name, file=self.logger)

  # ----------------------------------------------------------------------------

  def get_results(self):
    return group_args(
     output_file_name=self.output_file_name)
