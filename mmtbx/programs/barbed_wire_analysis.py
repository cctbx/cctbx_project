from __future__ import absolute_import, division, print_function
import sys
from libtbx.program_template import ProgramTemplate

from mmtbx.validation import barbed_wire_analysis

version = "1.0.0"

master_phil_str = '''
profile = False
  .type = bool
  .short_caption = Profile the run
  .help = Profile the performance of the entire run

output
  .style = menu_item auto_align
{
  type = *text kin selection_string selection_file
    .type = choice
    .help = """choose output type"""
  file_name = None
    .type = str
    .short_caption = Output file name
    .help = Output file name
}
''' #+ Helpers.probe_phil_parameters

class Program(ProgramTemplate):
  description = '''
barbed_wire_analysis

This is a tool for identifying barbed wire, near-folded, and other behaviors within AlphaFold2 predicted models

output.type = *text kin selection_string selection_file
output.file_name= (optional, prints to sys.stdout by default)

Examples:
phenix.barbed_wire_analysis your_prediction.pdb
phenix.barbed_wire_analysis your_prediction.pdb output.type=kin output.file_name=your_prediction_markup.kin
phenix.barbed_wire_analysis your_prediction.pdb output.type=selection_file output.file_name=your_prediction_pruned.pdb
'''
  datatypes = ['model', 'restraint', 'phil']
  master_phil_str = master_phil_str
  data_manager_options = ['model_skip_expand_with_mtrix',
                          'model_skip_ss_annotations']
  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    model = self.data_manager.get_model()
    bwa = barbed_wire_analysis.barbed_wire_analysis(model)
    if self.params.output.file_name:
      out = open(self.params.output.file_name, "w")
    else:
      out = sys.stdout
    if self.params.output.type == "text":
      bwa.as_text_chunks(out)
    elif self.params.output.type == "kin":
      bwa.as_kinemage(out)
    elif self.params.output.type == "selection_string":
      print(bwa.as_selection_string(), file=out)
    elif self.params.output.type == "selection_file":
      bwa.as_selection_file(model, out)
