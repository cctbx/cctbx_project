"""Tool for identifying barbed wire, near-folded, and other behaviors within AlphaFold2 predicted models"""
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
  type = *text kin selection_string selection_file json
    .type = choice
    .help = """choose output type"""
  file_name = None
    .type = str
    .short_caption = Output file name
    .help = Output file name
  filename = None
    .type = str
    .short_caption = Output file name
    .help = Output file name, optional spelling
  modes = *1 2 *3 4 5 6
    .type = choice(multi=True)
    .short_caption = Choose modes for selection string or file
    .help = Choose modes for selection string or file, subset of 123456, joined by +, eg 1+3
}
''' #+ Helpers.probe_phil_parameters

class Program(ProgramTemplate):
  description = '''
barbed_wire_analysis

This is a tool for identifying barbed wire, near-folded, and other behaviors within AlphaFold2 predicted models

output.type = *text kin selection_string selection_file json
output.file_name= (optional, prints to sys.stdout by default)

pdb or mmcif of residues matching selected behaviors can be printed with output.type=selection_file
output.modes = may be used to choose the printed behaviors with a string of numbers joined by +
1 = Predictive
2 = Unpacked high pLDDT
3 = Near-predictive
4 = Pseudostructure
5 = Barbed wire
6 = Unphysical
output.modes=1+3 is the default, and prints predictive + near-predictive

Examples:
phenix.barbed_wire_analysis your_prediction.pdb
phenix.barbed_wire_analysis your_prediction.pdb output.type=json output.filename=your_prediction_annotation.json
phenix.barbed_wire_analysis your_prediction.pdb output.type=kin output.file_name=your_prediction_markup.kin
phenix.barbed_wire_analysis your_prediction.pdb output.type=selection_file output.file_name=your_prediction_useful_parts.pdb
phenix.barbed_wire_analysis your_prediction.pdb output.type=selection_file modes=4+5+6 output.file_name=your_prediction_nonpredictive_only.pdb
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
    elif self.params.output.filename:
      out = open(self.params.output.filename, "w")
    else:
      out = sys.stdout
    if self.params.output.type == "text":
      bwa.as_text_chunks(out)
    elif self.params.output.type == "json":
      bwa.as_json(out)
    elif self.params.output.type == "kin":
      bwa.as_kinemage(out)
    elif self.params.output.type == "selection_string":
      print(bwa.as_selection_string(modes=self.params.output.modes), file=out)
    elif self.params.output.type == "selection_file":
      import os
      #Extension detection taken from mmtbx.pdbtools - output matches input
      input_file_name_base = os.path.basename(
        self.data_manager.get_default_model_name())[:-4]
      if(  model.input_model_format_cif()) or (
         not model.can_be_output_as_pdb()):
            extension = ".cif"
      elif(model.input_model_format_pdb()): extension = ".pdb"
      else:
        assert model.input_model_format_pdb() or \
          model.input_model_format_cif()
      bwa.as_selection_file(model, out, modes=self.params.output.modes, extension=extension)
