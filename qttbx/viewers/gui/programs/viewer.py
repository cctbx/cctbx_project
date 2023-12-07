from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from libtbx import group_args
import sys
from pathlib import Path


# =============================================================================

class Program(ProgramTemplate):

  description = """
  Demo program to visualize models in a QT gui connected to either Molstar or ChimeraX

  """

  # Define the data types that will be used
  datatypes = ['model', 'real_map']

  # Input parameters

  # Note: include of map_model_phil_str allows specification of
  #   full_map, half_map, another half_map, and model and Program Template
  #   produces input_files.map_model_manager containing these

  master_phil_str = """

  viewer {
    choice = None
      .type = str
      .help = 'molstar' or 'chimerax'
      .short_caption = The viewer to use
    }

  """

  # Define how to determine if inputs are ok
  def validate(self):
   pass


  # Run the program
  def run(self):
    viewer_name = self.params.viewer.choice
    main(viewer=viewer_name,dm=self.datamanager,log=self.logger)

  # Get the results
  def get_results(self):
    return group_args()




if __name__ == '__main__':
  from iotbx.cli_parser import run_program
  run_program(program_class=Program)