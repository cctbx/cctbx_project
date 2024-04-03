"""
Holton geometry validation score calculation.  Based on a
script by James Holton (untangle_score.csh).  See details in
mmtbx/validation/holton_geometry_validation.py
"""

from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.holton_geometry_validation import \
    holton_geometry_validation
from libtbx.program_template import ProgramTemplate
from datetime import datetime

try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  pass

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
%(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file          input PDB file

Example:

  %(prog)s model=1ubq.pdb
""" % locals()

  master_phil_str = """
    model = None
    .type = path
    .optional = False
    .help = '''input PDB file'''

"""


  datatypes = ['model','phil']

  def validate(self):
    self.set_defaults()
    self.data_manager.has_models(raise_sorry=True)

  def set_defaults(self):

    if not self.params.model:
      self.data_manager.has_models(raise_sorry=True)
      self.params.model = self.data_manager.get_default_model_name()
    self.model = self.data_manager.get_model(self.params.model)
    print("\nModel read from %s" %(self.params.model), file = self.logger)

  def run(self):


    self.results = holton_geometry_validation(
      dm = self.data_manager,
      filename = self.params.model,
      model = self.model,
      log =self.logger,)

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return self.results.as_JSON(self.info_json)
