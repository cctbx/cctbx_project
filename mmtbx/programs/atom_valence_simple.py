from __future__ import absolute_import, division, print_function

import os
# from datetime import datetime
from libtbx.program_template import ProgramTemplate
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  pass
#from libtbx.utils import Sorry
from mmtbx.ligands.hierarchy_utils import simple_valence_check

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s file.pdb [options ...]

Options:

Example:

  %(prog)s 1ubq.pdb
""" % locals()

  master_phil_str = """
"""
  datatypes = ['model','phil']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    model = self.data_manager.get_model()
    model.process(make_restraints=True)
    hierarchy = model.get_hierarchy()
    geometry_restraints_manager = model.get_restraints_manager()
    rc = simple_valence_check(hierarchy, geometry_restraints_manager)
    print(rc)

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return None
