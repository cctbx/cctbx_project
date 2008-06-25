import boost.python
ext = boost.python.import_ext("smtbx_refinement_ext")
from smtbx_refinement_ext import *

from libtbx import object_oriented_patterns as oop
import cctbx.xray

class smtbx_structure_enhancements(oop.injector, cctbx.xray.structure):

  def parameter_map(self):
    return parameter_map(self.scatterers())

