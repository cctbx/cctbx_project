from smtbx.refinement.ext import *

import smtbx.refinement.tests
import smtbx.refinement.minimization
import smtbx.refinement.barriers


from libtbx import object_oriented_patterns as oop
import cctbx.xray

class smtbx_structure_enhancements(oop.injector, cctbx.xray.structure):

  def parameter_map(self):
    return parameter_map(self.scatterers())


del oop
