from __future__ import division
from rstbx.apps.stills.simple_integration import IntegrationMetaProcedure
from rstbx.apps import simple_integration

class integrate_one_frame(IntegrationMetaProcedure):
  def __init__(self, triclinic_pairs):
    simple_integration.__init__(self)
    IntegrationMetaProcedure.__init__(self)
    self.triclinic_pairs = triclinic_pairs
