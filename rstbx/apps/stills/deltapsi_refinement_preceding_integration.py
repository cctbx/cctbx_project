from __future__ import absolute_import, division, print_function
from rstbx.apps.stills.simple_integration import IntegrationMetaProcedure
from rstbx.apps import simple_integration

class integrate_one_frame(IntegrationMetaProcedure):
  def __init__(self):
    simple_integration.__init__(self)
    IntegrationMetaProcedure.__init__(self)
