from __future__ import division

from libtbx import group_args
import libtbx.phil

master_phil_str = '''
cablam_idealization {
  nproc = 1
    .type = int
}
'''

# This is needed to import scope
master_phil = libtbx.phil.parse(master_phil_str)

class cablam_idealization(object):
  def __init__(self, model, params, log):
    """
    model is changed in place
    params - those in master_phil_str without scope name
    """
    self.model = model
    self.params = params

  def get_results(self):
    return group_args(model = self.model)
