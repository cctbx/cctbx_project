import scitbx.array_family.flex

import boost.python
ext = boost.python.import_ext("scitbx_rigid_body_ext")
from scitbx_rigid_body_ext import *

class _tardy_model(boost.python.injector, ext.tardy_model):

  def minimization(O, max_iterations=None, callback_after_step=None):
    from scitbx.rigid_body.essence.tardy import refinery
    return refinery(
      tardy_model=O,
      max_iterations=max_iterations,
      callback_after_step=callback_after_step)
