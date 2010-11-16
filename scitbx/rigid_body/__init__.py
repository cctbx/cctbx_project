import scitbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("scitbx_rigid_body_ext")
from scitbx_rigid_body_ext import *

class _tardy_model(boost.python.injector, ext.tardy_model):

  def __getinitargs__(O):
    return (
      O.labels,
      O.sites,
      O.masses,
      O.tardy_tree,
      O.potential_obj,
      O.near_singular_hinges_angular_tolerance_deg)

  def __getstate__(O):
    version = 1
    return (version, O.pack_q(), O.pack_qd())

  def __setstate__(O, state):
    assert len(state) == 3
    assert state[0] == 1 # version
    O.unpack_q(q_packed=state[1])
    O.unpack_qd(qd_packed=state[2])

  def minimization(O, max_iterations=None, callback_after_step=None):
    from scitbx.rigid_body.essence.tardy import refinery
    return refinery(
      tardy_model=O,
      max_iterations=max_iterations,
      callback_after_step=callback_after_step)
