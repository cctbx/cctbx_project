from __future__ import absolute_import, division, print_function
import inspect
from libtbx import adopt_init_args
from libtbx.utils import user_plus_sys_time

def broadcast(m, log):
  if log is None: return
  print("-"*79, file=log)
  print(m, file=log)
  print("*"*len(m), file=log)

class manager(object):
  '''
  Clsss to execute Xtrapol8 calculations.
  f_obs_reference and f_obs_triggered are experimental amplitudes and expecetd
  to have consistent indexing.
  '''

  def __init__(self,
               model_reference,
               f_obs_reference,
               f_obs_triggered,
               log=None):
    adopt_init_args(self, locals())
    self.total_time = 0
    #
    self._caller(self._initialize)

  def run(self):
    self._caller(self._scale)
    self._caller(self._compute_weights)
    self._caller(self._compute_fo_minus_fo_map)

  def _caller(self, func):
    timer = user_plus_sys_time()
    doc = inspect.getdoc(func)
    broadcast(m=doc, log=self.log)
    func()
    t = timer.elapsed()
    self.total_time += t

  def _initialize(self):
    '''
    Start: initializing and validating inputs
    '''
    self.f_obs_reference, self.f_obs_triggered = \
      self.f_obs_reference.common_sets(self.f_obs_triggered)
    # Compute Riso and CCiso

  def _scale(self):
    '''
    Scale f_obs_reference and f_obs_triggered
    '''
    pass

  def _compute_weights(self):
    '''
    Compute k- and q-weights
    '''
    pass

  def _compute_fo_minus_fo_map(self):
    '''
    Compute F_diff = {weight * (F_obs1 - F_obs2), Phase_reference} map
    '''
    pass

  def _compute_f_extrapolated(self, alpha):
    '''
    Compute F_extrapolated
    '''
    pass
