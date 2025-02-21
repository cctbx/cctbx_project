from __future__ import absolute_import, division, print_function
import platform
import sys
import time
from scitbx.array_family import flex
from libtbx import adopt_init_args
import scitbx.lbfgs
from scitbx import lbfgsb
from scitbx import minimizers
from libtbx.test_utils import approx_equal

class rosenbrock(object):
  def __init__(self, a, b, x,         # lbfgs and lbfgsb minimizer
               bound_flags    = None, # only lbfgsb
               lower_bound    = None, # only lbfgsb
               upper_bound    = None, # only lbfgsb
               initial_values = None  # only lbfgsb
               ):
    adopt_init_args(self, locals())
    assert self.x.size() == 2

  def update(self, x):
    self.x = x
    assert self.x.size() == 2

  def target(self):
    t = (self.a-self.x[0])**2+self.b*(self.x[1]-self.x[0]**2)**2
    return t

  def gradients(self):
    g1 = 2*(self.x[0]-self.a) + 4*self.b*(self.x[0]**3-self.x[0]*self.x[1])
    g2 = 2*self.b*(self.x[1]-self.x[0]**2)
    return flex.double([g1,g2])

  def curvatures(self):
    d1 = 2+4*self.b*(-self.x[1]+3*self.x[0]**2)
    d2 = 2*self.b
    d = flex.double([d1,d2])
    assert d.all_ne(0)
    return 1 / d


def run():
  tolerance = 1.e-16
  if sys.platform == 'win32':
    tolerance = 1.e-9
  elif (sys.platform == 'darwin' or sys.platform.startswith('linux')) \
    and ('arm64' in platform.machine() or 'arch64' in platform.machine()):
      tolerance = 1.e-6

  # Run L-BFGS (no boundaries)
  calculator = rosenbrock(a = 20, b = 10, x = flex.double([0,0]),
    bound_flags  = flex.int(2,2),
    lower_bound  = flex.double([-10000,-10000]),
    upper_bound  = flex.double([10000,10000]),
    initial_values = flex.double([0,0]))
  m_unbound = scitbx.minimizers.lbfgs(
    mode='lbfgs', max_iterations=100, calculator=calculator)
  #print('\tMinimum: ', list(m_unbound.x))
  res = (19.99999855596629, 399.99994289914525)
  assert approx_equal(m_unbound.x[0]-res[0],0, tolerance), (m_unbound.x[0], res[0], m_unbound.x[0]-res[0])
  assert approx_equal(m_unbound.x[1]-res[1],0, tolerance), (m_unbound.x[1], res[1], m_unbound.x[1]-res[1])

  # Run L-BFGS-B with boundaries
  calculator = rosenbrock(a = 20, b = 10, x = flex.double([0,0]),
    bound_flags  = flex.int(2,2),
    lower_bound  = flex.double([-10000,-10000]),
    upper_bound  = flex.double([10000,10000]),
    initial_values = flex.double([0,0]))
  m_bound = scitbx.minimizers.lbfgs(
    mode='lbfgsb', calculator=calculator)
  #print('\tMinimum: ', list(m_bound.x))
  res = (19.999999988074844, 399.99999950735986)
  assert approx_equal(m_bound.x[0]-res[0],0, tolerance), (m_bound.x[0], res[0], m_bound.x[0]-res[0])
  assert approx_equal(m_bound.x[1]-res[1],0, tolerance), (m_bound.x[1], res[1], m_bound.x[1]-res[1])

  # Run L-BFGS (no curvatures)
  calculator = rosenbrock(a = 1, b = 100, x = flex.double([-3,-4]))
  m_unbound = scitbx.minimizers.lbfgs(
    mode='lbfgs', calculator=calculator)
  #print('\tMinimum: ', list(m_unbound.x))
  res = (0.9999998308201578, 0.9999996829964546)
  assert approx_equal(m_unbound.x[0]-res[0],0, tolerance), (m_unbound.x[0], res[0], m_unbound.x[0]-res[0])
  assert approx_equal(m_unbound.x[1]-res[1],0, tolerance), (m_unbound.x[1], res[1], m_unbound.x[1]-res[1])


  # Run L-BFGS (with curvatures)
  calculator = rosenbrock(a = 1, b = 100, x = flex.double([-3,-4]))
  m_unbound2 = scitbx.minimizers.lbfgs(
    mode='lbfgs', calculator=calculator, diag_mode='always')
  #print('\tMinimum: ', list(m_unbound2.x))
  res = (1.0000002135019004, 1.000000406037043)
  assert approx_equal(m_unbound2.x[0]-res[0],0, tolerance), (m_unbound2.x[0], res[0], m_unbound2.x[0]-res[0])
  assert approx_equal(m_unbound2.x[1]-res[1],0, tolerance), (m_unbound2.x[1], res[1], m_unbound2.x[1]-res[1])


if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
