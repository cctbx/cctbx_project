from cctbx import uctbx, xray
from smtbx.refinement import constraints
from scitbx import sparse
from scitbx import matrix as mat
from libtbx.test_utils import Exception_expected, approx_equal


class test_case(object):

  eps = 1.e-15
  bond_length = 0.9

  def __init__(self):
    self.uc = uctbx.unit_cell((1, 2, 3))
    self.c0 = xray.scatterer("C0", site=(0.1, 0.2, 0.3))
    self.c0.flags.set_grad_site(True)
    self.c1 = xray.scatterer("C1", site=(0., 0.1, 0.2))
    self.h = xray.scatterer("H")
    self.reparam = constraints.reparametrisation(self.uc)
    x0 = self.reparam.add(constraints.independent_site_parameter, self.c0)
    x1 = self.reparam.add(constraints.independent_site_parameter, self.c1)
    l = self.reparam.add(constraints.independent_scalar_parameter,
                         self.bond_length)
    x_h = self.reparam.add(constraints.terminal_linear_ch_site,
                           pivot=x0,
                           pivot_neighbour=x1,
                           length=l,
                           hydrogen=self.h)
    self.reparam.finalise()

  def run(self):
    self.reparam.linearise()
    self.reparam.store()
    x_c0, x_c1, x_h = [ mat.col(sc.site)
                        for sc in (self.c0, self.c1, self.h) ]
    assert approx_equal(self.uc.angle(x_c1, x_c0, x_h), 180, self.eps)
    assert approx_equal(
      self.uc.distance(x_c0, x_h), self.bond_length, self.eps)


def exercise():
  t = test_case()
  t.run()

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
