from __future__ import division
from cctbx import uctbx, xray, sgtbx, crystal
from smtbx.refinement import constraints
from scitbx import sparse
from scitbx import matrix as mat
from libtbx.test_utils import Exception_expected, approx_equal


class terminal_linear_ch_site_test_case(object):

  eps = 1.e-15
  bond_length = 0.9

  def __init__(self, with_special_position_pivot):
    self.with_special_position_pivot = with_special_position_pivot
    self.uc = uctbx.unit_cell((1, 2, 3))
    self.sg = sgtbx.space_group("P 6")
    self.c0 = xray.scatterer("C0", site=(1e-5, 0., 0.1))
    self.site_symm = sgtbx.site_symmetry(self.uc, self.sg, self.c0.site)
    self.c0.flags.set_grad_site(True)
    self.c1 = xray.scatterer("C1", site=(0.09, 0.11, 0.))
    self.c1.flags.set_grad_site(True)
    self.h = xray.scatterer("H")
    self.reparam = constraints.reparametrisation(self.uc)
    if with_special_position_pivot:
      x0 = self.reparam.add(constraints.special_position_site,
                            self.site_symm, self.c0)
    else:
      x0 = self.reparam.add(constraints.independent_site_parameter, self.c0)
    x1 = self.reparam.add(constraints.independent_site_parameter, self.c1)
    l = self.reparam.add(constraints.independent_scalar_parameter,
                         self.bond_length, variable=False)
    x_h = self.reparam.add(constraints.terminal_linear_ch_site,
                           pivot=x0,
                           pivot_neighbour=x1,
                           length=l,
                           hydrogen=self.h)
    self.reparam.finalise()
    self.x0, self.x1, self.x_h, self.l = [ x.index for x in (x0, x1, x_h, l) ]
    if self.with_special_position_pivot:
      self.y0 = x0.independent_params.index

  def run(self):
    self.reparam.linearise()
    self.reparam.store()
    x_c0, x_c1, x_h = [ mat.col(sc.site)
                        for sc in (self.c0, self.c1, self.h) ]
    if self.with_special_position_pivot:
      assert approx_equal(x_c0, self.site_symm.exact_site(), self.eps)
    assert approx_equal(self.uc.angle(x_c1, x_c0, x_h), 180, self.eps)
    assert approx_equal(
      self.uc.distance(x_c0, x_h), self.bond_length, self.eps)

    if self.with_special_position_pivot:
      jt0 = sparse.matrix(1 + 3, # y0, x1
                          1 + 3 + 3 + 1 + 3) # y0, x0, x1, l, x_h
    else:
      jt0 = sparse.matrix(3 + 3, # x0, x1
                          3 + 3 + + 1 + 3) # x0, x1, l, x_h

    # Identity for independent parameters
    if self.with_special_position_pivot:
      jt0[self.y0, self.y0] = 1.
    for i in xrange(3): jt0[self.x0 + i, self.x0 + i] = 1.
    for i in xrange(3): jt0[self.x1 + i, self.x1 + i] = 1.

    # special position x0
    if self.with_special_position_pivot:
      jt0[self.y0, self.x0    ] = 0
      jt0[self.y0, self.x0 + 1] = 0
      jt0[self.y0, self.x0 + 2] = 1.

    # riding
    if self.with_special_position_pivot:
      jt0[self.y0, self.x_h + 2] = 1.
    else:
      for i in xrange(3): jt0[self.x0 + i, self.x_h + i] = 1.

    jt = self.reparam.jacobian_transpose
    assert sparse.approx_equal(self.eps)(jt, jt0)


class special_position_adp_test_case(object):

  eps = 1.e-15

  def __init__(self):
    cs = crystal.symmetry(uctbx.unit_cell((1, 1, 2, 90, 90, 120)), 'R3')
    sgi = sgtbx.space_group_info('R3(y+z, x+z, x+y+z)')
    op = sgi.change_of_basis_op_to_reference_setting()
    self.cs = cs.change_basis(op.inverse())
    self.sc = xray.scatterer('C',
                             site=(3/8,)*3,
                             u=(1/2, 1/4, 3/4, -3/2, -3/4, -1/4))
    self.sc.flags.set_grad_u_aniso(True)
    self.site_symm = sgtbx.site_symmetry(self.cs.unit_cell(),
                                         self.cs.space_group(),
                                         self.sc.site)
    self.reparam = constraints.reparametrisation(self.cs.unit_cell())
    u = self.reparam.add(constraints.special_position_cartesian_adp,
                         self.site_symm, self.cs.unit_cell(), self.sc)
    self.reparam.finalise()
    self.u, self.v = u.index, u.independent_params.index

  def run(self):
    self.reparam.linearise()
    self.reparam.store()
    assert approx_equal(self.sc.u_star, (19/6, 19/6, 17/2,
                                         11/6, 9/2, 9/2), self.eps)
    jt0 = sparse.matrix(2, 8)
    jt0[0, 0] = 1
    jt0[1, 1] = 1
    jac_u_cart_trans = self.site_symm.cartesian_adp_constraints(
      self.cs.unit_cell()).jacobian().transpose()
    for j in xrange(6):
      jt0[:, j + 2] = jac_u_cart_trans[:, j]
    jt = self.reparam.jacobian_transpose
    assert sparse.approx_equal(self.eps)(jt, jt0)


def exercise():
  terminal_linear_ch_site_test_case(with_special_position_pivot=False).run()
  terminal_linear_ch_site_test_case(with_special_position_pivot=True).run()
  special_position_adp_test_case().run()

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
