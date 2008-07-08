from __future__ import division

import math

from cctbx import crystal, xray
import smtbx.util
from smtbx import refinement
from smtbx.refinement import constraints

from scitbx import matrix as mat
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import random


class test_case(smtbx.util.test_case):

  def f(self):
    result = 0
    for sc in self.xs.scatterers():
      x = sc.site
      result += math.sin(x[0]*x[1]*x[2])
      if sc.flags.use_u_aniso():
        u = sc.u_star
        result += (u[0]**2 + 2*u[1]**2 + 3*u[2]**2
                   + 4*u[3]**2 + 5*u[4]**2 + 6*u[5]**2)
    return result

  def grad_f(self):
    result = flex.double()
    for sc in self.xs.scatterers():
      if sc.flags.grad_site():
        x = sc.site
        c = math.cos(x[0]*x[1]*x[2])
        result.extend(flex.double((x[1]*x[2]*c,
                                   x[0]*x[2]*c,
                                   x[0]*x[1]*c)))
      if sc.flags.use_u_aniso() and sc.flags.grad_u_aniso():
        u = sc.u_star
        result.extend(flex.double((2*u[0], 4*u[1], 6*u[2],
                                   8*u[3], 10*u[4], 12*u[5])))
    return result

def shift_site(sc, delta):
  sc.site = (mat.col(sc.site) + mat.col(delta)).elems

def shift_adp(sc, delta):
  sc.u_star = (mat.col(sc.u_star) + mat.col(delta)).elems


class special_position_test_case(test_case):

  def __init__(self):
    self.cs = crystal.symmetry((10, 10, 10, 90, 90, 90), "P432")
    x = 0.1
    u1, u2, u3, u4 = 0.04, -0.02, 0.03, -0.06
    self.xs = xray.structure(self.cs.special_position_settings())
    self.xs.add_scatterer(xray.scatterer("C1", site=(x, 1/2, 1/2),
                                               u=(u1, u2, u2, 0, 0, 0)))
    self.xs.add_scatterer(xray.scatterer("C2", site=(x, x, x),
                                               u=(u1, u1, u1, u2, u2, u2)))
    self.xs.add_scatterer(xray.scatterer("C3", site=(1/2, 1/2, 1/2),
                                               u=(u1, u1, u1, 0, 0, 0)))
    self.xs.add_scatterer(xray.scatterer("C4", site=(x, 1/2, 0),
                                               u=(u1, u2, u3, 0, 0, u4)))
    for i, sc in enumerate(self.xs.scatterers()):
      ops = self.xs.site_symmetry_table().get(i)
      assert ops.is_compatible_u_star(sc.u_star, tolerance=1e-6)
      sc.flags.set_grad_site(True)
      sc.flags.set_use_u_aniso(True)
      sc.flags.set_grad_u_aniso(True)

    self.constraint_flags = xray.scatterer_flags_array(
      len(self.xs.scatterers()))
    for f in self.constraint_flags:
      f.set_grad_site(True)
      f.set_grad_u_aniso(True)

  def shifted_structure(self, dx, du1, du2, du3, du4):
    result = self.xs.deep_copy_scatterers()
    sc0, sc1, sc2, sc3 = result.scatterers()
    f0, f1, f2, f3 = self.constraint_flags

    if f0.grad_site(): shift_site(sc0, (dx, 0, 0))
    if f0.grad_u_aniso(): shift_adp(sc0, (du1, du2, du2, 0, 0, 0))

    if f1.grad_site(): shift_site(sc1, (dx, dx, dx))
    if f1.grad_u_aniso(): shift_adp(sc1, (du1, du1, du1, du2, du2, du2))

    if f2.grad_u_aniso(): shift_adp(sc2, (du1, du1, du1, 0, 0, 0))

    if f3.grad_site(): shift_site(sc3, (dx, 0, 0))
    if f3.grad_u_aniso(): shift_adp(sc3, (du1, du2, du3, 0, 0, du4))

    return result

  def exercise(self):
    crystallographic_gradients = self.grad_f()
    parameter_map = self.xs.parameter_map()

    foo = (0,)*3

    # Reference gradients
    reparametrization_gradients_reference = flex.double(foo)
    for i,sc in enumerate(self.xs.scatterers()):
      ops = self.xs.site_symmetry_table().get(i)
      if self.constraint_flags[i].grad_site():
        j = parameter_map[i].site
        g = crystallographic_gradients[j:j+3]
        g1 = ops.site_constraints().independent_gradients(g)
        reparametrization_gradients_reference.extend(flex.double(g1))
    for i,sc in enumerate(self.xs.scatterers()):
      ops = self.xs.site_symmetry_table().get(i)
      if self.constraint_flags[i].grad_u_aniso():
        j = parameter_map[i].u_aniso
        g = crystallographic_gradients[j:j+6]
        g1 = ops.adp_constraints().independent_gradients(tuple(g))
        reparametrization_gradients_reference.extend(flex.double(g1))

    # Reference shifts
    crystallographic_shifts = flex.double()
    dx = 1e-4
    du1, du2, du3, du4 = -2e-4, -1e-4, 1e-4, 3e-4
    reference_shifted_xs = self.shifted_structure(dx, du1, du2, du3, du4)

    reparametrization_shifts = list(foo)
    f0, f1, f2, f3 = self.constraint_flags
    if f0.grad_site():
      reparametrization_shifts.append(dx)
    if f0.grad_u_aniso():
      reparametrization_shifts.extend((du1, du2))
    if f1.grad_site():
      reparametrization_shifts.append(dx)
    if f1.grad_u_aniso():
      reparametrization_shifts.extend((du1, du2))
    if f2.grad_u_aniso():
      reparametrization_shifts.append(du1)
    if f3.grad_site():
      reparametrization_shifts.append(dx)
    if f3.grad_u_aniso():
      reparametrization_shifts.extend((du1, du2, du3, du4))
    reparametrization_shifts = flex.double(reparametrization_shifts)

    # Testing...
    self.cts = constraints.special_positions(self.cs.unit_cell(),
                                             self.xs.site_symmetry_table(),
                                             self.xs.scatterers(),
                                             parameter_map,
                                             self.constraint_flags)

    for f in self.constraint_flags:
      assert not f.grad_site()
      assert not f.grad_u_aniso()

    # gradients
    reparametrization_gradients = flex.double(foo)
    self.cts.compute_gradients(crystallographic_gradients,
                               reparametrization_gradients)
    assert approx_equal(reparametrization_gradients,
                        reparametrization_gradients_reference)

    # apply_shifts
    original_xs = self.xs.deep_copy_scatterers()
    self.cts.apply_shifts(crystallographic_shifts,
                          reparametrization_shifts)
    assert approx_equal(self.xs.sites_frac(),
                        reference_shifted_xs.sites_frac())
    self.xs = original_xs

  def exercise_snap(self):
    cts = constraints.special_positions(self.cs.unit_cell(),
                                        self.xs.site_symmetry_table(),
                                        self.xs.scatterers(),
                                        self.xs.parameter_map(),
                                        self.constraint_flags)
    for sc in self.xs.scatterers():
      shift_site(sc, (1e-3,)*3)
      shift_adp(sc, (1e-3,)*6)
    cts.place_constrained_scatterers()
    xs = xray.structure(
      crystal.special_position_settings(crystal_symmetry=self.cs,
                                        min_distance_sym_equiv=0,
                                        u_star_tolerance=0,
                                        assert_min_distance_sym_equiv=True))
    xs.add_scatterers(self.xs.scatterers())
    assert (xs.site_symmetry_table().table()
            == self.xs.site_symmetry_table().table())


  def exercise_already_constrained(self):
    self.constraint_flags[1].set_grad_site(False)
    self.constraint_flags[2].set_grad_u_aniso(False)
    self.exercise()
    assert len(self.cts.already_constrained) == 2
    assert not self.cts.already_constrained[1].grad_site()
    assert not self.cts.already_constrained[2].grad_u_aniso()


class ch3_test_case(test_case):

  def __init__(self):
    self.cs = crystal.symmetry((8, 9, 10, 85, 95, 105), "P1")
    self.xs = xray.structure(self.cs.special_position_settings())
    pivot = xray.scatterer("C", site=(0.5, 0.5, 0.5),
                                u=(0.05, 0.04, 0.02,
                                   -0.01, -0.015, 0.005))
    self.xs.add_scatterer(pivot)
    self.i_pivot = 0
    v_CC = 1.54*(mat.col((-1, 2, 1.5)).normalize())
    v_CC = self.xs.unit_cell().fractionalize(v_CC)
    neighbour_site = mat.col(pivot.site) + mat.col(v_CC)
    self.xs.add_scatterer(xray.scatterer("C'", site=neighbour_site,
                                               u=(0,)*6))
    self.i_neighbour = 1
    for i in xrange(1,4):
      h = xray.scatterer("H%i" % i, u=(0,)*6)
      self.xs.add_scatterer(h)
    self.i_hydrogens = (2,3,4)
    for sc in self.xs.scatterers():
      sc.flags.set_grad_site(True)

    self.constraint_flags = xray.scatterer_flags_array(
      len(self.xs.scatterers()))
    for f in self.constraint_flags:
      f.set_grad_site(True)

    self.parameter_map = self.xs.parameter_map()

    self.cts = constraints.stretchable_rotatable_riding_terminal_X_Hn_array(
      self.cs.unit_cell(),
      self.xs.site_symmetry_table(),
      self.xs.scatterers(),
      self.parameter_map,
      self.constraint_flags)
    self.cts.append(constraints.stretchable_rotatable_riding_terminal_X_Hn(
      pivot=0,
      pivot_neighbour=1,
      hydrogens=(2,3,4),
      azimuth=50., # deg
      bond_length=1.))

  def pivot(self):
    return self.xs.scatterers()[self.i_pivot]
  pivot = property(pivot)

  def neighbour(self):
    return self.xs.scatterers()[self.i_neighbour]
  neighbour = property(neighbour)

  def hydrogens(self):
    return self.xs.scatterers()[self.i_hydrogens[0]:]
  hydrogens = property(hydrogens)

  def check_geometry(self, cartesian_frame=None):
    uc = self.cs.unit_cell()
    e0, e1, e2 = [ mat.col(e) for e in self.cts[0].local_cartesian_frame ]
    if cartesian_frame is not None:
      f0, f1, f2 = cartesian_frame
      for e,f in zip((e0,e1,e2), (f0,f1,f2)):
        assert approx_equal(e,f)
    else:
      assert approx_equal(abs(e0), 1)
      assert approx_equal(abs(e0), 1)
      assert approx_equal(abs(e0), 1)
      assert approx_equal(e0.cross(e1), e2)
    x_pivot = mat.col(uc.orthogonalize(self.pivot.site))
    x_neighbour = mat.col(uc.orthogonalize(self.neighbour.site))
    assert approx_equal(e2.cos_angle(x_pivot - x_neighbour), 1)
    x_h = [ mat.col(uc.orthogonalize(h.site)) for h in self.hydrogens ]
    dx = [x_neighbour - x_pivot] + [ x_h[i] - x_pivot for i in xrange(3) ]
    for i in xrange(4):
      for j in xrange(i+1,4):
        assert approx_equal(dx[i].cos_angle(dx[j]), -1/3)
    for i in xrange(1,4):
      assert approx_equal(abs(dx[i]), self.cts[0].bond_length)

  def exercise_riding(self):
    self.cts.place_constrained_scatterers()
    self.check_geometry()

    foo = (0,)*2

    crystallographic_gradients = self.grad_f()
    i_pivot_site = self.parameter_map[self.i_pivot].site
    grad_pivot_site = mat.col(crystallographic_gradients[i_pivot_site
                                                         :i_pivot_site+3])
    reparametrization_gradients = flex.double()
    self.cts.compute_gradients(crystallographic_gradients,
                               reparametrization_gradients)
    sum_grads_h = mat.col((0,0,0))
    for i_h in self.i_hydrogens:
      i_site = self.parameter_map[i_h].site
      grad_site = mat.col(crystallographic_gradients[i_site:i_site+3])
      sum_grads_h += grad_site
    assert approx_equal(
      crystallographic_gradients[i_pivot_site:i_pivot_site+3],
      grad_pivot_site + sum_grads_h)

    i_neigh_site = self.parameter_map[self.i_neighbour].site
    sum_grads = (mat.col(crystallographic_gradients[i_neigh_site
                                                    :i_neigh_site+3])
                 + grad_pivot_site + sum_grads_h)
    delta = mat.col((1e-3, -1.5e-3, 2e-3))
    shift_site(self.pivot, delta)
    shift_site(self.neighbour, delta)
    self.cts.place_constrained_scatterers()
    fp = self.f()
    shift_site(self.pivot, -2*delta)
    shift_site(self.neighbour, -2*delta)
    self.cts.place_constrained_scatterers()
    fm = self.f()
    true_delta_f = fp - fm
    riding_delta_f = sum_grads.dot(2*delta)
    assert approx_equal(true_delta_f, riding_delta_f)

  def exercise_rotate_stretch(self):
    foo = (0,)*3
    ct = self.cts[0]

    e0, e1, e2 = [ mat.col(e) for e in self.cts[0].local_cartesian_frame ]
    for rotating, stretching in [(True,True), (True,False),
                                 (False,True), (False,False)]:
      ct.rotating, ct.stretching = rotating, stretching
      assert self.cts[0].rotating == rotating
      assert self.cts[0].stretching == stretching

      self.cts.place_constrained_scatterers()
      crystallographic_gradients = self.grad_f()

      reparametrization_gradients = flex.double(foo)
      self.cts.compute_gradients(crystallographic_gradients,
                                 reparametrization_gradients)
      assert (len(reparametrization_gradients)
              == len(foo) + [rotating, stretching].count(True))
      if rotating:
        df_over_dphi = reparametrization_gradients[len(foo) + 0]
        if stretching:
          df_over_dl = reparametrization_gradients[len(foo) + 1]
      else:
        if stretching:
          df_over_dl = reparametrization_gradients[len(foo) + 0]

      if rotating:
        dphi = 1e-6 # deg
        phi = ct.azimuth
        ct.azimuth = phi + dphi
        self.cts.place_constrained_scatterers()
        self.check_geometry((e0,e1,e2))
        fp = self.f()
        ct.azimuth = phi - dphi
        self.cts.place_constrained_scatterers()
        self.check_geometry((e0,e1,e2))
        fm = self.f()
        df_over_dphi_approx = (fp - fm)/(2*dphi)
        assert approx_equal(df_over_dphi, df_over_dphi_approx)
        ct.azimuth = phi
      if stretching:
        h = 1e-6
        l = ct.bond_length
        ct.bond_length = l + h
        self.cts.place_constrained_scatterers()
        self.check_geometry((e0,e1,e2))
        fp = self.f()
        ct.bond_length = l - h
        self.cts.place_constrained_scatterers()
        self.check_geometry((e0,e1,e2))
        fm = self.f()
        df_over_dl_approx = (fp - fm)/(2*h)
        assert approx_equal(df_over_dl, df_over_dl_approx)
        ct.bond_length = l

def run():
  import sys
  verbose = '--verbose' in sys.argv[1:]
  ch3_test_case.run(verbose=verbose)
  special_position_test_case.run(verbose=verbose)
  print 'OK'

if __name__ == '__main__':
  run()
