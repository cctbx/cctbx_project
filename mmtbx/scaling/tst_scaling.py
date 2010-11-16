## Peter Zwart July 5, 2005
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import miller
from cctbx import xray
from cctbx import sgtbx
from cctbx import uctbx
from mmtbx import scaling
from libtbx.test_utils import approx_equal
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling import twin_analyses as t_a
import mmtbx.scaling

from scitbx.python_utils import random_transform
import random
import math
import time
from libtbx.str_utils import StringIO

random.seed(0)
flex.set_random_seed(0)

import scitbx.math as sm


##testing quick erf and quick eio
def test_luts():
  qerf = mmtbx.scaling.very_quick_erf(0.001)
  qeio = mmtbx.scaling.quick_ei0(5000)
  for i in xrange(-1000,1000):
    x=i/100.0
    assert approx_equal( qerf.erf(x), sm.erf(x), eps=1e-5 )
    if (x>=0):
      assert approx_equal( qeio.ei0(x), math.exp(-x)*sm.bessel_i0(x) , eps=1e-5 )
  number_of_iterations = 15000000
  for optimized in [False, True]:
    t0 = time.time()
    zero = qerf.loop_for_timings(number_of_iterations, optimized=optimized)
    print "very_quick_erf*%d optimized=%s: %.2f s" % (
      number_of_iterations, str(optimized), time.time()-t0)
    assert approx_equal(zero, 0)
  number_of_iterations = 5000000
  for optimized in [False, True]:
    t0 = time.time()
    zero = qeio.loop_for_timings(number_of_iterations, optimized=optimized)
    print "quick_ei0*%d optimized=%s: %.2f s" % (
      number_of_iterations, str(optimized), time.time()-t0)
    assert approx_equal(zero, 0)

## Testing Wilson parameters
def test_gamma_prot():
  gamma_prot_test = scaling.gamma_prot(0.011478)
  assert approx_equal(gamma_prot_test,-0.349085)

  gamma_prot_test = scaling.gamma_prot(0.028868)
  assert approx_equal(gamma_prot_test,-0.585563)

  d_star_sq = flex.double([0.011478,0.028868,1.0,0.0])
  gamma_array_test = scaling.get_gamma_prot(d_star_sq)
  assert approx_equal(gamma_array_test[0],-0.349085)
  assert approx_equal(gamma_array_test[1],-0.585563)
  assert approx_equal(gamma_array_test[2], 0.0)
  assert approx_equal(gamma_array_test[3], 0.0)

def test_sigma_prot():
  z_0 = scaling.sigma_prot_sq(0.0,1.0)
  z_0_theory = + 8.0*1.0*1.0 \
               + 5.0*6.0*6.0 \
               + 1.5*7.0*7.0 \
               + 1.2*8.0*8.0
  assert approx_equal(z_0,z_0_theory,eps=1e-0)

  d_star_sq = flex.double([0.0])
  z_0_array = scaling.get_sigma_prot_sq(d_star_sq,1.0)
  assert approx_equal(z_0_array[0],z_0)


## Testing isotropic wilson scaling
def finite_diffs_iso(p_scale=0.0,p_B_wilson=0.0,centric=False,h=0.0001):

  d_star_sq = flex.double(10,0.25)
  f_obs =  flex.double(10,1.0)
  centric_array = flex.bool(10,centric)
  sigma_f_obs = f_obs/10.0
  sigma_sq = flex.double(10,1.0)
  epsilon = flex.double(10,1.0)
  gamma =flex.double(10,0.0)

  stmp1 = scaling.wilson_total_nll(d_star_sq = d_star_sq,
                                    f_obs = f_obs,
                                    sigma_f_obs = sigma_f_obs,
                                    epsilon = epsilon,
                                    sigma_sq = sigma_sq,
                                    gamma_prot = gamma,
                                    centric = centric_array,
                                    p_scale = p_scale-h,
                                    p_B_wilson = p_B_wilson )

  stmp2 = scaling.wilson_total_nll(d_star_sq = d_star_sq,
                                    f_obs = f_obs,
                                    sigma_f_obs = sigma_f_obs,
                                    epsilon = epsilon,
                                    sigma_sq = sigma_sq,
                                    gamma_prot = gamma,
                                    centric = centric_array,
                                    p_scale = p_scale+h,
                                    p_B_wilson = p_B_wilson)

  s_grad_diff = (stmp1-stmp2)/(-2.0*h)

  btmp1 = scaling.wilson_total_nll(d_star_sq = d_star_sq,
                                    f_obs = f_obs,
                                    sigma_f_obs = sigma_f_obs,
                                    epsilon = epsilon,
                                    sigma_sq = sigma_sq,
                                    gamma_prot = gamma,
                                    centric = centric_array,
                                    p_scale = p_scale,
                                    p_B_wilson = p_B_wilson-h)

  btmp2 = scaling.wilson_total_nll(d_star_sq = d_star_sq,
                                    f_obs = f_obs,
                                    sigma_f_obs = sigma_f_obs,
                                    epsilon = epsilon,
                                    sigma_sq = sigma_sq,
                                    gamma_prot = gamma,
                                    centric = centric_array,
                                    p_scale = p_scale,
                                    p_B_wilson = p_B_wilson+h)

  b_grad_diff = (btmp1-btmp2)/(-2.0*h)

  grad  = scaling.wilson_total_nll_gradient(d_star_sq = d_star_sq,
                                             f_obs = f_obs,
                                             sigma_f_obs = sigma_f_obs,
                                             epsilon = epsilon,
                                             sigma_sq = sigma_sq,
                                             gamma_prot = gamma,
                                             centric = centric_array,
                                             p_scale = p_scale,
                                             p_B_wilson = p_B_wilson)
  assert approx_equal(s_grad_diff, grad[0])
  assert approx_equal(b_grad_diff, grad[1])


def test_likelihood_iso():
  d_star_sq = flex.double(10,0.250)
  f_obs = flex.double(10,1.0)
  sigma_f_obs = flex.double(10,0.0000)
  sigma_sq = flex.double(10,1.0)
  epsilon = flex.double(10,1.0)
  gamma = flex.double(10,0.0)
  centric = flex.bool(10,True)
  acentric = flex.bool(10,False)
  p_scale = 0.0
  p_B_wilson = 0.0



  centric_single_trans = scaling.wilson_single_nll(
    d_star_sq = d_star_sq[0],
    f_obs = f_obs[0],
    sigma_f_obs = sigma_f_obs[0],
    epsilon = epsilon[0],
    sigma_sq = sigma_sq[0],
    gamma_prot = gamma[0],
    centric = centric[0],
    p_scale = p_scale,
    p_B_wilson = p_B_wilson,
    transform = True)

  centric_single_no_trans = scaling.wilson_single_nll(
    d_star_sq = d_star_sq[0],
    f_obs = f_obs[0],
    sigma_f_obs = sigma_f_obs[0],
    epsilon = epsilon[0],
    sigma_sq = sigma_sq[0],
    gamma_prot = gamma[0],
    centric = centric[0],
    p_scale = 1.0,
    p_B_wilson = p_B_wilson,
    transform = False)

  assert approx_equal(centric_single_trans,  1.072364 ) ## from Mathematica
  assert approx_equal(centric_single_trans, centric_single_no_trans)

  acentric_single_trans = scaling.wilson_single_nll(
    d_star_sq = d_star_sq[0],
    f_obs = f_obs[0],
    sigma_f_obs = sigma_f_obs[0],
    epsilon = epsilon[0],
    sigma_sq = sigma_sq[0],
    gamma_prot = gamma[0],
    centric = acentric[0],
    p_scale = p_scale,
    p_B_wilson = p_B_wilson)

  acentric_single_no_trans  = scaling.wilson_single_nll(
    d_star_sq = d_star_sq[0],
    f_obs = f_obs[0],
    sigma_f_obs = sigma_f_obs[0],
    epsilon = epsilon[0],
    sigma_sq = sigma_sq[0],
    gamma_prot = gamma[0],
    centric = acentric[0],
    p_scale = 1.0,
    p_B_wilson =p_B_wilson,
    transform = False)

  assert approx_equal(acentric_single_trans, 0.306853) ## from Mathematica
  assert approx_equal(acentric_single_trans, acentric_single_no_trans)


  centric_total = scaling.wilson_total_nll(
    d_star_sq = d_star_sq,
    f_obs = f_obs,
    sigma_f_obs = sigma_f_obs,
    epsilon = epsilon,
    sigma_sq = sigma_sq,
    gamma_prot = gamma,
    centric = centric,
    p_scale = p_scale,
    p_B_wilson = p_B_wilson)

  acentric_total = scaling.wilson_total_nll(
    d_star_sq = d_star_sq,
    f_obs = f_obs,
    sigma_f_obs = sigma_f_obs,
    epsilon = epsilon,
    sigma_sq = sigma_sq,
    gamma_prot = gamma,
    centric = acentric,
    p_scale = p_scale,
    p_B_wilson = p_B_wilson)

  assert approx_equal(centric_total, centric_single_trans*10.0)
  assert approx_equal(acentric_total, acentric_single_trans*10.0)



def test_gradients_iso():
  ## Centrics
  finite_diffs_iso(p_scale=3.0,
               p_B_wilson=10.0,
               centric=True,h=0.000001)
  finite_diffs_iso(p_scale=-3.0,
               p_B_wilson=-10.0,
               centric=True,h=0.000001)
  finite_diffs_iso(p_scale=90.0,
               p_B_wilson=-10.0,
               centric=True,h=0.000001)
  finite_diffs_iso(p_scale=-90.0,
               p_B_wilson=10.0,
               centric=True,h=0.000001)
   ## Acentrics
  finite_diffs_iso(p_scale=3.0,
               p_B_wilson=10.0,
               centric=False,h=0.000001)
  finite_diffs_iso(p_scale=-3.0,
               p_B_wilson=-10.0,
               centric=False,h=0.000001)
  finite_diffs_iso(p_scale=90.0,
               p_B_wilson=-10.0,
               centric=True,h=0.000001)
  finite_diffs_iso(p_scale=-90.0,
               p_B_wilson=10.0,
               centric=True,h=0.000001)




## Testing anisotropic wilson scaling
def finite_diffs_aniso(p_scale,
                       u_star,
                       centric=False,
                       h=0.0001):
  d_star_sq = flex.double(2,0.25)
  f_obs =  flex.double(2,1.0)
  centric_array = flex.bool(2,centric)
  sigma_f_obs = f_obs/10.0
  sigma_sq = flex.double(2,1.0)
  epsilon = flex.double(2,1.0)
  gamma =flex.double(2,0.0)
  unit_cell = uctbx.unit_cell('20, 30, 40, 90.0, 90.0, 90.0')
  mi = flex.miller_index(((1,2,3), (1,2,3)))
  xs = crystal.symmetry((20,30,40), "P 2 2 2")
  ms = miller.set(xs, mi)

  nll_norm = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                             f_obs[0],
                                             sigma_f_obs[0],
                                             epsilon[0],
                                             sigma_sq[0],
                                             gamma[0],
                                             centric_array[0],
                                             p_scale,
                                             unit_cell,
                                             u_star)

  nll_scale = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                             f_obs[0],
                                             sigma_f_obs[0],
                                             epsilon[0],
                                             sigma_sq[0],
                                             gamma[0],
                                             centric_array[0],
                                             p_scale+h,
                                             unit_cell,
                                             u_star)
  u_star[0]+=h
  nll_u11 = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                             f_obs[0],
                                             sigma_f_obs[0],
                                             epsilon[0],
                                             sigma_sq[0],
                                             gamma[0],
                                             centric_array[0],
                                             p_scale,
                                             unit_cell,
                                             u_star)
  u_star[0]-=h
  u_star[1]+=h
  nll_u22 = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                             f_obs[0],
                                             sigma_f_obs[0],
                                             epsilon[0],
                                             sigma_sq[0],
                                             gamma[0],
                                             centric_array[0],
                                             p_scale,
                                             unit_cell,
                                             u_star)
  u_star[1]-=h
  u_star[2]+=h
  nll_u33 = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                             f_obs[0],
                                             sigma_f_obs[0],
                                             epsilon[0],
                                             sigma_sq[0],
                                             gamma[0],
                                             centric_array[0],
                                             p_scale,
                                             unit_cell,
                                             u_star)
  u_star[2]-=h
  u_star[3]+=h
  nll_u12 = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                             f_obs[0],
                                             sigma_f_obs[0],
                                             epsilon[0],
                                             sigma_sq[0],
                                             gamma[0],
                                             centric_array[0],
                                             p_scale,
                                             unit_cell,
                                             u_star)
  u_star[3]-=h
  u_star[4]+=h
  nll_u13 = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                             f_obs[0],
                                             sigma_f_obs[0],
                                             epsilon[0],
                                             sigma_sq[0],
                                             gamma[0],
                                             centric_array[0],
                                             p_scale,
                                             unit_cell,
                                             u_star)
  u_star[4]-=h
  u_star[5]+=h
  nll_u23 = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                            f_obs[0],
                                            sigma_f_obs[0],
                                            epsilon[0],
                                            sigma_sq[0],
                                            gamma[0],
                                            centric_array[0],
                                            p_scale,
                                            unit_cell,
                                            u_star)


  g = scaling.wilson_single_nll_aniso_gradient(ms.indices()[0],
                                               f_obs[0],
                                               sigma_f_obs[0],
                                               epsilon[0],
                                               sigma_sq[0],
                                               gamma[0],
                                               centric_array[0],
                                               p_scale,
                                               unit_cell,
                                               u_star)

  g2 = scaling.wilson_total_nll_aniso_gradient(ms.indices(),
                                               f_obs,
                                               sigma_f_obs,
                                               epsilon,
                                               sigma_sq,
                                               gamma,
                                               centric_array,
                                               p_scale,
                                               unit_cell,
                                               u_star)
  ds=(nll_norm-nll_scale)/-h
  du11=(nll_norm-nll_u11)/-h
  du22=(nll_norm-nll_u22)/-h
  du33=(nll_norm-nll_u33)/-h
  du12=(nll_norm-nll_u12)/-h
  du13=(nll_norm-nll_u13)/-h
  du23=(nll_norm-nll_u23)/-h
  assert approx_equal(ds,g[0]), (ds,g[0])
  assert approx_equal(du11,g[1]), (du11,g[1])
  assert approx_equal(du22,g[2])
  assert approx_equal(du33,g[3])
  assert approx_equal(du12,g[4])
  assert approx_equal(du13,g[5])
  assert approx_equal(du23,g[6])

  assert approx_equal(ds,g2[0]/2.0)
  assert approx_equal(du11,g2[1]/2.0)
  assert approx_equal(du22,g2[2]/2.0)
  assert approx_equal(du33,g2[3]/2.0)
  assert approx_equal(du12,g2[4]/2.0)
  assert approx_equal(du13,g2[5]/2.0)
  assert approx_equal(du23,g2[6]/2.0)


def test_likelihood_aniso():
  u_star = [0,0,0,0,0,0]
  d_star_sq = flex.double(2,0.25)
  f_obs =  flex.double(2,1.0)
  centric_array = flex.bool(2,True)
  sigma_f_obs = f_obs/10.0
  sigma_sq = flex.double(2,1.0)
  epsilon = flex.double(2,1.0)
  gamma =flex.double(2,0.0)
  unit_cell = uctbx.unit_cell('20, 30, 40, 90.0, 90.0, 90.0')
  mi = flex.miller_index(((1,2,3), (1,2,3)))
  xs = crystal.symmetry((20,30,40), "P 2 2 2")
  ms = miller.set(xs, mi)
  nll_centric_aniso = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                                      f_obs[0],
                                                      sigma_f_obs[0],
                                                      epsilon[0],
                                                      sigma_sq[0],
                                                      gamma[0],
                                                      centric_array[0],
                                                      0.0,
                                                      unit_cell,
                                                      u_star)


  assert approx_equal(nll_centric_aniso,  1.07239 ) ## from Mathematica
  nll_acentric_aniso = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                                       f_obs[0],
                                                       sigma_f_obs[0],
                                                       epsilon[0],
                                                       sigma_sq[0],
                                                       gamma[0],
                                                       centric_array[0],
                                                       0.0,
                                                       unit_cell,
                                                       u_star)
  centric_array = flex.bool(2,False)
  nll_acentric_aniso = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                                       f_obs[0],
                                                       sigma_f_obs[0],
                                                       epsilon[0],
                                                       sigma_sq[0],
                                                       gamma[0],
                                                       centric_array[0],
                                                       0.0,
                                                       unit_cell,
                                                       u_star)
  assert approx_equal(nll_acentric_aniso,0.306902 ) ## from Mathematica

  centric_array = flex.bool(2,True)
  u_star = [1,1,1,0,0,0]
  nll_centric_aniso = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                                      f_obs[0],
                                                      sigma_f_obs[0],
                                                      epsilon[0],
                                                      sigma_sq[0],
                                                      gamma[0],
                                                      centric_array[0],
                                                      0.0,
                                                      unit_cell,
                                                      u_star)
  assert approx_equal(nll_centric_aniso,  1.535008 ) ## from Mathematica
  centric_array = flex.bool(2,False)
  nll_acentric_aniso = scaling.wilson_single_nll_aniso(ms.indices()[0],
                                                      f_obs[0],
                                                      sigma_f_obs[0],
                                                      epsilon[0],
                                                      sigma_sq[0],
                                                      gamma[0],
                                                      centric_array[0],
                                                      0.0,
                                                      unit_cell,
                                                      u_star)
  assert approx_equal(nll_acentric_aniso,  0.900003 ) ## from Mathematica

  centric_array[1]=True
  nll_total_aniso = scaling.wilson_total_nll_aniso(ms.indices(),
                                                   f_obs,
                                                   sigma_f_obs,
                                                   epsilon,
                                                   sigma_sq,
                                                   gamma,
                                                   centric_array,
                                                   0.0,
                                                   unit_cell,
                                                   u_star)
  assert approx_equal(nll_total_aniso,  2.435011)


def test_grads_aniso():
  finite_diffs_aniso(0.0,[0.0,0.0,0.0,0.0,0.0,0.0],True, 0.0000001)
  finite_diffs_aniso(0.0,[0.0,0.0,0.0,2.0,0.0,0.0],False, 0.0000001)
  finite_diffs_aniso(0.0,[1.0,2.0,3.0,4.0,5.0,6.0],True, 0.0000001)
  finite_diffs_aniso(0.0,[1.0,2.0,3.0,4.0,5.0,6.0],False, 0.0000001)
  finite_diffs_aniso(-10.0,[1.0,2.0,3.0,4.0,5.0,6.0],True, 0.0000001)
  finite_diffs_aniso(-10.0,[1.0,2.0,3.0,4.0,5.0,6.0],False, 0.0000001)
  finite_diffs_aniso(10.0,[1.0,2.0,3.0,4.0,5.0,6.0],True, 0.0000001)
  finite_diffs_aniso(10.0,[1.0,2.0,3.0,4.0,5.0,6.0],False, 0.0000001)
  finite_diffs_aniso(10.0,[10.0,20.0,30.0,40.0,50.0,60.0],True, 0.0000001)
  finite_diffs_aniso(10.0,[10.0,20.0,30.0,40.0,50.0,60.0],False, 0.0000001)


## Testing relative scaling summats

class scaling_tester(object):
  def __init__(self):
    self.data_obs1 = flex.double(2,1.0)
    self.data_obs2 = flex.double(2,3.0)
    self.sigma_obs1 = flex.double(2,0.1)
    self.sigma_obs2 = flex.double(2,1)
    self.unit_cell = uctbx.unit_cell('20, 30, 40, 90.0, 90.0, 90.0')
    #mi = flex.miller_index(((1,2,3), (1,2,3)))
    self.mi = flex.miller_index(((1,2,3), (5,6,7)))
    self.xs = crystal.symmetry((20,30,40), "P 2 2 2")
    self.ms = miller.set(self.xs, self.mi)
    self.u = [1,2,3,4,5,6]
    self.p_scale = 0.40
    #self.u = [0,0,0,0,0,0]
    #self.p_scale = 0.00

    self.ls_i_wt = scaling.least_squares_on_i_wt(
      self.mi,
      self.data_obs1,
      self.sigma_obs1,
      self.data_obs2,
      self.sigma_obs2,
      self.p_scale,
      self.unit_cell,
      self.u)
    self.ls_i = scaling.least_squares_on_i(
      self.mi,
      self.data_obs1,
      self.sigma_obs1,
      self.data_obs2,
      self.sigma_obs2,
      self.p_scale,
      self.unit_cell,
      self.u)
    self.ls_f_wt = scaling.least_squares_on_f_wt(
      self.mi,
      self.data_obs1,
      self.sigma_obs1,
      self.data_obs2,
      self.sigma_obs2,
      self.p_scale,
      self.unit_cell,
      self.u)
    self.ls_f = scaling.least_squares_on_f(
      self.mi,
      self.data_obs1,
      self.sigma_obs1,
      self.data_obs2,
      self.sigma_obs2,
      self.p_scale,
      self.unit_cell,
      self.u)

    self.tst_ls_f_wt()
    self.tst_ls_i_wt()
    self.tst_ls_f()
    self.tst_ls_i()
    self.tst_hes_ls_i_wt()
    self.tst_hes_ls_f_wt()
    self.tst_hes_ls_i()
    self.tst_hes_ls_f()

  def tst_ls_i_wt(self, h=0.0000001):
    ## This function tests the gradients
    tmp = self.ls_i_wt.get_function()
    before = flex.double(7,tmp)
    after = flex.double(7,0)
    ## Test the pscale
    self.ls_i_wt.set_p_scale(self.p_scale+h)
    tmp = self.ls_i_wt.get_function()
    after[0]=tmp
    self.ls_i_wt.set_p_scale(self.p_scale)
    for ii in range(6):
      u_tmp=list(flex.double(self.u).deep_copy())
      u_tmp[ii]+=h
      self.ls_i_wt.set_u_rwgk(u_tmp)
      tmp = self.ls_i_wt.get_function()
      after[ii+1]=tmp
      self.ls_i_wt.set_u_rwgk(self.u)
      grads=self.ls_i_wt.get_gradient()
    f = max(1, flex.max(flex.abs(grads)))
    assert approx_equal(grads/f, (after-before)/h/f)

  def tst_ls_f_wt(self, h=0.0000001):
    ## This function tests the gradients
    tmp = self.ls_f_wt.get_function()
    before = flex.double(7,tmp)
    after = flex.double(7,0)
    ## Test the pscale
    self.ls_f_wt.set_p_scale(self.p_scale+h)
    tmp = self.ls_f_wt.get_function()
    after[0]=tmp
    self.ls_f_wt.set_p_scale(self.p_scale)
    for ii in range(6):
      u_tmp=list(flex.double(self.u).deep_copy())
      u_tmp[ii]+=h
      self.ls_f_wt.set_u_rwgk(u_tmp)
      tmp = self.ls_f_wt.get_function()
      after[ii+1]=tmp
      self.ls_f_wt.set_u_rwgk(self.u)
    grads=self.ls_f_wt.get_gradient()
    f = max(1, flex.max(flex.abs(grads)))
    assert approx_equal(grads/f, (after-before)/h/f)

  def tst_ls_f(self, h=0.0000001):
    ## This function tests the gradients
    tmp = self.ls_f.get_function()
    before = flex.double(7,tmp)
    after = flex.double(7,0)
    ## Test the pscale
    self.ls_f.set_p_scale(self.p_scale+h)
    tmp = self.ls_f.get_function()
    after[0]=tmp
    self.ls_f.set_p_scale(self.p_scale)
    for ii in range(6):
      u_tmp=list(flex.double(self.u).deep_copy())
      u_tmp[ii]+=h
      self.ls_f.set_u_rwgk(u_tmp)
      tmp = self.ls_f.get_function()
      after[ii+1]=tmp
      self.ls_f.set_u_rwgk(self.u)
    grads=self.ls_f.get_gradient()
    f = max(1, flex.max(flex.abs(grads)))
    assert approx_equal(grads/f, (after-before)/h/f)

  def tst_ls_i(self, h=0.0000001):
    ## This function tests the gradients
    tmp = self.ls_i.get_function()
    before = flex.double(7,tmp)
    after = flex.double(7,0)
    ## Test the pscale
    self.ls_i.set_p_scale(self.p_scale+h)
    tmp = self.ls_i.get_function()
    after[0]=tmp
    self.ls_i.set_p_scale(self.p_scale)
    #aniso tensor components
    for ii in range(6):
      u_tmp=list(flex.double(self.u).deep_copy())
      u_tmp[ii]+=h
      self.ls_i.set_u_rwgk(u_tmp)
      tmp = self.ls_i.get_function()
      after[ii+1]=tmp
      self.ls_i.set_u_rwgk(self.u)
    grads=self.ls_i.get_gradient()
    f = max(1, flex.max(flex.abs(grads)))
    assert approx_equal(grads/f, (after-before)/h/f)

  def tst_hes_ls_i_wt(self,h=0.0000001):

    hes_anal = self.ls_i_wt.hessian_as_packed_u()
    hes_anal=hes_anal.matrix_packed_u_as_symmetric()

    grads = self.ls_i_wt.get_gradient()

    self.ls_i_wt.set_p_scale(self.p_scale+h)
    tmp = self.ls_i_wt.get_gradient()
    tmp = list( (grads-tmp)/-h )
    tmp_hess=[]
    tmp_hess.append( tmp )
    self.ls_i_wt.set_p_scale(self.p_scale)

    for ii in range(6):
      u_tmp=list(flex.double(self.u).deep_copy())
      u_tmp[ii]+=h
      self.ls_i_wt.set_u_rwgk(u_tmp)
      tmp = self.ls_i_wt.get_gradient()
      tmp = (grads - tmp)/-h
      tmp_hess.append( list(tmp)  )
      self.ls_i_wt.set_u_rwgk(self.u)

    f = max(1, flex.max(flex.abs(hes_anal)))
    count=0
    for ii in range(7):
      for jj in range(7):
        assert approx_equal(tmp_hess[ii][jj]/f, hes_anal[count]/f)
        count+=1

  def tst_hes_ls_f_wt(self,h=0.0000001):

    hes_anal = self.ls_f_wt.hessian_as_packed_u()
    hes_anal=hes_anal.matrix_packed_u_as_symmetric()

    grads = self.ls_f_wt.get_gradient()

    self.ls_f_wt.set_p_scale(self.p_scale+h)
    tmp = self.ls_f_wt.get_gradient()
    tmp = list( (grads-tmp)/-h )
    tmp_hess=[]
    tmp_hess.append( tmp )
    self.ls_f_wt.set_p_scale(self.p_scale)

    for ii in range(6):
      u_tmp=list(flex.double(self.u).deep_copy())
      u_tmp[ii]+=h
      self.ls_f_wt.set_u_rwgk(u_tmp)
      tmp = self.ls_f_wt.get_gradient()
      tmp = (grads - tmp)/-h
      tmp_hess.append( list(tmp)  )
      self.ls_f_wt.set_u_rwgk(self.u)

    f = max(1, flex.max(flex.abs(hes_anal)))
    count=0
    for ii in range(7):
      for jj in range(7):
        assert approx_equal(tmp_hess[ii][jj]/f, hes_anal[count]/f)
        count+=1

  def tst_hes_ls_i(self,h=0.0000001):

    hes_anal = self.ls_i.hessian_as_packed_u()
    hes_anal=hes_anal.matrix_packed_u_as_symmetric()

    grads = self.ls_i.get_gradient()

    self.ls_i.set_p_scale(self.p_scale+h)
    tmp = self.ls_i.get_gradient()
    tmp = list( (grads-tmp)/-h )
    tmp_hess=[]
    tmp_hess.append( tmp )
    self.ls_i.set_p_scale(self.p_scale)

    for ii in range(6):
      u_tmp=list(flex.double(self.u).deep_copy())
      u_tmp[ii]+=h
      self.ls_i.set_u_rwgk(u_tmp)
      tmp = self.ls_i.get_gradient()
      tmp = (grads - tmp)/-h
      tmp_hess.append( list(tmp)  )
      self.ls_i.set_u_rwgk(self.u)

    count=0
    for ii in range(7):
      for jj in range(7):
        assert approx_equal(tmp_hess[ii][jj]/hes_anal[count], 1 , eps=1e-5)
        count+=1

  def tst_hes_ls_f(self,h=0.0000001):

    hes_anal = self.ls_f.hessian_as_packed_u()
    hes_anal=hes_anal.matrix_packed_u_as_symmetric()

    grads = self.ls_f.get_gradient()

    self.ls_f.set_p_scale(self.p_scale+h)
    tmp = self.ls_f.get_gradient()
    tmp = list( (grads-tmp)/-h )
    tmp_hess=[]
    tmp_hess.append( tmp )
    self.ls_f.set_p_scale(self.p_scale)

    for ii in range(6):
      u_tmp=list(flex.double(self.u).deep_copy())
      u_tmp[ii]+=h
      self.ls_f.set_u_rwgk(u_tmp)
      tmp = self.ls_f.get_gradient()
      tmp = (grads - tmp)/-h
      tmp_hess.append( list(tmp)  )
      self.ls_f.set_u_rwgk(self.u)

    count=0
    for ii in range(7):
      for jj in range(7):
        assert approx_equal(tmp_hess[ii][jj]/hes_anal[count], 1 , eps=1e-5)
        count+=1




def random_data(B_add=35,
                n_residues=585.0,
                d_min=3.5):
  unit_cell = uctbx.unit_cell( (81.0,  81.0,  61.0,  90.0,  90.0, 120.0) )
  xtal = crystal.symmetry(unit_cell, " P 3 ")
  ## In P3 I do not have to worry about centrics or reflections with different
  ## epsilons.
  miller_set = miller.build_set(
    crystal_symmetry = xtal,
    anomalous_flag = False,
    d_min = d_min)
  ## Now make an array with d_star_sq values
  d_star_sq = miller_set.d_spacings().data()
  d_star_sq = 1.0/(d_star_sq*d_star_sq)
  asu = {"H":8.0*n_residues*1.0,
         "C":5.0*n_residues*1.0,
         "N":1.5*n_residues*1.0,
         "O":1.2*n_residues*1.0}
  scat_info = absolute_scaling.scattering_information(
    asu_contents = asu,
    fraction_protein=1.0,
    fraction_nucleic=0.0)
  scat_info.scat_data(d_star_sq)
  gamma_prot = scat_info.gamma_tot
  sigma_prot = scat_info.sigma_tot_sq
  ## The number of residues is multriplied by the Z of the spacegroup
  protein_total = sigma_prot * (1.0+gamma_prot)
  ## add a B-value of 35 please
  protein_total = protein_total*flex.exp(-B_add*d_star_sq/2.0)
  ## Now that has been done,
  ## We can make random structure factors
  normalised_random_intensities = \
     random_transform.wilson_intensity_variate(protein_total.size())
  random_intensities = normalised_random_intensities*protein_total*math.exp(6)
  std_dev = random_intensities*5.0/100.0
  noise = random_transform.normal_variate(N=protein_total.size())
  noise = noise*std_dev
  random_intensities=noise+random_intensities
  ## STuff the arrays in the miller array
  miller_array = miller.array(miller_set,
                              data=random_intensities,
                              sigmas=std_dev)
  miller_array=miller_array.set_observation_type(
    xray.observation_types.intensity())
  miller_array = miller_array.f_sq_as_f()
  return (miller_array)


def test_scaling_on_random_data(B_add):
  miller_array = random_data(B_add,n_residues=100.0)
  scale_object_iso = absolute_scaling.ml_iso_absolute_scaling(
    miller_array,
    n_residues=100.0)

  ## compare the results please
  assert approx_equal(B_add, scale_object_iso.b_wilson, eps=5)

  scale_object_aniso = absolute_scaling.ml_aniso_absolute_scaling(
    miller_array,
    n_residues=100.0)

  assert approx_equal(B_add, scale_object_aniso.b_cart[0], eps=5)
  assert approx_equal(B_add, scale_object_aniso.b_cart[1], eps=5)
  assert approx_equal(B_add, scale_object_aniso.b_cart[2], eps=5)



def test_scattering_info():
  miller_array = random_data(35.0, d_min=2.5 )
  d_star_sq = miller_array.d_spacings().data()
  d_star_sq = 1.0/(d_star_sq*d_star_sq)

  asu = {"H":8.0*585.0,"C":5.0*585.0,"N":1.5*585.0, "O":1.2*585.0}
  scat_info = absolute_scaling.scattering_information(
    asu_contents = asu,
    fraction_protein=1.0,
    fraction_nucleic=0.0)
  scat_info.scat_data(d_star_sq)

  scat_info2 = absolute_scaling.scattering_information(
    n_residues=585.0)
  scat_info2.scat_data(d_star_sq)

  sigma_prot = scaling.get_sigma_prot_sq(d_star_sq,195.0*3.0)
  # Testing for consistency
  for ii in range(d_star_sq.size()):
    assert approx_equal(scat_info.sigma_tot_sq[ii],
                        scat_info2.sigma_tot_sq[ii],
                        eps=1e-03)
    assert approx_equal(scat_info.sigma_tot_sq[ii],
                        sigma_prot[ii],
                        eps=0.5)


def twin_the_data_and_analyse(twin_operator,twin_fraction=0.2):
  out_string = StringIO()

  miller_array = random_data(35).map_to_asu()
  miller_array = miller_array.f_as_f_sq()

  cb_op =  sgtbx.change_of_basis_op( twin_operator )

  miller_array_mod, miller_array_twin = miller_array.common_sets(
    miller_array.change_basis( cb_op ).map_to_asu() )
  twinned_miller = miller_array_mod.customized_copy(
    data = (1.0-twin_fraction)*miller_array_mod.data()
    + twin_fraction*miller_array_twin.data(),
    sigmas = flex.sqrt(
    flex.pow( ((1.0-twin_fraction)*miller_array_mod.sigmas()),2.0)+\
    flex.pow( ((twin_fraction)*miller_array_twin.sigmas()),2.0))
    )

  twinned_miller.set_observation_type( miller_array.observation_type())
  twin_anal_object = t_a.twin_analyses(twinned_miller,
                                       out=out_string,
                                       verbose=-100)

  index = twin_anal_object.twin_summary.twin_results.most_worrysome_twin_law

  assert approx_equal(
    twin_anal_object.twin_summary.twin_results.britton_alpha[index],
    twin_fraction,eps=0.1)

  assert approx_equal(twin_anal_object.twin_law_dependent_analyses[index].ml_murray_rust.estimated_alpha,
                      twin_fraction, eps=0.1)

  ## Untwinned data standards
  if twin_fraction==0:
    ## L-test
    assert approx_equal(twin_anal_object.l_test.mean_l, 0.50,eps=0.1)
    ## Wilson ratios
    assert approx_equal(twin_anal_object.twin_summary.twin_results.i_ratio,
                        2.00,
                        eps=0.1)
    ## H-test
    assert approx_equal(
      twin_anal_object.twin_law_dependent_analyses[index].h_test.mean_h,
      0.50,eps=0.1)


  ## Perfect twin standards
  if twin_fraction==0.5:
    assert approx_equal(twin_anal_object.l_test.mean_l, 0.375,eps=0.1)
    assert approx_equal(twin_anal_object.twin_summary.twin_results.i_ratio,
                        1.50,eps=0.1)
    assert approx_equal(
      twin_anal_object.twin_law_dependent_analyses[index].h_test.mean_h,
      0.00,eps=0.1)
  ## Just make sure we actually detect significant twinning
  if twin_fraction > 0.10:
    assert (twin_anal_object.twin_summary.twin_results.maha_l > 3.0)
  ## The patterson origin peak should be smallish ...
  assert (twin_anal_object.twin_summary.twin_results.patterson_p_value > 0.01)
  # and the brief test should be passed as well
  answer = t_a.twin_analyses_brief( twinned_miller,out=out_string,verbose=-100 )
  if twin_fraction > 0.10:
    assert answer is True


def test_twin_r_value(twin_operator):
  miller_array = random_data(35).map_to_asu()
  miller_array = miller_array.f_as_f_sq()

  for twin_fraction, expected_r_abs,expected_r_sq in zip(
     [0,0.1,0.2,0.3,0.4,0.5],
     [0.50,0.40,0.30,0.20,0.10,0.0],
     [0.333,0.213,0.120,0.0533,0.0133,0.00]):
    cb_op =  sgtbx.change_of_basis_op( twin_operator )

    miller_array_mod, miller_array_twin = miller_array.common_sets(
      miller_array.change_basis( cb_op ).map_to_asu() )
    twinned_miller = miller_array_mod.customized_copy(
      data = (1.0-twin_fraction)*miller_array_mod.data()
      + twin_fraction*miller_array_twin.data(),
      sigmas = flex.sqrt(
      flex.pow( ((1.0-twin_fraction)*miller_array_mod.sigmas()),2.0)+\
      flex.pow( ((twin_fraction)*miller_array_twin.sigmas()),2.0))
      )

    twinned_miller.set_observation_type( miller_array.observation_type())

    twin_r = scaling.twin_r( twinned_miller.indices(),
                             twinned_miller.data(),
                             twinned_miller.space_group(),
                             twinned_miller.anomalous_flag(),
                             cb_op.c().r().as_double()[0:9] )
    assert approx_equal(twin_r.r_abs_value(), expected_r_abs, 0.08)
    assert approx_equal(twin_r.r_sq_value(), expected_r_sq, 0.08)


def test_constant():
  # this is to make sure that the tmp_const in the class  symmetry_issues
  # does not result in any overflow problems
  math.log(1e-250)


def test_kernel_based_normalisation():
  miller_array = random_data(35.0, d_min=2.5 )
  normalizer = absolute_scaling.kernel_normalisation(
    miller_array, auto_kernel=True)
  z_values = normalizer.normalised_miller.data()/\
             normalizer.normalised_miller.epsilons().data().as_double()
  z_values = flex.mean(z_values)
  assert approx_equal(1.0,z_values,eps=0.05)

def test_ml_murray_rust():
  miller_array = random_data(35.0, d_min=4.5 )
  ml_mr_object = scaling.ml_murray_rust(
    miller_array.data(),
    miller_array.data(),
    miller_array.indices(),
    miller_array.space_group(),
    miller_array.anomalous_flag(),
    (1,0,0,0,1,0,0,0,1),
    6 )

  for ii in xrange(5,30):
    for jj in xrange(5,30):
       p1 = ml_mr_object.p_raw(ii/3.0, jj/3.0, 0.25)
       p2 = ml_mr_object.num_int(ii/3.0, 1e-13, jj/3.0, 1e-13, -5, 5,0.25, 20)
       assert approx_equal( p1, p2, eps=0.01)


if (__name__ == "__main__"):
  test_luts()
  test_ml_murray_rust()
  test_likelihood_iso()
  test_gradients_iso()
  test_gamma_prot()
  test_sigma_prot()
  test_likelihood_aniso()
  test_grads_aniso()
  test_scaling_on_random_data(10)
  test_scaling_on_random_data(20)
  test_scaling_on_random_data(40)
  test_scaling_on_random_data(70)
  test_scaling_on_random_data(80)
  scaling_tester()
  twin_the_data_and_analyse('h+k,-k,-l',0.00)
  twin_the_data_and_analyse('h+k,-k,-l',0.10)
  twin_the_data_and_analyse('h+k,-k,-l',0.20)
  twin_the_data_and_analyse('h+k,-k,-l',0.30)
  twin_the_data_and_analyse('h+k,-k,-l',0.50)
  test_scattering_info()
  test_kernel_based_normalisation()
  test_twin_r_value('h+k,-k,-l')
  test_constant()
  print "OK"
