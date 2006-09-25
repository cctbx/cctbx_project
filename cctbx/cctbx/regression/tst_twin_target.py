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
from scitbx import differential_evolution as de

from scitbx.python_utils import random_transform
import random
import math
import sys
from cStringIO import StringIO
from cctbx.development import random_structure as rs
import scitbx.lbfgs


random.seed(0)
flex.set_random_seed(0)


def tst_ls_on_i():
  tmp = rs.xray_structure(sgtbx.space_group_info( 'P4' ),
                          elements=['C']*310,
                          n_scatterers=310)

  sfs = tmp.structure_factors( False, 3.5  ).f_calc()
  f_mod = xray.f_model_core_data( hkl = sfs.indices(),
                        f_atoms= sfs.data(),
                        f_mask = sfs.data(),
                        unit_cell = sfs.unit_cell(),
                        k_overall=1.0,
                        u_star=(0,0,0,0,0,0),
                        k_sol=0.0,
                        u_sol=0.0,
                        f_part=sfs.data(),
                        k_part=0.0,
                        u_part=0.0 )

  target_evaluator = xray.least_squares_hemihedral_twinning_on_i(
    hkl_obs=sfs.indices(),
    i_obs=( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.0,
    w_obs=None,
    hkl_calc=sfs.indices(),
    space_group=sfs.space_group(),
    anomalous_flag=False,
    alpha=0.0,
    twin_law=[-1,0,0, 0,1,0, 0,0,-1]
    )

  target = target_evaluator.target( sfs.data() )
  # The target vaslue should be zero
  assert approx_equal( target, 0, eps=1e-6)

  # the derivatives as well
  derivs_complex = target_evaluator.d_target_d_fmodel( sfs.data()  )
  derivs_ab = target_evaluator.d_target_d_ab( sfs.data()  )
  for cmplx,da,db in zip(derivs_complex,
                         derivs_ab[0],
                         derivs_ab[1]):
    assert approx_equal( cmplx.real, da, eps=1e-5)
    assert approx_equal( cmplx.imag, -db, eps=1e-5)

  for alpha in flex.double(range(50))/100.0:
    #----------------------------------------------------------------
    # use fin diffs to check the derivatives to a and b
    old_target_evaluator = xray.least_squares_hemihedral_twinning_on_i(
      hkl_obs=sfs.indices(),
      i_obs=( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.1,
      w_obs=None,
      hkl_calc=sfs.indices(),
      space_group=sfs.space_group(),
      anomalous_flag=False,
      alpha=alpha,
      twin_law=[-1,0,0, 0,1,0, 0,0,-1]
      )
    old_target_value = old_target_evaluator.target( sfs.data() )
    old_derivs = old_target_evaluator.d_target_d_ab( sfs.data() )

    new_data =  sfs.data()
    h=0.0000001

    for N_test in xrange(sfs.data().size() ):
      ori = complex( sfs.data()[N_test] )
      #print "----------------"
      #print alpha
      #print sfs.indices()[N_test]
      #print sfs.data()[N_test]
      new_data[N_test] = ori+complex(h,0)
      new_target_value = old_target_evaluator.target( new_data )
      fdif_real = float((new_target_value-old_target_value)/h)
      new_data[N_test] = ori+complex(0,h)
      new_target_value = old_target_evaluator.target( new_data )
      fdif_imag = float( (new_target_value-old_target_value)/h )
      # only use 'large' first derivatives.
      if old_derivs[0][N_test]>0:
        #print "real", N_test, fdif_real,old_derivs[0][N_test], (fdif_real-old_derivs[0][N_test])/old_derivs[0][N_test]
        if old_derivs[0][N_test]>2500:
          assert approx_equal( (fdif_real-old_derivs[0][N_test])/fdif_real,0, eps=5e-1)
      if old_derivs[1][N_test]>0:
        #print  "Imag", N_test, fdif_imag,old_derivs[1][N_test], (fdif_imag-old_derivs[1][N_test])/old_derivs[1][N_test]
        if old_derivs[1][N_test]>2500:
          assert approx_equal( (fdif_imag-old_derivs[1][N_test])/fdif_imag,0, eps=5e-1)
      new_data[N_test] = ori

  #-------------------------------------
  # use fin diffs to test derivatives wrst alpha, the twin fraction
  h=0.0000000001

  target_evaluator = xray.least_squares_hemihedral_twinning_on_i(
    hkl_obs=sfs.indices(),
    i_obs=( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.0,
    w_obs=None,
    hkl_calc=sfs.indices(),
    space_group=sfs.space_group(),
    anomalous_flag=False,
    alpha=0,
    twin_law=[-1,0,0, 0,1,0, 0,0,-1]
    )

  tst_alpha = [0.1, 0.2, 0.3, 0.4, 0.5]

  for ii in tst_alpha:
    target_evaluator.alpha(ii)
    old_target = target_evaluator.target( sfs.data()*1.0 )
    target_evaluator.alpha( ii + h )
    new_target = target_evaluator.target( sfs.data()*1.0 )
    fd = (new_target-old_target)/h
    target_evaluator.alpha(ii)
    an = target_evaluator.d_target_d_alpha(sfs.data()*1.0)
    assert approx_equal( fd/an , 1.0, eps=1e-2 )






def tst_ls_on_f():
  tmp = rs.xray_structure(sgtbx.space_group_info( 'P4' ),
                          elements=['C']*310,
                          n_scatterers=310)

  sfs = tmp.structure_factors( False, 3.5  ).f_calc()
  f_mod = xray.f_model_core_data( hkl = sfs.indices(),
                        f_atoms= sfs.data(),
                        f_mask = sfs.data(),
                        unit_cell = sfs.unit_cell(),
                        k_overall=1.0,
                        u_star=(0,0,0,0,0,0),
                        k_sol=0.0,
                        u_sol=0.0,
                        f_part=sfs.data(),
                        k_part=0.0,
                        u_part=0.0 )

  target_evaluator = xray.least_squares_hemihedral_twinning_on_f(
    hkl_obs=sfs.indices(),
    f_obs=flex.sqrt( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.0,
    w_obs=None,
    hkl_calc=sfs.indices(),
    space_group=sfs.space_group(),
    anomalous_flag=False,
    alpha=0.0,
    twin_law=[-1,0,0, 0,1,0, 0,0,-1]
    )

  target = target_evaluator.target( sfs.data() )
  # The target vaslue should be zero
  assert approx_equal( target, 0, eps=1e-6)

  # the derivatives as well
  derivs_complex = target_evaluator.d_target_d_fmodel( sfs.data()  )
  derivs_ab = target_evaluator.d_target_d_ab( sfs.data()  )
  for cmplx,da,db in zip(derivs_complex,
                         derivs_ab[0],
                         derivs_ab[1]):
    assert approx_equal( cmplx.real, da, eps=1e-5)
    assert approx_equal( cmplx.imag, -db, eps=1e-5)

  for alpha in flex.double(range(50))/100.0:
    #----------------------------------------------------------------
    # use fin diffs to check the derivatives to a and b
    old_target_evaluator = xray.least_squares_hemihedral_twinning_on_f(
      hkl_obs=sfs.indices(),
      f_obs=flex.sqrt( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.1,
      w_obs=None,
      hkl_calc=sfs.indices(),
      space_group=sfs.space_group(),
      anomalous_flag=False,
      alpha=alpha,
      twin_law=[-1,0,0, 0,1,0, 0,0,-1]
      )
    old_target_value = old_target_evaluator.target( sfs.data() )
    old_derivs = old_target_evaluator.d_target_d_ab( sfs.data() )

    new_data =  sfs.data()
    h=0.0000001

    for N_test in xrange(sfs.data().size() ):
      ori = complex( sfs.data()[N_test] )
      #print "----------------"
      #print alpha
      #print sfs.indices()[N_test]
      #print sfs.data()[N_test]
      new_data[N_test] = ori+complex(h,0)
      new_target_value = old_target_evaluator.target( new_data )
      fdif_real = float((new_target_value-old_target_value)/h)
      new_data[N_test] = ori+complex(0,h)
      new_target_value = old_target_evaluator.target( new_data )
      fdif_imag = float( (new_target_value-old_target_value)/h )
      # only use 'large' first derivatives.
      if old_derivs[0][N_test]>0:
        #print "real", N_test, fdif_real,old_derivs[0][N_test], (fdif_real-old_derivs[0][N_test])/old_derivs[0][N_test]
        if old_derivs[0][N_test]>2500:
          assert approx_equal( (fdif_real-old_derivs[0][N_test])/fdif_real,0, eps=1e-1)
      if old_derivs[1][N_test]>0:
        #print  "Imag", N_test, fdif_imag,old_derivs[1][N_test], (fdif_imag-old_derivs[1][N_test])/old_derivs[1][N_test]
        if old_derivs[1][N_test]>2500:
          assert approx_equal( (fdif_imag-old_derivs[1][N_test])/fdif_imag,0, eps=1e-1)
      new_data[N_test] = ori

  #-------------------------------------
  # use fin diffs to test derivatives wrst alpha, the twin fraction
  h=0.0000000001

  target_evaluator = xray.least_squares_hemihedral_twinning_on_f(
    hkl_obs=sfs.indices(),
    f_obs=flex.sqrt( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.0,
    w_obs=None,
    hkl_calc=sfs.indices(),
    space_group=sfs.space_group(),
    anomalous_flag=False,
    alpha=0,
    twin_law=[-1,0,0, 0,1,0, 0,0,-1]
    )

  tst_alpha = [0.1, 0.2, 0.3, 0.4, 0.5]

  for ii in tst_alpha:
    target_evaluator.alpha(ii)
    old_target = target_evaluator.target( sfs.data()*1.0 )
    target_evaluator.alpha( ii + h )
    new_target = target_evaluator.target( sfs.data()*1.0 )
    fd = (new_target-old_target)/h
    target_evaluator.alpha(ii)
    an = target_evaluator.d_target_d_alpha(sfs.data()*1.0)
    #print ii, fd, an
    assert approx_equal( fd/an , 1.0, eps=1e-2 )


class find_max(object):
  def __init__(self,
               tot):
    self.to=tot
    self.n=2
    self.x=[]
    self.domain= [ (0,5), (0,5) ]
    self.eps=1e-16
    self.optim =  de.differential_evolution_optimizer(
      self,
      population_size=80,
      f=0.5,
      cr=0.9,
      n_cross=1,
      show_progress=True,
      show_progress_nth_cycle=1,eps=1e-12)


  def print_status(self, mean, min, best,count):
    None
    #print count, " \t:::  ", mean, min, list(best), list(
    #  self.to.d_log_p_d_f( best[0],best[1] ) )


  def target(self, vector):
    f1=vector[0]
    if f1<self.eps:
      f1=self.eps

    f2=vector[1]
    if f2<self.eps:
      f2=self.eps

    return -self.to.log_p(f1,f2)


class minim(object):
  def __init__(self,
               tot):
    self.tot=tot
    self.x = flex.double([3.0,2.0])
    scitbx.lbfgs.run(target_evaluator=self)

  def g(self,f1,f2,h=0.00001):
    fo = self.tot.log_p( f1, f2 )
    fp1= self.tot.log_p( f1+h, f2 )
    fp2= self.tot.log_p( f1, f2+h )
    tmp = [fp1-fo, fp2-fo]
    tmp = flex.double( tmp )/h
    return tmp

  def compute_functional_and_gradients(self):
    f = self.tot.log_p( self.x[0], self.x[1] )
    g = self.tot.d_log_p_d_f( self.x[0], self.x[1] )
    g1= self.g(self.x[0], self.x[1] )
    # test gradients wrst fin diffs
    div = int(math.log(abs(g[0])))
    if div<0:
      div=0
    div = math.exp(div)
    assert approx_equal( g[0]/div, g1[0]/div, eps=5e-1)
    div = int(math.log(abs(g[1])))
    if div<0:
      div=0
    div = math.exp(div)
    assert approx_equal( g[1]/div, g1[1]/div, eps=5e-1)
    #print list(g), list(g1)
    return -f,-flex.double(g)



def tst_snd_der(tot, f1, f2 , h=0.01):
  # see AMS55, 25.3.22 and 25.3.26
  #
  snd_der = tot.dd_log_p_dd_f(f1,f2)
  tmp0 = tot.log_p( f1,   f2   )
  tmp1 = tot.log_p( f1-h, f2   )
  tmp2 = tot.log_p( f1+h, f2   )

  ddf1 = (tmp1-2.0*tmp0 + tmp2)/(h*h)

  tmp0 = tot.log_p( f1,   f2   )
  tmp1 = tot.log_p( f1, f2-h   )
  tmp2 = tot.log_p( f1, f2+h   )

  ddf2 = (tmp1-2.0*tmp0+tmp2)/(h*h)

  div = int(math.log(abs( snd_der[0] )))
  if div<0:
    div=0
  div = math.exp(div)
  assert approx_equal( ddf1/div, snd_der[0]/div, eps=5e-2 )

  div = int(math.log(abs( snd_der[1] )))
  if div<0:
    div=0
  div = math.exp(div)
  assert approx_equal( ddf2/div, snd_der[1]/div, eps=5e-2 )

  tmp0 = tot.log_p( f1-h,   f2-h )
  tmp1 = tot.log_p( f1-h,   f2+h )
  tmp2 = tot.log_p( f1+h,   f2-h )
  tmp3 = tot.log_p( f1+h,   f2+h )

  tmp = (tmp3-tmp2-tmp1+tmp0)
  tmp = tmp/(h*h*4.0)

  div = int(math.log(abs( snd_der[2] )))
  if div<0:
    div=0
  div = math.exp(div)
  assert approx_equal( tmp/div, snd_der[2]/div, eps=5e-2 )

def tst_single_likelihood():
  io1=4.0
  io2=9.0
  so1=5.0/10.0
  so2=5.0/10.0
  fc1=2.0
  fc2=2.0
  a  =0.7
  b  =1.0-a*a
  centric=False
  tf=0.50

  tmp_1 = xray.single_twin_likelihood(io1,so1,io2,so2,
                                    fc1,fc2,1.0,1.0,
                                    centric,centric,
                                    a,b,
                                    tf, 14)
  tmp_1_de = find_max( tmp_1 )

  tmp_1_lbfgs = minim( tmp_1 )
  #Check if similar optimia are found
  assert approx_equal( tmp_1_de.x[0], tmp_1_lbfgs.x[0], eps=1e-3 )
  assert approx_equal( tmp_1_de.x[1], tmp_1_lbfgs.x[1], eps=1e-3 )

  # test second derivatives
  tst_snd_der(tmp_1, 3.0, 2.0 )
  tst_snd_der(tmp_1, 2.0, 3.0 )


def run():
  tst_ls_on_i()
  tst_ls_on_f()
  tst_single_likelihood()

if (__name__ == "__main__"):
  run()
  print "OK"
