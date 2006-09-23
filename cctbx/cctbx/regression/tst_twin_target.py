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

from scitbx.python_utils import random_transform
import random
import math
import sys
from cStringIO import StringIO
from cctbx.development import random_structure as rs

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
          assert approx_equal( (fdif_real-old_derivs[0][N_test])/fdif_real,0, eps=1e-1)
      if old_derivs[1][N_test]>0:
        #print  "Imag", N_test, fdif_imag,old_derivs[1][N_test], (fdif_imag-old_derivs[1][N_test])/old_derivs[1][N_test]
        if old_derivs[1][N_test]>2500:
          assert approx_equal( (fdif_imag-old_derivs[1][N_test])/fdif_imag,0, eps=1e-1)
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







def run():
  #tst_ls_on_i()
  tst_ls_on_f()


if (__name__ == "__main__"):
  run()
  print "OK"
