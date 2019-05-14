from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import xray
from cctbx import sgtbx
from libtbx.test_utils import approx_equal

import random
import math
from cctbx.development import random_structure as rs

random.seed(0)
flex.set_random_seed(0)

def compare_all(tmp1,tmp2):
  assert approx_equal( tmp1.ksol(),tmp2.ksol() )
  assert approx_equal( tmp1.usol(),tmp2.usol() )
  assert approx_equal( tmp1.kpart(),tmp2.kpart() )
  assert approx_equal( tmp1.upart(),tmp2.upart() )
  assert approx_equal( tmp1.koverall(),tmp2.koverall() )
  assert approx_equal( tmp1.ustar(), tmp2.ustar() )

def tst_f_model_derivative_holder():
  dinfo = xray.f_model_core_data_derivative_holder()

  dinfo.ksol(1)
  dinfo.usol(1)
  dinfo.kpart(1)
  dinfo.upart(1)
  dinfo.koverall(1)
  dinfo.ustar( [1,1,1,1,1,1] )

  dinfo.accumulate( dinfo )
  assert( dinfo.ksol()==2 )
  assert( dinfo.usol()==2 )
  assert( dinfo.kpart()==2 )
  assert( dinfo.upart()==2 )
  assert( dinfo.koverall()==2 )
  assert( dinfo.ustar()==(2,2,2,2,2,2) )


def tst_f_model():
  tmp = rs.xray_structure(sgtbx.space_group_info( 'P1' ),
                          elements=['C']*310,
                          n_scatterers=310)

  sfs = tmp.structure_factors( False, 3.5,  ).f_calc()
  f_mod = xray.f_model_core_data( hkl = sfs.indices(),
                        f_atoms= sfs.data(),
                        f_mask = sfs.data(),
                        unit_cell = sfs.unit_cell(),
                        k_overall=1.0,
                        u_star=(0,0,0,0,0,0),
                        k_sol=1.0,
                        u_sol=0.1,
                        f_part=sfs.data(),
                        k_part=1.0,
                        u_part=0.1 )
  f_mod.refresh()

  data1 = sfs.data()[123]
  hkl1 =  sfs.indices()[123]

  fbulk = f_mod.f_bulk()[123]
  fatoms = f_mod.f_atoms()[123]
  fpart = f_mod.f_part()[123]
  fmod = f_mod.f_model()[123]
  # we have unit scale now
  assert approx_equal( fmod, fbulk+fatoms+fpart )
  # get derivatives please
  #
  # make a mock target function: T=A^2 + B^2
  # dT/dF=2F
  # dT/dA=2A
  # dT/dB=2B
  f_model_data_complex = f_mod.f_model()
  fm = flex.double()
  a = flex.double()
  b = flex.double()

  for cmplx in f_model_data_complex:
    tmp_a=cmplx.real
    tmp_b=cmplx.imag
    fm.append( math.sqrt(tmp_a*tmp_a + tmp_b*tmp_b) )
    a.append( tmp_a )
    b.append( tmp_b )

  dtdf = 2.0*fm
  dtda = 2.0*a
  dtdb = 2.0*b

  gradient_flags=flex.bool([True,True,
                            True,True,
                            True,True])

  grads_f = f_mod.d_target_d_all(dtdf,gradient_flags)
  grads_ab = f_mod.d_target_d_all(dtda, dtdb,gradient_flags)
  compare_all(grads_ab,grads_f)

  f_mod = xray.f_model_core_data( hkl = sfs.indices(),
                        f_atoms= sfs.data(),
                        f_mask = sfs.data()*0.0,
                        unit_cell = sfs.unit_cell(),
                        k_overall=1.0,
                        u_star=(0,0,0,0,0,0),
                        k_sol=0.0,
                        u_sol=0.0,
                        f_part=sfs.data()*0.0,
                        k_part=0.0,
                        u_part=0.0 )
  f_mod.refresh()
  grad_123 = f_mod.d_target_d_all(0,1,123, gradient_flags)
  assert approx_equal( grad_123.koverall(),
                       f_mod.f_model()[123].imag/f_mod.koverall() )
  tps=19.7392088
  h=hkl1[0]
  k=hkl1[1]
  l=hkl1[2]
  assert approx_equal( -2*tps*h*h*f_mod.f_model()[123].imag,
                       grad_123.ustar()[0] )
  assert approx_equal( -2*tps*k*k*f_mod.f_model()[123].imag,
                       grad_123.ustar()[1] )
  assert approx_equal( -2*tps*l*l*f_mod.f_model()[123].imag,
                       grad_123.ustar()[2] )
  assert approx_equal( -4*tps*h*k*f_mod.f_model()[123].imag,
                       grad_123.ustar()[3] )
  assert approx_equal( -4*tps*h*l*f_mod.f_model()[123].imag,
                       grad_123.ustar()[4] )
  assert approx_equal( -4*tps*k*l*f_mod.f_model()[123].imag,
                       grad_123.ustar()[5] )





  grad_123 = f_mod.d_target_d_all(1,0,123, gradient_flags)
  assert approx_equal( grad_123.koverall(),
                       f_mod.f_model()[123].real/f_mod.koverall() )
  tps=19.7392088
  h=hkl1[0]
  k=hkl1[1]
  l=hkl1[2]
  assert approx_equal( -2*tps*h*h*f_mod.f_model()[123].real,
                       grad_123.ustar()[0] )
  assert approx_equal( -2*tps*k*k*f_mod.f_model()[123].real,
                       grad_123.ustar()[1] )
  assert approx_equal( -2*tps*l*l*f_mod.f_model()[123].real,
                       grad_123.ustar()[2] )
  assert approx_equal( -4*tps*h*k*f_mod.f_model()[123].real,
                       grad_123.ustar()[3] )
  assert approx_equal( -4*tps*h*l*f_mod.f_model()[123].real,
                       grad_123.ustar()[4] )
  assert approx_equal( -4*tps*k*l*f_mod.f_model()[123].real,
                       grad_123.ustar()[5] )





  oldfm = xray.f_model_core_data( hkl = sfs.indices(),
                        f_atoms= sfs.data(),
                        f_mask = sfs.data()*1.0,
                        unit_cell = sfs.unit_cell(),
                        k_overall=1.0,
                        u_star=(0,0,0,0,0,0),
                        k_sol=1.0,
                        u_sol=0.1,
                        f_part=sfs.data()*1.0,
                        k_part=1.0,
                        u_part=0.1 )
  h=0.0001


  newfm = xray.f_model_core_data( hkl = sfs.indices(),
                        f_atoms= sfs.data(),
                        f_mask = sfs.data()*1.0,
                        unit_cell = sfs.unit_cell(),
                        k_overall=1.0,
                        u_star=(0,0,0,0,0,0),
                        k_sol=1.0,
                        u_sol=0.1,
                        f_part=sfs.data()*1.0,
                        k_part=1.0,
                        u_part=0.1 )


  newfm.renew_overall_scale_parameters(1.0+h, (0,0,0,0,0,0))
  newfm.refresh()

  tmp = oldfm.d_target_d_all(1,0,123,gradient_flags).koverall()
  tmp_d = ((oldfm.f_model()[123]-newfm.f_model()[123])/(-h)).real
  assert approx_equal( tmp,tmp_d,eps=1e-4 )
  newfm.renew_overall_scale_parameters(1, (0,0,0,0,0,0))


  newfm.renew_bulk_solvent_scale_parameters(1+h,0.1)
  tmp = oldfm.d_target_d_all(1,0,123,gradient_flags).ksol()
  tmp_d = ((oldfm.f_model()[123]-newfm.f_model()[123])/(-h)).real
  assert approx_equal( tmp,tmp_d,eps=1e-2 )
  newfm.renew_bulk_solvent_scale_parameters(1,0.1)

  newfm.renew_bulk_solvent_scale_parameters(1.0,0.1+h)
  tmp = oldfm.d_target_d_all(1,0,123,gradient_flags).usol()
  tmp_d = ((oldfm.f_model()[123]-newfm.f_model()[123])/(-h)).real
  assert approx_equal( tmp,tmp_d,eps=1e-2 )
  newfm.renew_bulk_solvent_scale_parameters(1.0,0.1)

  newfm.renew_bulk_solvent_scale_parameters(1+h,0.1)
  tmp = oldfm.d_target_d_all(0,1,123,gradient_flags).ksol()
  tmp_d = ((oldfm.f_model()[123]-newfm.f_model()[123])/(-h)).imag
  assert approx_equal( tmp,tmp_d,eps=1e-2 )
  newfm.renew_bulk_solvent_scale_parameters(1,0.1)

  newfm.renew_bulk_solvent_scale_parameters(1.0,0.1+h)
  tmp = oldfm.d_target_d_all(0,1,123,gradient_flags).usol()
  tmp_d = ((oldfm.f_model()[123]-newfm.f_model()[123])/(-h)).imag
  assert approx_equal( tmp,tmp_d,eps=1e-2 )
  newfm.renew_bulk_solvent_scale_parameters(1.0,0.1)


  newfm.renew_partial_structure_scale_parameters(1+h,0.1)
  tmp = oldfm.d_target_d_all(1,0,123,gradient_flags).kpart()
  tmp_d = ((oldfm.f_model()[123]-newfm.f_model()[123])/(-h)).real
  assert approx_equal( tmp,tmp_d,eps=1e-2 )
  newfm.renew_bulk_solvent_scale_parameters(1,0.1)

  newfm.renew_partial_structure_scale_parameters(1.0,0.1+h)
  tmp = oldfm.d_target_d_all(1,0,123,gradient_flags).upart()
  tmp_d = ((oldfm.f_model()[123]-newfm.f_model()[123])/(-h)).real
  assert approx_equal( tmp,tmp_d,eps=1e-2 )
  newfm.renew_bulk_solvent_scale_parameters(1.0,0.1)

  newfm.renew_partial_structure_scale_parameters(1+h,0.1)
  tmp = oldfm.d_target_d_all(0,1,123,gradient_flags).kpart()
  tmp_d = ((oldfm.f_model()[123]-newfm.f_model()[123])/(-h)).imag
  assert approx_equal( tmp,tmp_d,eps=1e-2 )
  newfm.renew_bulk_solvent_scale_parameters(1,0.1)

  newfm.renew_partial_structure_scale_parameters(1.0,0.1+h)
  tmp = oldfm.d_target_d_all(0,1,123,gradient_flags).upart()
  tmp_d = ((oldfm.f_model()[123]-newfm.f_model()[123])/(-h)).imag
  assert approx_equal( tmp,tmp_d,eps=1e-2 )
  newfm.renew_bulk_solvent_scale_parameters(1.0,0.1)



def run():
  tst_f_model()
  tst_f_model_derivative_holder()


if (__name__ == "__main__"):
  ()
  print("OK")
