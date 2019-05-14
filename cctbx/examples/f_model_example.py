from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import xray
from cctbx import sgtbx
from cctbx.development import random_structure as rs

def f_model_example():
  """ This example illustrates the use of the f_model class"""

  # make up some structure factors
  random_structure = rs.xray_structure(
    sgtbx.space_group_info( 'P1' ),
    elements=['C']*310,
    n_scatterers=310)

  # here f_atoms, f_mask and f_part are all the same
  # doens't make sense of course
  sfs = random_structure.structure_factors( False, 3.5,  ).f_calc()
  f_model_structure_factors =  xray.f_model_core_data( hkl = sfs.indices(),
                                             f_atoms= sfs.data(),
                                             f_mask = sfs.data(),
                                             unit_cell = sfs.unit_cell(),
                                             k_overall=1.0,
                                             u_star=(0,0,0,0,0,0),
                                             k_sol=1.0,
                                             u_sol=1.0,
                                             f_part=None,
                                             k_part=0,
                                             u_part=0 )

  #Resetting model parameters
  #
  # over all scale parameters; scale and aniso Ustar
  f_model_structure_factors.renew_overall_scale_parameters(
    3.0,(0,0,0,0,0,0) )
  # bulk solvent scale parameters; scale and B
  f_model_structure_factors.renew_bulk_solvent_scale_parameters(0.55,50.0)
  # partial structure scale parameters; scale and B
  f_model_structure_factors.renew_partial_structure_scale_parameters(0.05,25.0)

  # is is also possible to reset the values term by term ratrher then grouped
  f_model_structure_factors.ksol( 1.0 )
  f_model_structure_factors.usol( 3.0 )
  f_model_structure_factors.kpart( 1.0 )
  f_model_structure_factors.upart( 3.0 )
  f_model_structure_factors.koverall( 1.0 )
  f_model_structure_factors.ustar( (0,0,0,0,0,0) )
  # the cached arrays of various scale factor as updated automatically.

  # Obtaining the current parameters
  ksol = f_model_structure_factors.ksol()
  bsol = f_model_structure_factors.usol()
  kpart = f_model_structure_factors.kpart()
  bpart = f_model_structure_factors.upart()
  koverall = f_model_structure_factors.koverall()
  ustar = f_model_structure_factors.ustar()

  # Derivatives
  #
  #
  #derivates can be obtained accumalated over the full dataset
  # or on a struct term by term bases
  #what is needed in any case are one of the following arrays or floats
  # - d(target)/d(|F_model|)
  # - d(target)/d(Re[F_model]), d(target)/d(Im[F_model])
  #
  # which derivates are computed, is controlled by an array of gradient flags:
  # the order is (koverall, ustar, ksol, bsol, kpart, bpart )
  #
  gradient_flags = flex.bool([True,True,
                              True,True,
                              True,True])
  #
  # Use this function call is your target function returns
  # d(target)/d(|Fmodel|)
  dt_dabsfmodel = flex.double( sfs.data().size(), 1 )
  dtdall = f_model_structure_factors.d_target_d_all(
    dt_dabsfmodel,gradient_flags )

  # Use this function call is your target function returns
  # d(target)/d(Re[Fmodel]), d(target)/d(Im[Fmodel])
  dtda = dt_dabsfmodel
  dtdb = dt_dabsfmodel

  dtdall = f_model_structure_factors.d_target_d_all(
    dtda, dtdb, gradient_flags)
  # if desired (likely in c++, not in python) you can do it term by term.
  hkl_no=123
  dtdsingle = f_model_structure_factors.d_target_d_all(
    dtda[hkl_no], dtdb[hkl_no], hkl_no, gradient_flags)

  # the resulting gradients are delivered as a
  # 'f_model_derivative_holder'
  # it has methods both to set as well as to get items.
  #
  # getting the values
  tmp = dtdsingle.koverall()
  tmp = dtdsingle.ustar()
  tmp = dtdsingle.ksol()
  tmp = dtdsingle.usol()
  tmp = dtdsingle.kpart()
  tmp = dtdsingle.upart()
  #
  # if desired, you can set values as well
  dtdsingle.koverall(1)
  dtdsingle.ustar(  (1,1,1,1,1,1)  )
  dtdsingle.ksol(1)
  dtdsingle.usol(1)
  dtdsingle.kpart(1)
  dtdsingle.upart(1)



  #if desired, a selection can be made on the f_model object.
  #currently, only an integer selection is supported
  new_f_model_object = f_model_structure_factors.select( flex.int([1,2,3]) )
  assert  new_f_model_object.f_model().size()==3


if (__name__ == "__main__" ):
  f_model_example()
  print("OK")
