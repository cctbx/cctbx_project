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
  f_model_structure_factors =  xray.f_model( hkl = sfs.indices(),
                                             f_atoms= sfs.data(),
                                             f_mask = sfs.data(),
                                             unit_cell = sfs.unit_cell(),
                                             k_overall=1.0,
                                             u_star=(0,0,0,0,0,0),
                                             k_sol=1.0,
                                             b_sol=1.0,
                                             f_part=None,
                                             k_part=0,
                                             b_part=0 )

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
  # take care however that than a 'refresh' command is needed to redo some
  # computations. If this is not done, the subsequent gradients will not
  # be correct!
  f_model_structure_factors.ksol( 1.0 )
  f_model_structure_factors.bsol( 3.0 )
  f_model_structure_factors.kpart( 1.0 )
  f_model_structure_factors.bpart( 3.0 )
  f_model_structure_factors.koverall( 1.0 )
  f_model_structure_factors.ustar( (0,0,0,0,0,0) )
  # this is the refresh command
  f_model_structure_factors.refresh()

  # Obtaining the current parameters
  ksol = f_model_structure_factors.ksol()
  bsol = f_model_structure_factors.bsol()
  kpart = f_model_structure_factors.kpart()
  bpart = f_model_structure_factors.bpart()
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

  # if desired (not likely) you can do it term by term.
  hkl_no=123
  dtdsingle = f_model_structure_factors.d_target_d_all(
    dtda[hkl_no], dtdb[hkl_no], hkl_no, gradient_flags)

  # the resulting gradients are delivered as a
  # 'f_model_derivatiove_holder'
  # it has methods both to set as well as to get items.
  #
  # getting the values
  tmp = dtdsingle.koverall()
  tmp = dtdsingle.ustar()
  tmp = dtdsingle.ksol()
  tmp = dtdsingle.bsol()
  tmp = dtdsingle.kpart()
  tmp = dtdsingle.bpart()
  #
  # if desired, you can set values as well
  dtdsingle.koverall(1)
  dtdsingle.ustar(  (1,1,1,1,1,1)  )
  dtdsingle.ksol(1)
  dtdsingle.bsol(1)
  dtdsingle.kpart(1)
  dtdsingle.bpart(1)













if (__name__ == "__main__" ):
  f_model_example()
