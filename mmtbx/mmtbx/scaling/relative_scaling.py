from cctbx.array_family import flex
from mmtbx import scaling
from cctbx import uctbx
from cctbx import adptbx
from cctbx import sgtbx
from cctbx import eltbx
from scitbx.math import chebyshev_lsq
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
from libtbx.utils import Sorry
import scitbx.lbfgs
import math
import sys

class ls_relative_scaling(object):
  def __init__(self,
               miller_native,
               miller_derivative,
               start_values=None,
               mask=[1,1]):
    ## This mask allows one to refine only scale factor and only B values
    self.mask = mask ## multiplier for gradients of [scale factor, u tensor]

    ## make deep copies just to avoid any possible problems
    self.native = miller_native.deep_copy().set_observation_type(
      miller_native)

    if not self.native.is_real_array():
      raise Sorry("A real array is need for ls scaling")
    self.derivative = miller_derivative.deep_copy().set_observation_type(
      miller_derivative)
    if not self.derivative.is_real_array():
      raise Sorry("A real array is need for ls scaling")

    ## Convert into intensities
    if not self.native.is_xray_intensity_array():
      self.native = self.native.f_as_f_sq()
    if not self.derivative.is_xray_intensity_array():
      self.derivative = self.derivative.f_as_f_sq()

    ## Get the common sets
    self.native, self.derivative = self.native.map_to_asu().common_sets(
       self.derivative.map_to_asu() )

    ## Get the requiered information
    self.hkl = self.native.indices()

    self.i_nat =  self.native.data()
    self.sig_nat = self.native.sigmas()

    self.i_der = self.derivative.data()
    self.sig_der = self.derivative.sigmas()

    self.unit_cell = self.native.unit_cell()

    ## Symmetry related issues
    self.sg = self.native.space_group()
    self.adp_constraints = self.sg.adp_constraints(
      initialize_gradient_handling=True)
    self.dim_u = self.adp_constraints.n_independent_params
    ## Setup number of parameters
    assert self.dim_u()<=6
    ## Optimisation stuff
    self.x = flex.double(self.dim_u()+1, 0.0) ## B-values and scale factor!
    if start_values is not None:
      assert( start_values.size()==self.x.size() )
      self.x = start_values

    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)
    ## All done

    Vrwgk = math.pow(self.unit_cell.volume(),2.0/3.0)
    self.p_scale = self.x[0]
    self.u_star = self.unpack()
    self.u_star = list( flex.double(self.u_star) / Vrwgk )
    self.b_cart = adptbx.u_as_b(adptbx.u_star_as_u_cart(self.unit_cell,
                                                        self.u_star))
    self.u_cif = adptbx.u_star_as_u_cif(self.unit_cell,
                                        self.u_star)



  def pack(self,grad_tensor):
    grad_independent = [ grad_tensor[0]*float(self.mask[0]) ]+\
      list( float(self.mask[1])*
            flex.double(self.adp_constraints.independent_gradients(
              list(grad_tensor[1:])))
            )

    return grad_independent

  def unpack(self):
    u_tensor = self.adp_constraints.all_params( list(self.x[1:]) )
    return u_tensor

  def compute_functional_and_gradients(self):
    ## First please 'unpack' the U tensor
    u_full = self.unpack()
    ## Now we can compute the least squares score
    f = scaling.rel_scale_total_ls_target(
      self.hkl,
      self.i_nat,
      self.sig_nat,
      self.i_der,
      self.sig_der,
      self.x[0],
      self.unit_cell,
      u_full)
    ## Now compute the gradients
    g_full = scaling.rel_scale_total_ls_gradient(
      self.hkl,
      self.i_nat,
      self.sig_nat,
      self.i_der,
      self.sig_der,
      self.x[0],
      self.unit_cell,
      u_full)
    ## g_full has to be 'packed'
    g = self.pack( g_full )
    return f ,flex.double(g)


class ls_rel_scale_driver(object):
  def __init__(self,
               miller_native,
               miller_derivative):
    self.nat = miller_native
    self.der = miller_derivative

    ## Although this might look crude, (and probably is)
    ## It might be a goods idea to actually refine things in part
    ## for stability reasons. Although the datastes should be reasonably
    ## close together allready, due to absolute scaling performed idealy
    ## in advance, for more troublesome data, the following might be
    ## more robust.


    ## first refine the scale factor
    lsq_object = ls_relative_scaling(miller_native,
                                     miller_derivative,
                                     mask=[0,1] )

    ## fix the scale and refine the aniso U
    lsq_object = ls_relative_scaling(miller_native,
                                     miller_derivative,
                                     start_values=lsq_object.x,
                                     mask=[0,1] )
    ## refine everything at once
    lsq_object = ls_relative_scaling(miller_native,
                                     miller_derivative,
                                     start_values=lsq_object.x,
                                     mask=[1,1] )

    self.p_scale = lsq_object.p_scale
    self.b_cart = lsq_object.b_cart
    self.u_star = lsq_object.u_star


    ## very well, all done and set.
    ## apply the scaling on the data please and compute some r values

    tmp_nat, tmp_der = self.nat.common_sets(self.der)

    self.r_val_before = flex.sum( flex.abs(tmp_nat.data()-tmp_der.data()) )
    self.r_val_before /=flex.sum( flex.abs(tmp_nat.data()) )

    self.der = scaling.absolute_scaling.anisotropic_correction(
      self.der,self.p_scale,self.u_star )

    tmp_nat, tmp_der = self.nat.common_sets(self.der)

    self.r_val_after = flex.sum( flex.abs(tmp_nat.data()-tmp_der.data()) )
    self.r_val_after /=flex.sum( flex.abs(tmp_nat.data()) )

    ## All done

  def show(self, out=None):
    if out is None:
      out=sys.stdout

    print >> out
    print >> out, "R-value before LS scaling  : %5.3f"%(self.r_val_before)
    print >> out, "R-value after LS scaling   : %5.3f"%(self.r_val_after)
    print >> out
