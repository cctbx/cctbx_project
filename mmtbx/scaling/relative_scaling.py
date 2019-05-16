from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from mmtbx import scaling
from mmtbx.scaling import absolute_scaling
from cctbx import adptbx
from cctbx import crystal
from cctbx import miller
from libtbx.utils import Sorry
from scitbx.minimizers import newton_more_thuente_1994
from scitbx import matrix
import math
import sys
from six.moves import range

# 2009-04-15, cctbx svn rev. 8940:
#   there are no unit tests for this module, but is is used by these commands:
#     phenix.model_vs_data, phenix.refine, phenix.real_space_correlation,
#     phenix.twin_map_utils, phenix.xmanip

class refinery:

  def __init__(self,
               miller_native,
               miller_derivative,
               use_intensities=True,
               scale_weight=False,
               use_weights=False,
               mask=[1,1],
               start_values=None ):


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


    if use_intensities:
      if not self.native.is_xray_intensity_array():
        self.native = self.native.f_as_f_sq()
      if not self.derivative.is_xray_intensity_array():
        self.derivative = self.derivative.f_as_f_sq()
    if not use_intensities:
      if self.native.is_xray_intensity_array():
        self.native = self.native.f_sq_as_f()
      if self.derivative.is_xray_intensity_array():
        self.derivative = self.derivative.f_sq_as_f()

    ## Get the common sets
    self.native, self.derivative = self.native.map_to_asu().common_sets(
       self.derivative.map_to_asu() )

    ## Get the required information
    self.hkl = self.native.indices()

    self.i_or_f_nat =  self.native.data()
    self.sig_nat = self.native.sigmas()
    if self.sig_nat is None:
      self.sig_nat = self.i_or_f_nat*0 + 1

    self.i_or_f_der = self.derivative.data()
    self.sig_der = self.derivative.sigmas()
    if self.sig_der is None:
      self.sig_der = self.i_or_f_der*0+1

    self.unit_cell = self.native.unit_cell()

    # Modifiy the weights if required
    if not use_weights:
      self.sig_nat = self.sig_nat*0.0 + 1.0
      self.sig_der = self.sig_der*0.0


    ## Set up the minimiser 'cache'
    self.minimizer_object = None
    if use_intensities:
      if scale_weight:
        self.minimizer_object = scaling.least_squares_on_i_wt(
          self.hkl,
          self.i_or_f_nat,
          self.sig_nat,
          self.i_or_f_der,
          self.sig_der,
          0,
          self.unit_cell,
          [0,0,0,0,0,0])
      else :
        self.minimizer_object = scaling.least_squares_on_i(
          self.hkl,
          self.i_or_f_nat,
          self.sig_nat,
          self.i_or_f_der,
          self.sig_der,
          0,
          self.unit_cell,
          [0,0,0,0,0,0])
    else:
      if scale_weight:
        self.minimizer_object = scaling.least_squares_on_f_wt(
          self.hkl,
          self.i_or_f_nat,
          self.sig_nat,
          self.i_or_f_der,
          self.sig_der,
          0,
          self.unit_cell,
          [0,0,0,0,0,0])
      else :
        self.minimizer_object = scaling.least_squares_on_f(
          self.hkl,
          self.i_or_f_nat,
          self.sig_nat,
          self.i_or_f_der,
          self.sig_der,
          0,
          self.unit_cell,
          [0,0,0,0,0,0])

    ## Symmetry related issues
    self.sg = self.native.space_group()
    self.adp_constraints = self.sg.adp_constraints()
    self.dim_u = self.adp_constraints.n_independent_params
    ## Setup number of parameters
    assert self.dim_u()<=6
    ## Optimisation stuff
    x0 = flex.double(self.dim_u()+1, 0.0) ## B-values and scale factor!
    if start_values is not None:
      assert( start_values.size()==self.x.size() )
      x0 = start_values

    minimized = newton_more_thuente_1994(
      function=self, x0=x0, gtol=0.9e-6, eps_1=1.e-6, eps_2=1.e-6,
      matrix_symmetric_relative_epsilon=1.e-6)


    Vrwgk = math.pow(self.unit_cell.volume(),2.0/3.0)
    self.p_scale = minimized.x_star[0]
    self.u_star = self.unpack( minimized.x_star )
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
    return flex.double(grad_independent)

  def unpack(self,x):
    u_tensor = self.adp_constraints.all_params( list(x[1:]) )
    return u_tensor

  def functional(self, x):
    ## unpack the u-tensor
    u_full = self.unpack(x)
    ## place the params in the whatever
    self.minimizer_object.set_params(
      x[0],
      u_full)
    return self.minimizer_object.get_function()

  def gradients(self, x):
    u_full = self.unpack(x)
    self.minimizer_object.set_params(
      x[0],
      u_full)
    g_full = self.minimizer_object.get_gradient()
    g = self.pack( g_full )
    return g

  def hessian(self, x, eps=1.e-6):

    u_full = self.unpack(x)
    self.minimizer_object.set_params(
      x[0],
      u_full)
    result = self.minimizer_object.hessian_as_packed_u()
    result = result.matrix_packed_u_as_symmetric()
    result = self.hessian_transform(result,self.adp_constraints )
    return(result)

  ## This function is *only* for hessian with scale + utensor components
  def hessian_transform(self,
                        original_hessian,
                        adp_constraints):
    constraint_matrix_tensor = matrix.rec(
      adp_constraints.gradient_sum_matrix(),
      adp_constraints.gradient_sum_matrix().focus())

    hessian_matrix = matrix.rec( original_hessian,
                                 original_hessian.focus())
      ## now create an expanded matrix
    rows=adp_constraints.gradient_sum_matrix().focus()[0]+1
    columns=adp_constraints.gradient_sum_matrix().focus()[1]+1
    expanded_constraint_array = flex.double(rows*columns,0)
    count_new=0
    count_old=0
    for ii in range(rows):
      for jj in range(columns):
        if (ii>0):
          if (jj>0):
            expanded_constraint_array[count_new]=\
               constraint_matrix_tensor[count_old]
            count_old+=1
        count_new+=1
      ## place the first element please
    expanded_constraint_array[0]=1
    result=matrix.rec(  expanded_constraint_array,
                        (rows, columns) )
    #print result.mathematica_form()
    new_hessian = result *  hessian_matrix * result.transpose()
    result = flex.double(new_hessian)
    result.resize(flex.grid( new_hessian.n ) )
    return(result)



class ls_rel_scale_driver(object):
  def __init__(self,
               miller_native,
               miller_derivative,
               use_intensities=True,
               scale_weight=True,
               use_weights=True):
    self.native = miller_native.deep_copy().map_to_asu()
    self.derivative = miller_derivative.deep_copy().map_to_asu()

    lsq_object = refinery(self.native,
                          self.derivative,
                          use_intensities=use_intensities,
                          scale_weight=scale_weight,
                          use_weights=use_weights)


    self.p_scale = lsq_object.p_scale
    self.b_cart = lsq_object.b_cart
    self.u_star = lsq_object.u_star

    ## very well, all done and set.
    ## apply the scaling on the data please and compute some r values
    tmp_nat, tmp_der = self.native.common_sets(self.derivative)

    self.r_val_before = flex.sum( flex.abs(tmp_nat.data()-tmp_der.data()) )
    if flex.sum( flex.abs(tmp_nat.data()+tmp_der.data()) ) > 0:
      self.r_val_before /=flex.sum( flex.abs(tmp_nat.data()+tmp_der.data()) )/2.0

    self.derivative = absolute_scaling.anisotropic_correction(
      self.derivative,self.p_scale,self.u_star )

    self.scaled_original_derivative = self.derivative.deep_copy().set_observation_type(
      self.derivative ).map_to_asu()

    tmp_nat = self.native
    tmp_der = self.derivative

    tmp_nat, tmp_der = self.native.map_to_asu().common_sets(self.derivative.map_to_asu())
    self.r_val_after = flex.sum( flex.abs( tmp_nat.data()-
                                           tmp_der.data()   )
                               )
    if (flex.sum( flex.abs(tmp_nat.data()) ) +
                        flex.sum( flex.abs(tmp_der.data()) )) > 0:
      self.r_val_after /=(flex.sum( flex.abs(tmp_nat.data()) ) +
                          flex.sum( flex.abs(tmp_der.data()) ))/2.0

    self.native=tmp_nat
    self.derivative=tmp_der

    ## All done

  def show(self, out=None):
    if out is None:
      out=sys.stdout

    print(file=out)
    print("p_scale                    : %5.3f"%(self.p_scale), file=out)
    print("                            (%5.3f)"%(math.exp( self.p_scale ) ), file=out)
    print("B_cart trace               : %5.3f, %5.3f, %5.3f"%(
      self.b_cart[0],
      self.b_cart[1],
      self.b_cart[2]), file=out)
    print(file=out)
    print("R-value before LS scaling  : %5.3f"%(self.r_val_before), file=out)
    print("R-value after LS scaling   : %5.3f"%(self.r_val_after), file=out)
    print(file=out)




class local_scaling_driver(object):
  def __init__(self,
               miller_native,
               miller_derivative,
               local_scaling_dict,
               use_intensities=True,
               use_weights=False,
               max_depth=10,
               target_neighbours=1000,
               sphere=1,
               threshold=1.0,
               out=None):

    if out == None:
      out = sys.stdout


    self.native = miller_native.deep_copy().set_observation_type(
      miller_native).map_to_asu()
    self.derivative = miller_derivative.deep_copy().set_observation_type(
      miller_derivative).map_to_asu()

    assert self.native.observation_type() != None
    assert self.derivative.observation_type() != None

    self.native, self.derivative = self.native.common_sets(
      self.derivative)

    ## Here we change things into intensities or amplitudes as asked

    if use_intensities:
      if not self.native.is_xray_intensity_array():
        self.native = self.native.f_as_f_sq()
      if not self.derivative.is_xray_intensity_array():
        self.derivative = self.derivative.f_as_f_sq()

    if not use_intensities:
      if self.native.is_xray_intensity_array():
        self.native = self.native.f_sq_as_f()
      if self.derivative.is_xray_intensity_array():
        self.derivative = self.derivative.f_sq_as_f()


    ## In order to avoid problems with abseces due to lattice
    ## types, we transform the system to the primitive setting

    ## make new symmetry object
    self.nat_primset=self.native.change_basis(
     self.native.change_of_basis_op_to_niggli_cell()
     ).set_observation_type( self.native ).map_to_asu()

    self.der_primset=self.derivative.change_basis(
      self.derivative.change_of_basis_op_to_niggli_cell()
      ).set_observation_type(  self.derivative ).map_to_asu()

    ## Get the symmetry of the intensity group
    ## this is to make suree systematic absense are present in the hkl list
    intensity_group = self.nat_primset.space_group() \
      .build_derived_reflection_intensity_group(
      anomalous_flag=self.nat_primset.anomalous_flag() )

    ## We need to define a master set
    ## This is a miller array that encompasses
    ## both native and derivative *and* has systematic absences present.
    tmp_reso_lim = self.nat_primset.d_min()-0.1
    ## the 0.1 limit is a bit of a hack, but is needed
    ## to ensure that the master set is equal or slightkly larger than
    ## the working sets
    assert( tmp_reso_lim > 0 )

    self.master_set = miller.build_set(
      crystal_symmetry=crystal.symmetry(
      unit_cell=self.nat_primset.unit_cell(),
      space_group=intensity_group),
      anomalous_flag=self.nat_primset.anomalous_flag(),
      d_min=tmp_reso_lim)

    self.local_scaler=None
    self.sphere=sphere
    self.max_depth=max_depth
    self.target_neighbours=target_neighbours
    self.use_weights=use_weights
    self.threshold=threshold

    # Moment based scaling
    if local_scaling_dict['local_moment']:
      self.local_moment_scaling(out)

    # then it will be lsq based local scaling
    if local_scaling_dict['local_lsq']:
      self.local_lsq_scaling(out)

    # or nikonov scaling
    if local_scaling_dict['local_nikonov']:
      self.local_nikonov_scaling(out)


    scales=self.local_scaler.get_scales()
    stats=self.local_scaler.stats()

    print("Mean number of neighbours           : %8.3f"%(stats[2]), file=out)
    print("Minimum number of neighbours        : %8i"%(stats[0]), file=out)
    print("Maximum number of neighbours        : %8i"%(stats[1]), file=out)
    print(file=out)
    print("Mean local scale                    : %8.3f"%(
      flex.mean(scales) ), file=out)
    print("Standard deviation of local scale   : %8.3f"%(
      math.sqrt(   flex.mean(scales*scales)
                 - flex.mean(scales)*flex.mean(scales))), file=out)
    print("Minimum local scale                 : %8.3f"%(
      flex.min( scales ) ), file=out)
    print("Maximum local scale                 : %8.3f"%(
      flex.max( scales ) ), file=out)

    self.der_primset = self.der_primset.customized_copy(
       data = self.der_primset.data()*scales,
       sigmas = self.der_primset.sigmas()*scales
      ).set_observation_type( self.der_primset )

    ## We now have to transform the thing back please

    self.derivative = self.der_primset.change_basis(
        self.native.change_of_basis_op_to_niggli_cell().inverse()
      ).set_observation_type( self.der_primset ).map_to_asu()

    del self.der_primset
    del self.nat_primset

  def r_value(self,out):
    top = flex.abs(self.der_primset.data()-
                   self.nat_primset.data())
    bottom = flex.abs(self.der_primset.data() +
                      self.nat_primset.data())/2.0
    top=flex.sum(top)
    bottom=flex.sum(bottom)
    print("Current R value: %4.3f"%(top/bottom), file=out)


  def local_moment_scaling(self,out):
    print(file=out)
    print("Moment based local scaling", file=out)
    print("Maximum depth        : %8i"%(self.max_depth), file=out)
    print("Target neighbours    : %8i"%(self.target_neighbours), file=out)
    print("neighbourhood sphere : %8i"%(self.sphere), file=out)
    print(file=out)
    self.local_scaler = scaling.local_scaling_moment_based(
      hkl_master=self.master_set.indices(),
      hkl_sets=self.nat_primset.indices(),
      data_set_a=self.nat_primset.data(),
      sigma_set_a=self.nat_primset.sigmas(),
      data_set_b=self.der_primset.data(),
      sigma_set_b=self.der_primset.sigmas(),
      space_group=self.nat_primset.space_group(),
      anomalous_flag=self.nat_primset.anomalous_flag(),
      radius=self.sphere,
      depth=self.max_depth,
      target_ref=self.target_neighbours,
      use_experimental_sigmas=self.use_weights)


  def local_lsq_scaling(self,out):
    print(file=out)
    print("Least squares based local scaling", file=out)
    print("Maximum depth        : %8i"%(self.max_depth), file=out)
    print("Target neighbours    : %8i"%(self.target_neighbours), file=out)
    print("neighbourhood sphere : %8i"%(self.sphere), file=out)
    print(file=out)
    self.local_scaler = scaling.local_scaling_ls_based(
      hkl_master=self.master_set.indices(),
      hkl_sets=self.nat_primset.indices(),
      data_set_a=self.nat_primset.data(),
      sigma_set_a=self.nat_primset.sigmas(),
      data_set_b=self.der_primset.data(),
      sigma_set_b=self.der_primset.sigmas(),
      space_group=self.nat_primset.space_group(),
      anomalous_flag=self.nat_primset.anomalous_flag(),
      radius=self.sphere,
      depth=self.max_depth,
      target_ref=self.target_neighbours,
      use_experimental_sigmas=self.use_weights)

  def local_nikonov_scaling(self,out):
    print(file=out)
    print("Nikonev based local scaling", file=out)
    print("Maximum depth        : %8i"%(self.max_depth), file=out)
    print("Target neighbours    : %8i"%(self.target_neighbours), file=out)
    print("neighbourhood sphere : %8i"%(self.sphere), file=out)
    print(file=out)

    if self.der_primset.is_xray_intensity_array():
      raise Sorry(" For Nikonev target in local scaling, amplitudes must be used")
      assert (False)

    self.local_scaler = scaling.local_scaling_nikonov(
      hkl_master=self.master_set.indices(),
      hkl_sets=self.nat_primset.indices(),
      data_set_a=self.nat_primset.data(),
      data_set_b=self.der_primset.data(),
      epsilons=self.der_primset.epsilons().data().as_double(),
      centric=flex.bool(self.der_primset.centric_flags().data()),
      threshold=self.threshold,
      space_group=self.nat_primset.space_group(),
      anomalous_flag=self.nat_primset.anomalous_flag(),
      radius=self.sphere,
      depth=self.max_depth,
      target_ref=self.target_neighbours)
