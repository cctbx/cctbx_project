from cctbx import xray
from cctbx.array_family import flex
import scitbx.lbfgs
from scitbx.python_utils.misc import adopt_init_args

class options:
  def __init__(self, site=00000, u_iso=00000, occupancy=00000):
    adopt_init_args(self, locals())

class lbfgs:

  def __init__(self, target_functor, options, xray_structure):
    adopt_init_args(self, locals())
    self.pack_parameters()
    self.first_target_value = None
    self.minimizer = scitbx.lbfgs.run(self)
    self.unpack_parameters()
    self.compute_target(compute_derivatives=00000)
    self.final_target_value = self.target_result.target()

  def pack_parameters(self):
    self.x = flex.double()
    self.n = xray.pack_parameters(
      None,
      self.xray_structure.scatterers(),
      self.x,
      self.options.site,
      self.options.u_iso,
      self.options.occupancy)

  def unpack_parameters(self):
    xray.unpack_parameters(
      None,
      self.xray_structure.space_group().order_z(),
      self.x, 0,
      self.xray_structure.scatterers(),
      self.options.site,
      self.options.u_iso,
      self.options.occupancy)

  def compute_target(self, compute_derivatives):
    self.f_calc_array = xray.structure_factors_direct(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs_array()).f_calc_array()
    self.target_result = self.target_functor(
      self.f_calc_array,
      compute_derivatives)

  def __call__(self):
    if (self.first_target_value != None):
      self.unpack_parameters()
    self.compute_target(compute_derivatives=0001)
    sf = xray.structure_factors_direct(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs_array(),
      d_target_d_f_calc=self.target_result.derivatives(),
      d_site_flag=self.options.site,
      d_u_iso_flag=self.options.u_iso,
      d_occupancy_flag=self.options.occupancy)
    self.g = flex.double()
    if (self.options.site):
      d_target_d_site = sf.d_target_d_site()
      self.xray_structure.apply_special_position_ops_d_target_d_site(
        d_target_d_site)
      self.g.append(d_target_d_site.as_double())
    if (self.options.u_iso):
      self.g.append(sf.d_target_d_u_iso())
    if (self.options.occupancy):
      self.g.append(sf.d_target_d_occupancy())
    if (self.first_target_value == None):
      self.pack_parameters()
    self.f = self.target_result.target()
    if (self.first_target_value == None):
      self.first_target_value = self.f
    return self.x, self.f, self.g
