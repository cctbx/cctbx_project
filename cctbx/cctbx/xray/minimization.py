import cctbx.xray.structure_factors
from cctbx.xray import ext
from cctbx.array_family import flex
import scitbx.lbfgs
from scitbx.python_utils.misc import adopt_init_args

class options:
  def __init__(self, site=00000, u_iso=00000, occupancy=00000):
    adopt_init_args(self, locals())

class lbfgs:

  def __init__(self, target_functor, options, xray_structure,
                     min_iterations=10, max_iterations=None):
    adopt_init_args(self, locals())
    self.structure_factors_from_scatterers = \
      cctbx.xray.structure_factors.from_scatterers(
        miller_set=self.target_functor.f_obs())
    self.structure_factor_gradients = \
      cctbx.xray.structure_factors.gradients(
        miller_set=self.target_functor.f_obs())
    self.pack_parameters()
    self.first_target_value = None
    self.minimizer = scitbx.lbfgs.run(
      self, min_iterations=min_iterations, max_iterations=max_iterations)
    self.unpack_parameters()
    self.compute_target(compute_derivatives=00000)
    self.final_target_value = self.target_result.target()

  def pack_parameters(self):
    self.x = flex.double()
    self.n = ext.pack_parameters(
      None,
      self.xray_structure.scatterers(),
      self.x,
      self.options.site,
      self.options.u_iso,
      self.options.occupancy)

  def unpack_parameters(self):
    ext.unpack_parameters(
      None,
      self.xray_structure.space_group().order_z(),
      self.x, 0,
      self.xray_structure.scatterers(),
      self.options.site,
      self.options.u_iso,
      self.options.occupancy)

  def compute_target(self, compute_derivatives):
    self.f_calc = self.structure_factors_from_scatterers(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs()).f_calc()
    self.target_result = self.target_functor(
      self.f_calc,
      compute_derivatives)

  def __call__(self):
    if (self.first_target_value is not None):
      self.unpack_parameters()
    self.compute_target(compute_derivatives=0001)
    sf = self.structure_factor_gradients(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs(),
      d_target_d_f_calc=self.target_result.derivatives(),
      gradient_flags=cctbx.xray.structure_factors.gradient_flags(
        site=self.options.site,
        u_iso=self.options.u_iso,
        occupancy=self.options.occupancy),
      direct=0001)
    self.g = flex.double()
    if (self.options.site):
      d_target_d_site = sf.d_target_d_site_frac()
      self.xray_structure.apply_special_position_ops_d_target_d_site(
        d_target_d_site)
      self.g.append(d_target_d_site.as_double())
    if (self.options.u_iso):
      self.g.append(sf.d_target_d_u_iso())
    if (self.options.occupancy):
      self.g.append(sf.d_target_d_occupancy())
    if (self.first_target_value is None):
      self.pack_parameters()
    self.f = self.target_result.target()
    if (self.first_target_value is None):
      self.first_target_value = self.f
    return self.x, self.f, self.g
