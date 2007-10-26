import libtbx.forward_compatibility

import cctbx.array_family.flex
from mmtbx.refinement import print_statistics
from cctbx.array_family import flex


class fmodels(object):
  def __init__(self, fmodel_xray = None,
                     fmodel_neutron = None,
                     xray_scattering_dict = None,
                     neutron_scattering_dict = None,
                     log = None):
    self.fmodel_x = fmodel_xray
    self.fmodel_n = fmodel_neutron
    self.xray_scattering_dict = xray_scattering_dict
    self.neutron_scattering_dict = neutron_scattering_dict
    self.log = log
    self.create_target_functors()

  def pseudo_deep_copy(self):
    fmodel_n_dc = None
    if(self.fmodel_n is not None):
      fmodel_n_dc = self.fmodel_n.deep_copy()
    result = fmodels(fmodel_xray             = self.fmodel_x.deep_copy(),
                     fmodel_neutron          = fmodel_n_dc,
                     xray_scattering_dict    = self.xray_scattering_dict,
                     neutron_scattering_dict = self.neutron_scattering_dict,
                     log                     = self.log)
    result.update_xray_structure(xray_structure =
      self.fmodel_x.xray_structure.deep_copy_scatterers())
    return result

  def fmodel_xray(self, xray_structure = None):
    if(self.fmodel_x is not None):
      if(self.fmodel_n is not None):
        if(xray_structure is not None):
          self.fmodel_x.xray_structure = xray_structure
        self.fmodel_x.xray_structure.scattering_type_registry(custom_dict =
          self.xray_scattering_dict)
    return self.fmodel_x

  def fmodel_neutron(self, xray_structure = None):
    if(self.fmodel_n is not None):
      if(xray_structure is not None):
        self.fmodel_n.xray_structure = xray_structure
      self.fmodel_n.xray_structure.scattering_type_registry(custom_dict =
        self.neutron_scattering_dict)
    return self.fmodel_n

  def update_xray_structure(self, xray_structure = None,
                                  update_f_calc  = None,
                                  update_f_mask  = None,
                                  force_update_f_mask = False):
    if(self.fmodel_x is not None):
      self.fmodel_xray(xray_structure = xray_structure).update_xray_structure(
        xray_structure = xray_structure,
        update_f_calc  = update_f_calc,
        update_f_mask  = update_f_mask,
        force_update_f_mask = force_update_f_mask)
    if(self.fmodel_n is not None):
      self.fmodel_neutron(xray_structure=xray_structure).update_xray_structure(
        xray_structure = xray_structure,
        update_f_calc  = update_f_calc,
        update_f_mask  = update_f_mask,
        force_update_f_mask = force_update_f_mask)

  def show_short(self):
    if(self.fmodel_x is not None):
      prefix = ""
      if(self.fmodel_n is not None): prefix = "x-ray data"
      self.fmodel_xray().info().show_rfactors_targets_scales_overall(
        header = prefix, out = self.log)
    if(self.fmodel_n is not None):
      print >> self.log
      self.fmodel_neutron().info().show_rfactors_targets_scales_overall(
        header = "neutron data", out = self.log)

  def show_comprihensive(self, message = ""):
    print_statistics.make_sub_header("X-ray data", out = self.log)
    if(self.fmodel_x is not None):
      self.fmodel_xray().info().show_all(header = message, out = self.log)
    if(self.fmodel_n is not None):
      print_statistics.make_sub_header("Neutron data", out = self.log)
      self.fmodel_neutron().info().show_all(header = message, out = self.log)

  def update_bulk_solvent_and_scale(self, params = None, optimize_mask= False,
                                    force_update_f_mask = False):
    print_statistics.make_header("bulk solvent modeling and scaling",
      out = self.log)
    self.update_xray_structure(update_f_calc = True, update_f_mask = True,
      force_update_f_mask = force_update_f_mask)
    if(self.fmodel_x is not None):
      if(optimize_mask):
        self.fmodel_xray().optimize_mask_and_update_solvent_and_scale(
          params = params, out = self.log, verbose =-1)
      else:
        self.fmodel_xray().update_solvent_and_scale(params = params,
          out = self.log, verbose =-1)
    if(self.fmodel_n is not None):
      if(optimize_mask):
        self.fmodel_neutron().optimize_mask_and_update_solvent_and_scale(
          params = params, out = self.log, verbose =-1)
      else:
        self.fmodel_neutron().update_solvent_and_scale(params = params,
          out = self.log, verbose =-1)
    self.show_short()

  def remove_outliers(self):
    print_statistics.make_header("Outliers rejection", out = self.log)
    if(self.fmodel_x is not None):
      if(self.fmodel_n is not None):
        print_statistics.make_sub_header("x-ray data", out = self.log)
      self.fmodel_xray().remove_outliers(show = True, log = self.log)
    if(self.fmodel_n is not None):
      print_statistics.make_sub_header("neutron data", out = self.log)
      self.fmodel_neutron().remove_outliers(show = True, log = self.log)

  def create_target_functors(self):
    self.target_functor_xray = self.fmodel_xray().target_functor()
    self.target_functor_neutron = None
    if(self.fmodel_n is not None):
      self.target_functor_neutron = self.fmodel_neutron().target_functor()

  def target_functor_result_xray(self, compute_gradients):
    fmx = self.fmodel_xray()
    return self.target_functor_xray(compute_gradients = compute_gradients)

  def target_functor_result_neutron(self, compute_gradients):
    result = None
    if(self.fmodel_n is not None):
      fmn = self.fmodel_neutron()
      result = self.target_functor_neutron(compute_gradients=compute_gradients)
    return result

  def target_and_gradients(self, weights, compute_gradients, hd_selection=None,
        h_flag = None, u_iso_refinable_params = None):
    tfx = self.target_functor_result_xray
    tfn = self.target_functor_result_neutron
    class tg(object):
      def __init__(self, fmodels):
        self.fmodels = fmodels
        tfx_r = tfx(compute_gradients = compute_gradients)
        wx = weights.wx * weights.wx_scale
        self.target_work_xray = tfx_r.target_work()
        self.target_work_xray_weighted = self.target_work_xray * wx
        self.gradient_xray = None
        if(compute_gradients):
           sf = tfx_r.gradients_wrt_atomic_parameters(
             u_iso_refinable_params = u_iso_refinable_params).packed()
           # do not count grads for H or D:
           if(h_flag):
             sf_v3d = flex.vec3_double(sf)
             sf_v3d_sel = sf_v3d.set_selected(hd_selection, [0,0,0])
             sf = sf_v3d_sel.as_double()
           self.gradient_xray = sf
           self.gradient_xray_weighted = sf * wx
        if(fmodels.fmodel_neutron() is not None):
          wn = weights.wn * weights.wn_scale
          tfn_r = tfn(compute_gradients = compute_gradients)
          self.target_work_neutron = tfn_r.target_work()
          self.target_work_neutron_weighted = self.target_work_neutron * wn
          if(compute_gradients):
            sf = tfn_r.gradients_wrt_atomic_parameters(
              u_iso_refinable_params=u_iso_refinable_params).packed()
            self.gradient_neutron = sf
            self.gradient_neutron_weighted = sf * wn
      def target(self):
        if(self.fmodels.fmodel_neutron() is not None):
          result = (self.target_work_xray * weights.wx_scale + \
                    self.target_work_neutron * weights.wn_scale) * weights.wxn
        else: result = self.target_work_xray_weighted
        return result
      def gradients(self):
        result = None
        if(compute_gradients):
          if(self.fmodels.fmodel_neutron() is not None):
            result = (self.gradient_xray  * weights.wx_scale + \
                      self.gradient_neutron * weights.wn_scale) * weights.wxn
          else: result = self.gradient_xray_weighted
        return result
    result = tg(fmodels = self)
    return result
