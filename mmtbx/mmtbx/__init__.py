import libtbx.forward_compatibility

import cctbx.array_family.flex
from mmtbx.refinement import print_statistics


class fmodels(object):
  def __init__(self, fmodel_xray = None,
                     fmodel_neutron = None,
                     model = None,
                     xray_scattering_dict = None,
                     neutron_scattering_dict = None,
                     log = None):
    self.fmodel_x = fmodel_xray
    self.fmodel_n = fmodel_neutron
    self.model = model
    self.xray_scattering_dict = xray_scattering_dict
    self.neutron_scattering_dict = neutron_scattering_dict
    self.log = log

  def fmodel_xray(self):
    if(self.fmodel_x is not None):
      if(self.fmodel_n is not None):
        self.fmodel_x.xray_structure.scattering_type_registry(custom_dict =
          self.xray_scattering_dict)
      assert self.fmodel_x.xray_structure is self.model.xray_structure
    return self.fmodel_x

  def fmodel_neutron(self):
    if(self.fmodel_n is not None):
      self.fmodel_n.xray_structure.scattering_type_registry(custom_dict =
        self.neutron_scattering_dict)
      assert self.fmodel_n.xray_structure is self.model.xray_structure
    return self.fmodel_n

  def update_xray_structure(self, xray_structure = None,
                                  update_f_calc  = None,
                                  update_f_mask  = None,
                                  force_update_f_mask = False):
    if(self.fmodel_x is not None):
      self.fmodel_xray().update_xray_structure(
        xray_structure = xray_structure,
        update_f_calc  = update_f_calc,
        update_f_mask  = update_f_mask,
        force_update_f_mask = force_update_f_mask)
      assert self.fmodel_x.xray_structure is self.model.xray_structure
    if(self.fmodel_n is not None):
      self.fmodel_neutron().update_xray_structure(
        xray_structure = xray_structure,
        update_f_calc  = update_f_calc,
        update_f_mask  = update_f_mask,
        force_update_f_mask = force_update_f_mask)
      assert self.fmodel_n.xray_structure is self.model.xray_structure

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
