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

  def resolution_filter(self, d_min):
    fmodel_n_dc = None
    if(self.fmodel_n is not None):
      fmodel_n_dc = self.fmodel_n.resolution_filter(d_min = d_min)
    result = fmodels(
      fmodel_xray             = self.fmodel_x.resolution_filter(d_min = d_min),
      fmodel_neutron          = fmodel_n_dc,
      xray_scattering_dict    = self.xray_scattering_dict,
      neutron_scattering_dict = self.neutron_scattering_dict,
      log                     = self.log)
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
      self.fmodel_x = self.fmodel_xray().remove_outliers(
        show = True, log = self.log)
    if(self.fmodel_n is not None):
      print_statistics.make_sub_header("neutron data", out = self.log)
      self.fmodel_n = self.fmodel_neutron().remove_outliers(
        show = True, log = self.log)

  def show_targets(self, log, text=""):
    prefix_x = ""
    if(self.fmodel_n is not None):
      prefix_x = "xray"
    self.fmodel_xray().info().show_targets(out = log, text = prefix_x+" "+text)
    if(self.fmodel_n is not None):
      self.fmodel_neutron().info().show_targets(out= log, text="neutron "+text)

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

  def target_and_gradients(self, weights, compute_gradients,
        u_iso_refinable_params = None, occupancy = False):
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
           if(occupancy):
             sf = tfx_r.gradients_wrt_atomic_parameters(occupancy = occupancy)
           else:
             sf = tfx_r.gradients_wrt_atomic_parameters(
               u_iso_refinable_params = u_iso_refinable_params).packed()
           self.gradient_xray = sf
           self.gradient_xray_weighted = sf * wx
        if(fmodels.fmodel_neutron() is not None):
          wn = weights.wn * weights.wn_scale
          tfn_r = tfn(compute_gradients = compute_gradients)
          self.target_work_neutron = tfn_r.target_work()
          self.target_work_neutron_weighted = self.target_work_neutron * wn
          if(compute_gradients):
            if(occupancy):
              sf = tfn_r.gradients_wrt_atomic_parameters(occupancy = occupancy)
            else:
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

class map_names(object):
  def __init__(self, map_name_string):
    k,n = None,None
    ml_map = False
    s = map_name_string.lower()
    self.k = k
    self.n = n
    self.ml_map = ml_map
    self.anomalous = False
    if(s.count('anom')):
      self.anomalous = True
    else:
      s = s.replace(" ","")
      s = s.replace("*","")
      if(s.count('_')==1): self.error(map_name_string)
      found = False
      for item in ['fobs','fob','fo']:
        if(s.count(item)==1):
          s = s.replace(item,"_")
          found = True
      if(not found): self.error(map_name_string)
      found = False
      for item in ['fcalc','fcal','fc', 'fmodel','fmod','fm']:
        if(s.count(item)==1):
          if(not s.endswith(item)): self.error(map_name_string)
          s = s.replace(item,"_")
          found = True
      if(not found): self.error(map_name_string)
      if(s.count('m')==1):
        if(s.count('d')==0): self.error(map_name_string)
        s = s.replace('m',"_")
        s = s.replace('d',"_")
        ml_map = True
      if(s.count('d')==1): self.error(map_name_string)
      # at this point we can only have numbers or eventually "-"
      if(s.count('-')==1 and not s.startswith('-')):
        if(len(s)==1): self.error(map_name_string)
        if(s.endswith('-')):
          n = 1
        else:
          n = s[s.index('-')+1:]
        k = s[:s.index('-')]
        s = s.replace('-',"_")
      if(len(s) == s.count('_')):
        k, n = 1, 1
        s = s.replace('_','')
      if(len(s) > 0):
        first = []
        second = []
        start_first = True
        start_second = False
        for item in s:
          if(item != '_' and start_first): first.append(item)
          if(item == '_' and len(first) > 0): start_first = False
          if(item != '_' and not start_first): second.append(item)
        tmp_result = "".join(first)
        if(tmp_result == '-'): tmp_result = -1
        else: tmp_result = "".join(first)
        try: k = float(tmp_result)
        except: self.error(map_name_string)
        if(len(second)==0): n = 1
        else:
          tmp_result = "".join(second)
          if(tmp_result.startswith('--')): tmp_result = tmp_result[1:]
          try: n = float(tmp_result)
          except: self.error(map_name_string)
      if([k,n].count(None) > 0): self.error(map_name_string)
      assert ml_map is not None
      self.k = float(k)
      self.n = float(n)
      self.ml_map = ml_map

  def error(self, s):
    raise RuntimeError("Wrong map type requested: "+s)

  def format(self):
    if(self.ml_map):
      if(self.n >= 0.):
        result = str(self.k)+"mFobs"+"-"+str(self.n)+"DFmodel"
      else:
        result = str(self.k)+"mFobs"+"-("+str(self.n)+")DFmodel"
    else:
      if(self.n >= 0.):
        result = str(self.k)+"Fobs"+"-"+str(self.n)+"Fmodel"
      else:
        result = str(self.k)+"Fobs"+"-("+str(self.n)+")Fmodel"
    return result
