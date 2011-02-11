
class fmodels(object):
  def __init__(self, fmodel_xray = None,
                     fmodel_neutron = None,
                     xray_scattering_dict = None,
                     neutron_scattering_dict = None,
                     neutron_refinement = None,
                     twin_law = None, # XXX used below in ONE plase to avoid running into a BUG in twin_f_model
                     log = None):
    self.fmodel_x = fmodel_xray
    self.fmodel_n = fmodel_neutron
    self.xray_scattering_dict = xray_scattering_dict
    self.neutron_scattering_dict = neutron_scattering_dict
    self.neutron_refinement = neutron_refinement
    self.log = log
    # pre-scale
    if(self.fmodel_n is not None and twin_law is None): # XXX This is broken if twin_f_model is used
      scale_k1_x = self.fmodel_x.scale_k1()
      scale_k1_n = self.fmodel_n.scale_k1()
      xn_scale = scale_k1_x / scale_k1_n
      f_obs_n_new = self.fmodel_n.f_obs().array(
        data = self.fmodel_n.f_obs().data()*xn_scale)
      f_obs_n_new.set_observation_type_xray_amplitude()
      self.fmodel_n.update(f_obs = f_obs_n_new)
    #
    self.create_target_functors()

  def pseudo_deep_copy(self):
    fmodel_n_dc = None
    if(self.fmodel_n is not None):
      fmodel_n_dc = self.fmodel_n.deep_copy()
    result = fmodels(fmodel_xray             = self.fmodel_x.deep_copy(),
                     fmodel_neutron          = fmodel_n_dc,
                     xray_scattering_dict    = self.xray_scattering_dict,
                     neutron_scattering_dict = self.neutron_scattering_dict,
                     neutron_refinement      = self.neutron_refinement,
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
      neutron_refinement      = self.neutron_refinement,
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

  def show_comprehensive(self, message = ""):
    from mmtbx.refinement import print_statistics
    print_statistics.make_sub_header("X-ray data", out = self.log)
    if(self.fmodel_x is not None):
      self.fmodel_xray().info().show_all(header = message, out = self.log)
    if(self.fmodel_n is not None):
      print_statistics.make_sub_header("Neutron data", out = self.log)
      self.fmodel_neutron().info().show_all(header = message, out = self.log)

  def update_bulk_solvent_and_scale(self, params = None, optimize_mask= False,
                                    optimize_mask_thorough= False,
                                    force_update_f_mask = False):
    from mmtbx.refinement import print_statistics
    print_statistics.make_header("bulk solvent modeling and scaling",
      out = self.log)
    self.update_xray_structure(update_f_calc = True, update_f_mask = True,
      force_update_f_mask = force_update_f_mask)
    if(self.fmodel_x is not None):
      if(optimize_mask_thorough):
        self.fmodel_xray().optimize_mask_and_update_solvent_and_scale(
          params = params, out = self.log, verbose =-1)
      else:
        self.fmodel_xray().update_solvent_and_scale(params = params,
          out = self.log, verbose =-1, optimize_mask = optimize_mask)
    if(self.fmodel_n is not None):
      if(optimize_mask_thorough):
        self.fmodel_neutron().optimize_mask_and_update_solvent_and_scale(
          params = params, out = self.log, verbose =-1)
      else:
        self.fmodel_neutron().update_solvent_and_scale(params = params,
          out = self.log, verbose =-1, optimize_mask = optimize_mask)
    self.show_short()

  def remove_outliers(self):
    from mmtbx.refinement import print_statistics
    n_old = self.fmodel_x.f_obs().size()
    print_statistics.make_sub_header("Outliers rejection", out = self.log)
    if(self.fmodel_x is not None):
      if(self.fmodel_n is not None):
        print_statistics.make_sub_header("x-ray data", out = self.log)
      self.fmodel_x = self.fmodel_xray().remove_outliers(
        show = True, log = self.log)
    if(self.fmodel_n is not None):
      print_statistics.make_sub_header("neutron data", out = self.log)
      self.fmodel_n = self.fmodel_neutron().remove_outliers(
        show = True, log = self.log)
    print >> self.log
    n_new = self.fmodel_x.f_obs().size()
    if(n_old != n_new):
      self.create_target_functors() # XXX Cover neutrons

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

  def prepare_target_functors_for_minimization(self):
    self.target_functor_xray.prepare_for_minimization()
    if (self.target_functor_neutron is not None):
      self.target_functor_neutron.prepare_for_minimization()

  def target_functions_are_invariant_under_allowed_origin_shifts(self):
    for f in [
          self.target_functor_xray,
          self.target_functor_neutron]:
      if (f is None): continue
      if (not f.target_function_is_invariant_under_allowed_origin_shifts()):
        return False
    return True

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

  FC = ['fcalc','fcal','fc', 'fmodel','fmod','fm']
  DFC = ['dfcalc','dfcal','dfc', 'dfmodel','dfmod','dfm']
  FO = ['fobs','fob','fo']
  MFO = ['mfobs','mfob','mfo']
  KICK = ["-kicked","+kicked","kicked","-kick","+kick","kick"]
  FILLED = ['+filled','-filled','filled','+fill','-fill','fill']

  def __init__(self, map_name_string):
    s = map_name_string.lower()
    self.k = None
    self.n = None
    self.ml_map = None
    self.anomalous = False
    self.kicked = False
    self.f_obs_filled = False
    for sym in ["~","!","@","#","$","%","^","&","*","(",")","=","<",">","?","/",
                ":",";","|","[","]","{","}",",","_"," "]:
      s = s.replace(sym,"")
    for tmp in self.KICK:
      if tmp in s:
        s = s.replace(tmp,"")
        self.kicked = True
    for tmp in self.FILLED:
      if tmp in s:
        s = s.replace(tmp,"")
        self.f_obs_filled = True
    if(s.count('ano')):
      self.anomalous = True
    elif(s in self.FC):
      self.k = 0
      self.n = 1
      self.ml_map = False
    elif(s in self.DFC):
      self.k = 0
      self.n = 1
      self.ml_map = True
    elif(s in self.FO):
      self.k = 1
      self.n = 0
      self.ml_map = False
    elif(s in self.MFO):
      self.k = 1
      self.n = 0
      self.ml_map = True
    else:
      #
      found_D = False
      for tmp in self.DFC:
        if(tmp in s):
          found_D = True
          s = s.replace(tmp,"C")
      found_M = False
      for tmp in self.MFO:
        if(tmp in s):
          found_M = True
          s = s.replace(tmp,"O")
      if(not ([found_D,found_M].count(True) in [0,2])):
        self.error(map_name_string)
      if([found_D,found_M].count(True)==2): self.ml_map=True
      elif([found_D,found_M].count(True)==0): self.ml_map=False
      else: self.error(map_name_string)
      #
      if(not self.ml_map):
        for tmp in self.FC:
          if(tmp in s): s = s.replace(tmp,"C")
        for tmp in self.FO:
          if(tmp in s): s = s.replace(tmp,"O")
      #
      if(s.count("O")+s.count("C")!=2):
        self.error(map_name_string)
      for tmp in s:
        if(not (tmp in ["C","O"])):
          if(tmp.isalpha()): self.error(map_name_string)
      if(s.index("O")<s.index("C")):
        tmp = s[s.index("O")+1:s.index("C")]
        sign = None
        if(tmp.count("+")==1): sign="+"
        elif(tmp.count("-")==1): sign="-"
        else: self.error(map_name_string)
        po = s[:s.index("O")+tmp.index(sign)+1]
        pc = s[s.index("O")+tmp.index(sign)+1:]
        po = po.replace("O","")
        pc = pc.replace("C","")
        if(len(po)==0): self.k = 1.
        elif(len(po)==1 and po in ["+","-"]): self.k = float("%s1"%po)
        else: self.k = float(po)
        if(len(pc)==0): self.n = 1.
        elif(len(pc)==1 and pc in ["+","-"]): self.n = float("%s1"%pc)
        else: self.n = float(pc)
      elif(s.index("O")>s.index("C")):
        tmp = s[s.index("C")+1:s.index("O")]
        sign = None
        if(tmp.count("+")==1): sign="+"
        elif(tmp.count("-")==1): sign="-"
        else: self.error(map_name_string)
        po = s[:s.index("C")+tmp.index(sign)+1]
        pc = s[s.index("C")+tmp.index(sign)+1:]
        po = po.replace("C","")
        pc = pc.replace("O","")
        if(len(po)==0): self.n = 1.
        elif(len(po)==1 and po in ["+","-"]): self.n = float("%s1"%po)
        else: self.n = float(po)
        if(len(pc)==0): self.k = 1.
        elif(len(pc)==1 and pc in ["+","-"]): self.k = float("%s1"%pc)
        else: self.k = float(pc)
      else: raise RuntimeError
    if(self.k is not None):
      self.k = float(self.k)
      self.n = float(self.n)
      #if self.n < 0: self.n *= -1

  def error(self, s):
    msg="""\n
Wrong map type requested: %s
  Allowed format is: %s
    where [p] and [q] are any numbers (optional),
          [m] and [D] indicate if the requested map is sigmaa (optional),
          Fo and Fc are Fobs and Fcalc,
          [kick] is for Average Kick Map (optional),
          [filled] is for missing Fobs filled map.
  Examples: 2mFo-DFc, 3.2Fo-2.3Fc, mFobs-DFcalc_kick, 2mFobs-DFcalc_filled, Fc,
            2mFobs-DFcalc_kick_fill, anom, anom_diff, anomalous_difference, Fo
"""
    format = "[p][m]Fo+[q][D]Fc[kick][filled]"
    raise RuntimeError(msg%(s,format))

  def format(self):
    if(not self.anomalous):
      if(abs(int(self.k)-self.k)<1.e-6): k = str(int(self.k))
      else: k = str(self.k)
      if(abs(int(self.n)-self.n)<1.e-6):
        sign = ""
        if(self.n>0): sign = "+"
        if(self.n==0): sign = "-"
        n = sign+str(int(self.n))
      else: n = str(self.n)
      if(k=="1" or k=="+1"): k = ""
      if(n=="1" or n=="+1"): n = "+"
      if(k=="-1"): k = "-"
      if(n=="-1"): n = "-"
      if(self.ml_map):
        result = k+"mFobs"+n+"DFmodel"
      else:
        result = k+"Fobs"+n+"Fmodel"
      if(self.kicked): result += "_kick"
      if(self.f_obs_filled): result += "_filled"
      return result
    else:
      assert [self.k,self.n,self.ml_map].count(None) == 3
      assert [self.kicked,self.f_obs_filled].count(False)==2
      return "anomalous_difference"
