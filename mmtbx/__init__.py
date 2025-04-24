"""fmodels and map_names container classes"""

from __future__ import absolute_import, division, print_function

from libtbx.utils import Sorry
from libtbx.version import get_version

__version__ = get_version()

class fmodels(object):
  """
  Container object for F_model values used during refinement.

  Attributes
  ----------
  fmodel_xray :
  fmodel_neutron :
  xray_scattering_dict :
  neutron_scattering_dict :
  neutron_refinement :
  """
  def __init__(self, fmodel_xray = None,
                     fmodel_neutron = None,
                     xray_scattering_dict = None,
                     neutron_scattering_dict = None,
                     neutron_refinement = None,
                     log = None):
    self.fmodel_x = fmodel_xray
    self.fmodel_n = fmodel_neutron
    self.xray_scattering_dict = xray_scattering_dict
    self.neutron_scattering_dict = neutron_scattering_dict
    self.neutron_refinement = neutron_refinement
    self.log = log
    self.create_target_functors()

  def pseudo_deep_copy(self):
    """
    Makes a deep copy of self.

    Returns
    -------
    mmtbx.fmodels
    """
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
    """
    Returns a copy of self with a resolution filter applied to the x-ray and
    neutron maps above a given resolution.

    Parameters
    ----------
    d_min : float
        Reflections with resolutions <= d_min are removed.

    Returns
    -------
    mmtbx.fmodels

    See Also
    --------
    mmtbx.f_model.manager.resolution_filter
    """
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
    """
    ...
    """
    if(self.fmodel_x is not None):
      if(self.fmodel_n is not None):
        if(xray_structure is not None):
          self.fmodel_x.xray_structure = xray_structure
        # mrt: discard old scattering dictionary to avoid mixing xray and neutron
        self.fmodel_x.xray_structure._scattering_type_registry = None
        # XXX: xray tables could mix here
        # update possibly changed xray dictionary
        self.xray_scattering_dict = \
          self.fmodel_x.xray_structure.scattering_type_registry(custom_dict =
          self.xray_scattering_dict).as_type_gaussian_dict()
      #assert not self.fmodel_x.xray_structure.guess_scattering_type_neutron()
    return self.fmodel_x

  def fmodel_neutron(self, xray_structure = None):
    """
    ...
    """
    if(self.fmodel_n is not None):
      if(xray_structure is not None):
        self.fmodel_n.xray_structure = xray_structure
      # mrt: discard old scattering dictionary to avoid mixing xray and neutron
      self.fmodel_n.xray_structure._scattering_type_registry = None
      # update possibly changed xray dictionary
      self.neutron_scattering_dict = \
        self.fmodel_n.xray_structure.scattering_type_registry(custom_dict =
        self.neutron_scattering_dict, table="neutron").as_type_gaussian_dict()
      assert self.fmodel_n.xray_structure.guess_scattering_type_neutron()
    return self.fmodel_n

  def update_xray_structure(self, xray_structure = None,
                                  update_f_calc  = None,
                                  update_f_mask  = None,
                                  force_update_f_mask = False):
    """
    ...
    """
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

  def show_short(self, log=None):
    """
    ...
    """
    if log is None: log = self.log
    if(self.fmodel_x is not None):
      prefix = ""
      if(self.fmodel_n is not None): prefix = "x-ray data"
      self.fmodel_xray().info().show_rfactors_targets_scales_overall(
        header = prefix, out = log)
    if(self.fmodel_n is not None):
      print(file=self.log)
      self.fmodel_neutron().info().show_rfactors_targets_scales_overall(
        header = "neutron data", out = log)

  def show_comprehensive(self, message = ""):
    """
    ...
    """
    from mmtbx.refinement import print_statistics
    print_statistics.make_sub_header("X-ray data", out = self.log)
    if(self.fmodel_x is not None):
      self.fmodel_xray().info().show_all(header = message, out = self.log)
    if(self.fmodel_n is not None):
      print_statistics.make_sub_header("Neutron data", out = self.log)
      self.fmodel_neutron().info().show_all(header = message, out = self.log)

  def update_all_scales(
        self,
        update_f_part1,
        remove_outliers     = True,
        params              = None,
        optimize_mask       = False,
        force_update_f_mask = False,
        nproc               = 1,
        log                 = None,
        apply_back_trace    = False,
        refine_hd_scattering=None):
    """
    ...
    """
    if log is None: log = self.log
    fast=True
    if(params.mode=="slow"): fast=False
    from mmtbx.refinement import print_statistics
    print_statistics.make_header("updating all scales", out = log)
    self.update_xray_structure(update_f_calc = True, update_f_mask = True,
      force_update_f_mask = force_update_f_mask)
    if([self.fmodel_x, self.fmodel_n].count(None)==0): apply_back_trace=False
    if(self.fmodel_x is not None):
      msg = None
      if(self.fmodel_n is not None):
        msg = "X-ray:"
        print(msg, file=log)
      self.fmodel_xray().update_all_scales(
        update_f_part1       = update_f_part1,
        remove_outliers      = remove_outliers,
        params               = params,
        fast                 = fast,
        log                  = log,
        show                 = True,
        optimize_mask        = optimize_mask,
        nproc                = nproc,
        apply_back_trace     = apply_back_trace,
        refine_hd_scattering = refine_hd_scattering)
      self.fmodel_x.show(log = log, suffix = msg)
    if(self.fmodel_n is not None):
      msg = "Neutron:"
      print(msg, file=log)
      self.fmodel_neutron().update_all_scales(
        update_f_part1       = update_f_part1,
        remove_outliers      = remove_outliers,
        params               = params,
        fast                 = fast,
        log                  = log,
        show                 = True,
        optimize_mask        = optimize_mask,
        nproc                = nproc,
        apply_back_trace     = apply_back_trace,
        refine_hd_scattering = refine_hd_scattering)
      self.fmodel_n.show(log = log, suffix = msg)

  def show_targets(self, log, text=""):
    """
    ...
    """
    prefix_x = ""
    if(self.fmodel_n is not None):
      prefix_x = "xray"
    self.fmodel_xray().info().show_targets(out = log, text = prefix_x+" "+text)
    if(self.fmodel_n is not None):
      self.fmodel_neutron().info().show_targets(out= log, text="neutron "+text)

  def update(self, target_name=None):
    """
    ...
    """
    if(self.fmodel_x is not None):
      self.fmodel_x.update(target_name=target_name)
    if(self.fmodel_n is not None):
      self.fmodel_n.update(target_name=target_name)

  def create_target_functors(self, alpha_beta=None):
    """
    ...
    """
    self.target_functor_xray = self.fmodel_xray().target_functor(
      alpha_beta = alpha_beta)
    self.target_functor_neutron = None
    if(self.fmodel_n is not None):
      self.target_functor_neutron = self.fmodel_neutron().target_functor()

  def prepare_target_functors_for_minimization(self):
    """
    ...
    """
    self.target_functor_xray.prepare_for_minimization()
    if (self.target_functor_neutron is not None):
      self.target_functor_neutron.prepare_for_minimization()

  def target_functions_are_invariant_under_allowed_origin_shifts(self):
    """
    ...
    """
    for f in [
          self.target_functor_xray,
          self.target_functor_neutron]:
      if (f is None): continue
      if (not f.target_function_is_invariant_under_allowed_origin_shifts()):
        return False
    return True

  def target_functor_result_xray(self, compute_gradients):
    """
    ...
    """
    fmx = self.fmodel_xray()
    return self.target_functor_xray(compute_gradients = compute_gradients)

  def target_functor_result_neutron(self, compute_gradients):
    """
    ...
    """
    result = None
    if(self.fmodel_n is not None):
      fmn = self.fmodel_neutron()
      result = self.target_functor_neutron(compute_gradients=compute_gradients)
    return result

  def target_and_gradients(self, compute_gradients, weights=None,
        u_iso_refinable_params = None, occupancy = False):
    """
    ...
    """
    tfx = self.target_functor_result_xray
    tfn = self.target_functor_result_neutron
    class tg(object):
      def __init__(self, fmodels):
        self.fmodels = fmodels
        tfx_r = tfx(compute_gradients = compute_gradients)
        wx=1.
        if(weights is not None): wx = weights.wx * weights.wx_scale
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
          wn=1
          if(weights is not None): wn = weights.wn * weights.wn_scale
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
          if(weights is not None):
            result = (self.target_work_xray * weights.wx_scale + \
                      self.target_work_neutron * weights.wn_scale) * weights.wxn
          else:
            result = self.target_work_xray + self.target_work_neutron
        else: result = self.target_work_xray_weighted
        return result
      def gradients(self):
        result = None
        if(compute_gradients):
          if(self.fmodels.fmodel_neutron() is not None):
            if(weights is not None):
              result = (self.gradient_xray  * weights.wx_scale + \
                        self.gradient_neutron * weights.wn_scale) * weights.wxn
            else:
              result = self.gradient_xray + self.gradient_neutron
          else: result = self.gradient_xray_weighted
        return result
    result = tg(fmodels = self)
    return result

class map_names(object):
  """
  Class used for parsing and external display of map's name.

  Attributes
  ----------
  k : float
     Scale for F_obs.
  n : float
      Scale for F_model.
  ml_map : bool
  anomalous : bool
  anomalous_residual : bool
  phaser_sad_llg : bool
  f_obs_filled : bool
  """

  FC = ['fcalc','fcal','fc', 'fmodel','fmod','fm']
  DFC = ['dfcalc','dfcal','dfc', 'dfmodel','dfmod','dfm']
  FO = ['fobs','fob','fo']
  MFO = ['mfobs','mfob','mfo']
  FILLED = ['+filled','-filled','filled','+fill','-fill','fill']

  def __init__(self, map_name_string):
    s = map_name_string.lower()
    self.k = None
    self.n = None
    self.ml_map = None
    self.anomalous = False
    self.anomalous_residual = False
    self.phaser_sad_llg = False
    self.f_obs_filled = False
    for sym in ["~","!","@","#","$","%","^","&","*","(",")","=","<",">","?","/",
                ":",";","|","[","]","{","}",",","_"," "]:
      s = s.replace(sym,"")
    for tmp in self.FILLED:
      if tmp in s:
        s = s.replace(tmp,"")
        self.f_obs_filled = True
    if(s.count('ano')):
      if (s.count('resid')):
        self.anomalous_residual = True
      else :
        self.anomalous = True
    elif (s.count("sad") or s.count("llg")):
      self.phaser_sad_llg = True
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
      else: raise RuntimeError("Error attempting to decode map name string "+
        "'%s'" % map_name_string)
    if(self.k is not None):
      self.k = float(self.k)
      self.n = float(self.n)
      #if self.n < 0: self.n *= -1

  def error(self, s):
    """
    Raises an exception for bad map names.
    """
    msg="""\n
Wrong map type requested: %s
  Allowed format is: %s
    where [p] and [q] are any numbers (optional),
          [m] and [D] indicate if the requested map is sigmaa (optional),
          Fo and Fc are Fobs and Fcalc,
          [filled] is for missing Fobs filled map.
  Examples: 2mFo-DFc, 3.2Fo-2.3Fc, mFobs-DFcalc, 2mFobs-DFcalc_filled, Fc,
            2mFobs-DFcalc_fill, anom, anom_diff, anomalous_difference, Fo
"""
    format = "[p][m]Fo+[q][D]Fc[filled]"
    raise Sorry(msg%(s,format))

  def format(self):
    """
    Formats a map name for external display.

    Examples
    --------
    >>> r = mmtbx.map_names(map_name_string="mFo-DFc ")
    >>> print r.format()
    mFobs-DFmodel
    """
    if (not self.anomalous) and (not self.phaser_sad_llg):
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
      if(self.f_obs_filled): result += "_filled"
      return result
    else:
      assert [self.k,self.n,self.ml_map].count(None) == 3
      assert [self.f_obs_filled].count(False)==1
      if (self.phaser_sad_llg):
        return "phaser_sad_llg"
      return "anomalous_difference"
