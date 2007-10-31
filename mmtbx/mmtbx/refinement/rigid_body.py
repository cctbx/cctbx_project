from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx import lbfgs
from mmtbx.refinement import print_statistics
import copy, time
from libtbx.utils import Sorry
from cctbx import xray
from libtbx.utils import user_plus_sys_time

time_rigid_body_total = 0.0

def euler(phi, psi, the, convention):
  if(convention == "zyz"):
     result = rb_mat_euler(the=the, psi=psi, phi=phi)
  elif(convention == "xyz"):
     result = rb_mat(the=the, psi=psi, phi=phi)
  else:
     raise Sorry("\nWrong rotation convention\n")
  return result

class rb_mat_euler(object):

   def __init__(self, the, psi, phi):
     the = the * math.pi/180
     psi = psi * math.pi/180
     phi = phi * math.pi/180
     self.c_psi = math.cos(psi)
     self.c_the = math.cos(the)
     self.c_phi = math.cos(phi)
     self.s_psi = math.sin(psi)
     self.s_the = math.sin(the)
     self.s_phi = math.sin(phi)

   def rot_mat(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 =  c_the*c_psi*c_phi - s_the*s_phi
     r12 = -c_the*c_psi*s_phi - s_the*c_phi
     r13 =  c_the*s_psi
     r21 =  s_the*c_psi*c_phi + c_the*s_phi
     r22 = -s_the*c_psi*s_phi + c_the*c_phi
     r23 =  s_the*s_psi
     r31 = -s_psi*c_phi
     r32 =  s_psi*s_phi
     r33 =  c_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_the(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 = -s_the*c_psi*c_phi - c_the*s_phi
     r12 =  s_the*c_psi*s_phi - c_the*c_phi
     r13 = -s_the*s_psi
     r21 =  c_the*c_psi*c_phi - s_the*s_phi
     r22 = -c_the*c_psi*s_phi - s_the*c_phi
     r23 =  c_the*s_psi
     r31 = 0.0
     r32 = 0.0
     r33 = 0.0
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_psi(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 = -c_the*s_psi*c_phi
     r12 =  c_the*s_psi*s_phi
     r13 =  c_the*c_psi
     r21 = -s_the*s_psi*c_phi
     r22 =  s_the*s_psi*s_phi
     r23 =  s_the*c_psi
     r31 = -c_psi*c_phi
     r32 =  c_psi*s_phi
     r33 = -s_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_phi(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 = -c_the*c_psi*s_phi - s_the*c_phi
     r12 = -c_the*c_psi*c_phi + s_the*s_phi
     r13 =  0.0
     r21 = -s_the*c_psi*s_phi + c_the*c_phi
     r22 = -s_the*c_psi*c_phi - c_the*s_phi
     r23 =  0.0
     r31 =  s_psi*s_phi
     r32 =  s_psi*c_phi
     r33 =  0.0
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

class rigid_body_shift_accumulator(object):

   def __init__(self, euler_angle_convention):
     self.euler_angle_convention = euler_angle_convention
     self.rotations = []
     self.translations = []

   def add(self, rotations, translations):
     assert len(rotations) == len(translations)
     new_rotations = []
     new_translations = []
     if(len(self.rotations) > 0):
        for rn, tn, r, t in zip(rotations, translations, self.rotations,
                                                            self.translations):
            new_rotations.append(rn + r)
            new_translations.append(tn + t)
     else:
        for rn, tn in zip(rotations, translations):
            new_rotations.append(rn)
            new_translations.append(tn)
     self.rotations = new_rotations
     self.translations = new_translations

   def show(self, out = None):
     if (out is None): out = sys.stdout
     print >> out, "|-rigid body shift (total)------------------------------"\
                   "----------------------|"
     print_statistics.show_rigid_body_rotations_and_translations(
       out=out,
       prefix="",
       frame="|",
       euler_angle_convention=self.euler_angle_convention,
       rotations=self.rotations,
       translations=self.translations)
     print >> out, "|"+"-"*77+"|"
     print >> out

class rb_mat(object):

   def __init__(self, phi, psi, the):
     phi = phi * math.pi/180
     psi = psi * math.pi/180
     the = the * math.pi/180
     self.c_psi = math.cos(psi)
     self.c_phi = math.cos(phi)
     self.c_the = math.cos(the)
     self.s_psi = math.sin(psi)
     self.s_phi = math.sin(phi)
     self.s_the = math.sin(the)

   def rot_mat(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 =  c_psi*c_phi
     r12 = -c_psi*s_phi
     r13 =  s_psi
     r21 =  c_the*s_phi + s_the*s_psi*c_phi
     r22 =  c_the*c_phi - s_the*s_psi*s_phi
     r23 = -s_the*c_psi
     r31 =  s_the*s_phi - c_the*s_psi*c_phi
     r32 =  s_the*c_phi + c_the*s_psi*s_phi
     r33 =  c_the*c_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_phi(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 = -c_psi*s_phi
     r12 = -c_psi*c_phi
     r13 =  0.0
     r21 =  c_the*c_phi - s_the*s_psi*s_phi
     r22 = -c_the*s_phi - s_the*s_psi*c_phi
     r23 =  0.0
     r31 =  s_the*c_phi + c_the*s_psi*s_phi
     r32 = -s_the*s_phi + c_the*s_psi*c_phi
     r33 =  0.0
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_psi(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 = -s_psi*c_phi
     r12 =  s_psi*s_phi
     r13 =  c_psi
     r21 =  s_the*c_psi*c_phi
     r22 = -s_the*c_psi*s_phi
     r23 =  s_the*s_psi
     r31 = -c_the*c_psi*c_phi
     r32 =  c_the*c_psi*s_phi
     r33 = -c_the*s_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_the(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 =  0.0
     r12 =  0.0
     r13 =  0.0
     r21 = -s_the*s_phi+c_the*s_psi*c_phi
     r22 = -s_the*c_phi-c_the*s_psi*s_phi
     r23 = -c_the*c_psi
     r31 =  c_the*s_phi+s_the*s_psi*c_phi
     r32 =  c_the*c_phi-s_the*s_psi*s_phi
     r33 = -s_the*c_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

def split_resolution_range(miller_array,
                           nref_first_low,
                           d_low = 8.0,
                           d_high = 2.0,
                           number_of_zones = 0,
                           log = None):
  if(log is None): log = sys.stdout
  assert number_of_zones > 0
  assert d_low > d_high
  assert nref_first_low > 0
  d_spacings = miller_array.d_spacings().data()
  d_max, d_min = flex.max(d_spacings), flex.min(d_spacings)
  d_high = max(d_min, d_high)
  sel = flex.sort_permutation(d_spacings, reverse = True)
  d_spacings = d_spacings.select(sel)
  d_nref_first_low = d_spacings[:nref_first_low]
  if((d_nref_first_low < d_high).count(True) > 0):
    sel_d_high = d_spacings >= d_high
    d_nref_first_low = d_spacings.select(sel_d_high)
    try: d_nref_first_low = d_nref_first_low[:nref_first_low]
    except: pass
  d_low_first = d_nref_first_low[len(d_nref_first_low)-1]
  if(d_low_first > d_low): d_low_first = d_low
  d_mids = []
  if(number_of_zones != 1):
     d_step = (d_low_first - d_high)/(number_of_zones - 1)
     if(d_step > 0.1):
       for i_zone in range(0, number_of_zones):
         d_mids.append(d_low_first - d_step*i_zone)
     else:
       d_mids.append(d_high)
  else:
     d_mids.append(d_high)
  print >> log, \
    "Information about resolution zones and data for rigid body refinement:"
  print >> log, \
    "  Data resolution: %.2f - %.2f"%(d_max,d_min)
  print >> log, \
    "  Resolution for rigid body refinement: %.2f - %.2f"%(d_max,d_high)
  print >> log, \
    "  Number of zones for Multi-Zone (MZ) protocol: %d"%(number_of_zones)
  print >> log, \
    "  Requested number of reflections for first low resolution zone: %d"%(
                                                                nref_first_low)
  print >> log, "  High resolution cutoffs for MZ protocol: "
  print >> log, "    zone  resolution  number of reflections"
  for i, d_i in enumerate(d_mids):
    nref = (d_spacings >= d_i).count(True)
    print >> log, \
      "     %d:   %.2f-%.2f                  %d" % (i+1, d_max, d_i, nref)
  return d_mids

class manager(object):
  def __init__(self, fmodel,
                     selections              = None,
                     refine_r                = True,
                     refine_t                = True,
                     r_initial               = None,
                     t_initial               = None,
                     nref_min                = 1000,
                     max_iterations          = 50,
                     bulk_solvent_and_scale  = True,
                     high_resolution         = 2.0,
                     max_low_high_res_limit  = 8.0,
                     bss                     = None,
                     euler_angle_convention  = "xyz",
                     lbfgs_maxfev            = 10,
                     number_of_zones         = 5,
                     protocol                = None,
                     log                     = None,
                     monitors = None):
    global time_rigid_body_total
    timer_rigid_body_total = user_plus_sys_time()
    save_r_work = fmodel.r_work()
    save_r_free = fmodel.r_free()
    save_xray_structure = fmodel.xray_structure.deep_copy_scatterers()
    xray.set_scatterer_grad_flags(
                               scatterers = fmodel.xray_structure.scatterers(),
                               site       = True)
    self.euler_angle_convention = euler_angle_convention
    if(protocol not in ["one_zone","multiple_zones"]):
       raise Sorry("Wrong rigid body refinement protocol: %s"%str(protocol))
    if(log is None): log = sys.stdout
    if(selections is None):
       selections = []
       selections.append(flex.bool(fmodel.xray_structure.scatterers().size(),
                                                                         True))
    else: assert len(selections) > 0
    self.total_rotation = []
    self.total_translation = []
    for item in selections:
        self.total_rotation.append(flex.double(3,0))
        self.total_translation.append(flex.double(3,0))
    if(r_initial is None):
       r_initial = []
       for item in selections:
           r_initial.append(flex.double(3,0))
    if(t_initial is None):
       t_initial = []
       for item in selections:
           t_initial.append(flex.double(3,0))
    fmodel_copy = fmodel.deep_copy()
    if(fmodel_copy.mask_params is not None):
       fmodel_copy.mask_params.verbose = -1
    if(protocol == "one_zone"):
      number_of_zones = 1
    d_mins = split_resolution_range(miller_array    = fmodel_copy.f_obs_w,
                                    nref_first_low  = nref_min,
                                    d_low           = max_low_high_res_limit,
                                    d_high          = high_resolution,
                                    number_of_zones = number_of_zones,
                                    log             = log)
    print >> log
    self.show(fmodel = fmodel,
              r_mat  = self.total_rotation,
              t_vec  = self.total_translation,
              header = "Start",
              out    = log)
    if (protocol == "one_zone" or monitors is None):
      monitors_call_back_after_collect = None
    else:
      monitors_call_back_after_collect = monitors.call_back_after_collect
      if (monitors_call_back_after_collect is not None):
        monitors_call_back_after_collect(
          monitor=None, model=None, fmodel=fmodel)
    for res in d_mins:
        xrs = fmodel_copy.xray_structure.deep_copy_scatterers()
        fmodel_copy = fmodel.resolution_filter(d_min = res)
        d_max_min = fmodel_copy.f_obs_w.d_max_min()
        line = "Refinement at resolution: "+\
                 str("%7.2f"%d_max_min[0]).strip()+" - "+\
                 str("%6.2f"%d_max_min[1]).strip()
        print_statistics.make_sub_header(line, out = log)
        fmodel_copy.update_xray_structure(xray_structure = xrs,
                                          update_f_calc  = True)
        rworks = flex.double()
        if(len(d_mins) == 1):
           n_rigid_body_minimizer_cycles = 1
        else:
           n_rigid_body_minimizer_cycles = min(int(res),4)
        for i_macro_cycle in xrange(n_rigid_body_minimizer_cycles):
            if(bss is not None and bulk_solvent_and_scale):
               if(fmodel_copy.f_obs.d_min() > 3.0):
                  save_bss_anisotropic_scaling = bss.anisotropic_scaling
                  bss.anisotropic_scaling=False
               fmodel_copy.update_solvent_and_scale(params  = bss,
                                                    out     = log,
                                                    verbose = -1)
               if(fmodel_copy.f_obs.d_min() > 3.0):
                  bss.anisotropic_scaling=save_bss_anisotropic_scaling
            minimized = rigid_body_minimizer(
                          fmodel                 = fmodel_copy,
                          selections             = selections,
                          r_initial              = r_initial,
                          t_initial              = t_initial,
                          refine_r               = refine_r,
                          refine_t               = refine_t,
                          max_iterations         = max_iterations,
                          euler_angle_convention = self.euler_angle_convention,
                          lbfgs_maxfev           = lbfgs_maxfev)
            rotation_matrices = []
            translation_vectors = []
            for i in xrange(len(selections)):
                self.total_rotation[i] += flex.double(minimized.r_min[i])
                self.total_translation[i] += flex.double(minimized.t_min[i])
                rot_obj = euler(phi        = minimized.r_min[i][0],
                                psi        = minimized.r_min[i][1],
                                the        = minimized.r_min[i][2],
                                convention = self.euler_angle_convention)
                rotation_matrices.append(rot_obj.rot_mat())
                translation_vectors.append(minimized.t_min[i])
            new_xrs = apply_transformation(
                         xray_structure      = minimized.fmodel.xray_structure,
                         rotation_matrices   = rotation_matrices,
                         translation_vectors = translation_vectors,
                         selections          = selections)
            fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                              update_f_calc  = True,
                                              update_f_mask  = True,
                                              out            = log)
            rwork = minimized.fmodel.r_work()
            rfree = minimized.fmodel.r_free()
            assert approx_equal(rwork, fmodel_copy.r_work())
        fmodel.update_xray_structure(xray_structure = fmodel_copy.xray_structure,
                                     update_f_calc  = True,
                                     update_f_mask  = True,
                                     out            = log)
        if(bss is not None and bulk_solvent_and_scale):
          fmodel.update_solvent_and_scale(params = bss, out = log, verbose= -1)
        self.show(fmodel = fmodel,
              r_mat  = self.total_rotation,
              t_vec  = self.total_translation,
              header = "Rigid body refinement",
              out    = log)
        if (monitors_call_back_after_collect is not None):
          monitors_call_back_after_collect(
            monitor=None, model=None, fmodel=fmodel)
    if(bss is not None and bulk_solvent_and_scale):
      fmodel.update_solvent_and_scale(out = log, verbose = -1)
    print >> log
    self.show(fmodel = fmodel,
              r_mat  = self.total_rotation,
              t_vec  = self.total_translation,
              header = "Rigid body end",
              out    = log)
    print >> log
    self.evaluate_after_end(fmodel, save_r_work, save_r_free,
                                                      save_xray_structure, log)
    self.fmodel = fmodel
    time_rigid_body_total += timer_rigid_body_total.elapsed()

  def evaluate_after_end(self, fmodel, save_r_work, save_r_free,
                                                     save_xray_structure, log):
    r_work = fmodel.r_work()
    r_free = fmodel.r_free()
    if((r_work > save_r_work and abs(r_work-save_r_work) > 0.005) or
       (r_free > save_r_free and abs(r_free-save_r_free) > 0.005)):
       print >> log
       print >> log, "The model after this rigid-body refinement step is not accepted."
       print >> log, "Reason: increase in R-factors after refinement."
       print >> log, "Start/final R-work: %6.4f/%-6.4f"%(save_r_work, r_work)
       print >> log, "Start/final R-free: %6.4f/%-6.4f"%(save_r_free, r_free)
       print >> log, "Return back to the previous model."
       print >> log
       fmodel.update_xray_structure(xray_structure = save_xray_structure,
                                    update_f_calc  = True,
                                    update_f_mask  = True)
       fmodel.info().show_rfactors_targets_scales_overall(
         header = "rigid body after step back", out = log)
       print >> log

  def rotation(self):
    return self.total_rotation

  def translation(self):
    return self.total_translation

  def show(self, fmodel, r_mat, t_vec, header = "", out=None):
    if(out is None): out = sys.stdout
    fmodel_info = fmodel.info()
    fmodel_info._header_resolutions_nreflections(header=header, out=out)
    print >> out, "| "+"  "*38+"|"
    fmodel_info._rfactors_and_bulk_solvent_and_scale_params(out=out)
    print >> out, "| "+"  "*38+"|"
    print >> out,"| Rigid body shift (Euler angles %s):"% \
                                        self.euler_angle_convention+" "*40 +"|"
    print_statistics.show_rigid_body_rotations_and_translations(
      out=out,
      prefix="",
      frame="|",
      euler_angle_convention=self.euler_angle_convention,
      rotations=r_mat,
      translations=t_vec)
    print >> out, "|" +"-"*77+"|"

class rigid_body_minimizer(object):
  def __init__(self,
               fmodel,
               selections,
               r_initial,
               t_initial,
               refine_r,
               refine_t,
               max_iterations,
               euler_angle_convention = "xyz",
               lbfgs_maxfev = 10):
    adopt_init_args(self, locals())
    self.fmodel_copy = self.fmodel.deep_copy()
    self.target_functor = self.fmodel_copy.target_functor()
    self.atomic_weights = self.fmodel.xray_structure.atomic_weights()
    self.sites_cart = self.fmodel.xray_structure.sites_cart()
    self.sites_frac = self.fmodel.xray_structure.sites_frac()
    self.n_groups = len(self.selections)
    assert self.n_groups > 0
    self.counter=0
    assert len(self.r_initial)  == len(self.t_initial)
    assert len(self.selections) == len(self.t_initial)
    self.dim_r = 3
    self.dim_t = 3
    self.r_min = copy.deepcopy(self.r_initial)
    self.t_min = copy.deepcopy(self.t_initial)
    for i in xrange(len(self.r_min)):
        self.r_min[i] = tuple(self.r_min[i])
        self.t_min[i] = tuple(self.t_min[i])
    self.x = self.pack(self.r_min, self.t_min)
    self.n = self.x.size()
    self.minimizer = lbfgs.run(
               target_evaluator = self,
               core_params = lbfgs.core_parameters(
                    maxfev = lbfgs_maxfev),
               termination_params = lbfgs.termination_parameters(
                    max_iterations = max_iterations),
               exception_handling_params = lbfgs.exception_handling_parameters(
                    ignore_line_search_failed_step_at_lower_bound = True,
                    ignore_line_search_failed_step_at_upper_bound = True)
                              )
    self.compute_functional_and_gradients(suppress_gradients=True)
    del self.x

  def pack(self, r, t):
    v = []
    for ri,ti in zip(r,t):
        if(self.refine_r): v += list(ri)
        if(self.refine_t): v += list(ti)
    return flex.double(tuple(v))

  def unpack_x(self):
    i = 0
    for j in xrange(self.n_groups):
        if(self.refine_r):
           self.r_min[j] = tuple(self.x)[i:i+self.dim_r]
           i += self.dim_r
        if(self.refine_t):
           self.t_min[j] = tuple(self.x)[i:i+self.dim_t]
           i += self.dim_t

  def compute_functional_and_gradients(self, suppress_gradients=False):
    self.unpack_x()
    self.counter += 1
    rotation_matrices   = []
    translation_vectors = []
    rot_objs = []
    for i in xrange(self.n_groups):
        rot_obj = euler(phi        = self.r_min[i][0],
                        psi        = self.r_min[i][1],
                        the        = self.r_min[i][2],
                        convention = self.euler_angle_convention)
        rotation_matrices.append(rot_obj.rot_mat())
        translation_vectors.append(self.t_min[i])
        rot_objs.append(rot_obj)
    new_sites_frac, new_sites_cart, centers_of_mass = apply_transformation_(
                              xray_structure      = self.fmodel.xray_structure,
                              sites_cart          = self.sites_cart,
                              sites_frac          = self.sites_frac,
                              rotation_matrices   = rotation_matrices,
                              translation_vectors = translation_vectors,
                              selections          = self.selections,
                              atomic_weights      = self.atomic_weights)
    self.fmodel_copy.xray_structure.set_sites_frac(new_sites_frac)
    new_xrs = self.fmodel_copy.xray_structure
    self.fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                           update_f_calc  = True)
    tg_obj = target_and_grads(
                   centers_of_mass = centers_of_mass,
                   sites_cart      = new_sites_cart,
                   target_functor  = self.target_functor,
                   rot_objs        = rot_objs,
                   selections      = self.selections,
                   suppress_gradients = suppress_gradients)
    self.f = tg_obj.target()
    if (suppress_gradients):
      self.g = None
    else:
      self.g = self.pack( tg_obj.gradients_wrt_r(), tg_obj.gradients_wrt_t() )
    return self.f, self.g

def apply_transformation_(xray_structure,
                          sites_cart,
                          sites_frac,
                          rotation_matrices,
                          translation_vectors,
                          selections,
                          atomic_weights):
  assert len(selections) == len(rotation_matrices)
  assert len(selections) == len(translation_vectors)
  centers_of_mass = []
  sites_cart = sites_cart.deep_copy()
  sites_frac = sites_frac.deep_copy()
  for sel,rot,trans in zip(selections,rotation_matrices,translation_vectors):
      apply_rigid_body_shift_obj = xray_structure.apply_rigid_body_shift_obj(
                                   sites_cart     = sites_cart,
                                   sites_frac     = sites_frac,
                                   rot            = rot.as_mat3(),
                                   trans          = trans,
                                   atomic_weights = atomic_weights,
                                   unit_cell      = xray_structure.unit_cell(),
                                   selection      = sel)
      sites_cart = apply_rigid_body_shift_obj.sites_cart
      sites_frac = apply_rigid_body_shift_obj.sites_frac
      centers_of_mass.append(apply_rigid_body_shift_obj.center_of_mass)
  return sites_frac, sites_cart, centers_of_mass

def apply_transformation(xray_structure,
                         rotation_matrices,
                         translation_vectors,
                         selections):
   assert len(selections) == len(rotation_matrices)
   assert len(selections) == len(translation_vectors)
   new_sites = xray_structure.sites_cart()
   for sel,rot,trans in zip(selections,rotation_matrices,translation_vectors):
       xrs = xray_structure.select(sel)
       cm_cart = xrs.center_of_mass()
       sites_cart = xrs.sites_cart()
       sites_cart_cm = sites_cart - cm_cart
       tmp = list(rot) * sites_cart_cm + trans + cm_cart
       new_sites.set_selected(sel, tmp)
   new_xrs = xray_structure.replace_sites_cart(new_sites = new_sites)
   return new_xrs

class target_and_grads(object):
  def __init__(self, centers_of_mass,
                     sites_cart,
                     target_functor,
                     rot_objs,
                     selections,
                     suppress_gradients):
    t_r = target_functor(compute_gradients=not suppress_gradients)
    self.f = t_r.target_work()
    if (suppress_gradients):
      self.grads_wrt_r = None
      self.grads_wrt_t = None
      return
    target_grads_wrt_xyz = t_r.gradients_wrt_atomic_parameters(site=True)
    self.grads_wrt_r = []
    self.grads_wrt_t = []
    target_grads_wrt_xyz = flex.vec3_double(target_grads_wrt_xyz.packed())
    for sel,rot_obj, cm in zip(selections, rot_objs, centers_of_mass):
        sites_cart_cm = sites_cart.select(sel) - cm
        target_grads_wrt_xyz_sel = target_grads_wrt_xyz.select(sel)
        target_grads_wrt_r = matrix.sqr(
                    sites_cart_cm.transpose_multiply(target_grads_wrt_xyz_sel))
        self.grads_wrt_t.append(flex.double(target_grads_wrt_xyz_sel.sum()))
        g_phi = (rot_obj.r_phi() * target_grads_wrt_r).trace()
        g_psi = (rot_obj.r_psi() * target_grads_wrt_r).trace()
        g_the = (rot_obj.r_the() * target_grads_wrt_r).trace()
        self.grads_wrt_r.append(flex.double([g_phi, g_psi, g_the]))

  def target(self):
    return self.f

  def gradients_wrt_r(self):
    return self.grads_wrt_r

  def gradients_wrt_t(self):
    return self.grads_wrt_t
