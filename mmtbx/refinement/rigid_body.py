from __future__ import division
from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx import lbfgs
from mmtbx.refinement import print_statistics
import iotbx.phil
import copy
from libtbx.utils import Sorry
from libtbx.utils import user_plus_sys_time
from libtbx.math_utils import iround

time_rigid_body_total = 0.0

def euler(phi, psi, the, convention):
  if(convention == "zyz"):
     result = rb_mat_zyz(the=the, psi=psi, phi=phi)
  elif(convention == "xyz"):
     result = rb_mat_xyz(the=the, psi=psi, phi=phi)
  else:
     raise Sorry("\nWrong rotation convention\n")
  return result

class rb_mat_zyz(object):

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

class rb_mat_xyz(object):

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

multiple_zones_params_str = """\
  min_number_of_reflections = 200
    .type = int
    .help = Number of reflections that defines the first lowest resolution \
            zone for the multiple_zones protocol. If very large \
            displacements are expected, decreasing this parameter to 100 \
            may lead to a larger convergence radius.
  multi_body_factor = 1
    .type = float
  zone_exponent = 3.0
    .type = float
  high_resolution = 3.0
    .type = float
    .help = High resolution cutoff (used for rigid body refinement only)
  max_low_high_res_limit = None
    .type = float
    .expert_level=2
    .help = Maximum value for high resolution cutoff for the first lowest \
            resolution zone
  number_of_zones = 5
    .type = int
    .help = Number of resolution zones for MZ protocol
"""

master_params = iotbx.phil.parse("""\
    mode = *first_macro_cycle_only every_macro_cycle
      .type = choice
      .help = Defines how many times the rigid body refinement is performed \
              during refinement run. first_macro_cycle_only to run only once \
              at first macrocycle, every_macro_cycle to do rigid body \
              refinement main.number_of_macro_cycles times
    target = ls_wunit_k1 ml *auto
      .type = choice
      .help = Rigid body refinement target function: least-squares or \
              maximum-likelihood
    target_auto_switch_resolution = 6.0
      .type = float
      .help = Used if target=auto, use optimal target for given working \
              resolution.
    refine_rotation = True
      .type = bool
      .help = Only rotation is refined (translation is fixed).
    refine_translation = True
      .type = bool
      .help = Only translation is refined (rotation is fixed).
    max_iterations = 25
      .type = int
      .help = Number of LBFGS minimization iterations
    bulk_solvent_and_scale = True
      .type = bool
      .help = Bulk-solvent and scaling within rigid body refinement (needed \
              since large rigid body shifts invalidate the mask).
    euler_angle_convention = *xyz zyz
      .type = choice
      .expert_level=2
      .help = Euler angles convention
    lbfgs_line_search_max_function_evaluations = 10
      .type = int
      .expert_level=2
    %s
""" % multiple_zones_params_str)

def split_resolution_range(
      d_spacings,
      n_bodies,
      target,
      target_auto_switch_resolution,
      n_ref_first,
      multi_body_factor_n_ref_first,
      d_low,
      d_high,
      number_of_zones,
      zone_exponent,
      log = None):
  assert n_bodies > 0
  assert target_auto_switch_resolution > 0
  assert multi_body_factor_n_ref_first is None or multi_body_factor_n_ref_first > 0
  assert n_ref_first is None or n_ref_first > 0
  assert d_low is None or d_low > 0
  assert d_high is None or d_high > 0
  assert n_ref_first is not None or d_low is not None
  if (d_low is not None and d_high is not None): assert d_low > d_high
  assert number_of_zones is None or number_of_zones > 0
  assert zone_exponent > 0
  if (log is None): log = sys.stdout
  n_refl_data = d_spacings.size()
  d_spacings = d_spacings.select(
    flex.sort_permutation(d_spacings, reverse = True))
  d_max, d_min = d_spacings[0], d_spacings[-1]
  if (d_high is not None and d_min < d_high):
    d_spacings = d_spacings.select(d_spacings >= d_high)
  d_high = d_spacings[-1]
  m_ref_first = n_ref_first
  if (n_ref_first is None):
    final_n_ref_first = (d_spacings >= d_low).count(True)
  else:
    if (multi_body_factor_n_ref_first is not None):
      m_ref_first += iround(
        m_ref_first * (n_bodies-1) * multi_body_factor_n_ref_first)
    assert m_ref_first > 0
    if (    d_low is not None
        and m_ref_first <= d_spacings.size()
        and d_spacings[m_ref_first-1] > d_low):
      final_n_ref_first = (d_spacings >= d_low).count(True)
    else:
      final_n_ref_first = m_ref_first
  d_mins = []
  if (number_of_zones is not None and number_of_zones > 1):
    degenerate = final_n_ref_first > d_spacings.size()
    if (degenerate):
      d_mins.append(d_high)
    else:
      zone_factor = (d_spacings.size() - final_n_ref_first) \
                  / ((number_of_zones-1)**zone_exponent)
      d_mins.append(d_spacings[final_n_ref_first-1])
      for i_zone in range(1, number_of_zones):
        n = iround(final_n_ref_first + zone_factor * i_zone**zone_exponent)
        if (i_zone == number_of_zones - 1):
          assert n == d_spacings.size() # sanity check
        d_mins.append(d_spacings[n-1])
  else:
    final_n_ref_first = min(final_n_ref_first, d_spacings.size())
    degenerate = False
    d_mins.append(d_high)
  print >> log, "Rigid body refinement:"
  print >> log, "  Requested number of resolution zones: %d" % number_of_zones
  if (len(d_mins) != 1 or degenerate):
    print >> log, "  Calculation for first resolution zone:"
    print >> log, "    Requested number of reflections per body:", n_ref_first
    print >> log, "    Requested factor per body:", \
      multi_body_factor_n_ref_first
    print >> log, "    Number of bodies:", n_bodies
    print >> log, "    Resulting number of reflections:", m_ref_first
    print >> log, "    Requested low-resolution limit:", d_low,
    if (final_n_ref_first != m_ref_first):
      print >> log, "(determines final number)",
    print >> log
    print >> log, "    Final number of reflections:", final_n_ref_first
  print >> log, "  Data resolution:                      %6.2f - %6.2f" \
    " (%d reflections)" % (d_max, d_min, n_refl_data)
  print >> log, "  Resolution for rigid body refinement: %6.2f - %6.2f" \
    " (%d reflections)" % (d_max, d_high, d_spacings.size())
  if (degenerate):
    print >> log, """\
  WARNING: Final number of reflections for first resolution zone is greater
           than the number of available reflections (%d > %d).
  INFO: Number of resolution zones reset to 1.""" % (
    final_n_ref_first, d_spacings.size())
  target_names = []
  for d_min in d_mins:
    if (target == "auto"):
      if (d_min > target_auto_switch_resolution):
        target_names.append("ls_wunit_k1")
      else:
        target_names.append("ml")
    else:
      target_names.append(target)
  if (len(d_mins) > 1):
    print >> log, "  Resolution cutoffs for multiple zones: "
    print >> log, "                          number of"
    print >> log, "    zone     resolution  reflections  target"
    for i, d_i in enumerate(d_mins):
      n_ref = (d_spacings >= d_i).count(True)
      print >> log, "    %3d  %6.2f -%6.2f %11d    %s" % (
        i+1, d_max, d_i, n_ref, target_names[i])
    print >> log, "    zone number of reflections =" \
      " %d + %.6g * (zone-1)**%.6g" % (
        final_n_ref_first, zone_factor, zone_exponent)
  return d_mins, target_names

class manager(object):
  def __init__(self, fmodel,
                     selections = None,
                     params     = None,
                     r_initial  = None,
                     t_initial  = None,
                     bss        = None,
                     log        = None,
                     monitors   = None):
    global time_rigid_body_total
    self.params = params
    save_original_target_name = fmodel.target_name
    save_bss_anisotropic_scaling = None
    if(bss is not None):
      save_bss_anisotropic_scaling = bss.anisotropic_scaling
    timer_rigid_body_total = user_plus_sys_time()
    save_r_work = fmodel.r_work()
    save_r_free = fmodel.r_free()
    save_xray_structure = fmodel.xray_structure.deep_copy_scatterers()
    if(log is None): log = sys.stdout
    if(selections is None):
      selections = []
      selections.append(flex.bool(
        fmodel.xray_structure.scatterers().size(), True).iselection())
    else: assert len(selections) > 0
    fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    fmodel.xray_structure.scatterers().flags_set_grad_site(
      iselection = flex.bool(fmodel.xray_structure.scatterers().size(), True
      ).iselection())
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
    d_mins, target_names = split_resolution_range(
      d_spacings                    = fmodel_copy.f_obs_work().d_spacings().data(),
      n_bodies                      = len(selections),
      target                        = params.target,
      target_auto_switch_resolution = params.target_auto_switch_resolution,
      n_ref_first                   = params.min_number_of_reflections,
      multi_body_factor_n_ref_first = params.multi_body_factor,
      d_low                         = params.max_low_high_res_limit,
      d_high                        = params.high_resolution,
      number_of_zones               = params.number_of_zones,
      zone_exponent                 = params.zone_exponent,
      log                           = log)
    print >> log
    if (fmodel.target_name != target_names[0]):
      fmodel.update(target_name=target_names[0])
    self.show(fmodel = fmodel,
              r_mat  = self.total_rotation,
              t_vec  = self.total_translation,
              header = "Start",
              out    = log)
    if (params.number_of_zones == 1 or monitors is None):
      monitors_call_back_handler = None
    else:
      monitors_call_back_handler = monitors.call_back_handler
      if (monitors_call_back_handler is not None):
        monitors_call_back_handler(
          monitor=None, model=None, fmodel=fmodel, method="rigid_body")
    for res,target_name in zip(d_mins, target_names):
        xrs = fmodel_copy.xray_structure.deep_copy_scatterers()
        fmodel_copy = fmodel.resolution_filter(d_min = res)
        if (fmodel_copy.target_name != target_name):
          fmodel_copy.update(target_name=target_name)
        d_max_min = fmodel_copy.f_obs_work().d_max_min()
        line = "Refinement at resolution: "+\
                 str("%7.2f"%d_max_min[0]).strip() + " - " \
               + str("%6.2f"%d_max_min[1]).strip() \
               + " target=" + fmodel_copy.target_name
        print_statistics.make_sub_header(line, out = log)
        fmodel_copy.update_xray_structure(xray_structure = xrs,
                                          update_f_calc  = True)
        rworks = flex.double()
        if(len(d_mins) == 1):
           n_rigid_body_minimizer_cycles = 1
        else:
           n_rigid_body_minimizer_cycles = min(int(res),4)
        for i_macro_cycle in xrange(n_rigid_body_minimizer_cycles):
            if(bss is not None and params.bulk_solvent_and_scale):
               if(fmodel_copy.f_obs().d_min() > 3.0):
                  bss.anisotropic_scaling=False
                  bss.apply_back_trace_of_b_cart=False
               fmodel_copy.update_solvent_and_scale(
                 params = bss, out = log, verbose = -1, optimize_mask=False)
               if(fmodel_copy.f_obs().d_min() > 3.0):
                  assert save_bss_anisotropic_scaling is not None
                  bss.anisotropic_scaling = save_bss_anisotropic_scaling
                  bss.minimization_b_cart = save_bss_anisotropic_scaling
            minimized = rigid_body_minimizer(
              fmodel                 = fmodel_copy,
              selections             = selections,
              r_initial              = r_initial,
              t_initial              = t_initial,
              refine_r               = params.refine_rotation,
              refine_t               = params.refine_translation,
              max_iterations         = params.max_iterations,
              euler_angle_convention = params.euler_angle_convention,
              lbfgs_maxfev = params.lbfgs_line_search_max_function_evaluations)
            rotation_matrices = []
            translation_vectors = []
            for i in xrange(len(selections)):
                self.total_rotation[i] += flex.double(minimized.r_min[i])
                self.total_translation[i] += flex.double(minimized.t_min[i])
                rot_obj = euler(phi        = minimized.r_min[i][0],
                                psi        = minimized.r_min[i][1],
                                the        = minimized.r_min[i][2],
                                convention = params.euler_angle_convention)
                rotation_matrices.append(rot_obj.rot_mat())
                translation_vectors.append(minimized.t_min[i])
            new_xrs = apply_transformation(
                         xray_structure      = minimized.fmodel.xray_structure,
                         rotation_matrices   = rotation_matrices,
                         translation_vectors = translation_vectors,
                         selections          = selections)
            fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                              update_f_calc  = True,
                                              update_f_mask  = True)
            rwork = minimized.fmodel.r_work()
            rfree = minimized.fmodel.r_free()
            assert approx_equal(rwork, fmodel_copy.r_work())
        fmodel.update_xray_structure(
          xray_structure = fmodel_copy.xray_structure,
          update_f_calc  = True,
          update_f_mask  = True)
        if(bss is not None and params.bulk_solvent_and_scale):
          fmodel.update_solvent_and_scale(params = bss, out = log, verbose= -1,
            optimize_mask=False)
        self.show(fmodel = fmodel,
                  r_mat  = self.total_rotation,
                  t_vec  = self.total_translation,
                  header = "Rigid body refinement",
                  out    = log)
        if (monitors_call_back_handler is not None):
          monitors_call_back_handler(
            monitor=None, model=None, fmodel=fmodel, method="rigid_body")
    if(bss is not None and params.bulk_solvent_and_scale):
      fmodel.update_solvent_and_scale(params = bss, out=log, verbose=-1,
        optimize_mask=False)
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
    self.fmodel.update(target_name = save_original_target_name)
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
      self.params.euler_angle_convention+" "*40 +"|"
    print_statistics.show_rigid_body_rotations_and_translations(
      out=out,
      prefix="",
      frame="|",
      euler_angle_convention= self.params.euler_angle_convention,
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
               euler_angle_convention,
               lbfgs_maxfev):
    adopt_init_args(self, locals())
    self.fmodel_copy = self.fmodel.deep_copy()
    self.target_functor = self.fmodel_copy.target_functor()
    self.target_functor.prepare_for_minimization()
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

# XXX this is a lot dumber than the method name implies - it simply looks for
# protein or nucleic acid chains, and extracts the atom selection for each
def identify_rigid_groups (pdb_hierarchy) :
  assert (not pdb_hierarchy.atoms().extract_i_seq().all_eq(0))
  model = pdb_hierarchy.models()[0]
  atom_labels = list(pdb_hierarchy.atoms_with_labels())
  segids = flex.std_string([ a.segid for a in atom_labels ])
  use_segid = not segids.all_eq('    ')
  selections = []
  for chain in model.chains() :
    main_conf = chain.conformers()[0]
    if (main_conf.is_protein() or main_conf.is_na()) :
      chain_sele = "(chain '%s'" % chain.id
      if (use_segid) :
        first_atom_labels = main_conf.atoms()[0].fetch_labels()
        chain_sele += " and segid '%s'" % first_atom_labels.segid
      resseq_first = main_conf.residues()[0].resseq.strip()
      resseq_last = main_conf.residues()[-1].resseq.strip()
      chain_sele += " and resseq %s:%s)" % (resseq_first, resseq_last)
      selections.append(chain_sele)
  return selections
