from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx import lbfgs
import copy

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
     r11 =  self.c_psi*self.c_phi
     r12 = -self.c_psi*self.s_phi
     r13 =  self.s_psi
     r21 =  self.c_the*self.s_phi + self.s_the*self.s_psi*self.c_phi
     r22 =  self.c_the*self.c_phi - self.s_the*self.s_psi*self.s_phi
     r23 = -self.s_the*self.c_psi
     r31 =  self.s_the*self.s_phi - self.c_the*self.s_psi*self.c_phi
     r32 =  self.s_the*self.c_phi + self.c_the*self.s_psi*self.s_phi
     r33 =  self.c_the*self.c_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_phi(self):
     r11 = -self.c_psi*self.s_phi
     r12 = -self.c_psi*self.c_phi
     r13 =  0.0
     r21 =  self.c_the*self.c_phi - self.s_the*self.s_psi*self.s_phi
     r22 = -self.c_the*self.s_phi - self.s_the*self.s_psi*self.c_phi
     r23 =  0.0
     r31 =  self.s_the*self.c_phi + self.c_the*self.s_psi*self.s_phi
     r32 = -self.s_the*self.s_phi + self.c_the*self.s_psi*self.c_phi
     r33 =  0.0
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_psi(self):
     r11 = -self.s_psi*self.c_phi
     r12 =  self.s_psi*self.s_phi
     r13 =  self.c_psi
     r21 =  self.s_the*self.c_psi*self.c_phi
     r22 = -self.s_the*self.c_psi*self.s_phi
     r23 =  self.s_the*self.s_psi
     r31 = -self.c_the*self.c_psi*self.c_phi
     r32 =  self.c_the*self.c_psi*self.s_phi
     r33 = -self.c_the*self.s_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_the(self):
     r11 =  0.0
     r12 =  0.0
     r13 =  0.0
     r21 = -self.s_the*self.s_phi+self.c_the*self.s_psi*self.c_phi
     r22 = -self.s_the*self.c_phi-self.c_the*self.s_psi*self.s_phi
     r23 = -self.c_the*self.c_psi
     r31 =  self.c_the*self.s_phi+self.s_the*self.s_psi*self.c_phi
     r32 =  self.c_the*self.c_phi-self.s_the*self.s_psi*self.s_phi
     r33 = -self.s_the*self.c_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

def setup_search_range(f, step, nref_min):
  d_max, d_min = f.d_max_min()
  r_maxs = [x/10. for x in range(int(d_min*10), int(d_max*10), int(step*10))]
  for r_max in r_maxs:
      nref = f.resolution_filter(d_max = d_max, d_min = r_max).data().size()
      if(nref < nref_min):
         sel = flex.double(r_maxs) < r_max
         r_maxs = list(flex.double(r_maxs).select(sel))
         break
  r_maxs.reverse()
  return r_maxs

class manager(object):
  def __init__(self, fmodel,
                     selections       = None,
                     refine_r         = True,
                     refine_t         = True,
                     r_initial        = None,
                     t_initial        = None,
                     resolution_step  = 1.5,
                     nref_min         = 100,
                     max_iterations   = 20,
                     convergence_test = True,
                     convergence_delta= 0.00001,
                     log              = None):
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
    if(nref_min < fmodel_copy.f_obs_w().data().size()):
       d_mins = setup_search_range(fmodel_copy.f_obs_w(),
                                   step     = resolution_step,
                                   nref_min = nref_min)
    else:
       d_mins = [fmodel_copy.f_obs_w().d_min()]
    print >> log, "High resolution cutoffs for grid search: ", d_mins
    for res in d_mins:
        print >> log
        xrs = fmodel_copy.xray_structure.deep_copy_scatterers()
        fmodel_copy = fmodel.resolution_filter(d_min = res)
        fmodel_copy.update_xray_structure(xray_structure = xrs,
                                          update_f_calc  = True)
        rworks = flex.double()
        for macro_cycle in xrange(min(int(res),5)):
            minimized = rigid_body_minimizer(fmodel         = fmodel_copy,
                                             selections     = selections,
                                             r_initial      = r_initial,
                                             t_initial      = t_initial,
                                             refine_r       = refine_r,
                                             refine_t       = refine_t,
                                             max_iterations = max_iterations)
            rotation_matrices = []
            translation_vectors = []
            for i in xrange(len(selections)):
                self.total_rotation[i] += flex.double(minimized.r_min[i])
                self.total_translation[i] += flex.double(minimized.t_min[i])
                rot_obj = rb_mat(phi = minimized.r_min[i][0],
                                 psi = minimized.r_min[i][1],
                                 the = minimized.r_min[i][2])
                rotation_matrices.append(rot_obj.rot_mat())
                translation_vectors.append(minimized.t_min[i])
            new_xrs = apply_transformation(
                         xray_structure      = minimized.fmodel.xray_structure,
                         rotation_matrices   = rotation_matrices,
                         translation_vectors = translation_vectors,
                         selections          = selections)
            fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                              update_f_calc  = True)
            rwork = minimized.fmodel.r_work()
            assert approx_equal(rwork, fmodel_copy.r_work())
            if(convergence_test):
               rworks.append(rwork)
               if(rworks.size() > 1):
                  size = rworks.size() - 1
                  if(abs(rworks[size]-rworks[size-1])<convergence_delta):
                     break
            self.show(f     = fmodel_copy.f_obs_w(),
                      r_mat = self.total_rotation,
                      t_vec = self.total_translation,
                      rw    = rwork,
                      tw    = minimized.fmodel.target_w(),
                      mc    = macro_cycle,
                      it    = minimized.counter,
                      ct    = convergence_test,
                      out   = log)
    fmodel.update_xray_structure(xray_structure = fmodel_copy.xray_structure,
                                 update_f_calc  = True)
    self.fmodel = fmodel

  def rotation(self):
    return self.total_rotation

  def translation(self):
    return self.total_translation

  def show(self, f,
                 r_mat,
                 t_vec,
                 rw,
                 tw,
                 mc,
                 it,
                 ct,
                 out = None):
    if(out is None): out = sys.stdout
    d_max, d_min = f.d_max_min()
    nref = f.data().size()
    mc = str(mc)
    it = str(it)
    part1 = "|-rigid body refinement (macro cycle = "
    part2 = "; iterations = "
    n = 77 - len(part1 + part2 + mc + it)
    part3 = ")"+"-"*n+"|"
    print >> out, part1 + mc + part2 + it + part3
    part1 = "| resolution range: "
    d_max = str("%.3f"%d_max)
    part2 = " - "
    d_min = str("%.3f"%d_min)
    part3 = " ("
    nref = str("%d"%nref)
    part4 = " reflections) r-work = "
    rw = str("%.6f"%rw)
    n = 78 - len(part1+d_max+part2+d_min+part3+nref+part4+rw)
    part5 = " "*n+"|"
    print >> out, part1+d_max+part2+d_min+part3+nref+part4+rw+part5
    for r,t in zip(r_mat,t_vec):
        part1 = "| rotation (deg.) = "
        part2 = str("%8.4f"%r[0])+" "+str("%8.4f"%r[1])+" "+str("%8.4f"%r[2])
        n = 53 - len(part1 + part2)
        part3 = " "*n + "target = "+str("%.6f"%tw)
        n = 78 - len(part1 + part2 + part3)
        part4 = " "*n+"|"
        print >> out, part1 + part2 + part3 + part4
        part1 = "| translation (A) = "
        part2 = str("%8.4f"%t[0])+" "+str("%8.4f"%t[1])+" "+str("%8.4f"%t[2])
        n = 53 - len(part1 + part2)
        if(ct): ct = "on"
        else:   ct = "off"
        part3 = " "*n + "convergence test = "+str("%s"%ct)
        n = 78 - len(part1 + part2 + part3)
        part4 = " "*n+"|"
        print >> out, part1 + part2 + part3 + part4
    print >> out, "|" +"-"*77+"|"

class rigid_body_minimizer(object):
  def __init__(self,
               fmodel,
               selections,
               r_initial,
               t_initial,
               refine_r,
               refine_t,
               max_iterations):
    adopt_init_args(self, locals())
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
               termination_params = lbfgs.termination_parameters(
                    max_iterations = max_iterations),
               exception_handling_params = lbfgs.exception_handling_parameters(
                    ignore_line_search_failed_step_at_lower_bound = True)
                              )
    self.compute_functional_and_gradients()
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

  def compute_functional_and_gradients(self):
    self.unpack_x()
    self.counter += 1
    rotation_matrices   = []
    translation_vectors = []
    rot_objs = []
    for i in xrange(self.n_groups):
        rot_obj = rb_mat(phi = self.r_min[i][0],
                         psi = self.r_min[i][1],
                         the = self.r_min[i][2])
        rotation_matrices.append(rot_obj.rot_mat())
        translation_vectors.append(self.t_min[i])
        rot_objs.append(rot_obj)
    new_xrs = apply_transformation(xray_structure = self.fmodel.xray_structure,
                                   rotation_matrices   = rotation_matrices,
                                   translation_vectors = translation_vectors,
                                   selections          = self.selections)
    fmodel_copy = self.fmodel.deep_copy()
    fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                      update_f_calc  = True)
    tg_obj = target_and_grads(fmodel  = fmodel_copy,
                              rot_objs = rot_objs,
                              selections = self.selections)
    del fmodel_copy
    self.f = tg_obj.target()
    self.g = self.pack( tg_obj.gradients_wrt_r(), tg_obj.gradients_wrt_t() )
    return self.f, self.g

def apply_transformation(xray_structure,
                         rotation_matrices,
                         translation_vectors,
                         selections):
  assert len(selections) == len(rotation_matrices)
  assert len(selections) == len(translation_vectors)
  new_sites = flex.vec3_double()
  for sel,rot,trans in zip(selections,rotation_matrices,translation_vectors):
      xrs = xray_structure.select(sel)
      cm_cart = xrs.center_of_mass()
      sites_cart = xrs.sites_cart()
      sites_cart_cm = sites_cart - cm_cart
      new_sites.extend(list(rot) * sites_cart_cm + trans + cm_cart)
  return xray_structure.replace_sites_cart(new_sites = new_sites)

class target_and_grads(object):
  def __init__(self, fmodel,
                     rot_objs,
                     selections):
    self.grads_wrt_r = []
    self.grads_wrt_t = []
    target_grads_wrt_xyz = fmodel.gradient_wrt_xyz()
    self.f = fmodel.target_w()
    for sel,rot_obj in zip(selections, rot_objs):
        xrs = fmodel.xray_structure.select(sel)
        cm_cart = xrs.center_of_mass()
        sites_cart = xrs.sites_cart()
        sites_cart_cm = sites_cart - cm_cart
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
