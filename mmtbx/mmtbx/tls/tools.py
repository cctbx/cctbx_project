from iotbx import pdb
import iotbx.pdb.interpretation
from cctbx.array_family import flex
import sys, math
from scitbx import lbfgs
from mmtbx_tls_ext import *
from libtbx import adopt_init_args

class show_tls(object):
   def __init__(self, tlsos):
     print "Total number of TLS groups %d: "% len(tlsos)
     counter = 0
     for item in tlsos:
         counter += 1
         T = item.t
         L = item.l
         S = item.s
         origin = item.origin
         print "   TLS group number %d: " % counter
         print "      Origin (x,y,z)   = %9.4f %9.4f %9.4f"%\
                                                (origin[0],origin[1],origin[2])
         print "      T:"
         print "      T11=%8.4f T12=%8.4f T13=%8.4f" % (T[0],T[3],T[4])
         print "      T21=%8.4f T22=%8.4f T23=%8.4f" % (T[3],T[1],T[5])
         print "      T31=%8.4f T32=%8.4f T33=%8.4f" % (T[4],T[5],T[2])
         print "      L:"
         print "      L11=%8.4f L12=%8.4f L13=%8.4f" % (L[0],L[3],L[4])
         print "      L21=%8.4f L22=%8.4f L23=%8.4f" % (L[3],L[1],L[5])
         print "      L31=%8.4f L32=%8.4f L33=%8.4f" % (L[4],L[5],L[2])
         print "      S:"
         print "      S11=%8.4f S12=%8.4f S13=%8.4f" % (S[0],S[1],S[2])
         print "      S21=%8.4f S22=%8.4f S23=%8.4f" % (S[3],S[4],S[5])
         print "      S31=%8.4f S32=%8.4f S33=%8.4f" % (S[6],S[7],S[8])

def uanisos_from_tls(sites_cart, selections, tlsos):
  uanisos = flex.sym_mat3_double()
  for selection, tlso in zip(selections, tlsos):
      uanisos.extend(
          uaniso_from_tls_one_group(tlso = tlso,
                                    sites_cart = sites_cart.select(selection)))
  return uanisos

def tls_from_uanisos(xray_structure, selections, tlsos_initial):
  assert xray_structure.scatterers().count_anisotropic() == \
                                             xray_structure.scatterers().size()
  u_cart=xray_structure.scatterers().extract_u_cart(xray_structure.unit_cell())
  T_min = []
  L_min = []
  S_min = []
  for tlso_initial, selection in zip(tlsos_initial, selections):
      T_initial = tlso_initial.t
      L_initial = tlso_initial.l
      S_initial = tlso_initial.s
      stop_flag = 0
      target_stop = -1.0
      sites_cart_selected = xray_structure.sites_cart().select(selection)
      u_cart_selected = u_cart.select(selection)
      for i in xrange(1000):
        target_start = target_stop
        minimized = tls_from_uaniso_minimizer(uaniso    = u_cart_selected,
                                              T_initial = T_initial,
                                              L_initial = L_initial,
                                              S_initial = S_initial,
                                              refine_T  = True,
                                              refine_L  = True,
                                              refine_S  = True,
                                              origin    = tlso_initial.origin,
                                              sites     = sites_cart_selected)
        T_initial = minimized.T_min
        L_initial = minimized.L_min
        S_initial = minimized.S_min
        target_stop = minimized.f
        if(abs(target_stop - target_start) < 1.e-30): stop_flag += 1
        if(stop_flag == 3):
          print "convergence at step ", i
          break
      T_min.append(minimized.T_min)
      L_min.append(minimized.L_min)
      S_min.append(minimized.S_min)
  tlsos_result = []
  for T,L,S,tlso_initial_ in zip(T_min,L_min,S_min, tlsos_initial):
      tlsos_result.append(tlso(t = T,
                                l = L,
                                s = S,
                                origin = tlso_initial_.origin))
  return tlsos_result

class tls_from_uaniso_minimizer(object):
  def __init__(self,
               uaniso,
               T_initial,
               L_initial,
               S_initial,
               refine_T,
               refine_L,
               refine_S,
               origin,
               sites,
               max_iterations=1000,
               max_calls=1000):
    adopt_init_args(self, locals())
    self.dim_T = len(self.T_initial)
    self.dim_L = len(self.L_initial)
    self.dim_S = len(self.S_initial)
    assert self.dim_T == 6 and self.dim_S == 9
    self.T_min = self.T_initial
    self.L_min = self.L_initial
    self.S_min = self.S_initial
    self.x = self.pack(self.T_min, self.L_min, self.S_min)
    self.n = self.x.size()
    self.minimizer = lbfgs.run(
                             target_evaluator = self,
                             termination_params = lbfgs.termination_parameters(
                                  max_iterations = max_iterations,
                                  max_calls = max_calls)
                              )
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, T, L, S):
    v = []
    if (self.refine_T): v += list(T)
    if (self.refine_L): v += list(L)
    if (self.refine_S): v += list(S)
    return flex.double(tuple(v))

  def unpack_x(self):
    i = 0
    if (self.refine_T):
      self.T_min = tuple(self.x)[i:self.dim_T]
      i = self.dim_T
    if (self.refine_L):
      self.L_min = tuple(self.x)[i:i+self.dim_L]
      i += self.dim_L
    if (self.refine_S):
      self.S_min = tuple(self.x)[i:i+self.dim_S]

  def compute_functional_and_gradients(self):
    self.unpack_x()
    manager = tls_from_uaniso_target_and_grads(self.T_min,
                                               self.L_min,
                                               self.S_min,
                                               self.origin,
                                               self.sites,
                                               self.uaniso)
    self.f = manager.target()
    self.g = self.pack(manager.grad_T(), manager.grad_L(), manager.grad_S())
    return self.f, self.g


class tls_xray_target_minimizer(object):
  def __init__(self,
               fmodel,
               T_initial,
               L_initial,
               S_initial,
               refine_T,
               refine_L,
               refine_S,
               origins,
               selections,
               min_iterations= 25,
               max_calls= 50):
    adopt_init_args(self, locals())
    self.fmodel_copy = self.fmodel.deep_copy()
    self.counter = 0
    assert len(self.T_initial) == len(self.L_initial) == len(self.S_initial)
    for Ti,Li,Si in zip(self.T_initial,self.L_initial,self.S_initial):
      assert len(Ti) == 6 and len(Li) == 6 and len(Si) == 9
    self.n_groups = len(self.T_initial)
    self.dim_T = len(self.T_initial[0])
    self.dim_L = len(self.L_initial[0])
    self.dim_S = len(self.S_initial[0])
    assert self.dim_T == 6 and self.dim_L == 6 and self.dim_S == 9
    self.T_min = self.T_initial
    self.L_min = self.L_initial
    self.S_min = self.S_initial
    self.x = self.pack(self.T_min, self.L_min, self.S_min)
    self.minimizer = lbfgs.run(
         target_evaluator = self,
         termination_params = lbfgs.termination_parameters(
              min_iterations = min_iterations,
              max_calls = max_calls))
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, T, L, S):
    v = []
    for Ti,Li,Si in zip(T,L,S):
      if (self.refine_T): v += list(Ti)
      if (self.refine_L): v += list(Li)
      if (self.refine_S): v += list(Si)
    return flex.double(tuple(v))

  def unpack_x(self):
    i = 0
    T_min = []
    L_min = []
    S_min = []
    for j in xrange(self.n_groups):
      if (self.refine_T):
        self.T_min[j] = tuple(self.x)[i:i+self.dim_T]
        i += self.dim_T
      if (self.refine_L):
        self.L_min[j] = tuple(self.x)[i:i+self.dim_L]
        i += self.dim_L
      if (self.refine_S):
        self.S_min[j] = tuple(self.x)[i:i+self.dim_S]
        i += self.dim_S


  def compute_functional_and_gradients(self):
    self.counter += 1
    self.unpack_x()

    S_new = []
    for item_1 in self.S_min:
      zz = list(item_1)
      zz[8] = -zz[0]-zz[4]
      S_new.append(zz)
    self.S_min = S_new

    tlsos = []
    for T,L,S,origin in zip(self.T_min, self.L_min, self.S_min, self.origins):
        tlsos.append(tlso(t      = T,
                          l      = L,
                          s      = S,
                          origin = origin))

    new_xrs = update_xray_structure_with_tls(
                              xray_structure = self.fmodel_copy.xray_structure,
                              selections     = self.selections,
                              tlsos          = tlsos)

    self.fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                           update_f_calc  = True)

    grad_manager = tls_xray_grads(fmodel     = self.fmodel_copy,
                                  selections = self.selections,
                                  tlsos      = tlsos)

    self.f = self.fmodel_copy.target_w()
    print self.f
    self.g = self.pack(grad_manager.grad_T,
                       grad_manager.grad_L,
                       grad_manager.grad_S)
    return self.f, self.g

class tls_xray_grads(object):
   def __init__(self, fmodel, selections, tlsos):
     self.grad_T = []
     self.grad_L = []
     self.grad_S = []
     d_target_d_uaniso = fmodel.gradient_wrt_u_aniso()
     for sel, tlso in zip(selections, tlsos):
         d_target_d_tls_manager = d_target_d_tls(
            sites             = fmodel.xray_structure.sites_cart().select(sel),
            origin            = tlso.origin,
            d_target_d_uaniso = d_target_d_uaniso.select(sel))
         self.grad_T.append(list(d_target_d_tls_manager.grad_T()))
         self.grad_L.append(list(d_target_d_tls_manager.grad_L()))
         self.grad_S.append(list(d_target_d_tls_manager.grad_S()))


def update_xray_structure_with_tls(xray_structure,
                                   selections,
                                   tlsos):
  u_cart_from_tls = uanisos_from_tls(sites_cart = xray_structure.sites_cart(),
                                     selections = selections,
                                     tlsos      = tlsos)
  xray_structure.scatterers().set_u_cart(xray_structure.unit_cell(),
                                         u_cart_from_tls)
  xray_structure.apply_symmetry_u_stars()
  return xray_structure



############################## OLD CODE: MUST BE REDONE
class tls_xray_target_minimizer__OLD(object):
  def __init__(self,
               fobs,
               T_initial,
               L_initial,
               S_initial,
               refine_T,
               refine_L,
               refine_S,
               origins,
               xray_structures,
               min_iterations=50,
               max_calls=100):
    adopt_init_args(self, locals())
    self.counter = 0
    assert len(self.T_initial) == len(self.L_initial) == len(self.S_initial)
    assert len(self.T_initial) == len(xray_structures) == len(origins)
    for Ti,Li,Si in zip(self.T_initial,self.L_initial,self.S_initial):
      assert len(Ti) == 6 and len(Li) == 6 and len(Si) == 9
    self.n_groups = len(self.T_initial)
    self.dim_T = len(self.T_initial[0])
    self.dim_L = len(self.L_initial[0])
    self.dim_S = len(self.S_initial[0])
    assert self.dim_T == 6 and self.dim_L == 6 and self.dim_S == 9
    self.T_min = self.T_initial
    self.L_min = self.L_initial
    self.S_min = self.S_initial
    self.x = self.pack(self.T_min, self.L_min, self.S_min)
    self.minimizer = lbfgs.run(
         target_evaluator = self,
         termination_params = lbfgs.termination_parameters(
              min_iterations = min_iterations,
              max_calls = max_calls))
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, T, L, S):
    v = []
    for Ti,Li,Si in zip(T,L,S):
      if (self.refine_T): v += list(Ti)
      if (self.refine_L): v += list(Li)
      if (self.refine_S): v += list(Si)
    return flex.double(tuple(v))

  def unpack_x(self):
    i = 0
    T_min = []
    L_min = []
    S_min = []
    for j in xrange(self.n_groups):
      if (self.refine_T):
        self.T_min[j] = tuple(self.x)[i:i+self.dim_T]
        i += self.dim_T
      if (self.refine_L):
        self.L_min[j] = tuple(self.x)[i:i+self.dim_L]
        i += self.dim_L
      if (self.refine_S):
        self.S_min[j] = tuple(self.x)[i:i+self.dim_S]
        i += self.dim_S


  def compute_functional_and_gradients(self):
    self.counter += 1
    self.unpack_x()

    S_new = []
    for item_1 in self.S_min:
      zz = list(item_1)
      zz[8] = -zz[0]-zz[4]
      S_new.append(zz)
    self.S_min = S_new

    manager = tls_fcalc_target_grads(self.xray_structures,
                                     self.fobs,
                                     self.origins,
                                     self.T_min,
                                     self.L_min,
                                     self.S_min)
    self.f = manager.target

    #print manager.grad_T
    #print manager.grad_L
    #print manager.grad_S


    #gt = flex.double()
    #for i in manager.grad_T:
    #  gt.extend(flex.double(i))
    #normT = math.sqrt(flex.sum(flex.pow2(gt)))
    #
    #gl = flex.double()
    #for i in manager.grad_L:
    #  gl.extend(flex.double(i))
    #normL = math.sqrt(flex.sum(flex.pow2(gl)))
    #
    #gs = flex.double()
    #for i in manager.grad_S:
    #  gs.extend(flex.double(i))
    #normS = math.sqrt(flex.sum(flex.pow2(gs)))
    #
    #print normT, normL, normS, normT/normL, normT/normS

    t = manager.grad_L
    #x = []
    #j = []
    #for i in t:
    #  x = tuple(flex.double(i) * (normT/normL) )
    #  j.append(x)
    #t = j

    s = manager.grad_S
    #x = []
    #j = []
    #for i in s:
    #  x = tuple(flex.double(i) * (normT/normS) )
    #  j.append(x)
    #s = j

    #print
    #print manager.grad_T
    #print tuple(t)
    #print tuple(s)


    self.g = self.pack(manager.grad_T, t, s)
    return self.f, self.g


class tls_fcalc_target_grads(object):
    def __init__(self, xray_structures,
                       fobs,
                       origins,
                       T,
                       L,
                       S):
      assert len(xray_structures) == len(T) == len(L) == len(S) == len(origins)
      unit_cell = xray_structures[0].unit_cell()
      xs_modified = []
      tls_group_flags = flex.int()
      utls = flex.sym_mat3_double()
      tmp = flex.vec3_double(3)
      counter = 1
      for Ti,Li,Si,oi,xsi in zip(T,L,S,origins,xray_structures):
        Si_ = list(Si)
        Si_[8] = -Si_[0]-Si_[4]
        Si = Si_
        assert len(Ti)==len(Li)==6 and len(Si)==9 and len(oi)==3
        u_from_tls = uaniso_from_tls_domain(Ti,Li,Si,oi,xsi.sites_cart())
        utls.extend(u_from_tls)
        xsi_ = xsi.deep_copy_scatterers()
        xsi_.scatterers().set_u_cart(unit_cell, u_from_tls)
        xs_modified.append(xsi_)
        tls_group_flags.extend( flex.int(xsi.sites_cart().size(),counter) )
        tmp[counter-1] = oi
        counter += 1

      xs = xs_modified[0]
      for i in range(1,len(xs_modified)):
        xs.add_scatterers(xs_modified[i].scatterers())
      xs.apply_symmetry_u_stars()
      fcalc = fobs.structure_factors_from_scatterers(
                                 xray_structure = xs,
                                 algorithm      = "direct").f_calc()
      manager_fast = tls_xray_target_grads(
                      fobs.indices(),
                      flex.abs(fobs.data()),
                      fcalc.data(),
                      unit_cell,
                      tmp,#origins,
                      utls,
                      xs.scattering_type_registry(),
                      xs.scatterers(),
                      tls_group_flags)
      self.target = manager_fast.target()
      self.grad_T = []
      self.grad_L = []
      self.grad_S = []
      gradTLS = manager_fast.gradTLS()
      zz=0
      for i in xrange(len(origins)):
        self.grad_T.append(list(gradTLS[i].t))
        self.grad_L.append(list(gradTLS[i].l))
        self.grad_S.append(list(gradTLS[i].s))
        if(0):
          Si = S[i]
          sum = (Si[0]+Si[4]+Si[8])
          delta_s =  sum**2
          delta_sg = 2.*sum
          self.target += delta_s
          self.grad_S[zz][0] += delta_sg
          self.grad_S[zz][4] += delta_sg
          self.grad_S[zz][8] += delta_sg
          zz+=1

def finite_differences_grads_of_xray_target_wrt_tls(xray_structure,
                                                    T,
                                                    L,
                                                    S,
                                                    fobs,
                                                    origin,
                                                    delta=0.00001):
  derivative_T = []
  for j in xrange(len(T)):
    dT = []
    for i in xrange(6):
      target_values = []
      for d_sign in (-1, 1):
        T_ = []
        for item in T:
          T_.append(list(item))
        d = d_sign*delta
        T_[j][i] += d
        target_result = tls_fcalc_target_grads(
                                 xray_structure, fobs, origin, T_, L, S).target
        target_values.append(target_result)
      derivative = (target_values[1] - target_values[0]) / (2 * delta)
      dT.append(derivative)
    derivative_T.append(dT)

  derivative_L = []
  for j in xrange(len(L)):
    dL = []
    for i in xrange(6):
      target_values = []
      for d_sign in (-1, 1):
        L_ = []
        for item in L:
          L_.append(list(item))
        d = d_sign*delta
        L_[j][i] += d
        target_result = tls_fcalc_target_grads(
                                 xray_structure, fobs, origin, T, L_, S).target
        target_values.append(target_result)
      derivative = (target_values[1] - target_values[0]) / (2 * delta)
      dL.append(derivative)
    derivative_L.append(dL)

  derivative_S = []
  for j in xrange(len(L)):
    dS = []
    for i in xrange(9):
      target_values = []
      for d_sign in (-1, 1):
        S_ = []
        for item in S:
          S_.append(list(item))
        d = d_sign*delta
        S_[j][i] += d
        target_result = tls_fcalc_target_grads(
                                 xray_structure, fobs, origin, T, L, S_).target
        target_values.append(target_result)
      derivative = (target_values[1] - target_values[0]) / (2 * delta)
      dS.append(derivative)
    derivative_S.append(dS)

  return derivative_T,derivative_L,derivative_S
