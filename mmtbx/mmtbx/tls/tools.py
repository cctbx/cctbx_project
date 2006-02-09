from iotbx import pdb
import iotbx.pdb.interpretation
from cctbx.array_family import flex
import sys, math, time
from scitbx import lbfgs
from mmtbx_tls_ext import *
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal


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
               min_iterations= 50,
               max_calls= 100,
               run_finite_differences_test = True):
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
    #scale = 8*math.pi**2
    #qqq = new_xrs.scatterers().extract_u_iso_or_u_equiv(new_xrs.unit_cell())*scale
    #ipd = new_xrs.is_positive_definite_u()
    #print "%10.2f %10.2f %10.2f %5d %5d" % (flex.max(qqq), flex.min(qqq), flex.mean(qqq), \
    #      ipd.count(True), ipd.count(False)),self.f


    self.g = self.pack(grad_manager.grad_T,
                       grad_manager.grad_L,
                       grad_manager.grad_S)
    if(self.run_finite_differences_test):
       GT,GL,GS = finite_differences_grads_of_xray_target_wrt_tls(
                                                 fmodel     = self.fmodel_copy,
                                                 T          = self.T_min,
                                                 L          = self.L_min,
                                                 S          = self.S_min,
                                                 origins    = self.origins,
                                                 selections = self.selections,
                                                 delta      = 0.00001)
       format   = "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f"
       formats="%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f"
       print "T"
       for m1,m2 in zip(grad_manager.grad_T, GT):
           print "1=" + format % (m1[0],m1[1],m1[2],m1[3],m1[4],m1[5])
           print "2=" + format % (m2[0],m2[1],m2[2],m2[3],m2[4],m2[5])
           assert approx_equal(m1,m2)
       print "L"
       for m1,m2 in zip(grad_manager.grad_L, GL):
           print "1=" + format % (m1[0],m1[1],m1[2],m1[3],m1[4],m1[5])
           print "2=" + format % (m2[0],m2[1],m2[2],m2[3],m2[4],m2[5])
           assert approx_equal(m1,m2)
       print "S"
       for m1,m2 in zip(grad_manager.grad_S, GS):
           print "1=" + formats %\
                        (m1[0],m1[1],m1[2],m1[3],m1[4],m1[5],m1[6],m1[7],m1[8])
           print "2=" + formats %\
                        (m2[0],m2[1],m2[2],m2[3],m2[4],m2[5],m2[6],m2[7],m2[8])
           assert approx_equal(m1,m2)
       print
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
            d_target_d_uaniso = d_target_d_uaniso.select(sel),
            scale_l_and_s     = True,
            use_trace_s_zero_constraint = True)
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
  #XXX refinement of only S or L does not work with this:
  #xray_structure.tidy_us(u_min = 1.e-6)
  return xray_structure

class tls_refinement(object):
   def __init__(self, fmodel,
                      selections,
                      refine_T,
                      refine_L,
                      refine_S,
                      number_of_macro_cycles,
                      max_number_of_iterations,
                      start_tls_value = None):
     fmodel.xray_structure.convert_to_anisotropic()
     fmodel.xray_structure.tidy_us(u_min = 1.e-6)
     origins = []
     for selection in selections:
         xrs = fmodel.xray_structure.select(selection)
         origins.append(xrs.center_of_mass())
     if(start_tls_value is not None):
        T_start, L_start,S_start = fill_tls(start_tls_value, len(selections))
     else:
        tlsos_ = []
        t_, l_, s_ = fill_tls(0.0, len(selections))
        for T_,L_,S_,o_ in zip(t_, l_, s_, origins):
            tlsos_.append(tlso(t = T_, l = L_, s = S_, origin = o_))
        tlsos_start = tls_from_uanisos(xray_structure = fmodel.xray_structure,
                                       selections     = selections,
                                       tlsos_initial  = tlsos_)
        print "*"*80
        show_tls(tlsos = tlsos_start)

def fill_tls(value=0.0, n=1):
  T = []
  L = []
  S = []
  v   = value
  vtl = [v,v,v,v,v,v]
  vs  = [v,v,v,v,v,v,v,v,v]
  for dummy in xrange(n):
      T.append(vtl)
      L.append(vtl)
      S.append(vs)
  return T,L,S

def finite_differences_grads_of_xray_target_wrt_tls(fmodel,
                                                    T,
                                                    L,
                                                    S,
                                                    origins,
                                                    selections,
                                                    delta=0.00001):
  fmodel_copy = fmodel.deep_copy()
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
        #
        tlsos = []
        for Ti,Li,Si,origini in zip(T_, L, S, origins):
            tlsos.append(tlso(t      = Ti,
                              l      = Li,
                              s      = Si,
                              origin = origini))
        new_xrs = update_xray_structure_with_tls(
                                  xray_structure = fmodel_copy.xray_structure,
                                  selections     = selections,
                                  tlsos          = tlsos)
        fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                          update_f_calc  = True)
        target_result = fmodel_copy.target_w()
        #
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
        #
        tlsos = []
        for Ti,Li,Si,origini in zip(T, L_, S, origins):
            tlsos.append(tlso(t      = Ti,
                              l      = Li,
                              s      = Si,
                              origin = origini))
        new_xrs = update_xray_structure_with_tls(
                                  xray_structure = fmodel_copy.xray_structure,
                                  selections     = selections,
                                  tlsos          = tlsos)
        fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                          update_f_calc  = True)
        target_result = fmodel_copy.target_w()
        #
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
        #
        tlsos = []
        for Ti,Li,Si,origini in zip(T, L, S_, origins):
            tlsos.append(tlso(t      = Ti,
                              l      = Li,
                              s      = Si,
                              origin = origini))
        new_xrs = update_xray_structure_with_tls(
                                  xray_structure = fmodel_copy.xray_structure,
                                  selections     = selections,
                                  tlsos          = tlsos)
        fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                          update_f_calc  = True)
        target_result = fmodel_copy.target_w()
        #
        target_values.append(target_result)
      derivative = (target_values[1] - target_values[0]) / (2 * delta)
      dS.append(derivative)
    derivative_S.append(dS)
  return derivative_T,derivative_L,derivative_S
