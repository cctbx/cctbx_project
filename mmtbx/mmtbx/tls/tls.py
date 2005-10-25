from iotbx import pdb
import iotbx.pdb.interpretation
from cctbx.array_family import flex
import sys, math
from scitbx import lbfgs
from mmtbx.tls.tls import *
from mmtbx_tls_ext import *
from libtbx import adopt_init_args
import boost.python
ext = boost.python.import_ext("mmtbx_tls_ext")
from mmtbx_tls_ext import *

class tls_parameters(object):
   def __init__(self, T,
                      L,
                      S,
                      origin,
                      number_of_tls_groups,
                      residue_range):
     adopt_init_args(self, locals())
     assert len(self.residue_range) == self.number_of_tls_groups == len(self.T)
     assert len(self.T) == len(self.L) == len(self.S) == len(self.origin)
     self.n_groups = len(self.origin)

   def print_as_rows(self):
     print "Total Number of TLS Groups Extracted = ", self.number_of_tls_groups
     for i in xrange(self.number_of_tls_groups):
       T = self.T[i]
       L = self.L[i]
       S = self.S[i]
       print "TLS Group Number = ", i
       print "Residue Range    = ", self.residue_range[i]
       print "Origin (x,y,z)   = %9.4f %9.4f %9.4f"%\
             (self.origin[i][0],self.origin[i][1],self.origin[i][2])
       print "T11=%8.4f T22=%8.4f T33=%8.4f T12=%8.4f T13=%8.4f T23=%8.4f"%\
             (T[0],T[1],T[2],T[3],T[4],T[5])
       print "L11=%8.4f L22=%8.4f L33=%8.4f L12=%8.4f L13=%8.4f L23=%8.4f"%\
             (L[0],L[1],L[2],L[3],L[4],L[5])
       format1 = "S11=%8.4f S22=%8.4f S33=%8.4f S12=%8.4f S13=%8.4f "
       format2 = "S23=%8.4f S21=%8.4f S31=%8.4f S32=%8.4f"
       print (format1+format2)% (S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[8])

   def print_as_matrices(self):
     print "Total number of TLS Groups %d: "% self.number_of_tls_groups
     for i in xrange(self.number_of_tls_groups):
       T = self.T[i]
       L = self.L[i]
       S = self.S[i]
       lims = self.residue_range[i]
       print
       print "   TLS Group Number %d: " % i
       print "      Residue Range    = ",lims[0]," ",lims[1]," : ",lims[2],\
                                                                    " ",lims[3]
       print "      Origin (x,y,z)   = %9.4f %9.4f %9.4f"%\
             (self.origin[i][0],self.origin[i][1],self.origin[i][2])
       print "      Tensor T:                             "
       print "      T11=%8.4f T12=%8.4f T13=%8.4f" % (T[0],T[3],T[4])
       print "      T21=%8.4f T22=%8.4f T23=%8.4f" % (T[3],T[1],T[5])
       print "      T31=%8.4f T32=%8.4f T33=%8.4f" % (T[4],T[5],T[2])
       print "      Tensor L:"
       print "      L11=%8.4f L12=%8.4f L13=%8.4f" % (L[0],L[3],L[4])
       print "      L21=%8.4f L22=%8.4f L23=%8.4f" % (L[3],L[1],L[5])
       print "      L31=%8.4f L32=%8.4f L33=%8.4f" % (L[4],L[5],L[2])
       print "      Tensor S:"
       print "      S11=%8.4f S12=%8.4f S13=%8.4f" % (S[0],S[1],S[2])
       print "      S21=%8.4f S22=%8.4f S23=%8.4f" % (S[3],S[4],S[5])
       print "      S31=%8.4f S32=%8.4f S33=%8.4f" % (S[6],S[7],S[8])


class tls(object):
   def __init__(self, tls_parameters,
                      xray_structure,
                      stage_1):
     T = tls_parameters.T
     L = tls_parameters.L
     S = tls_parameters.S
     origin = tls_parameters.origin
     number_of_tls_groups = tls_parameters.number_of_tls_groups
     residue_range = tls_parameters.residue_range
     adopt_init_args(self, locals())
     self.u_from_pdb_grouped = []
     self.xray_structures = []
     self.u_aniso_from_tls = []
     self.tls_groups_as_xray_structures()

   def tls_groups_as_xray_structures(self):
     for rrange in self.residue_range:
       assert rrange[0] == rrange[2]
       chainID = rrange[0]
       selection_cache = self.stage_1.selection_cache()
       chain_sel = selection_cache.get_chainID(pattern = chainID)
       if(len(chain_sel) != 0):
         resNumStart = int(rrange[1])
         resNumEnd = int(rrange[3])
         res_sel = []
         for i in xrange(resNumStart,resNumEnd+1):
           res_sel.extend(selection_cache.get_resSeq(i=i))
         bool_sel = selection_cache.intersection(chain_sel) & \
                    selection_cache.union(res_sel)
         sel_xray_structure = self.xray_structure.select(bool_sel.iselection())
         self.xray_structures.append(sel_xray_structure)
         sel_atom_attributes = flex.select(self.stage_1.atom_attributes_list,
                                           flags = bool_sel)
         u_pdb = []
         for atom in sel_atom_attributes:
           if(atom.Ucart is not None):
             u_pdb.append(atom.Ucart)
           else:
             u_pdb.append([0.,0.,0.,0.,0.,0.])
         self.u_from_pdb_grouped.append(u_pdb)
     assert len(self.T) == len(self.u_from_pdb_grouped)
     for xray_structure,updb in zip(self.xray_structures,self.u_from_pdb_grouped):
       assert len(updb) == xray_structure.scatterers().size()
     print
     print "Number of Selected Xray Structures = ", len(self.xray_structures)
     assert len(self.xray_structures) == len(self.T)
     for xray_structure in self.xray_structures:
       print
       xray_structure.show_summary()
       print "Center of Mass = ", xray_structure.center_of_mass()

   def u_from_tls_grouped(self):
     for xray_structure,T,L,S,origin in zip(self.xray_structures,self.T,self.L,self.S,self.origin):
       u_aniso = uaniso_from_tls_domain(T,L,S,origin,xray_structure.sites_cart())
       self.u_aniso_from_tls.append(u_aniso)
     return self.u_aniso_from_tls

   def u_from_tls(self):
     self.u_from_tls_grouped()
     u_from_tls = []
     for u_tls in self.u_aniso_from_tls:
       for u in u_tls:
         u_from_tls.append(u)
     return u_from_tls

   def tls_from_u(self):
     T_min = []
     L_min = []
     S_min = []
     for u,origin,xray_structure in zip(self.u_from_pdb_grouped,
                                        self.origin,
                                        self.xray_structures):
       T_initial = [0.5,0.5,0.5,0.5,0.5,0.5]
       L_initial = [1.0,1.0,1.0,1.0,1.0,1.0]
       S_initial = [0.0,0.5,0.5,0.5,0.0,0.5,0.5,0.5,0.0]
       stop_flag = 0
       target_stop = -1.0
       for i in xrange(10000):
         target_start = target_stop
         minimized = tls_from_uaniso_minimizer(
                  uaniso = u,
                  T_initial = T_initial,
                  L_initial = L_initial,
                  S_initial = S_initial,
                  refine_T = True,
                  refine_L = True,
                  refine_S = True,
                  origin = origin,
                  sites = flex.to_list(xray_structure.sites_cart()))
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
     tls_obj = tls_parameters(T = T_min,
                              L = L_min,
                              S = S_min,
                              origin = self.origin,
                              number_of_tls_groups = self.number_of_tls_groups,
                              residue_range = self.residue_range)
     tls_obj.print_as_matrices()
     return tls_obj

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
               min_iterations=1000,
               max_calls=10000):
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
                                  min_iterations = min_iterations,
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
                      xs.scattering_dict(),
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
