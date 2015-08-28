from __future__ import division
from scitbx import matrix
import math
from scitbx.linalg import eigensystem
from libtbx.test_utils import approx_equal
from libtbx import group_args
from scitbx.array_family import flex
import iotbx.pdb
from libtbx.utils import Sorry
import sys, os
from copy import deepcopy
import mmtbx.tls.tools
from libtbx.utils import multi_out

def print_step(s, log):
  n = 79-len(s)
  print >> log, s, "*"*n

def set_log(prefix, i_current, i_total):
  log = multi_out()
  fo = open("%s.tls.%s_of_%s"%(prefix,str(i_current),i_total),"w")
  log.register(
    label       = "log_buffer",
    file_object = fo)
  sys.stderr = log
  return log

def cmd_driver(pdb_file_name):
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  xrs = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry=pdb_inp.crystal_symmetry_from_cryst1())
  tls_extract = mmtbx.tls.tools.tls_from_pdb_inp(
    remark_3_records = pdb_inp.extract_remark_iii_records(3),
    pdb_hierarchy    = pdb_hierarchy)
  tlsos = tls_extract.tls_params
  for i_seq, tlso in enumerate(tlsos, 1):
    log = set_log(prefix=os.path.basename(pdb_file_name), i_current=i_seq,
      i_total=len(tlsos))
    deg_to_rad_scale = math.pi/180
    # Units: T[A], L[deg**2], S[A*deg]
    T = matrix.sym(sym_mat3=tlso.t)
    L = matrix.sym(sym_mat3=tlso.l)*(deg_to_rad_scale**2)
    S = matrix.sqr(tlso.s)*deg_to_rad_scale
    try:
      r = run(T=T, L=L, S=S, log=log)
    except Exception, e:
      print >> log, str(e)
    log.close()

def truncate(m, eps_string="%.8f"):
  if(type(m) is float):
    return float(eps_string%m)
  elif(type(m) is flex.double):
    return flex.double([float(eps_string%i) for i in m])
  else:
    if(not type(m) is matrix.sqr):
      raise Sorry("truncate: arg must be matrix.sqr")
    x = [m[0],m[1],m[2], m[3],m[4],m[5], m[6],m[7],m[8]]
    x_ = []
    for xi in x: x_.append(float(eps_string%xi))
    return matrix.sqr(
      [x_[0], x_[1], x_[2],
       x_[3], x_[4], x_[5],
       x_[6], x_[7], x_[8]])

class run(object):
  def __init__(self, T, L, S, log=sys.stdout, eps=1.e-6, self_check_eps=1.e-5,
               force_t_S=None):
    """
    Decompose TLS matrices into physically iterpretable components.
    """
    self.log = log
    self.ff = "%12.9f"
    self.eps = eps
    self.self_check_eps = self_check_eps
    self.force_t_S = force_t_S
    print >> self.log, "Small is defined as:", self.eps
    self.T_M, self.L_M, self.S_M = T, L, S
    print_step("Input TLS matrices:", self.log)
    self.show_matrix(x=self.T_M, prefix="  ", title="T_M")
    self.show_matrix(x=self.L_M, prefix="  ", title="L_M")
    self.show_matrix(x=self.S_M, prefix="  ", title="S_M")
    # set at step A
    self.T_L, self.L_L, self.S_L = None, None, None
    self.l_x, self.l_y, self.l_z = None, None, None
    self.R_ML = None
    self.Lxx, self.Lyy, self.Lzz = None, None, None
    self.Sxx, self.Syy, self.Szz = None, None, None
    # set at step B
    self.w_15 = None # formula (15)
    self.w = None # formula (16)
    self.T_CL = None
    self.D_WL = None
    # set at step C
    self.t_S = None
    self.sx = 0
    self.sy = 0
    self.sz = 0
    self.C_L_t_S = None
    self.V_L = None
    self.S_C = None
    # D
    self.V_V = None
    self.R_MV = None
    self.v_x, self.v_y, self.v_z = None, None, None
    self.v_x_M, self.v_y_M, self.v_z_M = None, None, None # formula (40)
    self.tx, self.ty, self.tz = None, None, None
    # all steps
    self.step_A()
    self.step_B()
    self.step_C()
    self.step_D()
    # compose result-object
    self.result = self.finalize()
    self.show_summary()
    # consistency check: restore input T_M, L_M, S_M from base elements
    self.self_check()

  def step_A(self):
    """
    Diagonalization of L matrix.
    """
    # check input matrices T_M>=0 and L_M>=0:
    if(not self.is_pd(self.T_M.as_sym_mat3())):
      raise Sorry("Step A: Input matrix T[M] is not positive semidefinite.")
    if(not self.is_pd(self.L_M.as_sym_mat3())):
      raise Sorry("Step A: Input matrix L[M] is not positive semidefinite.")
    print_step("Step A:", self.log)
    es = self.eigen_system_default_handler(m=self.L_M, suffix="L_M")
    self.l_x, self.l_y, self.l_z = es.x, es.y, es.z
    self.show_vector(x=self.l_x, title="l_x")
    self.show_vector(x=self.l_y, title="l_y")
    self.show_vector(x=self.l_z, title="l_z")
    self.R_ML = matrix.sqr(
      [self.l_x[0], self.l_y[0], self.l_z[0],
       self.l_x[1], self.l_y[1], self.l_z[1],
       self.l_x[2], self.l_y[2], self.l_z[2]])
    self.show_matrix(x=self.R_ML, title="Matrix R_ML, eq.(14)")
    R_ML_transpose = self.R_ML.transpose()
    assert approx_equal(R_ML_transpose, self.R_ML.inverse(), 1.e-5)
    self.T_L = truncate(m=R_ML_transpose*self.T_M*self.R_ML)
    self.L_L = truncate(m=R_ML_transpose*self.L_M*self.R_ML)
    self.S_L = truncate(m=R_ML_transpose*self.S_M*self.R_ML)
    self.show_matrix(x=self.T_L, title="T_L, eq.(13)")
    self.show_matrix(x=self.L_L, title="L_L, eq.(13)")
    self.show_matrix(x=self.S_L, title="S_L, eq.(13)")
    L_ = self.L_L.as_sym_mat3()
    self.Lxx, self.Lyy, self.Lzz = L_[0], L_[1], L_[2]
    self.Sxx, self.Syy, self.Szz = self.S_L[0], self.S_L[4], self.S_L[8]

  def step_B(self):
    """
    Position of the libration axes in the L-basis.
    """
    print_step("Step B:", self.log)
    wy_lx=0
    wz_lx=0
    wx_ly=0
    wz_ly=0
    wx_lz=0
    wy_lz=0
    if(self.is_zero(self.Lxx)):
      if(not (self.is_zero(self.S_L[2]) and self.is_zero(self.S_L[1]))):
        raise Sorry("Step B: Non-zero off-diagonal S[L] and zero L[L] elements.")
    else:
      wy_lx =-self.S_L[2]/self.Lxx
      wz_lx = self.S_L[1]/self.Lxx
    if(self.is_zero(self.Lyy)):
      if(not (self.is_zero(self.S_L[5]) and self.is_zero(self.S_L[3]))):
        raise Sorry("Step B: Non-zero off-diagonal S[L] and zero L[L] elements.")
    else:
      wx_ly = self.S_L[5]/self.Lyy
      wz_ly =-self.S_L[3]/self.Lyy
    if(self.is_zero(self.Lzz)):
      if(not (self.is_zero(self.S_L[7]) and self.is_zero(self.S_L[6]))):
        raise Sorry("Step B: Non-zero off-diagonal S[L] and zero L[L] elements.")
    else:
      wx_lz =-self.S_L[7]/self.Lzz
      wy_lz = self.S_L[6]/self.Lzz
    # formuala (15)...
    w_lx = matrix.col((0,     wy_lx, wz_lx))
    w_ly = matrix.col((wx_ly,     0, wz_ly))
    w_lz = matrix.col((wx_lz, wy_lz,     0))
    self.w_15 = group_args(
      wy_lx = wy_lx,
      wz_lx = wz_lx,
      wx_ly = wx_ly,
      wz_ly = wz_ly,
      wx_lz = wx_lz,
      wy_lz = wy_lz,
      w_lx = w_lx,
      w_ly = w_ly,
      w_lz = w_lz)
    self.show_vector(x=w_lx, title="w_lx, eq.(15)")
    self.show_vector(x=w_ly, title="w_ly, eq.(15)")
    self.show_vector(x=w_lz, title="w_lz, eq.(15)")
    # ... but really we want (16)
    w_lx = matrix.col(((wx_ly+wx_lz)/2.,            wy_lx,            wz_lx)) # overwrite inplace
    w_ly = matrix.col((           wx_ly, (wy_lx+wy_lz)/2.,            wz_ly)) # overwrite inplace
    w_lz = matrix.col((           wx_lz,            wy_lz, (wz_lx+wz_ly)/2.)) # overwrite inplace
    self.w = group_args(
      wy_lx = wy_lx,
      wz_lx = wz_lx,
      wx_ly = wx_ly,
      wz_ly = wz_ly,
      wx_lz = wx_lz,
      wy_lz = wy_lz,
      w_lx = w_lx,
      w_ly = w_ly,
      w_lz = w_lz)
    self.show_vector(x=w_lx, title="w_lx, eq.(16)")
    self.show_vector(x=w_ly, title="w_ly, eq.(16)")
    self.show_vector(x=w_lz, title="w_lz, eq.(16)")
    #
    d11 = wz_ly**2*self.Lyy + wy_lz**2*self.Lzz
    d22 = wz_lx**2*self.Lxx + wx_lz**2*self.Lzz
    d33 = wy_lx**2*self.Lxx + wx_ly**2*self.Lyy
    d12 = -wx_lz*wy_lz*self.Lzz
    d13 = -wx_ly*wz_ly*self.Lyy
    d23 = -wy_lx*wz_lx*self.Lxx
    self.D_WL = matrix.sqr(
      [d11, d12, d13,
       d12, d22, d23,
       d13, d23, d33])
    self.show_matrix(x=self.D_WL, title="D_WL, eq.(10)")
    self.T_CL = self.T_L - self.D_WL
    self.T_CL = truncate(matrix.sqr(self.T_CL))
    self.show_matrix(x=self.T_CL, title="T_CL")
    if(not self.is_pd(self.T_CL.as_sym_mat3())):
      raise Sorry("Step B: Matrix T_C[L] is not positive semidefinite.")

  def step_C(self):
    """
    Determination of screw components.
    """
    print_step("Step C:", self.log)
    T_ = self.T_CL.as_sym_mat3()
    self.T_CLxx, self.T_CLyy, self.T_CLzz = T_[0], T_[1], T_[2]
    tlxx = self.T_CLxx*self.Lxx
    tlyy = self.T_CLyy*self.Lyy
    tlzz = self.T_CLzz*self.Lzz
    if(tlxx < 0 or tlyy < 0 or tlzz < 0):
      raise Sorry("Step C: Arguments in eq.(24) must be non-negative.")
    rx, ry, rz = math.sqrt(tlxx), math.sqrt(tlyy), math.sqrt(tlzz)
    #
    t_min_C = truncate(m=max(self.Sxx-rx, self.Syy-ry, self.Szz-rz))
    t_max_C = truncate(m=min(self.Sxx+rx, self.Syy+ry, self.Szz+rz))
    self.show_number(x=[t_min_C, t_max_C], title="t_min_C,t_max_C eq.(24):")
    t11, t22, t33 = tlxx, tlyy, tlzz
    t12 = self.T_CL[1] * math.sqrt(self.Lxx*self.Lyy)
    t13 = self.T_CL[2] * math.sqrt(self.Lxx*self.Lzz)
    t23 = self.T_CL[5] * math.sqrt(self.Lyy*self.Lzz)
    #
    # Left branch
    #
    if(not (self.is_zero(self.Lxx) or
            self.is_zero(self.Lyy) or
            self.is_zero(self.Lzz))):
      if(t_min_C > t_max_C):
        raise Sorry("Step C (left branch): Empty (tmin_c,tmax_c) interval.")
      t_0 = self.S_L.trace()/3.
      self.show_number(x=t_0, title="t_0 eq.(20):")
      # compose T_lambda and find tau_max (30)
      T_lambda = truncate( m=matrix.sqr(
        [t11, t12, t13,
         t12, t22, t23,
         t13, t23, t33]) )
      self.show_matrix(x=T_lambda, title="T_lambda eq.(29)")
      es = eigensystem.real_symmetric(T_lambda.as_sym_mat3())
      vals = es.values()
      assert vals[0]>=vals[1]>=vals[2]
      tau_max = vals[0]
      #
      if(tau_max < 0):
        raise Sorry("Step C (left branch): Eq.(32): tau_max<0.")
      t_min_tau = max(self.Sxx,self.Syy,self.Szz)-math.sqrt(tau_max)
      t_max_tau = min(self.Sxx,self.Syy,self.Szz)+math.sqrt(tau_max)
      self.show_number(x=[t_min_tau, t_max_tau],
        title="t_min_tau, t_max_tau eq.(31):")
      if(t_min_tau > t_max_tau):
        raise Sorry("Step C (left branch): Empty (tmin_t,tmax_t) interval.")
      # (38):
      arg = t_0**2 + (t11+t22+t33)/3. - (self.Sxx**2+self.Syy**2+self.Szz**2)/3.
      if(arg < 0):
        raise Sorry("Step C (left branch): Negative argument when estimating tmin_a.")
      t_a = math.sqrt(arg)
      self.show_number(x=t_a, title="t_a eq.(38):")
      t_min_a = t_0-t_a
      t_max_a = t_0+t_a
      self.show_number(x=[t_min_a, t_max_a], title="t_min_a, t_max_a eq.(37):")
      # compute t_min, t_max - this is step b)
      t_min = truncate(m=max(t_min_C, t_min_tau, t_min_a), eps_string="%.7f")
      t_max = truncate(m=min(t_max_C, t_max_tau, t_max_a), eps_string="%.7f")
      if(t_min > t_max):
        raise Sorry("Step C (left branch): Intersection of the intervals for t_S is empty.")
      elif(self.is_zero(t_min-t_max)):
        _, b_s, c_s = self.as_bs_cs(t=t_min, txx=t11,tyy=t22,tzz=t33,
          txy=t12,tyz=t23,tzx=t13)
        b_s, c_s = truncate(b_s), truncate(c_s)
        if(b_s >= 0 and c_s <= 0): self.t_S = t_min
        else:
          raise Sorry("Step C (left branch): t_min=t_max gives non positive semidefinite V_lambda.")
      elif(t_min < t_max):
        step = (t_max-t_min)/100000.
        target = 1.e+9
        t_S_best = t_min
        while t_S_best <= t_max:
          _, b_s, c_s = self.as_bs_cs(t=t_S_best, txx=t11,tyy=t22,tzz=t33,
            txy=t12,tyz=t23,tzx=t13)
          if(b_s >= 0 and c_s <= 0):
            target_ = abs(t_0-t_S_best)
            if(target_ < target):
              target = target_
              self.t_S = t_S_best
          t_S_best += step
        assert self.t_S <= t_max
        if(self.t_S is None):
          raise Sorry("Step C (left branch): Interval (t_min,t_max) has no t giving positive semidefinite V.")
        self.t_S = truncate(self.t_S)
    #
    # Right branch, Section 4.4
    #
    else:
      # helper-function to check Cauchy conditions
      def cauchy_conditions(i,j,k, tSs): # diagonals: 0,4,8; 4,0,8; 8,0,4
        if(self.is_zero(self.L_L[i])):
          t_S = self.S_L[i]
          cp1 = (self.S_L[j] - t_S)**2 - self.T_CL[j]*self.L_L[j] # *0.5 is missing!!! ?
          cp2 = (self.S_L[k] - t_S)**2 - self.T_CL[k]*self.L_L[k] # *0.5 is missing!!! ?
          if( not ((cp1 < 0 or self.is_zero(cp1)) and
                   (cp2 < 0 or self.is_zero(cp2))) ):
            raise Sorry("Step C (right branch): Cauchy condition failed (23).")
          a_s, b_s, c_s = self.as_bs_cs(t=t_S, txx=t11,tyy=t22,tzz=t33,
            txy=t12,tyz=t23,tzx=t13)
          self.check_33_34_35(a_s=a_s, b_s=b_s, c_s=c_s)
          tSs.append(t_S)
      #
      tSs = []
      cauchy_conditions(i=0,j=4,k=8, tSs=tSs)
      cauchy_conditions(i=4,j=0,k=8, tSs=tSs)
      cauchy_conditions(i=8,j=0,k=4, tSs=tSs)
      if(len(tSs)==1): self.t_S = tSs[0]
      elif(len(tSs)==0): raise RuntimeError
      else:
        self.t_S = tSs[0]
        for tSs_ in tSs[1:]:
          if(not self.is_zero(x = self.t_S-tSs_)):
            assert 0
    # end-of-procedure check, then truncate
    if(self.t_S is None):
      raise RuntimeError
    self.t_S = truncate(m=self.t_S)
    # override with provided value
    if(self.force_t_S is not None):
      self.t_S = self.force_t_S
    #
    # At this point t_S is found or procedure terminated earlier.
    #
    self.show_number(x=self.t_S, title="t_S:")
    # compute S_C(t_S_), (19)
    self.S_C = self.S_L - matrix.sqr(
      [self.t_S,        0,        0,
              0, self.t_S,        0,
              0,        0, self.t_S])
    self.S_C = truncate(m=matrix.sqr(self.S_C))
    self.show_matrix(x=self.S_C, title="S_C, (26)")
    # find sx, sy, sz
    if(self.is_zero(self.Lxx)):
      if(not self.is_zero(self.S_C[0])):
        raise Sorry("Step C: incompatible L_L and S_C matrices.")
    else:
      self.sx = truncate(m=self.S_C[0]/self.Lxx)
    if(self.is_zero(self.Lyy)):
      if(not self.is_zero(self.S_C[4])):
        raise Sorry("Step C: incompatible L_L and S_C matrices.")
    else:
      self.sy = truncate(m=self.S_C[4]/self.Lyy)
    if(self.is_zero(self.Lzz)):
      if(not self.is_zero(self.S_C[8])):
        raise Sorry("Step C: incompatible L_L and S_C matrices.")
    else:
      self.sz = truncate(m=self.S_C[8]/self.Lzz)
    self.show_number(x=[self.sx,self.sy,self.sz],
      title="Screw parameters (section 4.5): sx,sy,sz:")
    # compose C_L_t_S (26), and V_L (27)
    self.C_L_t_S = matrix.sqr(
      [self.sx*self.S_C[0],                   0,                   0,
                         0, self.sy*self.S_C[4],                   0,
                         0,                   0, self.sz*self.S_C[8]])
    self.C_L_t_S = truncate(m=self.C_L_t_S)
    self.show_matrix(x=self.C_L_t_S, title="C_L(t_S) (26)")
    self.V_L = truncate(m=matrix.sqr(self.T_CL - self.C_L_t_S))
    self.show_matrix(x=self.V_L, title="V_L (26-27)")
    if(not self.is_pd(self.V_L.as_sym_mat3())):
      raise Sorry("Step C: Matrix V[L] is not positive semidefinite.")

  def step_D(self):
    """
    Determination of vibration components (Step D).
    """
    print_step("Step D:", self.log)
    es = self.eigen_system_default_handler(m=self.V_L, suffix="V_L")
    self.v_x, self.v_y, self.v_z = es.x, es.y, es.z
    self.tx, self.ty, self.tz = es.vals[0]**0.5,es.vals[1]**0.5,es.vals[2]**0.5
    self.show_vector(x=self.v_x, title="v_x")
    self.show_vector(x=self.v_y, title="v_y")
    self.show_vector(x=self.v_z, title="v_z")
    if(min(es.vals)<0): raise RuntimeError # checked with Sorry at Step C.
    R = matrix.sqr(
      [self.v_x[0], self.v_y[0], self.v_z[0],
       self.v_x[1], self.v_y[1], self.v_z[1],
       self.v_x[2], self.v_y[2], self.v_z[2]])
    self.V_V = truncate(m=R.transpose()*self.V_L*R)
    self.show_matrix(x=self.V_V, title="V_V")
    self.v_x_M = self.R_ML*self.v_x
    self.v_y_M = self.R_ML*self.v_y
    self.v_z_M = self.R_ML*self.v_z
    self.R_MV = matrix.sqr(
      [self.v_x_M[0], self.v_y_M[0], self.v_z_M[0],
       self.v_x_M[1], self.v_y_M[1], self.v_z_M[1],
       self.v_x_M[2], self.v_y_M[2], self.v_z_M[2]])

  def check_33_34_35(self, a_s, b_s, c_s):
    if(not ((a_s<0 or self.is_zero(a_s)) and
            (b_s>0 or self.is_zero(a_s)) and
            (c_s<0 or self.is_zero(c_s)))):
      raise Sorry("Step C (right branch): Conditions 33-35 failed.")

  def as_bs_cs(self, t, txx,tyy,tzz, txy,tyz,tzx):
    xx = (t-self.Sxx)**2-txx
    yy = (t-self.Syy)**2-tyy
    zz = (t-self.Szz)**2-tzz
    a_s = xx + yy + zz
    b_s = xx*yy + yy*zz + zz*xx - (txy**2+tyz**2+tzx**2)
    c_s = xx*yy*zz - tyz**2*xx - tzx**2*yy - txy**2*zz - 2*txy*tyz*tzx
    return a_s, b_s, c_s

  def is_pd(self, m):
    es = eigensystem.real_symmetric(deepcopy(m))
    r = flex.min(truncate(m=es.values()))
    if(r > 0 or self.is_zero(r)): return True
    else:                         return False

  def is_zero(self, x):
    if(abs(x)<self.eps): return True
    else: return False

  def show_number(self, x, title, prefix="  "):
    print >> self.log, prefix, title
    if(type(x) in [float, int]):
      print >> self.log, "   ", ("%s"%self.ff)%x
    else:
      print >> self.log, "   ", " ".join([str(("%s"%self.ff)%i) for i in x])
    print >> self.log

  def show_matrix(self, x, title, prefix="  "):
    print >> self.log, prefix, title
    f = "%s %s %s"%(self.ff,self.ff,self.ff)
    print >> self.log, prefix, f%(x[0], x[1], x[2])
    print >> self.log, prefix, f%(x[3], x[4], x[5])
    print >> self.log, prefix, f%(x[6], x[7], x[8])
    print >> self.log

  def show_vector(self, x, title, prefix="  "):
    print >> self.log, prefix, title
    print >> self.log, prefix, self.ff%x[0]
    print >> self.log, prefix, self.ff%x[1]
    print >> self.log, prefix, self.ff%x[2]
    print >> self.log

  def eigen_system_default_handler(self, m, suffix):
    es = eigensystem.real_symmetric(m.as_sym_mat3())
    vals, vecs = truncate(m=es.values()), truncate(m=es.vectors())
    print >> self.log, "  eigen values  (%s):"%suffix, " ".join([self.ff%i for i in vals])
    print >> self.log, "  eigen vectors (%s):"%suffix, " ".join([self.ff%i for i in vecs])
    assert vals[0]>=vals[1]>=vals[2]
    ###
    def zero(x, e):
      for i in xrange(len(x)):
        if(abs(x[i])<e): x[i]=0
      return x
    vals = zero(vals, self.eps)
    vecs = zero(vecs, self.eps)
    ###
    # case 1: all different
    if(abs(vals[0]-vals[1])>=self.eps and
       abs(vals[1]-vals[2])>=self.eps and
       abs(vals[0]-vals[2])>=self.eps):
      l_z = matrix.col((vecs[0], vecs[1], vecs[2]))
      l_y = matrix.col((vecs[3], vecs[4], vecs[5]))
      l_x = l_y.cross(l_z)
      vals = [vals[2], vals[1], vals[0]]
    # case 2: all three coincide
    elif((abs(vals[0]-vals[1])<self.eps and
         abs(vals[1]-vals[2])<self.eps and
         abs(vals[0]-vals[2])<self.eps)):
      print >> self.log, "  three eigenvalues are equal: make eigenvectors unit."
      l_x = matrix.col((1, 0, 0))
      l_y = matrix.col((0, 1, 0))
      l_z = matrix.col((0, 0, 1))
    elif([abs(vals[0]-vals[1])<self.eps,
          abs(vals[1]-vals[2])<self.eps,
          abs(vals[0]-vals[2])<self.eps].count(True)==1):
      print >> self.log, "  two eigenvalues are equal."
      #
      l_z = matrix.col((vecs[0], vecs[1], vecs[2]))
      l_y = matrix.col((vecs[3], vecs[4], vecs[5]))
      l_x = l_y.cross(l_z)
      vals = [vals[2], vals[1], vals[0]]
    return group_args(x=l_x, y=l_y, z=l_z, vals=vals)

  def finalize(self):
    return group_args(
      # Libration rms around L-axes
      dx = math.sqrt(self.Lxx),
      dy = math.sqrt(self.Lyy),
      dz = math.sqrt(self.Lzz),
      # Unit vectors defining three Libration axes
      l_x = self.l_x,
      l_y = self.l_y,
      l_z = self.l_z,
      # Rotation axes pass through the points in the L base
      w_L_lx = self.w.w_lx,
      w_L_ly = self.w.w_ly,
      w_L_lz = self.w.w_lz,
      # Rotation axes pass through the points in the M base
      w_M_lx = self.R_ML*self.w.w_lx,
      w_M_ly = self.R_ML*self.w.w_ly,
      w_M_lz = self.R_ML*self.w.w_lz,
      # Correlation shifts sx,sy,sz for libration
      sx = self.sx,
      sy = self.sy,
      sz = self.sz,
      # Vectors defining three Vibration axes
      v_x_M = self.v_x_M,
      v_y_M = self.v_y_M,
      v_z_M = self.v_z_M,
      v_x = self.v_x,
      v_y = self.v_y,
      v_z = self.v_z,
      # Vibration rms along V-axes
      tx = self.tx,
      ty = self.ty,
      tz = self.tz)

  def show_summary(self):
    print_step("SUMMARY:", self.log)
    r = self.result
    self.show_number(x=[r.dx, r.dy, r.dz], title="Libration rms around L-axes")
    #
    self.show_vector(x=r.l_x, title="Unit vector defining libration axis Lx (13):")
    self.show_vector(x=r.l_y, title="Unit vector defining libration axis Ly (13):")
    self.show_vector(x=r.l_z, title="Unit vector defining libration axis Lz (13):")
    #
    self.show_number(x=r.w_L_lx, title="Rotation axis passes through the point in the L basis (15-16):")
    self.show_number(x=r.w_L_ly, title="Rotation axis passes through the point in the L basis (15-16):")
    self.show_number(x=r.w_L_lz, title="Rotation axis passes through the point in the L basis (15-16):")
    #
    self.show_number(x=r.w_M_lx, title="Rotation axis passes through the point in the M basis (40):")
    self.show_number(x=r.w_M_ly, title="Rotation axis passes through the point in the M basis (40):")
    self.show_number(x=r.w_M_lz, title="Rotation axis passes through the point in the M basis (40):")
    #
    self.show_number(x=[r.sx, r.sy, r.sz], title="Correlation shifts sx,sy,sz for libration (4.5):")
    #
    self.show_vector(x=r.v_x, title="Vector defining vibration axis in M basis (41):")
    self.show_vector(x=r.v_y, title="Vector defining vibration axis in M basis (41):")
    self.show_vector(x=r.v_z, title="Vector defining vibration axis in M basis (41):")
    #
    self.show_number(x=[r.tx, r.ty, r.tz], title="Vibration rms along V-axes")

  def self_check(self):
    print_step("Recover T_M, L_M,S_M from base elements:", self.log)
    r = self.result
    R_ML = matrix.sqr(
      [r.l_x[0], r.l_y[0], r.l_z[0],
       r.l_x[1], r.l_y[1], r.l_z[1],
       r.l_x[2], r.l_y[2], r.l_z[2]])
    # T_M
    C_S = matrix.sqr(
      [(r.sx*r.dx)**2,              0,              0,
                    0, (r.sy*r.dy)**2,              0,
                    0,              0, (r.sz*r.dz)**2])
    v_x_M = R_ML*r.v_x
    v_y_M = R_ML*r.v_y
    v_z_M = R_ML*r.v_z
    R_MV = matrix.sqr(
      [v_x_M[0], v_y_M[0], v_z_M[0],
       v_x_M[1], v_y_M[1], v_z_M[1],
       v_x_M[2], v_y_M[2], v_z_M[2]])
    V_V = matrix.sqr(
      [r.tx**2,       0,       0,
             0, r.ty**2,       0,
             0,       0, r.tz**2])
    # D_WL
    wy_lx = r.w_L_lx[1]
    wz_lx = r.w_L_lx[2]
    wx_ly = r.w_L_ly[0]
    wz_ly = r.w_L_ly[2]
    wx_lz = r.w_L_lz[0]
    wy_lz = r.w_L_lz[1]
    d11 = wz_ly**2*r.dy**2 + wy_lz**2*r.dz**2
    d22 = wz_lx**2*r.dx**2 + wx_lz**2*r.dz**2
    d33 = wy_lx**2*r.dx**2 + wx_ly**2*r.dy**2
    d12 = -wx_lz*wy_lz*r.dz**2
    d13 = -wx_ly*wz_ly*r.dy**2
    d23 = -wy_lx*wz_lx*r.dx**2
    D_WL = matrix.sqr(
      [d11, d12, d13,
       d12, d22, d23,
       d13, d23, d33])
    #
    T_M = R_ML * (C_S + D_WL) * R_ML.transpose() + R_MV * V_V * R_MV.transpose()
    self.show_matrix(x=self.T_M, title="Input T_M:")
    self.show_matrix(x=T_M, title="Recoverd T_M:")
    if(flex.max(flex.abs(flex.double(T_M - self.T_M))) > self.self_check_eps):
      raise Sorry("Cannot reconstruct T_M")
    # L_M
    L_L = matrix.sqr(
      [r.dx**2,       0,       0,
             0, r.dy**2,       0,
             0,       0, r.dz**2])
    L_M = R_ML * L_L * R_ML.transpose()
    self.show_matrix(x=self.L_M, title="Input L_M:")
    self.show_matrix(x=L_M, title="Recoverd L_M:")
    if(flex.max(flex.abs(flex.double(L_M - self.L_M))) > self.self_check_eps):
      raise Sorry("Cannot reconstruct L_M")
    # S_M
    S = matrix.sqr(
      [r.sx*r.dx**2,            0,            0,
                  0, r.sy*r.dy**2,            0,
                  0,            0, r.sz*r.dz**2])
    matrix_9 = matrix.sqr(
      [             0,  wz_lx*r.dx**2, -wy_lx*r.dx**2,
       -wz_ly*r.dy**2,              0,  wx_ly*r.dy**2,
        wy_lz*r.dz**2, -wx_lz*r.dz**2,              0])
    S_M = R_ML * (S + matrix_9) * R_ML.transpose()
    self.show_matrix(x=self.S_M, title="Input S_M:")
    self.show_matrix(x=S_M, title="Recoverd S_M:")
    d = truncate(m=matrix.sqr(self.S_M-S_M))
    # remark: diagonal does not have to match, can be different by a constant
    max_diff = flex.max(flex.abs(
      flex.double([d[1],d[2],d[3],d[5],d[6],d[7],d[0]-d[4],d[0]-d[8]])))
    if(max_diff > self.self_check_eps):
      raise Sorry("Cannot reconstruct S_M")
