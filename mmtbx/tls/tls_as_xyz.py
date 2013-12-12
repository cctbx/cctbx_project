from __future__ import division
from scitbx import matrix
import math
from scitbx.linalg import eigensystem
import random
from libtbx.test_utils import approx_equal
from libtbx import group_args
import mmtbx.utils
from scitbx.array_family import flex
import iotbx.pdb
from libtbx.utils import Sorry
from cctbx import adptbx

random.seed(2679941)

def print_step(s, log):
  n = 50-len(s)
  print >> log, s, "*"*n

def run(pdb_file_name,
        n_models,
        log,
        eps=1.e-6,
        output_file_name="ensemble.pdb"):
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  xrs = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry=pdb_inp.crystal_symmetry_from_cryst1())
  tls_extract = mmtbx.tls.tools.tls_from_pdb_inp(
    remark_3_records = pdb_inp.extract_remark_iii_records(3),
    pdb_hierarchy    = pdb_hierarchy)
  tlso = tls_extract.tls_params
  if(len(tlso)!=1):
    raise Sorry("Only one TLS group per PDB is currently supported.")
  tlso = tlso[0] # XXX one group only
  deg_to_rad_scale = math.pi/180
  # Units: T[A], L[deg**2], S[A*deg]
  T = matrix.sym(sym_mat3=tlso.t)
  L = matrix.sym(sym_mat3=tlso.l)*(deg_to_rad_scale**2)
  S = matrix.sqr(tlso.s)*deg_to_rad_scale
  # sanity check
  if(not adptbx.is_positive_definite(tlso.t, eps)):
    raise Sorry("T matrix is not positive definite.")
  if(not adptbx.is_positive_definite(tlso.l, eps)):
    raise Sorry("L matrix is not positive definite.")
  r = decompose_tls(T=T, L=L, S=S, log=log)
  ensemble_generator(
    decompose_tls_object = r,
    pdb_hierarchy        = pdb_hierarchy,
    xray_structure       = xrs,
    n_models             = n_models,
    log                  = log).write_pdb_file(file_name=output_file_name)

class ensemble_generator(object):
  def __init__(self,
               decompose_tls_object,
               pdb_hierarchy,
               xray_structure,
               n_models,
               log=None):
    if(log is None): log = sys.stdout
    xray_structure.convert_to_isotropic()
    xray_structure = xray_structure.set_b_iso(value=0)
    sites_cart = xray_structure.sites_cart()
    self.states = mmtbx.utils.states(
      xray_structure = xray_structure,
      pdb_hierarchy  = pdb_hierarchy)
    r = decompose_tls_object
    print >> log
    print >> log, "Generating ensemble of %d models:"%n_models
    for trial in xrange(n_models):
      print >> log, "model #%d"%trial
      dx0,dy0,dz0 = step_i__get_dxdydz(
        L_L=r.b_o.L_L, R_PL=r.b_o.R_PL, log = log)
      d_r_M_V  = step_j(h_o=r.h_o, log = log)
      sites_cart_new = flex.vec3_double()
      for site_cart in sites_cart:
        r_L = r.b_o.R_PL.transpose() * site_cart
        d_r_M_L = step_i__compute_delta_L_r_dp(
          r_L=r_L,c_o=r.c_o,e_o=r.e_o,dx0=dx0,dy0=dy0,dz0=dz0, R_PL=r.b_o.R_PL)
        d_r_M = step_k(d_r_M_L=d_r_M_L, d_r_M_V=d_r_M_V)
        sites_cart_new.append(matrix.col(site_cart) + d_r_M)
      self.states.add(sites_cart = sites_cart_new)

  def write_pdb_file(self, file_name):
    self.states.write(file_name = file_name)

def step_i__get_dxdydz(L_L, R_PL, log, eps=1.e-9):
  """
  Generation of shifts from screw rotations.
  """
  print_step("step_i__get_dxdydz:", log)
  L_ = L_L.as_sym_mat3()
  Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
  dx0, dy0, dz0 = 0, 0, 0
  if(Lxx != 0): dx0 = random.normalvariate(0,Lxx)
  if(Lyy != 0): dy0 = random.normalvariate(0,Lyy)
  if(Lzz != 0): dz0 = random.normalvariate(0,Lzz)
  print >> log, "  dx0, dy0, dz0:", dx0, dy0, dz0
  return dx0, dy0, dz0

def step_i__compute_delta_L_r_dp(r_L, c_o, e_o, dx0, dy0, dz0, R_PL):
  sxb,syb,szb = e_o.sx_bar, e_o.sy_bar, e_o.sz_bar
  x,y,z = r_L
  cos, sin = math.cos, math.sin
  d_lx_r_L = matrix.col((
    sxb*dx0,
    (y-c_o.wy_lx)*(cos(dx0)-1) - (z-c_o.wz_lx)*sin(dx0),
    (y-c_o.wy_lx)*sin(dx0)     + (z-c_o.wz_lx)*(cos(dx0)-1)
    ))
  d_ly_r_L = matrix.col((
    (z-c_o.wz_ly)*sin(dy0)     + (x-c_o.wx_ly)*(cos(dy0)-1),
    syb*dy0,
    (z-c_o.wz_ly)*(cos(dy0)-1) - (x-c_o.wx_ly)*sin(dy0)
    ))
  d_lz_r_L = matrix.col((
    (x-c_o.wx_lz)*(cos(dz0)-1) - (y-c_o.wy_lz)*sin(dz0),
    (x-c_o.wx_lz)*sin(dz0)     + (y-c_o.wy_lz)*(cos(dz0)-1),
    szb*dz0
    ))
  d_r_L = d_lx_r_L + d_ly_r_L + d_lz_r_L
  d_r_M = R_PL * d_r_L
  return d_r_M

def step_j(h_o, log):
  """
  Generate shifts from group translation.
  """
  print_step("step_j:", log)
  V_V_ = h_o.V_V.as_sym_mat3()
  tx0, ty0, tz0 = 0, 0, 0
  if(V_V_[0] != 0): tx0 = random.normalvariate(0,V_V_[0])
  if(V_V_[1] != 0): ty0 = random.normalvariate(0,V_V_[1])
  if(V_V_[2] != 0): tz0 = random.normalvariate(0,V_V_[2])
  print >> log, "  u0, v0, w0:", tx0, ty0, tz0
  d_r_V = tx0*h_o.v_x + ty0*h_o.v_y + tz0*h_o.v_z
  d_r_M = h_o.R_MV * d_r_V
  return d_r_M

def step_k(d_r_M_L, d_r_M_V):
  """
  Calculate the total shift in original coordinate system.
  """
  return d_r_M_L + d_r_M_V

class decompose_tls(object):
  def __init__(self, T, L, S, log):
    """
    Decompose TLS matrices into components necessary to compute ensemble:
    steps a) - h).
    """
    self.log = log
    self.T, self.L, self.S = T, L, S
    print_step("Input TLS matrices:", self.log)
    self.show_matrix(x=self.T, prefix="  ", title="T_M")
    self.show_matrix(x=self.L, prefix="  ", title="L_M")
    self.show_matrix(x=self.S, prefix="  ", title="S_M")
    self.a_o  = self.step_a(T  =self.T, L=self.L, S=self.S)
    self.b_o  = self.step_b(T_p=self.a_o.T_p, L_p=self.a_o.L_p, S_p=self.a_o.S_p)
    self.c_o  = self.step_c(T_L=self.b_o.T_L, L_L=self.b_o.L_L, S_L=self.b_o.S_L)
    self.S_C  = self.step_d(S_L=self.b_o.S_L, L_L=self.b_o.L_L, T_L=self.b_o.T_L)
    self.e_o  = self.step_e(S_C=self.S_C, L_L = self.b_o.L_L)
    self.C_LW = self.step_f(c_o=self.c_o, T_L = self.b_o.T_L, L_L=self.b_o.L_L)
    self.g_o  = self.step_g(T_L=self.b_o.T_L, C_LW=self.C_LW, e_o=self.e_o, S_C=self.S_C)
    self.h_o  = self.step_h(V_L=self.g_o.V_L, b_o=self.b_o)

  def show_matrix(self, x, title, prefix="  "):
    print >> self.log, prefix, title
    f = "%11.6f %11.6f %11.6f"
    print >> self.log, prefix, f%(x[0], x[1], x[2])
    print >> self.log, prefix, f%(x[3], x[4], x[5])
    print >> self.log, prefix, f%(x[6], x[7], x[8])
    print >> self.log

  def show_vector(self, x, title, prefix="  "):
    print >> self.log, prefix, title
    print >> self.log, prefix, "%10.6f"%x[0]
    print >> self.log, prefix, "%10.6f"%x[1]
    print >> self.log, prefix, "%10.6f"%x[2]
    print >> self.log

  def eigen_system_default_handler(self, m, eps=1.e-9):
    es = eigensystem.real_symmetric(m.as_sym_mat3())
    vals, vecs = es.values(), es.vectors()
    print >> self.log, "  eigen values :", " ".join(["%11.6f"%i for i in vals])
    print >> self.log, "  eigen vectors:", " ".join(["%11.6f"%i for i in vecs])
    assert vals[0]>=vals[1]>=vals[2]
    # case 1: all different
    if(abs(vals[0]-vals[1])>=eps and
       abs(vals[1]-vals[2])>=eps and
       abs(vals[0]-vals[2])>=eps):
      l_1 = matrix.col((vecs[0], vecs[1], vecs[2]))
      l_2 = matrix.col((vecs[3], vecs[4], vecs[5]))
      l_3 = matrix.col((vecs[6], vecs[7], vecs[8]))
      l_x, l_y, l_z = l_3, l_2, l_1
      vals = [vals[2], vals[1], vals[0]]
      print >> self.log, "re-assign: l_x, l_y, l_z = 3, 2, 1"
      tmp = l_x.cross(l_y) - l_z
      if(abs(tmp[0])>eps or abs(tmp[1])>eps or abs(tmp[2])>eps):
        print >> self.log, "  convert to right-handed"
        l_z = -1. * l_z
    # case 2: all three coincide
    elif(abs(vals[0]-vals[1])<eps and
         abs(vals[1]-vals[2])<eps and
         abs(vals[0]-vals[2])<eps):
      l_x = matrix.col((1, 0, 0))
      l_y = matrix.col((0, 1, 0))
      l_z = matrix.col((0, 0, 1))
    elif([abs(vals[0]-vals[1])<eps,
          abs(vals[1]-vals[2])<eps,
          abs(vals[0]-vals[2])<eps].count(True)==1):
      l_x = matrix.col((1, 0, 0))
      l_y = matrix.col((0, 1, 0))
      l_z = matrix.col((0, 0, 1))
      #print "LOOK"
      #l_1 = matrix.col((vecs[0], vecs[1], vecs[2]))
      #l_2 = matrix.col((vecs[3], vecs[4], vecs[5]))
      #l_3 = matrix.col((vecs[6], vecs[7], vecs[8]))
      #l_x, l_y, l_z = l_3, l_2, l_1
      ##
      #from scitbx.array_family import flex
      #tmp = flex.double([l_z[0], l_z[1], l_z[2]])
      #ma = flex.min(flex.abs(tmp))
      #zz = 999
      #ib = None
      #for i in xrange(3):
      #  zz_ = abs(l_z[i]-ma)
      #  if(zz_<zz):
      #    zz = zz_
      #    ib = i
      #l_y =  matrix.col((l_z[0], l_z[1], l_z[2]))
      #l_y = [l_y[0], l_y[1], l_y[2]]
      #l_y[ib]=0
      #l_y = matrix.col(l_y)
      #l_y = l_y/math.sqrt(l_y[0]**2+l_y[1]**2+l_y[2]**2)
      #l_x = l_y.cross(l_z)
    #
    print >> self.log, "  eigen values :", " ".join(["%11.6f"%i for i in vals])
    print >> self.log, "  eigen vectors:"
    self.show_vector(x=l_x, title="vector x")
    self.show_vector(x=l_y, title="vector y")
    self.show_vector(x=l_z, title="vector z")
    return group_args(x=l_x, y=l_y, z=l_z, vals=vals)

  def step_a(self, T, L, S, eps=1.e-6):
    """
    Shift origin into reaction center. New S' must be symmetric as result of
    origin shift.
    """
    print_step("Step a:", self.log)
    print >> self.log, "  system of equations: m * p = b, p = m_inverse * b"
    L_ = L.as_sym_mat3()
    m = matrix.sqr((
               L_[4],          L_[5], -(L_[0]+L_[1]),
               L_[3], -(L_[0]+L_[2]),         L_[5],
      -(L_[1]+L_[2]),          L_[3],         L_[4]))
    self.show_matrix(x=m, title="m")
    if(abs(m.determinant())<eps):
      print >> self.log, "  det(m)<", eps
      T_p = T
      L_p = L
      S_p = S
      p   = matrix.col((0, 0, 0))
      P = matrix.sqr((
            0, p[2], -p[1],
        -p[2],    0,  p[0],
         p[1], -p[0],    0))
    else:
      b = matrix.col((S[3]-S[1], S[2]-S[6], S[7]-S[5]))
      self.show_vector(x=b, title="b")
      m_inv = m.inverse()
      self.show_matrix(x=m_inv, title="m_inverse")
      p = m_inv * b
      P = matrix.sqr((
            0, p[2], -p[1],
        -p[2],    0,  p[0],
         p[1], -p[0],    0))
      T_p = T #XXX+ (P*L*P.transpose() + P*S + S.transpose()*P.transpose())
      L_p = L #XXX
      S_p = S #XXX- L*P
    #XXX   assert approx_equal(S-L*P, S+L*P.transpose())
    #XXX   # check S_p is symmetric
    #XXX   assert approx_equal(S_p[1], S_p[3])
    #XXX   assert approx_equal(S_p[2], S_p[6])
    #XXX   assert approx_equal(S_p[5], S_p[7])
    #XXX   ###### check system (3)
    #XXX   assert approx_equal(p[1]*L_[3] + p[2]*L_[4] - p[0]*(L_[1]+L_[2]), S[7]-S[5])
    #XXX   assert approx_equal(p[2]*L_[5] + p[0]*L_[3] - p[1]*(L_[2]+L_[0]), S[2]-S[6])
    #XXX   assert approx_equal(p[0]*L_[4] + p[1]*L_[5] - p[2]*(L_[0]+L_[1]), S[3]-S[1])
    self.show_vector(x=p,   title="p")
    self.show_matrix(x=P,   title="P")
    self.show_matrix(x=T_p, title="T_P")
    self.show_matrix(x=L_p, title="L_P")
    self.show_matrix(x=S_p, title="S_P")
    return group_args(
      T_p = T_p,
      L_p = L_p,
      S_p = S_p,
      p   = p)

  def step_b(self, T_p, L_p, S_p, eps = 1.e-6):
    """
    Principal libration axes and transition to L-base.
    """
    print_step("Step b:", self.log)
    es = self.eigen_system_default_handler(m=L_p)
    l_x, l_y, l_z = es.x, es.y, es.z
    self.show_vector(x=l_x, title="l_x")
    self.show_vector(x=l_y, title="l_y")
    self.show_vector(x=l_z, title="l_z")
    R_PL = matrix.sqr(
      [l_x[0], l_y[0], l_z[0],
       l_x[1], l_y[1], l_z[1],
       l_x[2], l_y[2], l_z[2]])
    self.show_matrix(x=R_PL, title="rotation matrix R_PL")
    assert approx_equal(R_PL.transpose(), R_PL.inverse())
    T_L = R_PL.transpose()*T_p*R_PL
    L_L = R_PL.transpose()*L_p*R_PL
    S_L = R_PL.transpose()*S_p*R_PL
    self.show_matrix(x=T_L, title="T_L")
    self.show_matrix(x=L_L, title="L_L")
    self.show_matrix(x=S_L, title="S_L")
    return group_args(
      l_x  = l_x,
      l_y  = l_y,
      l_z  = l_z,
      T_L  = T_L,
      L_L  = L_L,
      S_L  = S_L,
      R_PL = R_PL)

  def step_c(self, T_L, L_L, S_L, eps = 1.e-6):
    """
    Position of rotation axes.
    """
    print_step("Step c:", self.log)
    L_ = L_L.as_sym_mat3()
    Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
    wy_lx=0
    wz_lx=0
    wx_ly=0
    wz_ly=0
    wx_lz=0
    wy_lz=0
    if(Lxx != 0):
      if(abs(S_L[2])>eps): wy_lx =-S_L[2]/Lxx
      if(abs(S_L[1])>eps): wz_lx = S_L[1]/Lxx
    if(Lyy != 0):
      if(abs(S_L[5])>eps): wx_ly = S_L[5]/Lyy
      if(abs(S_L[3])>eps): wz_ly =-S_L[3]/Lyy
    if(Lzz != 0):
      if(abs(S_L[7])>eps): wx_lz =-S_L[7]/Lzz
      if(abs(S_L[6])>eps): wy_lz = S_L[6]/Lzz
    w_lx = matrix.col((0,wy_lx,wz_lx))
    w_ly = matrix.col((wx_ly,0,wz_ly))
    w_lz = matrix.col((wx_lz,wy_lz,0))
    result = group_args(
      wy_lx = wy_lx,
      wz_lx = wz_lx,
      wx_ly = wx_ly,
      wz_ly = wz_ly,
      wx_lz = wx_lz,
      wy_lz = wy_lz,
      w_lx = w_lx,
      w_ly = w_ly,
      w_lz = w_lz)
    self.show_vector(x=w_lx, title="w_lx")
    self.show_vector(x=w_ly, title="w_ly")
    self.show_vector(x=w_lz, title="w_lz")
    return result

  def step_d(self, S_L, L_L, T_L, eps=1.e-6):
    """
    Minimization of correlation between translation and rotation.
    """
    print_step("Step d:", self.log)
    L_ = L_L.as_sym_mat3()
    Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
    T_ = T_L.as_sym_mat3()
    Txx, Tyy, Tzz = T_[0], T_[1], T_[2]
    Sxx, Syy, Szz = S_L[0], S_L[4], S_L[8]
    if(abs(Lxx)<eps): Lxx=0
    if(abs(Lyy)<eps): Lyy=0
    if(abs(Lzz)<eps): Lzz=0
    sqrt = math.sqrt
    TxxLxx = Txx*Lxx
    TyyLyy = Tyy*Lyy
    TzzLzz = Tzz*Lzz
    assert TxxLxx >= 0 and TyyLyy >= 0 and TzzLzz >= 0
    tmin = max(Sxx-sqrt(TxxLxx), Syy-sqrt(TyyLyy), Szz-sqrt(TzzLzz))
    tmax = min(Sxx+sqrt(TxxLxx), Syy+sqrt(TyyLyy), Szz+sqrt(TzzLzz))
    print >> self.log, "  tmin, tmax:", tmin,tmax
    if(abs(tmin)<eps): tmin=0
    if(abs(tmax)<eps): tmax=0
    tS=None
    if(tmin>tmax): assert 0 # must never happen, otherwise S is wrong
    if(abs(tmin-tmax)<eps): tS=0
    if(tmin<tmax):
      trS_L_over3 = S_L.trace()/3.
      if(tmin<=trS_L_over3 and tmax>=trS_L_over3):
        tS = trS_L_over3
      else:
        tmp1 = abs(trS_L_over3-tmin)
        tmp2 = abs(trS_L_over3-tmax)
        if(tmp1<tmp2): tS = tmin
        else:          tS = tmax
    assert tS is not None
    S_C = S_L - matrix.sqr((tS,0,0, 0,tS,0, 0,0,tS))
    print >> self.log, "  tS:", tS
    self.show_matrix(x=S_C, title="S_C")
    return S_C

  def step_e(self, S_C, L_L, eps=1.e-9):
    """
    Find screw components of libration.
    """
    print_step("Step e:", self.log)
    L_ = L_L.as_sym_mat3()
    Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
    sx_bar,sy_bar,sz_bar=0,0,0
    if(abs(Lxx)>eps): sx_bar = S_C[0]/Lxx
    if(abs(Lyy)>eps): sy_bar = S_C[4]/Lyy
    if(abs(Lzz)>eps): sz_bar = S_C[8]/Lzz
    result = group_args(sx_bar = sx_bar, sy_bar = sy_bar, sz_bar = sz_bar)
    print >> self.log, "  sx_bar:", result.sx_bar
    print >> self.log, "  sy_bar:", result.sy_bar
    print >> self.log, "  sz_bar:", result.sz_bar
    return result

  def step_f(self, c_o, T_L, L_L):
    """
    Calculate translational contribution C_LW of rotations due to axes dispalcement
    """
    print_step("Step f:", self.log)
    L_ = L_L.as_sym_mat3()
    Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
    wy_lx = c_o.wy_lx
    wz_lx = c_o.wz_lx
    wx_ly = c_o.wx_ly
    wz_ly = c_o.wz_ly
    wx_lz = c_o.wx_lz
    wy_lz = c_o.wy_lz
    C_LW = [
      wz_ly**2*Lyy+wy_lz**2*Lzz, wz_lx**2*Lxx+wx_lz**2*Lzz, wy_lx**2*Lxx+wx_ly**2*Lyy,
               -wx_lz*wy_lz*Lzz,          -wx_ly*wz_ly*Lyy,          -wy_lx*wz_lx*Lxx]
    C_LW = matrix.sym(sym_mat3=C_LW)
    self.show_matrix(x=C_LW, title="C_LW")
    return C_LW

  def step_g(self, T_L, C_LW, e_o, S_C):
    """
    Calculate translation contribution C_LS of rotation due to screw components
    and resulting V_L matrix.
    """
    print_step("Step g:", self.log)
    C_LS = matrix.sym(sym_mat3=[
      e_o.sx_bar*S_C[0], e_o.sy_bar*S_C[4], e_o.sz_bar*S_C[8], 0,0,0])
    C_L = C_LW + C_LS
    V_L = T_L - C_L
    self.show_matrix(x=C_LS, title="C_LS")
    self.show_matrix(x=C_L , title="C_L ")
    self.show_matrix(x=V_L , title="V_L ")
    return group_args(
      V_L  = V_L,
      C_L  = C_L,
      C_LS = C_LS)

  def step_h(self, V_L, b_o, eps=1.e-6):
    """
    Three uncorrelated translations.
    """
    print_step("Step h:", self.log)
    V_M = b_o.R_PL * V_L * b_o.R_PL.transpose()
    self.show_matrix(x=V_M, title="V_M ")
    es = self.eigen_system_default_handler(m=V_M)
    v_x, v_y, v_z = es.x, es.y, es.z
    lam_u,lam_v,lam_w = es.vals
    self.show_vector(x=v_x, title="v_x")
    self.show_vector(x=v_y, title="v_y")
    self.show_vector(x=v_z, title="v_z")
    assert approx_equal(v_x.dot(v_y), 0)
    assert approx_equal(v_y.dot(v_z), 0)
    assert approx_equal(v_z.dot(v_x), 0)
    R_MV = matrix.sqr([
      v_x[0], v_y[0], v_z[0],
      v_x[1], v_y[1], v_z[1],
      v_x[2], v_y[2], v_z[2]])
    self.show_matrix(x=R_MV, title="R_MV")
    V_V = matrix.sym(sym_mat3=[lam_u, lam_v, lam_w, 0,0,0])
    self.show_matrix(x=V_V, title="V_V")
    assert approx_equal(V_V, R_MV.transpose() * V_M * R_MV) # formula (20)
    return group_args(
      v_x   = v_x,
      v_y   = v_y,
      v_z   = v_z,
      V_M   = V_M,
      V_V   = V_V,
      R_MV  = R_MV)
