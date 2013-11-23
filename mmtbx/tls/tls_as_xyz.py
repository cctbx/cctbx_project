from __future__ import division
from scitbx import matrix
import math
from scitbx.linalg import eigensystem
import random
from libtbx.test_utils import approx_equal
from libtbx import group_args

random.seed(2679941)

def print_step(s):
  """
  Helper program - to be removed for final version.
  """
  n = 50-len(s)
  print s, "*"*n

def step_a(T,L,S):
  """
  Shift origin into reaction center. New S' must be symmetric as result of
  origin shift.
  """
  print_step("Step a:")
  L_ = L.as_sym_mat3()
  m = matrix.sqr((
             L_[4],          L_[5], -(L_[0]+L_[1]),
             L_[3], -(L_[0]+L_[2]),         L_[5],
    -(L_[1]+L_[2]),          L_[3],         L_[4]))
  print "  system of equations: m * p = b, p = m_inverse * b"
  print "  m:\n", m
  b = matrix.col((S[3]-S[1], S[2]-S[6], S[7]-S[5]))
  print "  b:\n", b
  m_inv = m.inverse()
  print "m_inverse:\n", m_inv
  p = m_inv * b
  P = matrix.sqr((
        0, p[2], -p[1],
    -p[2],    0,  p[0],
     p[1], -p[0],    0))
  print "  vector p:\n", p
  print "  matrix P:\n", P
  T_p = T #XXX+ (P*L*P.transpose() + P*S + S.transpose()*P.transpose())
  L_p = L #XXX
  S_p = S #XXX- L*P
#XXX   assert approx_equal(S-L*P, S+L*P.transpose())
  print "  T':\n", T_p
  print "  L':\n", L_p
  print "  S':\n", S_p
#XXX   # check S_p is symmetric
#XXX   assert approx_equal(S_p[1], S_p[3])
#XXX   assert approx_equal(S_p[2], S_p[6])
#XXX   assert approx_equal(S_p[5], S_p[7])
#XXX   ###### check system (3)
#XXX   assert approx_equal(p[1]*L_[3] + p[2]*L_[4] - p[0]*(L_[1]+L_[2]), S[7]-S[5])
#XXX   assert approx_equal(p[2]*L_[5] + p[0]*L_[3] - p[1]*(L_[2]+L_[0]), S[2]-S[6])
#XXX   assert approx_equal(p[0]*L_[4] + p[1]*L_[5] - p[2]*(L_[0]+L_[1]), S[3]-S[1])
  return group_args(
    T_p = T_p,
    L_p = L_p,
    S_p = S_p,
    p   = p)

def step_b(T_p, L_p, S_p, eps = 1.e-6):
  """
  Principal libration axes and transition to L-base.
  """
  print_step("Step b:")
  es = eigensystem.real_symmetric(L_p.as_sym_mat3())
  vals = list(es.values())
  #vals.reverse()
  v = list(es.vectors())
  vecs = v#[v[6],v[7],v[8], v[3],v[4],v[5], v[0],v[1],v[2]]
  print "  lam_1, lam_2, lam_3:\n",vals
  l_1 = matrix.col((vecs[0], vecs[1], vecs[2]))
  l_2 = matrix.col((vecs[3], vecs[4], vecs[5]))
  l_3 = matrix.col((vecs[6], vecs[7], vecs[8]))
  l_x, l_y, l_z = l_3, l_2, l_1
  print "  l_x:\n", l_x
  print "  l_y:\n", l_y
  print "  l_z:\n", l_z
  assert approx_equal(l_x.cross(l_y), l_z)
  tmp = l_x.cross(l_y) - l_z
  if(abs(tmp[0])>eps or abs(tmp[1])>eps or abs(tmp[2])>eps):
    print "  inverting l_z..."
    l_z = -1. * l_z
  R_PL = matrix.sqr(
    [l_x[0], l_y[0], l_z[0],
     l_x[1], l_y[1], l_z[1],
     l_x[2], l_y[2], l_z[2]])
  print "  rotation matrix R_PL:\n", R_PL
  assert approx_equal(R_PL.transpose(), R_PL.inverse())
  T_L = R_PL.transpose()*T_p*R_PL
  L_L = R_PL.transpose()*L_p*R_PL
  S_L = R_PL.transpose()*S_p*R_PL
  print "  T_L:\n", T_L
  print "  L_L:\n", L_L
  print "  S_L:\n", S_L
  return group_args(
    l_x  = l_x,
    l_y  = l_y,
    l_z  = l_z,
    T_L  = T_L,
    L_L  = L_L,
    S_L  = S_L,
    R_PL = R_PL)

def step_c(T_L, L_L, S_L):
  """
  Position of rotation axes.
  """
  print_step("Step c:")
  L_ = L_L.as_sym_mat3()
  Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
  assert Lxx != 0 and Lyy != 0 and Lzz != 0
  wy_lx =-S_L[2]/Lxx
  wz_lx = S_L[1]/Lxx
  wx_ly = S_L[5]/Lyy
  wz_ly =-S_L[3]/Lyy
  wx_lz =-S_L[7]/Lzz
  wy_lz = S_L[6]/Lzz
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
  print "  w_lx:\n", w_lx
  print "  w_ly:\n", w_ly
  print "  w_lz:\n", w_lz
  return result

def step_d(S_L):
  """
  Minimization of correlation between translation and rotation.
  """
  print_step("Step d:")
  trace_over_3 = S_L.trace()/3.
  print "  1/3*traceS'':", trace_over_3
  S_C = S_L - matrix.sqr(
    (trace_over_3,0,0, 0,trace_over_3,0, 0,0,trace_over_3))
  print "  S_C:\n", S_C
  return S_C

def step_e(S_C, L_L):
  """
  Find screw components of libration.
  """
  print_step("Step e:")
  L_ = L_L.as_sym_mat3()
  Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
  assert Lxx != 0 and Lyy != 0 and Lzz != 0
  result = group_args(
    sx_bar = S_C[0]/Lxx,
    sy_bar = S_C[4]/Lyy,
    sz_bar = S_C[8]/Lzz)
  print "  sx_bar:", result.sx_bar
  print "  sy_bar:", result.sy_bar
  print "  sz_bar:", result.sz_bar
  return result

def step_f(c_o, T_L, L_L):
  """
  Calculate translational contribution C_LW of rotations due to axes dispalcement
  """
  print_step("Step f:")
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
  print "  C_LW:\n", C_LW
  return C_LW

def step_g(T_L, C_LW, e_o, S_C):
  """
  Calculate translation contribution C_LS of rotation due to screw components
  and resulting V_L matrix.
  """
  print_step("Step g:")
  C_LS = matrix.sym(sym_mat3=[
    e_o.sx_bar*S_C[0], e_o.sy_bar*S_C[4], e_o.sz_bar*S_C[8], 0,0,0])
  C_L = C_LW + C_LS
  V_L = T_L - C_L
  print "  C_LS:\n", C_LS
  print "  C_L :\n", C_L
  print "  V_L :\n", V_L
  return group_args(
    V_L  = V_L,
    C_L  = C_L,
    C_LS = C_LS)

def step_h(V_L, b_o, eps=1.e-6):
  """
  Three uncorrelated translations.
  """
  print_step("Step h:")
  V_M = b_o.R_PL * V_L * b_o.R_PL.transpose()
  print "  V_M:\n", V_M
  es = eigensystem.real_symmetric(V_M.as_sym_mat3())
  lam_u,lam_v,lam_w = es.values()
  assert lam_u >= lam_v and lam_v >= lam_w
  vecs = es.vectors()
  t_u = matrix.col((vecs[0], vecs[1], vecs[2]))
  t_v = matrix.col((vecs[3], vecs[4], vecs[5]))
  t_w = matrix.col((vecs[6], vecs[7], vecs[8]))
  v_x, v_y, v_z = t_w, t_v, t_u
  tmp = v_x.cross(v_y) - v_z
  if(abs(tmp[0])>eps or abs(tmp[1])>eps or abs(tmp[2])>eps):
    print "  inverting t_z..."
    v_z = -1. * v_z
  print "  v_x:\n", v_x
  print "  v_y:\n", v_y
  print "  v_z:\n", v_z
  assert approx_equal(t_u.dot(t_v), 0)
  assert approx_equal(t_v.dot(t_w), 0)
  assert approx_equal(t_u.dot(t_w), 0)
  R_MV = matrix.sqr([
    v_x[0], v_y[0], v_z[0],
    v_x[1], v_y[1], v_z[1],
    v_x[2], v_y[2], v_z[2]])
  print "  R_MV:\n", R_MV
  V_V = matrix.sym(sym_mat3=[lam_w, lam_v, lam_u, 0,0,0])
  print "  V_V:\n", V_V
  #
  assert approx_equal(V_V, R_MV.transpose() * V_M * R_MV) # formula (20)
  #
  return group_args(
    lam_u = lam_u,
    lam_v = lam_v,
    lam_w = lam_w,
    v_x   = v_x,
    v_y   = v_y,
    v_z   = v_z,
    V_M   = V_M,
    V_V   = V_V)

#def step_i__get_dxdydz(L_L):
#  """
#  Generation of shifts from screw rotations.
#  """
#  print_step("step_i__get_dxdydz:")
#  L_ = L_L.as_sym_mat3()
#  Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
#  #max_ampl = (8 * math.pi / 180)**2
#  dx0 = random.normalvariate(0,Lxx)
#  dy0 = random.normalvariate(0,Lyy)
#  dz0 = random.normalvariate(0,Lzz)
#  print "  dx0, dy0, dz0:", dx0, dy0, dz0
#  return dx0, dy0, dz0
#
#def step_i__compute_delta_L_r_dp(r_dp, c_o, e_o, dx0, dy0, dz0):
#  sxb,syb,szb = e_o.sx_bar, e_o.sy_bar, e_o.sz_bar
#  x,y,z = r_dp
#  cos, sin = math.cos, math.sin
#  d_Li_r_dp = matrix.col((
#    sxb*dx0,
#    (y-c_o.wy_lx)*(cos(dx0)-1) - (z-c_o.wz_lx)*sin(dx0),
#    (y-c_o.wy_lx)*sin(dx0)     + (z-c_o.wz_lx)*(cos(dx0)-1)
#    ))
#  d_Lj_r_dp = matrix.col((
#    (z-c_o.wz_ly)*sin(dy0)     + (x-c_o.wx_ly)*(cos(dy0)-1),
#    syb*dy0,
#    (z-c_o.wz_ly)*(cos(dy0)-1) - (x-c_o.wx_ly)*sin(dy0)
#    ))
#  d_Lk_r_dp = matrix.col((
#    (x-c_o.wx_lz)*(cos(dz0)-1) - (y-c_o.wy_lz)*sin(dz0),
#    (x-c_o.wx_lz)*sin(dz0)     + (y-c_o.wy_lz)*(cos(dz0)-1),
#    szb*dz0
#    ))
#  d_L_rdp = d_Li_r_dp + d_Lj_r_dp + d_Lk_r_dp
#  #print d_L_rdp
#  return d_L_rdp
#
#def step_j(h_o):
#  """
#  Generate shifts from group translation.
#  """
#  print_step("step_j:")
#  max_ampl = 8 * math.pi / 180
#  u0 = random.normalvariate(0,h_o.lam_u)
#  v0 = random.normalvariate(0,h_o.lam_v)
#  w0 = random.normalvariate(0,h_o.lam_w)
#  print "  u0, v0, w0:", u0, v0, w0
#  d_T_r_dp = u0*h_o.t_u + v0*h_o.t_v + w0*h_o.t_w
#  print "  d_T_r_dp:", d_T_r_dp
#  return d_T_r_dp
#
#def step_k(d_T_r_dp, d_L_rdp, b_o):
#  """
#  Calculate shift in original coordinate system.
#  """
#  d_TL_r_dp = b_o.R * (d_T_r_dp + d_L_rdp)
#  return d_TL_r_dp
#
#def compute_r_dp(sites_cart, p, R):
#  """
#  Compute r".
#  """
#  print_step("Compute r'':")
#  r_p = sites_cart - p
#  r_dp = R.transpose().elems*r_p
#  return r_dp

class decompose_tls(object):
  def __init__(self, T, L, S):
    """
    Decompose TLS matrices into components necessary to compute ensemble:
    steps a) - h).
    """
    self.T, self.L, self.S = T, L, S
    self.a_o  = step_a(T  =self.T, L=self.L, S=self.S)
    self.b_o  = step_b(T_p=self.a_o.T_p, L_p=self.a_o.L_p, S_p=self.a_o.S_p)
    self.c_o  = step_c(T_L=self.b_o.T_L, L_L=self.b_o.L_L, S_L=self.b_o.S_L)
    self.S_C  = step_d(S_L=self.b_o.S_L)
    self.e_o  = step_e(S_C=self.S_C, L_L = self.b_o.L_L)
    self.C_LW = step_f(c_o=self.c_o, T_L = self.b_o.T_L, L_L=self.b_o.L_L)
    self.g_o  = step_g(T_L=self.b_o.T_L, C_LW=self.C_LW, e_o=self.e_o, S_C=self.S_C)
    self.h_o  = step_h(V_L=self.g_o.V_L, b_o=self.b_o)
