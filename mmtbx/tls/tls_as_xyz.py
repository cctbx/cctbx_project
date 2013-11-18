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

def generate_TLS_fs_55_56_57(wyi, wzi, wxj, wzj, wxk, wyk, dxs, dys, dzs):
  """
  Helper program - to be removed for final version.
  """
  T = [wzj**2*dys+wyk**2*dzs,
       wzi**2*dxs+wxk**2*dzs,
       wyi**2*dxs+wxj**2*dys,
       -wxk*wyk*dzs,
       -wxj*wzj*dys,
       -wyi*wzi*dxs]
  L = [dxs,dys,dzs,0,0,0]
  S = [       0,  wzi*dxs, -wyi*dxs,
       -wzj*dys,        0,  wxj*dys,
        wyk*dzs, -wxk*dzs,        0]
  return T, L, S

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
  T_p = T + (P*L*P.transpose() + P*S + S.transpose()*P.transpose())
  L_p = L
  S_p = S -L*P#+ L*P.transpose() # same thing
  print "  T':\n", T_p
  print "  L':\n", L_p
  print "  S':\n", S_p
  assert approx_equal(S_p[1], S_p[3])
  assert approx_equal(S_p[2], S_p[6])
  assert approx_equal(S_p[5], S_p[7])
  return group_args(
    T_p = T_p,
    L_p = L_p,
    S_p = S_p,
    p   = p)

def step_b(T_p, L_p, S_p, eps = 1.e-6):
  """
  Diagonalize L'.
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
  print "  l_1:\n", l_1
  print "  l_2:\n", l_2
  print "  l_3:\n", l_3
  l_x, l_y, l_z = l_3, l_2, l_1
  print "  l_x:\n", l_x
  print "  l_y:\n", l_y
  print "  l_z:\n", l_z
  assert approx_equal(l_x.cross(l_y), l_z)
  tmp = l_x.cross(l_y) - l_z
  if(abs(tmp[0])>eps or abs(tmp[1])>eps or abs(tmp[2])>eps):
    print "  inverting l_z..."
    l_z = -1. * l_z
  Rd = matrix.sqr([l_x[0], l_y[0], l_z[0],
                   l_x[1], l_y[1], l_z[1],
                   l_x[2], l_y[2], l_z[2]])
  print "  rotation matrix Rd:\n", Rd
  T_dp = Rd*T_p*Rd.transpose()
  L_dp = Rd*L_p*Rd.transpose()
  S_dp = Rd*S_p*Rd.transpose()
  print "  T'':\n", T_dp
  print "  L'':\n", L_dp
  print "  S'':\n", S_dp
  return group_args(
    T_dp = T_dp,
    L_dp = L_dp,
    S_dp = S_dp,
    Rd   = Rd)

def step_c(T_dp, L_dp, S_dp):
  """
  Position of rotation axes.
  """
  print_step("Step c:")
  L_ = L_dp.as_sym_mat3()
  Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
  assert Lxx != 0 and Lyy != 0 and Lzz != 0
  result = group_args(
    wyi = S_dp[2]/Lxx,
    wzi = S_dp[1]/Lxx,
    wxj = S_dp[5]/Lyy,
    wzj = S_dp[3]/Lyy,
    wxk = S_dp[7]/Lzz,
    wyk = S_dp[6]/Lzz)
  print "  wyi:", result.wyi
  print "  wzi:", result.wzi
  print "  wxj:", result.wxj
  print "  wzj:", result.wzj
  print "  wxk:", result.wxk
  print "  wyk:", result.wyk
  return result

def step_d(S_dp):
  """
  Minimize correlation between translation and rotation.
  """
  print_step("Step d:")
  trace_over_3 = S_dp.trace()/3.
  print "  1/3*traceS'':", trace_over_3
  S_t = S_dp - matrix.sqr(
    (trace_over_3,0,0, 0,trace_over_3,0, 0,0,trace_over_3))
  print "  S_t:\n", S_t
  return S_t

def step_e(S_t, L_dp):
  """
  Find screw components.
  """
  print_step("Step e:")
  L_ = L_dp.as_sym_mat3()
  Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
  assert Lxx != 0 and Lyy != 0 and Lzz != 0
  result = group_args(
    sx_bar = S_t[0]/Lxx,
    sy_bar = S_t[4]/Lyy,
    sz_bar = S_t[8]/Lzz)
  print "  sx_bar:", result.sx_bar
  print "  sy_bar:", result.sy_bar
  print "  sz_bar:", result.sz_bar
  return result

def step_f(c_o, T_dp, L_dp):
  """
  Remove contribution of rotations to T.
  """
  print_step("Step f:")
  L_ = L_dp.as_sym_mat3()
  Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
  wyi = c_o.wyi
  wzi = c_o.wzi
  wxj = c_o.wxj
  wzj = c_o.wzj
  wxk = c_o.wxk
  wyk = c_o.wyk
  delta = [wzj**2*Lxx+wyk**2*Lzz, wzi**2*Lxx+wxk**2*Lzz, wyi**2*Lxx+wxj**2*Lyy,
                    -wxk*wyk*Lzz,          -wxj*wzj*Lyy,          -wyi*wzi*Lxx]
  delta = matrix.sym(sym_mat3=delta)
  print "  delta:\n", delta
  T_tp = T_dp - delta
  print "  T''':\n", T_tp
  return T_tp

def step_g(T_tp, e_o, S_t):
  """
  Remove contribution of screw motion to T.
  """
  print_step("Step g:")
  T_t = T_tp - matrix.sym(sym_mat3=[e_o.sx_bar,e_o.sy_bar,e_o.sz_bar,0,0,0])
  print "  T_t:\n", T_t
  return T_t

def step_h(T_t, b_o, eps=1.e-6):
  """
  Getting the direction of uncorrelated translations.
  """
  print_step("Step h:")
  T_p = b_o.Rd.transpose() * T_t * b_o.Rd
  print "  T_p:\n", T_p
  es = eigensystem.real_symmetric(T_p.as_sym_mat3())
  lam_u,lam_v,lam_w = es.values()
  assert lam_u >= lam_v and lam_v >= lam_w
  vecs = es.vectors()
  t_u = matrix.col((vecs[0], vecs[1], vecs[2]))
  t_v = matrix.col((vecs[3], vecs[4], vecs[5]))
  t_w = matrix.col((vecs[6], vecs[7], vecs[8]))
  t_x = t_w
  t_y = t_v
  t_z = t_u
  tmp = t_x.cross(t_y) - t_z
  if(abs(tmp[0])>eps or abs(tmp[1])>eps or abs(tmp[2])>eps):
    print "  inverting t_z..."
    t_z = -1. * t_z
  print "  lam_u,lam_v,lam_w:", lam_u,lam_v,lam_w
  print "  t_u:\n", t_u
  print "  t_v:\n", t_v
  print "  t_w:\n", t_w
  print "  t_x:\n", t_x
  print "  t_y:\n", t_y
  print "  t_z:\n", t_z
  assert approx_equal(t_u.dot(t_v), 0)
  assert approx_equal(t_v.dot(t_w), 0)
  assert approx_equal(t_u.dot(t_w), 0)
  ####
  Rt = matrix.sqr([t_x[0], t_y[0], t_z[0],
                   t_x[1], t_y[1], t_z[1],
                   t_x[2], t_y[2], t_z[2]])
  print "  Rt:\n", Rt
  T_T = matrix.sym(sym_mat3=[lam_w, lam_v, lam_u, 0,0,0])
  print "  T_T:\n", T_T
  #
  T_T_ = Rt.transpose() * T_p * Rt
  print "  T_T_ = Rt.transpose() * T_p * Rt:\n", T_T_
  #
  return group_args(
    lam_u = lam_u,
    lam_v = lam_v,
    lam_w = lam_w,
    t_u   = t_u,
    t_v   = t_v,
    t_w   = t_w,
    t_x   = t_x,
    t_y   = t_y,
    t_z   = t_z)

#def step_i__get_dxdydz(L_dp):
#  """
#  Generation of shifts from screw rotations.
#  """
#  print_step("step_i__get_dxdydz:")
#  L_ = L_dp.as_sym_mat3()
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
#    (y-c_o.wyi)*(cos(dx0)-1) - (z-c_o.wzi)*sin(dx0),
#    (y-c_o.wyi)*sin(dx0)     + (z-c_o.wzi)*(cos(dx0)-1)
#    ))
#  d_Lj_r_dp = matrix.col((
#    (z-c_o.wzj)*sin(dy0)     + (x-c_o.wxj)*(cos(dy0)-1),
#    syb*dy0,
#    (z-c_o.wzj)*(cos(dy0)-1) - (x-c_o.wxj)*sin(dy0)
#    ))
#  d_Lk_r_dp = matrix.col((
#    (x-c_o.wxk)*(cos(dz0)-1) - (y-c_o.wyk)*sin(dz0),
#    (x-c_o.wxk)*sin(dz0)     + (y-c_o.wyk)*(cos(dz0)-1),
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

def decompose_tls(T, L, S):
  """
  Decompose TLS matrices into components necessary to compute ensemble:
  steps a) - h).
  """
  a_o  = step_a(T=T, L=L, S=S)
  b_o  = step_b(T_p=a_o.T_p, L_p=a_o.L_p, S_p=a_o.S_p)
  c_o  = step_c(T_dp=b_o.T_dp, L_dp=b_o.L_dp, S_dp=b_o.S_dp)
  S_t  = step_d(S_dp=b_o.S_dp)
  e_o  = step_e(S_t=S_t, L_dp=b_o.L_dp)
  T_tp = step_f(c_o=c_o, T_dp = b_o.T_dp, L_dp=b_o.L_dp)
  T_t  = step_g(T_tp=T_tp, e_o=e_o, S_t=S_t)
  h_o  = step_h(T_t = T_t, b_o=b_o)
