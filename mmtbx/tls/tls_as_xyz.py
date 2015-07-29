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
import sys
from copy import deepcopy

random.seed(2679941)

def print_step(s, log):
  n = 50-len(s)
  print >> log, s, "*"*n

def run(pdb_file_name,
        n_models,
        log,
        eps=1.e-7,
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

def step_i__get_dxdydz(L_L, R_PL, log, eps=1.e-7):
  """
  Generation of shifts from screw rotations.
  """
#XXX  print_step("step_i__get_dxdydz:", log)
  L_ = L_L.as_sym_mat3()
  Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
  dx0, dy0, dz0 = 0, 0, 0
  if(abs(Lxx)>eps): dx0 = random.normalvariate(0,Lxx)
  if(abs(Lyy)>eps): dy0 = random.normalvariate(0,Lyy)
  if(abs(Lzz)>eps): dz0 = random.normalvariate(0,Lzz)
#XXX  print >> log, "  dx0, dy0, dz0:", dx0, dy0, dz0
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
#XXX  print_step("step_j:", log)
  V_V_ = h_o.V_V.as_sym_mat3()
  tx0, ty0, tz0 = 0, 0, 0
  if(V_V_[0] != 0): tx0 = random.normalvariate(0,V_V_[0])
  if(V_V_[1] != 0): ty0 = random.normalvariate(0,V_V_[1])
  if(V_V_[2] != 0): tz0 = random.normalvariate(0,V_V_[2])
#XXX  print >> log, "  u0, v0, w0:", tx0, ty0, tz0
  d_r_V = tx0*h_o.v_x + ty0*h_o.v_y + tz0*h_o.v_z
  d_r_M = h_o.R_MV * d_r_V
  return d_r_M

def step_k(d_r_M_L, d_r_M_V):
  """
  Calculate the total shift in original coordinate system.
  """
  return d_r_M_L + d_r_M_V

def truncate(m, eps_string="%.6f"):
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

