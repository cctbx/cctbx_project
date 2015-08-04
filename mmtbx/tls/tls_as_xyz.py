from __future__ import division
from scitbx import matrix
import math
import random
import mmtbx.utils
from scitbx.array_family import flex
import iotbx.pdb
from libtbx.utils import Sorry
from cctbx import adptbx
import sys
from mmtbx.tls import analysis

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
  r = analysis.run(T=T, L=L, S=S, log=log)
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
        L_L=r.L_L, R_PL=r.R_ML, log = log)
      d_r_M_V  = step_j(r=r, log = log)
      sites_cart_new = flex.vec3_double()
      for site_cart in sites_cart:
        r_L = r.R_ML.transpose() * site_cart
        d_r_M_L = step_i__compute_delta_L_r_dp(r=r,
          r_L=r_L, dx0=dx0,dy0=dy0,dz0=dz0)
        d_r_M = step_k(d_r_M_L=d_r_M_L, d_r_M_V=d_r_M_V)
        sites_cart_new.append(matrix.col(site_cart) + d_r_M)
      self.states.add(sites_cart = sites_cart_new)

  def write_pdb_file(self, file_name):
    self.states.write(file_name = file_name)

def step_i__get_dxdydz(L_L, R_PL, log, eps=1.e-7):
  """
  Generation of shifts from screw rotations.
  """
  L_ = L_L.as_sym_mat3()
  Lxx, Lyy, Lzz = L_[0], L_[1], L_[2]
  dx0, dy0, dz0 = 0, 0, 0
  if(abs(Lxx)>eps): dx0 = random.normalvariate(0,(2*Lxx)**0.5) # ???
  if(abs(Lyy)>eps): dy0 = random.normalvariate(0,(2*Lyy)**0.5) # ???
  if(abs(Lzz)>eps): dz0 = random.normalvariate(0,(2*Lzz)**0.5) # ???
  return dx0, dy0, dz0

def step_i__compute_delta_L_r_dp(r, r_L, dx0, dy0, dz0):
  sxb,syb,szb = r.sx, r.sy, r.sz
  x,y,z = r_L
  cos, sin = math.cos, math.sin
  d_lx_r_L = matrix.col((
    sxb*dx0,
    (y-r.w.wy_lx)*(cos(dx0)-1) - (z-r.w.wz_lx)*sin(dx0),
    (y-r.w.wy_lx)*sin(dx0)     + (z-r.w.wz_lx)*(cos(dx0)-1)
    ))
  d_ly_r_L = matrix.col((
    (z-r.w.wz_ly)*sin(dy0)     + (x-r.w.wx_ly)*(cos(dy0)-1),
    syb*dy0,
    (z-r.w.wz_ly)*(cos(dy0)-1) - (x-r.w.wx_ly)*sin(dy0)
    ))
  d_lz_r_L = matrix.col((
    (x-r.w.wx_lz)*(cos(dz0)-1) - (y-r.w.wy_lz)*sin(dz0),
    (x-r.w.wx_lz)*sin(dz0)     + (y-r.w.wy_lz)*(cos(dz0)-1),
    szb*dz0
    ))
  d_r_L = d_lx_r_L + d_ly_r_L + d_lz_r_L
  d_r_M = r.R_ML * d_r_L
  return d_r_M

def step_j(r, log):
  """
  Generate shifts from group translation.
  """
#XXX  print_step("step_j:", log)
  V_V_ = r.V_V.as_sym_mat3()
  tx0, ty0, tz0 = 0, 0, 0
  if(V_V_[0] != 0): tx0 = random.gauss(0,(2*V_V_[0]**0.5)) # ???
  if(V_V_[1] != 0): ty0 = random.gauss(0,(2*V_V_[1]**0.5)) # ???
  if(V_V_[2] != 0): tz0 = random.gauss(0,(2*V_V_[2]**0.5)) # ???
  #print >> log, "  u0, v0, w0:, LOOK", tx0, ty0, tz0
  d_r_V = tx0*r.v_x + ty0*r.v_y + tz0*r.v_z
  d_r_M = r.R_MV * d_r_V
  return d_r_M

def step_k(d_r_M_L, d_r_M_V):
  """
  Calculate the total shift in original coordinate system.
  """
  return d_r_M_L + d_r_M_V # (42)
