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
from mmtbx.tls.decompose import decompose_tls_matrices
from libtbx.utils import multi_out
from libtbx import adopt_init_args

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
    # Reduce TLS matrices to fundamental motions
    try:
      r = run(T=T, L=L, S=S, log=log)
    except Exception, e:
      print >> log, str(e)
    log.close()

def show_matrix(x, title, prefix="  ", log=None):
  if(log is None): log = sys.stdout
  print >> log, prefix, title
  ff = "%12.9f"
  f = "%s %s %s"%(ff,ff,ff)
  print >> log, prefix, f%(x[0], x[1], x[2])
  print >> log, prefix, f%(x[3], x[4], x[5])
  print >> log, prefix, f%(x[6], x[7], x[8])
  print >> log

def show_vector(x, title, prefix="  ", log=None):
  if(log is None): log = sys.stdout
  ff = "%12.9f"
  print >> log, prefix, title
  print >> log, prefix, ff%x[0]
  print >> log, prefix, ff%x[1]
  print >> log, prefix, ff%x[2]
  print >> log

def show_number(x, title, prefix="  ", log=None):
  if(log is None): log = sys.stdout
  ff = "%12.9f"
  print >> log, prefix, title
  if(type(x) in [float, int]):
    print >> log, "   ", ("%s"%ff)%x
  else:
    print >> log, "   ", " ".join([str(("%s"%ff)%i) for i in x])
  print >> log

# remark: this class is a legacy of the original python implementation
# remark: and is now a wrapper around the c++ code. It is intended to
# remark: allow backwards compatibility with other code...
class run(object):
  def __init__(self, T, L, S, log=sys.stdout, eps=1.e-6, self_check_eps=1.e-5,
               force_t_S=None, find_t_S_using_formula="11"):
    """
    Decompose TLS matrices into physically iterpretable components.
    Input L & S matrices must be in radians.
    """
    self.log = log
    self.ff = "%12.9f"
    self.eps = eps
    self.self_check_eps = self_check_eps

    # Check that choices for finding t_S are compatible
    assert find_t_S_using_formula in ["10","11",None]
    assert [find_t_S_using_formula, force_t_S].count(None) == 1, \
            'cannot provide both find_t_S_using_formula and force_t_S'
    # Possible values of find_t_S_using_formula given force_t_S
    if force_t_S is None:
      assert find_t_S_using_formula is not None
      # Need to make double to pass to c++
      force_t_S = -1.0
    else:
      assert find_t_S_using_formula is None
      # Needs to be an string to pass to c++
      find_t_S_using_formula = "Force"

    # Choose how to deal with t_S END
    print >> self.log, "Small is defined as:", self.eps
    self.T_M, self.L_M, self.S_M = T, L, S
    print_step("Input TLS matrices:", self.log)
    show_matrix(x=self.T_M, prefix="  ", title="T_M", log=self.log)
    show_matrix(x=self.L_M, prefix="  ", title="L_M", log=self.log)
    show_matrix(x=self.S_M, prefix="  ", title="S_M", log=self.log)

    # apply decomposition routine
    self.decomposition = decompose_tls_matrices(T.as_sym_mat3(),
                                                L.as_sym_mat3(),
                                                S.as_mat3(),
                                                l_and_s_in_degrees=False,
                                                verbose=False,
                                                tol=eps,
                                                t_S_formula=find_t_S_using_formula,
                                                t_S_value=force_t_S)
    # Raise error if decomposition fails to mirror previous implementation
    if not self.decomposition.is_valid():
        raise Sorry(self.decomposition.error())

    # compose result-object
    self.result = self.finalize()
    self.show_summary()
    # consistency check: restore input T_M, L_M, S_M from base elements
    self.self_check()

  def finalize(self):
    """Return a summary object of the decomposition"""
    # unpack
    d = self.decomposition
    tx, ty, tz = d.v_amplitudes
    dx, dy, dz = d.l_amplitudes
    sx, sy, sz = d.s_amplitudes
    l_x, l_y, l_z = d.l_axis_directions
    w_x, w_y, w_z = d.l_axis_intersections
    v_x, v_y, v_z = d.v_axis_directions
    # Extract rotation matrix (tls_motions need input in the L-basis)
    R_ML = matrix.sqr(d.R_ML)
    return group_args(
      # Libration rms around L-axes
      dx = dx, dy = dy, dz = dz,
      # Unit vectors defining three Libration axes
      l_x = l_x, l_y = l_y, l_z = l_z,
      # Rotation axes pass through the points in the L base
      w_L_lx = R_ML*l_x,
      w_L_ly = R_ML*l_y,
      w_L_lz = R_ML*l_z,
      # Rotation axes pass through the points in the M base
      w_M_lx = w_x,
      w_M_ly = w_y,
      w_M_lz = w_z,
      # Correlation shifts sx,sy,sz for libration
      sx = sx, sy = sy, sz = sz,
      # Vectors defining three Vibration axes
      v_x_M = v_x,
      v_y_M = v_y,
      v_z_M = v_z,
      # Recalculate these in the L basis
      v_x = R_ML*v_x,
      v_y = R_ML*v_y,
      v_z = R_ML*v_z,
      # Vibration rms along V-axes
      tx = tx, ty = ty, tz = tz)

  def show_summary(self):
    print_step("SUMMARY:", self.log)
    r = self.result
    show_number(x=[r.dx, r.dy, r.dz], title="Libration rms around L-axes", log=self.log)
    #
    show_vector(x=r.l_x, title="Unit vector defining libration axis Lx (13):", log=self.log)
    show_vector(x=r.l_y, title="Unit vector defining libration axis Ly (13):", log=self.log)
    show_vector(x=r.l_z, title="Unit vector defining libration axis Lz (13):", log=self.log)
    #
    show_number(x=r.w_L_lx, title="Rotation axis passes through the point in the L basis (15-16):", log=self.log)
    show_number(x=r.w_L_ly, title="Rotation axis passes through the point in the L basis (15-16):", log=self.log)
    show_number(x=r.w_L_lz, title="Rotation axis passes through the point in the L basis (15-16):", log=self.log)
    #
    show_number(x=r.w_M_lx, title="Rotation axis passes through the point in the M basis (40):", log=self.log)
    show_number(x=r.w_M_ly, title="Rotation axis passes through the point in the M basis (40):", log=self.log)
    show_number(x=r.w_M_lz, title="Rotation axis passes through the point in the M basis (40):", log=self.log)
    #
    show_number(x=[r.sx, r.sy, r.sz], title="Correlation shifts sx,sy,sz for libration (8):", log=self.log)
    #
    show_vector(x=r.v_x_M, title="Vector defining vibration axis in M basis (41):", log=self.log)
    show_vector(x=r.v_y_M, title="Vector defining vibration axis in M basis (41):", log=self.log)
    show_vector(x=r.v_z_M, title="Vector defining vibration axis in M basis (41):", log=self.log)
    #
    show_number(x=[r.tx, r.ty, r.tz], title="Vibration rms along V-axes", log=self.log)
    return self

  def self_check(self, show=True):
    if(show):
      print_step("Recover T_M, L_M,S_M from base elements:", self.log)
    # i = "input"
    i = self
    # Extract result from decomposition
    r = self.result
    # Intentionally overwrite variable: r = "result"
    r = tls_from_motions(
      dx=r.dx, dy=r.dy, dz=r.dz,
      l_x=r.l_x, l_y=r.l_y, l_z=r.l_z,
      sx=r.sx, sy=r.sy, sz=r.sz,
      tx=r.tx, ty=r.ty, tz=r.tz,
      v_x=r.v_x, v_y=r.v_y, v_z=r.v_z,
      w_M_lx=r.w_M_lx,
      w_M_ly=r.w_M_ly,
      w_M_lz=r.w_M_lz)
    #
    if(show):
      show_matrix(x=i.T_M, title="Input T_M:", log=self.log)
      show_matrix(x=r.T_M, title="Recovered T_M:", log=self.log)
    if(flex.max(flex.abs(flex.double(r.T_M - i.T_M))) > self.self_check_eps):
      raise Sorry("Cannot reconstruct T_M")
    # L_M
    if(show):
      show_matrix(x=i.L_M, title="Input L_M:", log=self.log)
      show_matrix(x=r.L_M, title="Recovered L_M:", log=self.log)
    if(flex.max(flex.abs(flex.double(r.L_M - i.L_M))) > self.self_check_eps):
      raise Sorry("Cannot reconstruct L_M")
    # S_M
    if(show):
      show_matrix(x=i.S_M, title="Input S_M:", log=self.log)
      show_matrix(x=r.S_M, title="Recovered S_M:", log=self.log)
    d = matrix.sqr(r.S_M - i.S_M)
    # remark: diagonal does not have to match, can be different by a constant
    max_diff = flex.max(flex.abs(
      flex.double([d[1],d[2],d[3],d[5],d[6],d[7],d[0]-d[4],d[0]-d[8]])))
    if(max_diff > self.self_check_eps):
      raise Sorry("Cannot reconstruct S_M")
    return r

class tls_from_motions(object):
  def __init__(self,
               dx,dy,dz,
               l_x,l_y,l_z,
               sx,sy,sz,
               tx,ty,tz,
               v_x,v_y,v_z,
               w_M_lx,
               w_M_ly,
               w_M_lz):
    adopt_init_args(self, locals())
    self.R_ML = matrix.sqr(
      [self.l_x[0], self.l_y[0], self.l_z[0],
       self.l_x[1], self.l_y[1], self.l_z[1],
       self.l_x[2], self.l_y[2], self.l_z[2]])
    # T_M
    C_S = matrix.sqr(
      [(self.sx*self.dx)**2,                    0,                    0,
                          0, (self.sy*self.dy)**2,                    0,
                          0,                    0, (self.sz*self.dz)**2])
    self.v_x_M = self.R_ML*self.v_x
    self.v_y_M = self.R_ML*self.v_y
    self.v_z_M = self.R_ML*self.v_z
    self.R_MV = matrix.sqr(
      [self.v_x_M[0], self.v_y_M[0], self.v_z_M[0],
       self.v_x_M[1], self.v_y_M[1], self.v_z_M[1],
       self.v_x_M[2], self.v_y_M[2], self.v_z_M[2]])
    V_V = matrix.sqr(
      [self.tx**2,       0,       0,
             0, self.ty**2,       0,
             0,       0, self.tz**2])
    # D_WL
    w_L_lx = self.R_ML.transpose()*w_M_lx
    w_L_ly = self.R_ML.transpose()*w_M_ly
    w_L_lz = self.R_ML.transpose()*w_M_lz
    self.wy_lx = w_L_lx[1]
    self.wz_lx = w_L_lx[2]
    self.wx_ly = w_L_ly[0]
    self.wz_ly = w_L_ly[2]
    self.wx_lz = w_L_lz[0]
    self.wy_lz = w_L_lz[1]
    d11 =  self.wz_ly**2*self.dy**2 + self.wy_lz**2*self.dz**2
    d22 =  self.wz_lx**2*self.dx**2 + self.wx_lz**2*self.dz**2
    d33 =  self.wy_lx**2*self.dx**2 + self.wx_ly**2*self.dy**2
    d12 = -self.wx_lz*self.wy_lz*self.dz**2
    d13 = -self.wx_ly*self.wz_ly*self.dy**2
    d23 = -self.wy_lx*self.wz_lx*self.dx**2
    D_WL = matrix.sqr(
      [d11, d12, d13,
       d12, d22, d23,
       d13, d23, d33])
    # T_M
    self.T_M = self.R_ML*(C_S+D_WL)*self.R_ML.transpose() + \
      self.R_MV * V_V * self.R_MV.transpose()
    # L_M
    L_L = matrix.sqr(
      [self.dx**2,       0,       0,
             0, self.dy**2,       0,
             0,       0, self.dz**2])
    self.L_M = self.R_ML * L_L * self.R_ML.transpose()
    # S_M
    S = matrix.sqr(
      [self.sx*self.dx**2,            0,            0,
                  0, self.sy*self.dy**2,            0,
                  0,            0, self.sz*self.dz**2])
    matrix_9 = matrix.sqr(
      [                     0,  self.wz_lx*self.dx**2, -self.wy_lx*self.dx**2,
       -self.wz_ly*self.dy**2,                      0,  self.wx_ly*self.dy**2,
        self.wy_lz*self.dz**2, -self.wx_lz*self.dz**2,                      0])
    self.S_M = self.R_ML * (S + matrix_9) * self.R_ML.transpose()

  def show(self, log=None):
    if(log is None): log = sys.stdout
    show_number(x=[self.tx, self.ty, self.tz], title="tx, ty, tz (4):", log=log)
    show_number(x=[self.dx, self.dy, self.dz], title="dx, dy, dz (5):", log=log)
    show_number(x=[self.sx, self.sy, self.sz], title="sx, sy, sz (8):", log=log)
    show_vector(x=self.l_x, title="lx (13):", log=log)
    show_vector(x=self.l_y, title="ly (13):", log=log)
    show_vector(x=self.l_z, title="lz (13):", log=log)
    show_vector(x=self.v_x, title="vx (41):", log=log)
    show_vector(x=self.v_y, title="vy (41):", log=log)
    show_vector(x=self.v_z, title="vz (41):", log=log)
    show_number(x=self.w_M_lx, title="w_M_lx (40):", log=log)
    show_number(x=self.w_M_ly, title="w_M_ly (40):", log=log)
    show_number(x=self.w_M_lz, title="w_M_lz (40):", log=log)
    show_matrix(x=self.T_M, title="T_M (computed from inputs):", log=log)
    show_matrix(x=self.L_M, title="L_M (computed from inputs):", log=log)
    show_matrix(x=self.S_M, title="S_M (computed from inputs):", log=log)
