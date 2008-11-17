from scitbx.rigid_body_dynamics import featherstone
from scitbx import matrix
import math

class six_dof_euler_params(object):

  def __init__(O, qE, qr):
    O.qE = qE
    O.qr = qr
    #
    O.E = RBDA_Eq_4_12(qE)
    O.r = O.E.transpose() * qr # RBDA Tab. 4.1
    #
    O.T = matrix.rt((O.E, -O.E * O.r)) # RBDA Eq. 2.28
    O.T_inv = matrix.rt((O.E.transpose(), O.r))
    #
    O.Xj = featherstone.Xrot(O.E) \
         * featherstone.Xtrans(O.r) # RBDA Tab. 4.1 footnote
    O.S = None
    O.S_ring = None

  def Xj_S_S_ring(O, q, qd):
    return O.Xj, O.S, O.S_ring

  def time_step_position(O, v_spatial, delta_t):
    w_body_frame, v_body_frame = matrix.col_list([
      v_spatial.elems[:3], v_spatial.elems[3:]])
    qEd = RBDA_Eq_4_13(O.qE.elems) * (O.E * w_body_frame)
    qrd = O.E * v_body_frame
    new_qE = (O.qE + qEd * delta_t).normalize() # RBDA, bottom of p. 86
    new_qr = O.qr + qrd * delta_t
    return six_dof_euler_params(new_qE, new_qr)

  def time_step_velocity(O, v_spatial, a_spatial, delta_t):
    return v_spatial + a_spatial * delta_t

class six_dof_euler_angles_xyz(object):

  def __init__(O, qE, qr):
    if (len(qE.elems) == 4):
      qE = euler_params_qE_as_euler_angles_xyz_qE(qE=qE)
    O.qE = qE
    O.qr = qr
    #
    O.E = RBDA_Eq_4_7(qE)
    O.r = O.E.transpose() * qr # RBDA Tab. 4.1
    #
    O.T = matrix.rt((O.E, -O.E * O.r)) # RBDA Eq. 2.28
    O.T_inv = matrix.rt((O.E.transpose(), O.r))
    #
    O.Xj = featherstone.Xrot(O.E) \
         * featherstone.Xtrans(O.r) # RBDA Tab. 4.1 footnote
    O.S = None
    O.S_ring = None

  def Xj_S_S_ring(O, q, qd):
    return O.Xj, O.S, O.S_ring

  def time_step_position(O, v_spatial, delta_t):
    w_body_frame, v_body_frame = matrix.col_list([
      v_spatial.elems[:3], v_spatial.elems[3:]])
    qEd = RBDA_Eq_4_8_inv(q=O.qE.elems) * (O.E * w_body_frame)
    qrd = O.E * v_body_frame
    new_qE = O.qE + qEd * delta_t
    new_qr = O.qr + qrd * delta_t
    return six_dof_euler_angles_xyz(new_qE, new_qr)

  def time_step_velocity(O, v_spatial, a_spatial, delta_t):
    return v_spatial + a_spatial * delta_t

def RBDA_Eq_4_7(q):
  q1,q2,q3 = q
  def cs(a): return math.cos(a), math.sin(a)
  c1,s1 = cs(q1)
  c2,s2 = cs(q2)
  c3,s3 = cs(q3)
  return matrix.sqr((
              c1*c2,          s1*c2,   -s2,
     c1*s2*s3-s1*c3, s1*s2*s3+c1*c3, c2*s3,
     c1*s2*c3+s1*s3, s1*s2*c3-c1*s3, c2*c3))

def RBDA_Eq_4_8_inv(q):
  """
  Inverse of RBDA Eq. 4.8 (corresponding to Shabana G^-1)
  Mathematica code:
    ew = {{-s2, 0, 1}, {c2*s3, c3, 0}, {c2*c3, -s3, 0}}
    -Det[ew]
    Simplify[Inverse[ew] * (-Det[ew])]
  """
  q1,q2,q3 = q
  def cs(a): return math.cos(a), math.sin(a)
  c2,s2 = cs(q2)
  c3,s3 = cs(q3)
  if (abs(c2) < 1.e-6):
    raise RuntimeError("Euler angle near singular position.")
  return matrix.sqr((
    0,    s3/c2,    c3/c2,
    0,       c3,      -s3,
    1, s2*s3/c2, c3*s2/c2))

def RBDA_Eq_4_12(q):
  p0, p1, p2, p3 = q
  return matrix.sqr((
    p0**2+p1**2-0.5,   p1*p2+p0*p3,     p1*p3-p0*p2,
      p1*p2-p0*p3,   p0**2+p2**2-0.5,   p2*p3+p0*p1,
      p1*p3+p0*p2,     p2*p3-p0*p1,   p0**2+p3**2-0.5)) * 2

def RBDA_Eq_4_13(q):
  p0, p1, p2, p3 = q
  return matrix.rec((
    -p1, -p2, -p3,
    p0, -p3, p2,
    p3, p0, -p1,
    -p2, p1, p0), n=(4,3)) * 0.5

def safe_acos(a):
  return math.acos(max(-1, min(1, a)))

def safe_asin(a):
  return math.asin(max(-1, min(1, a)))

def safe_atan2(y, x):
  if (y == 0. and x == 0.): return 0.
  return math.atan2(y, x)

def euler_params_qE_as_euler_angles_xyz_qE(qE):
  p0, p1, p2, p3 = qE
  q2 = safe_asin(-2*(p1*p3-p0*p2))
  c1c2 = p0*p0+p1*p1-0.5
  s1c2 = p1*p2+p0*p3
  q1 = safe_atan2(s1c2, c1c2)
  c2c3 = p0*p0+p3*p3-0.5
  c2s3 = p2*p3+p0*p1
  q3 = safe_atan2(c2s3, c2c3)
  return matrix.col((q1,q2,q3))
