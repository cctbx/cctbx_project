from scitbx.rigid_body_dynamics import featherstone
from scitbx.rigid_body_dynamics.utils import center_of_mass_from_sites
from scitbx import matrix
import math

class six_dof_alignment(object):

  def __init__(O, sites):
    O.pivot = center_of_mass_from_sites(sites=sites)
    O.normal = None
    O.T0b = matrix.rt(((1,0,0,0,1,0,0,0,1), -O.pivot))
    O.Tb0 = matrix.rt(((1,0,0,0,1,0,0,0,1), O.pivot))

class six_dof(object):

  def __init__(O, type, qE, qr, r_is_qr=False):
    assert type in ["euler_params", "euler_angles_xyz"]
    if (type == "euler_params"):
      if (len(qE.elems) == 3):
        qE = euler_angles_xyz_qE_as_euler_params_qE(qE=qE)
    else:
      if (len(qE.elems) == 4):
        qE = euler_params_qE_as_euler_angles_xyz_qE(qE=qE)
    O.type = type
    O.qE = qE
    O.qr = qr
    O.r_is_qr = r_is_qr
    #
    if (type == "euler_params"):
      O.E = RBDA_Eq_4_12(q=qE)
    else:
      O.E = RBDA_Eq_4_7(q=qE)
    if (r_is_qr):
      O.r = qr
    else:
      O.r = O.E.transpose() * qr # RBDA Tab. 4.1
    #
    O.Tps = matrix.rt((O.E, -O.E * O.r)) # RBDA Eq. 2.28
    O.Tsp = matrix.rt((O.E.transpose(), O.r))
    O.Xj = T_as_X(O.Tps)
    O.S = None

  def Xj_S(O, q, qd):
    return O.Xj, O.S

  def time_step_position(O, qd, delta_t):
    w_body_frame, v_body_frame = matrix.col_list([qd.elems[:3], qd.elems[3:]])
    if (O.type == "euler_params"):
      qEd = RBDA_Eq_4_13(q=O.qE.elems) * w_body_frame
    else:
      qEd = RBDA_Eq_4_8_inv(q=O.qE.elems) * w_body_frame
    if (O.r_is_qr):
      qrd = O.E.transpose() * v_body_frame
    else:
      qrd = v_body_frame - w_body_frame.cross(O.qr) # RBDA Eq. 2.38 p. 27
    if (O.type == "euler_params"):
      new_qE = (O.qE + qEd * delta_t).normalize() # RBDA, bottom of p. 86
    else:
      new_qE = O.qE + qEd * delta_t
    new_qr = O.qr + qrd * delta_t
    return six_dof(O.type, new_qE, new_qr, O.r_is_qr)

  def time_step_velocity(O, qd, qdd, delta_t):
    return qd + qdd * delta_t

class five_dof_alignment(object):

  def __init__(O, sites):
    assert len(sites) == 2
    O.pivot = center_of_mass_from_sites(sites=sites)
    O.normal = (sites[1] - sites[0]).normalize()
    r = O.normal.vector_to_001_rotation()
    O.T0b = matrix.rt((r, -r * O.pivot))
    O.Tb0 = matrix.rt((r.transpose(), O.pivot))

class five_dof(object):

  def __init__(O, qE, qr, r_is_qr=False):
    assert len(qE) == 2
    O.qE = qE
    O.qr = qr
    O.r_is_qr = r_is_qr
    #
    O.E = RBDA_Eq_4_7(q=(0, qE[0], qE[1]))
    if (r_is_qr):
      O.r = qr
    else:
      O.r = O.E.transpose() * qr # RBDA Tab. 4.1
    #
    O.Tps = matrix.rt((O.E, -O.E * O.r)) # RBDA Eq. 2.28
    O.Tsp = matrix.rt((O.E.transpose(), O.r))
    O.Xj = T_as_X(O.Tps)
    O.S = matrix.rec((
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1), n=(6,5))

  def Xj_S(O, q, qd):
    return O.Xj, O.S

  def time_step_position(O, qd, delta_t):
    w_body_frame, v_body_frame = matrix.col_list([
      qd.elems[:2]+(0,), qd.elems[2:]])
    qEd = matrix.col(
      (RBDA_Eq_4_8_inv(q=(0, O.qE[0], O.qE[1])) * w_body_frame).elems[1:])
    if (O.r_is_qr):
      qrd = O.E.transpose() * v_body_frame
    else:
      qrd = v_body_frame - w_body_frame.cross(O.qr) # RBDA Eq. 2.38 p. 27
    new_qE = O.qE + qEd * delta_t
    new_qr = O.qr + qrd * delta_t
    return five_dof(new_qE, new_qr, O.r_is_qr)

  def time_step_velocity(O, qd, qdd, delta_t):
    return qd + qdd * delta_t

class revolute_alignment(object):

  def __init__(O, pivot, normal):
    O.pivot = pivot
    O.normal = normal
    r = normal.vector_to_001_rotation()
    O.T0b = matrix.rt((r, -r * pivot))
    O.Tb0 = matrix.rt((r.transpose(), pivot))

class revolute(object):

  def __init__(O, qE):
    O.qE = qE
    #
    c, s = math.cos(qE[0]), math.sin(qE[0])
    O.E = matrix.sqr((c, s, 0, -s, c, 0, 0, 0, 1)) # RBDA Tab. 2.2
    O.r = matrix.col((0,0,0))
    #
    O.Tps = matrix.rt((O.E, (0,0,0)))
    O.Tsp = matrix.rt((O.E.transpose(), (0,0,0)))
    O.Xj = T_as_X(O.Tps)
    O.S = matrix.col((0,0,1,0,0,0))

  def Xj_S(O, q, qd):
    return O.Xj, O.S

  def time_step_position(O, qd, delta_t):
    new_qE = O.qE + qd * delta_t
    return revolute(qE=new_qE)

  def time_step_velocity(O, qd, qdd, delta_t):
    return qd + qdd * delta_t

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

def euler_angles_xyz_qE_as_euler_params_qE(qE):
  s1,s2,s3 = [math.sin(q/2) for q in qE]
  c1,c2,c3 = [math.cos(q/2) for q in qE]
  return matrix.col((
    c1*c2*c3+s1*s2*s3,
    c1*c2*s3-s1*s2*c3,
    c1*s2*c3+s1*c2*s3,
    s1*c2*c3-c1*s2*s3))

def T_as_X(Tps):
  return featherstone.Xrot(Tps.r) \
       * featherstone.Xtrans(-Tps.r.transpose() * Tps.t)
