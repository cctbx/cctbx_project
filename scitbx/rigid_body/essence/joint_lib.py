from featherstone import matrix, T_as_X
from utils import center_of_mass_from_sites
import math

class six_dof_alignment(object):

  def __init__(O, sites):
    O.pivot = center_of_mass_from_sites(sites=sites)
    O.normal = None
    O.T0b = matrix.rt(((1,0,0,0,1,0,0,0,1), -O.pivot))
    O.Tb0 = matrix.rt(((1,0,0,0,1,0,0,0,1), O.pivot))

class six_dof(object):

  qd_zero = matrix.zeros(n=6)
  qdd_zero = matrix.zeros(n=6)

  def __init__(O, qE, qr):
    O.qE = qE
    O.qr = qr
    O.q_size = 7
    O.unit_quaternion = qE.normalize() # RBDA, bottom of p. 86
    O.E = RBDA_Eq_4_12(q=O.unit_quaternion)
    O.r = qr
    O.Tps = matrix.rt((O.E, -O.E * O.r)) # RBDA Eq. 2.28
    O.Tsp = matrix.rt((O.E.transpose(), O.r))
    O.Xj = T_as_X(O.Tps)
    O.S = None

  def time_step_position(O, qd, delta_t):
    w_body_frame, v_body_frame = matrix.col_list([qd.elems[:3], qd.elems[3:]])
    qEd = RBDA_Eq_4_13(q=O.unit_quaternion) * w_body_frame
    new_qE = (O.qE + qEd * delta_t).normalize()
    qrd = O.E.transpose() * v_body_frame
    new_qr = O.qr + qrd * delta_t
    return six_dof(qE=new_qE, qr=new_qr)

  def time_step_velocity(O, qd, qdd, delta_t):
    return qd + qdd * delta_t

  def tau_as_d_pot_d_q(O, tau):
    d = d_unit_quaternion_d_qE_matrix(q=O.qE)
    c = d * 4 * RBDA_Eq_4_13(q=O.unit_quaternion)
    n, f = matrix.col_list([tau.elems[:3], tau.elems[3:]])
    return matrix.col((c * n, O.E.transpose() * f)).resolve_partitions()

  def get_q(O):
    return O.qE.elems + O.qr.elems

  def new_q(O, q):
    new_qE, new_qr = matrix.col_list((q[:4], q[4:]))
    return six_dof(qE=new_qE, qr=new_qr)

class spherical_alignment(object):

  def __init__(O, sites):
    O.pivot = center_of_mass_from_sites(sites=sites)
    O.normal = None
    O.T0b = matrix.rt(((1,0,0,0,1,0,0,0,1), -O.pivot))
    O.Tb0 = matrix.rt(((1,0,0,0,1,0,0,0,1), O.pivot))

class spherical(object):

  qd_zero = matrix.zeros(n=3)
  qdd_zero = matrix.zeros(n=3)

  def __init__(O, qE):
    O.qE = qE
    O.q_size = 4
    O.unit_quaternion = qE.normalize() # RBDA, bottom of p. 86
    O.E = RBDA_Eq_4_12(q=O.unit_quaternion)
    O.Tps = matrix.rt((O.E, (0,0,0)))
    O.Tsp = matrix.rt((O.E.transpose(), (0,0,0)))
    O.Xj = T_as_X(O.Tps)
    O.S = matrix.rec((
      1,0,0,
      0,1,0,
      0,0,1,
      0,0,0,
      0,0,0,
      0,0,0), n=(6,3))

  def time_step_position(O, qd, delta_t):
    w_body_frame = qd
    d = d_unit_quaternion_d_qE_matrix(q=O.qE)
    qEd = d * RBDA_Eq_4_13(q=O.unit_quaternion) * w_body_frame
    new_qE = O.qE + qEd * delta_t
    return spherical(qE=new_qE)

  def time_step_velocity(O, qd, qdd, delta_t):
    return qd + qdd * delta_t

  def tau_as_d_pot_d_q(O, tau):
    d = d_unit_quaternion_d_qE_matrix(q=O.qE)
    c = d * 4 * RBDA_Eq_4_13(q=O.unit_quaternion)
    n = tau
    return c * n

  def get_q(O):
    return O.qE.elems

  def new_q(O, q):
    return spherical(qE=matrix.col(q))

class revolute_alignment(object):

  def __init__(O, pivot, normal):
    O.pivot = pivot
    O.normal = normal
    r = normal.vector_to_001_rotation()
    O.T0b = matrix.rt((r, -r * pivot))
    O.Tb0 = matrix.rt((r.transpose(), pivot))

class revolute(object):

  qd_zero = matrix.zeros(n=1)
  qdd_zero = matrix.zeros(n=1)

  def __init__(O, qE):
    O.qE = qE
    O.q_size = len(qE)
    #
    c, s = math.cos(qE[0]), math.sin(qE[0])
    O.E = matrix.sqr((c, s, 0, -s, c, 0, 0, 0, 1)) # RBDA Tab. 2.2
    O.r = matrix.col((0,0,0))
    #
    O.Tps = matrix.rt((O.E, (0,0,0)))
    O.Tsp = matrix.rt((O.E.transpose(), (0,0,0)))
    O.Xj = T_as_X(O.Tps)
    O.S = matrix.col((0,0,1,0,0,0))

  def time_step_position(O, qd, delta_t):
    new_qE = O.qE + qd * delta_t
    return revolute(qE=new_qE)

  def time_step_velocity(O, qd, qdd, delta_t):
    return qd + qdd * delta_t

  def tau_as_d_pot_d_q(O, tau):
    return tau

  def get_q(O):
    return O.qE.elems

  def new_q(O, q):
    return revolute(qE=matrix.col(q))

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

def d_unit_quaternion_d_qE_matrix(q):
  """
  Coefficent matrix for converting gradients w.r.t. normalized Euler
  parameters to gradients w.r.t. non-normalized parameters, as produced
  e.g. by a minimizer in the line search.
  Mathematica code:
    nsq = p0^2+p1^2+p2^2+p3^2
    p0p = p0 / Sqrt[nsq]
    p1p = p1 / Sqrt[nsq]
    p2p = p2 / Sqrt[nsq]
    p3p = p3 / Sqrt[nsq]
    n3 = (p0^2+p1^2+p2^2+p3^2)^(3/2)
    FortranForm[FullSimplify[D[p0p,p0]*n3]]
    FortranForm[FullSimplify[D[p1p,p0]*n3]]
    FortranForm[FullSimplify[D[p2p,p0]*n3]]
    FortranForm[FullSimplify[D[p3p,p0]*n3]]
    FortranForm[FullSimplify[D[p0p,p1]*n3]]
    FortranForm[FullSimplify[D[p1p,p1]*n3]]
    FortranForm[FullSimplify[D[p2p,p1]*n3]]
    FortranForm[FullSimplify[D[p3p,p1]*n3]]
    FortranForm[FullSimplify[D[p0p,p2]*n3]]
    FortranForm[FullSimplify[D[p1p,p2]*n3]]
    FortranForm[FullSimplify[D[p2p,p2]*n3]]
    FortranForm[FullSimplify[D[p3p,p2]*n3]]
    FortranForm[FullSimplify[D[p0p,p3]*n3]]
    FortranForm[FullSimplify[D[p1p,p3]*n3]]
    FortranForm[FullSimplify[D[p2p,p3]*n3]]
    FortranForm[FullSimplify[D[p3p,p3]*n3]]
  """
  p0,p1,p2,p3 = q
  p0s,p1s,p2s,p3s = p0**2, p1**2, p2**2, p3**2
  n3 = (p0s+p1s+p2s+p3s)**(3/2.)
  c00 = p1s+p2s+p3s
  c11 = p0s+p2s+p3s
  c22 = p0s+p1s+p3s
  c33 = p0s+p1s+p2s
  c01 = -p0*p1
  c02 = -p0*p2
  c03 = -p0*p3
  c12 = -p1*p2
  c13 = -p1*p3
  c23 = -p2*p3
  return matrix.sqr((
    c00, c01, c02, c03,
    c01, c11, c12, c13,
    c02, c12, c22, c23,
    c03, c13, c23, c33)) / n3
