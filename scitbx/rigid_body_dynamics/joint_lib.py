from scitbx.rigid_body_dynamics import featherstone
from scitbx import matrix

class six_dof_joint_euler_params(object):

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
    O.S = matrix.sqr(( # RBDA Tab. 4.1
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1))
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
    return six_dof_joint_euler_params(new_qE, new_qr)

  def time_step_velocity(O, v_spatial, a_spatial, delta_t):
    return v_spatial + a_spatial * delta_t

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
