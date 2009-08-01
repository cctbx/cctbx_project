"""\
Based in part on Roy Featherstone's spatial_v1 matlab code:

  http://axiom.anu.edu.au/~roy/spatial/

  Version 1: January 2008 (latest bug fix: 7 October 2008)

See also: RBDA:
  Rigid Body Dynamics Algorithms.
  Roy Featherstone,
  Springer, New York, 2007.
  ISBN-10: 0387743146
"""

from spatial_lib import \
  matrix, cb_as_spatial_transform, crm, crf, kinetic_energy

try:
  import scitbx
except ImportError:
  scitbx = None
  def generalized_inverse(m):
    return m.inverse()
else:
  import scitbx.math
  def generalized_inverse(m):
    # assumption to achieve stability: order of magnitude of masses is around 1
    return matrix.sqr(
      scitbx.linalg.eigensystem.real_symmetric(
        m=m.as_flex_double_matrix(),
        relative_epsilon=1e-6,
        absolute_epsilon=1e-6)
          .generalized_inverse_as_packed_u()
          .matrix_packed_u_as_symmetric())

class system_model(object):
  "RBDA Tab. 4.3, p. 87"

  def __init__(O, bodies):
    "Stores bodies and computes Xtree (RBDA Fig. 4.7, p. 74)"
    O.bodies = bodies
    for body in bodies:
      if (body.parent == -1):
        cb_tree = body.alignment.cb_0b
      else:
        cb_tree = body.alignment.cb_0b * bodies[body.parent].alignment.cb_b0
      body.cb_tree = cb_tree
    O.__cb_up_array = None
    O.__xup_array = None

  def cb_up_array(O):
    "RBDA Example 4.4, p. 80"
    if (O.__cb_up_array is None):
      O.__cb_up_array = [body.joint.cb_ps * body.cb_tree for body in O.bodies]
    return O.__cb_up_array

  def xup_array(O):
    if (O.__xup_array is None):
      O.__xup_array = [cb_as_spatial_transform(cb) for cb in O.cb_up_array()]
    return O.__xup_array

  def spatial_velocities(O):
    "RBDA Example 4.4, p. 80"
    result = []
    cb_up_array = O.cb_up_array()
    for ib in xrange(len(O.bodies)):
      body = O.bodies[ib]
      s = body.joint.motion_subspace
      if (s is None): vj = body.qd
      else:           vj = s * body.qd
      if (body.parent == -1): result.append(vj)
      else:
        cb_up = cb_up_array[ib]
        vp = result[body.parent].elems
        r_va = cb_up.r * matrix.col(vp[:3])
        vp = matrix.col((
          r_va,
          cb_up.r * matrix.col(vp[3:]) + cb_up.t.cross(r_va))) \
            .resolve_partitions()
        result.append(vp + vj)
    return result

  def e_kin(O):
    "RBDA Eq. 2.67, p. 35"
    result = 0
    for body,v in zip(O.bodies, O.spatial_velocities()):
      result += kinetic_energy(i_spatial=body.i_spatial, v_spatial=v)
    return result

  def accumulated_spatial_inertia(O):
    result = [body.i_spatial for body in O.bodies]
    xup_array = O.xup_array()
    for i in xrange(len(result)-1,-1,-1):
      body = O.bodies[i]
      if (body.parent != -1):
        result[body.parent] \
          += xup_array[i].transpose() * result[i] * xup_array[i]
    return result

  def qd_e_kin_scales(O, e_kin_epsilon=1.e-12):
    result = []
    rap = result.append
    for body,asi in zip(O.bodies, O.accumulated_spatial_inertia()):
      s = body.joint.motion_subspace
      j_dof = body.joint.degrees_of_freedom
      qd = [0] * j_dof
      for i in xrange(j_dof):
        qd[i] = 1
        vj = matrix.col(qd)
        qd[i] = 0
        if (s is not None): vj = s * vj
        e_kin = kinetic_energy(i_spatial=asi, v_spatial=vj)
        if (e_kin < e_kin_epsilon):
          rap(1)
        else:
          rap(1 / e_kin**0.5)
    return result

  def inverse_dynamics(O, qdd_array=None, f_ext_array=None, grav_accn=None):
    """RBDA Tab. 5.1, p. 96:
Inverse Dynamics of a kinematic tree via Recursive Newton-Euler Algorithm.
qdd_array is a vector of joint acceleration variables.
The return value (tau) is a vector of joint force variables.
f_ext_array specifies external forces acting on the bodies. If f_ext_array
is None then there are no external forces; otherwise, f_ext[i] is a spatial
force vector giving the force acting on body i, expressed in body i
coordinates.
grav_accn is a 6D vector expressing the linear acceleration due to gravity.
    """
    assert qdd_array is None or len(qdd_array) == len(O.bodies)
    assert f_ext_array is None or len(f_ext_array) == len(O.bodies)
    nb = len(O.bodies)
    xup_array = O.xup_array()
    v = O.spatial_velocities()
    a = [None] * nb
    f = [None] * nb
    for ib in xrange(nb):
      body = O.bodies[ib]
      s = body.joint.motion_subspace
      if (s is None):
        vj = body.qd
        if (qdd_array is None or qdd_array[ib] is None):
          aj = matrix.zeros(n=6)
        else:
          aj = qdd_array[ib]
      else:
        vj = s * body.qd
        if (qdd_array is None or qdd_array[ib] is None):
          aj = matrix.zeros(n=6)
        else:
          aj = s * qdd_array[ib]
      if (body.parent == -1):
        a[ib] = aj
        if (grav_accn is not None):
          a[ib] += xup_array[ib] * (-grav_accn)
      else:
        a[ib] = xup_array[ib] * a[body.parent] + aj + crm(v[ib]) * vj
      f[ib] = body.i_spatial * a[ib] + crf(v[ib]) * (body.i_spatial * v[ib])
      if (f_ext_array is not None and f_ext_array[ib] is not None):
        f[ib] -= f_ext_array[ib]
    #
    tau_array = [None] * nb
    for ib in xrange(nb-1,-1,-1):
      body = O.bodies[ib]
      s = body.joint.motion_subspace
      if (s is None):
        tau_array[ib] = f[ib]
      else:
        tau_array[ib] = s.transpose() * f[ib]
      if (body.parent != -1):
        f[body.parent] += xup_array[ib].transpose() * f[ib]
    #
    return tau_array

  def f_ext_as_tau(O, f_ext_array):
    """
Simplified version of Inverse Dynamics via Recursive Newton-Euler Algorithm,
with all qd, qdd zero, but non-zero external forces.
    """
    assert len(f_ext_array) == len(O.bodies)
    nb = len(O.bodies)
    xup_array = O.xup_array()
    f = [-e for e in f_ext_array]
    tau_array = [None] * nb
    for ib in xrange(nb-1,-1,-1):
      body = O.bodies[ib]
      s = body.joint.motion_subspace
      if (s is None): tau_array[ib] = f[ib]
      else:           tau_array[ib] = s.transpose() * f[ib]
      if (body.parent != -1):
        f[body.parent] += xup_array[ib].transpose() * f[ib]
    return tau_array

  def d_e_pot_d_q(O, f_ext_array):
    """
Gradients of potential energy (defined via f_ext_array) w.r.t. positional
coordinates q. Uses f_ext_as_tau().
    """
    return [body.joint.tau_as_d_e_pot_d_q(tau=tau)
      for body,tau in zip(O.bodies, O.f_ext_as_tau(f_ext_array=f_ext_array))]

  def forward_dynamics_ab(O, tau_array=None, f_ext_array=None, grav_accn=None):
    """RBDA Tab. 7.1, p. 132:
Forward Dynamics of a kinematic tree via the Articulated-Body Algorithm.
tau_array is a vector of force variables.
The return value (qdd_array) is a vector of joint acceleration variables.
f_ext_array specifies external forces acting on the bodies. If f_ext_array
is None then there are no external forces; otherwise, f_ext_array[i] is a
spatial force vector giving the force acting on body i, expressed in body i
coordinates.
grav_accn is a 6D vector expressing the linear acceleration due to gravity.
    """
    assert tau_array is None or len(tau_array) == len(O.bodies)
    assert f_ext_array is None or len(f_ext_array) == len(O.bodies)
    nb = len(O.bodies)
    xup_array = O.xup_array()
    v = O.spatial_velocities()
    c = [None] * nb
    ia = [None] * nb
    pa = [None] * nb
    for ib in xrange(nb):
      body = O.bodies[ib]
      s = body.joint.motion_subspace
      if (s is None): vj = body.qd
      else:           vj = s * body.qd
      if (body.parent == -1): c[ib] = matrix.col([0,0,0,0,0,0])
      else:                   c[ib] = crm(v[ib]) * vj
      ia[ib] = body.i_spatial
      pa[ib] = crf(v[ib]) * (body.i_spatial * v[ib])
      if (f_ext_array is not None and f_ext_array[ib] is not None):
        pa[ib] -= f_ext_array[ib]
    #
    u = [None] * nb
    d_inv = [None] * nb
    u_ = [None] * nb
    for ib in xrange(nb-1,-1,-1):
      body = O.bodies[ib]
      s = body.joint.motion_subspace
      if (s is None):
        u[ib] = ia[ib]
        d = u[ib]
        u_[ib] = -pa[ib]
      else:
        u[ib] = ia[ib] * s
        d = s.transpose() * u[ib]
        u_[ib] = -s.transpose() * pa[ib]
      if (tau_array is not None and tau_array[ib] is not None):
        u_[ib] += tau_array[ib]
      d_inv[ib] = generalized_inverse(d)
      if (body.parent != -1):
        u_d_inv = u[ib] * d_inv[ib];
        ia_ = ia[ib] - u_d_inv * u[ib].transpose()
        pa_ = pa[ib] + ia_*c[ib] + u_d_inv * u_[ib]
        ia[body.parent] += xup_array[ib].transpose() * ia_ * xup_array[ib]
        pa[body.parent] += xup_array[ib].transpose() * pa_
    #
    a = [None] * nb
    qdd_array = [None] * nb
    for ib in xrange(nb):
      body = O.bodies[ib]
      s = body.joint.motion_subspace
      if (body.parent == -1):
        a[ib] = c[ib]
        if (grav_accn is not None):
          a[ib] += xup_array[ib] * (-grav_accn)
      else:
        a[ib] = xup_array[ib] * a[body.parent] + c[ib]
      qdd_array[ib] = d_inv[ib] * (u_[ib] - u[ib].transpose()*a[ib])
      if (s is None):
        a[ib] += qdd_array[ib]
      else:
        a[ib] += s * qdd_array[ib]
    #
    return qdd_array
