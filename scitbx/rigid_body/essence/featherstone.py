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

from __future__ import division
from spatial_lib import \
  matrix, cb_as_spatial_transform, crm, crf, kinetic_energy

try:
  import scitbx
except ImportError:
  scitbx = None
  def generalized_inverse(m):
    return m.inverse()
else:
  import scitbx.linalg.eigensystem
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
    O.bodies = bodies
    O.number_of_trees = 0
    O.degrees_of_freedom = 0
    O.q_packed_size = 0
    for body in O.bodies:
      if (body.parent == -1): O.number_of_trees += 1
      O.degrees_of_freedom += body.joint.degrees_of_freedom
      O.q_packed_size += body.joint.q_size
    O.flag_positions_as_changed()

  def flag_positions_as_changed(O):
    O.__aja_array = None
    O.__jar_array = None
    O.__cb_up_array = None
    O.__xup_array = None
    O.flag_velocities_as_changed()

  def flag_velocities_as_changed(O):
    O.__spatial_velocities = None
    O.__e_kin = None

  def root_indices(O):
    result = []
    for ib,body in enumerate(O.bodies):
      if (body.parent == -1):
        result.append(ib)
    return result

  def pack_q(O):
    from scitbx.array_family import flex
    result = flex.double()
    result.reserve(O.q_packed_size)
    for body in O.bodies:
      result.extend(flex.double(body.joint.get_q()))
    assert result.size() == O.q_packed_size
    return result

  def unpack_q(O, q_packed):
    assert q_packed.size() == O.q_packed_size
    i = 0
    for body in O.bodies:
      n = body.joint.q_size
      body.joint = body.joint.new_q(q=q_packed[i:i+n])
      i += n
    assert i == O.q_packed_size
    O.flag_positions_as_changed()

  def pack_qd(O):
    from scitbx.array_family import flex
    result = flex.double()
    result.reserve(O.degrees_of_freedom)
    for body in O.bodies:
      result.extend(flex.double(body.qd))
    assert result.size() == O.degrees_of_freedom
    return result

  def unpack_qd(O, qd_packed):
    assert qd_packed.size() == O.degrees_of_freedom
    i = 0
    for body in O.bodies:
      n = body.joint.degrees_of_freedom
      body.qd = matrix.col(qd_packed[i:i+n])
      i += n
    assert i == O.degrees_of_freedom
    O.flag_velocities_as_changed()

  def _accumulate_in_each_tree(O, attr):
    result = []
    accu = [0] * len(O.bodies)
    for ib in xrange(len(O.bodies)-1,-1,-1):
      body = O.bodies[ib]
      accu[ib] += getattr(body, attr)
      if (body.parent == -1):
        result.append((ib, accu[ib]))
      else:
        accu[body.parent] += accu[ib]
    return result

  def number_of_sites_in_each_tree(O):
    return O._accumulate_in_each_tree(attr="number_of_sites")

  def sum_of_masses_in_each_tree(O):
    return O._accumulate_in_each_tree(attr="sum_of_masses")

  def mean_linear_velocity(O, number_of_sites_in_each_tree):
    if (number_of_sites_in_each_tree is None):
      number_of_sites_in_each_tree = O.number_of_sites_in_each_tree()
    sum_v = matrix.col((0,0,0))
    sum_n = 0
    for ib,n in number_of_sites_in_each_tree:
      body = O.bodies[ib]
      v = body.joint.get_linear_velocity(qd=body.qd)
      if (v is None): continue
      sum_v += v * n
      sum_n += n
    if (sum_n == 0):
      return None
    return sum_v / sum_n

  def subtract_from_linear_velocities(O, number_of_sites_in_each_tree, value):
    if (number_of_sites_in_each_tree is None):
      number_of_sites_in_each_tree = O.number_of_sites_in_each_tree()
    for ib,n in number_of_sites_in_each_tree:
      body = O.bodies[ib]
      v = body.joint.get_linear_velocity(qd=body.qd)
      if (v is None): continue
      body.qd = body.joint.new_linear_velocity(qd=body.qd, value=v-value)

  def aja_array(O):
    if (O.__aja_array is None):
      O.__aja_array = []
      for body in O.bodies:
        aja = body.alignment.cb_b0 * body.joint.cb_sp * body.alignment.cb_0b
        if (body.parent != -1):
          aja = O.__aja_array[body.parent] * aja
        O.__aja_array.append(aja)
    return O.__aja_array

  def jar_array(O):
    if (O.__jar_array is None):
      O_aja = O.aja_array()
      O.__jar_array = []
      for body in O.bodies:
        jar = body.joint.cb_ps.r * body.alignment.cb_0b.r
        if (body.parent != -1):
          jar *= O_aja[body.parent].r.transpose()
        O.__jar_array.append(jar)
    return O.__jar_array

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
    if (O.__spatial_velocities is None):
      O.__spatial_velocities = []
      cb_up_array = O.cb_up_array()
      for ib in xrange(len(O.bodies)):
        body = O.bodies[ib]
        s = body.joint.motion_subspace
        if (s is None): vj = body.qd
        else:           vj = s * body.qd
        if (body.parent == -1):
          O.__spatial_velocities.append(vj)
        else:
          cb_up = cb_up_array[ib]
          vp = O.__spatial_velocities[body.parent].elems
          r_va = cb_up.r * matrix.col(vp[:3])
          vp = matrix.col((
            r_va,
            cb_up.r * matrix.col(vp[3:]) + cb_up.t.cross(r_va))) \
              .resolve_partitions()
          O.__spatial_velocities.append(vp + vj)
    return O.__spatial_velocities

  def e_kin(O):
    "RBDA Eq. 2.67, p. 35"
    if (O.__e_kin is None):
      result = 0
      for body,v in zip(O.bodies, O.spatial_velocities()):
        result += kinetic_energy(i_spatial=body.i_spatial, v_spatial=v)
      O.__e_kin = result
    return O.__e_kin

  def reset_e_kin(O, e_kin_target, e_kin_epsilon=1e-12):
    assert e_kin_target >= 0
    assert e_kin_epsilon > 0
    O_e_kin = O.e_kin()
    if (O_e_kin >= e_kin_epsilon):
      factor = (e_kin_target / O_e_kin)**0.5
      for body in O.bodies:
        body.qd *= factor
    O.flag_velocities_as_changed()

  def assign_zero_velocities(O):
    for body in O.bodies:
      body.qd = body.joint.qd_zero
    O.flag_velocities_as_changed()

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

  def assign_random_velocities(O,
        e_kin_target=None,
        e_kin_epsilon=1e-12,
        random_gauss=None):
    if (e_kin_target is None):
      work_e_kin_target = 1
    elif (e_kin_target == 0):
      O.assign_zero_velocities()
      return
    else:
      assert e_kin_target >= 0
      work_e_kin_target = e_kin_target
    from scitbx.array_family import flex
    qd_e_kin_scales = flex.double(
      O.qd_e_kin_scales(e_kin_epsilon=e_kin_epsilon))
    if (O.degrees_of_freedom != 0):
      qd_e_kin_scales *= (work_e_kin_target / O.degrees_of_freedom)**0.5
    if (random_gauss is None):
      import random
      random_gauss = random.gauss
    i_qd = 0
    for body in O.bodies:
      qd_new = []
      for qd in body.joint.qd_zero:
        qd_new.append(qd + random_gauss(mu=0, sigma=qd_e_kin_scales[i_qd]))
        i_qd += 1
      body.qd = matrix.col(qd_new)
    assert i_qd == O.degrees_of_freedom
    O.flag_velocities_as_changed()
    if (e_kin_target is not None):
      O.reset_e_kin(e_kin_target=e_kin_target, e_kin_epsilon=e_kin_epsilon)
    return qd_e_kin_scales

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
