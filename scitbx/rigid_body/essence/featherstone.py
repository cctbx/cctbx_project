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

try:
  import scitbx
except ImportError:
  scitbx = None

if (scitbx is not None):
  import scitbx.math
  from scitbx import matrix

  def generalized_inverse(m):
    # assumption to achieve stability: order of magnitude of masses is around 1
    return matrix.sqr(
      scitbx.linalg.eigensystem.real_symmetric(
        m=m.as_flex_double_matrix(),
        relative_epsilon=1e-6,
        absolute_epsilon=1e-6)
          .generalized_inverse_as_packed_u()
          .matrix_packed_u_as_symmetric())

else:
  import scitbx_matrix as matrix

  def generalized_inverse(m):
    return m.inverse()

def xrot(E):
  """RBDA Tab. 2.2, p. 23:
Spatial coordinate transform (rotation around origin).
Calculates the coordinate transform matrix from A to B coordinates
for spatial motion vectors, in which frame B is rotated relative to
frame A.
  """
  a,b,c,d,e,f,g,h,i = E
  return matrix.sqr((
     a,  b,  c,  0,  0,  0,
     d,  e,  f,  0,  0,  0,
     g,  h,  i,  0,  0,  0,
     0,  0,  0,  a,  b,  c,
     0,  0,  0,  d,  e,  f,
     0,  0,  0,  g,  h,  i))

def xtrans(r):
  """RBDA Tab. 2.2, p. 23:
Spatial coordinate transform (translation of origin).
Calculates the coordinate transform matrix from A to B coordinates
for spatial motion vectors, in which frame B is translated by an
amount r (3D vector) relative to frame A.
  """
  r1,r2,r3 = r
  return matrix.sqr((
      1,   0,   0, 0, 0, 0,
      0,   1,   0, 0, 0, 0,
      0,   0,   1, 0, 0, 0,
      0,  r3, -r2, 1, 0, 0,
    -r3,   0,  r1, 0, 1, 0,
     r2, -r1,   0, 0, 0, 1))

def cb_as_spatial_transform(cb):
  """RBDA Eq. 2.28, p. 22:
Conversion of matrix.rt object cb to spatial transform.
"""
  return xrot(cb.r) * xtrans(-cb.r.transpose() * cb.t)

def crm(v):
  """RBDA Eq. 2.31, p. 25:
Spatial cross-product operator (motion).
Calculates the 6x6 matrix such that the expression crm(v)*m is the
cross product of the spatial motion vectors v and m.
  """
  v1,v2,v3,v4,v5,v6 = v
  return matrix.sqr((
      0, -v3,  v2,   0,   0,   0,
     v3,   0, -v1,   0,   0,   0,
    -v2,  v1,   0,   0,   0,   0,
      0, -v6,  v5,   0, -v3,  v2,
     v6,   0, -v4,  v3,   0, -v1,
    -v5,  v4,   0, -v2,  v1,   0))

def crf(v):
  """RBDA Eq. 2.32, p. 25:
Spatial cross-product operator (force).
Calculates the 6x6 matrix such that the expression crf(v)*f is the
cross product of the spatial motion vector v with the spatial force
vector f.
  """
  return -crm(v).transpose()

def mcI(m, c, I):
  """RBDA Eq. 2.63, p. 33:
Spatial rigid-body inertia from mass, CoM and rotational inertia.
Calculates the spatial inertia matrix of a rigid body from its
mass, centre of mass (3D vector) and rotational inertia (3x3 matrix)
about its centre of mass.
  """
  c1,c2,c3 = c
  C = matrix.sqr((
      0, -c3,  c2,
     c3,   0, -c1,
    -c2,  c1,  0))
  return matrix.sqr((
    I + m*C*C.transpose(), m*C,
    m*C.transpose(), m*matrix.identity(3))).resolve_partitions()

def kinetic_energy(i_spatial, v_spatial):
  "RBDA Eq. 2.67, p. 35"
  return 0.5 * v_spatial.dot(i_spatial * v_spatial)

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
    for i,body in enumerate(O.bodies):
      s = body.joint.S
      if (s is None): vj = body.qd
      else:           vj = s * body.qd
      if (body.parent == -1): result.append(vj)
      else:
        cb_up = cb_up_array[i]
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
      s = body.joint.S
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

  def inverse_dynamics(O, qdd_array, f_ext_array=None, grav_accn=None):
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

    nb = len(O.bodies)
    xup_array = O.xup_array()
    v = O.spatial_velocities()
    a = [None] * nb
    f = [None] * nb
    for ib in xrange(nb):
      body = O.bodies[ib]
      s = body.joint.S
      if (s is None):
        vj = body.qd
        aj = qdd_array[ib]
      else:
        vj = s * body.qd
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

    tau_array = [None] * nb
    for ib in xrange(nb-1,-1,-1):
      body = O.bodies[ib]
      s = body.joint.S
      if (s is None):
        tau_array[ib] = f[ib]
      else:
        tau_array[ib] = s.transpose() * f[ib]
      if (body.parent != -1):
        f[body.parent] += xup_array[ib].transpose() * f[ib]

    return tau_array

  def ID0(O, f_ext):
    """
Simplified version of Inverse Dynamics via Recursive Newton-Euler Algorithm,
with all qd, qdd zero, but non-zero external forces.
    """
    nb = len(O.bodies)
    xup_array = O.xup_array()
    f = [-e for e in f_ext]
    tau = [None] * nb
    for i in xrange(nb-1,-1,-1):
      body = O.bodies[i]
      if (body.joint.S is None): tau[i] = f[i]
      else:                      tau[i] = body.joint.S.transpose() * f[i]
      if (body.parent != -1):
        f[body.parent] += xup_array[i].transpose() * f[i]
    return tau

  def d_pot_d_q(O, f_ext):
    """
Gradients of potential energy (defined via f_ext) w.r.t. positional
coordinates q. Uses ID0().
    """
    return [body.joint.tau_as_d_pot_d_q(tau=tau)
      for body,tau in zip(O.bodies, O.ID0(f_ext=f_ext))]

  def FDab(O, tau=None, f_ext=None, grav_accn=None):
    """RBDA Tab. 7.1, p. 132:
Forward Dynamics of a kinematic tree via the Articulated-Body Algorithm.
tau is a vector of force variables.
The return value (qdd) is a vector of joint acceleration variables.
f_ext specifies external forces acting on the bodies. If f_ext is None
then there are no external forces; otherwise, f_ext[i] is a spatial force
vector giving the force acting on body i, expressed in body i coordinates.
grav_accn is a 6D vector expressing the linear acceleration due to gravity.
    """

    nb = len(O.bodies)
    xup_array = O.xup_array()
    v = O.spatial_velocities()
    c = [None] * nb
    IA = [None] * nb
    pA = [None] * nb
    for i in xrange(nb):
      body = O.bodies[i]
      if (body.joint.S is None): vJ = body.qd
      else:                      vJ = body.joint.S * body.qd
      if (body.parent == -1): c[i] = matrix.col([0,0,0,0,0,0])
      else:                c[i] = crm(v[i]) * vJ
      IA[i] = body.i_spatial
      pA[i] = crf(v[i]) * body.i_spatial * v[i]
      if (f_ext is not None and f_ext[i] is not None):
        pA[i] -= f_ext[i]

    U = [None] * nb
    d_inv = [None] * nb
    u = [None] * nb
    for i in xrange(nb-1,-1,-1):
      body = O.bodies[i]
      if (body.joint.S is None):
        U[i] = IA[i]
        d = U[i]
        if (tau is None or tau[i] is None):
          u[i] =        - pA[i]
        else:
          u[i] = tau[i] - pA[i]
      else:
        U[i] = IA[i] * body.joint.S
        d = body.joint.S.transpose() * U[i]
        if (tau is None or tau[i] is None):
          u[i] =        - body.joint.S.transpose() * pA[i]
        else:
          u[i] = tau[i] - body.joint.S.transpose() * pA[i]
      d_inv[i] = generalized_inverse(d)
      if (body.parent != -1):
        Ia = IA[i] - U[i] * d_inv[i] * U[i].transpose()
        pa = pA[i] + Ia*c[i] + U[i] * d_inv[i] * u[i]
        IA[body.parent] += xup_array[i].transpose() * Ia * xup_array[i]
        pA[body.parent] += xup_array[i].transpose() * pa

    a = [None] * nb
    qdd = [None] * nb
    for i in xrange(nb):
      body = O.bodies[i]
      if (body.parent == -1):
        a[i] = c[i]
        if (grav_accn is not None):
          a[i] += xup_array[i] * (-grav_accn)
      else:
        a[i] = xup_array[i] * a[body.parent] + c[i]
      qdd[i] = d_inv[i] * (u[i] - U[i].transpose()*a[i])
      if (body.joint.S is None):
        a[i] += qdd[i]
      else:
        a[i] += body.joint.S * qdd[i]

    return qdd
