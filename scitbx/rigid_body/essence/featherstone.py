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
      scitbx.math.eigensystem.real_symmetric(
        m=m.as_flex_double_matrix(),
        relative_epsilon=1e-6,
        absolute_epsilon=1e-6)
          .generalized_inverse_as_packed_u()
          .matrix_packed_u_as_symmetric())

else:
  import scitbx_matrix as matrix

  def generalized_inverse(m):
    return m.inverse()

def Xrot(E):
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

def Xtrans(r):
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

def T_as_X(T):
  """RBDA Eq. 2.28, p. 22:
Conversion of matrix.rt object T to spatial transform X.
"""
  return Xrot(T.r) * Xtrans(-T.r.transpose() * T.t)

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

def kinetic_energy(I_spatial, v_spatial):
  "RBDA Eq. 2.67, p. 35"
  return 0.5 * v_spatial.dot(I_spatial * v_spatial)

class system_model(object):
  "RBDA Tab. 4.3, p. 87"

  def __init__(O, bodies):
    "Stores bodies and computes Xtree (RBDA Fig. 4.7, p. 74)"
    O.bodies = bodies
    for B in bodies:
      if (B.parent == -1):
        Ttree = B.A.T0b
      else:
        Ttree = B.A.T0b * bodies[B.parent].A.Tb0
      B.Xtree = T_as_X(Ttree)

  def Xup(O):
    "RBDA Example 4.4, p. 80"
    return [B.J.Xj * B.Xtree for B in O.bodies]

  def spatial_velocities(O, Xup):
    "RBDA Example 4.4, p. 80"
    result = []
    if (Xup is None): Xup = O.Xup()
    for B,Xup_i in zip(O.bodies, O.Xup()):
      if (B.J.S is None): vJ = B.qd
      else:               vJ = B.J.S * B.qd
      if (B.parent == -1): result.append(vJ)
      else:                result.append(Xup_i * result[B.parent] + vJ)
    return result

  def e_kin(O, Xup=None):
    "RBDA Eq. 2.67, p. 35"
    result = 0
    for B,v in zip(O.bodies, O.spatial_velocities(Xup=Xup)):
      result += kinetic_energy(I_spatial=B.I, v_spatial=v)
    return result

  def accumulated_spatial_inertia(O, Xup=None):
    result = [B.I for B in O.bodies]
    if (Xup is None): Xup = O.Xup()
    for i in xrange(len(result)-1,-1,-1):
      B = O.bodies[i]
      if (B.parent != -1):
        result[B.parent] += Xup[i].transpose() * result[i] * Xup[i]
    return result

  def qd_e_kin_scales(O, Xup=None, e_kin_epsilon=1.e-12):
    result = []
    rap = result.append
    for B,asi in zip(O.bodies, O.accumulated_spatial_inertia(Xup=Xup)):
      S = B.J.S
      j_dof = B.J.degrees_of_freedom
      qd = [0] * j_dof
      for i in xrange(j_dof):
        qd[i] = 1
        vJ = matrix.col(qd)
        qd[i] = 0
        if (S is not None): vJ = S * vJ
        e_kin = kinetic_energy(I_spatial=asi, v_spatial=vJ)
        if (e_kin < e_kin_epsilon):
          rap(1)
        else:
          rap(1 / e_kin**0.5)
    return result

  def ID(O, qdd, f_ext=None, grav_accn=None):
    """RBDA Tab. 5.1, p. 96:
Inverse Dynamics of a kinematic tree via Recursive Newton-Euler Algorithm.
qdd is a vector of joint acceleration variables.
The return value (tau) is a vector of joint force variables.
f_ext specifies external forces acting on the bodies. If f_ext is None
then there are no external forces; otherwise, f_ext[i] is a spatial force
vector giving the force acting on body i, expressed in body i coordinates.
grav_accn is a 6D vector expressing the linear acceleration due to gravity.
    """

    Xup = O.Xup()
    v = O.spatial_velocities(Xup=Xup)
    a = [None] * len(Xup)
    f = [None] * len(Xup)
    for i in xrange(len(O.bodies)):
      B = O.bodies[i]
      if (B.J.S is None):
        vJ = B.qd
        aJ = qdd[i]
      else:
        vJ = B.J.S * B.qd
        aJ = B.J.S * qdd[i]
      if (B.parent == -1):
        a[i] = aJ
        if (grav_accn is not None):
          a[i] += Xup[i] * (-grav_accn)
      else:
        a[i] = Xup[i] * a[B.parent] + aJ + crm(v[i]) * vJ
      f[i] = B.I * a[i] + crf(v[i]) * B.I * v[i]
      if (f_ext is not None and f_ext[i] is not None):
        f[i] -= f_ext[i]

    tau = [None] * len(Xup)
    for i in xrange(len(Xup)-1,-1,-1):
      B = O.bodies[i]
      if (B.J.S is None): tau[i] = f[i]
      else:               tau[i] = B.J.S.transpose() * f[i]
      if (B.parent != -1):
        f[B.parent] += Xup[i].transpose() * f[i]

    return tau

  def ID0(O, f_ext):
    """
Simplified version of Inverse Dynamics via Recursive Newton-Euler Algorithm,
with all qd, qdd zero, but non-zero external forces.
    """
    Xup = O.Xup()
    f = [-e for e in f_ext]
    tau = [None] * len(f)
    for i in xrange(len(f)-1,-1,-1):
      B = O.bodies[i]
      if (B.J.S is None): tau[i] = f[i]
      else:               tau[i] = B.J.S.transpose() * f[i]
      if (B.parent != -1):
        f[B.parent] += Xup[i].transpose() * f[i]
    return tau

  def d_pot_d_q(O, f_ext):
    """
Gradients of potential energy (defined via f_ext) w.r.t. positional
coordinates q. Uses ID0().
    """
    return [B.J.tau_as_d_pot_d_q(tau=tau)
      for B,tau in zip(O.bodies, O.ID0(f_ext=f_ext))]

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

    Xup = O.Xup()
    v = O.spatial_velocities(Xup=Xup)
    c = [None] * len(Xup)
    IA = [None] * len(Xup)
    pA = [None] * len(Xup)
    for i in xrange(len(O.bodies)):
      B = O.bodies[i]
      if (B.J.S is None): vJ = B.qd
      else:               vJ = B.J.S * B.qd
      if (B.parent == -1): c[i] = matrix.col([0,0,0,0,0,0])
      else:                c[i] = crm(v[i]) * vJ
      IA[i] = B.I
      pA[i] = crf(v[i]) * B.I * v[i]
      if (f_ext is not None and f_ext[i] is not None):
        pA[i] -= f_ext[i]

    U = [None] * len(Xup)
    d_inv = [None] * len(Xup)
    u = [None] * len(Xup)
    for i in xrange(len(Xup)-1,-1,-1):
      B = O.bodies[i]
      if (B.J.S is None):
        U[i] = IA[i]
        d = U[i]
        if (tau is None or tau[i] is None):
          u[i] =        - pA[i]
        else:
          u[i] = tau[i] - pA[i]
      else:
        U[i] = IA[i] * B.J.S
        d = B.J.S.transpose() * U[i]
        if (tau is None or tau[i] is None):
          u[i] =        - B.J.S.transpose() * pA[i]
        else:
          u[i] = tau[i] - B.J.S.transpose() * pA[i]
      d_inv[i] = generalized_inverse(d)
      if (B.parent != -1):
        Ia = IA[i] - U[i] * d_inv[i] * U[i].transpose()
        pa = pA[i] + Ia*c[i] + U[i] * d_inv[i] * u[i]
        IA[B.parent] += Xup[i].transpose() * Ia * Xup[i]
        pA[B.parent] += Xup[i].transpose() * pa

    a = [None] * len(Xup)
    qdd = [None] * len(Xup)
    for i in xrange(len(O.bodies)):
      B = O.bodies[i]
      if (B.parent == -1):
        a[i] = c[i]
        if (grav_accn is not None):
          a[i] += Xup[i] * (-grav_accn)
      else:
        a[i] = Xup[i] * a[B.parent] + c[i]
      qdd[i] = d_inv[i] * (u[i] - U[i].transpose()*a[i])
      if (B.J.S is None):
        a[i] += qdd[i]
      else:
        a[i] += B.J.S * qdd[i]

    return qdd
