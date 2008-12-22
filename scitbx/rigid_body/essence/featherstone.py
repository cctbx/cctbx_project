"""\
Python version of a subset of Roy Featherstone's spatial_v1 matlab code:

  http://axiom.anu.edu.au/~roy/spatial/

  Version 1: January 2008 (latest bug fix: 7 October 2008)

The subset of converted files covers all dependencies of:
  ID.m
  FDab.m

The original matlab comments are preserved as Python docstrings.

See also:
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
        relative_epsilon=1e-12,
        absolute_epsilon=1e-12)
          .generalized_inverse_as_packed_u()
          .matrix_packed_u_as_symmetric())

else:
  import scitbx_matrix as matrix

  def generalized_inverse(m):
    return m.inverse()

import math

class InfType(object): pass
Inf = InfType()

def Xrot(E):
  """
  Featherstone (2007) Tab. 2.2
  Added in Python version.
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
  """
% Xtrans  spatial coordinate transform (translation of origin).
% Xtrans(r) calculates the coordinate transform matrix from A to B
% coordinates for spatial motion vectors, in which frame B is translated by
% an amount r (3D vector) relative to frame A.
  """
  r1,r2,r3 = r
  return matrix.sqr((
      1,   0,   0, 0, 0, 0,
      0,   1,   0, 0, 0, 0,
      0,   0,   1, 0, 0, 0,
      0,  r3, -r2, 1, 0, 0,
    -r3,   0,  r1, 0, 1, 0,
     r2, -r1,   0, 0, 0, 1))

def crm(v):
  """
% crm  spatial cross-product operator (motion).
% crm(v) calculates the 6x6 matrix such that the expression crm(v)*m is the
% cross product of the spatial motion vectors v and m.
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
  """
% crf  spatial cross-product operator (force).
% crf(v) calculates the 6x6 matrix such that the expression crf(v)*f is the
% cross product of the spatial motion vector v with the spatial force
% vector f.
  """
  return -crm(v).transpose()

def mcI(m, c, I):
  """
% mcI  spatial rigid-body inertia from mass, CoM and rotational inertia.
% mcI(m,c,I) calculates the spatial inertia matrix of a rigid body from its
% mass, centre of mass (3D vector) and rotational inertia (3x3 matrix)
% about its centre of mass.
  """
  c1,c2,c3 = c
  C = matrix.sqr((
      0, -c3,  c2,
     c3,   0, -c1,
    -c2,  c1,  0))
  return matrix.sqr((
    I + m*C*C.transpose(), m*C,
    m*C.transpose(), m*matrix.identity(3))).resolve_partitions()

def ID(model, q, qd, qdd, f_ext=None, grav_accn=None):
  """
% ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
% ID(model,q,qd,qdd,f_ext,grav_accn) calculates the inverse dynamics of a
% kinematic tree via the recursive Newton-Euler algorithm.  q, qd and qdd
% are vectors of joint position, velocity and acceleration variables; and
% the return value is a vector of joint force variables.  f_ext is a cell
% array specifying external forces acting on the bodies.  If f_ext == {}
% then there are no external forces; otherwise, f_ext{i} is a spatial force
% vector giving the force acting on body i, expressed in body i
% coordinates.  Empty cells in f_ext are interpreted as zero forces.
% grav_accn is a 3D vector expressing the linear acceleration due to
% gravity.  The arguments f_ext and grav_accn are optional, and default to
% the values {} and [0,0,0], respectively, if omitted.
  """

  S = [None] * model.NB
  Xup = [None] * model.NB
  v = [None] * model.NB
  a = [None] * model.NB
  f = [None] * model.NB
  for i in xrange(model.NB):
    XJ, S[i] = model.pitch[i].Xj_S(q=q[i])
    if (S[i] is None):
      vJ = qd[i]
      aJ = qdd[i]
    else:
      vJ = S[i]*qd[i]
      aJ = S[i]*qdd[i]
    Xup[i] = XJ * model.Xtree[i]
    if model.parent[i] == -1:
      v[i] = vJ
      a[i] = aJ
      if (grav_accn is not None):
        a[i] += Xup[i] * -grav_accn
    else:
      v[i] = Xup[i]*v[model.parent[i]] + vJ
      a[i] = Xup[i]*a[model.parent[i]] + aJ + crm(v[i])*vJ
    f[i] = model.I[i]*a[i] + crf(v[i])*model.I[i]*v[i]
    if (f_ext is not None and f_ext[i] is not None):
      f[i] = f[i] - f_ext[i]

  tau = [None] * model.NB
  for i in xrange(model.NB-1,-1,-1):
    if (S[i] is None):
      tau[i] = f[i]
    else:
      tau[i] = S[i].transpose() * f[i]
    if model.parent[i] != -1:
      f[model.parent[i]] = f[model.parent[i]] + Xup[i].transpose()*f[i]

  return tau

def FDab(model, q, qd, tau=None, f_ext=None, grav_accn=None):
  """
% FDab  Forward Dynamics via Articulated-Body Algorithm
% FDab(model,q,qd,tau,f_ext,grav_accn) calculates the forward dynamics of a
% kinematic tree via the articulated-body algorithm.  q, qd and tau are
% vectors of joint position, velocity and force variables; and the return
% value is a vector of joint acceleration variables.  f_ext is a cell array
% specifying external forces acting on the bodies.  If f_ext == {} then
% there are no external forces; otherwise, f_ext{i} is a spatial force
% vector giving the force acting on body i, expressed in body i
% coordinates.  Empty cells in f_ext are interpreted as zero forces.
% grav_accn is a 3D vector expressing the linear acceleration due to
% gravity.  The arguments f_ext and grav_accn are optional, and default to
% the values {} and [0,0,0], respectively, if omitted.
  """

  S = [None] * model.NB
  Xup = [None] * model.NB
  v = [None] * model.NB
  c = [None] * model.NB
  IA = [None] * model.NB
  pA = [None] * model.NB
  for i in xrange(model.NB):
    XJ, S[i] = model.pitch[i].Xj_S(q=q[i])
    if (S[i] is None):
      vJ = qd[i]
    else:
      vJ = S[i]*qd[i]
    Xup[i] = XJ * model.Xtree[i]
    if model.parent[i] == -1:
      v[i] = vJ
      c[i] = matrix.col([0,0,0,0,0,0])
    else:
      v[i] = Xup[i]*v[model.parent[i]] + vJ
      c[i] = crm(v[i]) * vJ
    IA[i] = model.I[i]
    pA[i] = crf(v[i]) * model.I[i] * v[i]
    if (f_ext is not None and f_ext[i] is not None):
      pA[i] = pA[i] - f_ext[i]

  U = [None] * model.NB
  d_inv = [None] * model.NB
  u = [None] * model.NB
  for i in xrange(model.NB-1,-1,-1):
    if (S[i] is None):
      U[i] = IA[i]
      d = U[i]
      if (tau is None or tau[i] is None):
        u[i] =        - pA[i]
      else:
        u[i] = tau[i] - pA[i]
    else:
      U[i] = IA[i] * S[i]
      d = S[i].transpose() * U[i]
      if (tau is None or tau[i] is None):
        u[i] =        - S[i].transpose()*pA[i]
      else:
        u[i] = tau[i] - S[i].transpose()*pA[i]
    d_inv[i] = generalized_inverse(d)
    if model.parent[i] != -1:
      Ia = IA[i] - U[i] * d_inv[i] * U[i].transpose()
      pa = pA[i] + Ia*c[i] + U[i] * d_inv[i] * u[i]
      IA[model.parent[i]] = IA[model.parent[i]] \
                          + Xup[i].transpose() * Ia * Xup[i]
      pA[model.parent[i]] = pA[model.parent[i]] \
                          + Xup[i].transpose() * pa

  a = [None] * model.NB
  qdd = [None] * model.NB
  for i in xrange(model.NB):
    if model.parent[i] == -1:
      a[i] = c[i]
      if (grav_accn is not None):
        a[i] += Xup[i] * -grav_accn
    else:
      a[i] = Xup[i] * a[model.parent[i]] + c[i]
    qdd[i] = d_inv[i] * (u[i] - U[i].transpose()*a[i])
    if (S[i] is None):
      a[i] = a[i] + qdd[i]
    else:
      a[i] = a[i] + S[i]*qdd[i]

  return qdd
