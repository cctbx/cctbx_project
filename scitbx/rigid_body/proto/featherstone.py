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
from __future__ import absolute_import, division, print_function
from six.moves import range

try:
  import scitbx
except ImportError:
  scitbx = None

if (scitbx is not None):
  import scitbx.math
  from scitbx import matrix
  from libtbx.math_utils import ifloor

  def generalized_inverse(m):
    # assumption to achieve stability: order of magnitude of masses is around 1
    return matrix.sqr(
      scitbx.linalg.eigensystem.real_symmetric(
        m=m.as_flex_double_matrix(),
        relative_epsilon=1e-12,
        absolute_epsilon=1e-12)
          .generalized_inverse_as_packed_u()
          .matrix_packed_u_as_symmetric())

else:
  import scitbx_matrix as matrix

  def ifloor(x):
    def iround(x):
      if (x < 0): return int(x-0.5)
      return int(x+.5)
    return iround(math.floor(x))

  def generalized_inverse(m):
    return m.inverse()

import math

class InfType(object): pass
Inf = InfType()

def Xrotx(theta):
  """
% Xrotx  spatial coordinate transform (X-axis rotation).
% Xrotx(theta) calculates the coordinate transform matrix from A to B
% coordinates for spatial motion vectors, where coordinate frame B is
% rotated by an angle theta (radians) relative to frame A about their
% common X axis.
  """
  c = math.cos(theta)
  s = math.sin(theta)
  return matrix.sqr((
     1,  0,  0,  0,  0,  0,
     0,  c,  s,  0,  0,  0,
     0, -s,  c,  0,  0,  0,
     0,  0,  0,  1,  0,  0,
     0,  0,  0,  0,  c,  s,
     0,  0,  0,  0, -s,  c))

def Xroty(theta):
  """
% Xroty  spatial coordinate transform (Y-axis rotation).
% Xroty(theta) calculates the coordinate transform matrix from A to B
% coordinates for spatial motion vectors, where coordinate frame B is
% rotated by an angle theta (radians) relative to frame A about their
% common Y axis.
  """
  c = math.cos(theta)
  s = math.sin(theta)
  return matrix.sqr((
     c,  0, -s,  0,  0,  0,
     0,  1,  0,  0,  0,  0,
     s,  0,  c,  0,  0,  0,
     0,  0,  0,  c,  0, -s,
     0,  0,  0,  0,  1,  0,
     0,  0,  0,  s,  0,  c))

def Xrotz(theta):
  """
% Xrotz  spatial coordinate transform (Z-axis rotation).
% Xrotz(theta) calculates the coordinate transform matrix from A to B
% coordinates for spatial motion vectors, where coordinate frame B is
% rotated by an angle theta (radians) relative to frame A about their
% common Z axis.
  """
  c = math.cos(theta)
  s = math.sin(theta)
  return matrix.sqr((
     c,  s,  0,  0,  0,  0,
    -s,  c,  0,  0,  0,  0,
     0,  0,  1,  0,  0,  0,
     0,  0,  0,  c,  s,  0,
     0,  0,  0, -s,  c,  0,
     0,  0,  0,  0,  0,  1))

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

class autoTree(object):

  def __init__(self, nb, bf=1, skew=0, taper=1):
    """
% autoTree  Create System Models of Kinematic Trees
% autoTree(nb,bf,skew,taper) creates system models of kinematic trees
% having revolute joints.  nb and bf specify the number of bodies in the
% tree, and the branching factor, respectively.  The latter is the average
% number of children of a nonterminal node, and must be >=1.  bf=1 produces
% an unbranched tree; bf=2 produces a binary tree; and non-integer values
% produce trees in which the number of children alternates between
% floor(bf) and ceil(bf) in such a way that the average is bf.  Trees are
% constructed (and numbered) breadth-first.  Link i is a thin-walled
% cylindrical tube of length l(i), radius l(i)/20, and mass m(i), lying
% between 0 and l(i) on the x axis of its local coordinate system.  The
% values of l(i) and m(i) are determined by the tapering coefficient:
% l(i)=taper^(i-1) and m(i)=taper^(3*(i-1)).  Thus, if taper=1 then
% m(i)=l(i)=1 for all i.  The inboard joint axis of link i lies on the
% local z axis, and its outboard axis passes through the point (l(i),0,0)
% and is rotated about the x axis by an angle of skew radians relative to
% the inboard axis.  If the link has more than one outboard joint then they
% all have the same axis.  If skew=0 then the mechanism is planar.  The
% final one, two or three arguments can be omitted, in which case they
% assume default values of taper=1, skew=0 and bf=1.
    """
    self.NB = nb
    self.pitch = [0] * nb
    self.parent = [None] * nb
    self.Xtree = []
    self.I = []
    len_ = []
    for i in range(nb):
      self.parent[i] = ifloor((i-1+math.ceil(bf))/bf)-1
      if (self.parent[i] == -1):
        self.Xtree.append(Xtrans([0,0,0]))
      else:
        self.Xtree.append(Xrotx(skew) * Xtrans([len_[self.parent[i]],0,0]))
      len_.append(taper**i)
      mass = taper**(3*i)
      CoM = len_[i] * matrix.col([0.5,0,0])
      Icm = mass * len_[i]**2 * matrix.diag([0.0025,1.015/12,1.015/12])
      self.I.append(mcI(mass, CoM, Icm))

def jcalc(pitch, q):
  """
% jcalc  Calculate joint transform and motion subspace.
% [Xj,S]=jcalc(pitch,q) calculates the joint transform and motion subspace
% matrices for a revolute (pitch==0), prismatic (pitch==inf) or helical
% (pitch==any other value) joint.  For revolute and helical joints, q is
% the joint angle.  For prismatic joints, q is the linear displacement.
  """
  if (not isinstance(pitch, (int, float, InfType))):
    return pitch.Xj_S(q=q)
  if pitch == 0:                          # revolute joint
    Xj = Xrotz(q)
    S = matrix.col([0,0,1,0,0,0])
  elif pitch == Inf:                      # prismatic joint
    Xj = Xtrans([0,0,q])
    S = matrix.col([0,0,0,0,0,1])
  else:                                   # helical joint
    Xj = Xrotz(q) * Xtrans([0,0,q*pitch])
    S = matrix.col([0,0,1,0,0,pitch])
  return Xj, S

def grav_accn_as_a_grav(grav_accn):
  if grav_accn is None:
    return matrix.col([0,0,0,0,0,-9.81])
  grav_accn = list(grav_accn)
  assert len(grav_accn) == 3
  return matrix.col([0,0,0]+grav_accn)

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
% the values {} and [0,0,-9.81], respectively, if omitted.
  """

  a_grav = grav_accn_as_a_grav(grav_accn)

  S = [None] * model.NB
  Xup = [None] * model.NB
  v = [None] * model.NB
  a = [None] * model.NB
  f = [None] * model.NB
  for i in range(model.NB):
    XJ, S[i] = jcalc( model.pitch[i], q[i] )
    if (S[i] is None):
      vJ = qd[i]
      aJ = qdd[i]
    else:
      vJ = S[i]*qd[i]
      aJ = S[i]*qdd[i]
    Xup[i] = XJ * model.Xtree[i]
    if model.parent[i] == -1:
      v[i] = vJ
      a[i] = Xup[i] * -a_grav + aJ
    else:
      v[i] = Xup[i]*v[model.parent[i]] + vJ
      a[i] = Xup[i]*a[model.parent[i]] + aJ + crm(v[i])*vJ
    f[i] = model.I[i]*a[i] + crf(v[i])*model.I[i]*v[i]
    if (f_ext is not None and f_ext[i] is not None):
      f[i] = f[i] - f_ext[i]

  tau = [None] * model.NB
  for i in range(model.NB-1,-1,-1):
    if (S[i] is None):
      tau[i] = f[i]
    else:
      tau[i] = S[i].transpose() * f[i]
    if model.parent[i] != -1:
      f[model.parent[i]] = f[model.parent[i]] + Xup[i].transpose()*f[i]

  return tau

def FDab(model, q, qd, tau=None, f_ext=None, grav_accn=None, f_ext_in_ff=False):
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
% the values {} and [0,0,-9.81], respectively, if omitted.
  """

  a_grav = grav_accn_as_a_grav(grav_accn)

  S = [None] * model.NB
  Xup = [None] * model.NB
  X0 = [None] * model.NB
  v = [None] * model.NB
  c = [None] * model.NB
  IA = [None] * model.NB
  pA = [None] * model.NB
  for i in range(model.NB):
    XJ, S[i] = jcalc( model.pitch[i], q[i] )
    if (S[i] is None):
      vJ = qd[i]
    else:
      vJ = S[i]*qd[i]
    Xup[i] = XJ * model.Xtree[i]
    if model.parent[i] == -1:
      X0[i] = Xup[i]
      v[i] = vJ
      c[i] = matrix.col([0,0,0,0,0,0])
    else:
      X0[i] = Xup[i] * X0[model.parent[i]]
      v[i] = Xup[i]*v[model.parent[i]] + vJ
      c[i] = crm(v[i]) * vJ
    IA[i] = model.I[i]
    pA[i] = crf(v[i]) * model.I[i] * v[i]
    if (f_ext is not None and f_ext[i] is not None):
      if (not f_ext_in_ff):
        pA[i] = pA[i] - f_ext[i]
      else:
        pA[i] = pA[i] - X0[i].inverse().transpose() * f_ext[i]

  U = [None] * model.NB
  d_inv = [None] * model.NB
  u = [None] * model.NB
  for i in range(model.NB-1,-1,-1):
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
  for i in range(model.NB):
    if model.parent[i] == -1:
      a[i] = Xup[i] * -a_grav + c[i]
    else:
      a[i] = Xup[i] * a[model.parent[i]] + c[i]
    qdd[i] = d_inv[i] * (u[i] - U[i].transpose()*a[i])
    if (S[i] is None):
      a[i] = a[i] + qdd[i]
    else:
      a[i] = a[i] + S[i]*qdd[i]

  return qdd
