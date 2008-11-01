"""\
Python version of a subset of Roy Featherstone's spatial_v1 matlab code:

  http://axiom.anu.edu.au/~roy/spatial/

  Version 1: January 2008 (latest bug fix: 7 October 2008)

The subset of converted files covers all dependencies of:
  ID.m
  FDab.m
  IDf.m
  FDf.m

The original matlab comments are preserved as Python docstrings.

See also:
  Rigid Body Dynamics Algorithms.
  R. Featherstone,
  Springer, New York, 2007.
  ISBN-10: 0387743146
"""

try:
  import scitbx
except ImportError:
  scitbx = None

if (scitbx is not None):
  from scitbx import matrix
  from libtbx.math_utils import ifloor
else:
  import matrix

  def ifloor(x):
    def iround(x):
      if (x < 0): return int(x-0.5)
      return int(x+.5)
    return iround(math.floor(x))

import math

class InfType(object): pass
Inf = InfType()

def mldivide(A, B):
  "http://www.mathworks.com/access/helpdesk/help/techdoc/ref/mldivide.html"
  return A.inverse() * B

def mrdivide(B, A):
  "http://www.mathworks.com/access/helpdesk/help/techdoc/ref/mrdivide.html"
  return (A.transpose().inverse() * B.transpose()).transpose()

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
    for i in xrange(nb):
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

class floatbase(object):

  def __init__(self, model):
    """
% floatbase  construct the floating-base equivalent of a fixed-base model
% floatbase(model)  converts a fixed-base kinematic tree to a floating-base
% kinematic tree as follows: old body 1 becomes new body 6, and is regarded
% as the floating base; old joint 1 is discarded; six new joints are added
% (three prismatic and three revolute, in that order, arranged along the x,
% y and z axes, in that order); and five new zero-mass bodies are added
% (numbered 1 to 5), to connect the new joints.  All other bodies and
% joints in the given model are preserved, but their identification numbers
% are incremented by 5.  The end result is that body 6 has a full 6 degrees
% of motion freedom relative to the fixed base, with joint variables 1 to 6
% serving as a set of 3 Cartesian coordinates and 3 rotation angles (about
% x, y and z axes) specifying the position and orientation of body 6's
% coordinate frame relative to the base coordinate frame.  CAUTION: A
% singularity occurs when q(5)=+-pi/2, and these orientations must be
% avoided.  Note: this function requires that old joint 1 is the only joint
% connected to the fixed base, and that Xlink{1} is the identity.
    """
    assert model.parent[0] == -1
    for ip in model.parent[1:]:
      if (ip == -1):
        raise RuntimeError("only one connection to a fixed base allowed")
    if ((model.Xtree[0] - Xtrans([0,0,0])).dot() > 1.e-8):
      raise RuntimeError("model.Xtree[0] not identity")

    self.NB = model.NB + 5
    self.pitch = [Inf,Inf,Inf,0,0,0] + model.pitch[1:]
    self.parent = [-1,0,1,2,3,4] + [i+5 for i in model.parent[1:]]
    pi = math.pi
    self.Xtree = [
      Xroty(pi/2),
      Xrotx(-pi/2) * Xroty(-pi/2),
      Xrotx(pi/2),
      Xroty(pi/2),
      Xrotx(-pi/2) * Xroty(-pi/2),
      Xrotx(pi/2)] + model.Xtree[1:]
    zeros3 = matrix.diag([0,0,0])
    self.I = [mcI(0, [0,0,0], zeros3) for i in xrange(5)] + model.I

def jcalc(pitch, q):
  """
% jcalc  Calculate joint transform and motion subspace.
% [Xj,S]=jcalc(pitch,q) calculates the joint transform and motion subspace
% matrices for a revolute (pitch==0), prismatic (pitch==inf) or helical
% (pitch==any other value) joint.  For revolute and helical joints, q is
% the joint angle.  For prismatic joints, q is the linear displacement.
  """
  if (not isinstance(pitch, (int, float, InfType))):
    return pitch(q)
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
  for i in xrange(model.NB):
    XJ, S[i] = jcalc( model.pitch[i], q[i] )
    vJ = S[i]*qd[i]
    Xup[i] = XJ * model.Xtree[i]
    if model.parent[i] == -1:
      v[i] = vJ
      a[i] = Xup[i] * -a_grav + S[i]*qdd[i]
    else:
      v[i] = Xup[i]*v[model.parent[i]] + vJ
      a[i] = Xup[i]*a[model.parent[i]] + S[i]*qdd[i] + crm(v[i])*vJ
    f[i] = model.I[i]*a[i] + crf(v[i])*model.I[i]*v[i]
    if (f_ext is not None and f_ext[i] is not None):
      f[i] = f[i] - f_ext[i]

  tau = [None] * model.NB
  for i in xrange(model.NB-1,-1,-1):
    tau[i] = S[i].transpose() * f[i]
    if model.parent[i] != -1:
      f[model.parent[i]] = f[model.parent[i]] + Xup[i].transpose()*f[i]

  return tau

def FDab(model, q, qd, tau, f_ext=None, grav_accn=None):
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
  v = [None] * model.NB
  c = [None] * model.NB
  IA = [None] * model.NB
  pA = [None] * model.NB
  for i in xrange(model.NB):
    XJ, S[i] = jcalc( model.pitch[i], q[i] )
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
  d = [None] * model.NB
  u = [None] * model.NB
  for i in xrange(model.NB-1,-1,-1):
    U[i] = IA[i] * S[i]
    d[i] = S[i].transpose() * U[i]
    u[i] = tau[i] - S[i].transpose()*pA[i]
    if model.parent[i] != -1:
      Ia = IA[i] - mrdivide(U[i],d[i])*U[i].transpose()
      pa = pA[i] + Ia*c[i] + mrdivide(U[i] * u[i],d[i])
      IA[model.parent[i]] = IA[model.parent[i]] \
                          + Xup[i].transpose() * Ia * Xup[i]
      pA[model.parent[i]] = pA[model.parent[i]] \
                          + Xup[i].transpose() * pa

  a = [None] * model.NB
  qdd = [None] * model.NB
  for i in xrange(model.NB):
    if model.parent[i] == -1:
      a[i] = Xup[i] * -a_grav + c[i]
    else:
      a[i] = Xup[i] * a[model.parent[i]] + c[i]
    qdd[i] = mldivide(d[i], u[i] - U[i].transpose()*a[i])
    a[i] = a[i] + S[i]*qdd[i]

  return qdd

def IDf(model, Xfb, vfb, q, qd, qdd, f_ext=None, grav_accn=None):
  """
% IDf  Floating-Base Inverse Dynamics
% [afb,tau]=IDf(model,Xfb,vfb,q,qd,qdd,f_ext,grav_accn) calculates the
% inverse dynamics of a floating-base kinematic tree, such as that created
% by the function floatbase (i.e., body 6 in the system model is the
% floating base).  Xfb is the coordinate transform from fixed to floating
% base coordinates; vfb is the spatial velocity of the floating base,
% expressed in fixed-base coordinates; and q, qd and qdd contain the
% position, velocity and acceleration variables for the real joints in the
% system (i.e., joints 7 onwards in the system model).  The return values
% are the spatial acceleration of the floating base, expressed in
% fixed-base coordinates, and the joint force vector for the real joints.
% f_ext is a cell array specifying external forces acting on the bodies.
% If f_ext == {} then there are no external forces; otherwise, f_ext{i} is
% a spatial force vector giving the force acting on body i, expressed in
% body i coordinates.  Thus, f_ext{6} is the force acting on the floating
% base, and f_ext{1} to f_ext{5} are ignored.  Empty cells in f_ext are
% interpreted as zero forces.  grav_accn is a 3D vector expressing the
% linear acceleration due to gravity.  The arguments f_ext and grav_accn
% are optional, and default to the values {} and [0,0,-9.81], respectively,
% if omitted.  (Caution: vfb is expressed in fixed-base coordinates, but
% f_ext{6} is expressed in floating-base coordinates.)
  """

  a_grav = grav_accn_as_a_grav(grav_accn)

  vfb = Xfb * vfb
  afb = Xfb * -a_grav                     # fictitious base acc for calc pC
  # alternative: afb = zeros(6,1);

  NBR = model.NB - 6                      # NB & parent array for Rest of model
  parentR = [i-6 for i in model.parent[6:]]

  S = [None] * NBR
  Xup = [None] * NBR
  v = [None] * NBR
  a = [None] * NBR
  IC = [None] * NBR
  pC = [None] * NBR
  for i in xrange(NBR):
    XJ, S[i] = jcalc( model.pitch[i+6], q[i] )
    vJ = S[i]*qd[i]
    Xup[i] = XJ * model.Xtree[i+6]
    if parentR[i] == -1:
      v[i] = Xup[i]*vfb + vJ
      a[i] = Xup[i]*afb + S[i]*qdd[i] + crm(v[i])*vJ
    else:
      v[i] = Xup[i]*v[parentR[i]] + vJ
      a[i] = Xup[i]*a[parentR[i]] + S[i]*qdd[i] + crm(v[i])*vJ
    IC[i] = model.I[i+6]
    pC[i] = IC[i]*a[i] + crf(v[i])*IC[i]*v[i]
    if (f_ext is not None and f_ext[i+6] is not None):
      pC[i] = pC[i] - f_ext[i+6]

  ICfb = model.I[5]
  pCfb = ICfb*afb + crf(vfb)*ICfb*vfb
  if (f_ext is not None and f_ext[5] is not None):
    pCfb = pCfb - f_ext[5]

  for i in xrange(NBR-1,-1,-1):
    if parentR[i] == -1:
      ICfb = ICfb + Xup[i].transpose()*IC[i]*Xup[i]
      pCfb = pCfb + Xup[i].transpose()*pC[i]
    else:
      IC[parentR[i]] = IC[parentR[i]] + Xup[i].transpose()*IC[i]*Xup[i]
      pC[parentR[i]] = pC[parentR[i]] + Xup[i].transpose()*pC[i]

  afb = - mldivide(ICfb, pCfb)       # true base acceleration

  tau = [None] * NBR
  for i in xrange(NBR):
    if parentR[i] == -1:
      a[i] = Xup[i]*afb
    else:
      a[i] = Xup[i]*a[parentR[i]]
    tau[i] = S[i].transpose()*(IC[i]*a[i] + pC[i])

  afb = Xfb.inverse() * afb
  # alternative: afb = inv(Xfb) * afb + a_grav;

  return afb, tau

def FDf(model, Xfb, vfb, q, qd, tau, f_ext=None, grav_accn=None):
  """
% FDf  Floating-Base Forward Dynamics via Articulated-Body Algorithm
% [afb,qdd]=FDf(model,Xfb,vfb,q,qd,tau,f_ext,grav_accn) calculates the
% forward dynamics of a floating-base kinematic tree, such as that created
% by the function floatbase (i.e., body 6 in the system model is the
% floating base), via the articulated-body algorithm.  Xfb is the
% coordinate transform from fixed to floating base coordinates; vfb is the
% spatial velocity of the floating base, expressed in fixed-base
% coordinates; and q, qd and tau contain the position, velocity and force
% variables for the real joints in the system (i.e., joints 7 onwards in
% the system model).  The return values are the spatial acceleration of the
% floating base, expressed in fixed-base coordinates, and the joint
% acceleration vector for the real joints.  f_ext is a cell array
% specifying external forces acting on the bodies.  If f_ext == {} then
% there are no external forces; otherwise, f_ext{i} is a spatial force
% vector giving the force acting on body i, expressed in body i
% coordinates.  Thus, f_ext{6} is the force acting on the floating base,
% and f_ext{1} to f_ext{5} are ignored.  Empty cells in f_ext are
% interpreted as zero forces.  grav_accn is a 3D vector expressing the
% linear acceleration due to gravity.  The arguments f_ext and grav_accn
% are optional, and default to the values {} and [0,0,-9.81], respectively,
% if omitted.  (Caution: vfb is expressed in fixed-base coordinates, but
% f_ext{6} is expressed in floating-base coordinates.)
  """

  a_grav = grav_accn_as_a_grav(grav_accn)

  vfb = Xfb * vfb

  NBR = model.NB - 6                      # NB & parent array for Rest of model
  parentR = [i-6 for i in model.parent[6:]]

  S = [None] * NBR
  Xup = [None] * NBR
  v = [None] * NBR
  c = [None] * NBR
  IA = [None] * NBR
  pA = [None] * NBR
  for i in xrange(NBR):
    XJ, S[i] = jcalc(model.pitch[i+6], q[i])
    vJ = S[i]*qd[i]
    Xup[i] = XJ * model.Xtree[i+6]
    if parentR[i] == -1:
      v[i] = Xup[i]*vfb + vJ
    else:
      v[i] = Xup[i]*v[parentR[i]] + vJ
    c[i] = crm(v[i]) * vJ
    IA[i] = model.I[i+6]
    pA[i] = crf(v[i]) * IA[i] * v[i]
    if (f_ext is not None and f_ext[i+6] is not None):
      pA[i] = pA[i] - f_ext[i+6]

  IAfb = model.I[5]
  pAfb = crf(vfb) * IAfb * vfb
  if (f_ext is not None and f_ext[5] is not None):
    pAfb = pAfb - f_ext[5]

  U = [None] * NBR
  d = [None] * NBR
  u = [None] * NBR
  for i in xrange(NBR-1,-1,-1):
    U[i] = IA[i] * S[i]
    d[i] = S[i].transpose() * U[i]
    u[i] = tau[i] - S[i].transpose()*pA[i]
    Ia = IA[i] - mrdivide(U[i],d[i])*U[i].transpose()
    pa = pA[i] + Ia*c[i] + U[i] * mrdivide(u[i],d[i])
    if parentR[i] == -1:
      IAfb = IAfb + Xup[i].transpose() * Ia * Xup[i]
      pAfb = pAfb + Xup[i].transpose() * pa
    else:
      IA[parentR[i]] = IA[parentR[i]] + Xup[i].transpose() * Ia * Xup[i]
      pA[parentR[i]] = pA[parentR[i]] + Xup[i].transpose() * pa

  afb = - mldivide(IAfb, pAfb)       # floating base accn without gravity

  a = [None] * NBR
  qdd = [None] * NBR
  for i in xrange(NBR):
    if parentR[i] == -1:
      a[i] = Xup[i] * afb + c[i]
    else:
      a[i] = Xup[i] * a[parentR[i]] + c[i]
    qdd[i] = mrdivide((u[i] - U[i].transpose()*a[i]), d[i])
    a[i] = a[i] + S[i]*qdd[i]

  afb = mldivide(Xfb, afb) + a_grav  # true flt base accn in ref coords

  return afb, qdd
