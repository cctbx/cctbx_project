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

try: import scitbx
except ImportError: import scitbx_matrix as matrix
else: from scitbx import matrix

def xrot(e):
  """RBDA Tab. 2.2, p. 23:
Spatial coordinate transform (rotation around origin).
Calculates the coordinate transform matrix from A to B coordinates
for spatial motion vectors, in which frame B is rotated relative to
frame A.
  """
  a,b,c,d,e,f,g,h,i = e
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

def mci(m, c, i):
  """RBDA Eq. 2.63, p. 33:
Spatial rigid-body inertia from mass, CoM and rotational inertia.
Calculates the spatial inertia matrix of a rigid body from its
mass, centre of mass (3D vector) and rotational inertia (3x3 matrix)
about its centre of mass.
  """
  cx = matrix.cross_product_matrix(c)
  return matrix.sqr((
    i + m*cx*cx.transpose(), m*cx,
    m*cx.transpose(), m*matrix.identity(3))).resolve_partitions()

def kinetic_energy(i_spatial, v_spatial):
  "RBDA Eq. 2.67, p. 35"
  return 0.5 * v_spatial.dot(i_spatial * v_spatial)
