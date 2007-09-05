from __future__ import division
from scitbx.array_family import flex
from cctbx.maptbx import iso_surface
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
from math import sin, cos
from scitbx import matrix
import sys

def exercise(f, err, nx, ny, nz, iso_level, verbose):
  # construct a map of f
  dx, dy, dz = 1/nx, 1/ny, 1/nz
  map = flex.double(flex.grid((nx, ny, nz)))
  loop = flex.nested_loop(end=(nx,ny,nz))
  while not loop.over():
    p = loop()
    map[p] = f(p[0]/nx, p[1]/ny, p[2]/nz)
    loop.incr()

  # iso-surface of the map
  s = iso_surface(map, iso_level, (dx, dy, dz))

  # the value of f on the vertices v shall not differ from iso_level by more
  # than ||1/2 f''(v).h|| where h=(dx,dy,dz)
  bad_vertices = []
  for i,v in enumerate(s.vertices):
    val = f(v[0], v[1], v[2])
    tol = err(v[0], v[1], v[2], dx, dy, dz)
    if not approx_equal(val, iso_level, tol):
      bad_vertices.append((v,val))
  assert not bad_vertices, bad_vertices_err(bad_vertices)

  # consistency check on the triangulation
  # Note: a comprehensive check of the graph being triangulated is
  # not attempted
  degenerates = []
  for v1, v2, v3 in s.triangles:
    assert v1 != v2 and v2 != v3 and v3 != v1
    for a,b in ((v1,v2), (v2,v3), (v3,v1)):
      if s.vertices[a] == s.vertices[b]: degenerates.append((a,b))
  if verbose:
    if degenerates:
      print "Degenerate edges for the isosurface of %s:" % f
      print degenerates

  vertices = {}
  edges = {}
  for v1,v2,v3 in s.triangles:
    print v1,v2,v3
    vertices[v1] = vertices[v2] = vertices[v3] = 1
    for a,b in ((v1,v2), (v2,v3), (v3,v1)):
      if a < b: e = (a,b)
      else: e = (b,a)
      edges[e] = edges.setdefault(e,0) + 1
  assert len(vertices) == len(s.vertices)
  for e,p in edges.iteritems():
    assert p in (1,2)
    if not is_near_boundary(s.vertices[e[0]], s.vertices[e[1]], (nx,ny,nz)):
      assert p == 2

  # consistency check on the normals
  assert len(s.normals) == len(s.vertices)
  for v,n in zip(s.vertices, s.normals):
    v, n = matrix.col(v), matrix.col(n)
    assert approx_equal(abs(n), 1, 1e-12)
    outward = v + 0.01*n
    inward  = v - 0.01*n
    assert f(*outward) > iso_level > f(*inward)

def is_near_boundary(v1,v2,n):
  for c1, c2, m in zip(v1, v2, n):
    if c1 != c2: continue
    if (   approx_equal(c1, 0)   or approx_equal(c1, 1)
        or approx_equal(c1, 1/m) or approx_equal(c1, 1 - 1/m) ): return True
  return False

def bad_vertices_err(seq):
  msg = ["Vertices with value not close enough to threshold %f:" % iso_level]
  for vertex, val in bad_vertices:
    msg.append( "\t(%f, %f, %f) -> %f" % (v+(val,)) )
  return "\n".join(msg)

def elliptic(x,y,z):
  return x*x + 2*y*y + 3*z*z

def elliptic_error(x,y,z, hx,hy,hz):
  return 0.5*abs(elliptic(hx,hy,hz))

def hyperbolic(x,y,z):
  return x*x - y*y - z*z

def hyperbolic_error(x,y,z, hx,hy,hz):
  return 0.5*abs(hyperbolic(hx,hy,hz))

def sinusoidal(x,y,z):
  return sin(x*y + y*z + z*x)

def sinusoidal_error(x,y,z, hx,hy,hz):
  s = sin(x*y + y*z + z*x)
  c = cos(x*y + y*z + z*x)
  return 0.5*abs( ( -(y+z)**2 * hx*hx -(z+x)**2 * hy*hy -(x+y)**2 * hz*hz )*s
                  + 2*(c - (z+x)*(y+z)*s) * hx*hy
                  + 2*(c - (x+y)*(z+x)*s) * hy*hz
                  + 2*(c - (x+y)*(y+z)*s) * hz*hx
                )

def run(args):
  verbose = "--verbose" in args
  exercise(elliptic, elliptic_error, nx=50, ny=40, nz=30, iso_level=3,
    verbose=verbose)
  exercise(hyperbolic, hyperbolic_error, nx=50, ny=40, nz=30, iso_level=3,
    verbose=verbose)
  exercise(sinusoidal, sinusoidal_error, nx=50, ny=40, nz=30, iso_level=3,
    verbose=verbose)
  print format_cpu_times()

if __name__ == '__main__':
  run(sys.argv[1:])
