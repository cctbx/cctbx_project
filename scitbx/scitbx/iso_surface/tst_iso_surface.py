from __future__ import division
from scitbx.array_family import flex
from scitbx import iso_surface
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
from math import sin, cos
from scitbx import matrix
import sys


class triangulation_test_case(object):

  def __init__(self, func, grid_size):
    """ Construct a map of func with the given grid_size """
    self.func = func
    self.grid_size = grid_size

    nx, ny, nz = grid_size
    self.grid_cell = 1/nx, 1/ny, 1/nz

    map = flex.double(flex.grid(grid_size))
    loop = flex.nested_loop(end=grid_size)
    while not loop.over():
      p = loop()
      map[p] = func((p[0]/nx, p[1]/ny, p[2]/nz))
      loop.incr()
    self.map = map

  def run(self, iso_level, verbose):
    """ Test triangulation of the iso-surface at the given iso-level """
    f = self.func

    # triangulation of the iso-surface of the map
    s = iso_surface.triangulation(self.map, iso_level, self.grid_cell)

    # the value of f on the vertices v shall not differ from iso_level by more
    # than ||1/2 f''(v).h|| where h=(dx,dy,dz)
    bad_vertices = []
    for i,v in enumerate(s.vertices):
      val = f(v)
      tol = f.second_order_error(v, self.grid_cell)
      if not approx_equal(val, iso_level, tol):
        bad_vertices.append((v,val))
    assert not bad_vertices, self.bad_vertices_err(bad_vertices, iso_level)

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
    self.degenerate_edges = degenerates

    vertices = {}
    edges = {}
    for v1,v2,v3 in s.triangles:
      vertices[v1] = vertices[v2] = vertices[v3] = 1
      for a,b in ((v1,v2), (v2,v3), (v3,v1)):
        if a < b: e = (a,b)
        else: e = (b,a)
        edges[e] = edges.setdefault(e,0) + 1
    assert len(vertices) == len(s.vertices)
    for e,p in edges.iteritems():
      assert p in (1,2)
      if not self.is_near_boundary(s.vertices[e[0]], s.vertices[e[1]]):
        assert p == 2

    # consistency check on the normals
    assert len(s.normals) == len(s.vertices)
    for v,n in zip(s.vertices, s.normals):
      v, n = matrix.col(v), matrix.col(n)
      assert approx_equal(abs(n), 1, 1e-12)
      outward = v + 0.01*n
      inward  = v - 0.01*n
      assert f(outward) > iso_level > f(inward)

  def is_near_boundary(self, v1, v2):
    delta = matrix.col(self.grid_cell)/2
    for c1, c2, m, eps in zip(v1, v2, self.grid_size, delta):
      if c1 != c2: continue
      if (   abs(c1) < eps       or abs(c1 - 1) < eps
          or abs(c1 - 1/m) < eps or abs(c1 - (1 - 1/m)) < eps ): return True
    return False

  def bad_vertices_err(self, bad_vertices, iso_level):
    msg = ["Vertices with value not close enough to threshold %f:" % iso_level]
    for vertex, val in bad_vertices:
      msg.append( "\t(%f, %f, %f) -> %f" % (v+(val,)) )
    return "\n".join(msg)


class elliptic(object):

  def __call__(self, p):
    x,y,z = p
    return x*x + 2*y*y + 3*z*z

  def second_order_error(self, p, h):
    return 0.5*abs(self(h))

class hyperbolic(object):

  def __call__(self, p):
    x,y,z = p
    return x*x - y*y - z*z

  def second_order_error(self, p, h):
    return 0.5*abs(self(h))

class sinusoidal(object):

  def __call__(self, p):
    x,y,z = p
    return sin(x*y + y*z + z*x)

  def second_order_error(self, p, h):
    x,y,z = p
    hx,hy,hz = h
    s = sin(x*y + y*z + z*x)
    c = cos(x*y + y*z + z*x)
    return 0.5*abs( ( -(y+z)**2 * hx*hx -(z+x)**2 * hy*hy -(x+y)**2 * hz*hz )*s
                    + 2*(c - (z+x)*(y+z)*s) * hx*hy
                    + 2*(c - (x+y)*(z+x)*s) * hy*hz
                    + 2*(c - (x+y)*(y+z)*s) * hz*hx
                  )

def run(args):
  verbose = "--verbose" in args
  grid_size = (50, 40, 30)

  test = triangulation_test_case(elliptic(), grid_size)
  test.run(iso_level=3, verbose=verbose)
  assert test.degenerate_edges == [(303, 418), (418, 407), (407, 303),
                                   (303, 407), (418, 303), (407, 418),
                                   (654, 830), (830, 822), (822, 654),
                                   (654, 822), (830, 654), (715, 909),
                                   (909, 890), (890, 715), (715, 890),
                                   (909, 715), (822, 830), (890, 909),
                                   (1898, 2185), (2185, 2177), (2177, 1898),
                                   (1898, 2177), (2185, 1898), (1966, 2278),
                                   (2278, 2261), (2261, 1966), (1966, 2261),
                                   (2278, 1966), (2177, 2185), (2261, 2278)]

  test = triangulation_test_case(hyperbolic(), grid_size)
  test.run(iso_level=3, verbose=verbose)
  assert test.degenerate_edges == []

  test = triangulation_test_case(sinusoidal(), grid_size)
  test.run(iso_level=3, verbose=verbose)
  assert test.degenerate_edges == []

  print format_cpu_times()

if __name__ == '__main__':
  run(sys.argv[1:])
