from __future__ import division
from scitbx.array_family import flex
from scitbx import iso_surface
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
from math import sin, cos
from scitbx import matrix
import sys
import time


class triangulation_test_case(object):

  def __init__(self, func, grid_size, lazy_normals, descending_normals):
    """ Construct a map of func with the given grid_size """
    self.func = func
    nx, ny, nz = self.grid_size = grid_size


    hx, hy, hz = self.grid_cell = 1/(nx-1), 1/(ny-1), 1/(nz-1)

    map = flex.double(flex.grid(grid_size))
    loop = flex.nested_loop(end=grid_size)
    while not loop.over():
      p = loop()
      map[p] = func((p[0]*hx, p[1]*hy, p[2]*hz))
      loop.incr()
    self.map = map
    self.lazy_normals = lazy_normals
    self.descending_normals = descending_normals

  def run(self, iso_level, verbose):
    """ Test triangulation of the iso-surface at the given iso-level """
    f = self.func

    # triangulation of the iso-surface of the map
    t0 = time.time()
    s = iso_surface.triangulation(
      self.map, iso_level,
      map_extent=(1,1,1),
      lazy_normals=self.lazy_normals,
      ascending_normal_direction = not self.descending_normals)
    t1 = time.time()
    if verbose:
      print "iso-surface triangulation per se: %f s" % (t1-t0)

    # make sure there is something to test!!
    assert s.vertices.size() > 0

    # the value of f on the vertices v shall not differ from iso_level by more
    # than ||1/2 f''(v).h|| where h=(dx,dy,dz)
    bad_vertices = []
    for i,v in enumerate(s.vertices):
      val = f(v)
      tol = f.second_order_error(v, self.grid_cell)
      if not abs(val - iso_level) < tol:
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

    # this maps each vertex to the set of its neighbours
    adjacencies = {}
    # triangle edges
    edges = {}
    for v1,v2,v3 in s.triangles:
      adjacencies.setdefault(v1,{}).update({v2:1, v3:1})
      adjacencies.setdefault(v2,{}).update({v1:1, v3:1})
      adjacencies.setdefault(v3,{}).update({v1:1, v2:1})
      for a,b in ((v1,v2), (v2,v3), (v3,v1)):
        if a < b: e = (a,b)
        else: e = (b,a)
        edges[e] = edges.setdefault(e,0) + 1
    assert len(adjacencies) == len(s.vertices)
    for e,p in edges.iteritems():
      assert p in (1,2)
      if not self.is_near_boundary(s.vertices[e[0]], s.vertices[e[1]]):
        assert p == 2

    # consistency check on the normals
    assert len(s.normals) == len(s.vertices)
    for i,v,n in zip(xrange(len(s.vertices)), s.vertices, s.normals):
      v, n = matrix.col(v), matrix.col(n)
      abs_n = abs(n)
      if abs_n == 0:
        for edge in degenerates:
          if i in edge: break
        else:
          raise "zero normal not on a vertex at the end of a degenerate edge"
      else:
        assert abs(abs(n) - 1) < 1e-12
        outward = v + 0.01*n
        inward  = v - 0.01*n
        if s.ascending_normal_direction:
          assert f(outward) > iso_level > f(inward), i
        else:
          assert f(outward) < iso_level < f(inward), i

  def is_near_boundary(self, v1, v2):
    delta = matrix.col(self.grid_cell)/2
    for c1, c2, m, eps in zip(v1, v2, self.grid_size, delta):
      if c1 != c2: continue
      if abs(c1) < eps or abs(c1 - 1) < eps: return True
    return False

  def bad_vertices_err(self, bad_vertices, iso_level):
    msg = ["Vertices with value not close enough to threshold %f:" % iso_level]
    for vertex, val in bad_vertices:
      msg.append( "\t(%f, %f, %f) -> %f" % (vertex+(val,)) )
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

  """ For this one, the iso-surface passes through points at corners of the
  map, e.g. (1, 1, 0). That makes it interesting for that corner vertex ends
  up being part of only one triangle which is degenerate and the normal
  associated to that vertex is therefore undefined """
  test = triangulation_test_case(elliptic(), grid_size,
                                 lazy_normals=False,
                                 descending_normals=False)
  test.run(iso_level=3, verbose=verbose)
  assert test.degenerate_edges == [(2973, 2912)]

  test = triangulation_test_case(elliptic(), grid_size,
                                 lazy_normals=True,
                                 descending_normals=False)
  test.run(iso_level=2.9, verbose=verbose)
  assert test.degenerate_edges == []

  test = triangulation_test_case(hyperbolic(), grid_size,
                                 lazy_normals=True,
                                 descending_normals=False)
  test.run(iso_level=0.2, verbose=verbose)
  assert test.degenerate_edges == []

  test = triangulation_test_case(sinusoidal(), grid_size,
                                 lazy_normals=False,
                                 descending_normals=True)
  test.run(iso_level=0.8, verbose=verbose)
  assert test.degenerate_edges == []

  print format_cpu_times()

if __name__ == '__main__':
  run(sys.argv[1:])
