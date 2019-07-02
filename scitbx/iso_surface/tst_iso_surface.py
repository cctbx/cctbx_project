from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx import iso_surface
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
from math import sin
from scitbx import matrix
import sys
import time
import six
from six.moves import range
from six.moves import zip


class triangulation_test_case(object):

  def __init__(self, func, grid_size, periodic,
               lazy_normals, descending_normals):
    """ Construct a map of func with the given grid_size """
    self.func = func
    nx, ny, nz = self.grid_size = grid_size
    self.periodic = periodic

    if periodic:
      self.grid_cell = 1/nx, 1/ny, 1/nz
    else:
      self.grid_cell = 1/(nx-1), 1/(ny-1), 1/(nz-1)
    hx, hy, hz = self.grid_cell

    map = flex.double(flex.grid(grid_size))
    loop = flex.nested_loop(end=grid_size)
    while not loop.over():
      p = loop()
      map[p] = func((p[0]*hx, p[1]*hy, p[2]*hz))
      loop.incr()
    self.map = map
    self.lazy_normals = lazy_normals
    self.descending_normals = descending_normals

  def run(self, iso_level, from_here, to_there,
          verbose=0):
    """ Test triangulation of the iso-surface at the given iso-level """
    f = self.func

    if verbose:
      print("Testing %s" % f.__class__.__name__)

    # triangulation of the iso-surface of the map
    t0 = time.time()
    s = iso_surface.triangulation(
      self.map, iso_level,
      map_extent=(1,1,1),
      from_here=from_here, to_there=to_there,
      periodic=self.periodic,
      lazy_normals=self.lazy_normals,
      ascending_normal_direction = not self.descending_normals)
    self.triangulation = s
    t1 = time.time()
    if verbose:
      print("iso-surface triangulation per se: %f s" % (t1-t0))

    # make sure there is something to test!!
    assert s.vertices.size() > 0

    outside = [ v for v in s.vertices
                if not(s.from_here <= v <=  s.to_there) ]
    assert not outside

    # the value of f on the vertices v shall be close to iso_level
    deltas = flex.double()
    for v in s.vertices:
      val = f(v)
      deltas.append(abs((val - iso_level)/iso_level))
    assert (deltas > 0.07).count(True) == 0

    # consistency check on the triangulation
    degenerates = []
    for v1, v2, v3 in s.triangles:
      assert v1 != v2 and v2 != v3 and v3 != v1
      for a,b in ((v1,v2), (v2,v3), (v3,v1)):
        if s.vertices[a] == s.vertices[b]: degenerates.append((a,b))
    if verbose:
      if degenerates:
        print("Degenerate edges for the isosurface of %s:" % f)
        print(degenerates)
    self.degenerate_edges = degenerates

    # triangle edges and vertices
    edges = {}
    vertices = {}
    for v1,v2,v3 in s.triangles:
      vertices.update({v1:1, v2:1, v3:1})
      for a,b in ((v1,v2), (v2,v3), (v3,v1)):
        if a < b: e = (a,b)
        else: e = (b,a)
        edges[e] = edges.setdefault(e,0) + 1
    assert len(vertices) == len(s.vertices)
    missing = [ i for i in range(len(s.vertices)) if i not in vertices ]
    assert not missing, missing
    d = abs(matrix.col(self.grid_cell))
    bad_edge_multiplicities = []
    for e,p in six.iteritems(edges):
      v0, v1 = s.vertices[e[0]], s.vertices[e[1]]
      # conservative bound: edges shall be inscribed on voxel faces
      # this is just to catch vertex indexing errors resulting in edges
      # spanning several voxels
      d1 = abs(matrix.col(v0) - matrix.col(v1))
      assert d1 <= d or approx_equal(d1, d, eps=1e-3)
      assert p in (1,2)
      if not self.is_near_boundary(v0, v1):
        if p == 1:
          bad_edge_multiplicities.append((e[0], e[1]))
    assert not bad_edge_multiplicities

    # consistency check on the normals
    assert len(s.normals) == len(s.vertices)
    for i,v,n in zip(range(len(s.vertices)), s.vertices, s.normals):
      v, n = matrix.col(v), matrix.col(n)
      abs_n = abs(n)
      if abs_n == 0:
        for edge in degenerates:
          if i in edge: break
        else:
          raise RuntimeError("zero normal not on a vertex at the end of a degenerate edge")
      else:
        assert abs(abs(n) - 1) < 1e-12
        outward = v + 0.05*n
        inward  = v - 0.05*n
        if s.ascending_normal_direction:
          assert f(outward) > iso_level > f(inward), i
        else:
          assert f(outward) < iso_level < f(inward), i

  def is_near_boundary(self, v1, v2):
    delta = matrix.col(self.grid_cell)
    lower, higher = self.triangulation.bounds
    for c1, c2, l, h, eps in zip(v1, v2, lower, higher, delta):
      if c1 != c2: continue
      eps += 1e-15
      if abs(c1 - l) < eps or abs(c1 - h) < eps: return True
    return False


class elliptic(object):

  def __call__(self, p):
    x,y,z = p
    return x*x + 2*y*y + 3*z*z

class hyperbolic(object):

  def __call__(self, p):
    x,y,z = p
    return x*x - y*y - z*z

class sinusoidal(object):

  def __call__(self, p):
    x,y,z = p
    return sin(x*y + y*z + z*x)

class periodic(object):

  def __call__(self, p):
    x,y,z = p
    x %= 1
    y %= 1
    z %= 1
    return 8/3*(x**2*(1-x)**2 + 2*y**2*(1-y)**2 + 3*z**2*(1-z)**2)


def run(args):
  verbose = "--verbose" in args
  grid_size = (50, 40, 30)

  test = triangulation_test_case(periodic(), grid_size, periodic=True,
                                 lazy_normals=False,
                                 descending_normals=False)
  test.run(iso_level=0.194,
           from_here=(-0.5, -0.5, -0.5), to_there=(1.5, 1.5, 1.5),
           verbose=verbose)

  """ For this one, the iso-surface passes through points at corners of the
  map, e.g. (1, 1, 0). That makes it interesting for that corner vertex ends
  up being part of only one triangle which is degenerate and the normal
  associated to that vertex is therefore undefined """
  test = triangulation_test_case(elliptic(), grid_size, periodic=False,
                                 lazy_normals=False,
                                 descending_normals=False)
  test.run(iso_level=3, from_here=None, to_there=None, verbose=verbose)
  assert test.degenerate_edges == [(2973, 2912)]
  test.run(iso_level=1.3,
           from_here=(0.3, 0.2, 0.4), to_there=(0.7, 0.8, 0.6),
           verbose=verbose)
  assert test.degenerate_edges == []
  test.run(iso_level=1.4,
           from_here=(-0.3, 0.2, 0.4), to_there=(0.7, 0.8, 1.6),
           verbose=verbose)
  assert test.degenerate_edges == []

  test = triangulation_test_case(elliptic(), grid_size, periodic=False,
                                 lazy_normals=True,
                                 descending_normals=False)
  test.run(iso_level=0.8, from_here=None, to_there=None, verbose=verbose)
  assert test.degenerate_edges == []

  test = triangulation_test_case(hyperbolic(), grid_size, periodic=False,
                                 lazy_normals=True,
                                 descending_normals=False)
  test.run(iso_level=0.2, from_here=None, to_there=None, verbose=verbose)
  assert test.degenerate_edges == []

  test = triangulation_test_case(sinusoidal(), grid_size, periodic=False,
                                 lazy_normals=False,
                                 descending_normals=True)
  test.run(iso_level=0.8, from_here=None, to_there=None, verbose=verbose)
  assert test.degenerate_edges == []

  print(format_cpu_times())

if __name__ == '__main__':
  run(sys.argv[1:])
