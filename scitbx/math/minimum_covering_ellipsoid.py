from __future__ import absolute_import, division, print_function
from libtbx import slots_getstate_setstate

class compute(slots_getstate_setstate):

  __slots__ = ["center", "radii", "rotation"]

  def __init__(O, points, epsilon=1e-2):
    """\
Computation of Minimum-Volume Covering Ellipsoid using the Khachiyan Algorithm.
Based on a Python implementation by Raj Rajashankar (ANL, Nov 2011),
which in turn was based on a Matlab script by Nima Moshtagh (2009,
http://stackoverflow.com/questions/1768197/bounding-ellipse/1768440#1768440).

Caveats:
  - center and radii are correct, but rotation may permute axes
  - scales with the square of the number of points
"""
    d = 3 # d is the dimension
    n = points.size()
    assert n > 0
    #
    from scitbx.array_family import flex
    p = points.as_double()
    p.reshape(flex.grid(n, 3))
    p = p.matrix_transpose()
    q = p.deep_copy()
    q.resize(flex.grid(4, n), 1)
    #
    u = flex.double(n, 1/n)
    umx = flex.double(flex.grid(n, n), 0)
    #
    err = epsilon + 1
    while (err > epsilon):
      umx.matrix_diagonal_set_in_place(u)
      x_inv = q.matrix_multiply(umx).matrix_multiply_transpose(
        q).matrix_inversion()
      m = q.matrix_transpose_multiply(
        x_inv).matrix_multiply(q).matrix_diagonal()
      j = flex.max_index(m)
      maximum = m[j]
      ascent_step_size = (maximum-d-1)/((d+1)*(maximum-1))
      new_u = (1 - ascent_step_size) * u
      new_u[j] += ascent_step_size
      err = flex.sum_sq(new_u - u)**0.5
      u = new_u
    #
    center = p.matrix_multiply(u)
    umx.matrix_diagonal_set_in_place(u)
    t1 = p.matrix_multiply(umx).matrix_multiply_transpose(p)
    t2 = center.matrix_outer_product(center)
    a = (t1 - t2).matrix_inversion() / d
    #
    import scitbx.linalg.svd
    svd = scitbx.linalg.svd.real(a, accumulate_u=False, accumulate_v=True)
    size = 1.0/flex.sqrt(svd.sigma)
    from scitbx import matrix
    O.center = matrix.col(center)
    O.radii = matrix.col(size)
    O.rotation = matrix.sqr(svd.v)
