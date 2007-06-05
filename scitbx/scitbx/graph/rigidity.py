from __future__ import division

def gcd(a, b):
  ri = a
  if (ri < 0): ri = -ri
  rj = b
  if (rj != 0):
    while True:
      rk = ri % rj
      if (rk == 0):
        ri = rj
        break
      ri = rj % rk
      if (ri == 0):
        ri = rk
        break
      rj = rk % ri
      if (rj == 0):
        break
    if (ri < 0): ri = -ri
  return ri

def integer_row_echelon_form(m):
  free_vars = []
  n_rows = len(m)
  n_cols = len(m[0])
  piv_r = 0
  piv_c = 0
  while (piv_c < n_cols):
    best_r = None
    min_v = 0
    for i_row in xrange(piv_r, n_rows):
      v = abs(m[i_row][piv_c])
      if (v != 0 and (min_v == 0 or v < min_v)):
        min_v = v
        best_r = i_row
    if (best_r is not None):
      if (best_r != piv_r):
        m[piv_r], m[best_r] = m[best_r], m[piv_r]
      fp = m[piv_r][piv_c]
      for r in xrange(piv_r+1, n_rows):
        fr = m[r][piv_c]
        if (fr == 0): continue
        g = 0
        for c in xrange(piv_c, n_cols):
          m[r][c] = m[r][c] * fp - m[piv_r][c] * fr
          g = gcd(m[r][c], g)
        if (g > 1):
          for c in xrange(piv_c, n_cols):
            m[r][c] //= g
      piv_r += 1
    else:
      free_vars.append(piv_c)
    piv_c += 1
  return free_vars

def create_linearly_independent_integer_vertices(n_dim, n_vertices):
  # Idea due to Neil Sloane
  assert n_vertices != 0
  result = [range(2,2+n_dim)]
  while (len(result) != n_vertices):
    vertex = []
    for i,v0 in enumerate(result[0]):
      vertex.append(v0 * result[-1][i])
    result.append(vertex)
  return result

def construct_integer_rigidity_matrix(n_dim, n_vertices, edge_list):
  n_edges = len(edge_list)
  if (n_edges == 0): return None
  vertices = create_linearly_independent_integer_vertices(
    n_vertices=n_vertices, n_dim=n_dim)
  n_columns = n_vertices * n_dim
  result = []
  for i,j in edge_list:
    assert i != j
    row = [0] * n_columns
    dij = [vi-vj for vi,vj in zip(vertices[i], vertices[j])]
    def copy_to_row(c, sign):
      c *= n_dim
      for d in xrange(n_dim):
        row[c+d] = sign * dij[d]
    copy_to_row(c=i, sign= 1)
    copy_to_row(c=j, sign=-1)
    result.append(row)
  return result

def determine_degrees_of_freedom(n_dim, n_vertices, edge_list):
  if (n_vertices == 0): return None
  m = construct_integer_rigidity_matrix(
    n_dim=n_dim, n_vertices=n_vertices, edge_list=edge_list)
  if (m is None): return None
  free_vars = integer_row_echelon_form(m=m)
  return len(free_vars)

def example():
  # Donald Jacobs
  # Generic rigidity in three-dimensional bond-bending networks
  # Math. Gen. 31 (1998) 6653-6668)
  # Figure 1(d)
  n_vertices = 8
  edge_list = [
    (0,1), (0,2), (1,2),
    (0,3), (1,3), (2,3),
    (0,4), (1,4), (2,4),
    (5,6), (5,7), (6,7),
    (5,3), (6,3), (7,3),
    (5,4), (6,4), (7,4)]
  print "double banana 3D dof:", determine_degrees_of_freedom(
    n_dim=3, n_vertices=n_vertices, edge_list=edge_list)

if (__name__ == "__main__"):
  example()
