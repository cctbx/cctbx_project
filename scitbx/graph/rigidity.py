from __future__ import absolute_import, division, print_function
import random
from six.moves import range
from six.moves import zip

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
    for i_row in range(piv_r, n_rows):
      v = abs(m[i_row][piv_c])
      if (v != 0 and (min_v == 0 or v < min_v)):
        min_v = v
        best_r = i_row
    if (best_r is not None):
      if (best_r != piv_r):
        m[piv_r], m[best_r] = m[best_r], m[piv_r]
      fp = m[piv_r][piv_c]
      for r in range(piv_r+1, n_rows):
        fr = m[r][piv_c]
        if (fr == 0): continue
        g = 0
        for c in range(piv_c, n_cols):
          m[r][c] = m[r][c] * fp - m[piv_r][c] * fr
          g = gcd(m[r][c], g)
        if (g > 1):
          for c in range(piv_c, n_cols):
            m[r][c] //= g
      piv_r += 1
    else:
      free_vars.append(piv_c)
    piv_c += 1
  return free_vars

def float_row_echelon_form(
      m,
      zero_pivot_tolerance=1.e-8,
      min_non_zero_pivot=1.e-2):
  free_vars = []
  n_rows = len(m)
  n_cols = len(m[0])
  piv_r = 0
  piv_c = 0
  while (piv_c < n_cols):
    # search for best pivot
    best_r = None
    max_v = zero_pivot_tolerance
    for i_row in range(piv_r, n_rows):
      v = abs(m[i_row][piv_c])
      if (v > max_v):
        max_v = v
        best_r = i_row
    if (best_r is not None):
      if (max_v < min_non_zero_pivot):
        return None
      if (best_r != piv_r):
        m[piv_r], m[best_r] = m[best_r], m[piv_r]
      fp = m[piv_r][piv_c]
      for r in range(piv_r+1, n_rows):
        fr = m[r][piv_c]
        if (fr == 0): continue
        for c in range(piv_c, n_cols):
          m[r][c] -= m[piv_r][c] * fr / fp
      piv_r += 1
    else:
      free_vars.append(piv_c)
    piv_c += 1
  return free_vars

def float_row_echelon_form_is_redundant(
      m,
      free_vars,
      addl_row,
      zero_pivot_tolerance=1.e-6,
      min_non_zero_pivot=1.e-3):
  n_rows = len(m)
  n_cols = len(m[0])
  assert len(addl_row) == n_cols
  free_flags = [False] * n_cols
  for c in free_vars:
    free_flags[c] = True
  def approx_zero(c):
    v = abs(addl_row[c])
    if (v < zero_pivot_tolerance):
      return True
    if (v < min_non_zero_pivot):
      return None
    return False
  piv_c = 0
  for piv_r in range(n_cols-len(free_vars)):
    while (free_flags[piv_c]):
      az = approx_zero(c=piv_c)
      if (not az): return az
      piv_c += 1
    fp = m[piv_r][piv_c]
    fr = addl_row[piv_c]
    if (fr != 0):
      for c in range(piv_c, n_cols):
        addl_row[c] -= m[piv_r][c] * fr / fp
    piv_c += 1
  for c in range(piv_c, n_cols):
    az = approx_zero(c=c)
    if (not az): return az
  return True

def float_row_echelon_form_back_substitution(m, free_vars, sol):
  n_rows = len(m)
  n_cols = len(m[0])
  assert len(sol) == n_cols
  free_flags = [False] * n_cols
  for c in free_vars:
    free_flags[c] = True
  piv_cols = []
  for c,f in enumerate(free_flags):
    if (not f): piv_cols.append(c)
  for r in range(len(piv_cols)-1,-1,-1):
    piv_c = piv_cols[r]
    s = 0
    for c in range(piv_c+1, n_cols):
      s += m[r][c] * sol[c]
    sol[piv_c] = -s / m[r][piv_c]

def create_fake_integer_vertices(n_dim, n_vertices):
  # Idea due to Neil Sloane
  assert n_vertices != 0
  v0 = list(range(2,2+n_dim))
  result = [v0]
  while (len(result) != n_vertices):
    vertex = []
    for i in range(n_dim):
      vertex.append(v0[i] * result[-1][i])
    result.append(vertex)
  # Shuffle coordinates. Required to obtain correct result for
  # "K6,6 minus six parallel edges" (Figure 3.23 of J.E. Graver,
  # Counting on Frameworks, 2001).
  for i in range(len(result)):
    vertex = result[i]
    j = i % n_dim
    result[i] = vertex[j:] + vertex[:j]
  return result

def create_fake_float_vertices(n_dim, n_vertices, max_coordinate=1.e4):
  assert n_vertices != 0
  result = []
  while (len(result) != n_vertices):
    vertex = []
    for i in range(n_dim):
      vertex.append(random.random()*max_coordinate)
    result.append(vertex)
  return result

def construct_numeric_rigidity_matrix(n_dim, vertices, edge_list):
  if (len(edge_list) == 0): return None
  n_columns = len(vertices) * n_dim
  result = []
  for i,j in edge_list:
    assert i != j
    row = [0] * n_columns
    dij = [vi-vj for vi,vj in zip(vertices[i], vertices[j])]
    def copy_to_row(c, sign):
      c *= n_dim
      for d in range(n_dim):
        row[c+d] = sign * dij[d]
    copy_to_row(c=i, sign= 1)
    copy_to_row(c=j, sign=-1)
    result.append(row)
  return result

def construct_integer_rigidity_matrix(n_dim, n_vertices, edge_list):
  vertices = create_fake_integer_vertices(
    n_vertices=n_vertices, n_dim=n_dim)
  return construct_numeric_rigidity_matrix(
    n_dim=n_dim, vertices=vertices, edge_list=edge_list)

def construct_float_rigidity_matrix(n_dim, n_vertices, edge_list):
  vertices = create_fake_float_vertices(
    n_vertices=n_vertices, n_dim=n_dim)
  return construct_numeric_rigidity_matrix(
    n_dim=n_dim, vertices=vertices, edge_list=edge_list)

def determine_degrees_of_freedom_integer(
      n_dim,
      n_vertices,
      edge_list,
      also_return_repeats=False):
  if (n_vertices == 0):
    if (also_return_repeats): return 0, 0
    return 0
  m = construct_integer_rigidity_matrix(
    n_dim=n_dim, n_vertices=n_vertices, edge_list=edge_list)
  if (m is None):
    if (also_return_repeats): return n_dim * n_vertices, 0
    return n_dim * n_vertices
  free_vars = integer_row_echelon_form(m=m)
  if (also_return_repeats): return len(free_vars), 0
  return len(free_vars)

class row_reduced_float_rigidity_matrix(object):

  def __init__(self, n_dim, n_vertices, edge_list):
    self.n_dim = n_dim
    self.n_vertices = n_vertices
    self.edge_list = edge_list
    self.m = None
    self.vertices = None
    self.free_vars = None
    self.repeats = -1
    if (n_vertices == 0): return
    self.construct_m()

  def construct_m(self):
    while True:
      self.repeats += 1
      self.vertices = create_fake_float_vertices(
        n_vertices=self.n_vertices, n_dim=self.n_dim)
      self.m = construct_numeric_rigidity_matrix(
        n_dim=self.n_dim, vertices=self.vertices, edge_list=self.edge_list)
      if (self.m is None): return
      self.free_vars = float_row_echelon_form(m=self.m)
      if (self.free_vars is not None):
        break

  def dof(self):
    if (self.free_vars is None):
      return self.n_dim * self.n_vertices
    return len(self.free_vars)

  def is_redundant(self, edge):
    assert self.m is not None
    while True:
      addl_row = construct_numeric_rigidity_matrix(
        n_dim=self.n_dim, vertices=self.vertices, edge_list=[edge])[0]
      result = float_row_echelon_form_is_redundant(
        m=self.m, free_vars=self.free_vars, addl_row=addl_row)
      if (result is not None):
        break
      self.construct_m()
    return result

def determine_degrees_of_freedom_float(
      n_dim,
      n_vertices,
      edge_list,
      also_return_repeats=False):
  if (also_return_repeats): return_zero = 0, 0
  else:                     return_zero = 0
  if (n_vertices == 0): return return_zero
  r = row_reduced_float_rigidity_matrix(
    n_dim=n_dim, n_vertices=n_vertices, edge_list=edge_list)
  if (also_return_repeats): return r.dof(), r.repeats
  return r.dof()

def determine_degrees_of_freedom(
      n_dim,
      n_vertices,
      edge_list,
      method="integer",
      also_return_repeats=False):
  assert method in ["integer", "float"]
  if (method == "integer"):
    return determine_degrees_of_freedom_integer(
      n_dim=n_dim, n_vertices=n_vertices, edge_list=edge_list,
      also_return_repeats=also_return_repeats)
  return determine_degrees_of_freedom_float(
    n_dim=n_dim, n_vertices=n_vertices, edge_list=edge_list,
    also_return_repeats=also_return_repeats)

# Donald Jacobs
# Generic rigidity in three-dimensional bond-bending networks
# Math. Gen. 31 (1998) 6653-6668)
# Figure 1(d)
double_banana_edge_list = [
  (0,1), (0,2), (1,2),
  (0,3), (1,3), (2,3),
  (0,4), (1,4), (2,4),
  (5,6), (5,7), (6,7),
  (5,3), (6,3), (7,3),
  (5,4), (6,4), (7,4)]

def example():
  n_vertices = 8
  edge_list = double_banana_edge_list
  for method in ["integer", "float"]:
    print("double banana 3D dof (method=%s):" % method, \
      determine_degrees_of_freedom(
        n_dim=3, n_vertices=n_vertices, edge_list=edge_list, method=method))

if (__name__ == "__main__"):
  example()
