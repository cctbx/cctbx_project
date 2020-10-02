from __future__ import absolute_import, division, print_function

class dict_with_default_0(dict):
  def __missing__(self, key):
    return 0

def iiexps_mul(a, b):
  result = dict_with_default_0(a)
  for ii,exp in b.items():
    result[ii] += exp
  return result

def iiexps_as_tuples(iiexps):
  return tuple(sorted(iiexps.items()))

class term(object):

  def __init__(self, f, iiexps):
    self.f, self.iiexps = f, iiexps

  def __neg__(self):
    return term(-self.f, self.iiexps)

  def __mul__(a, b):
    return term(a.f * b.f, iiexps_mul(a.iiexps, b.iiexps))

  def __str__(self):
    result = []
    if (self.f != 1 and self.f != -1):
      result.append("%d" % self.f)
    for ii,exp in iiexps_as_tuples(self.iiexps):
      s = "v%d%s" % (ii//3, "xyz"[ii%3])
      if (exp != 1): s += "**%d" % exp
      result.append(s)
    result = "*".join(result)
    if (self.f == -1):
      result = "-" + result
    return result

elem_common_factor_check = True

class elem(object):

  def __init__(self, terms):
    n = dict_with_default_0()
    for t in terms:
      n[iiexps_as_tuples(t.iiexps)] += t.f
    terms = []
    for ies,f in n.items():
      if (f != 0):
        terms.append(term(f=f, iiexps=dict_with_default_0(ies)))
    if (elem_common_factor_check):
      # assert that there are no common factors
      count_i = {}
      min_e = {}
      for it in range(len(terms)):
        t = terms[it]
        for i,e in t.iiexps.items():
          ci = count_i.get(i, 0)
          if (ci == it):
            if (ci == 0):
              count_i[i] = 1
              min_e[i] = e
            else:
              count_i[i] += 1
              min_e[i] = min(min_e[i], e)
          elif (ci != 0):
            del count_i[i]
            del min_e[i]
      for i,c in count_i.items():
        assert c != len(terms)
    self.terms = terms

  def __eq__(a, b):
    assert b == 0
    return len(a.terms) == 0

  def __ne__(a, b):
    assert b == 0
    return len(a.terms) != 0

  def __mul__(a, b):
    terms = []
    for ta in a.terms:
      for tb in b.terms:
        terms.append(ta * tb)
    return elem(terms=terms)

  def __sub__(a, b):
    return elem(terms=a.terms + [-tb for tb in b.terms])

  def __str__(self):
    result = "+".join([str(t) for t in self.terms]).replace("+-","-")
    if (len(result) != 0): return result
    return "0"

def symbolic_row_echelon_form(m):
  free_vars = []
  n_rows = len(m)
  n_cols = len(m[0])
  piv_r = 0
  piv_c = 0
  while (piv_c < n_cols):
    best_r = None
    for i_row in range(piv_r, n_rows):
      if (m[i_row][piv_c] != 0):
        best_r = i_row
        break
    if (best_r is not None):
      if (best_r != piv_r):
        m[piv_r], m[best_r] = m[best_r], m[piv_r]
      fp = m[piv_r][piv_c]
      for r in range(piv_r+1, n_rows):
        fr = m[r][piv_c]
        if (fr == 0): continue
        for c in range(piv_c, n_cols):
          m[r][c] = m[r][c] * fp - m[piv_r][c] * fr
      piv_r += 1
    else:
      free_vars.append(piv_c)
    piv_c += 1
  return free_vars

def construct_symbolic_rigidity_matrix(n_dim, n_vertices, edge_list):
  n_edges = len(edge_list)
  if (n_edges == 0): return None
  n_columns = n_vertices * n_dim
  result = []
  for i,j in edge_list:
    assert i != j
    row = []
    for ic in range(n_columns):
      row.append(elem([]))
    ci = i * n_dim
    cj = j * n_dim
    for d in range(n_dim):
      ei = {ci+d: 1}
      ej = {cj+d: 1}
      row[ci+d] = elem(terms=[term(1, ei), term(-1, ej)])
      row[cj+d] = elem(terms=[term(1, ej), term(-1, ei)])
    result.append(row)
  return result

def determine_degrees_of_freedom(n_dim, n_vertices, edge_list):
  if (n_vertices == 0): return None
  m = construct_symbolic_rigidity_matrix(
    n_dim=n_dim, n_vertices=n_vertices, edge_list=edge_list)
  if (m is None): return None
  free_vars = symbolic_row_echelon_form(m=m)
  return len(free_vars)

def compare_int_symb(n_vertices, edge_list):
  from scitbx.graph.rigidity import determine_degrees_of_freedom as dof_int
  di = dof_int(
    n_dim=3, n_vertices=n_vertices, edge_list=edge_list)
  ds = determine_degrees_of_freedom(
    n_dim=3, n_vertices=n_vertices, edge_list=edge_list)
  if (0):
    print("dof_int:", di, "dof_symb:", ds)
  assert ds == di
  if (len(edge_list) > 1):
    for i_delete in range(len(edge_list)):
      el = list(edge_list)
      del el[i_delete]
      compare_int_symb(n_vertices, el)

def exercise():
  if (elem_common_factor_check):
    try: elem([term(1,{0:2}), term(1,{0:1})])
    except AssertionError: pass
    else: raise RuntimeError("Exception expected.")
  n_vertices = 2
  edge_list = [(0,1)]
  compare_int_symb(n_vertices, edge_list)
  n_vertices = 3
  edge_list.extend([(0,2),(1,2)])
  compare_int_symb(n_vertices=n_vertices, edge_list=edge_list)
  n_vertices = 4
  edge_list.extend([(0,3), (1,3), (2,3)])
  compare_int_symb(n_vertices=n_vertices, edge_list=edge_list)
  if (0):
    n_vertices = 5
    edge_list.extend([(0,4), (1,4), (2,4), (3,4)])
    compare_int_symb(n_vertices=n_vertices, edge_list=edge_list)

if (__name__ == "__main__"):
  import os
  try:
    exercise()
  finally:
    t = os.times()
    print("u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1]))
