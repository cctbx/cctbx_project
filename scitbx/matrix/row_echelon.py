def form_rational(m, t=None):
  free_vars = []
  n_rows = len(m)
  assert n_rows != 0
  n_cols = len(m[0])
  assert t is None or len(t) == n_rows
  piv_r = 0
  for piv_c in xrange(0, n_cols):
    for i_row in xrange(piv_r, n_rows):
      if (m[i_row][piv_c] != 0):
        break
    else:
      free_vars.append(piv_c)
      continue
    if (i_row != piv_r):
      m[piv_r], m[i_row] = m[i_row], m[piv_r]
      if (t != None):
        t[piv_r], t[i_row] = t[i_row], t[piv_r]
    fp = m[piv_r][piv_c]
    for r in xrange(piv_r+1, n_rows):
      fr = m[r][piv_c]
      if (fr == 0): continue
      frp = fr / fp
      for c in xrange(piv_c, n_cols):
        m[r][c] -= m[piv_r][c] * frp
      if (t is not None):
        t[r] -= t[piv_r] * frp
    piv_r += 1
  return free_vars

def back_substitution_rational(m, t, free_vars, sol):
  n_rows = len(m)
  n_cols = len(m[0])
  assert t is None or len(t) == n_rows
  assert len(sol) == n_cols
  if (t is not None):
    rank = n_cols - len(free_vars)
    for r in xrange(rank, n_rows):
      if (t[r] != 0):
        return None
  free_flags = [False] * n_cols
  for c in free_vars:
    free_flags[c] = True
  piv_cols = []
  for c,f in enumerate(free_flags):
    if (not f): piv_cols.append(c)
  for r in xrange(len(piv_cols)-1,-1,-1):
    piv_c = piv_cols[r]
    if (t is None): s = 0
    else:           s = -t[r]
    for c in xrange(piv_c+1, n_cols):
      s += m[r][c] * sol[c]
    sol[piv_c] = -s / m[r][piv_c]
  return sol
