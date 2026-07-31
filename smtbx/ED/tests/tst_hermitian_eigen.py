from __future__ import absolute_import, division, print_function

"""The small Hermitian eigensolver the dynamical ED code runs on every beam group.

It replaced a LAPACK driver, which is a thing worth having a test for on its
own: nothing else in the suite would notice it going wrong, and what it would
produce instead is not a crash but a slightly wrong refinement.

No reference decomposition is needed to check it. If V is unitary and
A V = V diag(w) with w real, that *is* the eigendecomposition of A, by
definition -- so those two properties, plus the ascending order the callers
depend on, are the whole test. The matrices are chosen to reach the parts of
the algorithm a random one would not: repeated eigenvalues, a matrix already
split into blocks, one already diagonal, and one whose elements span eight
orders of magnitude.
"""

from smtbx_ed_data_ext import hermitian_eigen
from scitbx.array_family import flex

import random


def as_flex(rows):
  n = len(rows)
  a = flex.complex_double(flex.grid(n, n))
  for i in range(n):
    for j in range(n):
      a[i, j] = rows[i][j]
  return a


def hermitian(n, rng, scale_spread=0):
  rows = [[0j]*n for _ in range(n)]
  for i in range(n):
    rows[i][i] = complex(rng.uniform(-1, 1), 0)
    for j in range(i + 1, n):
      v = complex(rng.uniform(-1, 1), rng.uniform(-1, 1))
      if scale_spread:
        v *= 10**rng.uniform(-scale_spread, scale_spread)
      rows[i][j] = v
      rows[j][i] = v.conjugate()
  return rows


def diagonal(n, rng):
  rows = [[0j]*n for _ in range(n)]
  for i in range(n):
    rows[i][i] = complex(rng.uniform(-1, 1), 0)
  return rows


def block_split(n, rng):
  """Two blocks with nothing between them: exactly zero subdiagonals."""
  rows = hermitian(n, rng)
  for i in range(n):
    for j in range(n):
      if (i < n//2) != (j < n//2):
        rows[i][j] = 0j
  return rows


def degenerate(n, rng):
  """U diag(1,1,1,2,2,2,...) U^H -- eigenvalues repeated in threes."""
  u = [[complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(n)]
       for _ in range(n)]
  for k in range(n):                      # Gram-Schmidt, by column
    for j in range(k):
      d = sum(u[i][j].conjugate()*u[i][k] for i in range(n))
      for i in range(n):
        u[i][k] -= d*u[i][j]
    nrm = sum(abs(u[i][k])**2 for i in range(n))**0.5
    for i in range(n):
      u[i][k] /= nrm
  lam = [float(k//3 + 1) for k in range(n)]
  return [[sum(u[i][k]*lam[k]*u[j][k].conjugate() for k in range(n))
           for j in range(n)] for i in range(n)]


def check(rows, tol=1e-11):
  n = len(rows)
  a = as_flex(rows)
  w, v = hermitian_eigen(a)
  assert len(w) == n

  scale = max([abs(x) for x in a] + [1e-300])

  # A V = V diag(w)
  for k in range(n):
    for i in range(n):
      lhs = sum(rows[i][j]*v[j, k] for j in range(n))
      assert abs(lhs - w[k]*v[i, k]) <= tol*scale, (n, i, k)

  # V^H V = I
  for p in range(n):
    for q in range(n):
      s = sum(v[i, p].conjugate()*v[i, q] for i in range(n))
      assert abs(s - (1 if p == q else 0)) <= tol, (n, p, q)

  for k in range(n - 1):
    assert w[k] <= w[k + 1], (n, k)

  # the input must come back untouched
  for i in range(n):
    for j in range(n):
      assert a[i, j] == rows[i][j]


def run():
  rng = random.Random(0)
  for n in range(1, 13):
    for _ in range(4):
      check(hermitian(n, rng))
      check(diagonal(n, rng))
      check(block_split(n, rng))
      check(degenerate(n, rng))
      check(hermitian(n, rng, scale_spread=8))
  print('OK')


if __name__ == '__main__':
  run()
