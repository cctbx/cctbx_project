from __future__ import division
import math
import scitbx.math
try: import tntbx
except ImportError: tntbx = None
import libtbx.utils
import sys
import random

from scitbx.array_family import flex


class test_case(object):

  n_matrices_per_dimension = 1
  full_coverage = True
  show_progress = False
  exercise_tntbx = False
  eps=scitbx.math.double_numeric_limits.epsilon

  def __init__(self, **kwds):
    libtbx.adopt_optional_init_args(self, kwds)
    self.exercise_tntbx = self.exercise_tntbx and tntbx is not None
    if self.exercise_tntbx:
      self.accumulate_u = self.accumulate_v = True
    else:
      self.accumulate_u = self.accumulate_v = False

  def matrices(self):
    for k in (1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
              2, 3, 4, 5, 6, 7, 8, 9, 10,
              20, 30, 40, 50, 60, 70, 80, 90, 100):
      for n in (10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
        if not self.full_coverage and random.random() < 0.9: continue
        m = int(k*n)
        gen = scitbx.math.random_normal_matrix_generator(m, n)
        for p in xrange(self.n_matrices_per_dimension):
          sigma = self.sigma(m,n)
          yield sigma, gen.matrix_with_singular_values(sigma)

  def exercise_increasing_dimensions(self):
    for sigma, a in self.matrices():
      m, n = a.focus()
      if self.show_progress:
        print (m,n),
        sys.stdout.flush()
      svd = scitbx.math.svd.real(a, self.accumulate_u, self.accumulate_v)
      if self.show_progress:
        print '!',
        sys.stdout.flush()
      sigma = sigma.select(flex.sort_permutation(sigma, reverse=True))
      delta = (svd.sigma - sigma).norm()/sigma.norm()/min(m,n)/self.eps
      assert delta < 10

      if not self.exercise_tntbx:
        if self.show_progress: print
        continue
      svd = tntbx.svd_m_ge_n_double(a)
      if self.show_progress:
        print '!'
        sys.stdout.flush()
      sigma = sigma.select(flex.sort_permutation(sigma, reverse=True))
      delta = ((svd.singular_values() - sigma).norm()
               /sigma.norm()/min(m,n)/self.eps)
      assert delta < 10

  def time(self):
    from libtbx.easy_profile import easy_profile
    self.scitbx_report = []
    self.tntbx_report = []
    prof_scitbx = easy_profile(scitbx.math.svd.real,
                                file_name='svd.py', func_name='__init__',
                                line=None)
    if tntbx is not None:
      prof_tntbx = easy_profile(lambda m: tntbx.svd_m_ge_n_double(m),
                                  file_name='tst_svd.py', func_name='<lambda>',
                                  line=None)
    for sigma, a in self.matrices():
      m, n = a.focus()
      if self.show_progress:
        print (m,n),
        sys.stdout.flush()
      flops = min(2*m*n**2 - 2*n**3/3, m*n**2 + n**3)
      runs = max(int(1e4/flops), 1)
      if tntbx:
        prof_tntbx.runs = runs
        sys.stdout.flush()
        t = prof_tntbx.time(a)
        if self.show_progress:
          print '!',
        self.tntbx_report.append((m, n, t))
      prof_scitbx.runs = runs
      t = prof_scitbx.time(a, accumulate_u=True, accumulate_v=True)
      if self.show_progress:
        print '!'
      self.scitbx_report.append((m, n, t))
    print 'scitbx = {'
    print ', '.join([ "{%i, %i, %f}" % (m, n, t)
                      for m, n, t in self.scitbx_report ])
    print '};'
    if tntbx:
      print '(***************************************************)'
      print '(***************************************************)'
      print 'tntbx = {'
      print ', '.join([ "{%i, %i, %f}" % (m, n, t)
                        for m, n, t in self.tntbx_report ])
      print '};'


class chunks_of_small_and_big_singular_values_case(test_case):
  """ Top left half of singular values small, lower right big """

  big = 1e50

  def sigma(self, m, n):
    p = min(m,n)
    sigma = flex.random_double(p - p//2, factor=1/self.big)
    sigma.extend(flex.random_double(p//2, factor=self.big) + self.big)
    return sigma


def run(show_progress, exercise_tntbx, full_coverage):
  t = chunks_of_small_and_big_singular_values_case(
    show_progress=show_progress,
    exercise_tntbx=exercise_tntbx,
    full_coverage=full_coverage,
  )
  t.exercise_increasing_dimensions()
  print libtbx.utils.format_cpu_times()

def time(show_progress, exercise_tntbx):
  t = chunks_of_small_and_big_singular_values_case(
    show_progress=show_progress,
    exercise_tntbx=exercise_tntbx)
  t.time()

if __name__ == '__main__':
  show_progress = '--show-progress' in sys.argv[1:]
  exercise_tntbx = '--exercise-tntbx' in sys.argv[1:]
  full_coverage = '--full-coverage' in sys.argv[1:]
  if '--time' not in sys.argv[1:]:
    run(show_progress=show_progress, exercise_tntbx=exercise_tntbx,
        full_coverage=full_coverage)
  else:
    time(show_progress=show_progress, exercise_tntbx=exercise_tntbx)
