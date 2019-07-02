from __future__ import absolute_import, division, print_function
import itertools
from smtbx.development import random_xray_structure
import smtbx.utils
from smtbx.refinement import constraints, least_squares
from scitbx.array_family import flex
import timeit
import pprint
from six.moves import range

def run(sizes, d_mins):
  timing = {
    'Serial': {
      'BLAS 2': [],
      'BLAS 3': []
    },
    'Parallel': {
      'BLAS 2': [],
      'BLAS 3': []
    },
  }
  print("%16s%11s%7s%9s%8s" % ('Algorithms', 'resolution', '#atoms',
                               'building', 'solving'))
  for observations, model in test_structures(sizes, d_mins):
    for blas in ('BLAS 2', 'BLAS 3'):
      for may_parallelise in (False, True):
        concurrency = ('Serial', 'Parallel')[may_parallelise]
        result = benchmark(observations, model, may_parallelise, blas)
        print("%9s%7s%11.2f%7i%9.3f%8.3f" % ((concurrency, blas) + result))
        timing[concurrency][blas].append(result)
  with open('smtbx_refinement.m', 'w') as f:
    pprint.pprint(timing, stream=f)

non_linear_ls_engine = {
  'BLAS 2':
  least_squares.normal_eqns.non_linear_ls_with_separable_scale_factor_BLAS_2,
  'BLAS 3':
  least_squares.normal_eqns.non_linear_ls_with_separable_scale_factor_BLAS_3,
}
def benchmark(observations, model, may_parallelise, blas):
  time = timeit.default_timer
  ls = least_squares.crystallographic_ls(
    observations, model, non_linear_ls_engine[blas], may_parallelise,
    weighting_scheme=least_squares.mainstream_shelx_weighting(a=0))
  m = ls.observations.fo_sq.size()
  n = ls.reparametrisation.n_independents
  # let's do worth of 5 Gflops at least
  n_trials = max(int(5e9/(0.5*m*n**2)), 1)
  building = 0
  solving = 0
  for i in range(n_trials):
    t0 = time()
    ls.build_up()
    t1 = time()
    rls = ls.reduced_problem()
    neqns = rls.step_equations()
    neqns.solve()
    t2 = time()
    building += t1 - t0
    solving += t2 - t1
  return (ls.observations.fo_sq.d_min(), len(ls.xray_structure.scatterers()),
         building/n_trials, solving/n_trials)

def test_structures(sizes, d_mins):
  for na in sizes:
    xs = random_xray_structure(
      space_group_symbol='hall: -P 2ybc',
      n_scatterers=na,
      proportion_of_elements={'C':5, 'O':2, 'N':1},
      use_u_iso=False,
      use_u_aniso=True,
    )
    for d_min in d_mins:
      mi = xs.build_miller_set(anomalous_flag=False, d_min=d_min)
      ma = mi.structure_factors_from_scatterers(xs, algorithm='direct').f_calc()
      fo_sq = ma.norm().customized_copy(sigmas=flex.double(ma.size(), 1.))
      xs.shake_sites_in_place(rms_difference=0.1)
      xs.shake_adp()
      for sc in xs.scatterers():
        sc.flags.set_use_u_iso(False).set_use_u_aniso(True)
        sc.flags.set_grad_site(True).set_grad_u_aniso(True)
      connectivity_table = smtbx.utils.connectivity_table(xs)
      reparametrisation = constraints.reparametrisation(
        structure=xs,
        constraints=[],
        connectivity_table=connectivity_table)
      yield (fo_sq.as_xray_observations(), reparametrisation)

if __name__ == '__main__':
  import argparse
  def arange(s, atype):
    parsed = tuple(atype(x) for x in s.split(':'))
    start, stop = parsed[:2] if len(parsed) >= 2 else parsed*2
    step = parsed[2] if len(parsed) == 3 else 1
    if atype == int:
      for i in range(start, stop, step):
        yield i
      yield stop
    elif atype == float:
      i = start
      yield i
      i += step
      while i < stop:
        yield i
        i += step
      yield stop
  def sizerange(s):
    for i in arange(s, int):
      yield i
  def floatrange(s):
    for i in  arange(s, float):
      yield i
  p = argparse.ArgumentParser(description='Benchmark smtbx refinement')
  p.add_argument('--atoms', metavar='FIRST:LAST:STEP',
                 type=sizerange, action='append',
                 help='A range of number of atoms in the ASU (Python style)')
  p.add_argument('--dmin', metavar='FIRST:LAST:STEP',
                 type=floatrange, action='append',
                 help='A range of d_min')
  cmd = p.parse_args()
  sizes, d_mins = (
    sorted(set(itertools.chain.from_iterable(opt)))
    for opt in (cmd.atoms, cmd.dmin))
  run(sizes, d_mins)
