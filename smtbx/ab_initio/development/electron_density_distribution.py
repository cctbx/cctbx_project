from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex

def find_delta(rho_map, tol):
  """ Find delta as hinted on fig. 1 of ref. [1] in module charge_flipping """
  rho = rho_map.real_map_unpadded().as_1d()
  max_rho = flex.max(rho)
  rho /= max_rho
  sorting = flex.sort_permutation(rho)
  sorted_rho = rho.select(sorting)
  n = len(sorted_rho)
  p,q = n//4, 3*n//4
  indexes = flex.double_range(p,q)
  values = sorted_rho[p:q]
  c = flex.linear_correlation(indexes, values)
  assert c.is_well_defined() and c.coefficient() > 0.99
  r = flex.linear_regression(indexes, values)
  a,b = r.y_intercept(), r.slope()
  deviation = flex.abs(a + b*flex.double_range(n) - sorted_rho)
  non_linear_sel = deviation > tol
  low = flex.first_index(non_linear_sel, False)
  high = flex.last_index(non_linear_sel, False)
  assert non_linear_sel[low:high].count(False)/(high-low+1) > 0.99
  assert sorted_rho[low] < 0 and sorted_rho[high] > 0
  return min(sorted_rho[high], -sorted_rho[low]), max_rho

def write_sorted_moduli_as_mathematica_plot(f, filename):
  """ To obtain fig. 1 in ref [2] in module charge_flipping """
  abs_f = flex.abs(f.data())
  sorted = abs_f.select(flex.sort_permutation(abs_f))
  sorted /= flex.max(sorted)
  mf = open(os.path.expanduser(filename), 'w')
  print('fp1 = {', file=mf)
  for f in sorted:
    print("%f, " % f, file=mf)
  print("1 };", file=mf)
  print("ListPlot[fp1]", file=mf)
  mf.close()
