from cctbx.eltbx import xray_scattering
from cctbx.eltbx import gaussian_fit
from cctbx.array_family import flex
import sys

def run(args):
  if ("--help" in args or len(args) not in [0,1]):
    print "usage: python wk1995_it1992.py [d_min]"
    return
  d_min = 1/4.
  if (len(args) == 1):
    d_min = float(args[0])
    assert d_min > 0
  sampling_points = \
    gaussian_fit.international_tables_sampling_points_and_value_sigmas_up_to(
      d_min=d_min)
  labels = flex.std_string()
  max_errors = flex.double()
  for wk in xray_scattering.wk1995_iterator():
    it = xray_scattering.it1992(wk.label(), 1)
    diff_gaussian = xray_scattering.difference_gaussian(wk.fetch(), it.fetch())
    labels.append(wk.label())
    max_errors.append(flex.max(gaussian_fit.get_significant_relative_errors(
      diff_gaussian, sampling_points.d_star_sq, sampling_points.sigmas)))
    gaussian_fit.write_plots(
      plots_dir="wk1995_it1992_plots",
      label=wk.label(),
      gaussians=[wk.fetch(), it.fetch()],
      d_star_sq=sampling_points.d_star_sq)
  perm = flex.sort_permutation(max_errors, 0001)
  labels = labels.select(perm)
  max_errors = max_errors.select(perm)
  for l,e in zip(labels, max_errors):
    print "%.6s %.4f" % (l,e)

if (__name__ == "__main__"):
  run(sys.argv[1:])
