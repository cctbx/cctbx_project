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
  stol_max = 1/(2*d_min)
  n_points = gaussian_fit.n_less_than(
    sorted_array=gaussian_fit.international_tables_sampling_stols,
    cutoff=stol_max)
  labels = flex.std_string()
  max_errors = flex.double()
  for wk in xray_scattering.wk1995_iterator():
    it = xray_scattering.it1992(wk.label(), 1)
    fit_object = xray_scattering.gaussian_fit(
      gaussian_fit.international_tables_sampling_stols[:n_points],
      wk.fetch(),
      gaussian_fit.international_tables_sampled_value_sigmas[:n_points],
      it.fetch())
    labels.append(wk.label())
    max_errors.append(flex.max(gaussian_fit.get_significant_relative_errors(
      fit_object=fit_object)))
    gaussian_fit.write_plots(
      plots_dir="wk1995_it1992_plots",
      label=wk.label(),
      fit_object=fit_object)
  perm = flex.sort_permutation(max_errors, 0001)
  labels = labels.select(perm)
  max_errors = max_errors.select(perm)
  for l,e in zip(labels, max_errors):
    print "%.6s %.4f" % (l,e)

if (__name__ == "__main__"):
  run(sys.argv[1:])
