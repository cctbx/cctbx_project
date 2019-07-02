from __future__ import absolute_import, division, print_function
from cctbx.eltbx.development import itvc_section61_io
from cctbx.eltbx import xray_scattering
import scitbx.math.gaussian_fit
import cctbx.eltbx.gaussian_fit
from cctbx.array_family import flex
from libtbx.option_parser import OptionParser
import sys, os
from six.moves import range
from six.moves import zip

def run(file_name, table_of_gaussians, cutoff,
        low_resolution_only=False,
        high_resolution_only=False,
        significant_errors_only=False,
        plots_dir=None,
        quiet=0,
        verbose=0):
  assert not (low_resolution_only and high_resolution_only)
  tab = itvc_section61_io.read_table6111(file_name)
  for wk in xray_scattering.wk1995_iterator():
    label = wk.label()
    if (not label in tab.entries):
      print("Warning: missing scatterer:", label)
  stols = cctbx.eltbx.gaussian_fit.international_tables_stols
  sel = stols <= cutoff + 1.e-6
  stols = stols.select(sel)
  if (low_resolution_only):
    sel = stols <= 2
    stols = stols.select(sel)
    assert stols.size() == 56
  elif (high_resolution_only):
    sel = stols > 2
    stols = stols.select(sel)
    assert stols.size() == 6
  range_62 = flex.size_t(range(62))
  labels = flex.std_string()
  errors = []
  correlations = flex.double()
  max_errors = flex.double()
  cmp_plots = flex.std_string()
  for element in tab.elements:
    entry = tab.entries[element]
    wk = table_of_gaussians(element, 1)
    assert entry.table_y.size() == 62
    if (not flex.sort_permutation(data=entry.table_y, reverse=True)
            .all_eq(range_62)):
      print("Increasing: %s (%d)" % (element, entry.atomic_number))
      prev_y = entry.table_y[0]
      for y in entry.table_y:
        if (y > prev_y):
          print("higher:", y, "before:", prev_y)
        prev_y = y
      raise RuntimeError("Data values are not increasing.")
    if (low_resolution_only):
      gaussian_fit = scitbx.math.gaussian.fit(
        stols,
        entry.table_y[:-6],
        entry.table_sigmas[:-6],
        wk.fetch())
    elif (high_resolution_only):
      gaussian_fit = scitbx.math.gaussian.fit(
        stols,
        entry.table_y[-6:],
        entry.table_sigmas[-6:],
        wk.fetch())
    elif (    entry.element != entry.atomic_symbol
          and entry.table_y[-6:].all_eq(0)):
      atom_entry = tab.entries[entry.atomic_symbol]
      patched_table_y = entry.table_y[:-6]
      patched_table_y.append(atom_entry.table_y[-6:])
      patched_table_sigmas = entry.table_sigmas[:-6]
      patched_table_sigmas.append(atom_entry.table_sigmas[-6:])
      gaussian_fit = scitbx.math.gaussian.fit(
        stols,
        patched_table_y,
        patched_table_sigmas,
        wk.fetch())
    else:
      gaussian_fit = scitbx.math.gaussian.fit(
        stols,
        entry.table_y[:stols.size()],
        entry.table_sigmas[:stols.size()],
        wk.fetch())
    labels.append(element)
    errors.append(gaussian_fit.significant_relative_errors())
    max_errors.append(flex.max(errors[-1]))
    correlations.append(flex.linear_correlation(
      gaussian_fit.table_y(), gaussian_fit.fitted_values()).coefficient())
    if (plots_dir is not None):
      if (not os.path.isdir(plots_dir)):
        print("No plots because the directory %s does not exist." % plots_dir)
        plots_dir = None
      else:
        cmp_plots.append(cctbx.eltbx.gaussian_fit.write_plots(
          plots_dir=plots_dir,
          label=element,
          gaussian_fit=gaussian_fit))
  perm = flex.sort_permutation(data=max_errors, reverse=True)
  labels = labels.select(perm)
  errors = flex.select(errors, perm)
  correlations = correlations.select(perm)
  if (plots_dir is None):
    cmp_plots = [None] * len(labels)
  else:
    cmp_plots = cmp_plots.select(perm)
  for l,e,cc,p in zip(labels, errors, correlations, cmp_plots):
    entry = tab.entries[l]
    y = entry.table_y
    perm = flex.sort_permutation(data=e, reverse=True)[:3]
    high = []
    for i in perm:
      if (significant_errors_only and e[i] < 0.01): break
      s = stols[i]
      a = ""
      if (not quiet and s < 2.1): a = "@%.3f" % y[i]
      high.append("%7.4f %4.2f%s" % (e[i],s,a))
      if (high_resolution_only): break
    if (verbose or len(high) > 0):
      print("Element %-5s(%2d) cc=%.4f:" % (
        l, entry.atomic_number, cc), ", ".join(high))
    if (verbose and p is not None):
      print(p)
      sys.stdout.write(open(p).read())
      print()

def main():
  parser = OptionParser(
    usage="usage: python %prog [options] file_name")
  parser.add_option("-t", "--table", default="wk1995",
    action="store", dest="table_of_gaussians",
    help="wk1995 or it1992")
  parser.add_option("-q", "--quiet",
    action="store_true", default=0,
    help="do not show values for large errors, only stol")
  parser.add_option("-v", "--verbose",
    action="store_true", default=0,
    help="show comparison table for each element")
  parser.add_option("-l", "--low_resolution_only",
    action="store_true", default=0,
    help="analyze points up to sin(theta)/lambda=2A-1 only")
  parser.add_option("-g", "--high_resolution_only",
    action="store_true", default=0,
    help="analyze points beyond sin(theta)/lambda=2A-1 only")
  parser.add_option("-s", "--significant_errors_only",
    action="store_true", default=0,
    help="show errors greater than 1% only")
  (options, args) = parser.parse_args()
  if (len(args) != 1):
    parser.print_help()
    return
  assert options.table_of_gaussians in ["wk1995", "it1992"]
  if (options.table_of_gaussians == "wk1995"):
    table_of_gaussians = xray_scattering.wk1995
    plots_dir = "itvc_wk1995_plots"
    cutoff = 6
  else:
    table_of_gaussians = xray_scattering.it1992
    plots_dir = "itvc_it1992_plots"
    cutoff = 2
    assert not options.low_resolution_only
    assert not options.high_resolution_only
  print("table_of_gaussians:", options.table_of_gaussians)
  run(
    file_name=args[0],
    table_of_gaussians=table_of_gaussians,
    cutoff=cutoff,
    low_resolution_only=options.low_resolution_only,
    high_resolution_only=options.high_resolution_only,
    significant_errors_only=options.significant_errors_only,
    plots_dir=plots_dir,
    quiet=options.quiet,
    verbose=options.verbose)

if (__name__ == "__main__"):
  main()
