from cctbx.eltbx.development import itvc_section61_io
from cctbx.eltbx import xray_scattering
import scitbx.math.gaussian_fit
import cctbx.eltbx.gaussian_fit
from cctbx.array_family import flex
from libtbx.optparse_wrapper import OptionParser
import sys, os

def run(file_name, plots_dir="itvc_wk1995_plots", quiet=0, verbose=0):
  tab = itvc_section61_io.read_table6111(file_name)
  for wk in xray_scattering.wk1995_iterator():
    label = wk.label()
    if (label in ["H'", "D"]): continue
    if (label == "Siv"):
      label = "Sival"
    for sign in ["+", "-"]:
      i = label.find(sign)
      if (i > 0):
        label = label[:i-1] + sign + label[i-1] + label[i+1:]
        break
    if (not label in tab.entries):
      print "Warning: missing scatterer:", label
  labels = flex.std_string()
  errors = []
  max_errors = flex.double()
  cmp_plots = flex.std_string()
  for element in tab.elements:
    entry = tab.entries[element]
    element_of_ion = None
    label = element
    for sign in ["+", "-"]:
      i = label.find(sign)
      if (i > 0):
        element_of_ion = label[:i]
        label = label.replace(sign,"") + sign
        break
    wk = xray_scattering.wk1995(label, 1)
    if (element_of_ion is None or not entry.table_y[-6:].all_eq(0)):
      gaussian_fit = scitbx.math.gaussian.fit(
        cctbx.eltbx.gaussian_fit.international_tables_stols,
        entry.table_y,
        cctbx.eltbx.gaussian_fit.international_tables_sigmas,
        wk.fetch())
    else:
      element_of_ion_entry = tab.entries[element_of_ion]
      patched_table_y = entry.table_y[:-6]
      patched_table_y.append(element_of_ion_entry.table_y[-6:])
      gaussian_fit = scitbx.math.gaussian.fit(
        cctbx.eltbx.gaussian_fit.international_tables_stols,
        patched_table_y,
        cctbx.eltbx.gaussian_fit.international_tables_sigmas,
        wk.fetch())
    labels.append(element)
    errors.append(gaussian_fit.significant_relative_errors())
    max_errors.append(flex.max(errors[-1]))
    if (plots_dir is not None):
      if (not os.path.isdir(plots_dir)):
        print "No plots because the directory %s does not exist." % plots_dir
        plots_dir = None
      else:
        cmp_plots.append(cctbx.eltbx.gaussian_fit.write_plots(
          plots_dir=plots_dir,
          label=element,
          gaussian_fit=gaussian_fit))
  perm = flex.sort_permutation(max_errors, 0001)
  labels = labels.select(perm)
  errors = flex.select(errors, perm)
  if (plots_dir is None):
    cmp_plots = [None] * len(labels)
  else:
    cmp_plots = cmp_plots.select(perm)
  for l,e,p in zip(labels, errors, cmp_plots):
    entry = tab.entries[l]
    y = entry.table_y
    perm = flex.sort_permutation(e, 0001)[:3]
    high = []
    for i in perm:
      s = cctbx.eltbx.gaussian_fit.international_tables_stols[i]
      a = ""
      if (not quiet and s < 2.1): a = "@%.3f" % y[i]
      high.append("%7.4f %4.2f%s" % (e[i],s,a))
    print "Element %-2s(%2d):" % (l, entry.atomic_number), ", ".join(high)
    if (verbose and p is not None):
      print p
      sys.stdout.write(open(p).read())
      print

def main():
  parser = OptionParser(
    usage="usage: python %prog [options] file_name")
  parser.add_option("-q", "--quiet",
                    action="store_true", dest="quiet", default=0,
                    help="do not show values for large errors, only stol")
  parser.add_option("-v", "--verbose",
                    action="store_true", dest="verbose", default=0,
                    help="show comparison table for each element")
  (options, args) = parser.parse_args()
  if (len(args) != 1):
    parser.print_help()
    return
  run(file_name=args[0], quiet=options.quiet, verbose=options.verbose)

if (__name__ == "__main__"):
  main()
