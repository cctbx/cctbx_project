from cctbx.eltbx import xray_scattering
from cctbx.array_family import flex
import sys, os

def it1974_vol_4_table_2_2a_points_stol():
  return flex.double(
    [0.00, 0.01, 0.02, 0.03, 0.04, 0.05,
     0.06, 0.07, 0.08, 0.09, 0.10,
     0.11, 0.12, 0.13, 0.14, 0.15,
     0.16, 0.17, 0.18, 0.19, 0.20,
     0.22, 0.24, 0.25, 0.26, 0.28, 0.30,
     0.32, 0.34, 0.35, 0.36, 0.38, 0.40,
     0.42, 0.44, 0.45, 0.46, 0.48, 0.50,
     0.55, 0.60, 0.65, 0.70, 0.80, 0.90, 1.00,
     1.10, 1.20, 1.30, 1.40, 1.50,
     1.60, 1.70, 1.80, 1.90, 2.00])

def get_max_error(diff_gaussian, d_star_sq):
  reference_values = diff_gaussian.reference_gaussian().at_d_star_sq(d_star_sq)
  errors = diff_gaussian.target_terms_at_points(d_star_sq)
  assert flex.min(reference_values) > 0
  return flex.max(flex.abs(errors/reference_values))

def write_plot(f, gaussian, d_star_sq):
  for d_star in flex.sqrt(d_star_sq):
    print >> f, d_star, gaussian.at_d_star_sq(d_star**2)
  print >> f, "&"

def write_plots(plots_dir, label, diff_gaussian, d_star_sq):
  label = label.replace("'", "prime")
  file_name = os.path.join(plots_dir, label+".xy")
  f = open(file_name, "w")
  write_plot(f, diff_gaussian.reference_gaussian(), d_star_sq)
  write_plot(f, diff_gaussian, d_star_sq)
  f.close()

def run(args):
  if ("--help" in args or len(args) not in [0,1]):
    print "usage: python wk1995_it1992.py [d_min]"
    return
  d_star_sq = flex.pow2(it1974_vol_4_table_2_2a_points_stol() * 2)
  if (len(args) == 1):
    d_star_sq_max = 1/float(args[0])**2
    d_star_sq = d_star_sq.select(d_star_sq <= d_star_sq_max)
  labels = flex.std_string()
  max_errors = flex.double()
  for wk in xray_scattering.wk1995_iterator():
    it = xray_scattering.it1992(wk.label(), 1)
    diff_gaussian = xray_scattering.difference_gaussian(wk.fetch(), it.fetch())
    labels.append(wk.label())
    max_errors.append(get_max_error(diff_gaussian, d_star_sq))
    write_plots("wk1995_it1992_plots", wk.label(), diff_gaussian, d_star_sq)
  perm = flex.sort_permutation(max_errors, 0001)
  labels = labels.select(perm)
  max_errors = max_errors.select(perm)
  for l,e in zip(labels, max_errors):
    print "%.6s %.4f" % (l,e)

if (__name__ == "__main__"):
  run(sys.argv[1:])
