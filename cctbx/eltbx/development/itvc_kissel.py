from __future__ import absolute_import, division, print_function
from cctbx.eltbx.development import itvc_section61_io
from cctbx.eltbx.development import kissel_io
import cctbx.eltbx.gaussian_fit
from cctbx.array_family import flex
from libtbx.option_parser import OptionParser
import cStringIO as StringIO
import sys, os
from six.moves import zip

def run(args, cutoff, high_resolution_only,
        plots_dir="itvc_kissel_plots", verbose=0):
  itab = itvc_section61_io.read_table6111(args[0])
  itab_x = cctbx.eltbx.gaussian_fit.international_tables_stols
  isel = itab_x <= cutoff + 1.e-6
  if (high_resolution_only):
    isel &= itab_x > 2 + 1.e-6
  itab_x = itab_x.select(isel)
  for file_name in args[1:]:
    ktab = kissel_io.read_table(file_name)
    if (ktab.element == "Es"): continue
    sel = ktab.x <= cutoff + 1
    ktab_x = ktab.x.select(sel)
    ktab_y = ktab.y.select(sel)
    ktabs_sigmas = ktab.sigmas.select(sel)
    itab_entry = itab.entries[ktab.element]
    itab_y = itab_entry.table_y.select(isel)
    itab_sigmas = itab_entry.table_sigmas.select(isel)
    f = open(os.path.join(plots_dir, ktab.element + ".xy"), "w")
    cctbx.eltbx.gaussian_fit.write_plot(f, ktab_x, ktab_y)
    cctbx.eltbx.gaussian_fit.write_plot(f, itab_x, itab_y)
    f.close()
    ktab_y_i = flex.linear_interpolation(ktab.x, ktab.y, itab_x)
    ktab_sigmas_i = flex.linear_interpolation(ktab.x, ktab.sigmas, itab_x)
    assert ktab_y_i.all_gt(0)
    s = StringIO.StringIO()
    print("stol  kissel    itvc   delta sig_itvc rel_sig rel_del tol_del", file=s)
    max_delta = 0
    max_tol_del = 0
    for x,ky,ksig,iy,isig in zip(itab_x, ktab_y_i, ktab_sigmas_i,
                                         itab_y, itab_sigmas):
      if (iy > 0): ie = "%7.4f" % abs(isig/iy)
      else: ie = " ******"
      delta = iy - ky
      rel_del = abs(delta) / ky
      tol_del = max(0, abs(delta)-ksig-isig) / ky
      print("%4.2f %7.4f %7.4f %7.4f %8.5f %-s %7.4f %7.4f" % (
        x,ky,iy,delta,isig,ie,rel_del,tol_del), file=s)
      max_delta = max(max_delta, abs(delta))
      max_tol_del = max(max_tol_del, tol_del)
    print("Element:", ktab.element, "max_delta=%.4f, max_tol_del=%.4f" % (
      max_delta, max_tol_del))
    sys.stdout.write(s.getvalue())
    print()

def main():
  parser = OptionParser(
    usage="usage: python %prog [options] itvc_table kissel_files ...")
  parser.add_option("-c", "--cutoff",
    type="float", default=6, metavar="FLOAT",
    help="maximum sin(theta)/lambda")
  parser.add_option("-g", "--high_resolution_only",
    action="store_true", default=0,
    help="analyze points beyond sin(theta)/lambda=2A-1 only")
  (options, args) = parser.parse_args()
  if (len(args) < 2):
    parser.print_help()
    return
  run(
    args=args,
    cutoff=options.cutoff,
    high_resolution_only=options.high_resolution_only)

if (__name__ == "__main__"):
  main()
