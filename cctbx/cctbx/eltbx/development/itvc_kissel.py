from cctbx.eltbx.development import itvc_section61_io
from cctbx.eltbx.development import kissel_io
import cctbx.eltbx.gaussian_fit
from cctbx.array_family import flex
from libtbx.optparse_wrapper import OptionParser
import os

def run(args, cutoff, high_resolution_only,
        plots_dir="itvc_kissel_plots", verbose=0):
  itab = itvc_section61_io.read_table6111(args[0])
  itab_x = cctbx.eltbx.gaussian_fit.international_tables_stols
  isel = itab_x <= cutoff
  if (high_resolution_only):
    isel &= itab_x > 2 + 1.e-6
  itab_x = itab_x.select(isel)
  for file_name in args[1:]:
    ktab = kissel_io.read_table(file_name)
    if (ktab.element == "Es"): continue
    sel = ktab.x <= cutoff
    ktab_x = ktab.x.select(sel)
    ktab_y = ktab.y.select(sel)
    itab_y = itab.entries[ktab.element].table_y.select(isel)
    f = open(os.path.join(plots_dir, ktab.element + ".xy"), "w")
    cctbx.eltbx.gaussian_fit.write_plot(f, ktab_x, ktab_y)
    cctbx.eltbx.gaussian_fit.write_plot(f, itab_x, itab_y)
    f.close()
    ktab_yi = flex.linear_interpolation(ktab.x, ktab.y, itab_x)
    assert ktab_yi.all_gt(0)
    dyr = flex.double([int(abs(d)*1000)/1000. for d in list(ktab_yi-itab_y)]) \
        / ktab_yi
    print ktab.element, "max_delta=%.4f, relative=%.4f" % (
      flex.max(flex.abs(ktab_yi-itab_y)), flex.max(dyr))
    for x,ky,iy,d in zip(itab_x, ktab_yi, itab_y, dyr):
      print "%4.2f %7.4f %7.4f %7.4f %7.4f" % (
        x,ky,iy,iy-ky, d)
    print

def main():
  parser = OptionParser(
    usage="usage: python %prog [options] itvc_table kissel_files ...")
  parser.add_option("-c", "--cutoff",
    type="float", dest="cutoff", default=6.05, metavar="FLOAT",
    help="maximum sin(theta)/lambda")
  parser.add_option("-g", "--high_resolution_only",
    action="store_true", dest="high_resolution_only", default=0,
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
