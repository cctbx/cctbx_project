from cctbx.eltbx import xray_scattering
from cctbx.eltbx import tiny_pse
from cctbx.array_family import flex
import cctbx.eltbx.gaussian_fit
import scitbx.math.gaussian
from libtbx.str_utils import line_feeder
from libtbx.option_parser import OptionParser
from libtbx import adopt_init_args

class table(object):

  def __init__(self, atomic_number, x, y, sigmas):
    adopt_init_args(self, locals())
    self.element = tiny_pse.table(atomic_number).symbol()

  def itvc_sampling_selection(self):
    xi = cctbx.eltbx.gaussian_fit.international_tables_stols
    xk = self.x
    selection = flex.bool(xk.size(), False)
    i_kissel = 0
    for i_itvc in xrange(xi.size()):
      while (xk[i_kissel] < xi[i_itvc]):
        i_kissel += 1
      if (xk[i_kissel] == xi[i_itvc]):
        selection[i_kissel] = True
      elif (i_kissel > 0 and xk[i_kissel-1] < xi[i_itvc] < xk[i_kissel]):
        if (xi[i_itvc] - xk[i_kissel-1] < xk[i_kissel] - xi[i_itvc]):
          selection[i_kissel-1] = True
        else:
          selection[i_kissel] = True
    return selection

def read_table(file_name):
  atomic_number = None
  number_of_electrons = None
  x = flex.double()
  y = flex.double()
  sigmas = flex.double()
  lf = line_feeder(open(file_name))
  while 1:
    line = lf.next()
    if (lf.eof): break
    if (line.startswith("   FORM: ATOMIC NUMBER=")):
      atomic_number = float(line.split("=")[1])
      assert int(atomic_number) == atomic_number
    elif (line.startswith("   FORM: # ELECTRONS=")):
      number_of_electrons = float(line.split("=")[1])
      assert int(number_of_electrons) == number_of_electrons
    elif (line.startswith("        X (1/A)")):
      assert atomic_number == number_of_electrons
      while 1:
        line = lf.next()
        assert not lf.eof
        if (line == " *** END OF DATA ***"):
          lf.eof = True
          break
        vals_str = line.split()
        for val_str in vals_str: assert len(val_str) == 13
        x.append(float(vals_str[0]))
        y.append(float(vals_str[1]))
        assert vals_str[1][-4] == "E"
        sigmas.append(float("0.00000005"+vals_str[1][-4:]))
  return table(int(atomic_number), x, y*atomic_number, sigmas)

def main():
  parser = OptionParser(
    usage="usage: python %prog [options] file_name ...")
  parser.add_option("-c", "--cutoff",
    type="float", default=6.05, metavar="FLOAT",
    help="maximum sin(theta)/lambda")
  (options, args) = parser.parse_args()
  if (len(args) < 1):
    parser.print_help()
    return
  cutoff = options.cutoff
  for file_name in args:
    tab = read_table(file_name)
    if (tab.element == "Es"): continue
    wk = xray_scattering.wk1995(tab.element, True).fetch()
    sel = tab.x <= cutoff
    tab_x = tab.x.select(sel)
    tab_y = tab.y.select(sel)
    sigmas = flex.double(tab_x.size(), 0.0005)
    wky = wk.at_x(tab_x)
    errors_abs = flex.abs(wky-tab_y)
    fit = scitbx.math.gaussian.fit(tab_x, tab_y, sigmas, wk)
    errors_rel = fit.significant_relative_errors(1.e-6)
    print tab.element, tab.atomic_number,
    print "max error < %.1fA-1 abs, rel: %7.4f %7.4f" % (
      cutoff, flex.max(errors_abs), flex.max(errors_rel))
    for x,y,f,ea,er in zip(tab_x,tab_y,wky,errors_abs,errors_rel):
      print "%7.4f %7.4f %7.4f %7.4f %7.4f" % (x, y, f, ea, er)
    print

if (__name__ == "__main__"):
  main()
