from cctbx.eltbx import xray_scattering
from cctbx.eltbx.e_scattering import ito_vol_c_table_4_3_2_2
import sys

def run(args):
  assert len(args) == 0
  assert len(ito_vol_c_table_4_3_2_2) == 98
  gaussians = []
  labels_found = set()
  for line in ito_vol_c_table_4_3_2_2:
    flds = line.split()
    assert len(flds) == 12
    std_lbl = xray_scattering.get_standard_label(flds[0], exact=True)
    assert flds[0] == std_lbl
    assert std_lbl not in labels_found
    labels_found.add(std_lbl)
    assert flds[1] == str(len(labels_found))
    def vals(i,j): return [float(s) for s in flds[i:j]]
    array_of_a = vals(2,7)
    array_of_b = vals(7,12)
    g = xray_scattering.gaussian(array_of_a, array_of_b)
    gaussians.append((std_lbl, g))
  #
  n_samples = 1000
  from matplotlib.backends.backend_pdf import PdfPages
  all_pdf = PdfPages("all.pdf")
  from libtbx import pyplot
  from scitbx.array_family import flex
  for std_lbl,g in gaussians:
    fig = pyplot.figure()
    fig.set_size_inches(11, 8.5)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(std_lbl, fontsize=12)
    def one_curv(g, code):
      x = flex.double()
      y = flex.double()
      for i_stol in xrange(n_samples+1):
        stol = 6 * i_stol / n_samples
        x.append(stol)
        y.append(g.at_stol(stol))
      ax.plot(x.as_numpy_array(), y.as_numpy_array(), code)
    one_curv(g, "b-")
    one_curv(xray_scattering.gaussian(
      g.array_of_a()[:4],
      g.array_of_b()[:4]), "r-")
    all_pdf.savefig(fig, bbox_inches="tight")
  all_pdf.close()
  print "plots written to file: all.pdf"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
