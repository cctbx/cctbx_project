from __future__ import absolute_import, division, print_function
from cctbx.eltbx import e_scattering
import sys
from six.moves import range

def run(args):
  assert len(args) == 0
  n_samples = 1000
  from matplotlib.backends.backend_pdf import PdfPages
  all_pdf = PdfPages("all.pdf")
  from libtbx import pyplot
  from scitbx.array_family import flex
  for element in e_scattering.ito_vol_c_2011_table_4_3_2_2_elements():
    fig = pyplot.figure()
    fig.set_size_inches(11, 8.5)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(element, fontsize=12)
    def one_curv(g, code):
      x = flex.double()
      y = flex.double()
      for i_stol in range(n_samples+1):
        stol = 6 * i_stol / n_samples
        x.append(stol)
        y.append(g.at_stol(stol))
      ax.plot(x.as_numpy_array(), y.as_numpy_array(), code)
    g = e_scattering.ito_vol_c_2011_table_4_3_2_2_entry_as_gaussian(
      label=element, exact=True)
    one_curv(g, "b-")
    one_curv(e_scattering.gaussian(
      g.array_of_a()[:4],
      g.array_of_b()[:4]), "r-")
    all_pdf.savefig(fig, bbox_inches="tight")
  all_pdf.close()
  print("plots written to file: all.pdf")

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
