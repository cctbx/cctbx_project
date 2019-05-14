from __future__ import absolute_import, division, print_function

import iotbx.pdb
import mmtbx.model
from mmtbx.validation.ramalyze import ramalyze
from libtbx.utils import null_out
from mmtbx.validation.comparama import add_arrows_on_plot
import hashlib

pdb_str = """\
ATOM     13  CB  PHE A   1      11.914  10.410  11.811  1.00  2.00           C
ATOM     14  CG  PHE A   1      11.204   9.472  12.746  1.00  2.00           C
ATOM     15  CD1 PHE A   1      10.636   8.301  12.273  1.00  2.00           C
ATOM     16  CD2 PHE A   1      11.105   9.762  14.096  1.00  2.00           C
ATOM     17  CE1 PHE A   1       9.982   7.436  13.131  1.00  2.00           C
ATOM     18  CE2 PHE A   1      10.452   8.901  14.958  1.00  2.00           C
ATOM     19  CZ  PHE A   1       9.890   7.737  14.475  1.00  2.00           C
ATOM     20  C   PHE A   1      11.828  12.443  10.351  1.00  2.00           C
ATOM     21  O   PHE A   1      11.808  12.365   9.123  1.00  2.00           O
ATOM     22  OXT PHE A   1      12.531  13.314  10.864  1.00  2.00           O
ATOM     23  N   PHE A   1       9.929  10.880  10.444  1.00  2.00           N
ATOM     24  CA  PHE A   1      11.008  11.488  11.213  1.00  2.00           C
"""

reference_md5 = "10bd3c50dc5857a45df9ed698fe69eb4"

def exercise_1(prefix="tst_add_arrows_on_plot_1"):
  model = mmtbx.model.manager(
      model_input = iotbx.pdb.input(source_info=None, lines=pdb_str))
  rama = ramalyze(model.get_hierarchy(), out=null_out())
  plots = rama.get_plots(
      show_labels=True,
      point_style='bo',
      markersize=1,
      markeredgecolor="black",
      dpi=300,
      markerfacecolor="white")

  ad_testing = []
  ad_testing.append( ((60,-120), (120, -120)) )
  ad_testing.append( ((-125, 120), (-125,  179)) )
  ad_testing.append( ((-120, 120), (-120, -120)) ) # wrapping up
  ad_testing.append( ((-115, -120), (-115, 120)) ) # wrapping down
  ad_testing.append( ((120, -60), (-120, -60)) ) # wrapping right
  ad_testing.append( ((-120, -65), (120, -65)) ) # wrapping left
  ad_testing.append( ((120, 0), (-120, 60)) ) # diag right
  ad_testing.append( ((-120, 55), (120, -5)) )# diag left
  ad_testing.append( ((-60, 120), (0, -120)) ) # diag up
  ad_testing.append( ((5, -120), (-55, 120)) ) # diag up
  ad_testing.append( ((150, 150), (-150, -150)) ) # going to top right corner straight
  ad_testing.append( ((140, 155), (-130, -140)) ) # going to top right corner not straight
  ad_testing.append( ((150, -150), (-150, 150)) ) # going to bottom right corner straight
  ad_testing.append( ((140, -155), (-130, 140)) ) # going to bottom right corner not straight
  ad_testing.append( ((-150, 150), (150, -150)) ) # going to top left corner straight
  ad_testing.append( ((-140, 155), (130, -140)) ) # going to top left corner not straight
  ad_testing.append( ((-150, -150), (150, 150)) ) # going to bottom left corner straight
  ad_testing.append( ((-140, -155), (130, 140)) ) # going to bottom left corner not straight

  plot = plots[0]
  add_arrows_on_plot(
      plot,
      ad_testing,
      color="black")
  plot_file_name = "%s.png" % prefix
  plot.save_image(plot_file_name, dpi=300)

  hasher = hashlib.md5()
  with open(plot_file_name, 'rb') as afile:
    buf = afile.read()
    hasher.update(buf)
  fhash = hasher.hexdigest()
  assert fhash == reference_md5, "%s != %s" % (fhash, reference_md5)

if __name__ == '__main__':
  exercise_1()
