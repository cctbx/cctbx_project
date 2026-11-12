from __future__ import absolute_import, division, print_function

import iotbx.pdb
import mmtbx.model
from mmtbx.validation.ramalyze import ramalyze
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from mmtbx.validation.comparama import add_arrows_on_plot
from mmtbx.validation.comparama import breake_arrow_if_needed

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

# Each case: an arrow (start, end) in (phi, psi) degrees and the segments it
# must be drawn as. An arrow takes the shortest way between its ends on the
# periodic plot; leaving through one edge it comes back through the opposite
# one. A crossing point lies on the straight line from the start to the end
# shifted by 360: (120, 0) -> (-120, 60) runs as (120, 0) -> (240, 60), with
# slope 1/2, so it leaves at phi = 180 with psi = 30. An arrow that crosses
# both edges is drawn through the corner instead.
cases = [
  # no crossing
  (((60, -120), (120, -120)),
    [((60, -120), (120, -120))]),
  (((-125, 120), (-125, 179)),
    [((-125, 120), (-125, 179))]),
  # straight up and down
  (((-120, 120), (-120, -120)),
    [((-120, 120), (-120, 180)), ((-120, -180), (-120, -120))]),
  (((-115, -120), (-115, 120)),
    [((-115, -120), (-115, -180)), ((-115, 180), (-115, 120))]),
  # straight right and left
  (((120, -60), (-120, -60)),
    [((120, -60), (180, -60)), ((-180, -60), (-120, -60))]),
  (((-120, -65), (120, -65)),
    [((-120, -65), (-180, -65)), ((180, -65), (120, -65))]),
  # diagonal right: slope 1/2, psi = 0 + 60/2 at phi = 180
  (((120, 0), (-120, 60)),
    [((120, 0), (180, 30)), ((-180, 30), (-120, 60))]),
  # diagonal left: slope 1/2, psi = 55 - 60/2 at phi = -180
  (((-120, 55), (120, -5)),
    [((-120, 55), (-180, 25)), ((180, 25), (120, -5))]),
  # diagonal up: slope 2, phi = -60 + 60/2 at psi = 180
  (((-60, 120), (0, -120)),
    [((-60, 120), (-30, 180)), ((-30, -180), (0, -120))]),
  # diagonal down: slope 2, phi = 5 - 60/2 at psi = -180
  (((5, -120), (-55, 120)),
    [((5, -120), (-25, -180)), ((-25, 180), (-55, 120))]),
  # through a corner, straight and not straight
  (((150, 150), (-150, -150)),
    [((150, 150), (180, 180)), ((-180, -180), (-150, -150))]),
  (((140, 155), (-130, -140)),
    [((140, 155), (180, 180)), ((-180, -180), (-130, -140))]),
  (((150, -150), (-150, 150)),
    [((150, -150), (180, -180)), ((-180, 180), (-150, 150))]),
  (((140, -155), (-130, 140)),
    [((140, -155), (180, -180)), ((-180, 180), (-130, 140))]),
  (((-150, 150), (150, -150)),
    [((-150, 150), (-180, 180)), ((180, -180), (150, -150))]),
  (((-140, 155), (130, -140)),
    [((-140, 155), (-180, 180)), ((180, -180), (130, -140))]),
  (((-150, -150), (150, 150)),
    [((-150, -150), (-180, -180)), ((180, 180), (150, 150))]),
  (((-140, -155), (130, 140)),
    [((-140, -155), (-180, -180)), ((180, 180), (130, 140))]),
]

def exercise_breake_arrow_if_needed():
  for arrow, expected in cases:
    result = breake_arrow_if_needed(
        arrow[0], arrow[1], plot_ranges=[(-180, 180), (-180, 180)])
    assert approx_equal(result, expected), arrow

def exercise_add_arrows_on_plot():
  import matplotlib.colors
  model = mmtbx.model.manager(
      model_input = iotbx.pdb.input(source_info=None, lines=pdb_str))
  rama = ramalyze(model.get_hierarchy(), out=null_out())
  plot = rama.get_plots()[0]
  n_lines = len(plot.plot.lines)
  n_patches = len(plot.plot.patches)
  add_arrows_on_plot(plot, [arrow for arrow, expected in cases], color="red")
  lines = list(plot.plot.lines)[n_lines:]
  arrows = list(plot.plot.patches)[n_patches:]
  # all segments of an arrow but the last are lines; the last has the head
  segments = []
  for arrow, expected in cases:
    segments.extend(expected[:-1])
  assert len(lines) == len(segments)
  assert len(arrows) == len(cases)
  for line, segment in zip(lines, segments):
    assert approx_equal(line.get_xdata(), [segment[0][0], segment[1][0]])
    assert approx_equal(line.get_ydata(), [segment[0][1], segment[1][1]])
  red = matplotlib.colors.to_rgba("red")
  for line in lines:
    assert matplotlib.colors.to_rgba(line.get_color()) == red
  for a in arrows:
    assert matplotlib.colors.to_rgba(a.get_facecolor()) == red
  # over the filling and contours of the plot
  for artist in lines + arrows:
    assert artist.get_zorder() == 10
  plot.save_image("tst_add_arrows_on_plot_1.png", dpi=100)

if __name__ == '__main__':
  exercise_breake_arrow_if_needed()
  exercise_add_arrows_on_plot()
  print("OK")
