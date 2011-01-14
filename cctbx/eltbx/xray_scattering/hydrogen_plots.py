from __future__ import division
from cctbx.eltbx import xray_scattering
import sys

# Table 6.1.1.2. Spherical bonded hydrogen-atom scattering
# factors from Stewart, Davidson & Simpson (1965)
# Manually padded with itvc/v5.txt HF values (similar to what
# the ITC Table 6.1.1.1 authors did in some other cases).
itc_tab_6112_padded = [
  (0.0000, 1.0000),
  (0.0215, 0.9924),
  (0.0429, 0.9704),
  (0.0644, 0.9352),
  (0.0859, 0.8892),
  (0.1073, 0.8350),
  (0.1288, 0.7752),
  (0.1503, 0.7125),
  (0.1718, 0.6492),
  (0.1932, 0.5871),
  (0.2147, 0.5277),
  (0.2576, 0.4201),
  (0.3006, 0.3301),
  (0.3435, 0.2573),
  (0.3864, 0.1998),
  (0.4294, 0.1552),
  (0.4723, 0.1208),
  (0.5153, 0.0945),
  (0.5582, 0.0744),
  (0.6011, 0.0592),
  (0.6441, 0.0474),
  (0.6870, 0.0383),
  (0.7300, 0.0311),
  (0.7729, 0.0254),
  (0.8158, 0.0208),
  (0.8588, 0.0171),
  (0.9017, 0.0140),
  (0.9447, 0.0116),
  (0.9876, 0.0096),
  (1.0305, 0.0080),
  (1.0735, 0.0066),
  (1.1164, 0.0056),
  (1.1593, 0.0047),
  (1.2023, 0.0040),
  (1.2452, 0.0035),
  (1.2882, 0.0031),
  (1.3311, 0.0027),
  (1.3740, 0.0025),
  (1.4170, 0.0022),
  (1.4599, 0.0020),
  (1.5029, 0.0018),
  (1.5458, 0.0016),
  (1.5887, 0.0015),
  (1.6317, 0.0013),
  (1.6746, 0.0011),
  (1.7176, 0.0010),
  # Manually padded with itvc/v5.txt HF values.
  (1.80, 0.0008),
  (1.90, 0.0006),
  (2.00, 0.0005),
  (2.50, 0.0002),
  (3.00, 0.0001),
  (3.50, 0.0001),
  (4.00, 0.0000),
  (5.00, 0.0000),
  (6.00, 0.0000)]

def run(args):
  assert len(args) == 0
  sds = xray_scattering.gaussian( # ITC 1992 Table 6.1.1.4
    (0.493002, 0.322912, 0.140191, 0.040810),
    (10.5109, 26.1257, 3.14236, 57.7997),
    0.003038)
  hf = xray_scattering.it1992("H").fetch()
  wk = xray_scattering.wk1995("H").fetch()
  ng = xray_scattering.n_gaussian_table_entry("H", 6).gaussian()
  print "@with g0"
  print '@ s0 legend "ITC Tab 6.1.1.2"'
  for i,lbl in enumerate(["SDS", "HF", "WK", "NG"]):
    print '@ s%d legend "%s"' % (i+1, lbl)
  print "@ s0 symbol 1"
  print "@ s0 line linestyle 0"
  for x,y in itc_tab_6112_padded:
    print x, y
  print "&"
  n_samples = 100
  for g in [sds, hf, wk, ng]:
    for i_stol in xrange(n_samples+1):
      stol = 6 * i_stol / n_samples
      print stol, g.at_stol(stol)
    print "&"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
