from __future__ import absolute_import, division, print_function
import iotbx.pdb
from mmtbx import real_space

pdb_str = """\
CRYST1  115.062  115.062   90.199  90.00  90.00 120.00 P 62
HETATM 3333 ZN    ZN E   1     -26.232  76.952   9.703  1.00 49.37          ZN
HETATM 3334 ZN    ZN E   2     -31.336  69.295  -1.147  0.00 11.59          ZN
TER
"""

def run():
  xrs = iotbx.pdb.input(source_info=None, lines=pdb_str).xray_structure_simple()
  smd = real_space.sampled_model_density(xray_structure=xrs, grid_step=0.3)

if __name__ == '__main__':
  run()
