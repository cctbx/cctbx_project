"""Creates d*trek <-> cctbx space group symbol dictionary based on
   d*trek symmetry file (DTREK_SPACEGROUP_FILE, DTREK73 distribution)

   usage: python dtrek_symmetry_dict.py DTREK_SPACEGROUP_FILE
"""
from __future__ import absolute_import, division, print_function

from libtbx.str_utils import line_feeder
from cctbx import sgtbx
import sys
import pprint
from six.moves import range

class dtrek_symmetry_entry(object):

  def __init__(self, lf):
    self.symbol = None
    l = next(lf)
    print(l.strip())
    if (lf.eof): return
    no, number, symbol, m = l.split()
    assert no == "NO."
    number = int(number)
    assert int(float(m)) == float(m)
    m = int(float(m))
    space_group = sgtbx.space_group()
    if (m < 0):
      space_group.expand_inv(sgtbx.tr_vec((0,0,0)))
    space_group.expand_conventional_centring_type(symbol[0])
    t_den = space_group.t_den()
    matrices = []
    for i in range(abs(m)):
      l = next(lf)
      assert not lf.eof
      flds = l.split()
      assert len(flds) == 12
      r = [int(e) for e in flds[1:4] + flds[5:8] + flds[9:12]]
      t = [int(round(float(e)*t_den)) for e in (flds[0],flds[4],flds[8])]
      try:
        s = sgtbx.rt_mx(sgtbx.rot_mx(r), sgtbx.tr_vec(t))
      except RuntimeError as e:
        print(e)
        print(l)
      else:
        try:
          matrices.append(s)
          space_group.expand_smx(s)
        except RuntimeError as e:
          print(e)
          print(l)
          print(s)
    space_group_info = sgtbx.space_group_info(group=space_group)
    if (space_group_info.type().number() != number):
      print("Space group number mismatch:")
      print("   from file:", number)
      print("  operations:", space_group_info.type().number())
      for s in matrices:
        space_group = sgtbx.space_group_info(symbol=number).group()
        order_z = space_group.order_z()
        space_group.expand_smx(s)
        if (space_group.order_z() != order_z):
          print("  misfit:", s)
      space_group = sgtbx.space_group_info(symbol=number).group()
      for i_smx in range(space_group.n_smx()):
        OK = False
        for i_inv in range(space_group.f_inv()):
          for i_ltr in range(space_group.n_ltr()):
            sg = space_group(i_ltr, i_inv, i_smx).mod_positive()
            for sm in matrices:
              sm = sm.mod_positive()
              if (sm == sg):
                OK = True
                break
            if (OK): break
          if (OK): break
        if (not OK):
          print("  missing:", sg)
    self.number = number
    self.symbol = symbol
    self.space_group_info = space_group_info

def run():
  assert len(sys.argv) == 2
  f = open(sys.argv[1], "r")
  lf = line_feeder(f)
  lookup_dict = {}
  while 1:
    entry = dtrek_symmetry_entry(lf)
    if (entry.symbol is None):
      break
    cctbx_lookup_symbol = str(entry.space_group_info).replace(" ", "")
    if (cctbx_lookup_symbol != entry.symbol):
      lookup_dict[entry.symbol] = cctbx_lookup_symbol
  pprint.pprint(lookup_dict)

if (__name__ == "__main__"):
  run()
