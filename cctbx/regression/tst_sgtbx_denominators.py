from __future__ import absolute_import, division, print_function
from cctbx import sgtbx
from cctbx.development import debug_utils
import sys

def run_call_back(flags, space_group_info):
  t_den = sgtbx.sg_t_den * 2
  g = sgtbx.space_group()
  g.reset(t_den)
  for s in space_group_info.group():
    g.expand_smx(s.new_denominators(1,t_den))
  info = sgtbx.space_group_info(group=g)
  for fr in (1,2):
    for ft in (1,2):
      r = sgtbx.cb_r_den * fr
      t = sgtbx.cb_t_den * ft
      type = info.type(r_den=r, t_den=t)
      assert type.number() == space_group_info.type().number()
      assert type.lookup_symbol() == space_group_info.type().lookup_symbol()
      assert type.cb_op().c().r().den() == r
      assert type.cb_op().c().t().den() == t

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print("OK")

if (__name__ == "__main__"):
  run()
