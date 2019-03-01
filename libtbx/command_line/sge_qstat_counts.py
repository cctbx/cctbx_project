from __future__ import absolute_import, division, print_function
from operator import itemgetter
from libtbx.queuing_system_utils.sge_utils import qstat_parse
from libtbx import dict_with_default_0
import sys

def run(args):
  assert len(args) == 0
  qstat_info = qstat_parse()
  user_states = {}
  for items in qstat_info:
    counts = user_states.setdefault(items.user, dict_with_default_0())
    counts[items.state] += items.counts()
  sum_counts = []
  for user,counts in user_states.items():
    sum_counts.append((user, sum(counts.values())))
  sum_counts.sort(key=itemgetter(1, 0))
  cpus=0
  for user,sc in sum_counts:
    counts = user_states[user]
    print("%-10s   %5d r   %5d qw" % (user, counts["r"], counts["qw"]), end=' ')
    cpus+=counts["r"]
    for state,c in counts.items():
      if (state in ["r", "qw"]): continue
      print("  %5d %s" % (c, state), end=' ')
    print("  %5d total" % sc)
  print("-"*45)
  print("%-10s   %5d r" % ("total", cpus))
  print("-"*45)

if (__name__ == "__main__"):
  run(sys.argv[1:])
