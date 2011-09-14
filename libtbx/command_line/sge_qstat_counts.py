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
  def cmp_sum_counts(a,b):
    result = cmp(b[1], a[1])
    if (result != 0): return result
    return cmp(b[0], a[0])
  sum_counts.sort(cmp_sum_counts)
  for user,sc in sum_counts:
    counts = user_states[user]
    print "%-10s   %5d r   %5d qw" % (user, counts["r"], counts["qw"]),
    for state,c in counts.items():
      if (state in ["r", "qw"]): continue
      print "  %5d %s" % (c, state),
    print "  %5d total" % sc

if (__name__ == "__main__"):
  run(sys.argv[1:])
