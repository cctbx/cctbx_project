from __future__ import absolute_import, division, print_function
import sys
from six.moves import zip

items = [
  "from_scatterers_direct",
  "from_scatterers_fft",
  "gradients_direct",
  "gradients_fft"]

class calls_and_accumulated_time(object):

  def __init__(self):
    self.calls = 0
    self.time = 0

  def process(self, time):
    self.calls += 1
    self.time += time

for item in items:
  exec("%s = calls_and_accumulated_time()" % item)
del item

def show(out=None, prefix="", show_zero_calls=False):
  if (out is None): out = sys.stdout
  g = globals()
  calls = ["%d" % c for c in [g[item].calls for item in items]]
  times = ["%.2f" % t for t in [g[item].time for item in items]]
  fmt = "%%-23s %%%ds calls, %%%ds s" % (
    max([len(s) for s in calls]),
    max([len(s) for s in times]))
  for i,c,t in zip(items, calls, times):
    if (not show_zero_calls and c == "0"): continue
    print(prefix + (fmt % (i+":", c, t)).replace(
      " 1 calls,", " 1 call, "), file=out)
