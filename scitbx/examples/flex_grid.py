"""\
grid = grid over N-dim box (but there is no array with data!)

origin = N-dim indices of lower-left corner of box
 focus = N-dim indices of upper-right corner of box as we want to use it
  last = N-dim indices of upper-right corner of box as actually allocated
   all = actual number or grid points to be allocated in each dimension

Motivation for focus/last distinction:
  padding required by real-to-complex FFT algorithms
"""
from __future__ import absolute_import, division, print_function

from scitbx.array_family.flex import grid

def show(g):
  print("origin:", g.origin())
  print(" focus:", g.focus())
  print("  last:", g.last())
  print("   all:", g.all())
  print()

def run():
  print(__doc__)
  g = grid((3,4,6))
  show(g)
  g.set_focus((3,4,5))
  show(g)
  g = grid((-2,-3,-4), (3,4,6))
  show(g)
  g.set_focus((3,4,5))
  show(g)
  print("OK")

if (__name__ == "__main__"):
  run()
