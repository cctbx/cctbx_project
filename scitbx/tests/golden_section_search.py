from __future__ import absolute_import, division, print_function
import scitbx.golden_section_search

def function(x):
  y = (x-3.5)**2.0
  return y

def tst_it():
  opti = scitbx.golden_section_search.gss(function,-10,10)
  assert abs(opti-3.5)<1e-4

if __name__ == "__main__":
  tst_it()
