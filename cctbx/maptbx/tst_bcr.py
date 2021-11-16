from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from scitbx.array_family import flex

def rfactor(a,b):
  n = flex.sum(flex.abs(a-b))
  d = flex.sum(flex.abs(a+b))
  return n/d*100*2

def run():
  for i, mxp in enumerate([0,5]):
    o = maptbx.atom_curves(scattering_type="C", scattering_table="wk1995")
    b = o.bcr_approx(
      d_min       = 2.0,
      b_iso       = 0,
      radius_max  = 5,
      radius_step = 0.01,
      mxp=mxp, epsc=0.001, kpres=0 # BCR params
      )
    r = rfactor(b.image_values,b.bcr_approx_values)
    if(i == 0): assert r>15
    else:       assert r<0.2

if (__name__ == "__main__"):
  run()
