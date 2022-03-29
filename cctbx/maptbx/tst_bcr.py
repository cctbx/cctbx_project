from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from scitbx.array_family import flex
from cctbx.maptbx.bcr import bcr

def rfactor(a,b):
  n = flex.sum(flex.abs(a-b))
  d = flex.sum(flex.abs(a+b))
  return n/d*100*2

def run(d_min       = 2.0,
        radius_max  = 5,
        radius_step = 0.01):
  for i, mxp in enumerate([0,6]):
    o = maptbx.atom_curves(scattering_type="C", scattering_table="wk1995")
    b = o.bcr_approx(
      d_min       = d_min,
      radius_max  = radius_max,
      radius_step = radius_step,
      mxp=mxp, epsc=0.001, kpres=0 # BCR params
      )
    r = rfactor(b.image_values, b.bcr_approx_values)
    if(i == 0): assert r>10, r
    else:       assert r<0.4, r
    #
    # b_iso is not 0
    #
    b_iso = 50
    im = o.image(
      d_min=d_min, b_iso=b_iso, radius_max=radius_max, radius_step=radius_step)
    bcr_approx_values = flex.double()
    for r in im.radii:
      first = 0
      second = 0
      for B, C, R in zip(b.B, b.C, b.R):
        if(abs(R)<1.e-6):
          first += bcr.gauss(B=B, C=C, r=r, b_iso=b_iso)
        else:
          second += C*bcr.chi(B=B, R=R, r=r, b_iso=b_iso)
      bcr_approx_values.append(first + second)
    r = rfactor(im.image_values, bcr_approx_values)
    if(i == 0): assert r>10, r
    else:       assert r<0.4, r

if (__name__ == "__main__"):
  run()
