from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from cctbx import maptbx
from cctbx.maptbx.bcr import qmap

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

def rfactor(a,b):
  n = flex.sum(flex.abs(a-b))
  d = flex.sum(flex.abs(a+b))
  return n/d*100*2, flex.max(flex.abs(a-b))

def run():
  if 0:
    o = maptbx.atom_curves(scattering_type="S", scattering_table="wk1995")
    b = o.bcr_approx(
      d_min       = 0.41,
      radius_max  = 14.,
      radius_step = 0.01,
      mxp   = 1000,
      epsc  = 0.001,
      epsp  = 0.000,
      edist = 1.0E-13,
      kpres = 1,
      kprot = 112)
    #print("B, C, R")
    #for B,C,R in zip(b.B, b.C, b.R):
    #  print("%12.6f %12.6f %12.6f"%(B,C,R))
    #print()
    #print("dist, image, approx, image-approx")
    for a,b,c in zip(
      b.radii, b.image_values, b.bcr_approx_values):
      print("%5.3f %15.9f %15.9f %12.6f"%(a,b,c, b-c))
  #
  if 1:
    tH = qmap.load_table(element="H", table="wk1995")
    tC = qmap.load_table(element="C", table="wk1995")
    tO = qmap.load_table(element="O", table="wk1995")
    tN = qmap.load_table(element="N", table="wk1995")
    tS = qmap.load_table(element="S", table="wk1995")
    f = "%4s %8.6f %8.6f %8.6f %8.6f %8.6f"
    for k in tC.keys():
      print(f%(k, tH[k]['err'], tC[k]['err'], tO[k]['err'], tN[k]['err'], tS[k]['err']))

if (__name__ == "__main__"):
  run()
