from cctbx_boost import adptbx
from cctbx_boost.eltbx.caasf_it1992 import CAASF_IT1992

def topology(sites):
  cns_input = []
  l = cns_input.append
  l("topology")
  for site in sites:
    lbl = site.Label()
    l("  residue %s" % (lbl,))
    l("    atom %s mass=1 charge=0 {chemical}type=%s end" % (lbl, lbl))
    l("  end")
  l("end")
  l("")
  l("segment")
  l("  name=A")
  for site in sites:
    lbl = site.Label()
    l("  molecule {res}name=%s number=1 end" % (lbl,))
  l("end")
  l("")
  return cns_input

def coordinates(Sites):
  cns_input = []
  l = cns_input.append
  resid = 1
  for Site in Sites:
    X = Site.Coordinates()
    q = Site.Occ()
    b = adptbx.U_as_B(Site.Uiso())
    sf = CAASF_IT1992(Site.CAASF().Label())
    fp = Site.fpfdp().real
    fdp = Site.fpfdp().imag
    for i in xrange(3):
      l("do (%s=%.12g) (resid=%d)" % ("xyz"[i], X[i], resid))
    l("do (q=%.12g) (resid=%d)" % (q, resid))
    l("do (b=%.12g) (resid=%d)" % (b, resid))
    l("xray")
    l("  scatter (chemical=%s)" % (Site.Label(),))
    for i in xrange(4):
      l("    %.6g %.6g" % (sf.a(i), sf.b(i)))
    l("    %.6g" % (sf.c(),))
    l("end")
    l("do (scatter_fp=%.12g) (resid=%d)" % (fp, resid))
    l("do (scatter_fdp=%.12g) (resid=%d)" % (fdp, resid))
    l("")
    resid += 1
  l("")
  l("show (name) (all)")
  l("show (chemical) (all)")
  l("show (scatter_fp) (all)")
  l("show (scatter_fdp) (all)")
  l("")
  return cns_input

def unit_cell(unit_cell):
  cns_input = []
  l = cns_input.append
  l("xray")
  l("  a=%.12g b=%.12g c=%.12g alpha=%.12g beta=%.12g gamma=%.12g"
    % unit_cell.getParameters())
  l("end")
  l("")
  return cns_input

def symmetry(space_group_ops):
  cns_input = []
  l = cns_input.append
  l("xray")
  for m in space_group_ops:
    mp = m.modPositive()
    l("  symm=(" + mp.as_xyz() + ")")
  l("end")
  l("")
  return cns_input

def predict(reciprocal_space_array_name, method):
  cns_input = []
  l = cns_input.append
  l("""xray
  declare name=%s domain=reciprocal type=complex end
  method=%s
  predict
    mode=reciprocal
    to=%s
    selection=(all)
    atomselection=(all)
  end
end
"""
% (reciprocal_space_array_name, method, reciprocal_space_array_name))
  return cns_input
