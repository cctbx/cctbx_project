from cctbx import adptbx
from cctbx.eltbx.caasf import it1992
import sys, os

def check_cns_availability():
  if (not os.environ.has_key("CNS_INST")):
    print "CNS not available."
    sys.exit(1)

def topology(structure):
  cns_input = []
  l = cns_input.append
  l("topology")
  for atom in structure:
    lbl = atom.label()
    l("  residue %s" % (lbl,))
    l("    atom %s mass=1 charge=0 {chemical}type=%s end" % (lbl, lbl))
    l("  end")
  l("end")
  l("")
  l("segment")
  l("  name=A")
  for atom in structure:
    lbl = atom.label()
    l("  molecule {res}name=%s number=1 end" % (lbl,))
  l("end")
  l("")
  return cns_input

def coordinates(xray_structure):
  cns_input = []
  l = cns_input.append
  resid = 1
  for scatterer in xray_structure.scatterers():
    x = scatterer.site
    q = scatterer.occupancy
    assert not scatterer.anisotropic_flag
    b = adptbx.u_as_b(scatterer.u_iso)
    caasf = it1992(scatterer.caasf.label())
    fp = scatterer.fp_fdp.real
    fdp = scatterer.fp_fdp.imag
    for i in xrange(3):
      l("do (%s=%.12g) (resid=%d)" % ("xyz"[i], x[i], resid))
    l("do (q=%.12g) (resid=%d)" % (q, resid))
    l("do (b=%.12g) (resid=%d)" % (b, resid))
    l("xray")
    l("  scatter (chemical=%s)" % (scatterer.label,))
    for i in xrange(4):
      l("    %.6g %.6g" % (caasf.a(i), caasf.b(i)))
    l("    %.6g" % (caasf.c(),))
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

def xray_unit_cell(unit_cell):
  cns_input = []
  l = cns_input.append
  l("xray")
  l("  a=%.12g b=%.12g c=%.12g alpha=%.12g beta=%.12g gamma=%.12g"
    % unit_cell.parameters())
  l("end")
  l("")
  return cns_input

def xray_symmetry(space_group):
  cns_input = []
  l = cns_input.append
  l("xray")
  for m in space_group:
    mp = m.mod_positive()
    l("  symm=(" + mp.as_xyz() + ")")
  l("end")
  l("")
  return cns_input

def xray_anomalous(flag):
  cns_input = []
  l = cns_input.append
  l("xray")
  if (flag):
    l("  anomalous=true")
  else:
    l("  anomalous=false")
  l("end")
  l("")
  return cns_input

def xray_generate(d_max, d_min):
  cns_input = []
  l = cns_input.append
  l("xray")
  l("  generate %.6g %.6g" % (d_max, d_min))
  l("end")
  l("")
  return cns_input

def xray_predict(reciprocal_space_array_name, method):
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
