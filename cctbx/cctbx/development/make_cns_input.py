from cctbx_boost import adptbx

def topology(elements):
  cns_input = []
  l = cns_input.append
  ElemDict = {}
  for Elem in elements: ElemDict[Elem] = 0
  l("topology")
  for Elem in ElemDict.keys():
    l("  residue %s" % (Elem,))
    l("    atom %s mass=1 charge=0 {chemical}type=%s end" % (Elem, Elem))
    l("  end")
  l("end")
  l("")
  l("segment")
  l("  name=A")
  for Elem in elements:
    l("  molecule {res}name=%s number=1 end" % (Elem,))
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
    for i in xrange(3):
      l("do (%s=%.12g) (resid=%d)" % ("xyz"[i], X[i], resid))
    l("do (q=%.12g) (resid=%d)" % (q, resid))
    l("do (b=%.12g) (resid=%d)" % (b, resid))
    resid += 1
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
