from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx.stdlib import math
from six.moves import range

def rosca(m=9, hemisphere=True):
  """
  Regular grid on the unit sphere, Rosca (2010).
  """
  def truncate(x):
    if(abs(x)<1.e-6): return 0
    else:             return x
  def add(result, new):
    foud = False
    for s in result:
      d = math.sqrt((s[0]-new[0])**2+(s[1]-new[1])**2+(s[2]-new[2])**2)
      if(d<1.e-3): return
    result.append(new)
    return
  d_l = math.sqrt(2)/m
  result = flex.vec3_double()
  for m_ in range(m+1):
    if(m_==0):
      add(result, [0,0,1])
    else:
      l_m = m_ * d_l
      d_phi = math.pi/(4*m_)
      assert l_m>=0 and l_m<=math.sqrt(2)
      for k in range(m_+1):
        arg1 = k*d_phi
        arg2 = l_m*math.sqrt(1-l_m**2/4)
        x_km = truncate(math.cos(arg1)*arg2)
        y_km = truncate(math.sin(arg1)*arg2)
        z_m  = truncate(1.-l_m**2/2                                )
        add(result, [ x_km, y_km,z_m])
        add(result, [ y_km, x_km,z_m])
        add(result, [-y_km, x_km,z_m])
        add(result, [-x_km, y_km,z_m])
        add(result, [ x_km,-y_km,z_m])
        add(result, [ y_km,-x_km,z_m])
        add(result, [-x_km,-y_km,z_m])
        add(result, [-y_km,-x_km,z_m])
        if(not hemisphere):
          add(result, [ x_km, y_km,-z_m])
          add(result, [ y_km, x_km,-z_m])
          add(result, [-y_km, x_km,-z_m])
          add(result, [-x_km, y_km,-z_m])
          add(result, [ x_km,-y_km,-z_m])
          add(result, [ y_km,-x_km,-z_m])
          add(result, [-x_km,-y_km,-z_m])
          add(result, [-y_km,-x_km,-z_m])
  for r in result:
    assert abs(1.-math.sqrt(r[0]**2+r[1]**2+r[2]**2))<1.e-6
  # XXX for debugging
  if(0):
    f = "HETATM%5d  O   HOH A%4d    %8.3f%8.3f%8.3f  1.00 23.99           O "
    o = open("junk.pdb", "w")
    for i, r in enumerate(result):
      print(f%(i, i, r[0],r[1],r[2]), file=o)
    o.close()
  return result
