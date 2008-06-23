from __future__ import division
from cctbx import crystal, xray
from smtbx import refinement

def exercise_parameter_map():
  cs = crystal.symmetry((8,9,10, 85, 95, 105), "P1")
  xs = xray.structure(cs.special_position_settings())
  for i in xrange(5): xs.add_scatterer(xray.scatterer("C%i" % i))
  grad_site      = (True , False, False, True , False)
  grad_u_iso     = (False, True , True , False, True )
  grad_u_aniso   = (True , False, False, True , True )
  grad_occupancy = (False, True , True , False, False)
  grad_fp        = (True , True , False, False, False)
  grad_fdp       = (True , False, True , False, False)
  for sc, site, u_iso, u_aniso, occ, fp, fdp in zip(xs.scatterers(),
    grad_site, grad_u_iso, grad_u_aniso, grad_occupancy, grad_fp, grad_fdp):
    f = sc.flags
    f.set_grad_site(site)
    f.set_grad_u_iso(u_iso)
    f.set_grad_u_aniso(u_aniso)
    f.set_grad_occupancy(occ)
    f.set_grad_fp(fp)
    f.set_grad_fdp(fdp)
  m = xs.parameter_map()
  assert m.n_parameters == xs.n_parameters()

  indices = m[0]
  assert indices.site == 0
  assert indices.u_iso == -1
  assert indices.u_aniso == 3
  assert indices.occupancy == -1
  assert indices.fp == 9
  assert indices.fdp == 10

  indices = m[1]
  assert indices.site == -1
  assert indices.u_iso == 11
  assert indices.u_aniso == -1
  assert indices.occupancy == 12
  assert indices.fp == 13
  assert indices.fdp == -1

def run():

  print 'OK'

if __name__ == '__main__':
  run()
