from __future__ import absolute_import, division, print_function
import time
from cctbx import xray, crystal, adptbx
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx import sgtbx

"""
xray.structure.is_similar: what makes two structures different, and that the
comparison is fast for large structures (it is used as a consistency check
after every step of real-space refinement).
"""

def get(n = 50, seed = 0):
  import random
  random.seed(seed)
  flex.set_random_seed(seed)
  return random_structure.xray_structure(
    space_group_info = sgtbx.space_group_info("P212121"),
    elements         = ["C","N","O","S"]*(n//4),
    unit_cell        = (30, 40, 50, 90, 90, 90),
    use_u_iso        = True,
    random_u_iso     = True,
    random_occupancy = True)

def exercise_what_differs():
  xrs = get()
  assert xrs.is_similar(xrs)
  assert xrs.is_similar(xrs.deep_copy_scatterers())
  eps = 1.e-6
  # sites
  for shift, expected in [(0.5*eps, True), (2*eps, False)]:
    o = xrs.deep_copy_scatterers()
    sites = o.sites_frac()
    x = list(sites[7]); x[1] += shift; sites[7] = x
    o.set_sites_frac(sites)
    assert xrs.is_similar(o) is expected, shift
    assert o.is_similar(xrs) is expected, shift
  # u_iso
  for shift, expected in [(0.5*eps, True), (2*eps, False)]:
    o = xrs.deep_copy_scatterers()
    u = o.extract_u_iso_or_u_equiv(); u[3] += shift
    o.set_u_iso(values = u)
    assert xrs.is_similar(o) is expected, shift
  # occupancy
  for shift, expected in [(0.5*eps, True), (2*eps, False)]:
    o = xrs.deep_copy_scatterers()
    q = o.scatterers().extract_occupancies(); q[11] -= shift
    o.scatterers().set_occupancies(q)
    assert xrs.is_similar(o) is expected, shift
  # larger eps
  assert xrs.is_similar(o, eps = 1.e-3)
  # scattering type
  o = xrs.deep_copy_scatterers()
  assert o.scatterers()[0].scattering_type == "C"
  o.scatterers()[0].scattering_type = "N"
  assert not xrs.is_similar(o)
  # number of scatterers
  o = xrs.select(flex.bool([i > 0 for i in range(xrs.scatterers().size())]))
  assert not xrs.is_similar(o) and not o.is_similar(xrs)
  # symmetry
  o = xray.structure(
    crystal_symmetry = crystal.symmetry((30, 40, 55, 90, 90, 90), "P212121"),
    scatterers       = xrs.scatterers())
  assert not xrs.is_similar(o)
  # no scatterers
  e = xray.structure(crystal_symmetry = xrs.crystal_symmetry())
  assert e.is_similar(e.deep_copy_scatterers())

def exercise_time():
  n = 200000
  flex.set_random_seed(0)
  sites = flex.vec3_double(flex.random_double(n), flex.random_double(n),
    flex.random_double(n))
  u = flex.random_double(n)*0.5
  scatterers = flex.xray_scatterer(n, xray.scatterer(scattering_type = "C"))
  xrs = xray.structure(
    crystal_symmetry = crystal.symmetry((300, 400, 500, 90, 90, 90), "P1"),
    scatterers       = scatterers)
  xrs.set_sites_frac(sites)
  xrs.set_u_iso(values = u)
  o = xrs.deep_copy_scatterers()
  t0 = time.time()
  assert xrs.is_similar(o)
  t = time.time()-t0
  print("is_similar, %d scatterers: %6.3f s"%(xrs.scatterers().size(), t))
  assert t < 1.0, t # was 4-5 s (element by element in Python)

if(__name__ == "__main__"):
  t0 = time.time()
  exercise_what_differs()
  exercise_time()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
