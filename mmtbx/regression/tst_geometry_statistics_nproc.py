from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
from mmtbx.secondary_structure.build import ss_idealization as ssb
from cctbx import uctbx
from scitbx.array_family import flex
from libtbx.utils import null_out
from six.moves import cStringIO as StringIO

"""
model.geometry_statistics(nproc=...): validation of geometry (clashscore,
Ramachandran plot, rotamers, C-beta deviations, CaBLAM, omega, Rama-Z) is
done in parallel, one part per process; the result is the same as with one
process.
"""

def get_model():
  ph = ssb.secondary_structure_from_sequence(
    ssb.alpha_helix_str, "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY")
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart = ph.atoms().extract_xyz(), buffer_layer = 5)
  ph.atoms().set_xyz(box.sites_cart)
  model = mmtbx.model.manager(
    model_input      = None,
    pdb_hierarchy    = ph,
    crystal_symmetry = box.crystal_symmetry(),
    log              = null_out())
  model.process(make_restraints = True)
  # distorted, so that statistics are not all zeros
  xrs = model.get_xray_structure()
  flex.set_random_seed(0)
  xrs.shake_sites_in_place(mean_distance = 0.4)
  model.set_sites_cart(sites_cart = xrs.sites_cart())
  return model

def numbers(gs):
  r = gs.result()
  return [
    r.bond.mean, r.angle.mean, r.chirality.mean, r.dihedral.mean,
    r.planarity.mean, r.nonbonded.mean,
    r.ramachandran.outliers, r.ramachandran.allowed, r.ramachandran.favored,
    r.rotamer.outliers, r.rotamer.allowed, r.rotamer.favored,
    r.c_beta.outliers,
    r.clash.score,
    r.molprobity_score,
    r.cablam.outliers, r.cablam.disfavored, r.cablam.ca_outliers,
    r.omega.cis_proline, r.omega.cis_general, r.omega.twisted_general,
    r.omega.twisted_proline,
    r.rama_z.whole.value, r.rama_z.whole.std, r.rama_z.whole.n]

def text(gs):
  s = StringIO()
  gs.show(log = s)
  return s.getvalue()

def exercise():
  model = get_model()
  g1 = model.geometry_statistics()
  assert len(g1.computed_in_parallel()) == 0
  n1, t1 = numbers(g1), text(g1)
  assert n1[6] is not None and n1[9] is not None
  assert n1[0] > 0.01 and n1[1] > 1 # distorted model
  for nproc in [2, 7]:
    g = model.geometry_statistics(nproc = nproc)
    parts = g.computed_in_parallel()
    assert sorted(parts) == sorted(["clash", "ramachandran", "rotamer",
      "c_beta", "cablam", "omega", "rama_z_score"]), parts
    n, t = numbers(g), text(g)
    assert n == n1, [nproc, n, n1]
    assert t == t1
    # objects that come with the numbers are there too
    r = g.result()
    assert r.ramachandran.ramalyze.percent_favored == n1[8]
    assert r.rotamer.rotalyze.percent_outliers == n1[9]
    assert len(r.cablam.gui_table.data) == len(g1.result().cablam.gui_table.data)
  # nproc=1 or None: nothing is done in other processes
  for nproc in [1, None]:
    g = model.geometry_statistics(nproc = nproc)
    assert len(g.computed_in_parallel()) == 0
    assert numbers(g) == n1

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
