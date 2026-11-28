from __future__ import absolute_import, division, print_function
import time
import iotbx.phil
import mmtbx.model
import mmtbx.f_model
import mmtbx.refinement.real_space.adp as adp
import mmtbx.refinement.real_space.occupancy as occupancy
from mmtbx.secondary_structure.build import ss_idealization as ssb
from cctbx import uctbx
from iotbx import map_manager, map_model_manager
from scitbx.array_family import flex
from libtbx import group_args
from libtbx.utils import null_out

"""
Refinement of B-factors and occupancies against data calculated from a known
model: occupancies only, B-factors only, and both (occupancies, then
B-factors, twice: what phenix.real_space_refine does in its last
macro-cycle).

Two kinds of data:
- "exact": amplitudes calculated from the answer model, given to the
  refinement in place of those it makes from the map. Nothing is lost on the
  way, so the answer must come back (nearly) exactly; this checks targets,
  gradients, minimizers and the logic around them.
- "map": the map calculated from the answer model, used the way the program
  uses it (box around the model, mask around atoms). The answer comes back
  within what this route allows.
"""

d_min = 2.0
q_answer, q_start = 0.6, 0.9
b_start = 60.

occupancies_phil_str = """
individual = None
  .type = str
  .multiple = True
remove_selection = None
  .type = str
  .multiple = True
constrained_group
  .multiple = True
{
  selection = None
    .type = str
    .multiple = True
}
"""

def get_inputs():
  """
  Helix; B-factors vary along the chain and between main and side chains;
  side-chain atoms of one residue (ARG 15) have occupancy 0.6. Returns the
  model (with the answer in it), the answer and the map calculated from it.
  """
  ph = ssb.secondary_structure_from_sequence(
    ssb.alpha_helix_str, "ACDEFGHIKLMNPQRSTVWY")
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart = ph.atoms().extract_xyz(), buffer_layer = 7)
  ph.atoms().set_xyz(box.sites_cart)
  b, partial = flex.double(), flex.bool()
  for i, rg in enumerate(ph.residue_groups()):
    for a in rg.atoms():
      main = a.name.strip() in ["N", "CA", "C", "O"]
      b.append(40 + 3*abs(i-10) + (0 if main else 10))
      partial.append(i == 14 and not main and a.name.strip() != "CB")
  assert partial.count(True) > 4
  q = flex.double(b.size(), 1).set_selected(partial, q_answer)
  ph.atoms().set_b(b)
  ph.atoms().set_occ(q)
  model = mmtbx.model.manager(
    model_input      = None,
    pdb_hierarchy    = ph,
    crystal_symmetry = box.crystal_symmetry(),
    log              = null_out())
  model.setup_scattering_dictionaries(scattering_table = "electron")
  fc = model.get_xray_structure().structure_factors(d_min = d_min).f_calc()
  map_data = fc.fft_map(resolution_factor = 0.25).real_map_unpadded()
  mm = map_manager.map_manager(
    map_data                   = map_data,
    unit_cell_grid             = map_data.all(),
    unit_cell_crystal_symmetry = model.crystal_symmetry(),
    wrapping                   = False)
  return group_args(model = model, mm = mm, b = b, q = q, partial = partial,
    xrs_answer = model.get_xray_structure().deep_copy_scatterers())

class exact_data(object):
  """
  While active, B-factor and occupancy refinement get amplitudes calculated
  from the answer model instead of those made from the map.
  """
  def __init__(self, xrs_answer):
    self.xrs_answer = xrs_answer
  def __enter__(self):
    self.saved = adp.map_and_model_to_fmodel
    xrs_answer = self.xrs_answer
    def f(map_data, xray_structure, atom_radius, d_min, reset_adp=True):
      f_obs = abs(xrs_answer.structure_factors(d_min = d_min).f_calc())
      fmodel = mmtbx.f_model.manager(
        f_obs          = f_obs,
        xray_structure = xray_structure.deep_copy_scatterers())
      fmodel.update(target_name = "ls_wunit_k1")
      return fmodel
    adp.map_and_model_to_fmodel = f
    assert occupancy.dependency is adp
  def __exit__(self, *args):
    adp.map_and_model_to_fmodel = self.saved

def refine(inp, what, exact, b_restraints, rounds_both = 2):
  """
  what: "occupancy", "adp" or "both". Start: wrong values of what is
  refined, the answer for the rest.
  """
  model = inp.model.deep_copy()
  if(what in ["adp", "both"]):
    model.set_b_iso(values = flex.double(inp.b.size(), b_start))
  if(what in ["occupancy", "both"]):
    model.set_occupancies(
      values = flex.double(inp.q.size(), 1).set_selected(inp.partial, q_start))
  mmm = map_model_manager.map_model_manager(
    model = model, map_manager = inp.mm.deep_copy())
  occ_params = iotbx.phil.parse(occupancies_phil_str).extract()
  def run():
    rounds = rounds_both if what == "both" else 1
    for i_round in range(rounds):
      if(what in ["occupancy", "both"]):
        occupancy.ncs_aware_refinement(
          map_model_manager = mmm,
          d_min             = d_min,
          atom_radius       = 3,
          occupancy_params  = occ_params,
          log               = null_out())
      if(what in ["adp", "both"]):
        adp.ncs_aware_refinement(
          map_model_manager = mmm,
          d_min             = d_min,
          atom_radius       = 3,
          individual        = True,
          restraints_weight = 1 if b_restraints else None,
          group_mode        = "one_per_residue" if b_restraints else None,
          log               = null_out())
  if(exact):
    with exact_data(xrs_answer = inp.xrs_answer): run()
  else: run()
  b = mmm.model().get_b_iso()
  q = mmm.model().get_occ()
  qp = q.select(inp.partial)
  return group_args(
    b_diff     = flex.mean(flex.abs(b - inp.b)),
    b_max      = flex.max(flex.abs(b - inp.b)),
    q_diff     = flex.max(flex.abs(qp - q_answer)),
    q_mean     = flex.mean(qp),
    q_others   = flex.max(flex.abs(q.select(~inp.partial) - 1)))

# (what, exact data, restraints on B): limits on mean |B - answer| and
# largest |occupancy - answer|
cases = [
  ("occupancy", True,  False, None, 0.01),
  ("adp",       True,  False, 0.05, None),
  ("adp",       True,  True,  4.0,  None), # restraints: not the exact answer
  ("both",      True,  False, 0.5,  0.06),
  ("occupancy", False, False, None, 0.02),
  ("adp",       False, False, 5.0,  None), # no restraints: noisy
  ("adp",       False, True,  1.5,  None),
  ("both",      False, True,  1.5,  0.10),
]

def exercise(report_only = False):
  inp = get_inputs()
  for what, exact, b_restraints, b_limit, q_limit in cases:
    r = refine(inp = inp, what = what, exact = exact,
      b_restraints = b_restraints)
    print("%-9s %-5s data, B %-12s: mean |B-answer| %5.2f (max %5.2f) "
      "occupancy %4.2f +- %4.2f (answer %3.1f)"%(what,
      "exact" if exact else "map",
      "restrained" if b_restraints else "unrestrained",
      r.b_diff, r.b_max, r.q_mean, r.q_diff, q_answer))
    if(report_only): continue
    assert r.q_others < 1.e-6
    if(b_limit is None): assert r.b_diff < 1.e-6, r.b_diff # not refined
    else:                assert r.b_diff < b_limit, r.b_diff
    if(q_limit is None): assert abs(r.q_mean - q_answer) < 1.e-6
    else:                assert r.q_diff < q_limit, r.q_diff

def exercise_rounds():
  """
  Both refined, exact data: B-factors and occupancies of the same atoms are
  correlated, so the pair (occupancies, B-factors) gets to the answer over
  several rounds; every round is closer than the previous one.
  """
  inp = get_inputs()
  q_diff, b_max = [], []
  for rounds in [1, 2, 5]:
    r = refine(inp = inp, what = "both", exact = True, b_restraints = False,
      rounds_both = rounds)
    print("both, exact data, rounds %d: largest |occupancy-answer| %5.3f "
      "largest |B-answer| %5.2f"%(rounds, r.q_diff, r.b_max))
    q_diff.append(r.q_diff)
    b_max.append(r.b_max)
  assert q_diff[0] > q_diff[1] > q_diff[2], q_diff
  assert b_max[0] > b_max[1] > b_max[2], b_max
  assert q_diff[2] < 0.025 and b_max[2] < 2, [q_diff, b_max]

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  exercise_rounds()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
