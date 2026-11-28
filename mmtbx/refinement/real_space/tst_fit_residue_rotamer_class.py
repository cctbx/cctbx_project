from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space
from mmtbx.refinement.real_space import fit_residue
from mmtbx.rotamer.rotamer_eval import RotamerEval
from mmtbx.secondary_structure.build import ss_idealization as ssb
from scitbx.array_family import flex
from libtbx import group_args

"""
Side-chain fitting keeps a trial conformation (rotamer fitting, then
tune-up) instead of the start one by map values; it does not do so if the
rotamer class (favoured, allowed, outlier) gets worse, unless the gain in map
value is large. It used to accept any gain in map value as long as the
result was not an outlier.
"""

def get_residue(resname = "LEU"):
  ph = ssb.secondary_structure_from_sequence(ssb.alpha_helix_str, "ALA")
  for r in ph.only_model().only_chain().residues():
    if(r.resname == resname): return r

def exercise_rotamer_rank():
  re = RotamerEval()
  rank = mmtbx.refinement.real_space.rotamer_rank
  r = get_residue()
  assert re.evaluate_residue_2(residue = r) == "Favored"
  assert rank(re, r) == 2
  assert rank(None, r) is None
  assert rank(re, get_residue("ALA")) is None
  # all classes are met when chi1 goes round
  cl = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
    residue         = r,
    mon_lib_srv     = re.mon_lib_srv,
    backbone_sample = False).clusters[0]
  from scitbx.matrix import rotate_point_around_axis
  sites = r.atoms().extract_xyz()
  found = set()
  for angle in range(0, 360, 5):
    s = sites.deep_copy()
    for i in cl.atoms_to_rotate:
      s[i] = rotate_point_around_axis(
        axis_point_1 = sites[cl.axis[0]], axis_point_2 = sites[cl.axis[1]],
        point = sites[i], angle = angle, deg = True)
    r.atoms().set_xyz(s)
    found.add(rank(re, r))
  r.atoms().set_xyz(sites)
  assert found == set([0,1,2]), found

def state(target, rank, shift):
  r = get_residue()
  return group_args(
    vals = {}, sites_cart = r.atoms().extract_xyz() + (shift,0,0),
    target = target, target_neg = 0, exceed_map_max_value = False,
    rot = "OUTLIER" if rank == 0 else "x", rank = rank)

def kept(start, fitting, tuneup = None):
  """
  Which of the states monitor.finalize leaves in the residue.
  """
  r = get_residue()
  m = fit_residue.monitor(id_str = "", selection = flex.size_t([0]),
    map_data = None, unit_cell = None, weights = None, pairs = None,
    cmv = None, rotamer_evaluator = None, log = None)
  m.add = lambda residue, state: None # no map in this test
  m.states["start"]   = state(*start,   shift = 0)
  m.states["fitting"] = state(*fitting, shift = 1)
  if(tuneup is not None): m.states["tuneup"] = state(*tuneup, shift = 2)
  x0 = r.atoms().extract_xyz()[0][0]
  m.finalize(residue = r)
  return ["start","fitting","tuneup"][
    int(round(r.atoms().extract_xyz()[0][0] - x0))]

def exercise_finalize():
  F, A, O = 2, 1, 0
  # better map value, same or better class: taken
  assert kept((10,F), (11,F)) == "fitting"
  assert kept((10,A), (11,F)) == "fitting"
  assert kept((10,O), (11,A)) == "fitting"
  assert kept((10,O), ( 9,F)) == "fitting" # outlier fixed
  # better map value, worse class: not taken
  assert kept((10,F), (11,A)) == "start"
  assert kept((10,F), (11,A), (12,A)) == "start"
  assert kept((10,A), (11,O)) == "start"
  # ... unless the gain in map value is large (the side chain was not in
  # density at all): 1.5 times or more
  assert fit_residue.rotamer_class_min_gain == 1.5
  assert kept((10,F), (14.9,A)) == "start"
  assert kept((10,F), (15.1,A)) == "fitting"
  assert kept((10,F), (30,A), (31,A)) == "tuneup"
  assert kept((10,F), (11,F), (16.4,A)) == "fitting"
  assert kept((10,F), (11,F), (16.6,A)) == "tuneup"
  assert kept((-1,F), (5,A)) == "fitting"
  assert kept((10,A), (20,O)) == "start" # never to an outlier
  # worse map value: not taken
  assert kept((10,F), ( 9,F)) == "start"
  # tune-up after fitting: same rules
  assert kept((10,F), (11,F), (12,F)) == "tuneup"
  assert kept((10,F), (11,F), (12,A)) == "fitting"
  assert kept((10,F), ( 9,F), (12,A)) == "start"
  assert kept((10,F), ( 9,F), (12,F)) == "tuneup"
  # class unknown: map values decide, as before
  assert kept((10,None), (11,A)) == "fitting"

if(__name__ == "__main__"):
  t0 = time.time()
  exercise_rotamer_rank()
  exercise_finalize()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
