from __future__ import absolute_import, division, print_function
import time
from libtbx.test_utils import approx_equal
from scitbx import matrix
import scitbx.lbfgs
from cctbx import geometry_restraints
import mmtbx.solvent.ordered_solvent as ordered_solvent
from mmtbx.refinement import geometry_minimization
from mmtbx.regression import tst_ordered_solvent_altlocs as t
from mmtbx.regression import tst_ordered_solvent_chains as tc

"""
Water coordination restraints of mmtbx.solvent.ordered_solvent: distance
between water and the atom that coordinates it is restrained to its current
value.
"""

W_O = t.O + matrix.col((0, 2.8, 0)) # 2.8 A from O of SER

def get_model(lines, cs=t.cryst1):
  model = t.get_model(lines, cs=cs)
  model.process(make_restraints=True)
  return model

def coordination(model):
  """ [(i_seq, j_seq, distance_ideal, weight, delta)] """
  g = model.get_restraints_manager().geometry
  origin_id = ordered_solvent.coordination_origin_id()
  sites_cart = model.get_sites_cart()
  pp = g.pair_proxies(sites_cart = sites_cart)
  result = []
  for p in pp.bond_proxies.simple:
    if(p.origin_id != origin_id): continue
    b = geometry_restraints.bond(sites_cart=sites_cart, proxy=p)
    result.append((p.i_seqs[0], p.i_seqs[1], p.distance_ideal, p.weight,
      b.delta))
  asu_mappings = pp.bond_proxies.asu_mappings()
  for p in pp.bond_proxies.asu:
    if(p.origin_id != origin_id): continue
    b = geometry_restraints.bond(sites_cart=sites_cart,
      asu_mappings=asu_mappings, proxy=p)
    result.append((p.i_seq, p.j_seq, p.distance_ideal, p.weight, b.delta))
  return result

def n_bonds(model):
  g = model.get_restraints_manager().geometry
  pp = g.pair_proxies(sites_cart = model.get_sites_cart())
  return pp.bond_proxies.simple.size() + pp.bond_proxies.asu.size()

def nonbonded_pairs(model):
  g = model.get_restraints_manager().geometry
  pp = g.pair_proxies(sites_cart = model.get_sites_cart())
  result = set()
  for p in pp.nonbonded_proxies.simple:
    result.add(tuple(sorted(p.i_seqs)))
  for p in pp.nonbonded_proxies.asu:
    result.add(tuple(sorted((p.i_seq, p.j_seq))))
  return result

def i_seq_of(model, resname, name, altloc=""):
  for a in model.get_hierarchy().atoms():
    if(a.parent().resname == resname and a.name.strip() == name and
       a.parent().altloc == altloc): return a.i_seq

def exercise_basic():
  """
  One restraint per water, to the nearest N/O/S; the pair is not subject to
  nonbonded repulsion; rebuilt for current distances; removed if disabled.
  """
  model = get_model(t.ser_lines(altlocs=False) + [t.hoh(1, W_O)])
  params = ordered_solvent.master_params().extract()
  assert params.coordination_restraints.enabled
  assert approx_equal(params.coordination_restraints.sigma, 0.5)
  i_o, i_w = i_seq_of(model, "SER", "O"), i_seq_of(model, "HOH", "O")
  n_covalent = n_bonds(model)
  assert (i_o, i_w) in nonbonded_pairs(model)
  i_c = i_seq_of(model, "SER", "C")
  assert (i_c, i_w) in nonbonded_pairs(model)
  stats = model.geometry_statistics(use_hydrogens=False)
  bond_before = (stats.bond().mean, stats.bond().n)
  n = ordered_solvent.update_coordination_restraints(model=model,
    params=params)
  assert n == 1
  c = coordination(model)
  assert len(c) == 1
  assert sorted(c[0][:2]) == [i_o, i_w]
  assert approx_equal(c[0][2], 2.8, 1.e-3)       # ideal: current distance
  assert approx_equal(c[0][3], 1/0.5**2)         # sigma 0.5
  assert approx_equal(c[0][4], 0, 1.e-6)
  assert n_bonds(model) == n_covalent+1
  assert (i_o, i_w) not in nonbonded_pairs(model)
  # only this pair: atom bonded to O still repels water
  assert (i_c, i_w) in nonbonded_pairs(model)
  # not in bond statistics
  stats = model.geometry_statistics(use_hydrogens=False)
  assert approx_equal((stats.bond().mean, stats.bond().n), bond_before)
  # in .geo file
  geo = model.restraints_as_geo()
  assert "Water coordination" in geo
  # water moved: restraint follows after update, still one
  sites_cart = model.get_sites_cart()
  sites_cart[i_w] = tuple(matrix.col(sites_cart[i_w])+matrix.col((0,0.3,0)))
  model.set_sites_cart(sites_cart)
  assert approx_equal(abs(coordination(model)[0][4]), 0.3, 1.e-3)
  assert ordered_solvent.update_coordination_restraints(model=model,
    params=params) == 1
  c = coordination(model)
  assert len(c) == 1 and approx_equal(c[0][2], 3.1, 1.e-3)
  assert n_bonds(model) == n_covalent+1
  # survives selection of atoms
  m2 = model.select(~model.selection("name CB"))
  assert len(coordination(m2)) == 1
  # removed
  assert ordered_solvent.remove_coordination_restraints(model=model) == 1
  assert coordination(model) == [] and n_bonds(model) == n_covalent
  assert (i_o, i_w) in nonbonded_pairs(model)
  # disabled: none is built and existing ones are removed
  ordered_solvent.update_coordination_restraints(model=model, params=params)
  params.coordination_restraints.enabled = False
  assert ordered_solvent.update_coordination_restraints(model=model,
    params=params) == 0
  assert coordination(model) == [] and n_bonds(model) == n_covalent
  # model without restraints: nothing to do
  m3 = t.get_model(t.ser_lines(altlocs=False) + [t.hoh(1, W_O)])
  params.coordination_restraints.enabled = True
  assert ordered_solvent.update_coordination_restraints(model=m3,
    params=params) == 0

def exercise_partner():
  """ Which atom water is restrained to """
  params = ordered_solvent.master_params().extract()
  # conformers: water A next to OG B (2.6 A) and N (2.9 A): N
  w = t.N + matrix.col((0, -2.9, 0))
  lines = t.ser_lines() + [t.hoh(1, w, altloc="A", occ=0.7)]
  model = get_model(lines)
  assert ordered_solvent.update_coordination_restraints(model=model,
    params=params) == 1
  c = coordination(model)[0]
  assert sorted(c[:2]) == [i_seq_of(model, "SER", "N"),
                           i_seq_of(model, "HOH", "O", "A")], c
  # water of the same conformer as the atom: that atom if it is the nearest
  w = t.OGA + (t.OGA-t.CB).normalize()*2.7
  model = get_model(t.ser_lines() + [t.hoh(1, w, altloc="A", occ=0.7)])
  assert ordered_solvent.update_coordination_restraints(model=model,
    params=params) == 1
  c = coordination(model)[0]
  assert sorted(c[:2]) == [i_seq_of(model, "SER", "OG", "A"),
                           i_seq_of(model, "HOH", "O", "A")], c
  # protein atom is preferred to a closer water; water that has only water
  # is restrained to the nearest one
  w2 = W_O + matrix.col((0, 2.6, 0))
  w3 = w2 + matrix.col((0, 2.7, 0))
  lines = t.ser_lines(altlocs=False) + [
    t.hoh(1, W_O), t.hoh(2, w2), t.hoh(3, w3)]
  model = get_model(lines)
  assert ordered_solvent.update_coordination_restraints(model=model,
    params=params) == 3
  pairs = sorted([tuple(sorted(c[:2])) for c in coordination(model)])
  i_o = i_seq_of(model, "SER", "O")
  iw = [a.i_seq for a in t.waters(model)]
  assert pairs == sorted([(i_o, iw[0]), (iw[0], iw[1]), (iw[1], iw[2])]), pairs
  # two waters that have only each other: one restraint
  lines = t.ser_lines(altlocs=False) + [
    t.hoh(1, t.FREE), t.hoh(2, t.FREE+matrix.col((2.8, 0, 0)))]
  model = get_model(lines)
  assert ordered_solvent.update_coordination_restraints(model=model,
    params=params) == 1
  # two waters coordinated by the same atom still repel each other
  w_a = t.O + matrix.col((0, 2.8, 0))
  w_b = t.O + matrix.col((-2.0, 1.9, 0.5))
  lines = t.ser_lines(altlocs=False) + [t.hoh(1, w_a), t.hoh(2, w_b)]
  model = get_model(lines)
  iw = [a.i_seq for a in t.waters(model)]
  assert tuple(iw) in nonbonded_pairs(model)
  assert ordered_solvent.update_coordination_restraints(model=model,
    params=params) == 2
  i_o = i_seq_of(model, "SER", "O")
  assert sorted([tuple(sorted(c[:2])) for c in coordination(model)]) == \
    sorted([(i_o, iw[0]), (i_o, iw[1])])
  assert tuple(iw) in nonbonded_pairs(model)
  # nothing within dist_max: no restraint
  model = get_model(t.ser_lines(altlocs=False) + [
    t.hoh(1, t.O + matrix.col((0, 4.0, 0)))])
  assert ordered_solvent.update_coordination_restraints(model=model,
    params=params) == 0

def exercise_symmetry():
  """ Coordinating atom is a symmetry copy; restraint can be removed """
  params = ordered_solvent.master_params().extract()
  w = tc.symmetry_copy(tc.A_O + matrix.col((0, 2.8, 0)))
  model = get_model(tc.gly("A", (0,0,0)) + [tc.hoh(1, w)], cs=tc.CS)
  n_covalent = n_bonds(model)
  assert ordered_solvent.update_coordination_restraints(model=model,
    params=params) == 1
  c = coordination(model)
  assert len(c) == 1
  assert approx_equal(c[0][2], 2.8, 1.e-3) and approx_equal(c[0][4], 0, 1.e-6)
  # again and again
  for i in range(2):
    assert ordered_solvent.update_coordination_restraints(model=model,
      params=params) == 1
    assert len(coordination(model)) == 1 and n_bonds(model) == n_covalent+1
  assert ordered_solvent.remove_coordination_restraints(model=model) == 1
  assert n_bonds(model) == n_covalent

def exercise_two_copies():
  """
  Atom and its symmetry copy are both next to the water: no restraint to
  that atom (it is not known to which copy it would be made); another atom
  is used.
  """
  params = ordered_solvent.master_params().extract()
  # P 1 2 1: O of chain A is 2.0 A from the axis, water is 0.6 A off the
  # axis on the other side: 2.6 A from O and 3.0 A from its copy... so both
  # copies of O are in reach
  cs = "CRYST1   25.000   25.000   25.000  90.00  90.00  90.00 P 1 2 1"
  o = matrix.col((6.3, 7.4, 5.0))
  shift = matrix.col((12.5-2.0, 0, 12.5)) - matrix.col((o[0], 0, o[2]))
  w = matrix.col((12.5+0.7, 7.4, 12.5))
  model = get_model(tc.gly("A", shift) + [tc.hoh(1, w)], cs=cs)
  nb = dict(ordered_solvent.get_water_neighbors(model, 3.3))
  i_w = model.size()-1
  n = ordered_solvent.update_coordination_restraints(model=model,
    params=params)
  i_o = 3
  for c in coordination(model):
    assert i_o not in c[:2], c
    assert approx_equal(c[4], 0, 1.e-6)

def minimize(model):
  g = model.get_restraints_manager().geometry
  sites_cart = model.get_sites_cart()
  geometry_minimization.lbfgs(
    sites_cart                         = sites_cart,
    geometry_restraints_manager        = g,
    geometry_restraints_flags          = None,
    lbfgs_termination_params           = scitbx.lbfgs.termination_parameters(
      max_iterations = 50),
    correct_special_position_tolerance = 1.0)
  return sites_cart

def exercise_minimization():
  """
  Water 2.2 A from O is pushed away by nonbonded repulsion; not so with the
  restraint.
  """
  params = ordered_solvent.master_params().extract()
  w = t.O + matrix.col((0, 2.2, 0))
  shifts = []
  for restrain in [False, True]:
    model = get_model(t.ser_lines(altlocs=False) + [t.hoh(1, w)])
    i_o, i_w = i_seq_of(model, "SER", "O"), i_seq_of(model, "HOH", "O")
    if(restrain):
      ordered_solvent.update_coordination_restraints(model=model,
        params=params)
    sites_cart = minimize(model)
    d = (matrix.col(sites_cart[i_o])-matrix.col(sites_cart[i_w])).length()
    shifts.append(d-2.2)
  assert shifts[0] > 0.2, shifts
  assert abs(shifts[1]) < 0.05, shifts

def exercise_manager():
  """ Water update leaves restraints for water present at its end """
  w_o = W_O
  w_n = t.N + matrix.col((0, 0, -2.8))
  lines = t.ser_lines(altlocs=False) + [t.hoh(1, w_o)]
  lines_true = lines + [t.hoh(2, w_n)]
  o = t.run_manager(lines_true, lines, d_min=1.2, map_filter=False)
  m = o.model
  c = coordination(m)
  i_ws = [a.i_seq for a in t.waters(m)]
  assert len(i_ws) >= 2
  assert len(c) >= 2
  for i, j, d_ideal, weight, delta in c:
    assert i in i_ws or j in i_ws
    assert approx_equal(delta, 0, 1.e-6)
    assert approx_equal(weight, 4)
  # second update starts with them present
  fmodel = o.fmodel
  o2 = t.manager_no_map_filter(fmodel=fmodel, fmodels=None, model=m,
    is_neutron_scat_table=False, params=o.params, log=o.log)
  assert len(coordination(o2.model)) == len(c)

def run():
  exercise_basic()
  exercise_partner()
  exercise_symmetry()
  exercise_two_copies()
  exercise_minimization()
  exercise_manager()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %6.2f"%(time.time()-t0))
