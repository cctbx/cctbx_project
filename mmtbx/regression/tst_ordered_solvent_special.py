from __future__ import absolute_import, division, print_function
import time
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from scitbx import matrix
import scitbx.lbfgs
import mmtbx.solvent.ordered_solvent as ordered_solvent
from mmtbx.refinement import geometry_minimization
from mmtbx.regression import tst_ordered_solvent_altlocs as t

"""
Water on and next to special positions in mmtbx.solvent.ordered_solvent.
"""

# 2-fold axes are parallel to b and go through (0,0), (10,0), (0,10) and
# (10,10) in (x,z)
P2 = "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1 2 1"
P222 = "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 2 2 2"

def ser(shift):
  s = matrix.col(shift)
  return [
    t.atom_line(1, " N  ", t.N+s),
    t.atom_line(2, " CA ", t.CA+s),
    t.atom_line(3, " C  ", t.C+s),
    t.atom_line(4, " O  ", t.O+s),
    t.atom_line(5, " CB ", t.CB+s),
    t.atom_line(6, " OG ", t.OGA+s)]

# O of the residue is 2.8 A from the axis at (10, y, 10)
SHIFT = (-4.1, 0, 0)
W = matrix.col((10., 12.4, 10.)) # on the axis, 2.8 A from O and its copy

def run_manager(lines_true, lines, cs, d_min=1.2, params=None):
  model_true = t.get_model(lines_true, cs=cs)
  model = t.get_model(lines, cs=cs)
  for m in [model_true, model]:
    m.setup_scattering_dictionaries(scattering_table="n_gaussian")
  model.process(make_restraints=True)
  fmodel = t.get_fmodel(model_true, model, d_min)
  if(params is None):
    params = ordered_solvent.master_params().extract()
  params.add_hydrogens.d_min = None
  o = t.manager_no_map_filter(fmodel=fmodel, fmodels=None, model=model,
    is_neutron_scat_table=False, params=params, log=null_out())
  assert approx_equal(o.model.get_xray_structure().sites_cart(),
    o.model.get_hierarchy().atoms().extract_xyz())
  t.no_relabel_needed(o.model)
  return o

def n_coordination_restraints(model, i_seq):
  g = model.get_restraints_manager().geometry
  origin_id = ordered_solvent.coordination_origin_id()
  pp = g.pair_proxies(sites_cart = model.get_sites_cart())
  n = 0
  for p in pp.bond_proxies.simple:
    if(p.origin_id == origin_id and i_seq in p.i_seqs): n += 1
  for p in pp.bond_proxies.asu:
    if(p.origin_id == origin_id and i_seq in [p.i_seq, p.j_seq]): n += 1
  return n

def check_on_axis(model, xyz=W, y_tolerance=0.3):
  a = t.find_water(model, xyz, tolerance=0.4)
  assert approx_equal(a.xyz[0], xyz[0], 1.e-3), a.xyz
  assert approx_equal(a.xyz[2], xyz[2], 1.e-3), a.xyz
  assert abs(a.xyz[1]-xyz[1]) < y_tolerance, a.xyz
  xrs = model.get_xray_structure()
  assert a.i_seq in list(xrs.special_position_indices())
  g = model.get_restraints_manager().geometry
  assert g.site_symmetry_table.is_special_position(i_seq = a.i_seq)
  return a

def exercise_on_axis():
  """
  Water on 2-fold axis is built on the axis, has no coordination restraint,
  stays on the axis in coordinate minimization.
  """
  lines = ser(SHIFT)
  lines_true = lines + [t.hoh(1, W)]
  o = run_manager(lines_true, lines, cs=P2)
  m = o.model
  assert len(t.waters(m)) == 1, len(t.waters(m))
  a = check_on_axis(m)
  p = ordered_solvent.master_params().extract()
  assert p.occupancy_min <= a.occ <= 1, a.occ
  assert n_coordination_restraints(m, a.i_seq) == 0
  # coordinate minimization
  g = m.get_restraints_manager().geometry
  sites_cart = m.get_sites_cart()
  sites_cart_start = sites_cart.deep_copy()
  geometry_minimization.lbfgs(
    sites_cart                         = sites_cart,
    geometry_restraints_manager        = g,
    geometry_restraints_flags          = None,
    lbfgs_termination_params           = scitbx.lbfgs.termination_parameters(
      max_iterations = 10),
    correct_special_position_tolerance = 1.0)
  assert sites_cart.max_distance(sites_cart_start) > 0.01 # something moved
  xyz = sites_cart[a.i_seq]
  assert approx_equal(xyz[0], 10, 1.e-6) and approx_equal(xyz[2], 10, 1.e-6)

def exercise_input_water_on_axis():
  """ Water of input model on the axis stays as it is """
  lines = ser(SHIFT) + [t.hoh(1, W, b=15)]
  o = run_manager(lines, lines, cs=P2)
  m = o.model
  a = check_on_axis(m, y_tolerance=1.e-3)
  assert approx_equal(a.occ, 1) and approx_equal(a.b, 15, 1.e-2)
  # 0.1 A off the axis: on the axis in the result
  lines = ser(SHIFT) + [t.hoh(1, W+matrix.col((0.1,0,0)), b=15)]
  o = run_manager(ser(SHIFT) + [t.hoh(1, W, b=15)], lines, cs=P2)
  check_on_axis(o.model, y_tolerance=1.e-3)

def exercise_near_axis():
  """
  Water 0.5 A off the axis is 1 A from its own symmetry copy: removed,
  unless it is to be kept.
  """
  lines = ser(SHIFT) + [t.hoh(1, W+matrix.col((0.5,0,0)))]
  m = t.get_model(lines, cs=P2)
  assert m.get_xray_structure().special_position_indices().size() == 0
  nb = ordered_solvent.get_water_neighbors(m, distance_cutoff=3.2)
  i_w = m.size()-1
  assert approx_equal(dict(nb[i_w])[i_w], 1.0, 1.e-3) # own copy
  r = ordered_solvent.filter_by_distance(m.deep_copy())
  assert len(t.waters(r)) == 0
  fixed = flex.bool(m.size(), False)
  fixed[i_w] = True
  r = ordered_solvent.filter_by_distance(m.deep_copy(),
    fixed_water_selection=fixed)
  assert len(t.waters(r)) == 1
  # 2.8 A from own copy: this is its partner
  lines = [t.hoh(1, matrix.col((1.4, 5., 0.)))]
  m = t.get_model(lines, cs=P2)
  r = ordered_solvent.filter_by_distance(m.deep_copy())
  assert len(t.waters(r)) == 1
  # site next to its own copy is not accepted
  m = t.get_model(ser(SHIFT), cs=P2)
  r = t.assign(m, [W+matrix.col((0.6,0,0)), W+matrix.col((0,3.,0))])
  assert list(r.keep) == [False, True]

def exercise_anisotropic():
  """ Anisotropic water on the axis: ADP has symmetry of the site """
  lines = ser(SHIFT)
  lines_true = lines + [t.hoh(1, W)]
  params = ordered_solvent.master_params().extract()
  params.new_solvent = "anisotropic"
  o = run_manager(lines_true, lines, cs=P2, params=params)
  m = o.model
  a = check_on_axis(m)
  xrs = m.get_xray_structure()
  sc = xrs.scatterers()[a.i_seq]
  assert sc.flags.use_u_aniso()
  assert xrs.site_symmetry_table().get(a.i_seq).is_compatible_u_star(
    sc.u_star)
  assert a.uij != (-1,-1,-1,-1,-1,-1)

def exercise_p222():
  """ Water at the intersection of three 2-fold axes """
  w = matrix.col((10., 10., 10.))
  lines = ser((-3.28, -4.38, 0)) # O is 2.8 A from w, not on an axis
  assert approx_equal((t.O+matrix.col((-3.28, -4.38, 0))-w).length(), 2.8,
    0.01)
  lines_true = lines + [t.hoh(1, w)]
  o = run_manager(lines_true, lines, cs=P222)
  m = o.model
  a = t.find_water(m, w, tolerance=1.e-3)
  xrs = m.get_xray_structure()
  assert a.i_seq in list(xrs.special_position_indices())
  assert xrs.scatterers()[a.i_seq].multiplicity() == 1
  assert n_coordination_restraints(m, a.i_seq) == 0

def run():
  exercise_near_axis()
  exercise_on_axis()
  exercise_input_water_on_axis()
  exercise_anisotropic()
  exercise_p222()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %6.2f"%(time.time()-t0))
