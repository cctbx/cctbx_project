from __future__ import absolute_import, division, print_function
import time
import mmtbx.idealized_aa_residues.rotamer_manager as rm
import mmtbx.refinement.real_space
import scitbx.math
from scitbx.array_family import flex
from mmtbx.secondary_structure.build import ss_idealization as ssb
from cctbx import uctbx
from libtbx.test_utils import approx_equal
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_rotamer_fit_ext")

"""
Search for the best rotamer (mmtbx_rotamer_fit_ext.fit).

1. The table of chi angles is taken as a list of lists or as one array of
   numbers: same result. Tables have up to 5.e5 sets of angles; converting
   the list of lists for every residue took longer than the search.
2. rotamer_manager.get_chi_angles_flat gives the table sorted by angles. The
   search does not repeat rotations for the leading angles that a set shares
   with the one before it, and it skips all sets that share leading angles
   that put an atom into a clash or out of the allowed range of map values.
   The result does not depend on the order of the sets (the same table in
   random order, where nearly nothing is shared, gives the same answer), the
   time does.
"""

def exercise():
  m = rm.load(rotamers = "favored_allowed")
  sct = scitbx.math.sin_cos_table(n=10000)
  ph = ssb.secondary_structure_from_sequence(ssb.alpha_helix_str, "AEKRLA")
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart = ph.atoms().extract_xyz(), buffer_layer = 5)
  ph.atoms().set_xyz(box.sites_cart)
  ph.atoms().reset_i_seq()
  xrs = ph.extract_xray_structure(crystal_symmetry = box.crystal_symmetry())
  md = xrs.structure_factors(d_min=3).f_calc().fft_map(
    resolution_factor=0.25).real_map_unpadded()
  md = md / flex.max(md)
  sites_all = ph.atoms().extract_xyz()
  t_list, t_flat, t_sorted, t_random = 0, 0, 0, 0
  n_real_differs = 0
  for r in ph.only_model().only_chain().residues():
    if(r.resname == "ALA"): continue
    co = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
      residue         = r,
      mon_lib_srv     = m.rotamer_evaluator.mon_lib_srv,
      backbone_sample = True)
    chi = m.get_chi_angles(resname = r.resname)
    n_chi = len(chi[0])
    flat = flex.double(chi).as_1d()
    flat_sorted = m.get_chi_angles_flat(resname = r.resname)
    assert flat_sorted is m.get_chi_angles_flat(resname = r.resname) # once
    assert flat_sorted.size() == flat.size() == len(chi)*n_chi
    assert approx_equal(flat_sorted[:n_chi], min(chi))
    assert approx_equal(flat_sorted[-n_chi:], max(chi))
    perm = flex.random_permutation(len(chi))
    random_order = [chi[i] for i in perm]
    assert sorted(random_order) == sorted(chi)
    flat_random = flex.double(random_order).as_1d()
    cl = co.clusters[1:]
    axes = [flex.size_t(c.axis) for c in cl[:n_chi]]
    atr  = [flex.size_t(c.atoms_to_rotate) for c in cl[:n_chi]]
    sites = r.atoms().extract_xyz()
    start = sites.deep_copy()
    common = dict(axes = axes, rotatable_points_indices = atr,
      sin_table = sct.sin_table, cos_table = sct.cos_table, step = sct.step,
      n = sct.n)
    def moving(ref_map_min, ref_map_max):
      return ext.moving(
        sites_cart       = start,
        sites_cart_start = start,
        radii            = flex.double(sites.size(), 1.5),
        weights          = flex.double(sites.size(), 6),
        bonded_pairs     = [],
        ref_map_max      = flex.double(sites.size(), ref_map_max),
        ref_map_min      = flex.double(sites.size(), ref_map_min))
    # nothing to bump into, any map value
    fixed_far = ext.fixed(
      sites_cart = flex.vec3_double([(0,0,0), sites[0]]) + (50,50,50),
      radii      = flex.double([1.5, 1.5]))
    # all other atoms of the model, map values limited
    others = flex.bool(sites_all.size(), True).set_selected(
      r.atoms().extract_i_seq(), False)
    fixed_real = ext.fixed(
      sites_cart = sites_all.select(others),
      radii      = flex.double(others.count(True), 1.5))
    def search(angles):
      t0 = time.time()
      o1 = ext.fit(target_value = -1.e9, angles_array = angles,
        density_map = md, all_points = start, unit_cell = xrs.unit_cell(),
        selection = co.rsr_eval_selection, **common)
      o2 = ext.fit(fixed = fixed_far, angles_array = angles, density_map = md,
        moving = moving(-1.e9, 1.e9), unit_cell = xrs.unit_cell(),
        selection_clash = co.clash_eval_selection,
        selection_rsr = co.rsr_eval_selection, **common)
      o3 = ext.fit(sites_cart_start = sites, angles_array = angles,
        all_points = start, **common)
      t1 = time.time()
      o4 = ext.fit(fixed = fixed_real, angles_array = angles,
        density_map = md, moving = moving(0.02, 0.9),
        unit_cell = xrs.unit_cell(),
        selection_clash = co.clash_eval_selection,
        selection_rsr = co.rsr_eval_selection, **common)
      t2 = time.time()
      return group(
        results = [o1.result(), o1.score(), o2.result(), o2.score(),
                   o3.result(), o4.result(), o4.score()],
        t_all   = t2-t0,
        t_real  = t2-t1)
    def same(a, b):
      for x, y in zip(a.results, b.results):
        if(isinstance(x, float)): assert x == y, [x, y]
        else:
          assert x.size() == y.size()
          assert x.size() in [0, sites.size()]
          if(x.size() > 0): assert flex.max(flex.sqrt((x-y).dot())) == 0
    a = search(chi)
    b = search(flat)
    c = search(flat_sorted)
    d = search(flat_random)
    same(a, b)
    same(a, c)
    same(a, d)
    t_list   += a.t_all
    t_flat   += b.t_all
    t_sorted += c.t_real
    t_random += d.t_real
    # the search found something close to the answer
    dist = flex.max(flex.sqrt((a.results[0] - sites).dot()))
    # limits on map values and neighbours make a difference
    real = a.results[5]
    if(real.size() == 0 or
       flex.max(flex.sqrt((real - a.results[2]).dot())) > 0.1):
      n_real_differs += 1
    print("%s: sets of angles %6d, largest distance to answer %5.2f"%(
      r.resname, len(chi), dist), "restricted search: found",
      real.size() > 0, "time sorted %5.3f random order %5.3f"%(
      c.t_real, d.t_real))
    assert dist < 1.5, dist # carboxylate O are equivalent
  assert n_real_differs > 0
  print("time: list of lists %5.2f one array %5.2f"%(t_list, t_flat))
  assert t_flat < t_list*0.8, [t_list, t_flat]
  print("time, restricted search: sorted %5.2f random order %5.2f"%(
    t_sorted, t_random))
  assert t_sorted < t_random*0.5, [t_sorted, t_random]
  # empty table: nothing found
  o = ext.fit(sites_cart_start = sites, angles_array = flex.double(),
    all_points = start, **common)
  assert o.result().size() == 0

def group(**kw):
  from libtbx import group_args
  return group_args(**kw)

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
