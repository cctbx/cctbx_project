from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from mmtbx.refinement.real_space import individual_sites
from libtbx.test_utils import approx_equal
from cctbx import maptbx
from libtbx import adopt_init_args
import scitbx.rigid_body
import sys
from six.moves import range

def apply_rigid_body_shift(sites_cart, cm, x,y,z, the, psi, phi):
  rot_matrix = scitbx.rigid_body.rb_mat_zyz(
    the=the, psi=psi, phi=phi).rot_mat().as_mat3()
  transl = (x,y,z)
  return rot_matrix * (sites_cart-cm) + transl + cm

class scorer(object):
  def __init__(self, unit_cell, sites_frac, target_map, log=None):
    adopt_init_args(self, locals())
    self.target = maptbx.map_sum_at_sites_frac(
      map_data   = target_map,
      sites_frac = sites_frac)
    self.sites_frac_best = None
    if(log):
      print("Target (start):", self.target, file=log)

  def update(self, sites_cart):
    sites_frac = self.unit_cell.fractionalize(sites_cart)
    t = maptbx.map_sum_at_sites_frac(
      map_data   = self.target_map,
      sites_frac = sites_frac)
    if(t>self.target):
      self.target = t
      self.sites_frac_best = sites_frac.deep_copy()
      if(self.log):
        print("Target:", self.target, file=self.log)

def prepare_maps(fofc, two_fofc, fem, fofc_cutoff=2, two_fofc_cutoff=0.5,
                 fem_cutoff=0.5, connectivity_cutoff=0.5, local_average=True):
  """
  - This takes 3 maps: mFo-DFc, 2mFo-DFc and FEM and combines them into one map
    that is most suitable for real-space refinement.
  - Maps are the boxes extracted around region of interest from the whole unit
    cell map.
  - All maps are expected to be normalized by standard deviation (sigma-scaled)
    BEFORE extracting the box. There is no way to assert it at this point.
  - Map gridding equivalence is asserted.
  """
  m1,m2,m3 = fofc, two_fofc, fem
  # assert identical gridding
  for m_ in [m1,m2,m3]:
    for m__ in [m1,m2,m3]:
      assert m_.all()    == m__.all()
      assert m_.focus()  == m__.focus()
      assert m_.origin() == m__.origin()
  # binarize residual map
  sel = m1 <= fofc_cutoff
  mask = m1  .set_selected( sel, 0)
  mask = mask.set_selected(~sel, 1)
  del sel, m1
  assert approx_equal([flex.max(mask), flex.min(mask)], [1,0])
  def truncate_and_filter(m, cutoff, mask):
    return m.set_selected(m<=cutoff, 0)*mask
  # truncate and filter 2mFo-DFc map
  m2 = truncate_and_filter(m2, two_fofc_cutoff, mask)
  # truncate and filter FEM
  m3 = truncate_and_filter(m3, fem_cutoff, mask)
  del mask
  # combined maps
  def scale(m):
    sd = m.sample_standard_deviation()
    if(sd != 0): return m/sd
    else: return m
  m2 = scale(m2)
  m3 = scale(m3)
  m = (m2+m3)/2.
  del m2, m3
  m = scale(m)
  # connectivity analysis
  co = maptbx.connectivity(map_data=m, threshold=connectivity_cutoff)
  v_max=-1.e+9
  i_max=None
  for i, v in enumerate(co.regions()):
    if(i>0):
      if(v>v_max):
        v_max=v
        i_max=i
  mask2 = co.result()
  selection = mask2==i_max
  mask2 = mask2.set_selected(selection, 1)
  mask2 = mask2.set_selected(~selection, 0)
  assert mask2.count(1) == v_max
  # final filter
  m = m * mask2.as_double()
  if(local_average):
    maptbx.map_box_average(map_data=m, cutoff=0.5, index_span=1)
  return m

def shift_to_center_of_mass(xray_structure, target_map, cutoff=0.5):
  cm = xray_structure.center_of_mass()
  sites_cart = xray_structure.sites_cart()
  cmm = maptbx.center_of_mass(
    map_data=target_map, unit_cell=xray_structure.unit_cell(), cutoff=cutoff)
  x,y,z = cm[0]-cmm[0],cm[1]-cmm[1],cm[2]-cmm[2]
  sites_cart_new = apply_rigid_body_shift(sites_cart=sites_cart,
    cm=cm, x=-x,y=-y,z=-z, the=0, psi=0, phi=0)
  return xray_structure.replace_sites_cart(new_sites=sites_cart_new)

def run_refine(rsr_simple_refiner, xray_structure, scorer, log, weight):
  weight_best = None
  try:
    refined = individual_sites.refinery(
      refiner                  = rsr_simple_refiner,
      xray_structure           = xray_structure,
      start_trial_weight_value = 1.0,
      rms_bonds_limit          = 0.02,
      rms_angles_limit         = 2.0)
    scorer.update(sites_cart=refined.sites_cart_result)
    weight_best = refined.weight_final
  # except KeyboardInterrupt: raise
  except Exception:
    if(log): print("A trial failed: keep going...", file=log)
  return weight_best

def macro_cycle(
      xray_structure,
      target_map,
      geometry_restraints,
      max_iterations = 50,
      expload        = False,
      n_expload      = 1,
      log            = None):
  if(not log):
    if(log is None): log = sys.stdout
  d_min = maptbx.d_min_from_map(
      map_data  = target_map,
      unit_cell = xray_structure.unit_cell())
  all_selection = flex.bool(xray_structure.scatterers().size(),True)
  rsr_simple_refiner = individual_sites.simple(
    target_map                  = target_map,
    selection                   = all_selection,
    real_space_gradients_delta  = d_min/4,
    max_iterations              = max_iterations,
    geometry_restraints_manager = geometry_restraints)
  xray_structure = shift_to_center_of_mass(xray_structure=xray_structure,
    target_map=target_map)
  cm = xray_structure.center_of_mass()
  sites_cart = xray_structure.sites_cart()
  sc = scorer(
    unit_cell  = xray_structure.unit_cell(),
    sites_frac = xray_structure.sites_frac(),
    target_map = target_map,
    log        = log)
  weights = flex.double()
  sampling_range = range(0, 370, 50)   #[0, 90, 180, 270]
  for the in sampling_range:
    for psi in sampling_range:
      for phi in sampling_range:
        sites_cart_new = apply_rigid_body_shift(
          sites_cart=sites_cart,
          cm=cm, x=0,y=0,z=0, the=the, psi=psi, phi=phi)
        xray_structure = xray_structure.replace_sites_cart(
          new_sites=sites_cart_new)
        w = run_refine(
          rsr_simple_refiner = rsr_simple_refiner,
          xray_structure     = xray_structure,
          scorer             = sc,
          log                = log,
          weight             = flex.mean_default(weights, 1.0))
        weights.append(w)
        if(expload):
          for i in range(n_expload):
            xray_structure_ = xray_structure.deep_copy_scatterers()
            xray_structure_.shake_sites_in_place(mean_distance=1.0)
            run_refine(
              rsr_simple_refiner = rsr_simple_refiner,
              xray_structure     = xray_structure_,
              scorer             = sc,
              log                = log)
  if(log): print("Final target:", sc.target, file=log)
  return xray_structure.replace_sites_frac(new_sites=sc.sites_frac_best)
