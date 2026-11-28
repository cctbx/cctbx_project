from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
import mmtbx.model
from scitbx.array_family import flex
from libtbx.utils import null_out
from mmtbx.secondary_structure.build import ss_idealization as ssb
from cctbx import uctbx

"""
geometry_restraints.manager.plain_pair_sym_table (pairs of atoms within
plain_pairs_radius, used by ADP restraints) is made when it is asked for, from
the coordinates of the last update of nonbonded interactions, and not at
every such update: coordinate minimization of a large model updates
nonbonded interactions many times and never needs this table.
"""

def get_model():
  ph = ssb.secondary_structure_from_sequence(
    ssb.alpha_helix_str, "ACDEFGHIKLMNPQRSTVWY")
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart = ph.atoms().extract_xyz(), buffer_layer = 5)
  ph.atoms().set_xyz(box.sites_cart)
  model = mmtbx.model.manager(
    model_input      = None,
    pdb_hierarchy    = ph,
    crystal_symmetry = box.crystal_symmetry(),
    log              = null_out())
  model.process(make_restraints = True)
  return model

def as_list(table):
  result = []
  for i, d in enumerate(table):
    for j, ops in d.items():
      result.append((i, j, tuple([str(o) for o in ops])))
  return sorted(result)

def exercise():
  model = get_model()
  grm = model.get_restraints_manager().geometry
  assert grm.plain_pairs_radius == 5
  uc = model.crystal_symmetry().unit_cell()
  # count how many times the table is made
  n_made = []
  original = grm.update_plain_pair_sym_table
  def counting(sites_frac):
    n_made.append(1)
    return original(sites_frac = sites_frac)
  grm.update_plain_pair_sym_table = counting
  sites_cart = model.get_sites_cart().deep_copy()
  n_updates_start = grm.n_updates_pair_proxies
  cm = sites_cart.mean()
  for trial in range(3): # each move is larger than the nonbonded buffer
    sites_previous = sites_cart
    sites_cart = (sites_cart - cm)*1.08 + cm # model expands: fewer pairs
    assert sites_cart.max_distance(sites_previous) > 1.05
    grm.pair_proxies(sites_cart = sites_cart)
  assert grm.n_updates_pair_proxies == n_updates_start + 3
  assert len(n_made) == 0, len(n_made) # not made by updates
  table = grm.plain_pair_sym_table
  assert len(n_made) == 1, len(n_made) # made when asked for
  assert grm.plain_pair_sym_table is table and len(n_made) == 1 # and kept
  # it is the table for coordinates of the last update
  other = get_model().get_restraints_manager().geometry
  other.update_plain_pair_sym_table(sites_frac = uc.fractionalize(sites_cart))
  expected = as_list(other.plain_pair_sym_table)
  assert len(expected) > 500
  assert as_list(table) == expected
  other.update_plain_pair_sym_table(sites_frac=uc.fractionalize(sites_previous))
  assert len(as_list(other.plain_pair_sym_table)) > len(expected)
  # explicit update: made at once
  grm.update_plain_pair_sym_table(
    sites_frac = uc.fractionalize(model.get_sites_cart()))
  assert len(n_made) == 2
  assert as_list(grm.plain_pair_sym_table) != expected
  assert len(n_made) == 2
  # assignment, reset
  grm.plain_pair_sym_table = table
  assert grm.plain_pair_sym_table is table
  grm.pair_proxies(sites_cart = sites_cart + flex.vec3_double(
    sites_cart.size(), (3,3,3)))
  grm.reset_internals()
  assert grm.plain_pair_sym_table is None and len(n_made) == 2
  # no radius: no table
  grm.plain_pairs_radius = None
  grm.pair_proxies(sites_cart = sites_cart)
  assert grm.plain_pair_sym_table is None
  # selected manager, pickle
  model = get_model()
  grm = model.get_restraints_manager().geometry
  sites_cart = model.get_sites_cart()
  grm.pair_proxies(sites_cart = sites_cart)
  from libtbx import easy_pickle
  grm2 = easy_pickle.loads(easy_pickle.dumps(grm))
  assert as_list(grm2.plain_pair_sym_table) == as_list(grm.plain_pair_sym_table)
  sel = flex.bool([i%2 == 0 for i in range(sites_cart.size())])
  grm_sel = grm.select(selection = sel)
  grm_sel.pair_proxies(sites_cart = sites_cart.select(sel))
  n = len(as_list(grm_sel.plain_pair_sym_table))
  assert 0 < n < len(as_list(grm.plain_pair_sym_table))
  # ADP restraints use it
  rm = model.get_restraints_manager()
  e = rm.energies_adp_iso(
    xray_structure = model.get_xray_structure(),
    parameters     = mmtbx.refinement.adp_refinement.adp_restraints_master_params.extract().iso,
    use_u_local_only = False,
    use_hd = False,
    compute_gradients = True)
  assert e.gradients.size() == sites_cart.size()

if(__name__ == "__main__"):
  t0 = time.time()
  import mmtbx.refinement.adp_refinement
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
