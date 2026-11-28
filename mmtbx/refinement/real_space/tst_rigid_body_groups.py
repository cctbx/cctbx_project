from __future__ import absolute_import, division, print_function
import time, random
import iotbx.pdb
import mmtbx.model
import mmtbx.refinement.real_space.rigid_body as rigid_body
from mmtbx.secondary_structure.build import ss_idealization as ssb
from cctbx import uctbx
from iotbx import map_manager, map_model_manager
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import null_out
from six.moves import cStringIO as StringIO

chain_ids = "ABCDEFGHIJ"

def assembly(template, sequence, shifts, rotation=None):
  """
  Hierarchy made of copies of one secondary-structure element, one chain per
  shift (translation vector). rotation: applied to the element first.
  """
  ph = ssb.secondary_structure_from_sequence(template, sequence)
  sites = ph.atoms().extract_xyz()
  sites = sites - sites.mean()
  if(rotation is not None): sites = rotation.elems * sites
  lines = []
  for i, shift in enumerate(shifts):
    ph.atoms().set_xyz(sites + tuple(shift))
    for chain in ph.chains(): chain.id = chain_ids[i]
    lines.extend([l for l in ph.as_pdb_string().splitlines()
                  if l.startswith("ATOM")])
    lines.append("TER")
  return iotbx.pdb.input(
    source_info=None, lines=lines).construct_hierarchy(), sites

def get_mmm(ph, d_min, buffer_layer, n_chains_in_model=None, b_iso=None):
  """
  Map calculated from all chains of ph; model has the first n_chains_in_model
  of them (all by default). b_iso: B-factor of each chain used for the map.
  """
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart = ph.atoms().extract_xyz(), buffer_layer = buffer_layer)
  ph.atoms().set_xyz(box.sites_cart)
  model = mmtbx.model.manager(
    model_input      = None,
    pdb_hierarchy    = ph,
    crystal_symmetry = box.crystal_symmetry(),
    log              = null_out())
  model.setup_scattering_dictionaries(scattering_table = "electron")
  if(b_iso is not None):
    values = flex.double()
    for chain, b in zip(ph.chains(), b_iso):
      values.extend(flex.double(chain.atoms_size(), b))
    model.set_b_iso(values = values)
  fc = model.get_xray_structure().structure_factors(d_min = d_min).f_calc()
  fft_map = fc.fft_map(resolution_factor = 0.25)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  mm = map_manager.map_manager(
    map_data                   = map_data,
    unit_cell_grid             = map_data.all(),
    unit_cell_crystal_symmetry = model.crystal_symmetry(),
    wrapping                   = False)
  if(n_chains_in_model is not None):
    model = model.select(model.selection(" or ".join(
      ["chain %s"%c for c in chain_ids[:n_chains_in_model]])))
  return map_model_manager.map_model_manager(model = model, map_manager = mm)

def chain_selections(model):
  return [flex.bool(model.size(), c.atoms().extract_i_seq())
          for c in model.get_hierarchy().chains()]

def move(sites_cart, selections, angle, shift, seed=1, direction=None):
  """
  Rotate each selection by angle (degrees) about a random axis through its
  center and shift it by shift (A) along direction (random if None).
  """
  random.seed(seed)
  result = sites_cart.deep_copy()
  for sel in selections:
    s = sites_cart.select(sel)
    cm = matrix.col(s.mean())
    axis = matrix.col([random.gauss(0,1) for i in range(3)]).normalize()
    r = axis.axis_and_angle_as_r3_rotation_matrix(angle=angle, deg=True)
    t = direction
    if(t is None): t = [random.gauss(0,1) for i in range(3)]
    t = matrix.col(t).normalize()*shift
    result.set_selected(sel, r.elems*(s-cm) + cm + t)
  return result

def distances(sites_1, sites_2, selections):
  d = flex.sqrt((sites_1 - sites_2).dot())
  return flex.double([flex.mean(d.select(s)) for s in selections])

def exercise_bundle(d_min = 4.0):
  """
  Four helices side by side, each a rigid group. In place: not moved at all.
  Displaced by 3 degrees and 1.5 A: come back. Same result for nproc 1 and 2.
  """
  ph, junk = assembly(
    template = ssb.alpha_helix_str,
    sequence = "ACDEFGHIKLMNPQRSTVWY",
    shifts   = [(0,0,0), (12,0,0), (0,12,0), (12,12,0)])
  mmm = get_mmm(ph = ph, d_min = d_min, buffer_layer = 8)
  model = mmm.model()
  answer = model.get_sites_cart().deep_copy()
  sels = chain_selections(model)
  assert len(sels) == 4
  # in place
  o = rigid_body.refine_groups(
    map_model_manager = mmm, selections = sels, d_min = d_min)
  assert o.result == ["not moved",]*4, o.result
  assert flex.max(distances(model.get_sites_cart(), answer, sels)) < 1.e-6
  # displaced
  start = move(answer, sels, angle = 3, shift = 1.5)
  assert flex.min(distances(start, answer, sels)) > 1.4
  results = []
  for nproc in [1, 2]:
    model.set_sites_cart(start)
    log = StringIO()
    o = rigid_body.refine_groups(
      map_model_manager = mmm, selections = sels, d_min = d_min,
      nproc = nproc, log = log)
    d = distances(model.get_sites_cart(), answer, sels)
    print("bundle, nproc=%d: distance to answer (mean, max): %6.3f %6.3f"%(
      nproc, flex.mean(d), flex.max(d)), o.result)
    assert flex.max(d) < 0.3, list(d)
    assert o.result.count("not moved") == 0, o.result
    assert len([l for l in log.getvalue().splitlines()
                if l.count("moved (mean):")]) == 4
    results.append(model.get_sites_cart().deep_copy())
  assert flex.max(flex.sqrt((results[0]-results[1]).dot())) < 1.e-6

def exercise_far(d_min = 4.0):
  """
  Whole model as one group, 30 degrees and 12 A away. The local search alone
  does not get there; the long-range multi-zone search does, and its result
  is taken because it is much better.
  """
  ph, junk = assembly(
    template = ssb.alpha_helix_str,
    sequence = "ACDEFGHIKLMNPQRSTVWY",
    shifts   = [(0,0,0), (12,0,0), (0,12,0), (12,12,0)])
  mmm = get_mmm(ph = ph, d_min = d_min, buffer_layer = 25)
  model = mmm.model()
  answer = model.get_sites_cart().deep_copy()
  sels = [flex.bool(model.size(), True)]
  start = move(answer, sels, angle = 30, shift = 12)
  d_start = distances(start, answer, sels)[0]
  assert d_start > 11, d_start
  for max_shift, converged in [(None, False), (5., True)]:
    model.set_sites_cart(start)
    o = rigid_body.refine_groups(
      map_model_manager = mmm, selections = sels, d_min = d_min,
      max_shift = max_shift)
    d = distances(model.get_sites_cart(), answer, sels)[0]
    print("far, max_shift=%s: distance to answer %6.3f -> %6.3f"%(
      str(max_shift), d_start, d), o.result)
    if(converged):
      assert d < 0.5, d
      assert o.result == ["long-range"], o.result
    else:
      assert d > 3, d

def exercise_layers(d_min = 3.5, expect_failure = False):
  """
  Stack of identical strands 4.8 A apart (as layers of a fibril); the map has
  one more layer at each end than the model. A neighbouring layer, modelled
  or not, fits a strand as well as its own, or better: here the middle layer
  of the model has much stronger density than the others (low B), and the
  two layers that are not modelled have somewhat stronger density. Strands
  in place are not moved; strands displaced by 1.5 A towards the next layer
  (and 2 degrees) come back to their own layer. Without the rule that a
  long move must pay the end strands go to the layers that are not modelled;
  without the rejection of overlaps the neighbours of the middle strand go
  on top of it (checked by switching each rule off).
  """
  ph, sites = assembly(
    template = ssb.beta_pdb_str, sequence = "VKVTVEVKVT", shifts = [(0,0,0)])
  # stacking direction: along C=O of the first residue, perpendicular to the
  # strand
  atoms = ph.atoms()
  c = [a for a in atoms if a.name == " C  "][0]
  o = [a for a in atoms if a.name == " O  "][0]
  ca = [a for a in atoms if a.name == " CA "]
  axis = (matrix.col(ca[-1].xyz) - matrix.col(ca[0].xyz)).normalize()
  v = matrix.col(o.xyz) - matrix.col(c.xyz)
  v = (v - axis*v.dot(axis)).normalize()
  # turn the strand so that the stacking direction is z (as in maps of
  # fibrils, and a direction the translation grid search samples)
  z = matrix.col((0,0,1))
  rotation = v.cross(z).normalize().axis_and_angle_as_r3_rotation_matrix(
    angle = v.angle(z, deg=True), deg = True)
  v = z
  n_model = 5
  # model layers first, then the two that are in the map only
  layers = list(range(n_model)) + [-1, n_model]
  ph, junk = assembly(
    template = ssb.beta_pdb_str,
    sequence = "VKVTVEVKVT",
    shifts   = [(v*4.8*i).elems for i in layers],
    rotation = rotation)
  mmm = get_mmm(
    ph = ph, d_min = d_min, buffer_layer = 8, n_chains_in_model = n_model,
    b_iso = [40, 40, 5, 40, 40] + [36, 36])
  model = mmm.model()
  answer = model.get_sites_cart().deep_copy()
  sels = chain_selections(model)
  assert len(sels) == n_model
  # in place
  o = rigid_body.refine_groups(
    map_model_manager = mmm, selections = sels, d_min = d_min)
  d = distances(model.get_sites_cart(), answer, sels)
  print("layers, in place: moved (max): %6.3f"%flex.max(d), o.result)
  d_in_place = d
  if(not expect_failure): assert flex.max(d) < 0.1, list(d)
  # displaced towards the next layer
  start = move(answer, sels, angle = 2, shift = 1.5, direction = v.elems)
  assert flex.min(distances(start, answer, sels)) > 1.4
  model.set_sites_cart(start)
  o = rigid_body.refine_groups(
    map_model_manager = mmm, selections = sels, d_min = d_min)
  d = distances(model.get_sites_cart(), answer, sels)
  print("layers, displaced: distance to answer (mean, max): %6.3f %6.3f"%(
    flex.mean(d), flex.max(d)), o.result)
  if(not expect_failure): assert flex.max(d) < 0.5, list(d)
  return flex.max(d_in_place), flex.max(d)

def exercise_layers_rules_matter():
  """
  exercise_layers fails if either safety rule of refine_groups is off.
  """
  original = rigid_body.refine_groups
  for off in [dict(min_gain = -1.e9), dict(overlap_tolerance = 1.e9)]:
    class modified(original):
      def __init__(self, **kwds):
        kwds.update(off)
        original.__init__(self, **kwds)
    rigid_body.refine_groups = modified
    try:
      d_in_place, d_displaced = exercise_layers(expect_failure = True)
    finally:
      rigid_body.refine_groups = original
    assert max(d_in_place, d_displaced) > 4, [off, d_in_place, d_displaced]

def exercise_translation_search():
  """
  translation_search and protocol_1: same position found with all atoms and
  with a subset (max_sites); the search range is honoured.
  """
  ph, junk = assembly(
    template = ssb.alpha_helix_str,
    sequence = "ACDEFGHIKLMNPQRSTVWY",
    shifts   = [(0,0,0), (10,0,0)])
  mmm = get_mmm(ph = ph, d_min = 4.0, buffer_layer = 10)
  model = mmm.model()
  map_data = mmm.map_manager().map_data()
  answer = model.get_sites_cart().deep_copy()
  shifted = answer + (2.0, 0, -2.0)
  for max_sites in [None, 100]:
    model.set_sites_cart(shifted)
    shift = rigid_body.translation_search(
      model = model, map_data = map_data, max_sites = max_sites)
    assert [round(s, 6) for s in shift] == [-2.0, 0, 2.0], shift
    assert flex.max(flex.sqrt((model.get_sites_cart()-answer).dot())) < 1.e-6
  # out of range: best within range is taken
  model.set_sites_cart(shifted)
  shift = rigid_body.translation_search(
    model = model, map_data = map_data, max_shift = 1.0)
  assert max([abs(s) for s in shift]) <= 1.0 + 1.e-6, shift
  # protocol_1 with rotation, all atoms and subset
  sels = [flex.bool(model.size(), True)]
  start = move(answer, sels, angle = 5, shift = 2)
  for max_sites in [None, 100]:
    model.set_sites_cart(start)
    rigid_body.protocol_1(
      model = model, map_data = map_data, max_sites = max_sites)
    d = distances(model.get_sites_cart(), answer, sels)[0]
    print("protocol_1, max_sites=%s: distance to answer %6.3f"%(
      str(max_sites), d))
    assert d < 0.2, d

def exercise_empty_and_tiny_groups(d_min = 4.0):
  """
  A selection that selects nothing is skipped (it used to crash or mix up
  groups). A group of one atom is refined locally only: for so few atoms a
  strong peak a few A away always "pays" (it used to hop there).
  """
  ph, junk = assembly(
    template = ssb.alpha_helix_str,
    sequence = "ACDEFGHIKLMNPQRSTVWY",
    shifts   = [(0,0,0), (12,0,0)])
  mmm = get_mmm(ph = ph, d_min = d_min, buffer_layer = 8)
  model = mmm.model()
  answer = model.get_sites_cart().deep_copy()
  n = model.size()
  chains = chain_selections(model)
  # atom with the weakest density in chain B, alone in a group
  md = mmm.map_manager().map_data()
  uc = model.crystal_symmetry().unit_cell()
  values = [md.tricubic_interpolation(uc.fractionalize(s)) for s in answer]
  i_weak = min([i for i in range(n) if chains[1][i]], key=lambda i: values[i])
  one = flex.bool(n, False)
  one[i_weak] = True
  empty = flex.bool(n, False)
  sels = [chains[0], empty, chains[1] & ~one, one]
  start = move(answer, [chains[0], chains[1]], angle = 3, shift = 1.5)
  model.set_sites_cart(start)
  o = rigid_body.refine_groups(
    map_model_manager = mmm, selections = sels, d_min = d_min)
  assert len(o.result) == 4 and o.result[1] == "not moved", o.result
  d = flex.sqrt((model.get_sites_cart() - answer).dot())
  assert flex.mean(d.select(sels[0])) < 0.3 and flex.mean(d.select(sels[2])) < 0.3
  assert o.result[3] in ["local", "not moved"], o.result
  moved = flex.sqrt((model.get_sites_cart() - start).dot())[i_weak]
  assert moved < 2.5, moved

if(__name__ == "__main__"):
  t0 = time.time()
  exercise_empty_and_tiny_groups()
  exercise_translation_search()
  exercise_bundle()
  exercise_far()
  exercise_layers()
  exercise_layers_rules_matter()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
