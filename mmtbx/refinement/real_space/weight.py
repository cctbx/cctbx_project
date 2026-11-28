from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.refinement.real_space.individual_sites
import math

def search(refine, start, tolerance=0.2, max_trials=10, w_min=1.e-3,
           w_max=1.e+4):
  """
  Find the largest weight that satisfies geometry limits.

  refine(weight) does a trial refinement and returns score, the larger of
  (rmsd of bonds)/(its limit) and (rmsd of angles)/(its limit) of the refined
  model: limits are satisfied if score <= 1. The score grows with the weight.
  A weight is accepted if 1-tolerance <= score <= 1.

  The search is done on logarithms of weight and score. From the start value
  the weight is changed by a factor (2, or 4 if the score is more than twice
  off; from the third trial on the factor is extrapolated from the last two
  trials, between 2 and 100) until the acceptable range of scores is
  bracketed, then the bracket is narrowed by interpolation. A start value
  near the answer (the weight found in the previous macro-cycle) needs one or
  two trials.

  If the score is above 1 and lowering the weight does not lower it (less
  than 3 % for a twofold lower weight) the limits cannot be reached by
  changing the weight: the refinement did not converge in the iterations
  given (start model with poor geometry) or restraints conflict. The search
  stops; giving up the fit to the map would buy nothing.

  Returns list of trials [weight, score] in the order they were done and the
  index of the best one: the largest score not above 1; if there is none, the
  largest weight among trials with score within 5 % of the smallest score.
  """
  assert tolerance > 0 and tolerance < 1
  target = math.log(1 - tolerance/2.)
  log = lambda x: math.log(max(x, 1.e-9))
  trials = []
  lo, hi = None, None # [weight, score] closest to the target from each side
  w = min(w_max, max(w_min, start))
  while True:
    s = refine(w)
    trials.append([w, s])
    if(s <= 1 and s >= 1-tolerance): break
    if(s < 1):
      if(lo is None or w > lo[0]): lo = [w, s]
    else:
      if(hi is None or w < hi[0]): hi = [w, s]
    if(lo is not None and hi is not None and hi[0] <= lo[0]):
      # not monotonic (noise): keep the pair that brackets with larger weight
      hi = None
    if(len(trials) >= max_trials): break
    if(lo is not None and hi is not None):
      if(hi[0]/lo[0] < 1.1): break
      f = (target - log(lo[1]))/(log(hi[1]) - log(lo[1]))
      f = min(0.85, max(0.15, f))
      w = math.exp(log(lo[0]) + f*(log(hi[0]) - log(lo[0])))
    else:
      up = lo is not None
      if(up and w >= w_max): break
      if(not up and w <= w_min): break
      step = math.log(4. if (s < 0.5 or s > 2) else 2.)
      if(len(trials) > 1):
        (w1, s1), (w2, s2) = trials[-2], trials[-1]
        if(w1 != w2):
          slope = (log(s2) - log(s1))/(log(w2) - log(w1))
          if(slope > 0.1):
            step = abs(target - log(s2))/slope
            step = min(math.log(100.), max(math.log(2.), step))
          elif(up):
            step = math.log(8.) # weight too small to matter yet
          elif(s1 > 1 and abs(log(s2) - log(s1)) < 0.03):
            break # lowering the weight does not help
      if(not up): step = -step
      w = min(w_max, max(w_min, math.exp(log(w) + step)))
  ok = [i for i, t in enumerate(trials) if t[1] <= 1]
  if(len(ok) > 0):
    best = max(ok, key = lambda i: trials[i][1])
  else:
    s_min = min([t[1] for t in trials])
    best = max([i for i, t in enumerate(trials) if t[1] <= s_min*1.05],
      key = lambda i: trials[i][0])
  return trials, best

class run(object):
  def __init__(
      self,
      map_data,
      xray_structure,
      pdb_hierarchy,
      geometry_restraints_manager,
      gradients_method="fd",
      ncs_groups=None,
      rms_bonds_limit=0.0075,
      rms_angles_limit=1.4,
      real_space_gradients_delta=1./4,
      max_iterations = 100,
      range_size=10,
      n_ranges=10,
      default_weight=50,
      start_trial_weight_value=None,
      nproc=1,
      selection=None,
      selection_real_space=None,
      tolerance=0.2,
      max_trials=10,
      max_shift=mmtbx.refinement.real_space.individual_sites.default_max_shift,
      history=None,
      one_trial=False,
      must_satisfy=True,
      fragment_residues=10,
      n_fragments=10):
    """
Data/restraints weight for real-space refinement of individual sites: the
largest weight for which, after refinement of the whole model from its
current position, rmsd of bonds and angles from ideal values are within
rms_bonds_limit and rms_angles_limit (see search above). The rmsd are those
of non-hydrogen atoms, as in reported statistics; if only a part of the
model is refined (selection), of the refined atoms (restraints between
them): atoms that do not move tell nothing about the weight.

Trial refinements are the same refinement the caller is about to do
(max_iterations, gradients_method, selections, NCS constraints), so the
result of the trial with the chosen weight is the refined model: sites_cart,
rms_bonds, rms_angles.

start_trial_weight_value: where the search starts. If it is not given (and
  there is no history, see below) it is estimated with the same search done
  on n_fragments fragments of fragment_residues residues each, spread evenly
  over the model (always the same for a given model), each refined in a box
  with the atoms around it fixed. This is cheap, its cost does not depend on
  model size, and fragments are done in parallel (nproc; the result does not
  depend on nproc). Fragments do not predict the rmsd of the whole model
  well, so this is a first guess only.

one_trial: for a series of calls (macro-cycles). A trial refinement of the
  whole model costs as much as the refinement itself, so with one_trial only
  one is done per call - the one that has to be done anyway - and the weight
  for it is predicted from previous calls (history: list of [weight, score]
  of the refinements done before, most recent last; returned extended as
  self.history): score ~ weight**p, with p from the last two entries (0.7 if
  not available), aiming at score 0.9, the weight changed by a factor of
  four at most. The weight thus converges to its target over macro-cycles
  instead of within each. If must_satisfy is true (pass it for the last
  macro-cycle only) and the score is above 1, the search continues from the
  corrected weight until the limits are satisfied.

xray_structure is not modified. range_size, n_ranges and default_weight are
not used (they belonged to the previous algorithm).
    """
    self.msg_strings = []
    grm = geometry_restraints_manager.geometry
    sites_cart_start = xray_structure.sites_cart()
    if(selection is None):
      selection = flex.bool(xray_structure.scatterers().size(), True)
    self._refiner = mmtbx.refinement.real_space.individual_sites.simple(
      target_map                  = map_data,
      selection                   = selection,
      selection_real_space        = selection_real_space,
      real_space_gradients_delta  = real_space_gradients_delta,
      max_iterations              = max_iterations,
      geometry_restraints_manager = grm,
      gradients_method            = gradients_method,
      ncs_groups                  = ncs_groups,
      max_shift                   = max_shift)
    # rmsd over refined non-hydrogen atoms
    self._not_hd, self._grm_rmsd = None, grm
    use = ~xray_structure.hd_selection()
    if(use.count(True) == 0): use = ~use # hydrogen atoms only
    use = use & selection
    if(use.count(True) > 0 and use.count(False) > 0):
      self._not_hd = use
      self._grm_rmsd = grm.select(self._not_hd)
    self._results = {}
    def refine(weight):
      self._refiner.refine(
        weight         = weight,
        xray_structure = xray_structure.deep_copy_scatterers())
      sites_cart = self._refiner.sites_cart()
      b, a = self._rmsds(sites_cart = sites_cart)
      self._results[weight] = [sites_cart.deep_copy(), b, a,
        self._refiner.largest_shift]
      return max(b/rms_bonds_limit, a/rms_angles_limit)
    self.history = [] if history is None else [list(h) for h in history]
    start = start_trial_weight_value
    if(start is None and len(self.history) == 0):
      start = self._estimate_from_fragments(
        pdb_hierarchy        = pdb_hierarchy,
        fragment_residues    = fragment_residues,
        n_fragments          = n_fragments,
        nproc                = nproc,
        tolerance            = tolerance,
        max_trials           = max_trials,
        setup = dict(
          xray_structure       = xray_structure,
          map_data             = map_data,
          grm                  = grm,
          selection            = selection,
          selection_real_space = selection_real_space,
          limits               = (rms_bonds_limit, rms_angles_limit),
          refiner_args         = dict(
            real_space_gradients_delta = real_space_gradients_delta,
            max_iterations             = max_iterations,
            gradients_method           = gradients_method,
            max_shift                  = max_shift)))
    if(start is None): start = 1.
    b, a = self._rmsds(sites_cart = sites_cart_start)
    self.msg_strings.append(
      "start: rmsd(bonds)=%6.4f rmsd(angles)=%5.2f; limits: %6.4f %5.2f"%(
      b, a, rms_bonds_limit, rms_angles_limit))
    def predicted(w, s):
      # weight expected to give score 0.9, from [w, s] and the entry before
      p = 0.7
      if(len(self.history) > 1):
        (w1, s1), (w2, s2) = self.history[-2], self.history[-1]
        if(w1 > 0 and w2 > 0 and s1 > 0 and s2 > 0 and
           abs(math.log(w2/w1)) > 0.05):
          p = min(2., max(0.3, math.log(s2/s1)/math.log(w2/w1)))
      f = (0.9/max(s, 1.e-6))**(1./p)
      return min(1.e+4, max(1.e-3, w*min(4., max(0.25, f))))
    if(one_trial):
      w = start
      if(len(self.history) > 0):
        w = predicted(*self.history[-1])
      s = refine(w)
      self.history.append([w, s])
      trials, best = [[w, s]], 0
      if(s > 1 and must_satisfy):
        more, junk = search(
          refine     = refine,
          start      = predicted(w, s),
          tolerance  = tolerance,
          max_trials = max_trials)
        trials = trials + more
        ok = [i for i, t in enumerate(trials) if t[1] <= 1]
        if(len(ok) > 0): best = max(ok, key = lambda i: trials[i][1])
        else: best = min(range(len(trials)), key = lambda i: trials[i][1])
        self.history.append(list(trials[best]))
    else:
      trials, best = search(
        refine     = refine,
        start      = start,
        tolerance  = tolerance,
        max_trials = max_trials)
      self.history.append(list(trials[best]))
    for w, s in trials:
      self.msg_strings.append(
        "trial weight: %10.4f rmsd(bonds)=%6.4f rmsd(angles)=%5.2f"%(
        w, self._results[w][1], self._results[w][2])+
        " largest shift: %5.2f"%self._results[w][3])
    self.n_trials = len(trials)
    self.trials = trials
    self.weight = trials[best][0]
    self.sites_cart, self.rms_bonds, self.rms_angles, self.largest_shift = \
      self._results[self.weight]
    self.limits_satisfied = trials[best][1] <= 1
    self.msg_strings.append("trial refinements: %d"%self.n_trials)
    self.msg_strings.append("overall best weight: %9.4f"%self.weight)
    if(not self.limits_satisfied):
      self.msg_strings.append(
        "rmsd limits are not reached with any weight tried")
    del self._results, self._refiner, self._grm_rmsd
    self._fragments = None

  def _estimate_from_fragments(self, pdb_hierarchy, fragment_residues,
        n_fragments, nproc, tolerance, max_trials, setup):
    """
    First guess of the weight: median of weights found for fragments. None
    if there are no suitable fragments.
    """
    runs = []
    selection = setup["selection"]
    for chain in pdb_hierarchy.chains():
      if(not (chain.is_protein() or chain.is_na())): continue
      rgs = chain.residue_groups()
      for i in range(0, len(rgs), fragment_residues):
        sel = flex.size_t()
        for rg in rgs[i:i+fragment_residues]:
          sel.extend(rg.atoms().extract_i_seq())
        # whole fragment must be refinable, and not too small
        if(sel.size() >= 20 and selection.select(sel).all_eq(True)):
          runs.append(sel)
    if(len(runs) == 0): return None
    if(len(runs) > n_fragments):
      step = len(runs)/float(n_fragments)
      runs = [runs[int(i*step + step/2)] for i in range(n_fragments)]
    self._fragments = [self._setup_fragment(iselection = r, **setup)
                       for r in runs]
    argss = [[i, tolerance, max_trials] for i in range(len(runs))]
    if(nproc > 1 and len(argss) > 1):
      from libtbx import easy_mp
      results = easy_mp.pool_map(
        processes  = min(nproc, len(argss)),
        fixed_func = self._fragment_search,
        args       = argss)
    else:
      results = [self._fragment_search(a) for a in argss]
    self._fragments = None
    weights = sorted([r[0] for r in results])
    n = len(weights)
    w = math.exp(0.5*(math.log(weights[(n-1)//2]) + math.log(weights[n//2])))
    self.msg_strings.append(
      "first guess from %d fragments (%d trial refinements): %s -> %.4f"%(
      n, sum([r[1] for r in results]),
      " ".join(["%.3g"%x for x in weights]), w))
    return w

  def _setup_fragment(self, iselection, xray_structure, map_data, grm,
        selection, selection_real_space, limits, refiner_args,
        buffer_radius=5, box_cushion=2):
    """
    Box with the fragment and atoms within buffer_radius of it (fixed).
    """
    from mmtbx import utils
    n = xray_structure.scatterers().size()
    fragment = flex.bool(n, iselection)
    within = xray_structure.selection_within(
      radius = buffer_radius, selection = fragment)
    box = utils.extract_box_around_model_and_map(
      xray_structure = xray_structure.select(within),
      map_data       = map_data,
      box_cushion    = box_cushion)
    geo = grm.select(within).discard_symmetry(
      new_unit_cell = box.xray_structure_box.unit_cell())
    geo.shift_sites_cart(box.shift_cart)
    sel = fragment.select(within)
    sel_rs = None
    if(selection_real_space is not None):
      sel_rs = selection_real_space.select(within)
    refiner = mmtbx.refinement.real_space.individual_sites.simple(
      target_map                  = box.map_box,
      selection                   = sel,
      selection_real_space        = sel_rs,
      geometry_restraints_manager = geo,
      **refiner_args)
    use = sel & ~box.xray_structure_box.hd_selection()
    if(use.count(True) == 0): use = sel
    return dict(refiner = refiner, xrs = box.xray_structure_box, use = use,
      grm_rmsd = geo.select(use), limits = limits)

  def _fragment_search(self, args):
    """
    Search for the weight on one fragment: weight and number of trials.
    """
    i, tolerance, max_trials = args
    c = self._fragments[i]
    def refine(weight):
      c["refiner"].refine(
        weight = weight, xray_structure = c["xrs"].deep_copy_scatterers())
      es = c["grm_rmsd"].energies_sites(
        sites_cart = c["refiner"].sites_cart().select(c["use"]))
      return max(es.bond_deviations()[2]/c["limits"][0],
                 es.angle_deviations()[2]/c["limits"][1])
    trials, best = search(refine = refine, start = 1.,
      tolerance = tolerance, max_trials = max_trials)
    return trials[best][0], len(trials)

  def _rmsds(self, sites_cart):
    if(self._not_hd is not None):
      sites_cart = sites_cart.select(self._not_hd)
    es = self._grm_rmsd.energies_sites(sites_cart = sites_cart)
    return es.bond_deviations()[2], es.angle_deviations()[2]

  def show(self, log, prefix=""):
    for m in self.msg_strings:
      print("%s %s"%(prefix, m), file=log)
