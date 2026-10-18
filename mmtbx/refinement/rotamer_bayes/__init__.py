"""
rotamer_bayes: decide side-chain rotamer changes by Bayes' rule.

For each side chain, the hypotheses are the rotamer wells of the Top8000
distribution (priors.py), with prior probabilities from that distribution.
The evidence for moving the side chain from its current conformer to a
candidate in another well is the change in the log-likelihood of the
diffraction data (evidence.py) for the rigid move (rest of the model, mask,
scales and error model fixed), so the log posterior odds are

    ln[P(b|data)/P(a|data)] = dLL(a->b) + ln(prior_b) - ln(prior_a)

and a change is made when the best candidate's log odds exceed a threshold
(0 by default: the posterior favours it). Neighbours are not relaxed, so the
evidence for a change is a lower bound on what refinement would give: the
test is conservative.

Candidates start at each well's canonical chi angles and their chi angles
are then searched within tune_range of them, on a cheap approximation of
dLL, so that they are compared with the (refined) current conformer on a
fair footing. Candidates whose heavy atoms come closer than clash_distance
to other heavy atoms (crystal symmetry included) are left out.

Independent of phenix.refine: needs an mmtbx.f_model.manager (any target)
and the matching pdb hierarchy.

  result = mmtbx.refinement.rotamer_bayes.run(fmodel, pdb_hierarchy)
  result.show()
  n = result.apply(pdb_hierarchy, fmodel.xray_structure)
"""
from __future__ import absolute_import, division, print_function
import math, sys, time
import numpy as np
import iotbx.phil
from libtbx import group_args
from libtbx.utils import Sorry
from libtbx import Auto
from cctbx.array_family import flex
from mmtbx.rotamer.sidechain_angles import SidechainAngles
from mmtbx.refinement.rotamer_bayes import priors as priors_module
from mmtbx.refinement.rotamer_bayes import side_chain as side_chain_module
from mmtbx.refinement.rotamer_bayes import evidence as evidence_module
from mmtbx.refinement.rotamer_bayes import ncs as ncs_module

master_params_str = """
p_min = 0.001
  .type = float
  .help = Wells with prior probability below this are not tried
threshold = 0.
  .type = float
  .help = Change when the log posterior odds of the best candidate over \
          the current conformer exceed this
screen = 10.
  .type = float
  .help = Compute the exact dLL for candidates whose predicted log odds \
          exceed -screen
calibration_sample = 0
  .type = int
  .help = Also compute the exact dLL for this many other candidates, at \
          random (for calibrating the prediction)
random_seed = 0
  .type = int
exact = False
  .type = bool
  .help = Compute the exact dLL (direct summation, target, local mask) for \
          candidates passing the screen and decide on it. The prediction \
          (with its mask term) agrees with it closely (CC ~0.999 in tests) \
          and is much cheaper, so by default decisions use the prediction.
min_dll = 0.
  .type = float
  .help = A change is made only if the data themselves favour it (dLL > \
          min_dll): the prior can veto a change but not drive one
min_information = 10.
  .type = float
  .help = A change is made only if the expected separation of the two \
          hypotheses, I_eff = 4 alpha^2/beta I (alpha, beta from the \
          self-check, 1/2 and 1 if not estimated), is at least this. \
          Residues whose side chains the data cannot define are left alone.
mask_correction = True
  .type = bool
mask_predict_spacing = 0.5
  .type = float
mask_cell = 1.5
  .type = float
  .help = Size (A) of the cells into which changes of the mask are gathered \
          (as pseudo-atoms)
  .help = Grid spacing (A) of the local mask in the prediction
  .help = Include the change of the bulk-solvent mask in the exact dLL \
          (recomputed on a local grid). Without it every change is biased \
          against: new atoms land in solvent, the old sites stay excluded.
ncs_pooling = True
  .type = bool
  .help = Pool the evidence of NCS-related copies of a residue: the prior \
          over their wells is a mixture of "all copies share one well" \
          (probability rho, estimated from the structure) and "independent"
stubs = True
  .type = bool
  .help = Treat side chains whose moving atoms all have zero occupancy as \
          stubs to rebuild: every well is a candidate, the evidence is the \
          dLL for adding the side chain, and the most probable well is built \
          if its posterior probability exceeds stub_threshold
stub_threshold = 0.9
  .type = float
removal_exact = False
  .type = bool
  .help = Also compute the exact dLL for removing each tested side chain \
          (the atoms that the chi angles move). The prediction is always made
  .help = Decide on the exact dLL where computed (else on the prediction)
tune_range = 15.
  .type = float
tune_step = 5.
  .type = float
tune_passes = 2
  .type = int
clash_distance = 2.5
  .type = float
  .help = Minimum distance (A) between candidate heavy atoms and other \
          heavy atoms
bound_distance = 2.2
  .type = float
  .help = Residues with side-chain heavy atoms this close to other heavy \
          atoms (links, disulfides) are not tried
clash_with_water = False
  .type = bool
self_check = True
  .type = bool
  .help = Make no changes unless the evidence behaves as expected: for \
          wrong candidates dLL ~ N(-alpha I, beta I) with alpha ~ 1/2 and \
          beta <~ 1 (I the information of the swap). alpha and beta are \
          estimated robustly from all candidates (most are wrong). In \
          poor models beta is inflated (local changes absorb other errors).
self_check_alpha = 0.4 0.6
  .type = floats(size=2)
self_check_beta_max = Auto
  .type = float
  .help = Largest beta accepted by the self-check. Auto: 1.5 for llgi, 3 \
          for ml and mlhl, whose likelihood leaves out the measurement \
          errors, so that dLL varies about twice as much as the curvature \
          predicts even for a correct model (beta ~1.1-1.3 in tests)
self_check_min_candidates = 200
  .type = int
n_shells = 40
  .type = int
"""

def master_params():
  return iotbx.phil.parse(master_params_str)

SKIP_RESNAMES = ("ALA", "GLY", "PRO")

SUPPORTED_TARGETS = ("ml", "mlhl", "llgi")

def check_fmodel(fmodel):
  """ The evidence is a change of log-likelihood with F_model changing
  reflection by reflection: likelihood targets only, no twinning. """
  if(fmodel.is_twin_fmodel_manager()):
    raise Sorry("rotamer_bayes: twinned data are not supported (a change of "
      "the model affects twin-related reflections together).")
  name = fmodel.target_name
  if(name not in SUPPORTED_TARGETS):
    if(fmodel.target_attributes().family == "ls"):
      raise Sorry("rotamer_bayes: target %s is least squares, not a "
        "likelihood; the evidence would have no probabilistic scale. Use one "
        "of: %s." % (name, ", ".join(SUPPORTED_TARGETS)))
    raise Sorry("rotamer_bayes: target %s is not supported (tested: %s)."
      % (name, ", ".join(SUPPORTED_TARGETS)))

class candidate(object):
  def __init__(self, well, xyz):
    self.well = well.name
    self.prior = well.mass
    self.chis = list(well.chis)
    self.xyz = xyz            # all atoms of the residue
    self.clash = False
    self.search = None        # search approximation of dLL
    self.predicted = None     # second-order prediction of dLL
    self.information = None   # expected separation of dLL (curvature)
    self.information_size = None  # the same from the f^2 law (size only)
    self.exact = None         # rigid-move dLL
    self.log_odds = None
    self.posterior = None     # stubs: posterior probability of the well

  def dll(self, use_exact=True):
    if(use_exact and self.exact is not None): return self.exact
    return self.predicted + getattr(self, "offset", 0.)

class record(object):
  def __init__(self, label, sc, current_well, current_prior, chis):
    self.label = label
    self.sc = sc
    self.current_well = current_well
    self.current_prior = current_prior
    self.current_chis = chis
    self.candidates = []
    self.best = None
    self.accept = False
    self.posterior_current = None
    self.stub = False           # side chain absent (zero occupancy): rebuild
    self.posterior_best = None  # stubs: posterior of the best well

class run(object):

  def __init__(self, fmodel, pdb_hierarchy, selection=None, params=None,
               log=None, ncs_groups=None):
    """ ncs_groups: lists of NCS-equivalent residues, as (chain id, resseq,
    icode) (stripped); by default from the NCS search of the hierarchy """
    if(params is None): params = master_params().extract()
    if(log is None): log = sys.stdout
    self.params = params
    self.log = log
    check_fmodel(fmodel)
    self.target_name = fmodel.target_name
    t0 = time.time()
    self.ev = evidence_module.evidence(fmodel=fmodel,
      n_shells=params.n_shells)
    self.priors = priors_module.rotamer_priors()
    # Expected discrimination from size alone: LLG ~ K f^2 (f the fraction
    # of the scattering accounted for), so placing a part of size d right
    # instead of wrong changes the LLG by ~ 2 LLG_total d/f (d, f as sums of
    # (occ Z)^2 over non-hydrogen atoms).
    self.llg_total = self.ev.llg_total()
    self.model_scattering = self.ev.scattering_power()
    sa = SidechainAngles(False)
    atoms = pdb_hierarchy.atoms()
    sites = np.array(list(atoms.extract_xyz()))
    # clash partners: heavy atoms that are present (zero-occupancy atoms,
    # e.g. side chains of stubs, are not)
    keep = np.array([a.element.strip().upper() not in ("H", "D") and
      a.occ >= 0.01 for a in atoms])
    if(not params.clash_with_water):
      water = pdb_hierarchy.atom_selection_cache().selection("water")
      keep &= ~np.array(list(water))
    self.env = side_chain_module.environment(
      crystal_symmetry=fmodel.xray_structure.crystal_symmetry(),
      sites_cart=sites, keep=keep)
    if(params.mask_correction):
      self.ev.setup_mask(self.env, cell=params.mask_cell)
    self.records = []
    self.skipped = {}
    self._sample = set()
    for chain in pdb_hierarchy.only_model().chains():
      if(not chain.is_protein()): continue
      for rg in chain.residue_groups():
        if(len(rg.atom_groups()) != 1):
          self._skip("alternative conformations"); continue
        ag = rg.atom_groups()[0]
        resname = ag.resname.strip().upper()
        if(resname in SKIP_RESNAMES): continue
        if(selection is not None and
           not selection[ag.atoms()[0].i_seq]): continue
        wells = self.priors.wells(resname)
        if(wells is None):
          self._skip("no rotamer table"); continue
        sc = side_chain_module.side_chain(residue=ag, sidechain_angles=sa)
        if(not sc.ok):
          self._skip(sc.reason); continue
        if(sc.moving.sum() == 0): continue
        occ = np.array([a.occ for a in ag.atoms()])
        stub = bool(np.all(occ[sc.moving & ~sc.is_h] < 0.01))
        if(stub and not params.stubs):
          self._skip("zero occupancy"); continue
        if("CA" not in sc.names):
          self._skip("no CA"); continue
        sc.near = self.env.near(sc.xyz[sc.names.index("CA")], 13.,
          sc.i_seqs)
        if(not stub and
           self._min_distance(sc, sc.xyz) < params.bound_distance):
          self._skip("bound side chain"); continue
        chis = sc.chis()
        name = "STUB" if stub else self.priors.classify(resname, chis)
        label = "%s %s%s %s" % (chain.id, rg.resseq, rg.icode.strip(),
          resname)
        rec = record(label, sc, name, 1. if stub else (wells[name].mass
          if name in wells else priors_module.OUTLIER_THRESHOLD), chis)
        rec.stub = stub
        rec.key = (chain.id, rg.resseq.strip(), rg.icode.strip())
        rec.occ_old = occ[sc.moving]
        rec.occ_new = np.ones(int(sc.moving.sum())) if stub else rec.occ_old
        rec.moving_scattering = self.ev.scattering_power(
          sc.i_seqs[sc.moving])
        for w in wells.values():
          if(w.name in (priors_module.OUTLIER, priors_module.UNCLASSIFIED)):
            continue
          if(w.name == name or w.mass < params.p_min): continue
          rec.candidates.append(candidate(w, sc.set_chis(w.chis)))
        self.records.append(rec)
    self.time_setup = time.time() - t0
    t0 = time.time()
    self._fit_search_scale()
    for rec in self.records:
      self._tune(rec)
    self.time_search = time.time() - t0
    t0 = time.time()
    self._predict()
    self.time_predict = time.time() - t0
    self._self_check()
    t0 = time.time()
    self._exact_and_decide()
    self.time_exact = time.time() - t0
    self.ncs_rho = None
    self.ncs_n_groups = 0
    if(params.ncs_pooling):
      self._ncs_decide(pdb_hierarchy, ncs_groups)
    if(params.self_check and self.self_check_passed is False):
      for rec in self.records:
        if(not rec.stub): rec.accept = False

  def _self_check(self):
    """ Robust alpha, beta of dLL ~ N(-alpha I, beta I) over all candidates
    (median of dLL/I; MAD of the standardised residuals). """
    p = self.params
    d, inf = [], []
    for rec in self.records:
      if(rec.stub): continue
      for c in rec.candidates:
        if(c.clash or c.information is None or c.information <= 0): continue
        d.append(c.predicted); inf.append(c.information)
    self.alpha = self.beta = None
    self.self_check_passed = None
    if(len(d) < p.self_check_min_candidates): return
    d, inf = np.array(d), np.array(inf)
    self.alpha = float(-np.median(d / inf))
    z = (d + self.alpha * inf) / np.sqrt(inf)
    self.beta = float((np.median(np.abs(z - np.median(z))) * 1.4826) ** 2)
    lo, hi = p.self_check_alpha
    self.beta_max = p.self_check_beta_max
    if(self.beta_max is None or self.beta_max is Auto):
      self.beta_max = 1.5 if self.target_name == "llgi" else 3.
    self.self_check_passed = (lo <= self.alpha <= hi and
      self.beta <= self.beta_max)

  def _information_scale(self):
    """ I_eff / I = 4 alpha^2 / beta (theory 1) """
    if(self.alpha is None or self.beta is None or self.beta <= 0): return 1.
    return 4 * self.alpha ** 2 / self.beta

  def _gates(self, c):
    """ the data favour the change, with enough information to tell """
    p = self.params
    return (c.dll(p.exact) > p.min_dll and c.information is not None and
      self._information_scale() * c.information >= p.min_information)

  def _skip(self, reason):
    self.skipped[reason] = self.skipped.get(reason, 0) + 1

  def _min_distances(self, sc, xyz):
    """ _min_distance for T sets of residue sites (T, n_atoms, 3) """
    heavy = sc.moving & ~sc.is_h
    near = sc.near
    if(heavy.sum() == 0 or len(near) == 0): return np.full(len(xyz), 1.e9)
    x = xyz[:, heavy]                                            # (T, m, 3)
    d = ((x[:, :, None, :] - near[None, None, :, :]) ** 2).sum(axis=3)
    return np.sqrt(d.min(axis=(1, 2)))

  def _min_distance(self, sc, xyz):
    heavy = sc.moving & ~sc.is_h
    near = sc.near
    if(heavy.sum() == 0 or len(near) == 0): return 1.e9
    d = ((xyz[heavy][:, None, :] - near[None, :, :]) ** 2).sum(axis=2)
    return math.sqrt(d.min())

  # --- search ---------------------------------------------------------------
  def _fit_search_scale(self):
    """ Scale of the search maps: least squares against the exact occupancy
    gradients of the current moving atoms. """
    i_list, x_list = [], []
    for rec in self.records:
      sc = rec.sc
      i_list.append(sc.i_seqs[sc.moving]); x_list.append(sc.xyz[sc.moving])
    if(len(i_list) == 0):
      self.search_scale = 1.; self.search_cc = None; return
    exact = np.concatenate(self.ev.occupancy_gradients(i_list, x_list))
    approx = np.concatenate([self.ev.search_linear(x,
      self.ev.search_setup(i), 1.) for i, x in zip(i_list, x_list)])
    self.search_scale = float((exact * approx).sum() / (approx * approx).sum())
    self.search_cc = float(np.corrcoef(exact, approx)[0, 1])

  def _search_score(self, rec, xyz):
    sc = rec.sc
    setup, ff_old, ff_new, lin_old = rec._search
    new = xyz[sc.moving]
    if(self._min_distance(sc, xyz) < self.params.clash_distance):
      return -1.e30
    lin = float((rec.occ_new * self.ev.search_linear(new, setup,
      self.search_scale) - rec.occ_old * lin_old).sum())
    quad = self.ev.quadratic(np.concatenate([sc.xyz[sc.moving], new]),
      np.concatenate([-ff_old, ff_new]))
    return -self.ev.n_work * (lin + quad)

  def _search_scores(self, rec, xyzs):
    """ _search_score for T sets of residue sites (T, n_atoms, 3) """
    sc = rec.sc
    setup, ff_old, ff_new, lin_old = rec._search
    new = xyzs[:, sc.moving]                                     # (T, m, 3)
    t, m = new.shape[0], new.shape[1]
    lin_new = self.ev.search_linear(new.reshape(-1, 3), (np.tile(setup[0], t),
      np.tile(setup[1], t), np.tile(setup[2], t)), self.search_scale)
    lin = (rec.occ_new[None, :] * lin_new.reshape(t, m)).sum(axis=1) - \
      float((rec.occ_old * lin_old).sum())
    old = np.broadcast_to(sc.xyz[sc.moving], (t, m, 3))
    quad = self.ev.quadratic_batch(np.concatenate([old, new], axis=1),
      np.concatenate([-ff_old, ff_new]))
    score = -self.ev.n_work * (lin + quad)
    clash = self._min_distances(sc, xyzs) < self.params.clash_distance
    score[clash] = -1.e30
    return score

  def _tune(self, rec):
    p = self.params
    sc = rec.sc
    i_seqs = sc.i_seqs[sc.moving]
    setup = self.ev.search_setup(i_seqs)
    types, b, occ = self.ev.atom_params(i_seqs)
    ff_old = self.ev.form_factors(types, b, rec.occ_old)
    ff_new = self.ev.form_factors(types, b, rec.occ_new)
    lin_old = self.ev.search_linear(sc.xyz[sc.moving], setup,
      self.search_scale)
    rec._search = (setup, ff_old, ff_new, lin_old)
    n_steps = int(round(p.tune_range / p.tune_step)) if p.tune_step > 0 else 0
    offsets = [k * p.tune_step for k in range(-n_steps, n_steps + 1)]
    for cand in rec.candidates:
      centre = list(cand.chis)
      best_chis = list(centre)
      best_xyz = cand.xyz
      best = self._search_score(rec, best_xyz)
      for _ in range(p.tune_passes if n_steps > 0 else 0):
        for k in range(len(centre)):
          base_xyz, base_chi = best_xyz, best_chis[k]
          targets = [centre[k] + off for off in offsets
            if centre[k] + off != base_chi]
          if(not targets): continue
          trials = np.array([sc.rotate_chi(base_xyz, k, x - base_chi)
            for x in targets])
          scores = self._search_scores(rec, trials)
          i = int(np.argmax(scores))
          if(scores[i] > best):
            best, best_xyz = float(scores[i]), trials[i]
            best_chis = list(best_chis); best_chis[k] = targets[i]
      cand.chis = best_chis
      cand.xyz = best_xyz
      cand.search = best
      cand.clash = best <= -1.e29

  # --- prediction and exact evidence ------------------------------------------
  def _predict(self):
    i_list, x_list = [], []
    for rec in self.records:
      sc = rec.sc
      i_list.append(sc.i_seqs[sc.moving]); x_list.append(sc.xyz[sc.moving])
      for cand in rec.candidates:
        if(cand.clash): continue
        i_list.append(sc.i_seqs[sc.moving])
        x_list.append(cand.xyz[sc.moving])
    grads = self.ev.occupancy_gradients(i_list, x_list)
    if(self.params.mask_correction): self._predict_mask()
    k = 0
    for rec in self.records:
      sc = rec.sc
      i_seqs = sc.i_seqs[sc.moving]
      lin_old = grads[k]; k += 1
      ff = rec._search[1]
      ff_new = rec._search[2]
      # evidence for removing the moving atoms (current conformer wrong):
      # dLL ~ -I/2 if they are right
      occ = rec.occ_old
      q_rm = self.ev.quadratic(sc.xyz[sc.moving], -ff)
      rec.removal = -self.ev.n_work * (-float((occ * lin_old).sum()) + q_rm)
      rec.removal_information = 2 * self.ev.n_work * q_rm
      # per atom (self term only): dLL for removing each moving atom alone,
      # positive where the atom sits in no density
      q_atom = 0.5 * self.ev.n_sym * (ff ** 2).dot(self.ev.shell_c)
      rec.atom_removal = -self.ev.n_work * (-occ * lin_old + q_atom)
      rec.atom_removal_information = 2 * self.ev.n_work * q_atom
      rec.removal_exact = None
      if(self.params.removal_exact):
        rec.removal_exact = self.ev.exact_removal(i_seqs, sc.xyz[sc.moving])
      for cand in rec.candidates:
        if(cand.clash): continue
        lin_new = grads[k]; k += 1
        lin = float((rec.occ_new * lin_new - rec.occ_old * lin_old).sum())
        xyz_q = [sc.xyz[sc.moving], cand.xyz[sc.moving]]
        ff_q = [-ff, ff_new]
        m = getattr(cand, "_mask", None)
        if(m is not None and len(m[1])):
          # the solvent displaced or restored by the move
          lin += cand._mask_lin
          xyz_q.append(m[0]); ff_q.append(self.ev.mask_form_factors(m[1]))
        quad = self.ev.quadratic(np.concatenate(xyz_q), np.concatenate(ff_q))
        cand.predicted = -self.ev.n_work * (lin + quad)
        # information (Fisher, from the curvature) for telling a from b:
        # dLL ~ N(-I/2, I) if a is right, N(+I/2, I) if b is
        cand.information = 2 * self.ev.n_work * quad
        cand.information_size = 2 * self.llg_total * \
          rec.moving_scattering / self.model_scattering

  def _predict_mask(self):
    """ The change of the bulk-solvent mask for every candidate (local grid,
    once per residue for the atoms that stay): pseudo-atoms of solvent and
    their linear term, from one FFT gradient calculation. """
    ev = self.ev
    pos_list, cands = [], []
    for rec in self.records:
      sc = rec.sc; mv = sc.moving; ii = sc.i_seqs[mv]
      live = [c for c in rec.candidates if not c.clash]
      if(not live): continue
      ctx = ev.mask_context(ii, [sc.xyz[mv]] + [c.xyz[mv] for c in live],
        rec.occ_old, sc.i_seqs, spacing=self.params.mask_predict_spacing)
      if(ctx is None): continue
      s_old = ev.mask_solvent(ctx, sc.xyz[mv], rec.occ_old)
      for c in live:
        c._mask = ev.mask_cells(ctx, s_old, ev.mask_solvent(ctx, c.xyz[mv],
          rec.occ_new))
        pos_list.append(c._mask[0]); cands.append(c)
    for c, g in zip(cands, ev.mask_gradients(pos_list)):
      c._mask_lin = float((c._mask[1] * g).sum())

  def _exact_and_decide(self):
    p = self.params
    others = [c for rec in self.records for c in rec.candidates
      if not c.clash]
    rng = np.random.RandomState(p.random_seed)
    n = min(p.calibration_sample, len(others))
    self._sample = set(id(others[i]) for i in
      rng.choice(len(others), n, replace=False)) if n > 0 else set()
    for rec in self.records:
      sc = rec.sc
      i_seqs = sc.i_seqs[sc.moving]
      live = [c for c in rec.candidates if not c.clash]
      if(len(live) == 0): continue
      if(rec.stub):
        self._decide_stub(rec, live, i_seqs); continue
      lp_a = math.log(rec.current_prior)
      for c in live:
        c.log_odds = c.predicted + math.log(c.prior) - lp_a
      f_old = None
      for c in live:
        if((p.exact and c.log_odds > -p.screen) or id(c) in self._sample):
          if(f_old is None): f_old = self.ev.f_calc_of(i_seqs,
            sc.xyz[sc.moving])
          c.exact = self.ev.exact(i_seqs, sc.xyz[sc.moving],
            c.xyz[sc.moving], f_old=f_old, mask_exclude=sc.i_seqs
            if self.params.mask_correction else None)
      self._set_offsets(live)
      for c in live:
        c.log_odds = c.dll(p.exact) + math.log(c.prior) - lp_a
      rec.best = max(live, key=lambda c: c.log_odds)
      rec.accept = rec.best.log_odds > p.threshold and self._gates(rec.best)
      lo = np.array([0.] + [c.log_odds for c in live])
      w = np.exp(lo - lo.max())
      rec.posterior_current = float(w[0] / w.sum())

  def _set_offsets(self, live):
    """ Put the candidates without an exact dLL on the scale of those with
    one (median exact - predicted in the residue), so that they are not
    compared on the prediction's scale with exact values """
    d = [c.exact - c.predicted for c in live if c.exact is not None]
    off = float(np.median(d)) if d else 0.
    for c in live: c.offset = off

  def _decide_stub(self, rec, live, i_seqs):
    """ Posterior over wells for building a stub: P(b) ~ prior_b exp(dLL_add)
    (exact dLL for the wells within screen of the best prediction) """
    p = self.params
    sc = rec.sc
    lp = np.array([c.predicted + math.log(c.prior) for c in live])
    f_old = None
    for c, x in zip(live, lp):
      if((p.exact and x > lp.max() - p.screen) or id(c) in self._sample):
        if(f_old is None):
          f_old = flex.complex_double(self.ev.fmodel.f_obs().size(), 0)
        c.exact = self.ev.exact(i_seqs, sc.xyz[sc.moving], c.xyz[sc.moving],
          f_old=f_old, occ_new=rec.occ_new, occ_old=rec.occ_old,
          mask_exclude=sc.i_seqs if p.mask_correction else None)
    self._set_offsets(live)
    lp = np.array([c.dll(p.exact) + math.log(c.prior) for c in live])
    w = np.exp(lp - lp.max()); w /= w.sum()
    for c, x in zip(live, w):
      c.posterior = float(x)
      c.log_odds = math.log(max(x, 1.e-300) / max(1 - x, 1.e-300))
    rec.best = live[int(np.argmax(w))]
    rec.posterior_best = float(w.max())
    rec.accept = (rec.posterior_best > p.stub_threshold and
      self._gates(rec.best))

  def _ncs_decide(self, pdb_hierarchy, ncs_groups):
    """ decisions with the evidence of NCS-related copies pooled """
    p = self.params
    keys = ncs_groups if ncs_groups is not None else \
      ncs_module.groups_from_search(pdb_hierarchy)
    by_key = dict((r.key, r) for r in self.records if not r.stub)
    groups = []
    for g in keys:
      recs = [by_key[tuple(k)] for k in g if tuple(k) in by_key]
      if(len(recs) < 2): continue
      names = [r.sc.resname for r in recs]
      main = max(set(names), key=names.count)
      recs = [r for r in recs if r.sc.resname == main]
      if(len(recs) >= 2): groups.append(recs)
    if(not groups): return
    pool = ncs_module.pooled(groups, self.priors, p.exact)
    self.ncs_rho = pool.rho
    self.ncs_n_groups = len(pool.groups)
    scale = self._information_scale()
    for rec, lp, r_shared, others, e in pool.posteriors():
      rec.ncs_shared = r_shared
      cur = rec.current_well if rec.current_well in lp else "CURRENT"
      best = max(lp, key=lp.get)
      rec.accept = False
      if(best == cur): continue
      c = e.cands.get(best)
      if(c is None): continue
      lo = lp[best] - lp[cur]
      # evidence and information expected under the group's posterior
      d = c.dll(p.exact)
      info = c.information or 0.
      if(cur != "CURRENT"):
        for x in others:
          db = x.log_l.get(best, ncs_module.NEG)
          dc = x.log_l.get(cur, ncs_module.NEG)
          if(db > ncs_module.NEG / 2 and dc > ncs_module.NEG / 2):
            d += r_shared * (db - dc)
          xc = x.cands.get(best) or x.cands.get(cur)
          if(xc is not None and xc.information):
            info += r_shared * xc.information
      c.log_odds = lo
      rec.best = c
      rec.ncs_dll = d
      # the copy's own data must favour the change too: pooling adds
      # confidence, it does not override a copy (copies can differ, and
      # errors shared by the copies add up as if independent)
      rec.accept = (lo > p.threshold and d > p.min_dll and
        c.dll(p.exact) > p.min_dll and scale * info >= p.min_information)

  # --- output -----------------------------------------------------------------
  def accepted(self):
    """ Accepted records, without conflicts: of two accepted changes whose
    moved atoms come within 4 A of each other, only the one with the larger
    log odds is kept (the other is tested again at the next call). """
    acc = sorted([r for r in self.records if r.accept],
      key=lambda r: -r.best.log_odds)
    kept = []
    for r in acc:
      x = np.concatenate([r.best.xyz[r.sc.moving], r.sc.xyz[r.sc.moving]])
      ok = True
      for q in kept:
        y = np.concatenate([q.best.xyz[q.sc.moving], q.sc.xyz[q.sc.moving]])
        if(((x[:, None, :] - y[None, :, :]) ** 2).sum(axis=2).min() < 16.):
          ok = False; break
      if(ok): kept.append(r)
    return kept

  def apply(self, pdb_hierarchy, xray_structure=None):
    """ Make the accepted changes; returns their number. The caller updates
    any f_model (f_calc and mask) afterwards. """
    atoms = pdb_hierarchy.atoms()
    sites = atoms.extract_xyz()
    kept = self.accepted()
    for r in kept:
      for i, x in zip(r.sc.i_seqs, r.best.xyz):
        sites[int(i)] = tuple(x)
      if(r.stub):   # built: the moving atoms get their occupancy
        for i, o in zip(r.sc.i_seqs[r.sc.moving], r.occ_new):
          atoms[int(i)].set_occ(float(o))
        if(xray_structure is not None):
          occ = xray_structure.scatterers().extract_occupancies()
          for i, o in zip(r.sc.i_seqs[r.sc.moving], r.occ_new):
            occ[int(i)] = float(o)
          xray_structure.set_occupancies(occ)
    atoms.set_xyz(sites)
    if(xray_structure is not None):
      xray_structure.set_sites_cart(sites)
    return len(kept)

  def show(self, out=None, all_residues=False):
    if(out is None): out = self.log
    print("rotamer_bayes: %d residues tested" % len(self.records), file=out)
    for k, v in sorted(self.skipped.items()):
      print("  skipped (%s): %d" % (k, v), file=out)
    print("  LLG of the model %.1f; scattering power %.0f" % (
      self.llg_total, self.model_scattering), file=out)
    if(self.alpha is not None):
      print("  evidence self-check (%s target): alpha %.2f beta %.2f "
        "(expected ~0.5, <= %.1f): %s%s" % (self.target_name, self.alpha,
        self.beta, self.beta_max, "passed" if
        self.self_check_passed else "FAILED", "" if self.self_check_passed
        or not self.params.self_check else " -- no changes made"), file=out)
    else:
      print("  evidence self-check: too few candidates", file=out)
    if(self.ncs_rho is not None):
      print("  NCS pooling: %d groups, rho (copies share a well) %.3f" % (
        self.ncs_n_groups, self.ncs_rho), file=out)
    print("  search map scale %.4g (CC with exact gradients %s)" % (
      self.search_scale, "%.3f" % self.search_cc
      if self.search_cc is not None else "n/a"), file=out)
    print("  times: setup %.1f s, search %.1f s, predict %.1f s, exact %.1f s"
      % (self.time_setup, self.time_search, self.time_predict,
      self.time_exact), file=out)
    print("  %-16s %-8s %7s -> %-8s %7s %9s %9s %9s %s" % ("residue",
      "current", "prior", "best", "prior", "dLL_pred", "dLL_exact",
      "ln_odds", ""), file=out)
    n_acc = 0
    for r in self.records:
      if(r.best is None): continue
      if(not (all_residues or r.accept)): continue
      b = r.best
      print("  %-16s %-8s %7.4f -> %-8s %7.4f %9.2f %9s %9.2f %s" % (
        r.label, r.current_well, r.current_prior, b.well, b.prior,
        b.predicted, "%.2f" % b.exact if b.exact is not None else "-",
        b.log_odds, "ACCEPT" if r.accept else ""), file=out)
    print("  accepted: %d (%d without conflicts)" % (
      sum([1 for r in self.records if r.accept]), len(self.accepted())),
      file=out)
