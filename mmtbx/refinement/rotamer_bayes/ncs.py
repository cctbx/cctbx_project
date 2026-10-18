"""
Pooling the evidence of NCS-related copies of a residue.

For a group of NCS-equivalent residues (copies j), with L_j(w) = exp dLL_j(w)
the evidence for well w in copy j (relative to copy j's current conformer,
so L_j(current) = 1), the prior over the copies' wells is a mixture: with
probability rho all copies share one well w ~ pi (the library prior), with
probability 1 - rho the copies' wells are independent, each ~ pi. Then

  P(w_i = w | data) ~ rho pi_w prod_j L_j(w)
                    + (1 - rho) pi_w L_i(w) prod_{j != i} Z_j,
  Z_j = sum_v pi_v L_j(v),

and rho is estimated by maximising the marginal likelihood of all groups,
sum_g ln[rho S_g + (1 - rho) prod_j Z_j], S_g = sum_w pi_w prod_j L_j(w)
(empirical Bayes: if the copies differ, rho is small and pooling switches
itself off). A copy whose current conformer is not in a named well
(outlier, unclassified) keeps that conformer as a hypothesis of its own in
the independent part.

The gates of the decision use the evidence expected under the posterior of
the group: the copy's own dLL and information, plus those of the other
copies weighted by the posterior probability r that the group shares a
conformer; and the copy's own dLL must favour the change as well (in tests,
without it copies that really differ, and near-isosteric flips whose small
systematic errors are shared by the copies, were changed wrongly).
"""
from __future__ import absolute_import, division, print_function
import math
import numpy as np

NEG = -1.e300

def logsumexp(a):
  a = np.asarray(a, dtype=float)
  m = a.max() if a.size else NEG
  if(m <= NEG / 2): return NEG
  return float(m + math.log(np.exp(a - m).sum()))

def residue_key(rg, chain_id):
  return (chain_id, rg.resseq.strip(), rg.icode.strip())

def groups_from_search(pdb_hierarchy):
  """ NCS-equivalent residues: lists of (chain, resseq, icode), from the
  NCS search used by phenix.refine (CA atoms of master and copies) """
  import iotbx.ncs
  from libtbx.utils import null_out
  key_of = {}
  for chain in pdb_hierarchy.only_model().chains():
    for rg in chain.residue_groups():
      for a in rg.atoms():
        key_of[a.i_seq] = (residue_key(rg, chain.id), a.name.strip())
  try:
    ncs = iotbx.ncs.input(hierarchy=pdb_hierarchy,
      log=null_out()).get_ncs_restraints_group_list()
  except Exception:
    return []
  out = []
  for g in ncs:
    copies = [list(c.iselection) for c in g.copies]
    for k, i in enumerate(g.master_iselection):
      kr, name = key_of.get(int(i), (None, None))
      if(name != "CA"): continue
      grp = [kr]
      for c in copies:
        kc, nc = key_of.get(int(c[k]), (None, None))
        if(nc == "CA"): grp.append(kc)
      if(len(grp) > 1): out.append(grp)
  return out

class copy_evidence(object):
  """ log evidence of one copy for each named well (0 for its current well,
  -inf where there is no candidate), and its own unnamed current
  hypothesis if it has one """
  def __init__(self, rec, wells, exact):
    self.rec = rec
    self.cands = dict((c.well, c) for c in rec.candidates
      if not c.clash and c.predicted is not None)
    self.log_l = {}
    for w in wells:
      if(w == rec.current_well): self.log_l[w] = 0.
      elif(w in self.cands): self.log_l[w] = self.cands[w].dll(exact)
      else: self.log_l[w] = NEG
    self.own = None if rec.current_well in wells else \
      math.log(max(rec.current_prior, 1.e-300))
    self.log_z = None

  def logz(self, log_pi):
    terms = [log_pi[w] + self.log_l[w] for w in log_pi]
    if(self.own is not None): terms.append(self.own)
    return logsumexp(terms)

class pooled(object):
  """ rho, and the per-copy posteriors over wells for the NCS groups """

  def __init__(self, groups, priors, exact):
    self.groups = []
    for recs in groups:
      wells = priors.wells(recs[0].sc.resname)
      named = dict((w.name, w.mass) for w in wells.values() if w.name not in
        ("OUTLIER", "UNCLASSIFIED"))
      # the wells tried in every copy (same residue type, same p_min)
      tried = set(named)
      for r in recs:
        tried &= set([c.well for c in r.candidates] + [r.current_well])
      if(not tried): continue
      log_pi = dict((w, math.log(named[w])) for w in tried)
      ev = [copy_evidence(r, tried, exact) for r in recs]
      for e in ev: e.log_z = e.logz(log_pi)
      log_s = logsumexp([log_pi[w] + sum(e.log_l[w] for e in ev)
        for w in tried])
      log_p = sum(e.log_z for e in ev)
      self.groups.append((recs, log_pi, ev, log_s, log_p))
    self.rho = self._estimate_rho()

  def _estimate_rho(self):
    if(not self.groups): return 0.
    s = np.array([g[3] for g in self.groups])
    p = np.array([g[4] for g in self.groups])
    best = (None, NEG)
    for rho in np.linspace(0.001, 0.999, 999):
      ll = np.logaddexp(math.log(rho) + s, math.log(1 - rho) + p).sum()
      if(ll > best[1]): best = (float(rho), ll)
    return best[0]

  def posteriors(self):
    """ for every copy: (record, {hypothesis: log posterior (unnormalised)},
    r = posterior probability that its group shares a conformer, the other
    copies' evidence) """
    rho = self.rho
    out = []
    for recs, log_pi, ev, log_s, log_p in self.groups:
      a = math.log(rho) + log_s
      b = math.log(1 - rho) + log_p
      r_shared = 1. / (1. + math.exp(min(700., b - a)))
      for i, e in enumerate(ev):
        rest = sum(x.log_z for k, x in enumerate(ev) if k != i)
        lp = {}
        for w in log_pi:
          shared = math.log(rho) + log_pi[w] + sum(x.log_l[w] for x in ev)
          indep = math.log(1 - rho) + log_pi[w] + e.log_l[w] + rest
          lp[w] = float(np.logaddexp(shared, indep))
        if(e.own is not None):
          lp["CURRENT"] = math.log(1 - rho) + e.own + rest
        out.append((e.rec, lp, r_shared, [x for k, x in enumerate(ev)
          if k != i], e))
    return out
