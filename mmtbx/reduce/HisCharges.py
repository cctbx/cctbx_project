"""His charges across conformers: decision, rerun cost, trial planning and the
written-state check, as pure functions over plain records (no Optimizer object).

A His whose conformers end the Optimizer passes with different ring protonation
(charge +1 with HD1 and HE2, 0 with one of them) is a qmready error unless there
is a reason for it (a metal, a protonation set by the user, an ion at the
flipped ring). The caller (the
Optimizer) records each pass's His states and scores, calls decide()
per His, reruns the passes with the reconciled states fixed, and accepts or
releases them with trial() and plan(). check_written() checks the result.

Policy:
- Cost of a charge: per conformer, the chosen state's score minus the best
  state of that charge in the same (written) orientation, clamped at 0, summed.
  A flipped state is never used as the cost; a suggested flip is reported only.
- BUDGET (0.5) is a fixed provisional constant, not tied to non_flip_preference.
- No evidence (no conformer loses more than the both-H preference under either
  charge, or equal totals): neutral (charge 0) by policy. This is a stated
  tie-break, not a claim about protonation energetics.
- Left unchanged: partially split rings (no pass decides a whole conformer) and
  unscorable conformers (no Mover for a geometric reason). Kept with the reason:
  protonation set by the user, ion lock, metal-bound N, bond outside the ring,
  ion at the flipped ring.
- After the rerun the fresh scores decide, not the first-stage cost (a screen).
  Fixed His whose neighbourhoods overlap form one component with one BUDGET.
  Changes outside every neighbourhood count only as losses, never as credit.
  A His whose conformers matched in the free run and differ after the rerun (a
  new mismatch) fails its component, or the run if it lies outside them all.
- A His is kept only if it passes together with the His already accepted.

- Reruns are not assumed reproducible (the Optimizer's search depends on its
  deadline): a failure is blamed only after a fresh confirmation run, and an
  unpinned rerun that differs from the free run is never accepted silently.

Contract for the one-state Mover.
- A fixed His is built with one coarse state, so every index into the Mover's
  arrays (positions, extras, deleteMes, preference energies, fix-up positions)
  is 0 and stays valid for FixUp. The original state is recorded separately on
  the Mover: its index in the full table, whether it is flipped, which ring H it
  deletes, and its coarse and fix-up positions.
- The pose text (PoseDescription, _printPose) and the flip report come from the
  recorded original state, not from index 0 (which would always read
  "Unflipped HD1Placed HE2Placed"); consumers parsing the pose words get the
  original state's words. check_fixed() takes the recorded original index.
- check_written()'s expected positions come from the free run's full Mover at
  the original index, never from the one-state Mover (that would compare it
  with itself). fixup is always filled; an empty FixUp().positions means the
  coarse positions are retained, which is checked, not skipped.
- Mover ids must be the same in every run: the one-state Mover is built over
  the full Mover's atom list, and every rerun starts from the free run's
  pre-deletion model (same atoms, same i_seqs). Otherwise every id is unmatched
  and every trial fails (safe, but the feature does nothing).
- A His of the free run absent from a rerun (new_mismatches' missing; a His
  without a Mover has no id in unmatched()) is passed to trial() as missing.
- A separate test of search-time atom information (extras of the fixed state
  during the search, not only after FixUp) is needed: the written state can
  be right while the search scored wrong flags.
- Contract failures raise ContractError (explicit checks). The Optimizer
  handles them, and any other error after the free result is saved
  (MemoryError and a user's Stop, Abort, apart), by keeping the free result
  (restored if a rerun started) and warning; in strict mode (tests) it
  restores first, then re-raises. A MemoryError propagates as it is; an Abort
  after the restore. If the restore cannot be verified, it stops with a Sorry
  naming the error that led to it. Errors in the free run propagate as
  without reconciliation.
- Limits: at most run_cap(n) = min(max_runs(n), RUN_CAP_LIMIT) reruns per model
  (confirmations and retries count). A free run that reached its search
  deadline gets no rerun; each rerun may search as long as the free run was
  allowed to (the time left when it started, not its duration), and one that
  reaches that limit is not comparable. Reruns beyond WALL_BACKSTOP_FACTOR x the
  free run's wall time stop (backstop only). If another rerun is needed after a
  limit, the free result is restored. The rerun time is added to later models'
  deadline.
- The model outcome is a code (OUTCOME_TEXT); the report text derives from it.

Known limits: a state is scored against the His's own symmetry mates in the
chosen state (the mates follow only after the search), as in the suggested
flips.
"""
from __future__ import absolute_import, division, print_function

BUDGET = 0.5        # cost after the rerun allowed per component
NOISE = 0.005       # score change below this is not a change
MAX_ATTEMPTS = 4    # trial runs per His: joint, split, single, after a confirmation
XYZ_TOL = 0.01      # A, written coordinate against the fixed state
RUN_CAP_LIMIT = 20       # reruns per model; provisional, from validation
WALL_BACKSTOP_FACTOR = 50  # reruns' wall time, x the free run's; backstop only

# model outcome codes (Optimizer.getHisCharges()['code'])
NO_RERUN = 'no rerun'                     # nothing to rerun
RECONCILED = 'reconciled'                 # DONE: the last rerun is written
REJECTED = 'scientific rejection'         # plan() restored: no trial passed
RUN_CAP = 'run cap'                       # another rerun needed after run_cap(n)
SEARCH_DEADLINE = 'search deadline'       # the free run or a rerun reached its search time
WALL_BACKSTOP = 'wall-time backstop'      # reruns past WALL_BACKSTOP_FACTOR x the free run
INTERNAL = 'internal check'               # a ContractError or other error; the free result kept
OUTCOME_TEXT = {
  NO_RERUN: 'no rerun',
  RECONCILED: 'reconciled',
  REJECTED: 'first result restored: no fixed state passed',
  RUN_CAP: 'first result restored: rerun cap reached',
  SEARCH_DEADLINE: 'first result restored: search time limit reached',
  WALL_BACKSTOP: 'first result restored: wall-time backstop reached',
  INTERNAL: 'first result kept: internal error',
}

def run_cap(n):
  """Reruns allowed for n candidates."""
  return min(max_runs(n), RUN_CAP_LIMIT)

class ContractError(AssertionError):
  """A reconciliation check failed. Raised explicitly (kept under python -O);
  the Optimizer restores the free result on it."""

def _require(ok, what):
  if not ok:
    raise ContractError(what)

USER_PROTONATION = 'protonation set by the user'   # no Mover: protonation { } entry
METAL_REASONS = (USER_PROTONATION, 'bond outside the ring')
METAL_PREFIXES = ('ion lock', 'metal-bound')

def key_text(k):
  """'A 107' for the His key (model, chain, resseq, icode)."""
  return ' '.join(str(x) for x in k[1:]).strip()

def charge_text(q):
  """A His charge for reports: '+1', '0'."""
  return '%+d' % q if q else '0'

def decide(recs, budget=BUDGET, both_pref=0.05):
  """Decision for one His whose conformers end with different charges.
  recs: {conformer letter: record of the pass that decides its ring H}. A record
    has partial (bool), mover (None if no Mover; then noMover gives the reason),
    suggested (None or dict(gain, ion)), chosen (state index) and states
    [dict(index, score, charge, flipped)].
  both_pref: the both-H preference, 0.05 x the Optimizer's preferenceMagnitude
    (the Optimizer passes it; the default is for magnitude 1).
  Returns dict(kind='reconciled'|'kept'|'unchanged', decision=text, and, when
  costed, target, costs {'0'|'1': dict(total, per)}, per {conf: dict(cost, raw,
  to)}, policy)."""
  partial = sorted(c for c, r in recs.items() if r['partial'])
  if partial:
    return dict(kind='unchanged', decision='unchanged: partially split ring in %s (blank and split ring atoms; '
      'no pass decides a whole conformer)' % ','.join(c or 'blank' for c in partial))
  noMover = dict((c, r['noMover']) for c, r in recs.items() if r.get('mover') is None)
  reasons = ["%s %s" % (c or 'blank', w) for c, w in sorted(noMover.items())
             if w in METAL_REASONS or w.startswith(METAL_PREFIXES)]
  if reasons:
    return dict(kind='kept', decision='kept: ' + ', '.join(reasons))
  if noMover:
    return dict(kind='unchanged', decision='unchanged: unscorable conformer ' + ', '.join(
      "%s (%s)" % (c or 'blank', w) for c, w in sorted(noMover.items())))
  ion = sorted(c for c, r in recs.items() if r.get('suggested') and r['suggested']['ion'])
  if ion:
    return dict(kind='kept', decision='kept: an ion binds the flipped ring of ' + ','.join(c or 'blank' for c in ion))
  costs = {}
  for q in (0, 1):
    per = {}
    tot = 0.0
    for c, r in sorted(recs.items()):
      st = r['states']
      ch = st[r['chosen']]
      same = [s for s in st if s['flipped'] == ch['flipped'] and s['charge'] == q]
      if not same:
        per = None
        break
      best = ch if ch['charge'] == q else max(same, key=lambda s: s['score'])
      raw = ch['score'] - best['score']
      per[c] = dict(cost=round(max(0.0, raw), 3), raw=round(raw, 3), to=best['index'])
      tot += max(0.0, raw)
    costs[str(q)] = dict(total=round(tot, 3) if per is not None else None, per=per)
  valid = [q for q in (0, 1) if costs[str(q)]['total'] is not None]
  if not valid:
    return dict(kind='unchanged', decision='unchanged: no state of either charge', costs=costs)
  tot = dict((q, costs[str(q)]['total']) for q in valid)
  noEvidence = all(p['cost'] <= both_pref + 1e-3 for q in valid for p in costs[str(q)]['per'].values())
  policy = None
  if noEvidence or (len(valid) == 2 and abs(tot[0] - tot[1]) <= 1e-3):
    q = 0 if 0 in valid else valid[0]
    policy = 'tie without evidence: neutral by policy'
  else:
    q = min(valid, key=lambda x: tot[x])
  other = ''.join('; to %s: %.2f' % (charge_text(x), tot[x]) for x in valid if x != q)
  if tot[q] > budget:
    return dict(kind='kept', decision='kept: cost to %s is %.2f > %.2f%s' % (charge_text(q), tot[q], budget, other),
      target=q, costs=costs, per=costs[str(q)]['per'], policy=policy)
  return dict(kind='reconciled', decision='reconciled to charge %s (cost %.2f%s%s)' % (charge_text(q), tot[q], other,
    ('; ' + policy) if policy else ''), target=q, costs=costs, per=costs[str(q)]['per'], policy=policy)

def _owners(runs):
  # a Mover borrowed by several passes belongs to the last one (the pass that writes it)
  owner = {}
  for alt, sc in runs:
    for mid in sc:
      owner[mid] = alt
  return owner

def unmatched(before, after, ids=None):
  """Mover ids whose comparison between the free run and the rerun is not valid:
  present in only one run, or written by a different pass. ids: restrict to these
  (None: all). Returns [(id, reason)], sorted."""
  ob = _owners(before)
  oa = _owners(after)
  out = []
  for mid in set(ob) | set(oa):
    if ids is not None and mid not in ids:
      continue
    if mid not in oa:
      out.append((mid, 'absent from the rerun'))
    elif mid not in ob:
      out.append((mid, 'absent from the free run'))
    elif ob[mid] != oa[mid]:
      out.append((mid, 'written by pass %s, before by %s' % (oa[mid] or 'blank', ob[mid] or 'blank')))
  return sorted(out)

def rerun_cost(before, after, hood, his=()):
  """Loss of a group of Movers after the rerun.
  before, after: [(pass letter, {Mover id: fresh final score})] in pass order, from
  the free run and the rerun. A Mover id is the tuple of its atoms' i_seqs, so
  conformers with their own atoms (a blank ring with split H) stay distinct.
  hood: the ids counted; his: the fixed His ids among them (report split).
  A borrowed Mover counts once, in the last pass that has it (the writing pass).
  Each pass's net loss is clamped at 0: a gain in one pass does not cancel a loss
  in another. his and neighbours are the unclamped sums (they need not add up to
  total). unmatched: hood ids without a valid comparison (see unmatched()); the
  cost of a group with any is not valid.
  Returns dict(total, per {pass: loss}, changed [(pass, id, old, new)], his,
  neighbours, unmatched [(id, reason)])."""
  owner = _owners(after)
  first = dict(before)
  bad = unmatched(before, after, hood)
  skip = set(m for m, why in bad)
  per = {}
  changed = []
  hisLoss = 0.0
  nbLoss = 0.0
  for alt, sc in after:
    loss = 0.0
    for mid, score in sc.items():
      if owner[mid] != alt or mid not in hood or mid in skip:
        continue
      _require(mid in first.get(alt, {}), ('unmatched() missed', alt, mid))
      old = first[alt][mid]
      if abs(old - score) > NOISE:
        changed.append((alt, mid, old, score))
      loss += old - score
      if mid in his:
        hisLoss += old - score
      else:
        nbLoss += old - score
    per[alt] = round(max(0.0, loss), 3)
  return dict(total=round(sum(per.values()), 3), per=per, changed=changed, unmatched=bad,
    his=round(hisLoss, 3), neighbours=round(nbLoss, 3))

def components(his, cliques):
  """Group fixed His that must share one budget.
  his: {His key: set of its Mover ids}; cliques: [set of Mover ids], every
  Clique of every pass of the free run and the rerun.
  A His's neighbourhood is its Movers plus every Clique holding one of them.
  His whose neighbourhoods overlap are merged, transitively (also through a
  blank Mover that sits in different Cliques in different passes).
  Returns [dict(keys (sorted), mids, hood)], sorted by keys."""
  keys = sorted(his)
  hood = {}
  for k in keys:
    h = set(his[k])
    for c in cliques:
      if c & his[k]:
        h |= c
    hood[k] = h
  parent = dict((k, k) for k in keys)
  def find(k):
    while parent[k] != k:
      k = parent[k]
    return k
  for i, a in enumerate(keys):
    for b in keys[i + 1:]:
      if hood[a] & hood[b]:
        parent[find(b)] = find(a)
  groups = {}
  for k in keys:
    groups.setdefault(find(k), []).append(k)
  out = []
  for ks in groups.values():
    out.append(dict(keys=sorted(ks), mids=set().union(*(his[k] for k in ks)),
      hood=set().union(*(hood[k] for k in ks))))
  return sorted(out, key=lambda c: c['keys'])

def trial(before, after, his, cliques, trial_keys, mismatch, budget=BUDGET, missing=()):
  """Evaluate one rerun.
  his: {His key: set of Mover ids} for every fixed His of the run (empty for the
  unpinned rerun); trial_keys: those tested in this run (the others were accepted
  before); mismatch: {His key: Mover ids} for His whose conformers had equal
  charges in the free run and differ now. A component (components()) passes when
  the fresh cost of its neighbourhood, His included, is within budget, every id
  in it has a valid comparison, and it holds no new mismatch. missing: His keys
  of the free run absent from this run (new_mismatches()); unattributed.
  Outside every neighbourhood: losses only, never credit. All losses there are
  summed (a run of small losses counts); when the sum exceeds NOISE every losing
  Mover is listed. An unmatched id or a mismatch there is unattributed.
  Returns dict(components [dict(keys, trial, cost, mismatch, ok)],
  unattributed [(pass, id, old, new)] losses, unattributedUnmatched [(id, why)],
  unattributedMismatch [keys], unattributedMissing [keys])."""
  comps = []
  hoodAll = set()
  for g in components(his, cliques):
    hoodAll |= g['hood']
    rc = rerun_cost(before, after, g['hood'], g['mids'])
    mm = sorted(k for k, ms in mismatch.items() if ms & g['hood'])
    comps.append(dict(keys=g['keys'], trial=sorted(set(g['keys']) & set(trial_keys)), cost=rc,
      mismatch=mm, ok=rc['total'] <= budget and not mm and not rc['unmatched']))
  owner = _owners(after)
  first = dict(before)
  um = [(m, why) for m, why in unmatched(before, after) if m not in hoodAll]
  skip = set(m for m, why in um)
  losses = []
  for alt, sc in after:
    for mid, score in sc.items():
      if owner[mid] == alt and mid not in hoodAll and mid not in skip and first[alt][mid] - score > 0:
        losses.append((alt, mid, first[alt][mid], score))
  if sum(o - n for a, m, o, n in losses) <= NOISE:
    losses = []
  return dict(components=comps, unattributed=losses, unattributedUnmatched=um,
    unattributedMismatch=sorted(k for k, ms in mismatch.items() if not ms & hoodAll),
    unattributedMissing=sorted(missing))

RERUN = 'rerun'                 # run the passes again with pins = accepted + trial
DONE = 'done'                   # the last run is the written result
RESTORE_FREE = 'restore free'   # write the saved free-run result
FAILED = 'failed'               # no valid result; report failure

def new_plan(candidates):
  """Planner state for the round-0 joint trial of every reconcilable His.
  accepted: kept His; trial: His tested in the next run; pending: trials to run;
  single: His already tried alone; rejected: {key: reason}; attempts: {key: trial
  runs}; held: (trial keys, evaluation) of a failed trial waiting for its
  confirmation run; runs: evaluated reruns; requeued: held trials run again after a passing
  confirmation; n: candidates; reason: why the plan ended in RESTORE_FREE or
  FAILED."""
  _require(len(set(candidates)) == len(candidates), candidates)
  st = dict(accepted=set(), trial=sorted(candidates), pending=[], single=set(), rejected={},
    attempts={}, held=None, requeued=set(), runs=0, n=len(candidates), reason=None)
  _start(st)
  return st

def max_runs(n):
  """Bound on the reruns for n candidates: at most MAX_ATTEMPTS * n trial runs,
  each followed by at most one confirmation run, plus one written run after each
  release (at most n) and the last one."""
  return 2 * MAX_ATTEMPTS * n + n + 1

def _start(st):
  for k in st['trial']:
    st['attempts'][k] = st['attempts'].get(k, 0) + 1
    # backstop: the pending loop rejects a His at the limit before it starts
    _require(st['attempts'][k] <= MAX_ATTEMPTS, (k, st['attempts'][k]))

def _why(c):
  if c['mismatch']:
    return 'new mismatch at His %s' % ', '.join(key_text(k) for k in c['mismatch'])
  if c['cost']['unmatched']:
    return 'no valid comparison for %d Mover(s): %s' % (len(c['cost']['unmatched']), c['cost']['unmatched'][0][1])
  return 'cost after rerun %.2f > budget' % c['cost']['total']

def _unattributed(ev):
  ua = []
  if ev['unattributed']:
    ua.append('loss %.2f outside every neighbourhood' % sum(o - n for a, m, o, n in ev['unattributed']))
  if ev['unattributedUnmatched']:
    ua.append('no valid comparison for %d Mover(s) outside every neighbourhood' % len(ev['unattributedUnmatched']))
  if ev['unattributedMismatch']:
    ua.append('new mismatch at His %s' % ', '.join(key_text(k) for k in ev['unattributedMismatch']))
  if ev.get('unattributedMissing'):
    ua.append('His %s absent from the rerun' % ', '.join(key_text(k) for k in ev['unattributedMissing']))
  return ua

def plan(st, ev, first_cost, free_saved=True):
  """Next step after a rerun evaluated by trial(). st: from new_plan(), updated in
  place. first_cost: {key: first-stage cost}, used only to order retries.
  free_saved: the caller kept the free run's result and can write it.
  Returns RERUN (run with pins = accepted + trial), DONE (the last run is the
  result), RESTORE_FREE or FAILED (st['reason'] says why).
  Reruns are not assumed reproducible (the search depends on the deadline), so a
  failure is blamed only on fresh evidence:
  - Everything passes: the trial is accepted.
  - A trial fails while His are accepted: nothing is blamed yet. The next run
    has the accepted set alone (confirmation), and the trial is held.
    - The confirmation passes: the held run is judged, the accepted set now
      shown to pass alone (below). If only accepted components failed in it
      and every trial component passed, the failure may be noise: the held
      trial is run again, once.
    - It fails: the failing accepted component is released (the whole
      component); an unattributed loss, unmatched Mover or mismatch releases
      every accepted His. The held trial is retried.
  - Judging a failed trial (no His accepted, or after a passing confirmation):
    a failing component's trial His are rejected (one His) or retried one at a
    time (several, unless tried alone already), lowest first-stage cost first,
    on top of those accepted. An unattributed loss, unmatched Mover or mismatch,
    or an accepted component that fails, fails a one-component trial and splits
    a trial of several components into one trial per component.
  - Nothing accepted and nothing left to try: RESTORE_FREE when free_saved,
    without another run. Else (free_saved False) a run without pins, which is
    the free state again; if it differs from the free run (unattributed loss,
    unmatched Mover, new mismatch) the result is not accepted silently: FAILED.
  Known limits: blame within one run assumes components independent; a rejection
  is not revisited when a later release removes its reason; release-all after an
  unattributed change is conservative.
  Each His enters at most MAX_ATTEMPTS trial runs; a His at the limit is rejected
  (asserted in _start). The reruns are bounded by max_runs(n) (asserted). DONE
  only when the last run's pins equal the accepted set and it passed (asserted).
  ev must be this run's evaluation: its component keys equal the run's pins
  (asserted)."""
  st['runs'] += 1
  _require(st['runs'] <= max_runs(st['n']), (st['runs'], st['n']))
  trialKeys = list(st['trial'])
  runPins = frozenset(st['accepted'] | set(trialKeys))
  evKeys = frozenset(k for c in ev['components'] for k in c['keys'])
  _require(evKeys == runPins, ('evaluation of another run', sorted(evKeys), sorted(runPins)))
  ua = _unattributed(ev)
  bad = [c for c in ev['components'] if not c['ok']]
  order = lambda ks: sorted(ks, key=lambda k: (first_cost.get(k, 0.0), k))
  def reject(k, why):
    st['rejected'][k] = why
    st['accepted'].discard(k)
  def retry(keys, why):
    for k in reversed(order(keys)):
      if len(keys) > 1 and k not in st['single']:
        st['single'].add(k)
        st['pending'].insert(0, [k])
      else:
        reject(k, why)
  def judge(trialKeys, ev):
    comps = ev['components']
    ua = _unattributed(ev)
    tcomps = [c for c in comps if c['trial']]
    acc = ['accepted His %s fails' % ', '.join(key_text(k) for k in c['keys'])
           for c in comps if not c['ok'] and not c['trial']]
    if acc and not ua and all(c['ok'] for c in tcomps) and frozenset(trialKeys) not in st['requeued']:
      # only accepted His failed, and they passed alone since: run the trial again, once
      st['requeued'].add(frozenset(trialKeys))
      st['pending'].insert(0, list(trialKeys))
    elif ua or acc:
      why = '; '.join(ua + acc)
      if len(tcomps) > 1:
        for c in reversed(sorted(tcomps, key=lambda c: order(c['trial']))):
          st['pending'].insert(0, c['trial'])
      else:
        for c in tcomps:
          retry(c['trial'], why)
    else:
      for c in tcomps:
        if c['ok']:
          st['accepted'] |= set(c['trial'])
        else:
          retry(c['trial'], _why(c))
  if not runPins:
    if ua:
      st['trial'] = []
      st['reason'] = 'the unpinned rerun differs from the free run: ' + '; '.join(ua)
      return RESTORE_FREE if free_saved else FAILED
  elif not bad and not ua:
    st['accepted'] |= set(trialKeys)
    if st['held'] is not None:
      # passing confirmation: judge the held trial against the accepted set
      held, st['held'] = st['held'], None
      judge(*held)
  elif not trialKeys:
    # failing confirmation, or a failing written run
    if ua:
      for k in sorted(st['accepted']):
        reject(k, 'released: ' + '; '.join(ua))
    for c in bad:
      for k in c['keys']:
        reject(k, 'released: ' + _why(c))
    if st['held'] is not None:
      st['pending'].insert(0, st['held'][0])
      st['held'] = None
  elif st['accepted']:
    st['held'] = (trialKeys, ev)
    st['trial'] = []
    return RERUN
  else:
    judge(trialKeys, ev)
  while st['pending']:
    nxt = [k for k in st['pending'].pop(0) if k not in st['rejected'] and k not in st['accepted']]
    over = [k for k in nxt if st['attempts'].get(k, 0) >= MAX_ATTEMPTS]
    for k in over:
      reject(k, 'attempt limit %d' % MAX_ATTEMPTS)
    nxt = [k for k in nxt if k not in over]
    if nxt:
      st['trial'] = nxt
      _start(st)
      return RERUN
  st['trial'] = []
  if not st['accepted'] and runPins and free_saved:
    st['reason'] = 'nothing accepted'
    return RESTORE_FREE
  if frozenset(st['accepted']) != runPins:
    return RERUN
  # the written run: its pins are the accepted set, and it passed as a whole
  _require(st['held'] is None and not bad and not ua, (runPins, bad, ua))
  return DONE

RING_H = (('HE2', 'NE2'), ('HD1', 'ND1'))
RING = ('CG', 'ND1', 'CD2', 'CE1', 'NE2')

def check_written(fixed, written, skip_bond_fixup=False, tol=XYZ_TOL):
  """Check a reconciled His after the deletion step against its fixed states,
  driven by the expected inventory.
  fixed: {conformer: dict(dels={'HE2': bool, 'HD1': bool} (True: deleted),
    coarse={atom name: xyz}, fixup={atom name: xyz})}, from the free run's full
    Mover at the original state index (see the one-state Mover contract). The
    expected positions are coarse with skip_bond_fixup or an empty fixup (FixUp
    retains the coarse positions), else fixup over coarse.
  written: {conformer: dict(present={'HE2': bool, 'HD1': bool},
    xyz={atom name: xyz}, acceptor={'NE2': bool, 'ND1': bool},
    writer={atom name: conformer of the pass that writes it})}. Each conformer's
    record is checked; a blank (shared) atom against the fixed state of the
    conformer whose pass writes it (so two records of one atom must agree).
  Per fixed conformer: its record exists; each ring H is present exactly when
  the state keeps it; each expected atom has its coordinates within tol (ring,
  retained H; this also fixes the orientation: a flipped twin keeps the same H
  but moves the ring); each ring N is an acceptor exactly when its H is deleted.
  Anything expected but absent is an error, and so is an incomplete inventory
  (a fixed state without the five ring atoms or a kept ring H) or an atom
  written by a pass whose conformer is not fixed (every conformer of a
  reconciled His is fixed). fixup must be a dict. Counted skips: a deleted H's
  coordinates; a written conformer that is not fixed.
  Returns dict(errors [text], checked {'H', 'xyz', 'N'}, skipped [(conf, item,
  reason)])."""
  errors = []
  checked = dict(H=0, xyz=0, N=0)
  skipped = []
  for c in sorted(fixed):
    cn = c or 'blank'
    f0 = fixed[c]
    _require(isinstance(f0['fixup'], dict), ('fixup must be a dict', cn, f0['fixup']))
    need = list(RING) + [hn for hn, nn in RING_H if not f0['dels'][hn]]
    lack = [n for n in need if n not in f0['coarse']]
    if lack:
      errors.append('His conformer %s: fixed state lacks %s' % (cn, ' '.join(lack)))
    if c not in written:
      errors.append('His conformer %s: no written record' % cn)
      continue
    wr = written[c]
    def writer(name):
      t = wr.get('writer', {}).get(name, c)
      if t not in fixed:
        errors.append('His conformer %s %s written by pass %s, not fixed' % (cn, name, t or 'blank'))
        return None
      return fixed[t]
    for hn, nn in RING_H:
      f = writer(hn)
      if f is None:
        continue
      if hn not in wr.get('present', {}):
        errors.append('His conformer %s %s: presence not recorded' % (cn, hn))
        continue
      checked['H'] += 1
      if wr['present'][hn] == f['dels'][hn]:
        errors.append('His conformer %s %s %s, fixed state %s it' % (cn, hn,
          'kept' if wr['present'][hn] else 'deleted', 'deletes' if f['dels'][hn] else 'keeps'))
    names = set(f0['coarse']) | set(f0['fixup'])
    for name in sorted(names):
      f = writer(name)
      if f is None:
        continue
      if name in f['dels'] and f['dels'][name]:
        skipped.append((c, name, 'deleted H'))
        continue
      exp = dict(f['coarse'])
      if not skip_bond_fixup:
        exp.update(f['fixup'])
      if name not in exp:
        errors.append('His conformer %s %s: no expected position' % (cn, name))
        continue
      if name not in wr.get('xyz', {}):
        errors.append('His conformer %s %s: coordinates not recorded' % (cn, name))
        continue
      checked['xyz'] += 1
      d = sum((a - b) ** 2 for a, b in zip(wr['xyz'][name], exp[name])) ** 0.5
      if d > tol:
        errors.append('His conformer %s %s written %.3f A from the fixed state' % (cn, name, d))
    for hn, nn in RING_H:
      f = writer(nn)
      if f is None:
        continue
      if nn not in wr.get('acceptor', {}):
        errors.append('His conformer %s %s: acceptor flag not recorded' % (cn, nn))
        continue
      checked['N'] += 1
      if wr['acceptor'][nn] != f['dels'][hn]:
        errors.append('His conformer %s %s acceptor %s with %s %s in the writing state' % (
          cn, nn, wr['acceptor'][nn], hn, 'deleted' if f['dels'][hn] else 'kept'))
  for c in sorted(set(written) - set(fixed)):
    skipped.append((c, 'conformer', 'not fixed'))
  return dict(errors=errors, checked=checked, skipped=skipped)

def check_fixed(key, pins, chosen, final):
  """Assert that a run kept a fixed His as fixed. pins: {conformer: fixed state
  (original index)}; chosen: {conformer: state the run selected, as the recorded
  original index of the one-state Mover}; final: {conformer: written charge},
  which must hold every pinned conformer."""
  missing = sorted(set(pins) - set(final))
  if missing:
    raise ContractError("His %s fixed but conformer(s) %s missing from the final state: %s" % (
      key_text(key), ','.join(c or 'blank' for c in missing), final))
  for c, pin in sorted(pins.items()):
    if chosen.get(c) != pin:
      raise ContractError("His %s conformer %s state %s, fixed %s" % (key_text(key), c or 'blank', chosen.get(c), pin))
  if len(set(final.values())) != 1:
    raise ContractError("His %s fixed but charges differ: %s" % (key_text(key), final))

def new_mismatches(free_final, now_final, fixed):
  """His whose conformers had equal charges in the free run and differ now.
  free_final, now_final: {key: {conformer: charge}}; fixed: keys fixed in this
  run (excluded). Returns (mismatch keys, keys of the free run missing now)."""
  mm = sorted(k for k, f0 in free_final.items() if k not in fixed and k in now_final
    and len(set(f0.values())) == 1 and len(set(now_final[k].values())) > 1)
  missing = sorted(k for k in free_final if k not in now_final)
  return mm, missing
