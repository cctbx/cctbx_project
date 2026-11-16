"""mmtbx.reduce.HisCharges: the decision, cost, planning and written-state
check for His whose conformers end with different charges, on synthetic
records (no model runs). Mover ids are tuples of atom i_seqs. Every case was
checked once against a mutation of the function it guards."""
from __future__ import absolute_import, division, print_function
import random
from mmtbx.reduce import HisCharges as HC

K1 = (0, 'A', '107', '')
K2 = (0, 'A', '200', '')
K3 = (0, 'A', '300', '')

# ---------------------------------------------------------------- decide
def rec(scores, chosen, flipped=False, partial=False, sug=None):
  # unflipped states 0 both H (+1), 1 HD1 only (0), 2 HE2 only (0); 3-5 flipped
  st = [dict(index=i, score=s, charge=1 if i % 3 == 0 else 0, flipped=i >= 3) for i, s in enumerate(scores)]
  return dict(states=st, chosen=chosen, partial=partial, suggested=sug, mover='MoverHisFlip')

# id -> (kind, target, selected state per conformer)
DECIDED = dict(
  neutral_1exr_like=('reconciled', 0, {'A': 2, 'B': 1}), both_h_1us0_like=('reconciled', 1, {'A': 0, 'B': 0}),
  tie_neutral_by_policy=('reconciled', 0, {'A': 2, 'B': 1}), kept_cost_3i39_like=('kept', 0, {'A': 1, 'B': 1}),
  suggested_flip_warning_only=('reconciled', 0, {'A': 1, 'B': 1}), ion_at_flipped_ring=('kept', None, {}),
  negative_cost_clamped=('kept', 0, {'A': 1, 'B': 1, 'C': 1}), partly_split_ring=('unchanged', None, {}),
  unscorable_failed_mover=('unchanged', None, {}), ion_lock_kept=('kept', None, {}), user_protonation_kept=('kept', None, {}),
  flipped_ring_same_orientation=('reconciled', 0, {'A': 5, 'B': 1}), edge_cost_0_49=('reconciled', 0, {'A': 1, 'B': 1}),
  edge_cost_0_51=('kept', 0, {'A': 1, 'B': 1}), both_h_preference_edge=('reconciled', 1, {'A': 0, 'B': 0}))

DECIDE = [
  ('neutral_1exr_like', '1exr-like: A +1 by 0.29, B 0 by 3.71',
    dict(A=rec([-22.24, -22.57, -22.53], 0), B=rec([-6.05, -2.34, -12.43], 1)), 'reconciled to charge 0'),
  ('both_h_1us0_like', '1us0-like: B pays only the both-H preference',
    dict(A=rec([0.0, -5.28, -6.0], 0), B=rec([-0.05, 0.0, -3.0], 1)), 'reconciled to charge +1'),
  ('tie_neutral_by_policy', 'tie, no evidence: 0.05 each way',
    dict(A=rec([0.0, -0.05, -0.04], 0), B=rec([-0.05, 0.0, -1.0], 1)), 'neutral by policy'),
  ('kept_cost_3i39_like', '3i39-like: costs above the cutoff',
    dict(A=rec([0.0, -4.08, -5.0], 0), B=rec([-10.34, 0.0, -1.0], 1, sug=dict(gain=1.15, ion=False))), 'kept: cost'),
  ('suggested_flip_warning_only', 'suggested flip is a warning only',
    dict(A=rec([0.0, -0.2, -0.3], 0, sug=dict(gain=20.0, ion=False)), B=rec([-6.0, 0.0, -2.0], 1)), 'reconciled to charge 0'),
  ('ion_at_flipped_ring', 'ion at the flipped ring: metal reason',
    dict(A=rec([0.0, -0.2, -0.3], 0, sug=dict(gain=None, ion=True)), B=rec([-6.0, 0.0, -2.0], 1)), 'kept: an ion binds'),
  ('negative_cost_clamped', 'negative cost clamped, does not cancel',
    dict(A=rec([0.0, 0.4, -1.0], 0), C=rec([0.0, -0.8, -1.0], 0), B=rec([-2.0, 0.0, -1.0], 1)), 'kept: cost to 0 is 0.80'),
  ('partly_split_ring', 'partially split ring',
    dict(A=rec([0.0, -0.1, -0.1], 0, partial=True), B=rec([-6.0, 0.0, -2.0], 1, partial=True)), 'unchanged: partially split'),
  ('unscorable_failed_mover', 'failed Mover: unscorable',
    dict(A=rec([0.0, -0.1, -0.1], 0), C=dict(mover=None, noMover='failed: CG does not have three bonded neighbors', partial=False)),
    'unchanged: unscorable'),
  ('ion_lock_kept', 'ion lock: kept',
    dict(A=rec([0.0, -0.1, -0.1], 0), B=dict(mover=None, noMover='ion lock (ZN)', partial=False)), 'kept: B ion lock'),
  ('user_protonation_kept', 'protonation set by the user: kept',
    dict(A=rec([0.0, -0.1, -0.1], 0), B=dict(mover=None, noMover=HC.USER_PROTONATION, partial=False)),
    'kept: B protonation set by the user'),
  ('flipped_ring_same_orientation', 'flipped ring, same-orientation states',
    dict(A=rec([-9, 0.5, 0.5, 0.0, -0.4, -0.2], 3), B=rec([-6.0, 0.0, -2.0], 1)), 'reconciled to charge 0 (cost 0.20'),
  ('edge_cost_0_49', 'edge: cost 0.49 reconciled',
    dict(A=rec([0.0, -0.49, -0.6], 0), B=rec([-6.0, 0.0, -2.0], 1)), 'reconciled to charge 0 (cost 0.49'),
  ('edge_cost_0_51', 'edge: cost 0.51 kept',
    dict(A=rec([0.0, -0.51, -0.6], 0), B=rec([-6.0, 0.0, -2.0], 1)), 'kept: cost to 0 is 0.51 > 0.50'),
  ('both_h_preference_edge', 'both-H preference edge: 0.06 is evidence against 0.05',
    dict(A=rec([0.0, -0.06, -0.5], 0), B=rec([-0.05, 0.0, -1.0], 1)), 'reconciled to charge +1 (cost 0.05; to 0: 0.06)'),
]

def make_decide_test(n, want, recs):
  def t():
    d = HC.decide(recs)
    got = (d['kind'], d.get('target'), dict((c, p['to']) for c, p in (d.get('per') or {}).items()))
    assert want in d['decision'] and got == DECIDED[n], (d['decision'], got)
  return t

# ---------------------------------------------------------------- rerun_cost
HA = (1, 2, 3); HB = (4, 5, 6); N1 = (11,); N2 = (12,); N3 = (13,); BL = (20, 21)
BEFORE3 = [('A', {HA: -10.0, N1: 5.0, N2: 2.0, BL: 1.0}), ('B', {HB: -3.0, N3: 4.0, BL: 1.0})]
HOOD3 = {HA, HB, N1, N3, BL}
def after3(**ch):
  names = dict(HA=HA, HB=HB, N1=N1, N2=N2, N3=N3, BL=BL)
  r = [(a, dict(d)) for a, d in BEFORE3]
  for k, v in ch.items():
    alt, lab = k.split('_', 1)
    dict(r)[alt][names[lab]] = v
  return r

def t_cost_neighbour_loss():
  rc = HC.rerun_cost(BEFORE3, after3(A_N1=-5.0), HOOD3)
  assert rc['total'] == 10.0, rc
def t_cost_gain_does_not_offset():
  rc = HC.rerun_cost(BEFORE3, after3(A_N1=5.3, B_N3=3.0), HOOD3)
  assert rc['total'] == 1.0 and rc['per']['A'] == 0.0, rc
def t_cost_blank_mover_once():
  rc = HC.rerun_cost(BEFORE3, after3(A_BL=0.0, B_BL=0.0), HOOD3)
  assert rc['total'] == 1.0 and rc['per']['A'] == 0.0, rc
def t_cost_outside_neighbourhood():
  rc = HC.rerun_cost(BEFORE3, after3(A_N2=-8.0), HOOD3)
  assert rc['total'] == 0.0, rc
def t_unattributed_loss_listed():
  # unattributed: N2 outside every neighbourhood loses 10 -> listed as a loss
  ev = HC.trial(BEFORE3, after3(A_N2=-8.0), {K1: {HA, HB}}, [{HA, N1}, {HB, N3}], [K1], {})
  assert ev['unattributed'] == [('A', N2, 2.0, -8.0)] and ev['components'][0]['ok'], ev

# ---------------------------------------------------------------- simulation helper
def simulate(cands, mids, before, world, cliques, first_cost, free_saved=True):
  """world(pins) -> (after, mismatch). Returns (state, list of run pins); the
  outcome is in st['outcome']."""
  st = HC.new_plan(cands)
  runs = []
  while True:
    pins = frozenset(st['accepted'] | set(st['trial']))
    after, mm = world(pins)
    runs.append(pins)
    assert len(runs) < 50, runs
    ev = HC.trial(before, after, dict((k, mids[k]) for k in pins), cliques, st['trial'], mm)
    out = HC.plan(st, ev, first_cost, free_saved)
    if out != HC.RERUN:
      st['outcome'] = out
      return st, runs

# ---------------------------------------------------------------- plan: rounds of one or two His
HIS1 = {K1: {HA, HB}}
NB = (30,)
BEFORE2 = [('A', {HA: -10.0, NB: 5.0}), ('B', {HB: -3.0})]
CL2 = [{HA, NB}, {HB}]
def one_his(nloss, free_saved=True):
  def world(pins):
    if K1 in pins:
      return [('A', {HA: -10.3, NB: 5.0 - nloss}), ('B', {HB: -3.0})], {}
    return [(a, dict(d)) for a, d in BEFORE2], {}
  return simulate([K1], HIS1, BEFORE2, world, CL2, {K1: 0.3}, free_saved)

def t_one_his_accepted():
  st, runs = one_his(0.0)
  assert st['accepted'] == {K1} and len(runs) == 1, (st, runs)
def t_one_his_rejected():
  # nothing accepted: the saved free result, no unpinned rerun
  st, runs = one_his(10.0)
  assert K1 in st['rejected'] and not st['accepted'] and st['outcome'] == HC.RESTORE_FREE and len(runs) == 1, (st, runs)
  assert st['reason'] == 'nothing accepted', st
  # free result not saved: one unpinned rerun, written when it matches the free run
  st, runs = one_his(10.0, free_saved=False)
  assert st['outcome'] == HC.DONE and runs == [frozenset([K1]), frozenset()], (st, runs)
def t_one_his_loss_above_budget():
  st, runs = one_his(0.5)
  assert K1 in st['rejected'] and '0.80' in st['rejected'][K1], st['rejected']
def t_one_his_larger_loss():
  st, runs = one_his(0.6)
  assert K1 in st['rejected'], st
def t_one_his_neighbour_gain():
  st, runs = one_his(-0.4)
  assert st['accepted'] == {K1}, st
def t_new_mismatch_and_missing():
  # a mismatch first seen in a rerun is reported, never attempted
  # K2 and K3 differed in the free run: not new
  mm, missing = HC.new_mismatches({K1: {'A': 0, 'B': 0}, K2: {'A': 1, 'B': 0}, K3: {'A': 1, 'B': 0}},
    {K1: {'A': 1, 'B': 0}, K2: {'A': 0, 'B': 0}, K3: {'A': 1, 'B': 0}}, set())
  assert mm == [K1] and missing == [], (mm, missing)
  # a His of the free run absent now is reported
  mm, missing = HC.new_mismatches({K1: {'A': 0, 'B': 0}}, {}, set())
  assert mm == [] and missing == [K1], (mm, missing)
def t_two_his_one_rejected():
  H2A = (40,); H2B = (41,)
  mids = {K1: {HA, HB}, K2: {H2A, H2B}}
  before = [('A', {HA: -10.0, NB: 5.0, H2A: -1.0}), ('B', {HB: -3.0, H2B: -1.1})]
  def world(pins):
    return [('A', {HA: -10.3 if K1 in pins else -10.0, NB: -5.0 if K1 in pins else 5.0,
      H2A: -1.1 if K2 in pins else -1.0}), ('B', {HB: -3.0, H2B: -1.1})], {}
  st, runs = simulate([K1, K2], mids, before, world, [{HA, NB}, {H2A}, {HB}, {H2B}], {K1: 0.3, K2: 0.1})
  assert K1 in st['rejected'] and st['accepted'] == {K2} and runs[-1] == frozenset([K2]), (st, runs)

# ---------------------------------------------------------------- check_fixed
def t_check_fixed_wrong_state():
  try:
    HC.check_fixed(K1, {'A': 2}, {'A': 0, 'B': 2}, {'A': 1, 'B': 0})
  except AssertionError as e:
    assert 'fixed 2' in str(e), str(e)
  else:
    raise AssertionError('no exception')
def t_check_fixed_charges_differ():
  try:
    HC.check_fixed(K1, {'A': 2, 'B': 0}, {'A': 2, 'B': 0}, {'A': 0, 'B': 1})
  except AssertionError as e:
    assert 'charges differ' in str(e), str(e)
  else:
    raise AssertionError('no exception')

# ---------------------------------------------------------------- Mover identity
# 1us0 His A 83: blank ring (CG 701, ND1 702, CD2 703, CE1 704, NE2 705 in the file, i_seq =
# serial - 1), split CB (709 A, 713 B); MoverHisFlip order NE2 HE2 CE1 HE1 ND1 HD1 CD2 HD2
# CG CB HB2 HB3. H i_seqs follow the heavy atoms, as placement appends them (illustrative).
def his83(conf):
  cb = 708 if conf == 'A' else 712
  h = 3400 if conf == 'A' else 3410
  return (704, h, 703, h + 1, 701, h + 2, 702, h + 3, 700, cb, h + 4, h + 5)
M83A = his83('A'); M83B = his83('B')
# 1us0 His A 306: CB and ring blank (CB 2891, CG 2892, ND1 2893, CD2 2894, CE1 2895, NE2 2896);
# only N/CA split, so the A and B ids differ in the split ring H alone (illustrative i_seqs).
def his306(conf):
  h = 3500 if conf == 'A' else 3510
  return (2895, h, 2894, 3520, 2892, h + 2, 2893, 3521, 2891, 2890, 3522, 3523)
M306A = his306('A'); M306B = his306('B')
def t_mover_id_per_conformer():
  assert M83A != M83B and M83A[0] == M83B[0] and M83A[2] == M83B[2], (M83A, M83B)
  assert M306A != M306B and len(set(M306A) ^ set(M306B)) == 4, (M306A, M306B)
  rc = HC.rerun_cost([('A', {M306A: 1.0}), ('B', {M306B: 2.0})], [('A', {M306A: 1.0}), ('B', {M306B: 0.5})],
    {M306A, M306B})
  assert rc['per'] == {'A': 0.0, 'B': 1.5} and not rc['unmatched'], rc
  rc = HC.rerun_cost([('A', {M83A: 1.0}), ('B', {M83B: 2.0})], [('A', {M83A: 0.0}), ('B', {M83B: 2.0})], {M83A, M83B})
  assert rc['per'] == {'A': 1.0, 'B': 0.0}, rc
def t_cost_his_own_loss():
  rc = HC.rerun_cost([('A', {M83A: 5.0}), ('B', {M83B: 3.0})], [('A', {M83A: 5.0}), ('B', {M83B: -7.0})], {M83A, M83B}, {M83A, M83B})
  assert rc['total'] == 10.0 and rc['his'] == 10.0, rc
def t_cost_same_id_both_passes():
  bl = (1, 2, 3)
  rc = HC.rerun_cost([('A', {bl: 1.0}), ('B', {bl: 1.0})], [('A', {bl: 0.0}), ('B', {bl: 0.0})], {bl})
  assert rc['total'] == 1.0, rc

# ---------------------------------------------------------------- unattributed gains
h1 = (1,); h2 = (2,); nb = (3,); un = (4,)
BEF4 = [('A', {h1: -10.0, nb: 5.0, un: 0.0, h2: -1.0})]
def plan_state(accepted=(), trial=()):
  st = HC.new_plan(list(trial) + list(accepted))
  st['trial'] = list(trial)
  st['accepted'] = set(accepted)
  return st
def t_unattributed_gain_no_credit():
  ev = HC.trial(BEF4, [('A', {h1: -10.3, nb: -5.0, un: 10.0, h2: -1.0})], {K1: {h1}}, [{h1, nb}], [K1], {})
  c = ev['components'][0]
  assert c['cost']['total'] == 10.3 and not c['ok'] and not ev['unattributed'], c
def t_cost_his_and_neighbours():
  ev = HC.trial(BEF4, [('A', {h1: -10.3, nb: -5.0, un: 10.0, h2: -1.0})], {K1: {h1}}, [{h1, nb}], [K1], {})
  c = ev['components'][0]
  assert c['cost']['his'] == 0.3 and c['cost']['neighbours'] == 10.0, c['cost']
def t_neighbour_loss_restores_first():
  ev = HC.trial(BEF4, [('A', {h1: -10.3, nb: -5.0, un: 10.0, h2: -1.0})], {K1: {h1}}, [{h1, nb}], [K1], {})
  st = plan_state(trial=[K1]); again = HC.plan(st, ev, {K1: 0.3}) == HC.RESTORE_FREE
  assert again and K1 in st['rejected'] and not st['accepted'] and st['trial'] == [], st
def t_loss_outside_rejected():
  ev = HC.trial(BEF4, [('A', {h1: -10.1, nb: 5.0, un: -1.0, h2: -1.0})], {K1: {h1}}, [{h1, nb}], [K1], {})
  st = plan_state(trial=[K1]); again = HC.plan(st, ev, {K1: 0.1}) == HC.RESTORE_FREE
  assert again and K1 in st['rejected'] and 'outside' in st['rejected'][K1], st['rejected']
def t_joint_trial_split():
  ev = HC.trial(BEF4, [('A', {h1: -10.1, nb: 5.0, un: -1.0, h2: -1.1})], {K1: {h1}, K2: {h2}},
    [{h1, nb}, {h2}], [K1, K2], {})
  st = plan_state(trial=[K1, K2]); again = HC.plan(st, ev, {K1: 0.1, K2: 0.2}) == HC.RERUN
  assert again and st['trial'] == [K1] and st['pending'] == [[K2]] and not st['rejected'] and not st['accepted'], st

# ---------------------------------------------------------------- new mismatch
def mismatch_case(same_clique):
  H2 = (50,)
  mids = {K1: {h1}, K2: {H2}}
  before = [('A', {h1: -10.0, H2: -1.0})]
  free = {K2: {'A': 1, 'B': 1}}
  def world(pins):
    now = {K2: {'A': 1, 'B': 0 if K1 in pins else 1}}
    mm, missing = HC.new_mismatches(free, now, pins)
    return [('A', {h1: -10.3 if K1 in pins else -10.0, H2: -1.0})], dict((k, mids[k]) for k in mm)
  cl = [{h1, H2}] if same_clique else [{h1}, {H2}]
  return simulate([K1], mids, before, world, cl, {K1: 0.3})
def t_new_mismatch_same_clique():
  st, runs = mismatch_case(True)
  assert K1 in st['rejected'] and 'mismatch' in st['rejected'][K1] and st['outcome'] == HC.RESTORE_FREE and len(runs) == 1, (st, runs)
def t_new_mismatch_other_clique():
  st, runs = mismatch_case(False)
  assert K1 in st['rejected'] and 'mismatch at His A 200' in st['rejected'][K1] and len(runs) == 1, (st, runs)

# ---------------------------------------------------------------- component policy
def policy(costs, cliques, first=None):
  """costs: {key: (mid, loss when fixed)}; the neighbourhood loses nothing."""
  mids = dict((k, {m}) for k, (m, l) in costs.items())
  base = dict((m, 0.0) for k, (m, l) in costs.items())
  def world(pins):
    return [('A', dict((m, -l if k in pins else 0.0) for k, (m, l) in costs.items()))], {}
  return simulate(sorted(costs), mids, [('A', base)], world, cliques,
    first or dict((k, l) for k, (m, l) in costs.items()))
def t_component_one_budget():
  st, runs = policy({K1: (h1, 0.30), K2: (h2, 0.30)}, [{h1, h2}])
  assert st['accepted'] == {K1} and K2 in st['rejected'] and runs[-1] == frozenset([K1]), (st, runs)
def t_separate_components():
  st, runs = policy({K1: (h1, 0.30), K2: (h2, 0.30)}, [{h1}, {h2}])
  assert st['accepted'] == {K1, K2} and len(runs) == 1, (st, runs)
def t_component_cheaper_kept():
  st, runs = policy({K1: (h1, 0.40), K2: (h2, 0.20)}, [{h1, h2}])
  assert st['accepted'] == {K2} and K1 in st['rejected'] and runs[1] == frozenset([K2]), (st, runs)
def t_component_within_budget():
  st, runs = policy({K1: (h1, 0.20), K2: (h2, 0.25)}, [{h1, h2}])
  assert st['accepted'] == {K1, K2} and len(runs) == 1, (st, runs)

# ---------------------------------------------------------------- cost
def t_n_clamp():
  # gain 5 in A, loss 1 in B: the gain does not cancel the loss
  rc = HC.rerun_cost([('A', {HA: 0.0}), ('B', {HB: 0.0})], [('A', {HA: 5.0}), ('B', {HB: -1.0})], {HA, HB})
  assert rc['total'] == 1.0, rc
def t_n_overlap():
  # two passes: blank Mover M in His 107's Clique in A and in His 200's in B
  M = (60,); P = (61,); Q = (62,)
  before = [('A', {P: 0.0, Q: 0.0, M: 0.0}), ('B', {})]
  after = [('A', {P: -0.2, Q: -0.6, M: 0.2}), ('B', {})]
  ev = HC.trial(before, after, {K1: {P}, K2: {Q}}, [{P, M}, {Q, M}], [K1, K2], {})
  cs = ev['components']
  assert len(cs) == 1 and cs[0]['keys'] == [K1, K2] and cs[0]['cost']['total'] == 0.6 and not cs[0]['ok'], cs
def t_n_transitive():
  # K1-K2 overlap through M, K2-K3 through R: one component
  M = (60,); R = (63,); P = (61,); Q = (62,); S = (64,)
  g = HC.components({K1: {P}, K2: {Q}, K3: {S}}, [{P, M}, {Q, M}, {Q, R}, {S, R}])
  assert [c['keys'] for c in g] == [[K1, K2, K3]], g

def t_n_star():
  # K1-K3 and K2-K3 overlap, K1-K2 do not: one component
  P = (61,); Q = (62,); S = (64,); M = (60,); R = (63,)
  g = HC.components({K1: {P}, K2: {Q}, K3: {S}}, [{P, M}, {S, M}, {Q, R}, {S, R}])
  assert [c['keys'] for c in g] == [[K1, K2, K3]], g

# ---------------------------------------------------------------- new: both inventories
def t_n_pass_change():
  # a neighbour written by pass A in the free run, by B in the rerun, 1 -> -5
  before = [('A', {h1: -10.0, nb: 1.0}), ('B', {})]
  after = [('A', {h1: -10.1}), ('B', {nb: -5.0})]
  ev = HC.trial(before, after, {K1: {h1}}, [{h1, nb}], [K1], {})
  c = ev['components'][0]
  assert not c['ok'] and c['cost']['unmatched'] == [(nb, 'written by pass B, before by A')], c
  assert c['cost']['total'] == 0.1 and c['cost']['changed'] == [('A', h1, -10.0, -10.1)], c
def t_n_disappear():
  before = [('A', {h1: -10.0, nb: 1.0})]
  after = [('A', {h1: -10.1})]
  ev = HC.trial(before, after, {K1: {h1}}, [{h1, nb}], [K1], {})
  c = ev['components'][0]
  assert not c['ok'] and c['cost']['unmatched'] == [(nb, 'absent from the rerun')], c
  # outside every neighbourhood: unattributed, the trial fails
  ev = HC.trial([('A', {h1: -10.0, un: 1.0})], [('A', {h1: -10.1})], {K1: {h1}}, [{h1}], [K1], {})
  st = plan_state(trial=[K1])
  HC.plan(st, ev, {K1: 0.1})
  assert ev['unattributedUnmatched'] and 'no valid comparison' in st['rejected'].get(K1, ''), (ev, st)
def t_n_small_losses_summed():
  u1 = (71,); u2 = (72,); u3 = (73,)
  before = [('A', {h1: 0.0, u1: 0.0, u2: 0.0, u3: 0.0})]
  ev = HC.trial(before, [('A', {h1: 0.0, u1: -0.003, u2: -0.003, u3: -0.003})], {K1: {h1}}, [{h1}], [K1], {})
  assert len(ev['unattributed']) == 3, ev
  ev = HC.trial(before, [('A', {h1: 0.0, u1: -0.003, u2: 0.0, u3: 0.0})], {K1: {h1}}, [{h1}], [K1], {})
  assert ev['unattributed'] == [], ev

# ---------------------------------------------------------------- new: planner (fresh confirmation)
def ev_of(comps, ua=(), uam=(), uu=()):
  return dict(components=[dict(keys=sorted(k), trial=sorted(t), cost=dict(total=c, unmatched=[]), mismatch=[],
    ok=c <= HC.BUDGET) for k, t, c in comps], unattributed=list(ua), unattributedMismatch=list(uam),
    unattributedUnmatched=list(uu))
def t_p_accepted_fails():
  # K1 fails while K2 is tried: nothing blamed, confirmation run with K1 alone
  st = plan_state(accepted=[K1], trial=[K2])
  out = HC.plan(st, ev_of([([K1], [], 0.9), ([K2], [K2], 0.1)]), {})
  assert out == HC.RERUN and st['trial'] == [] and not st['rejected'] and st['held'][0] == [K2], st
def t_p_always_confirm():
  # even after K1 passed alone before, a failing trial gets a fresh confirmation
  st = plan_state(trial=[K1])
  assert HC.plan(st, ev_of([([K1], [K1], 0.1)]), {}) == HC.DONE and st['accepted'] == {K1}
  st['trial'] = [K2]
  assert HC.plan(st, ev_of([([K1, K2], [K2], 0.9)]), {}) == HC.RERUN and st['trial'] == [] and not st['rejected'], st
def t_p_confirm_then_combination_fails():
  st = plan_state(accepted=[K1], trial=[K2])
  assert HC.plan(st, ev_of([([K1, K2], [K2], 0.9)]), {}) == HC.RERUN and st['trial'] == [] and not st['rejected'], st
  # K1 alone passes: the held combination is judged, K2 rejected; this run is the written one
  assert HC.plan(st, ev_of([([K1], [], 0.3)]), {}) == HC.DONE and K2 in st['rejected'] and st['accepted'] == {K1}, st
def t_p_confirm_fails_release_retry():
  st = plan_state(accepted=[K1], trial=[K2])
  HC.plan(st, ev_of([([K1], [], 0.9), ([K2], [K2], 0.1)]), {})
  out = HC.plan(st, ev_of([([K1], [], 0.9)]), {})
  assert out == HC.RERUN and 'released' in st['rejected'].get(K1, '') and st['trial'] == [K2] and not st['accepted'], st
  assert HC.plan(st, ev_of([([K2], [K2], 0.1)]), {}) == HC.DONE and st['accepted'] == {K2}, st
def t_p_attempt_bound():
  st = plan_state(accepted=[K1], trial=[K2])
  st['attempts'][K2] = HC.MAX_ATTEMPTS
  HC.plan(st, ev_of([([K1], [], 0.9), ([K2], [K2], 0.1)]), {})
  HC.plan(st, ev_of([([K1], [], 0.9)]), {})
  assert 'attempt limit' in st['rejected'].get(K2, ''), st
def t_p_single_guard():
  # K1 was tried alone already: in a failing joint trial it is rejected, K2 retried alone
  st = plan_state(trial=[K1, K2])
  st['single'] = {K1}
  HC.plan(st, ev_of([([K1, K2], [K1, K2], 0.9)]), {})
  assert K1 in st['rejected'] and st['trial'] == [K2], st
def t_p_random_worlds():
  # fresh-confirmation rule on 1,200 seeded random worlds, noise inside components and in
  # unpinned runs, random free_saved: termination, attempts, run bound, every outcome reached
  rnd = random.Random(7)
  worst = dict(attempts=0, runs=0)
  outcomes = {}
  for i in range(1200):
    keys = [(0, 'A', str(j), '') for j in range(rnd.randint(1, 6))]
    mids = dict((k, {(j,)}) for j, k in enumerate(keys))
    loss = dict((k, rnd.choice([0.1, 0.3, 0.6])) for k in keys)
    noise = rnd.random() < 0.4
    saved = rnd.random() < 0.7
    cl = [set().union(*(mids[k] for k in keys[j:j + 2])) for j in range(0, len(keys), rnd.choice([1, 2]))]
    def world(pins):
      sc = dict(((j,), -loss[k] - (0.5 if noise and rnd.random() < 0.2 else 0.0) if k in pins else 0.0)
        for j, k in enumerate(keys))
      sc[(99,)] = -1.0 if (noise and rnd.random() < 0.3) else 0.0
      return [('A', sc)], {}
    base = dict(((j,), 0.0) for j in range(len(keys))); base[(99,)] = 0.0
    st, runs = simulate(keys, mids, [('A', base)], world, cl, loss, saved)
    out = st['outcome']
    outcomes[out] = outcomes.get(out, 0) + 1
    worst['attempts'] = max([worst['attempts']] + list(st['attempts'].values()))
    worst['runs'] = max(worst['runs'], len(runs))
    assert len(runs) <= HC.max_runs(len(keys)), (len(runs), len(keys))
    if out == HC.DONE:
      assert runs[-1] == frozenset(st['accepted']), (runs, st)
    else:
      assert not st['accepted'] and st['reason'], st
      assert out == HC.RESTORE_FREE or not saved, (out, saved)
  assert worst['attempts'] <= HC.MAX_ATTEMPTS, worst
  assert all(outcomes.get(o) for o in (HC.DONE, HC.RESTORE_FREE, HC.FAILED)), outcomes
  print("  random worlds: max attempts %d, max reruns %d (bound for 6 His %d); %s" % (
    worst['attempts'], worst['runs'], HC.max_runs(6), sorted(outcomes.items())))
def t_p_unpinned_differs():
  # the trial is rejected, the unpinned rerun loses 6 and has a new mismatch
  st = plan_state(trial=[K1])
  out = HC.plan(st, ev_of([([K1], [K1], 0.9)]), {}, free_saved=False)
  assert out == HC.RERUN and K1 in st['rejected'] and st['trial'] == [], st
  ev = ev_of([], ua=[('A', (9,), 1.0, -5.0)], uam=[K2])
  assert HC.plan(dict(st, runs=1), ev, {}) == HC.RESTORE_FREE
  out = HC.plan(st, ev, {}, free_saved=False)
  assert out == HC.FAILED and 'loss 6.00' in st['reason'] and 'mismatch at His A 200' in st['reason'], st
  # the same run without a difference is the written result
  st2 = plan_state(trial=[])
  assert HC.plan(st2, ev_of([]), {}) == HC.DONE
def t_p_release_all():
  st = plan_state(accepted=[K1, K2], trial=[K3])
  HC.plan(st, ev_of([([K1], [], 0.1), ([K2], [], 0.1), ([K3], [K3], 0.9)]), {})
  out = HC.plan(st, ev_of([([K1], [], 0.1), ([K2], [], 0.1)], ua=[('A', (9,), 1.0, 0.0)]), {})
  assert out == HC.RERUN and not st['accepted'] and K1 in st['rejected'] and K2 in st['rejected'] and st['trial'] == [K3], st
def t_p_pins_equal_accepted():
  # K1 passes and K2 fails in one run: one more run with K1 alone before the end
  st, runs = policy({K1: (h1, 0.30), K2: (h2, 0.60)}, [{h1}, {h2}])
  assert runs == [frozenset([K1, K2]), frozenset([K1])] and st['accepted'] == {K1} and st['outcome'] == HC.DONE, (st, runs)
def t_p_run_bound():
  st = plan_state(trial=[K1])
  st['runs'] = HC.max_runs(st['n'])
  try:
    HC.plan(st, ev_of([([K1], [K1], 0.1)]), {})
  except AssertionError:
    pass
  else:
    raise AssertionError('run bound not asserted')

# ---------------------------------------------------------------- plan: confirmation runs
def t_p_requeue_after_confirmation():
  # K1 accepted fails (0.9) while K2's own component passes (0.1); K1
  # alone passes: K2 is run again once, not rejected
  st = plan_state(accepted=[K1], trial=[K2])
  assert HC.plan(st, ev_of([([K1], [], 0.9), ([K2], [K2], 0.1)]), {}) == HC.RERUN and st['trial'] == [], st
  assert HC.plan(st, ev_of([([K1], [], 0.3)]), {}) == HC.RERUN and st['trial'] == [K2] and not st['rejected'], st
  assert HC.plan(st, ev_of([([K1], [], 0.3), ([K2], [K2], 0.1)]), {}) == HC.DONE and st['accepted'] == {K1, K2}, st
  # the same failure again after a second passing confirmation: K2 is rejected
  st = plan_state(accepted=[K1], trial=[K2])
  for i in range(2):
    HC.plan(st, ev_of([([K1], [], 0.9), ([K2], [K2], 0.1)]), {})
    out = HC.plan(st, ev_of([([K1], [], 0.3)]), {})
  assert out == HC.DONE and 'accepted His A 107 fails' in st['rejected'].get(K2, '') and st['accepted'] == {K1}, st
  # K2's own component failed too: judged at once after the passing confirmation, no rerun
  st = plan_state(accepted=[K1], trial=[K2])
  HC.plan(st, ev_of([([K1], [], 0.9), ([K2], [K2], 0.9)]), {})
  out = HC.plan(st, ev_of([([K1], [], 0.3)]), {})
  assert out == HC.DONE and K2 in st['rejected'] and st['runs'] == 2, st
def t_p_restore_keeps_pending():
  # nothing accepted, but trials pending: they run before the free result is restored
  st = plan_state(trial=[K1, K2])
  assert HC.plan(st, ev_of([([K1, K2], [K1, K2], 0.9)]), {K1: 0.1, K2: 0.2}) == HC.RERUN and st['trial'] == [K1], st
  assert HC.plan(st, ev_of([([K1], [K1], 0.9)]), {}) == HC.RERUN and st['trial'] == [K2], st
  out = HC.plan(st, ev_of([([K2], [K2], 0.9)]), {})
  assert out == HC.RESTORE_FREE and st['reason'] == 'nothing accepted' and st['runs'] == 3, st
def t_p_component_release():
  # a failing confirmation releases the whole accepted component, both His
  st = plan_state(accepted=[K1, K2], trial=[K3])
  HC.plan(st, ev_of([([K1, K2], [], 0.1), ([K3], [K3], 0.9)]), {})
  out = HC.plan(st, ev_of([([K1, K2], [], 0.9)]), {})
  assert out == HC.RERUN and 'released' in st['rejected'].get(K1, '') and 'released' in st['rejected'].get(K2, ''), st
  assert st['trial'] == [K3] and not st['accepted'], st
def t_p_eval_of_other_run():
  st = plan_state(accepted=[K1], trial=[K2])
  try:
    HC.plan(st, ev_of([([K1], [], 0.1)]), {})
  except AssertionError as e:
    assert 'evaluation of another run' in str(e), str(e)
  else:
    raise AssertionError('evaluation of another run accepted')
def t_p_duplicate_candidates():
  try:
    HC.new_plan([K1, K1])
  except AssertionError:
    pass
  else:
    raise AssertionError('duplicate candidates accepted')
def t_p_missing_his():
  # a His of the free run absent from the rerun (no Mover id): unattributed, the trial fails
  ev = HC.trial([('A', {h1: -10.0})], [('A', {h1: -10.1})], {K1: {h1}}, [{h1}], [K1], {}, missing=[K2])
  st = plan_state(trial=[K1])
  HC.plan(st, ev, {})
  assert ev['unattributedMissing'] == [K2] and 'His A 200 absent from the rerun' in st['rejected'].get(K1, ''), (ev, st)

# ---------------------------------------------------------------- check_written
RING = ('NE2', 'CE1', 'ND1', 'CD2', 'CG')
def xyz(flipped):
  # flipped twin: ND1/CD2 and CE1/NE2 swap places about CB-CG; the ring H follow their N
  p = dict(CG=(0.0, 0.0, 0.0), ND1=(1.0, 1.0, 0.0), CD2=(1.0, -1.0, 0.0), CE1=(2.0, 0.7, 0.0), NE2=(2.0, -0.7, 0.0),
    HD1=(1.2, 2.0, 0.0), HE2=(2.8, -1.4, 0.0))
  if flipped:
    p = dict(CG=p['CG'], ND1=p['CD2'], CD2=p['ND1'], CE1=p['NE2'], NE2=p['CE1'], HD1=(1.2, -2.0, 0.0), HE2=(2.8, 1.4, 0.0))
  return p
def fixed_state(dels, flipped=False):
  p = xyz(flipped)
  return dict(dels=dels, coarse=p, fixup=dict(p))
def written_from(st, writer=None):
  return dict(present=dict((h, not d) for h, d in st['dels'].items()),
    xyz=dict((n, v) for n, v in st['fixup'].items() if not st['dels'].get(n)),
    acceptor={'NE2': st['dels']['HE2'], 'ND1': st['dels']['HD1']}, writer=writer or {})
HE2_ONLY = {'HE2': False, 'HD1': True}
HD1_ONLY = {'HE2': True, 'HD1': False}
def t_w_pass():
  f = {'A': fixed_state(HE2_ONLY), 'B': fixed_state(HE2_ONLY)}
  r = HC.check_written(f, dict((c, written_from(s)) for c, s in f.items()))
  assert not r['errors'] and r['checked'] == dict(H=4, xyz=12, N=4), r
  assert r['skipped'] == [('A', 'HD1', 'deleted H'), ('B', 'HD1', 'deleted H')], r
def t_w_h_wrong():
  f = {'A': fixed_state(HE2_ONLY)}
  w = {'A': written_from(fixed_state(HD1_ONLY))}
  w['A']['acceptor'] = {'NE2': False, 'ND1': True}
  r = HC.check_written(f, w)
  assert len([e for e in r['errors'] if ' kept, ' in e or ' deleted, ' in e]) == 2, r
def t_w_retained_h_moved():
  f = {'A': fixed_state(HE2_ONLY)}
  w = {'A': written_from(f['A'])}
  w['A']['xyz']['HE2'] = (12.8, -1.4, 0.0)
  r = HC.check_written(f, w)
  assert r['errors'] == ['His conformer A HE2 written 10.000 A from the fixed state'], r
def t_w_flipped_twin():
  f = {'A': fixed_state(HE2_ONLY)}
  w = {'A': written_from(fixed_state(HE2_ONLY, flipped=True))}
  r = HC.check_written(f, w)
  assert len(r['errors']) == 5 and all('from the fixed state' in e for e in r['errors']), r
def t_w_missing_conformer():
  f = {'A': fixed_state(HE2_ONLY), 'B': fixed_state(HE2_ONLY)}
  r = HC.check_written(f, {'A': written_from(f['A'])})
  assert r['errors'] == ['His conformer B: no written record'], r
def t_w_empty_xyz():
  f = {'A': fixed_state(HE2_ONLY)}
  w = {'A': written_from(f['A'])}
  w['A']['xyz'] = {}
  r = HC.check_written(f, w)
  assert len(r['errors']) == 6 and all('coordinates not recorded' in e for e in r['errors']), r
def t_w_missing_present():
  f = {'A': fixed_state(HE2_ONLY)}
  w = {'A': written_from(f['A'])}
  del w['A']['present']['HE2']
  r = HC.check_written(f, w)
  assert r['errors'] == ['His conformer A HE2: presence not recorded'], r
def t_w_shared_blank_n():
  # 1us0 His 83 layout: blank ring, split HE2/HD1; A keeps HE2, B keeps HD1; pass B writes the ring
  f = {'A': fixed_state(HE2_ONLY), 'B': fixed_state(HD1_ONLY)}
  wr = dict((n, 'B') for n in RING)
  good = dict(A=written_from(f['A'], wr), B=written_from(f['B'], wr))
  good['A']['acceptor'] = dict(good['B']['acceptor'])       # one blank N, written by B
  r = HC.check_written(f, good)
  assert not r['errors'] and r['checked']['N'] == 4 and r['checked']['xyz'] == 12, r
  # each conformer's record of the shared N is checked: A's record with A's own flags fails
  bad = dict(A=written_from(f['A'], wr), B=written_from(f['B'], wr))
  r = HC.check_written(f, bad)
  assert len(r['errors']) == 2 and all(e.startswith('His conformer A N') for e in r['errors']), r
  # swapped records (A carries B's flags, B carries A's): B's record fails
  sw = dict(A=written_from(f['A'], wr), B=written_from(f['B'], wr))
  sw['A']['acceptor'], sw['B']['acceptor'] = dict(sw['B']['acceptor']), dict(sw['A']['acceptor'])
  r = HC.check_written(f, sw)
  assert len(r['errors']) == 2 and all(e.startswith('His conformer B N') for e in r['errors']), r
def t_w_unfixed_writer():
  # ring atoms written by an unfixed (blank) pass are errors for a reconciled His:
  # ring xyz 5 + ring N 2; skips: the deleted HD1 and the unfixed conformer C
  f = {'A': fixed_state(HE2_ONLY)}
  w = {'A': written_from(f['A'], dict((n, '') for n in RING)), 'C': written_from(f['A'])}
  r = HC.check_written(f, w)
  assert len(r['errors']) == 7 and all('written by pass blank, not fixed' in e for e in r['errors']), r
  assert r['skipped'] == [('A', 'HD1', 'deleted H'), ('C', 'conformer', 'not fixed')], r
  assert r['checked'] == dict(H=2, xyz=1, N=0), r
def t_w_small_offset():
  # 0.05 A off one ring atom (fix-up moves are small): an error
  f = {'A': fixed_state(HE2_ONLY)}
  w = {'A': written_from(f['A'])}
  w['A']['xyz']['CE1'] = (2.05, 0.7, 0.0)
  r = HC.check_written(f, w)
  assert r['errors'] == ['His conformer A CE1 written 0.050 A from the fixed state'], r
def t_w_incomplete_inventory():
  # a fixed state holding only HE2: the missing ring atoms are an error, not an empty check
  f = {'A': dict(dels=HE2_ONLY, coarse={'HE2': xyz(False)['HE2']}, fixup={})}
  w = {'A': written_from(fixed_state(HE2_ONLY))}
  r = HC.check_written(f, w)
  assert r['errors'] == ['His conformer A: fixed state lacks CG ND1 CD2 CE1 NE2'], r
def t_w_fixup_none():
  f = {'A': dict(fixed_state(HE2_ONLY), fixup=None)}
  try:
    HC.check_written(f, {'A': written_from(fixed_state(HE2_ONLY))})
  except AssertionError as e:
    assert 'fixup must be a dict' in str(e), str(e)
  else:
    raise AssertionError('fixup=None accepted')
def t_w_skip_bond_fixup():
  f = dict(dels=HE2_ONLY, coarse=xyz(False), fixup=dict(xyz(False), CE1=(2.0, 0.9, 0.0)))
  w = {'A': written_from(dict(dels=HE2_ONLY, fixup=f['coarse']))}
  assert not HC.check_written({'A': f}, w, skip_bond_fixup=True)['errors']
  assert HC.check_written({'A': f}, w, skip_bond_fixup=False)['errors']
def t_w_empty_fixup():
  # empty FixUp positions: the coarse positions are retained and checked, not skipped
  f = dict(dels=HE2_ONLY, coarse=xyz(False), fixup={})
  w = {'A': written_from(dict(dels=HE2_ONLY, fixup=f['coarse']))}
  r = HC.check_written({'A': f}, w)
  assert not r['errors'] and r['checked']['xyz'] == 6, r
  w['A']['xyz']['CG'] = (0.0, 0.0, 1.0)
  assert HC.check_written({'A': f}, w)['errors']

# ---------------------------------------------------------------- check_fixed
def t_n_fixed_missing_conformer():
  try:
    HC.check_fixed(K1, {'A': 2, 'B': 2}, {'A': 2, 'B': 2}, {'A': 0})
  except AssertionError as e:
    assert 'missing from the final state' in str(e), str(e)
  else:
    raise AssertionError('no exception')

def t_x_fixed_all_final():
  # every final conformer counts, not only the pinned ones
  try:
    HC.check_fixed(K1, {'A': 2, 'B': 2}, {'A': 2, 'B': 2}, {'A': 0, 'B': 0, 'C': 1})
  except AssertionError as e:
    assert 'charges differ' in str(e), str(e)
  else:
    raise AssertionError('no exception')

TESTS = [(n, make_decide_test(n, want, recs)) for n, name, recs, want in DECIDE] + [
  (n[2:], f) for n, f in sorted(globals().items()) if n.startswith('t_')]

def run():
  failed = []
  for name, f in TESTS:
    try:
      f()
      print("PASS", name)
    except Exception as e:
      print("FAIL", name, type(e).__name__, str(e)[:300])
      failed.append(name)
  assert not failed, failed

if __name__ == "__main__":
  run()
  print("OK")
