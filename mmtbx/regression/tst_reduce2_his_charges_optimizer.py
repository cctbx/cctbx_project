"""reduce2 with reconcile_his_charges=True: the Optimizer integration of
mmtbx.reduce.HisCharges. Fixtures: residues within 6 A of a His ring, the
symmetry mates written out as chains (no crystal symmetry; the fragments without
the mates lose the case):
- pdb_1exr: 1exr His A 107 (A +1, B 0 -> 0); mates (0,-1,1) as chain B and
  (1,0,0) as chain C.
- pdb_1us0_163: 1us0 His A 163 (A +1, B 0 -> +1); mate -x,y+1/2,-z (1,0,1) as B.
- pdb_1us0_83: 1us0 His A 83, a blank ring with ring H split A/B (no case).
- pdb_1exr_sym: 1exr His A 107 with the crystal symmetry kept: the residues
  whose symmetry mates come within 6 A of the ring (the Optimizer's symmetry
  copies supply the mates; for the copies' atom info).
Built from these: two candidates in one model, a Zn at conformer B's NE2 or
conformer B's ND1 (kept), two models.
"""
from __future__ import absolute_import, division, print_function
import os, re
from six.moves import cStringIO as StringIO
import iotbx.pdb
import mmtbx.model
from libtbx.utils import Sorry, Abort
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.reduce import Optimizers, HisCharges

LOG = [None]   # the last run's log

def run(pdb, prefix, *extra):
  """reduce2 on pdb with extra arguments. Returns (output ATOM records without
  serials, description, the Optimizers made); the log in LOG[0]."""
  fn_in, fn_out, fn_desc = prefix + ".pdb", prefix + "_out.pdb", prefix + "_out.txt"
  with open(fn_in, "w") as f:
    f.write(pdb)
  for fn in (fn_out, fn_desc):
    if os.path.exists(fn): os.remove(fn)
  made = []
  orig = Optimizers.Optimizer.__init__
  def init(self, *a, **k):
    made.append(self)
    orig(self, *a, **k)
  Optimizers.Optimizer.__init__ = init
  try:
    LOG[0] = StringIO()
    run_program(program_class=reduce2.Program, logger=LOG[0], args=[fn_in,
      "output.filename=" + fn_out, "output.description_file_name=" + fn_desc,
      "overwrite=True"] + list(extra))
  finally:
    Optimizers.Optimizer.__init__ = orig
  with open(fn_out) as f:
    atoms = [l[12:80] for l in f.read().splitlines() if l.startswith(("ATOM", "HETATM"))]
  with open(fn_desc) as f:
    desc = f.read()
  return atoms, desc, made

ON = "reconcile_his_charges=True"

def ring_h(atoms, chain, resseq):
  return sorted((l[0:4].strip(), l[4]) for l in atoms if l[9] == chain and
    l[10:14].strip() == resseq and l[0:4] in (" HD1", " HE2"))

def only(made):
  assert len(made) == 1, len(made)
  his = made[0].getHisCharges()
  assert len(his) == 1, his
  return his[0]

def case(m, resseq):
  cs = [c for c in m['cases'] if c['key'][2] == resseq]
  assert len(cs) == 1, m['cases']
  return cs[0]

def removed(off, on):
  return sorted(set(off) - set(on)), sorted(set(on) - set(off))

class patched(object):
  """Set attributes for the block; the old values are always put back."""
  def __init__(self, *triples):
    self.triples = triples
  def __enter__(self):
    self.old = [(o, n, getattr(o, n)) for o, n, v in self.triples]
    for o, n, v in self.triples:
      setattr(o, n, v)
  def __exit__(self, *exc):
    for o, n, v in reversed(self.old):
      setattr(o, n, v)

def failing_check(on_call=None):
  """check_written with an injected error (on every call, or on call on_call)."""
  orig = HisCharges.check_written
  calls = []
  def bad(*a, **k):
    r = orig(*a, **k)
    calls.append(1)
    if on_call is None or len(calls) == on_call:
      r['errors'].append('injected')
    return r
  return (HisCharges, 'check_written', bad)

def test_done_1exr():
  """1exr His A 107 (6 A, symmetry mates as chains B and C): A +1, B 0 -> 0 by
  removing HD1 A; one rerun; the description section and the marked flip."""
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_1")
  assert made[0].getHisCharges() == [], made[0].getHisCharges()
  assert "His charges across" not in desc_off
  on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_1on", ON)
  m = only(made)
  assert m['outcome'] == HisCharges.DONE and m['code'] == HisCharges.RECONCILED, m
  assert m['reruns'] == 1 and m['cap'] == 10 and m['maxRuns'] == 10, m
  assert m['backstopSeconds'] == HisCharges.WALL_BACKSTOP_FACTOR * m['freeSeconds'], m
  c = case(m, '107')
  assert c['free'] == {'A': 1, 'B': 0} and c['final'] == {'A': 0, 'B': 0} and c.get('fixed'), c
  assert c['category'] == 'fixed', c
  assert abs(c['first']['costs']['0']['total'] - 0.294) < 0.002, c['first']
  gone, new = removed(off, on)
  assert len(gone) == 1 and gone[0][:15] == " HD1AHIS A 107 " and new == [], (gone, new)
  assert ("His A 107: A +1 B 0 -> A 0 B 0; set to charge 0; screen cost 0.29 (A 0.29, B 0.00); "
    "rerun cost 0.29 (His 0.29, neighbours 0.00); suggested flip A 22.80") in desc, desc
  line = [l for l in desc.splitlines() if l.startswith("  model: reconciled; first run ")]
  assert len(line) == 1 and "; 1 rerun(s), cap 10 (max_runs 10), " in line[0], desc
  assert line[0].endswith("; rejected: none; unmatched Movers: 0; skipped checks: 2"), line
  flip = [l for l in desc.splitlines() if "chain A HIS 107 conformer A: flipped" in l]
  assert len(flip) == 1 and flip[0].endswith("(state fixed by charge reconciliation)"), flip

def test_check_failure_restores():
  """A failed written-state check: production restores the free result with a
  warning; strict mode restores, then raises."""
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_2")
  with patched(failing_check()):
    on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_2on", ON)
    m = only(made)
    assert m['outcome'] == HisCharges.RESTORE_FREE and m['code'] == HisCharges.INTERNAL, m
    assert m['detail'] == 'His A 107 written state: injected', m
    assert on == off, removed(off, on)
    assert "Warning: His charges: reconciliation stopped by an error, the first result is kept: ContractError: " \
      "His A 107 written state: injected" in LOG[0].getvalue(), LOG[0].getvalue()
    c = case(m, '107')
    assert c['final'] == c['free'] == {'A': 1, 'B': 0} and not c.get('fixed'), c
    assert c['category'] == 'restored', c
    with Optimizers.his_charges_strict():
      try:
        run(pdb_1exr, "tst_reduce2_his_charges_opt_2s", ON)
      except HisCharges.ContractError as e:
        assert 'injected' in str(e), e
      else:
        raise AssertionError("strict mode did not raise")
  assert Optimizers.Optimizer.hisChargesStrict is False

def strict_run(pdb, prefix, *patches):
  """Strict run expected to raise ContractError; returns the Optimizer."""
  made = []
  init = Optimizers.Optimizer.__init__
  def rec(self, *a, **k):
    made.append(self)
    init(self, *a, **k)
  with patched((Optimizers.Optimizer, '__init__', rec), *patches):
    with Optimizers.his_charges_strict():
      try:
        run(pdb, prefix, ON)
      except HisCharges.ContractError as e:
        return made[-1], str(e)
  raise AssertionError("strict mode did not raise")

def hd1_a_deleted(opt):
  return [h.id_str() for h in opt._hydrogensToDelete if h.name.strip() == 'HD1' and h.parent().altloc == 'A']

def test_strict_restores_first():
  """Strict mode: the free result is restored before the raise (HD1 A not deleted)."""
  opt, e = strict_run(pdb_1exr, "tst_reduce2_his_charges_opt_3", failing_check())
  assert hd1_a_deleted(opt) == [], hd1_a_deleted(opt)
  assert opt.getHisCharges()[0]['code'] == HisCharges.INTERNAL

def test_unpin_contract():
  """_HisUnpin's two checks (the one-state Mover's index, the pinned state's
  positions) raise ContractError: production restores with a warning, strict
  mode restores, then raises."""
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_u")
  orig = Optimizers.Optimizer._HisUnpin
  def moved(self, his):
    for m in his['pinned']:
      a = m.CoarsePositions().atoms[0]
      a.xyz = (a.xyz[0] + 0.5, a.xyz[1], a.xyz[2])
      break
    return orig(self, his)
  class Index(object):
    def __init__(self, o):
      self.o = o
    def GetCoarseLocation(self, m):
      return 1
  def index(self, his):
    for m in his['pinned']:
      self._optimizerOf[m] = Index(self._optimizerOf[m])
    return orig(self, his)
  for tag, f, what in (('p', moved, 'not at its positions'), ('i', index, 'not at index 0')):
    with patched((Optimizers.Optimizer, '_HisUnpin', f)):
      on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_u" + tag, ON)
    m = only(made)
    assert m['code'] == HisCharges.INTERNAL and what in m['detail'], m
    assert on == off, removed(off, on)
    assert "Warning: His charges: reconciliation stopped by an error" in LOG[0].getvalue()
    opt, e = strict_run(pdb_1exr, "tst_reduce2_his_charges_opt_us" + tag, (Optimizers.Optimizer, '_HisUnpin', f))
    assert what in e and hd1_a_deleted(opt) == [], (e, hd1_a_deleted(opt))

def test_restore_check_catches_omission():
  """Mutations: the save omits the deletions, the model atoms' atom info, or the
  symmetry copies' atom info; the independent restore check stops the run with
  a Sorry instead of returning the rerun's state. The copies need the crystal
  symmetry (pdb_1exr_sym)."""
  orig_save = Optimizers.Optimizer._HisSave
  for item, pdb, where in (('deleted', pdb_1exr, 'deleted'), ('info', pdb_1exr, 'model'),
                           ('copyInfo', pdb_1exr_sym, 'copies')):
    def save(self, myModel):
      s = orig_save(self, myModel)
      s[item] = []
      return s
    with patched(failing_check(), (Optimizers.Optimizer, '_HisSave', save)):
      try:
        run(pdb, "tst_reduce2_his_charges_opt_4" + item, ON)
      except Sorry as e:
        assert "restored first result differs from it (%s)" % where in str(e), (item, e)
        assert str(e).endswith("; restoring after ContractError: His A 107 written state: injected"), (item, e)
      else:
        raise AssertionError("omitted %s not caught" % item)
  # without the omission the symmetry fixture restores cleanly
  with patched(failing_check()):
    on, desc, made = run(pdb_1exr_sym, "tst_reduce2_his_charges_opt_4s", ON)
  assert only(made)['code'] == HisCharges.INTERNAL
  assert len(made[0]._passCopies) > 0

class _Clock(object):
  """time for Optimizers: the real clock plus an offset, switched on after a given
  call of an Optimizer method."""
  def __init__(self, after_call, offset):
    import time
    self._t, self.after, self.offset, self.calls, self.on = time, after_call, offset, 0, False
  def time(self):
    return self._t.time() + (self.offset if self.on else 0.0)
  def __getattr__(self, n):
    return getattr(self._t, n)

def clocked(method, after_call, offset, prefix, *extra):
  clock = _Clock(after_call, offset)
  orig = getattr(Optimizers.Optimizer, method)
  def f(self, *a, **k):
    r = orig(self, *a, **k)
    clock.calls += 1
    if clock.calls == clock.after:
      clock.on = True
    return r
  with patched((Optimizers, 'time', clock), (Optimizers.Optimizer, method, f)):
    return run(pdb_1exr, prefix, ON, *extra)

def test_deadline():
  """The search deadline (the free run reaches it; a rerun reaches its own) and
  the wall-time backstop restore the free result; a rerun gets its own search
  time, not the free run's deadline."""
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_5")
  # the free run ends past its deadline: no rerun
  on, desc, made = clocked('_OptimizeModel', 1, 1.0e11, "tst_reduce2_his_charges_opt_5f")
  m = only(made)
  assert (m['code'], m['detail'], m['reruns']) == (HisCharges.SEARCH_DEADLINE, 'the first run reached it', 0), m
  assert on == off, removed(off, on)
  # the first rerun ends past its own search time
  on, desc, made = clocked('_OptimizeModel', 2, 1.0e11, "tst_reduce2_his_charges_opt_5a")
  m = only(made)
  assert (m['code'], m['detail'], m['reruns']) == (HisCharges.SEARCH_DEADLINE, 'rerun 1 reached it', 1), m
  assert on == off, removed(off, on)
  # 1000 s pass after the free run: past 50 x its time
  on, desc, made = clocked('_HisByResidue', 1, 1.0e3, "tst_reduce2_his_charges_opt_5b")
  m = only(made)
  assert (m['code'], m['reruns']) == (HisCharges.WALL_BACKSTOP, 0), m
  assert m['reason'] == 'first result restored: wall-time backstop reached (50 x the first run)', m
  assert on == off, removed(off, on)
  # a 1000 s search limit; 2000 s pass before the rerun starts: the rerun has its
  # own 1000 s (the free run's deadline is past), so it is DONE
  on, desc, made = clocked('_HisRestore', 1, 2.0e3, "tst_reduce2_his_charges_opt_5r",
    "treat_cliques_as_singletons_after_seconds=1000")
  m = only(made)
  assert m['code'] == HisCharges.RECONCILED and m['reruns'] == 1, m

def failing_call(owner, name, on_call, exc, after=True):
  """owner.name raising exc on call on_call (after the original ran, or instead
  of it)."""
  orig = getattr(owner, name)
  calls = []
  def f(*a, **k):
    calls.append(1)
    if len(calls) == on_call and not after:
      raise exc
    r = orig(*a, **k)
    if len(calls) == on_call:
      raise exc
    return r
  return (owner, name, f)

def test_any_error_keeps_first():
  """Any error after the first result is saved, not only a ContractError: a
  RuntimeError or a Sorry at the end of the first rerun (the model then holds
  the rerun's state) and an error in decide() (no rerun) keep the first result
  (output equal to flag off), with a warning naming the exception, outcome code
  INTERNAL and no His fixed. Strict mode re-raises the original exception after
  the restore."""
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_e")
  for tag, p, cls, reruns in (
      ('r', lambda: failing_call(Optimizers.Optimizer, '_OptimizeModel', 2, RuntimeError('injected')), RuntimeError, 1),
      ('s', lambda: failing_call(Optimizers.Optimizer, '_OptimizeModel', 2, Sorry('injected')), Sorry, 1),
      ('d', lambda: failing_call(HisCharges, 'decide', 1, ValueError('injected'), after=False), ValueError, 0)):
    with patched(p()):
      on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_e" + tag, ON)
    m = only(made)
    name = cls.__name__
    assert (m['outcome'], m['code'], m['reruns']) == (HisCharges.RESTORE_FREE, HisCharges.INTERNAL, reruns), m
    assert m['detail'] == name + ': injected', m
    assert on == off, (tag, removed(off, on))
    assert ("Warning: His charges: reconciliation stopped by an error, the first result is kept: %s: injected"
      % name) in LOG[0].getvalue(), LOG[0].getvalue()
    c = case(m, '107')
    assert c['final'] == c['free'] == {'A': 1, 'B': 0} and not c.get('fixed'), c
    assert c['category'] == 'restored', c
    assert "His A 107: A +1 B 0; kept: first result kept: internal error" in desc, desc
    made = []
    init = Optimizers.Optimizer.__init__
    def rec(self, *a, **k):
      made.append(self)
      init(self, *a, **k)
    with patched((Optimizers.Optimizer, '__init__', rec), p()):
      with Optimizers.his_charges_strict():
        try:
          run(pdb_1exr, "tst_reduce2_his_charges_opt_es" + tag, ON)
        except cls as e:
          assert str(e) == 'injected', e
        else:
          raise AssertionError("strict mode did not raise")
    assert hd1_a_deleted(made[-1]) == [], hd1_a_deleted(made[-1])
    assert made[-1].getHisCharges()[0]['code'] == HisCharges.INTERNAL
  assert Optimizers.Optimizer.hisChargesStrict is False

def test_first_run_error_propagates():
  """An error in the first run propagates the same with the flag off and on."""
  for extra in ((), (ON,)):
    with patched(failing_call(Optimizers.Optimizer, '_OptimizeModel', 1, RuntimeError('injected first'))):
      try:
        run(pdb_1exr, "tst_reduce2_his_charges_opt_f", *extra)
      except RuntimeError as e:
        assert str(e) == 'injected first', e
      else:
        raise AssertionError("no error with %s" % (extra,))

def partial_record(on_call, exc):
  """_HisRecordPass failing partway on call on_call: each His Mover is left in a
  state it did not choose (as the record loop sets them in turn), then exc."""
  orig = Optimizers.Optimizer._HisRecordPass
  calls = []
  def f(self, *a, **k):
    calls.append(1)
    if len(calls) != on_call:
      return orig(self, *a, **k)
    for m in self._movers:
      if isinstance(m, Optimizers.Movers.MoverHisFlip):
        pr = m.CoarsePositions()
        chosen = self._optimizerOf[m].GetCoarseLocation(m)
        self._setMoverState(pr, (chosen + 1) % len(pr.positions))
    raise exc
  return (Optimizers.Optimizer, '_HisRecordPass', f)

def test_first_run_record_error():
  """An error in code that runs only with the flag, in the first run (the
  record pass, partway, in the first or the last pass; the save before or
  after it): the passes run again without reconciliation, so the output equals
  flag off; a warning naming the exception, code FIRST_RUN_ERROR, no His case,
  the description line. Strict mode re-raises it. A MemoryError propagates."""
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_p")
  calls = []
  orig = Optimizers.Optimizer._HisRecordPass
  def count(self, *a, **k):
    calls.append(1)
    return orig(self, *a, **k)
  with patched((Optimizers.Optimizer, '_HisRecordPass', count)):
    run(pdb_1exr, "tst_reduce2_his_charges_opt_pc", ON)
  assert len(calls) == 4, calls    # two passes: the first run (calls 1, 2) and one rerun
  for tag, p, cls in (
      ('1', lambda: partial_record(1, RuntimeError('injected')), RuntimeError),
      ('2', lambda: partial_record(2, Sorry('injected')), Sorry),
      ('s', lambda: failing_call(Optimizers.Optimizer, '_HisSave', 1, ValueError('injected'), after=False), ValueError),
      ('t', lambda: failing_call(Optimizers.Optimizer, '_HisSave', 2, KeyError('injected'), after=False), KeyError)):
    with patched(p()):
      on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_p" + tag, ON)
    assert len(made) == 1, made
    his = made[0].getHisCharges()
    assert len(his) == 1, his
    m = his[0]
    detail = '%s: %s' % (cls.__name__, cls('injected'))
    assert (m['outcome'], m['code'], m['reruns'], m['cases']) == (None, HisCharges.FIRST_RUN_ERROR, 0, []), m
    assert m['detail'] == detail, m
    assert on == off, (tag, removed(off, on))
    assert ("Warning: His charges: an error in the first run, the model was optimized again "
      "without reconciliation: " + detail) in LOG[0].getvalue(), LOG[0].getvalue()
    assert ("  model: not reconciled: error in the first run (%s); first run " % detail) in desc, desc
    assert "; 1 model(s) not checked (error in the first run)" in desc, desc
    with patched(p()):
      with Optimizers.his_charges_strict():
        try:
          run(pdb_1exr, "tst_reduce2_his_charges_opt_ps" + tag, ON)
        except cls as e:
          assert type(e) is cls, e
        else:
          raise AssertionError("strict mode did not raise")
  with patched(partial_record(1, MemoryError('injected'))):
    try:
      run(pdb_1exr, "tst_reduce2_his_charges_opt_pm", ON)
    except MemoryError as e:
      assert str(e) == 'injected', e
    else:
      raise AssertionError("MemoryError not re-raised")

def test_fallback_deadline():
  """The fallback run's search time starts after the restore: a first-pass
  error, then a restore that takes 2000 s (clock offset). The fallback still gets
  the first run's allowance (1000 s), the output equals flag off, and the
  failed run's seconds reported include the restore."""
  dl = "treat_cliques_as_singletons_after_seconds=1000"
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_q", "reconcile_his_charges=False", dl)
  clock = _Clock(1, 2000.0)
  O = Optimizers.Optimizer
  restore, optimize = O._HisRestore, O._OptimizeModel
  seen = []
  def slow_restore(self, *a, **k):
    r = restore(self, *a, **k)
    clock.on = True
    return r
  def budget(self, mi, model, ctx, his=None):
    seen.append((his is not None, ctx['stopOptimizingAt'] - clock.time()))
    return optimize(self, mi, model, ctx, his)
  with patched((Optimizers, 'time', clock), (O, '_HisRestore', slow_restore), (O, '_OptimizeModel', budget),
      partial_record(1, RuntimeError('injected'))):
    on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_qon", ON, dl)
  assert [x for x, b in seen] == [True, False] and all(990 < b <= 1000 for x, b in seen), seen
  m = only(made)
  assert m['code'] == HisCharges.FIRST_RUN_ERROR and m['freeSeconds'] >= 2000, m
  assert on == off, removed(off, on)
  # a failed first save that takes 2000 s: the same allowance, the save credited
  clock = _Clock(1, 2000.0)
  save, seen = O._HisSave, []
  def slow_failed_save(self, *a, **k):
    save(self, *a, **k)
    clock.on = True
    raise ValueError('injected')
  with patched((Optimizers, 'time', clock), (O, '_HisSave', slow_failed_save), (O, '_OptimizeModel', budget)):
    on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_qs", ON, dl)
  assert [x for x, b in seen] == [False] and all(990 < b <= 1000 for x, b in seen), seen
  m = only(made)
  assert m['code'] == HisCharges.FIRST_RUN_ERROR and m['freeSeconds'] >= 2000, m
  assert on == off, removed(off, on)

def test_report_errors():
  """An error in the His report (Optimizer._HisReport, reduce2's _HisChargesText)
  after the model is final: qmready (Auto) on the 1exr fixture writes the same
  model and total charge as without the error, with a warning and a one-line
  section; strict mode re-raises."""
  ok, desc_ok, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_w", "qmready=True")
  assert "Total charge: 0" in LOG[0].getvalue()
  for tag, owner, name, warning, line in (
      ('o', Optimizers.Optimizer, '_HisReport', "Warning: His charges: the report could not be written: "
        "RuntimeError: injected", "  model: His charges not reported: RuntimeError: injected"),
      ('r', reduce2.Program, '_HisChargesText', "Warning: His charges: the description section could not "
        "be written: RuntimeError: injected", "His charges across conformers (reconcile_his_charges=True): "
        "could not be written: RuntimeError: injected")):
    with patched(failing_call(owner, name, 1, RuntimeError('injected'), after=False)):
      atoms, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_w" + tag, "qmready=True")
    assert atoms == ok, (tag, removed(ok, atoms))
    assert "Total charge: 0" in LOG[0].getvalue(), tag
    assert warning in LOG[0].getvalue(), (tag, LOG[0].getvalue())
    assert line in desc.splitlines(), (tag, desc)
    if tag == 'o':
      assert ("His charges across conformers (reconcile_his_charges=True): counts not available "
        "(the report could not be written)") in desc.splitlines(), desc
    with patched(failing_call(owner, name, 1, RuntimeError('injected'), after=False)):
      with Optimizers.his_charges_strict():
        try:
          run(pdb_1exr, "tst_reduce2_his_charges_opt_ws" + tag, "qmready=True")
        except RuntimeError as e:
          assert str(e) == 'injected', e
        else:
          raise AssertionError("strict mode did not raise")

def test_memory_error_propagates():
  """A MemoryError in a rerun is re-raised as it is, not handled."""
  with patched(failing_call(Optimizers.Optimizer, '_OptimizeModel', 2, MemoryError('injected'))):
    try:
      run(pdb_1exr, "tst_reduce2_his_charges_opt_m", ON)
    except MemoryError as e:
      assert str(e) == 'injected', e
    else:
      raise AssertionError("MemoryError not re-raised")

def test_decide_error_second_his():
  """decide() fails on the second of two His (A 107 decided, D 163 not): both
  are reported with their first charges as restored, no rerun, output equal
  to flag off."""
  pdb = two_candidates()
  off, desc_off, made = run(pdb, "tst_reduce2_his_charges_opt_d2")
  with patched(failing_call(HisCharges, 'decide', 2, ValueError('injected'), after=False)):
    on, desc, made = run(pdb, "tst_reduce2_his_charges_opt_d2on", ON)
  m = only(made)
  assert (m['code'], m['detail'], m['reruns']) == (HisCharges.INTERNAL, 'ValueError: injected', 0), m
  assert on == off, removed(off, on)
  for r in ('107', '163'):
    c = case(m, r)
    assert c['category'] == 'restored' and c['final'] == c['free'] and not c.get('fixed'), c
  assert case(m, '107')['first'].get('target') == 0 and case(m, '163')['first']['decision'] == 'not decided', m

def test_abort_in_rerun():
  """A user's Stop (Abort) in a rerun propagates, after the first result is
  restored: the model then equals the first result."""
  first = []
  fp = Optimizers.Optimizer._HisFingerprint
  def fingerprint(self, myModel):
    r = fp(self, myModel)
    if not first:
      first.append((self, myModel, r))    # the first call: the first result
    return r
  with patched((Optimizers.Optimizer, '_HisFingerprint', fingerprint),
               failing_call(Optimizers.Optimizer, '_OptimizeModel', 2, Abort('injected stop'))):
    try:
      run(pdb_1exr, "tst_reduce2_his_charges_opt_ab", ON)
    except Abort as e:
      assert str(e) == 'injected stop', e
    else:
      raise AssertionError("Abort did not propagate")
  opt, myModel, want = first[0]
  assert fp(opt, myModel) == want
  assert hd1_a_deleted(opt) == [], hd1_a_deleted(opt)
  # a restore that cannot be verified after the Stop: the Sorry names it
  orig_save = Optimizers.Optimizer._HisSave
  def save(self, myModel):
    s = orig_save(self, myModel)
    s['deleted'] = []
    return s
  with patched((Optimizers.Optimizer, '_HisSave', save),
               failing_call(Optimizers.Optimizer, '_OptimizeModel', 2, Abort('injected stop'))):
    try:
      run(pdb_1exr, "tst_reduce2_his_charges_opt_ab2", ON)
    except Sorry as e:
      assert not isinstance(e, Abort), e
      assert str(e) == ("His charges: the restored first result differs from it (deleted); "
        "restoring after Abort: injected stop"), e
    else:
      raise AssertionError("unverified restore after the Stop not reported")

def test_diagnostic_errors():
  """The metal diagnostic for a His without a Mover (approach=optimize on the
  flag-off add output: B has one ring H) ignores its own errors, so the Mover's
  failure stands, but lets MemoryError and Abort through. The error is raised
  in the diagnostic's spatial-query call (a Python wrapper set only while
  _PlaceMovers runs)."""
  run(pdb_1exr, "tst_reduce2_his_charges_opt_g")
  with open("tst_reduce2_his_charges_opt_g_out.pdb") as f:
    added = f.read()
  base, desc, made = run(added, "tst_reduce2_his_charges_opt_g0", "approach=optimize", ON)
  place = Optimizers.Optimizer._PlaceMovers
  for tag, exc, raised in (('r', RuntimeError('injected'), False), ('m', MemoryError('injected'), True),
      ('a', Abort('injected'), True)):
    calls = []
    class Query(object):
      def __init__(self, real):
        self.real = real
      def __getattr__(self, n):
        return getattr(self.real, n)
      def neighbors(self, xyz, low, high):
        if (low, high) == (0.0, 3.0):
          calls.append(1)
          raise exc
        return self.real.neighbors(xyz, low, high)
    def wrapped(self, *a, **k):
      self._spatialQuery = Query(self._spatialQuery)
      try:
        return place(self, *a, **k)
      finally:
        if isinstance(self._spatialQuery, Query):
          self._spatialQuery = self._spatialQuery.real
    with patched((Optimizers.Optimizer, '_PlaceMovers', wrapped)):
      try:
        atoms, desc, made = run(added, "tst_reduce2_his_charges_opt_g" + tag, "approach=optimize", ON)
      except (MemoryError, Abort) as e:
        assert raised and e is exc, (tag, e)
      else:
        assert not raised, tag
        assert atoms == base, removed(base, atoms)
        c = case(only(made), '107')
        assert c['decision'].startswith('unchanged: unscorable conformer B (failed: MoverHisFlip()'), c
    assert calls, tag

def ring_xyz(atoms, alt):
  return dict((l[0:4], l[18:42]) for l in atoms if l[9] == 'A' and l[10:14].strip() == '107'
    and l[0:4] in (' CG ', ' ND1', ' CD2', ' CE1', ' NE2') and l[4] == alt)

def test_set_flip_states():
  """With flips, conformer A of His A 107 locked by set_flip_states. Unflipped:
  A takes both H (B one), reconciled to 0 within the locked orientation (ring
  of A as in the input; HisPlace line, pose word Unflipped). Flipped: A then
  takes HE2 only, both conformers are neutral, no case; output as flag off."""
  inp = ring_xyz([l[12:80] for l in pdb_1exr.splitlines() if l.startswith("ATOM")], 'A')
  lock = "set_flip_states=1 a A HIS 107 Unflipped"
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_14", "add_flip_movers=True", lock)
  assert ring_h(off, 'A', '107') == [('HD1', 'A'), ('HD1', 'B'), ('HE2', 'A')], ring_h(off, 'A', '107')
  on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_14on", "add_flip_movers=True", lock, ON)
  m = only(made)
  assert m['code'] == HisCharges.RECONCILED, m
  c = case(m, '107')
  assert c['final'] == {'A': 0, 'B': 0} and c.get('fixed'), c
  assert ring_h(on, 'A', '107') == [('HD1', 'B'), ('HE2', 'A')], ring_h(on, 'A', '107')
  assert ring_xyz(on, 'A') == inp, (ring_xyz(on, 'A'), inp)
  got = [x for x in reduce2._FindFlipsInOutputString(desc, 'HisPlace') if x.resName == 'AHIS']
  assert len(got) == 1 and not got[0].flipped, got
  lock = "set_flip_states=1 a A HIS 107 Flipped AnglesAdjusted"
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_15", "add_flip_movers=True", lock)
  on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_15on", "add_flip_movers=True", lock, ON)
  assert only(made)['cases'] == [] and on == off, removed(off, on)
  assert ring_xyz(on, 'A') != inp
  got = [x for x in reduce2._FindFlipsInOutputString(desc, 'HisPlace') if x.resName == 'AHIS']
  assert len(got) == 1 and got[0].flipped, got

def test_optimize():
  """approach=optimize. On reduce2's flag-off add output, B has one ring H, so
  no His Mover for B: unchanged with the reason, output as flag off. On the
  approach=place output (both ring H on both conformers): the decision of add,
  charge 0."""
  run(pdb_1exr, "tst_reduce2_his_charges_opt_16")
  with open("tst_reduce2_his_charges_opt_16_out.pdb") as f:
    added = f.read()
  off, desc_off, made = run(added, "tst_reduce2_his_charges_opt_16o", "approach=optimize")
  on, desc, made = run(added, "tst_reduce2_his_charges_opt_16on", "approach=optimize", ON)
  c = case(only(made), '107')
  assert c['category'] == 'unchanged' and c['decision'].startswith('unchanged: unscorable conformer B'), c
  assert on == off, removed(off, on)
  run(pdb_1exr, "tst_reduce2_his_charges_opt_17", "approach=place")
  with open("tst_reduce2_his_charges_opt_17_out.pdb") as f:
    placed = f.read()
  on, desc, made = run(placed, "tst_reduce2_his_charges_opt_17on", "approach=optimize", ON)
  m = only(made)
  c = case(m, '107')
  assert m['code'] == HisCharges.RECONCILED and c['free'] == {'A': 1, 'B': 0}, m
  assert c['final'] == {'A': 0, 'B': 0} and c.get('fixed'), c
  assert ring_h(on, 'A', '107') == [('HD1', 'B'), ('HE2', 'A')], ring_h(on, 'A', '107')

def test_real_deadline():
  """treat_cliques_as_singletons_after_seconds=0: the first run reaches its
  deadline, no rerun, output as flag off with the same setting."""
  dl = "treat_cliques_as_singletons_after_seconds=0"
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_18", dl)
  on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_18on", dl, ON)
  m = only(made)
  assert (m['code'], m['detail'], m['reruns']) == (HisCharges.SEARCH_DEADLINE, 'the first run reached it', 0), m
  assert not case(m, '107').get('fixed')
  assert on == off, removed(off, on)

def two_candidates():
  a = iotbx.pdb.input(source_info=None, lines=pdb_1exr).construct_hierarchy()
  b = iotbx.pdb.input(source_info=None, lines=pdb_1us0_163).construct_hierarchy()
  for c in b.chains():
    c.id = dict(A='D', B='E')[c.id]
  for x in b.atoms():
    x.xyz = (x.xyz[0] + 60, x.xyz[1], x.xyz[2])
  for c in b.models()[0].chains():
    a.models()[0].append_chain(c.detached_copy())
  a.atoms().reset_i_seq()
  return a.as_pdb_string()

D163 = (0, 'D', '163', '')

def loss_on_d163():
  orig = HisCharges.trial
  def trial(before, after, his, *a, **k):
    if D163 in his:
      after = [(alt, dict((m, s - 1.0 if m in his[D163] else s) for m, s in sc.items())) for alt, sc in after]
    return orig(before, after, his, *a, **k)
  return (HisCharges, 'trial', trial)

def test_two_candidates():
  """1exr His A 107 and 1us0 His A 163 (as chain D, 60 A away) in one model: one
  joint rerun fixes both. With a loss of 1.0 injected on His D 163 whenever it is
  fixed, the joint trial fails, His A 107 is accepted alone (confirmed fresh) and
  His D 163 is released. The rerun cap: with the cap at the reruns this plan
  needs it ends DONE; one less restores (run cap)."""
  pdb = two_candidates()
  off, desc_off, made = run(pdb, "tst_reduce2_his_charges_opt_6")
  on, desc, made = run(pdb, "tst_reduce2_his_charges_opt_6on", ON)
  m = only(made)
  assert m['code'] == HisCharges.RECONCILED and m['reruns'] == 1 and m['cap'] == 19, m
  assert case(m, '107')['final'] == {'A': 0, 'B': 0} and case(m, '163')['final'] == {'A': 1, 'B': 1}, m
  gone, new = removed(off, on)
  assert [g[:15] for g in gone] == [" HD1AHIS A 107 "] and [n[:15] for n in new] == [" HE2BHIS D 163 "], (gone, new)
  with patched(loss_on_d163()):
    on, desc, made = run(pdb, "tst_reduce2_his_charges_opt_6i", ON)
  m = only(made)
  k = m['reruns']
  assert m['code'] == HisCharges.RECONCILED and k > 1, m
  assert case(m, '107').get('fixed') and case(m, '107')['final'] == {'A': 0, 'B': 0}, m
  c = case(m, '163')
  assert not c.get('fixed') and c['decision'].startswith('kept: rerun failed') and c['category'] == 'rejected', c
  assert c['final'] == {'A': 1, 'B': 0}, c
  assert list(m['rejected']) == [D163], m['rejected']
  gone, new = removed(off, on)
  assert [g[:15] for g in gone] == [" HD1AHIS A 107 "] and new == [], (gone, new)
  # the cap boundary
  with patched(loss_on_d163(), (HisCharges, 'RUN_CAP_LIMIT', k)):
    on2, desc, made = run(pdb, "tst_reduce2_his_charges_opt_6k", ON)
  m = only(made)
  assert (m['code'], m['reruns'], m['cap']) == (HisCharges.RECONCILED, k, k), m
  assert on2 == on
  with patched(loss_on_d163(), (HisCharges, 'RUN_CAP_LIMIT', k - 1)):
    on2, desc, made = run(pdb, "tst_reduce2_his_charges_opt_6c", ON)
  m = only(made)
  assert (m['code'], m['reruns'], m['cap']) == (HisCharges.RUN_CAP, k - 1, k - 1), m
  assert m['reason'] == 'first result restored: rerun cap reached (%d rerun(s))' % (k - 1), m
  assert on2 == off, removed(off, on2)

def fixed_keys(made, pdb):
  class P(object):
    pass
  p = P()
  p.model = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None, lines=pdb))
  p._hisCharges = made[0].getHisCharges()
  return reduce2.Program._HisFixedKeys(p)

def test_failure_after_first_his():
  """Two accepted His; the written-state check fails for the
  second (His D 163). Nothing is marked: the written model is the free one, the
  header counts none set, each His line shows its free charges and the restore,
  no suggested flip is marked, and _HisFixedKeys() is empty."""
  pdb = two_candidates()
  off, desc_off, made = run(pdb, "tst_reduce2_his_charges_opt_m1")
  with patched(failing_check(on_call=2)):
    on, desc, made = run(pdb, "tst_reduce2_his_charges_opt_m1on", ON)
  m = only(made)
  assert m['code'] == HisCharges.INTERNAL and m['detail'] == 'His D 163 written state: injected', m
  assert on == off, removed(off, on)
  for c in m['cases']:
    assert c['final'] == c['free'] and c['category'] == 'restored', c
    assert not [x for x in ('fixed', 'checked', 'skipped') if x in c], c
  assert "2 His with different charges, 0 set to one charge" in desc, desc
  for r in ('A 107: A +1 B 0; kept: first result kept: internal error',
            'D 163: A +1 B 0; kept: first result kept: internal error'):
    assert "  His " + r in desc, desc
  assert "set to charge" not in desc and "(state fixed by charge reconciliation)" not in desc, desc
  assert fixed_keys(made, pdb) == set()

def test_blank_ring_rule_in_driver():
  """The blank-ring rule through its caller: with the helper reporting pass A as
  the last writer without the ring H, His A 107 is left unchanged, no rerun."""
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_r1")
  with patched((Optimizers, '_HisBlankRingLastPassLacksH', lambda recs, key: 'A')):
    on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_r1on", ON)
  m = only(made)
  assert m['code'] == HisCharges.NO_RERUN and m['reruns'] == 0, m
  c = case(m, '107')
  assert c['decision'] == 'unchanged: blank ring last written by pass A, which lacks the ring H', c
  assert c['category'] == 'unchanged' and c['final'] == {'A': 1, 'B': 0}, c
  assert on == off, removed(off, on)
  assert "  model: no rerun; first run " in desc, desc

def test_check_written_marks():
  """_HisCheckWritten on plain records: the accepted His is marked fixed with the
  written charges; a His that was not accepted gets its written charges too; a His
  that was no case but changed in the written run is added as 'changed'; with an
  error on any accepted His nothing is marked."""
  k1, k2, k3 = (0, 'A', '1', ''), (0, 'A', '2', ''), (0, 'A', '3', '')
  st = dict(accepted={k1})
  state = dict(dels={}, coarse={}, fixup={})
  freeRes = {k1: ({'A': ({'states': [state, state]}, {}), 'B': ({'states': [state, state]}, {})}, {'A': 1, 'B': 0}),
    k2: ({}, {'A': 0, 'B': 0}), k3: ({}, {'A': 1, 'B': 0})}
  now = {k1: ({'A': ({}, {}), 'B': ({}, {})}, {'A': 0, 'B': 0}), k2: ({}, {'A': 1, 'B': 1}),
    k3: ({}, {'A': 1, 'B': 1})}
  pins = {(k1, 'A'): 1, (k1, 'B'): 1}
  def cases():
    return [dict(key=k1, free={'A': 1, 'B': 0}, final={'A': 1, 'B': 0}, first=dict(kind='reconciled')),
      dict(key=k3, free={'A': 1, 'B': 0}, final={'A': 1, 'B': 0}, first=dict(kind='kept', target=1))]
  o = Optimizers.Optimizer.__new__(Optimizers.Optimizer)
  o._skipBondFixup = False
  o._HisWritten = lambda his, k, confs: {}
  ok = lambda fixed, written, skip: dict(errors=[], checked={'H': 1}, skipped=[])
  bad = lambda fixed, written, skip: dict(errors=['x'], checked={}, skipped=[])
  cs = cases()
  with patched((HisCharges, 'check_written', ok)):
    o._HisCheckWritten(st, None, now, pins, freeRes, cs)
  by = dict((c['key'], c) for c in cs)
  assert by[k1].get('fixed') and by[k1]['final'] == {'A': 0, 'B': 0}, by[k1]
  assert not by[k3].get('fixed') and by[k3]['final'] == {'A': 1, 'B': 1}, by[k3]
  assert by[k2]['first']['kind'] == 'changed' and by[k2]['final'] == {'A': 1, 'B': 1}, by[k2]
  assert by[k2]['free'] == {'A': 0, 'B': 0}, by[k2]
  cs = cases()
  with patched((HisCharges, 'check_written', bad)):
    try:
      o._HisCheckWritten(st, None, now, pins, freeRes, cs)
    except HisCharges.ContractError:
      pass
    else:
      raise AssertionError("no ContractError")
  assert cs == cases(), cs

def program_stub(his):
  class P(object):
    pass
  p = P()
  h = iotbx.pdb.hierarchy.root()
  h.append_model(iotbx.pdb.hierarchy.model())
  class M(object):
    def get_hierarchy(self):
      return h
  p.model = M()
  p._hisCharges = his
  p.params = P()
  p.params.qmready, p.params.reconcile_his_charges = True, True
  for n in ('_HisChargesText', '_HisChargesHint', '_WithHisChargesHint', '_HisFixedKeys', '_HisWithDifferentRingH'):
    setattr(p, n, getattr(reduce2.Program, n).__get__(p))
  return p

def test_hint_and_report():
  """The qmready hint and the report: classified from the structure (category,
  outcome code), never from the text; read from final, not free; attached only to
  the conformer-charge Sorry."""
  def his(code, detail, *cases):
    return [dict(model=0, code=code, detail=detail, cases=list(cases))]
  def c(n, free, final, category, decision):
    return dict(key=(0, 'A', str(n), ''), free=free, final=final, category=category, decision=decision,
      first=dict(kind='x'))
  two, one = {'A': 1, 'B': 0}, {'A': 1, 'B': 1}
  # internal failures whose text contains 'ion'
  for detail in ("('evaluation of another run', [], [])",
                 "His A 1 written state: His conformer A NE2: no expected position"):
    p = program_stub(his(HisCharges.INTERNAL, detail, c(1, two, two, 'restored',
      'kept: first result kept: internal error (%s)' % detail)))
    h = p._HisChargesHint()
    assert "expected)" not in h and "kept after an internal error" in h and "please report" in h, h
  p = program_stub(his(HisCharges.NO_RERUN, None,
    c(1, two, two, 'metal', 'kept: B metal-bound NE2 without H (ZN)'),
    c(2, two, two, 'cost', 'kept: cost to 0 is 0.80 > 0.50'),
    c(3, two, one, 'cost', 'kept: cost to +1 is 0.90 > 0.50'),
    c(4, one, two, 'changed', 'changed in the rerun (not fixed)')))
  h = p._HisChargesHint().splitlines()
  assert h == ["His A 1: conformer charges differ (B metal-bound NE2 without H (ZN); expected)",
    "reduce2 kept His A 2 with different charges (kept: cost to 0 is 0.80 > 0.50); set a common "
    "state with protonation { selection = ... state = ... }; a broader alternative (changes all "
    "flip decisions): add_flip_movers=True",
    "His A 4: conformer charges differ (changed in the rerun (not fixed))"], h
  # only the marked conformer-charge Sorry gets the hint; others pass unchanged
  def raise_(e):
    raise e
  marked = Sorry('charges')
  marked.conformer_charge_mismatch = True
  for e, hinted in ((Sorry('unknown residue XYZ'), False), (marked, True)):
    try:
      p._WithHisChargesHint(lambda: raise_(e))
    except Sorry as r:
      assert type(r) is Sorry, type(r)
      assert ("His A 1: conformer charges differ" in str(r)) == hinted, (e, r)
      assert str(r).startswith(str(e)), r
      assert (r is e) == (not hinted), (e, r)
      assert getattr(r, 'conformer_charge_mismatch', False) == hinted, r
  # the report: free -> final
  m = his(HisCharges.NO_RERUN, None, c(3, two, one, 'cost', 'kept: cost to +1 is 0.90 > 0.50'))
  m[0].update(reason='no rerun', freeSeconds=0.5, rejected={})
  text = program_stub(m)._HisChargesText(m, [])
  assert "  His A 3: A +1 B 0 -> A +1 B +1; kept: cost to +1 is 0.90 > 0.50" in text, text
  assert "  model: no rerun; first run 0.50 s" in text, text

def test_blank_ring_split_h():
  """1us0 His A 83: blank ring, ring H split A/B. Both conformers end +1, so
  there is no case: this checks only that the flag leaves this layout unchanged
  (a blank ring with split H and different charges goes through decide() in the
  module tests; a real blank-ring-last-pass fixture is still to come)."""
  off, desc_off, made = run(pdb_1us0_83, "tst_reduce2_his_charges_opt_7")
  on, desc, made = run(pdb_1us0_83, "tst_reduce2_his_charges_opt_7on", ON)
  m = only(made)
  assert m['code'] == HisCharges.NO_RERUN and m['cases'] == [] and m['outcome'] is None, m
  assert ring_h(on, 'A', '83') == [('HD1', 'A'), ('HD1', 'B'), ('HE2', 'A'), ('HE2', 'B')], ring_h(on, 'A', '83')
  assert on == off, removed(off, on)
  assert "  model: no rerun; first run " in desc, desc

def with_zn(alt, n='NE2'):
  from scitbx import matrix
  h = iotbx.pdb.input(source_info=None, lines=pdb_1exr).construct_hierarchy()
  at = dict(((a.name.strip(), a.parent().altloc), matrix.col(a.xyz)) for a in h.atoms()
    if a.parent().parent().resseq.strip() == '107')
  x = at[(n, alt)]
  ring = dict(NE2=('CD2', 'CE1'), ND1=('CG', 'CE1'))[n]
  zn = x + (x - (at[(ring[0], alt)] + at[(ring[1], alt)]) * 0.5).normalize() * 2.05
  return pdb_1exr + "HETATM  999 ZN   ZN%s Z   1    %8.3f%8.3f%8.3f  0.50 20.00          ZN\n" % (
    alt, zn[0], zn[1], zn[2])

def test_metal_kept():
  """A Zn at NE2 of conformer B only: B loses HE2, the His keeps A +1, B 0 with
  the reason; no rerun. qmready then stops with the hint. A Zn at conformer B's
  ND1 is labelled metal-bound ND1 (flag on only; here both conformers end at 0,
  so the label is checked where it is recorded)."""
  pdb = with_zn('B')
  off, desc_off, made = run(pdb, "tst_reduce2_his_charges_opt_8")
  on, desc, made = run(pdb, "tst_reduce2_his_charges_opt_8on", ON)
  m = only(made)
  assert m['outcome'] is None and m['reruns'] == 0 and m['code'] == HisCharges.NO_RERUN, m
  c = case(m, '107')
  assert c['decision'] == 'kept: B metal-bound NE2 without H (ZN)' and c['final'] == {'A': 1, 'B': 0}, c
  assert c['category'] == 'metal', c
  assert on == off, removed(off, on)
  try:
    run(pdb, "tst_reduce2_his_charges_opt_8q", ON, "qmready=True")
  except Sorry as e:
    assert "His A 107: conformer charges differ (B metal-bound NE2 without H (ZN); expected)" in str(e), e
  else:
    raise AssertionError("qmready ran with different conformer charges")
  pdb = with_zn('B', 'ND1')
  seen = []
  orig = Optimizers.Optimizer._HisRecordPass
  def rec(self, mi, alt, *a, **k):
    seen.extend(self._hisNoMover.values())
    return orig(self, mi, alt, *a, **k)
  with patched((Optimizers.Optimizer, '_HisRecordPass', rec)):
    on, desc, made = run(pdb, "tst_reduce2_his_charges_opt_8d", ON)
  assert 'metal-bound ND1 without H (ZN)' in seen, seen
  off, desc, made = run(pdb, "tst_reduce2_his_charges_opt_8doff")
  assert not [w for w in made[0]._hisNoMover.values() if w.startswith(('metal-bound', 'failed'))], made[0]._hisNoMover
  assert on == off, removed(off, on)

def test_multi_model():
  """Two models: each reconciled, keys and report per model; the reruns' time is
  added to the second model's search deadline."""
  h = iotbx.pdb.input(source_info=None, lines=pdb_1exr).construct_hierarchy()
  m2 = h.models()[0].detached_copy()
  m2.id, h.models()[0].id = '2', '1'
  h.append_model(m2)
  seen = []
  orig = Optimizers.Optimizer._ReconcileHisCharges
  def rec(self, mi, myModel, ctx):
    r = orig(self, mi, myModel, ctx)
    seen.append((ctx['stopOptimizingAt'], r))
    return r
  with patched((Optimizers.Optimizer, '_ReconcileHisCharges', rec)):
    on, desc, made = run(h.as_pdb_string(), "tst_reduce2_his_charges_opt_9", ON)
  his = made[0].getHisCharges()
  assert [m['model'] for m in his] == [0, 1] and all(m['code'] == HisCharges.RECONCILED for m in his), his
  assert [c['key'] for m in his for c in m['cases']] == [(0, 'A', '107', ''), (1, 'A', '107', '')], his
  assert len(seen) == 2 and seen[0][1] > 0 and seen[1][0] == seen[0][0] + seen[0][1], seen
  for n in (1, 2):
    assert "  model %d His A 107: A +1 B 0 -> A 0 B 0; set to charge 0" % n in desc, desc
    assert "  model %d reconciled; first run " % n in desc, desc

def test_qmready():
  """qmready on the 1exr fixture: stops with the conformer-charge Sorry (a plain
  Sorry marked conformer_charge_mismatch, which the hint is attached to) with the
  flag off, runs (total charge 0) with it on. get_total_charge() itself raises a
  plain Sorry, as every other caller sees it."""
  from mmtbx.ligands import electrons
  seen = []
  orig = electrons.electron_distribution.get_total_charge
  def total(self):
    try:
      return orig(self)
    except Exception as e:
      seen.append(e)
      raise
  try:
    with patched((electrons.electron_distribution, 'get_total_charge', total)):
      run(pdb_1exr, "tst_reduce2_his_charges_opt_10", "qmready=True", "reconcile_his_charges=False")
  except Sorry as e:
    assert type(e) is Sorry and e.conformer_charge_mismatch, type(e)
    assert "different total charges (A: 1, B: 0)" in str(e), e
  else:
    raise AssertionError("qmready ran without reconciliation")
  assert len(seen) == 1 and type(seen[0]) is Sorry and seen[0].conformer_charge_mismatch, seen
  on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_10on", ON, "qmready=True")
  assert only(made)['code'] == HisCharges.RECONCILED

def run_command(pdb, prefix, *extra):
  """mmtbx.reduce2 in a subprocess: (return code, stderr lines)."""
  from libtbx import easy_run
  with open(prefix + ".pdb", "w") as f:
    f.write(pdb)
  r = easy_run.fully_buffered(" ".join(["mmtbx.reduce2", prefix + ".pdb",
    "output.filename=%s_out.pdb" % prefix, "overwrite=True"] + list(extra)))
  return r.return_code, [l for l in r.stderr_lines if l.strip()]

def test_qmready_stderr():
  """What the user sees when qmready stops on conformer charges: flag off
  (explicit False; Auto is on with qmready), master's text (one Sorry: message,
  no class name, no traceback, exit 1) and the hint to turn it on; flag on with
  a kept His (Zn at conformer B's NE2), the same message followed by the hint,
  no "During handling"."""
  first = ("Sorry: The alternate conformations have different total charges "
    "(A: %d, B: %d). The charge of the model is the charge of one conformer; "
    "check the protonation of the residues with altlocs.")
  rows = "Residues whose charge differs between the conformers: A 107 (A: 1, B: 0)"
  hint = "His A 107: conformer charges differ (B metal-bound NE2 without H (ZN); expected)"
  turn_on = ("His A 107: conformer charges differ; reconcile_his_charges=True (the default with "
    "qmready=True) gives a His one charge where the scores allow it")
  for pdb, prefix, extra, want in (
      (pdb_1exr, "tst_reduce2_his_charges_opt_12", ("reconcile_his_charges=False",), [first % (1, 0), rows, turn_on]),
      (with_zn('B'), "tst_reduce2_his_charges_opt_12on", (ON,), [first % (3, 2), rows, hint])):
    rc, err = run_command(pdb, prefix, "qmready=True", *extra)
    assert rc == 1, (rc, err)
    assert err == want, err

def test_auto():
  """reconcile_his_charges=Auto (the default): on with qmready, off otherwise.
  1exr fixture: Auto with qmready reconciles (total charge 0, no stop); Auto
  without qmready writes what explicit False writes; explicit False with
  qmready stops with the hint to turn it on; explicit True without qmready
  reconciles."""
  on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_a1", "qmready=True")
  assert [o._reconcileHisCharges for o in made] == [True], made
  assert only(made)['code'] == HisCharges.RECONCILED
  assert "qmready=True: reconcile_his_charges set to True" in LOG[0].getvalue()
  assert "Total charge: 0" in LOG[0].getvalue()
  auto, desc_auto, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_a2")
  assert [(o._reconcileHisCharges, o.getHisCharges()) for o in made] == [(False, [])], made
  assert "reconcile_his_charges set to" not in LOG[0].getvalue()
  off, desc_off, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_a3", "reconcile_his_charges=False")
  assert auto == off
  mask = lambda d: [re.sub(r" in [0-9.]+ seconds", "", l) for l in d.splitlines()
    if not l.startswith("reduce2 v.") and "Time to " not in l]
  assert mask(desc_auto) == mask(desc_off)
  try:
    run(pdb_1exr, "tst_reduce2_his_charges_opt_a4", "qmready=True", "reconcile_his_charges=False")
  except Sorry as e:
    assert str(e).splitlines()[-1] == ("His A 107: conformer charges differ; reconcile_his_charges=True "
      "(the default with qmready=True) gives a His one charge where the scores allow it"), e
  else:
    raise AssertionError("qmready with reconcile_his_charges=False ran")
  on, desc, made = run(pdb_1exr, "tst_reduce2_his_charges_opt_a5", ON)
  assert [o._reconcileHisCharges for o in made] == [True], made
  assert only(made)['code'] == HisCharges.RECONCILED
  # an explicit value is not announced
  assert "reconcile_his_charges set to" not in LOG[0].getvalue()

def test_flip_report_words():
  """With flips, a fixed His keeps the Mover report format that reduce2's
  _FindFlipsInOutputString reads by word position (pose at word 13, the flip
  state at 14), so flipkin and set_flip_states read it right; the fixed state
  is noted at the end of the line. 1us0 His A 163, both conformers fixed."""
  on, desc, made = run(pdb_1us0_163, "tst_reduce2_his_charges_opt_13", ON, "add_flip_movers=True")
  assert only(made)['code'] == HisCharges.RECONCILED
  lines = [l for l in desc.splitlines() if l.split()[:1] == ['HisFlip'] and 'Initial score:' in l]
  fixed = [l for l in lines if l.endswith(' (state 0 fixed)')]
  assert len(fixed) == 2 and all('HIS 163' in l for l in fixed), lines
  for l in lines:
    w = l.split()
    assert w[13] == 'pose' and w[14] in ('Flipped', 'Unflipped'), l
  got = reduce2._FindFlipsInOutputString(desc, 'HisFlip')
  assert len(got) == len(lines) and [s.flipped for s in got] == [l.split()[14] == 'Flipped' for l in lines], got

def test_search_atom_info():
  """At the end of the search (in _HisUnpin, before the full Movers return),
  the atom info of every atom of a fixed His (ring N acceptors, donors, radii)
  is that of the fixed state, the one-state Mover's only extras: A state 2 (HE2
  only; in the free run A had both H), so ND1 accepts; B state 1 (HD1 only), so
  ND1 does not. A one-state Mover can apply no other extras during the search."""
  seen = []
  orig = Optimizers.Optimizer._HisUnpin
  def unpin(self, his):
    for m, pin in his['pinned'].items():
      assert len(m._extras) == 1, len(m._extras)
      for i, a in enumerate(m.CoarsePositions().atoms):
        e, f = self._extraAtomInfo.getMappingFor(a), m._extras[0][i]
        assert (e.isAcceptor, e.isDonor, e.vdwRadius) == (f.isAcceptor, f.isDonor, f.vdwRadius), (a.id_str(), pin)
      seen.append((m.CoarsePositions().atoms[4].parent().altloc, pin,
        self._extraAtomInfo.getMappingFor(m.CoarsePositions().atoms[4]).isAcceptor))
    return orig(self, his)
  with patched((Optimizers.Optimizer, '_HisUnpin', unpin)):
    run(pdb_1exr, "tst_reduce2_his_charges_opt_11", ON)
  assert sorted(seen) == [('A', 2, True), ('B', 1, False)], seen

def test_blank_ring_last_pass():
  """Ring H in B and C, blank ring, the model has A: passes C,
  B, A; pass A sees the bare ring and writes it last. The His is left unchanged
  (_HisBlankRingLastPassLacksH), and _HisWritten names pass A as the writer."""
  k = (0, 'A', '5', '')
  recs = [dict(key=k, passAlt='C', conf='C'), dict(key=k, passAlt='B', conf='B'),
    dict(key=k, passAlt='A', conf='', hNotInPass=True)]
  assert Optimizers._HisBlankRingLastPassLacksH(recs, k) == 'A'
  assert Optimizers._HisBlankRingLastPassLacksH(recs[:2], k) is None
  assert Optimizers._HisBlankRingLastPassLacksH(recs[:1] + recs[2:] + recs[1:2], k) is None
  assert Optimizers._HisBlankRingLastPassLacksH(recs, (0, 'A', '6', '')) is None
  class A(object):
    def __init__(self, i):
      self.i_seq, self.xyz = i, (0.0, 0.0, float(i))
  class Info(object):
    isAcceptor = True
  class Ex(object):
    def getMappingFor(self, a):
      return Info()
  ring = dict(NE2=A(1), ND1=A(2), CG=A(3), CD2=A(4), CE1=A(5))
  ringB = dict(ring, HE2=A(6))
  ringC = dict(ring, HD1=A(7))
  o = Optimizers.Optimizer.__new__(Optimizers.Optimizer)
  o._hydrogensToDelete, o._extraAtomInfo = [], Ex()
  his = dict(recs=[(dict(recs[0], ring=ringC), {}), (dict(recs[1], ring=ringB), {}),
    (dict(recs[2], ring=ring), {})])
  w = o._HisWritten(his, k, {'B': (dict(recs[1], ring=ringB), {}), 'C': (dict(recs[0], ring=ringC), {})})
  assert w['B']['writer']['NE2'] == 'pass A without the ring H' and w['B']['writer']['HE2'] == 'B', w


pdb_1exr = """\
ATOM    935  N   ARG A  94      34.136 -16.706  29.792  1.00 13.44           N
ATOM    936  CA  ARG A  94      35.261 -17.201  30.584  1.00 15.05           C
ATOM    937  C   ARG A  94      35.967 -18.368  29.930  1.00 15.88           C
ATOM    938  O   ARG A  94      37.204 -18.484  30.059  1.00 19.56           O
ATOM    939  CB  ARG A  94      34.779 -17.552  31.990  1.00 16.59           C
ATOM    940  CG  ARG A  94      34.216 -16.398  32.818  1.00 19.74           C
ATOM    941  CD  ARG A  94      33.862 -16.888  34.231  1.00 21.05           C
ATOM    942  NE  ARG A  94      35.044 -17.164  35.019  1.00 23.55           N
ATOM    943  CZ  ARG A  94      35.279 -18.159  35.849  1.00 24.35           C
ATOM    944  NH1 ARG A  94      34.386 -19.122  36.089  1.00 39.46           N
ATOM    945  NH2 ARG A  94      36.442 -18.178  36.481  1.00 31.24           N
ATOM   1005  N   ALA A 103      26.078 -19.251  33.669  1.00 11.12           N
ATOM   1006  CA  ALA A 103      27.044 -19.412  34.744  1.00 11.64           C
ATOM   1007  C   ALA A 103      28.186 -18.414  34.623  1.00  9.77           C
ATOM   1008  O   ALA A 103      28.620 -17.809  35.617  1.00 10.63           O
ATOM   1009  CB  ALA A 103      27.572 -20.847  34.756  1.00 14.26           C
ATOM   1010  N   GLU A 104      28.647 -18.214  33.408  1.00  9.68           N
ATOM   1011  CA  GLU A 104      29.733 -17.224  33.187  1.00  9.38           C
ATOM   1012  C   GLU A 104      29.282 -15.833  33.521  1.00  9.15           C
ATOM   1013  O   GLU A 104      29.995 -15.072  34.194  1.00 10.22           O
ATOM   1014  CB  GLU A 104      30.184 -17.298  31.725  1.00  9.64           C
ATOM   1015  CG  GLU A 104      30.993 -18.518  31.406  1.00 11.21           C
ATOM   1016  CD  GLU A 104      31.281 -18.734  29.941  1.00 12.08           C
ATOM   1017  OE1 GLU A 104      30.537 -18.227  29.084  1.00 11.43           O
ATOM   1018  OE2 GLU A 104      32.279 -19.456  29.640  1.00 14.06           O
ATOM   1027  N   ARG A 106      26.895 -14.883  35.508  1.00  9.22           N
ATOM   1028  CA  ARG A 106      26.719 -14.722  36.953  1.00  9.65           C
ATOM   1029  C   ARG A 106      28.056 -14.566  37.655  1.00  9.45           C
ATOM   1030  O   ARG A 106      28.210 -13.758  38.572  1.00 10.64           O
ATOM   1031  CB  ARG A 106      25.870 -15.881  37.488  1.00  9.92           C
ATOM   1032  CG  ARG A 106      25.440 -15.634  38.915  1.00 11.16           C
ATOM   1033  CD  ARG A 106      24.342 -16.590  39.346  1.00 11.19           C
ATOM   1034  NE  ARG A 106      24.571 -17.943  39.133  1.00 12.34           N
ATOM   1035  CZ  ARG A 106      25.209 -18.729  39.967  1.00 12.70           C
ATOM   1036  NH1 ARG A 106      25.783 -18.243  41.075  1.00 13.78           N
ATOM   1037  NH2 ARG A 106      25.363 -20.040  39.717  1.00 18.51           N
ATOM   1038  N  AHIS A 107      29.019 -15.368  37.207  0.34  8.71           N
ATOM   1039  CA AHIS A 107      30.341 -15.358  37.814  0.34 10.82           C
ATOM   1040  C  AHIS A 107      31.077 -14.031  37.693  0.34 12.08           C
ATOM   1041  O  AHIS A 107      31.706 -13.553  38.647  0.34 10.34           O
ATOM   1042  CB AHIS A 107      31.184 -16.457  37.142  0.34 11.82           C
ATOM   1043  CG AHIS A 107      32.441 -16.649  37.944  0.34 14.84           C
ATOM   1044  ND1AHIS A 107      33.504 -15.796  37.976  0.34 15.41           N
ATOM   1045  CD2AHIS A 107      32.745 -17.648  38.805  0.34 17.00           C
ATOM   1046  CE1AHIS A 107      34.432 -16.266  38.805  0.34 16.26           C
ATOM   1047  NE2AHIS A 107      33.988 -17.399  39.330  0.34 16.81           N
ATOM   1048  N  BHIS A 107      29.054 -15.396  37.312  0.66 10.87           N
ATOM   1049  CA BHIS A 107      30.352 -15.252  37.989  0.66 10.92           C
ATOM   1050  C  BHIS A 107      30.956 -13.870  37.775  0.66 11.22           C
ATOM   1051  O  BHIS A 107      31.380 -13.167  38.696  0.66 11.35           O
ATOM   1052  CB BHIS A 107      31.290 -16.355  37.492  0.66 11.21           C
ATOM   1053  CG BHIS A 107      32.739 -16.176  37.885  0.66 13.34           C
ATOM   1054  ND1BHIS A 107      33.159 -16.647  39.130  0.66 16.57           N
ATOM   1055  CD2BHIS A 107      33.826 -15.605  37.308  0.66 17.50           C
ATOM   1056  CE1BHIS A 107      34.453 -16.356  39.258  0.66 16.43           C
ATOM   1057  NE2BHIS A 107      34.900 -15.723  38.175  0.66 16.67           N
ATOM   1058  N   VAL A 108      30.977 -13.409  36.516  1.00 11.16           N
ATOM   1059  CA  VAL A 108      31.596 -12.123  36.256  1.00 11.53           C
ATOM   1060  C   VAL A 108      30.769 -11.023  36.930  1.00 11.56           C
ATOM   1061  O   VAL A 108      31.342 -10.062  37.513  1.00 12.76           O
ATOM   1062  CB  VAL A 108      31.797 -11.899  34.755  1.00 13.24           C
ATOM   1063  CG1 VAL A 108      32.393 -10.525  34.485  1.00 16.48           C
ATOM   1064  CG2 VAL A 108      32.707 -13.012  34.228  1.00 15.11           C
ATOM   1073  N   THR A 110      29.187 -11.175  39.754  1.00 11.08           N
ATOM   1074  CA  THR A 110      29.596 -11.157  41.171  1.00 11.80           C
ATOM   1075  C   THR A 110      30.881 -10.390  41.314  1.00 11.50           C
ATOM   1076  O   THR A 110      31.060  -9.583  42.258  1.00 12.59           O
ATOM   1077  CB  THR A 110      29.715 -12.619  41.665  1.00 14.17           C
ATOM   1078  OG1 THR A 110      28.415 -13.216  41.608  1.00 18.55           O
ATOM   1079  CG2 THR A 110      30.144 -12.649  43.113  1.00 21.50           C
ATOM   1080  N  AASN A 111      31.845 -10.634  40.457  0.72 11.72           N
ATOM   1081  CA AASN A 111      33.111  -9.834  40.563  0.72 10.91           C
ATOM   1082  C  AASN A 111      32.887  -8.368  40.312  0.72 10.52           C
ATOM   1083  O  AASN A 111      33.646  -7.556  40.907  0.72 12.14           O
ATOM   1084  CB AASN A 111      34.138 -10.489  39.672  0.72 12.73           C
ATOM   1085  CG AASN A 111      34.671 -11.785  40.267  0.72 15.99           C
ATOM   1086  OD1AASN A 111      35.190 -12.663  39.560  0.72 19.42           O
ATOM   1087  ND2AASN A 111      34.605 -12.004  41.573  0.72 15.05           N
ATOM   1088  N  BASN A 111      31.778 -10.642  40.375  0.28 11.41           N
ATOM   1089  CA BASN A 111      33.114 -10.070  40.212  0.28 11.93           C
ATOM   1090  C  BASN A 111      33.055  -8.580  39.879  0.28 11.88           C
ATOM   1091  O  BASN A 111      34.093  -7.890  39.997  0.28 15.67           O
ATOM   1092  CB BASN A 111      33.896 -10.846  39.145  0.28 15.73           C
ATOM   1093  CG BASN A 111      35.370 -10.896  39.477  0.28 19.65           C
ATOM   1094  OD1BASN A 111      35.791 -10.030  40.232  0.28 36.84           O
ATOM   1095  ND2BASN A 111      36.139 -11.837  38.963  0.28 29.48           N
HETATM 1504  O   HOH A2031      28.414 -19.030  38.108  1.00 21.58           O
HETATM 1527  O   HOH A2054      31.784 -19.349  35.227  1.00 24.96           O
HETATM 1538  O   HOH A2065      36.943 -15.021  40.078  1.00 30.53           O
HETATM 1542  O   HOH A2069      31.525 -17.270  41.211  1.00 31.52           O
HETATM 1565  O   HOH A2092      29.078 -17.394  40.211  1.00 24.67           O
HETATM 1575  O   HOH A2102      32.697 -14.477  41.218  1.00 30.96           O
HETATM 1588  O   HOH A2115      31.012 -20.152  38.014  1.00 34.73           O
TER
ATOM    518  N   GLY B  59      31.866 -15.926  47.231  1.00 14.69           N
ATOM    519  CA  GLY B  59      31.205 -16.582  46.117  1.00 16.07           C
ATOM    520  C   GLY B  59      32.022 -16.500  44.839  1.00 14.80           C
ATOM    521  O   GLY B  59      31.406 -16.714  43.776  1.00 17.34           O
ATOM    522  N   ASN B  60      33.309 -16.229  44.841  1.00 13.93           N
ATOM    523  CA  ASN B  60      34.088 -16.165  43.628  1.00 15.63           C
ATOM    524  C   ASN B  60      34.542 -17.532  43.211  1.00 14.22           C
ATOM    525  O   ASN B  60      35.190 -17.610  42.141  1.00 17.75           O
ATOM    526  CB  ASN B  60      35.226 -15.165  43.723  1.00 15.69           C
ATOM    527  CG  ASN B  60      36.309 -15.542  44.679  1.00  9.01           C
ATOM    528  OD1 ASN B  60      36.215 -16.541  45.441  1.00  9.98           O
ATOM    529  ND2 ASN B  60      37.331 -14.719  44.670  1.00 10.09           N
ATOM    530  N   GLY B  61      34.228 -18.609  43.959  1.00 11.70           N
ATOM    531  CA  GLY B  61      34.534 -19.950  43.508  1.00 11.82           C
ATOM    532  C   GLY B  61      35.771 -20.575  44.112  1.00  9.90           C
ATOM    533  O   GLY B  61      36.023 -21.737  43.892  1.00 11.11           O
HETATM 1633  O   HOH B2160      31.400 -20.960  42.613  1.00 38.85           O
TER
ATOM   1166  N  AASP C 119      42.595 -16.656  38.457  0.34  8.94           N
ATOM   1167  CA AASP C 119      41.394 -16.069  37.854  0.34  8.53           C
ATOM   1168  C  AASP C 119      41.665 -14.679  37.279  0.34  8.66           C
ATOM   1169  O  AASP C 119      41.105 -14.294  36.242  0.34  6.57           O
ATOM   1170  CB AASP C 119      40.242 -15.989  38.848  0.34  8.93           C
ATOM   1171  CG AASP C 119      39.452 -17.290  38.953  0.34  9.03           C
ATOM   1172  OD1AASP C 119      39.676 -18.205  38.150  0.34 10.97           O
ATOM   1173  OD2AASP C 119      38.611 -17.271  39.874  0.34 10.12           O
ATOM   1174  N  BASP C 119      42.571 -16.687  38.419  0.66  8.99           N
ATOM   1175  CA BASP C 119      41.443 -16.100  37.662  0.66 10.11           C
ATOM   1176  C  BASP C 119      41.852 -14.731  37.126  0.66  8.54           C
ATOM   1177  O  BASP C 119      41.558 -14.410  35.961  0.66  8.09           O
ATOM   1178  CB BASP C 119      40.195 -16.037  38.509  0.66 11.81           C
ATOM   1179  CG BASP C 119      38.948 -15.929  37.672  0.66 18.19           C
ATOM   1180  OD1BASP C 119      38.120 -15.097  37.998  0.66 32.94           O
ATOM   1181  OD2BASP C 119      38.767 -16.725  36.715  0.66 27.56           O
HETATM 1548  O   HOH C2075      36.863 -19.490  39.935  1.00 31.49           O
HETATM 1558  O   HOH C2085      39.161 -15.462  34.876  1.00 28.44           O
HETATM 1577  O   HOH C2104      36.493 -14.774  35.326  1.00 31.02           O
HETATM 1582  O   HOH C2109      35.170 -12.003  36.709  1.00 30.46           O
HETATM 1604  O   HOH C2131      39.303 -12.648  38.691  1.00 26.44           O
HETATM 1634  O   HOH C2161      37.476 -17.449  40.907  0.60 18.47           O
HETATM 1646  O   HOH C2173      36.745 -12.749  34.283  1.00 40.97           O
TER
"""

pdb_1us0_163 = """\
ATOM   1444  N   ASN A 162      15.491   6.955  11.774  1.00  4.01           N
ATOM   1445  CA  ASN A 162      16.514   7.956  11.555  1.00  4.29           C
ATOM   1446  C   ASN A 162      15.876   9.352  11.522  1.00  4.45           C
ATOM   1447  O   ASN A 162      14.666   9.517  11.615  1.00  4.83           O
ATOM   1448  CB  ASN A 162      17.352   7.610  10.335  1.00  4.77           C
ATOM   1449  CG  ASN A 162      16.603   7.743   9.029  1.00  5.60           C
ATOM   1450  OD1 ASN A 162      15.717   8.549   8.888  1.00  6.92           O
ATOM   1451  ND2 ASN A 162      17.012   6.946   8.041  1.00  6.01           N
ATOM   1452  N  AHIS A 163      16.783  10.331  11.429  0.60  4.57           N
ATOM   1453  CA AHIS A 163      16.333  11.697  11.471  0.60  5.57           C
ATOM   1454  C  AHIS A 163      15.236  11.990  10.437  0.60  5.73           C
ATOM   1455  O  AHIS A 163      14.211  12.606  10.704  0.60  6.47           O
ATOM   1456  CB AHIS A 163      17.576  12.593  11.290  0.60  6.31           C
ATOM   1457  CG AHIS A 163      17.372  14.039  11.401  0.60  7.51           C
ATOM   1458  ND1AHIS A 163      18.351  14.790  11.987  0.60 12.51           N
ATOM   1459  CD2AHIS A 163      16.461  14.915  11.072  0.60  8.19           C
ATOM   1460  CE1AHIS A 163      17.978  16.080  11.917  0.60 14.99           C
ATOM   1461  NE2AHIS A 163      16.822  16.141  11.289  0.60  9.46           N
ATOM   1462  N  BHIS A 163      16.653  10.488  11.436  0.40  7.46           N
ATOM   1463  CA BHIS A 163      16.005  11.780  11.701  0.40  7.75           C
ATOM   1464  C  BHIS A 163      15.042  12.105  10.533  0.40  8.52           C
ATOM   1465  O  BHIS A 163      14.036  12.715  10.882  0.40  8.72           O
ATOM   1466  CB BHIS A 163      16.986  12.891  11.987  0.40  9.69           C
ATOM   1467  CG BHIS A 163      17.482  13.617  10.803  0.40 12.73           C
ATOM   1468  ND1BHIS A 163      18.476  13.045  10.014  0.40 16.25           N
ATOM   1469  CD2BHIS A 163      17.249  14.846  10.325  0.40 16.59           C
ATOM   1470  CE1BHIS A 163      18.774  13.917   9.050  0.40 16.30           C
ATOM   1471  NE2BHIS A 163      18.056  15.017   9.206  0.40 17.83           N
ATOM   1472  N  ALEU A 164      15.497  11.558   9.215  0.60  5.79           N
ATOM   1473  CA ALEU A 164      14.485  11.866   8.181  0.60  7.69           C
ATOM   1474  C  ALEU A 164      13.189  11.126   8.408  0.60  5.71           C
ATOM   1475  O  ALEU A 164      12.108  11.628   8.035  0.60  6.21           O
ATOM   1476  CB ALEU A 164      15.130  11.571   6.834  0.60 10.89           C
ATOM   1477  CG ALEU A 164      16.275  12.542   6.551  0.60 14.93           C
ATOM   1478  CD1ALEU A 164      16.842  12.203   5.192  0.60 28.12           C
ATOM   1479  CD2ALEU A 164      15.917  14.038   6.554  0.60 23.68           C
ATOM   1480  N  BLEU A 164      15.340  11.718   9.302  0.40  7.83           N
ATOM   1481  CA BLEU A 164      14.466  11.842   8.126  0.40  7.56           C
ATOM   1482  C  BLEU A 164      13.107  11.214   8.470  0.40  6.53           C
ATOM   1483  O  BLEU A 164      12.056  11.781   8.195  0.40  7.00           O
ATOM   1484  CB BLEU A 164      15.055  11.262   6.839  0.40 11.14           C
ATOM   1485  CG BLEU A 164      14.246  11.145   5.556  0.40 15.07           C
ATOM   1486  CD1BLEU A 164      15.050  11.027   4.283  0.40 16.25           C
ATOM   1487  CD2BLEU A 164      13.400   9.806   5.554  0.40 13.62           C
ATOM   1488  N   GLN A 165      13.212   9.955   9.005  1.00  5.28           N
ATOM   1489  CA  GLN A 165      12.008   9.190   9.271  1.00  5.02           C
ATOM   1490  C   GLN A 165      11.208   9.799  10.438  1.00  5.16           C
ATOM   1491  O   GLN A 165       9.977   9.837  10.396  1.00  5.88           O
ATOM   1492  CB  GLN A 165      12.389   7.732   9.566  1.00  4.91           C
ATOM   1493  CG  GLN A 165      12.943   7.065   8.325  1.00  5.42           C
ATOM   1494  CD  GLN A 165      13.575   5.729   8.626  1.00  4.61           C
ATOM   1495  OE1 GLN A 165      14.140   5.478   9.685  1.00  5.06           O
ATOM   1496  NE2 GLN A 165      13.482   4.834   7.637  1.00  5.75           N
ATOM   1504  N   GLU A 167      11.286  12.975  11.154  1.00  6.08           N
ATOM   1505  CA  GLU A 167      10.736  14.158  10.506  1.00  6.83           C
ATOM   1506  C   GLU A 167       9.509  13.802   9.653  1.00  6.47           C
ATOM   1507  O   GLU A 167       8.526  14.576   9.649  1.00  7.28           O
ATOM   1508  CB  GLU A 167      11.836  14.797   9.685  1.00  7.87           C
ATOM   1509  CG  GLU A 167      11.398  16.100   9.067  1.00 10.85           C
ATOM   1510  CD  GLU A 167      12.591  16.674   8.344  1.00 15.72           C
ATOM   1511  OE1 GLU A 167      13.153  16.205   7.525  1.00 21.30           O
ATOM   1512  OE2 GLU A 167      12.959  17.862   8.907  1.00 21.59           O
ATOM   1782  CD  LYS A 194      25.051  14.502  13.003  0.70 10.95           C
ATOM   1783  CE  LYS A 194      26.444  14.923  12.742  0.70 19.43           C
ATOM   1784  NZ  LYS A 194      26.729  16.370  12.695  0.70 43.62           N
ATOM   1785  N  ALYS A 194      23.655  13.562  17.130  0.70  4.68           N
ATOM   1786  CA ALYS A 194      22.812  14.254  16.175  0.70  4.77           C
ATOM   1787  C  ALYS A 194      21.362  13.831  16.297  0.70  4.67           C
ATOM   1788  O  ALYS A 194      20.446  14.670  16.344  0.70  5.34           O
ATOM   1789  CB ALYS A 194      23.325  14.023  14.746  0.70  5.52           C
ATOM   1790  CG ALYS A 194      24.667  14.567  14.485  0.70  7.80           C
ATOM   1791  N  BLYS A 194      23.377  14.057  17.407  0.30  6.13           N
ATOM   1792  CA BLYS A 194      22.487  14.712  16.449  0.30  6.47           C
ATOM   1793  C  BLYS A 194      21.017  14.263  16.567  0.30  5.83           C
ATOM   1794  O  BLYS A 194      20.091  15.066  16.685  0.30  6.90           O
ATOM   1795  CB BLYS A 194      23.031  14.422  15.033  0.30  7.90           C
ATOM   1796  CG BLYS A 194      22.456  15.091  13.917  0.30 11.47           C
ATOM   1797  N  ALEU A 195      21.108  12.513  16.317  0.70  4.33           N
ATOM   1798  CA ALEU A 195      19.708  12.074  16.353  0.70  4.27           C
ATOM   1799  C  ALEU A 195      19.052  12.429  17.687  0.70  4.07           C
ATOM   1800  O  ALEU A 195      17.887  12.788  17.730  0.70  4.53           O
ATOM   1801  CB ALEU A 195      19.642  10.584  16.070  0.70  4.44           C
ATOM   1802  CG ALEU A 195      18.245   9.976  16.006  0.70  4.77           C
ATOM   1803  CD1ALEU A 195      17.336  10.679  15.013  0.70  5.85           C
ATOM   1804  CD2ALEU A 195      18.358   8.517  15.614  0.70  5.36           C
ATOM   1805  N  BLEU A 195      20.843  12.932  16.559  0.30  5.53           N
ATOM   1806  CA BLEU A 195      19.447  12.459  16.589  0.30  6.09           C
ATOM   1807  C  BLEU A 195      18.801  12.738  17.933  0.30  5.61           C
ATOM   1808  O  BLEU A 195      17.633  13.042  18.004  0.30  6.46           O
ATOM   1809  CB BLEU A 195      19.427  10.990  16.258  0.30  5.25           C
ATOM   1810  CG BLEU A 195      18.075  10.295  16.107  0.30  5.66           C
ATOM   1811  CD1BLEU A 195      17.280  11.004  15.013  0.30  9.49           C
ATOM   1812  CD2BLEU A 195      18.271   8.829  15.807  0.30  5.68           C
ATOM   1847  N  ATYR A 198      18.083  16.110  17.640  0.70  4.92           N
ATOM   1848  CA ATYR A 198      16.891  16.333  16.843  0.70  5.02           C
ATOM   1849  C  ATYR A 198      15.637  15.900  17.588  0.70  4.74           C
ATOM   1850  O  ATYR A 198      14.667  16.638  17.728  0.70  5.16           O
ATOM   1851  CB ATYR A 198      17.013  15.612  15.484  0.70  5.91           C
ATOM   1852  CG ATYR A 198      15.766  15.841  14.685  0.70  5.77           C
ATOM   1853  CD1ATYR A 198      15.447  17.095  14.152  0.70  6.94           C
ATOM   1854  CD2ATYR A 198      14.814  14.843  14.482  0.70  7.47           C
ATOM   1855  CE1ATYR A 198      14.299  17.342  13.429  0.70  7.49           C
ATOM   1856  CE2ATYR A 198      13.635  15.083  13.767  0.70  8.66           C
ATOM   1857  CZ ATYR A 198      13.394  16.328  13.239  0.70  7.40           C
ATOM   1858  OH ATYR A 198      12.229  16.482  12.569  0.70  9.50           O
ATOM   1859  N  BTYR A 198      17.781  16.304  18.148  0.30  5.77           N
ATOM   1860  CA BTYR A 198      16.641  16.583  17.299  0.30  6.00           C
ATOM   1861  C  BTYR A 198      15.365  16.109  18.013  0.30  4.99           C
ATOM   1862  O  BTYR A 198      14.378  16.812  18.164  0.30  6.53           O
ATOM   1863  CB BTYR A 198      16.790  15.938  15.927  0.30  6.23           C
ATOM   1864  CG BTYR A 198      15.512  16.184  15.143  0.30  6.73           C
ATOM   1865  CD1BTYR A 198      15.078  17.465  14.765  0.30  9.37           C
ATOM   1866  CD2BTYR A 198      14.715  15.091  14.811  0.30  8.99           C
ATOM   1867  CE1BTYR A 198      13.899  17.652  14.070  0.30  9.39           C
ATOM   1868  CE2BTYR A 198      13.535  15.282  14.108  0.30 10.11           C
ATOM   1869  CZ BTYR A 198      13.147  16.557  13.743  0.30  9.66           C
ATOM   1870  OH BTYR A 198      11.992  16.732  13.049  0.30 14.54           O
HETATM 3041  C1  CIT A 400      21.468  13.190  11.677  0.38  8.99           C
HETATM 3042  C2  CIT A 400      22.280  11.994  11.364  0.38  7.86           C
HETATM 3043  C3  CIT A 400      21.718  10.592  11.722  0.38  7.31           C
HETATM 3044  C4  CIT A 400      22.680   9.482  11.410  0.38  7.81           C
HETATM 3045  C5  CIT A 400      23.997   9.419  12.127  0.38  8.12           C
HETATM 3046  C6  CIT A 400      20.308  10.327  11.067  0.38  5.91           C
HETATM 3047  O1  CIT A 400      20.565  13.088  12.605  0.38  7.28           O
HETATM 3048  O2  CIT A 400      21.741  14.362  11.155  0.38  9.59           O
HETATM 3049  O3  CIT A 400      24.138   8.452  13.038  0.38 10.61           O
HETATM 3050  O4  CIT A 400      24.894  10.125  11.705  0.38 14.93           O
HETATM 3051  O5  CIT A 400      19.396   9.939  11.959  0.38  7.12           O
HETATM 3052  O6  CIT A 400      20.126  10.734   9.907  0.38  9.70           O
HETATM 3053  O7  CIT A 400      21.618  10.670  13.172  0.38  8.00           O
HETATM 3418  O   HOH A2352      15.451  18.525  10.386  1.00 18.97           O
HETATM 3420  O   HOH A2354      14.708  17.729   6.157  1.00 41.92           O
HETATM 3459  O   HOH A2393      21.027  17.116  15.310  0.70 13.37           O
HETATM 3460  O   HOH A2394      20.162  17.744  16.045  0.30  8.14           O
HETATM 3464  O   HOH A2398      18.884  18.805  15.332  0.70 11.32           O
HETATM 3468  O   HOH A2402      11.916  18.727  11.267  1.00 18.00           O
HETATM 3665  O   HOH A2599      18.218  10.669   8.288  1.00 12.65           O
HETATM 3666  O   HOH A2600      20.973  12.054  12.732  0.63 18.43           O
HETATM 3668  O   HOH A2602      22.987  10.982  12.801  0.63 31.62           O
HETATM 3669  O   HOH A2603      19.570   9.769  11.804  0.63  9.09           O
TER
HETATM 3279  O   HOH B2213      19.656  20.261   9.495  1.00 40.45           O
HETATM 3285  O   HOH B2219      15.870  20.471  12.491  1.00 38.88           O
HETATM 3550  O   HOH B2484      21.301  19.333  11.758  1.00 17.90           O
HETATM 3554  O   HOH B2488      18.501  19.093  12.751  1.00 17.74           O
HETATM 3659  O   HOH B2593      22.014  18.973  14.153  1.00 22.91           O
"""

pdb_1us0_83 = """\
ATOM    653  N   LEU A  78       5.549  -3.331  17.819  1.00  2.68           N
ATOM    654  CA  LEU A  78       4.907  -3.926  16.650  1.00  2.67           C
ATOM    655  C   LEU A  78       5.988  -4.455  15.713  1.00  2.67           C
ATOM    656  O   LEU A  78       6.827  -3.679  15.247  1.00  3.05           O
ATOM    657  CB  LEU A  78       4.149  -2.801  15.929  1.00  3.07           C
ATOM    658  CG  LEU A  78       3.452  -3.228  14.634  1.00  3.40           C
ATOM    659  CD1 LEU A  78       2.203  -4.043  14.930  1.00  4.38           C
ATOM    660  CD2 LEU A  78       3.134  -2.004  13.802  1.00  4.78           C
ATOM    661  N   TRP A  79       5.965  -5.757  15.424  1.00  2.80           N
ATOM    662  CA  TRP A  79       6.995  -6.310  14.562  1.00  2.81           C
ATOM    663  C   TRP A  79       6.718  -5.966  13.094  1.00  2.94           C
ATOM    664  O   TRP A  79       5.623  -5.594  12.681  1.00  3.41           O
ATOM    665  CB  TRP A  79       7.150  -7.833  14.821  1.00  3.01           C
ATOM    666  CG  TRP A  79       8.293  -8.378  14.033  1.00  2.97           C
ATOM    667  CD1 TRP A  79       8.225  -9.243  12.961  1.00  3.37           C
ATOM    668  CD2 TRP A  79       9.658  -7.934  14.100  1.00  2.97           C
ATOM    669  NE1 TRP A  79       9.449  -9.347  12.360  1.00  3.59           N
ATOM    670  CE2 TRP A  79      10.346  -8.568  13.041  1.00  3.36           C
ATOM    671  CE3 TRP A  79      10.358  -7.074  14.960  1.00  3.10           C
ATOM    672  CZ2 TRP A  79      11.707  -8.362  12.820  1.00  3.79           C
ATOM    673  CZ3 TRP A  79      11.704  -6.875  14.746  1.00  3.56           C
ATOM    674  CH2 TRP A  79      12.378  -7.514  13.685  1.00  3.91           C
ATOM    675  N   CYS A  80       7.794  -6.091  12.300  1.00  3.00           N
ATOM    676  CA  CYS A  80       7.844  -5.635  10.928  1.00  3.10           C
ATOM    677  C   CYS A  80       6.856  -6.331  10.003  1.00  3.15           C
ATOM    678  O   CYS A  80       6.469  -5.736   8.991  1.00  3.71           O
ATOM    679  CB  CYS A  80       9.252  -5.828  10.374  1.00  3.25           C
ATOM    680  SG  CYS A  80      10.508  -4.912  11.293  1.00  3.60           S
ATOM    681  N   THR A  81       6.477  -7.569  10.308  1.00  3.26           N
ATOM    682  CA  THR A  81       5.531  -8.319   9.493  1.00  3.35           C
ATOM    683  C   THR A  81       4.078  -7.978   9.799  1.00  3.63           C
ATOM    684  O   THR A  81       3.172  -8.623   9.260  1.00  4.38           O
ATOM    685  CB  THR A  81       5.784  -9.816   9.687  1.00  3.50           C
ATOM    686  OG1 THR A  81       5.750 -10.034  11.093  1.00  4.22           O
ATOM    687  CG2 THR A  81       7.101 -10.258   9.095  1.00  3.93           C
ATOM    688  N   TYR A  82       3.846  -6.932  10.615  1.00  3.67           N
ATOM    689  CA  TYR A  82       2.500  -6.530  10.978  1.00  4.09           C
ATOM    690  C   TYR A  82       2.256  -5.053  10.644  1.00  4.50           C
ATOM    691  O   TYR A  82       1.361  -4.451  11.209  1.00  5.48           O
ATOM    692  CB  TYR A  82       2.205  -6.782  12.465  1.00  4.34           C
ATOM    693  CG  TYR A  82       2.425  -8.204  12.894  1.00  4.48           C
ATOM    694  CD1 TYR A  82       1.484  -9.212  12.778  1.00  6.09           C
ATOM    695  CD2 TYR A  82       3.650  -8.597  13.450  1.00  5.43           C
ATOM    696  CE1 TYR A  82       1.739 -10.522  13.205  1.00  6.44           C
ATOM    697  CE2 TYR A  82       3.905  -9.883  13.872  1.00  6.09           C
ATOM    698  CZ  TYR A  82       2.934 -10.837  13.756  1.00  5.63           C
ATOM    699  OH  TYR A  82       3.255 -12.113  14.163  1.00  7.02           O
ATOM    700  N   HIS A  83       3.032  -4.481   9.706  1.00  4.55           N
ATOM    701  CG  HIS A  83       5.157  -2.189   9.548  1.00  4.23           C
ATOM    702  ND1 HIS A  83       6.061  -1.219   9.232  1.00  4.40           N
ATOM    703  CD2 HIS A  83       5.569  -2.723  10.728  1.00  4.46           C
ATOM    704  CE1 HIS A  83       7.009  -1.169  10.162  1.00  4.52           C
ATOM    705  NE2 HIS A  83       6.727  -2.079  11.087  1.00  4.33           N
ATOM    706  CA AHIS A  83       2.838  -3.163   9.157  0.56  4.11           C
ATOM    707  C  AHIS A  83       1.730  -3.173   8.083  0.56  4.43           C
ATOM    708  O  AHIS A  83       2.003  -2.962   6.902  0.56  4.82           O
ATOM    709  CB AHIS A  83       4.061  -2.481   8.569  0.56  4.90           C
ATOM    710  CA BHIS A  83       2.726  -3.051   9.515  0.44  4.76           C
ATOM    711  C  BHIS A  83       1.402  -2.766   8.778  0.44  6.93           C
ATOM    712  O  BHIS A  83       0.906  -1.664   8.864  0.44  8.11           O
ATOM    713  CB BHIS A  83       3.912  -2.461   8.743  0.44  5.47           C
ATOM    714  N  AGLU A  84       0.522  -3.439   8.546  0.56  5.45           N
ATOM    715  CA AGLU A  84      -0.535  -3.486   7.518  0.56  8.98           C
ATOM    716  C  AGLU A  84      -1.670  -2.953   8.349  0.56  7.61           C
ATOM    717  O  AGLU A  84      -1.770  -3.583   9.403  0.56 11.11           O
ATOM    718  CB AGLU A  84      -0.879  -4.948   7.188  0.56  6.08           C
ATOM    719  CG AGLU A  84       0.298  -5.676   6.580  0.56  5.58           C
ATOM    720  CD AGLU A  84       0.880  -5.117   5.294  0.56  5.05           C
ATOM    721  OE1AGLU A  84       0.183  -4.371   4.571  0.56  5.96           O
ATOM    722  OE2AGLU A  84       2.074  -5.462   4.988  0.56  5.42           O
ATOM    723  N  BGLU A  84       0.825  -3.716   8.077  0.44  6.98           N
ATOM    724  CA BGLU A  84      -0.568  -3.539   7.559  0.44  6.80           C
ATOM    725  C  BGLU A  84      -1.743  -3.132   8.442  0.44  6.97           C
ATOM    726  O  BGLU A  84      -1.949  -3.457   9.602  0.44  5.15           O
ATOM    727  CB BGLU A  84      -0.913  -4.748   6.682  0.44 11.19           C
ATOM    728  CG BGLU A  84      -0.239  -5.108   5.371  0.44 11.33           C
ATOM    729  CD BGLU A  84       1.142  -5.543   5.442  0.44  9.99           C
ATOM    730  OE1BGLU A  84       1.615  -5.746   6.565  0.44 11.23           O
ATOM    731  OE2BGLU A  84       1.756  -5.689   4.378  0.44 15.04           O
ATOM    962  N   ILE A 109       6.602   1.301  16.966  1.00  2.69           N
ATOM    963  CA  ILE A 109       6.932   0.659  15.695  1.00  2.64           C
ATOM    964  C   ILE A 109       8.385   0.182  15.829  1.00  2.81           C
ATOM    965  O   ILE A 109       9.302   0.985  15.979  1.00  3.38           O
ATOM    966  CB  ILE A 109       6.777   1.605  14.509  1.00  2.96           C
ATOM    967  CG1 ILE A 109       5.367   2.189  14.443  1.00  3.27           C
ATOM    968  CG2 ILE A 109       7.131   0.833  13.237  1.00  3.78           C
ATOM    969  CD1 ILE A 109       5.161   3.142  13.266  1.00  3.92           C
ATOM    970  N   HIS A 110       8.601  -1.149  15.821  1.00  2.63           N
ATOM    971  CA  HIS A 110       9.893  -1.695  16.224  1.00  2.60           C
ATOM    972  C   HIS A 110      11.039  -1.265  15.301  1.00  2.72           C
ATOM    973  O   HIS A 110      12.117  -0.915  15.771  1.00  2.92           O
ATOM    974  CB  HIS A 110       9.756  -3.226  16.304  1.00  2.67           C
ATOM    975  CG  HIS A 110      10.805  -3.903  17.126  1.00  2.63           C
ATOM    976  ND1 HIS A 110      10.509  -5.092  17.765  1.00  2.83           N
ATOM    977  CD2 HIS A 110      12.113  -3.568  17.353  1.00  2.89           C
ATOM    978  CE1 HIS A 110      11.640  -5.443  18.360  1.00  3.01           C
ATOM    979  NE2 HIS A 110      12.629  -4.560  18.145  1.00  2.99           N
ATOM    980  N   TRP A 111      10.779  -1.354  13.990  1.00  2.95           N
ATOM    981  CA  TRP A 111      11.705  -0.874  12.981  1.00  3.00           C
ATOM    982  C   TRP A 111      10.855  -0.174  11.928  1.00  3.18           C
ATOM    983  O   TRP A 111       9.694  -0.511  11.724  1.00  3.46           O
ATOM    984  CB  TRP A 111      12.481  -2.039  12.305  1.00  3.30           C
ATOM    985  CG  TRP A 111      13.473  -2.696  13.211  1.00  3.46           C
ATOM    986  CD1 TRP A 111      13.235  -3.667  14.155  1.00  4.01           C
ATOM    987  CD2 TRP A 111      14.857  -2.360  13.322  1.00  3.57           C
ATOM    988  NE1 TRP A 111      14.386  -3.940  14.853  1.00  4.36           N
ATOM    989  CE2 TRP A 111      15.396  -3.146  14.372  1.00  4.08           C
ATOM    990  CE3 TRP A 111      15.702  -1.470  12.638  1.00  3.76           C
ATOM    991  CZ2 TRP A 111      16.714  -3.030  14.776  1.00  4.75           C
ATOM    992  CZ3 TRP A 111      17.011  -1.368  13.023  1.00  4.30           C
ATOM    993  CH2 TRP A 111      17.522  -2.141  14.098  1.00  4.73           C
ATOM    994  N   PRO A 112      11.451   0.781  11.185  1.00  3.34           N
ATOM    995  CB  PRO A 112      11.540   2.695   9.826  1.00  4.36           C
ATOM    996  CG  PRO A 112      12.956   2.217  10.140  1.00  4.32           C
ATOM    997  CD  PRO A 112      12.778   1.375  11.391  1.00  3.71           C
ATOM    998  CA APRO A 112      10.633   1.520  10.225  0.70  3.82           C
ATOM    999  C  APRO A 112      10.177   0.651   9.056  0.70  3.76           C
ATOM   1000  O  APRO A 112       9.207   0.983   8.403  0.70  5.34           O
ATOM   1001  CA BPRO A 112      10.660   1.505  10.189  0.30  3.89           C
ATOM   1002  C  BPRO A 112      10.415   0.818   8.856  0.30  4.60           C
ATOM   1003  O  BPRO A 112       9.840   1.426   7.974  0.30  6.86           O
ATOM   1004  N   THR A 113      10.848  -0.442   8.785  1.00  3.81           N
ATOM   1005  CA  THR A 113      10.622  -1.289   7.644  1.00  3.64           C
ATOM   1006  C   THR A 113       9.485  -2.280   7.855  1.00  3.83           C
ATOM   1007  O   THR A 113       9.498  -3.061   8.801  1.00  4.78           O
ATOM   1008  CB  THR A 113      11.911  -2.120   7.432  1.00  3.93           C
ATOM   1009  OG1 THR A 113      12.275  -2.755   8.648  1.00  4.10           O
ATOM   1010  CG2 THR A 113      13.110  -1.281   7.018  1.00  4.76           C
ATOM   1011  N   GLY A 114       8.536  -2.297   6.909  1.00  3.88           N
ATOM   1012  CA  GLY A 114       7.533  -3.337   6.874  1.00  4.07           C
ATOM   1013  C   GLY A 114       7.988  -4.504   6.004  1.00  3.63           C
ATOM   1014  O   GLY A 114       8.598  -4.300   4.948  1.00  3.84           O
ATOM   1015  N   PHE A 115       7.686  -5.710   6.455  1.00  3.63           N
ATOM   1016  CA  PHE A 115       8.009  -6.932   5.733  1.00  3.51           C
ATOM   1017  C   PHE A 115       6.714  -7.651   5.363  1.00  3.64           C
ATOM   1018  O   PHE A 115       5.677  -7.485   6.005  1.00  4.19           O
ATOM   1019  CB  PHE A 115       8.824  -7.899   6.592  1.00  4.27           C
ATOM   1020  CG  PHE A 115      10.110  -7.405   7.194  1.00  3.71           C
ATOM   1021  CD1 PHE A 115      10.786  -6.287   6.752  1.00  4.16           C
ATOM   1022  CD2 PHE A 115      10.641  -8.136   8.253  1.00  4.02           C
ATOM   1023  CE1 PHE A 115      11.970  -5.886   7.380  1.00  4.68           C
ATOM   1024  CE2 PHE A 115      11.815  -7.753   8.867  1.00  4.55           C
ATOM   1025  CZ  PHE A 115      12.473  -6.608   8.438  1.00  4.82           C
ATOM   1176  N  ASER A 133       7.107  -5.240   2.629  0.56  3.99           N
ATOM   1177  CA ASER A 133       5.753  -4.756   2.436  0.56  4.16           C
ATOM   1178  C  ASER A 133       5.762  -3.397   1.732  0.56  4.26           C
ATOM   1179  O  ASER A 133       6.615  -2.567   1.920  0.56  6.01           O
ATOM   1180  CB ASER A 133       5.016  -4.611   3.762  0.56  4.36           C
ATOM   1181  OG ASER A 133       3.774  -3.954   3.584  0.56  5.13           O
ATOM   1182  N  BSER A 133       7.093  -5.190   2.622  0.44  5.10           N
ATOM   1183  CA BSER A 133       5.707  -4.870   2.265  0.44  5.43           C
ATOM   1184  C  BSER A 133       5.755  -3.421   1.807  0.44  6.60           C
ATOM   1185  O  BSER A 133       6.614  -2.637   2.163  0.44  5.18           O
ATOM   1186  CB BSER A 133       4.811  -5.153   3.467  0.44  6.15           C
ATOM   1187  OG BSER A 133       5.198  -4.245   4.464  0.44  6.90           O
ATOM   1204  N  ATHR A 135       3.595  -1.171   2.557  0.56  5.41           N
ATOM   1205  CA ATHR A 135       2.956  -0.212   3.452  0.56  5.79           C
ATOM   1206  C  ATHR A 135       3.941   0.780   4.055  0.56  7.30           C
ATOM   1207  O  ATHR A 135       4.893   0.317   4.669  0.56  7.99           O
ATOM   1208  CB ATHR A 135       2.295  -0.956   4.638  0.56  6.76           C
ATOM   1209  OG1ATHR A 135       1.393  -1.921   4.102  0.56  9.43           O
ATOM   1210  CG2ATHR A 135       1.572  -0.005   5.572  0.56  8.26           C
ATOM   1211  N  BTHR A 135       3.849  -1.125   2.503  0.44  6.93           N
ATOM   1212  CA BTHR A 135       3.420  -0.293   3.597  0.44  6.85           C
ATOM   1213  C  BTHR A 135       4.454   0.819   3.797  0.44  5.23           C
ATOM   1214  O  BTHR A 135       5.646   0.624   3.586  0.44  6.13           O
ATOM   1215  CB BTHR A 135       3.270  -1.070   4.916  0.44  8.44           C
ATOM   1216  OG1BTHR A 135       4.524  -1.646   5.301  0.44  9.29           O
ATOM   1217  CG2BTHR A 135       2.192  -2.117   4.691  0.44 11.07           C
ATOM   1218  N  AASN A 136       3.784   2.102   3.989  0.56  8.62           N
ATOM   1219  CA AASN A 136       4.821   3.055   4.502  0.56  8.87           C
ATOM   1220  C  AASN A 136       4.751   3.424   5.978  0.56  6.19           C
ATOM   1221  O  AASN A 136       3.768   2.953   6.625  0.56  5.47           O
ATOM   1222  CB AASN A 136       4.752   4.323   3.624  0.56  9.08           C
ATOM   1223  CG AASN A 136       3.377   4.872   3.713  0.56 12.60           C
ATOM   1224  OD1AASN A 136       2.982   5.295   4.733  0.56 12.09           O
ATOM   1225  ND2AASN A 136       2.690   4.978   2.570  0.56 22.23           N
ATOM   1226  N  BASN A 136       3.892   1.913   4.171  0.44  9.36           N
ATOM   1227  CA BASN A 136       4.746   3.055   4.378  0.44  7.36           C
ATOM   1228  C  BASN A 136       4.510   3.438   5.823  0.44  7.84           C
ATOM   1229  O  BASN A 136       3.484   3.108   6.449  0.44 13.28           O
ATOM   1230  CB BASN A 136       4.489   4.174   3.349  0.44  6.97           C
ATOM   1231  CG BASN A 136       3.301   5.003   3.646  0.44 10.89           C
ATOM   1232  OD1BASN A 136       2.564   4.969   4.572  0.44 13.72           O
ATOM   1233  ND2BASN A 136       3.078   5.887   2.669  0.44 18.99           N
ATOM   1234  N  AILE A 137       5.693   4.184   6.559  0.56  5.93           N
ATOM   1235  CA AILE A 137       5.664   4.512   7.976  0.56  6.32           C
ATOM   1236  C  AILE A 137       4.454   5.413   8.273  0.56  6.52           C
ATOM   1237  O  AILE A 137       4.015   5.416   9.428  0.56  8.36           O
ATOM   1238  CB AILE A 137       6.995   5.119   8.466  0.56  6.01           C
ATOM   1239  CG1AILE A 137       7.380   6.443   7.756  0.56  6.46           C
ATOM   1240  CG2AILE A 137       8.107   4.100   8.356  0.56  7.48           C
ATOM   1241  CD1AILE A 137       8.625   7.134   8.253  0.56  8.08           C
ATOM   1242  N  BILE A 137       5.542   4.177   6.142  0.44  6.94           N
ATOM   1243  CA BILE A 137       5.666   4.543   7.544  0.44  6.46           C
ATOM   1244  C  BILE A 137       4.435   5.276   8.103  0.44  7.93           C
ATOM   1245  O  BILE A 137       3.947   5.102   9.225  0.44 10.29           O
ATOM   1246  CB BILE A 137       6.956   5.358   7.699  0.44  7.51           C
ATOM   1247  CG1BILE A 137       8.193   4.490   7.608  0.44  8.75           C
ATOM   1248  CG2BILE A 137       6.917   6.145   8.994  0.44 10.05           C
ATOM   1249  CD1BILE A 137       9.510   5.239   7.560  0.44 12.80           C
ATOM   1266  N   THR A 140       1.444   3.151   8.697  1.00  5.66           N
ATOM   1267  CA  THR A 140       1.625   2.403   9.922  1.00  5.19           C
ATOM   1268  C   THR A 140       1.371   3.309  11.135  1.00  4.68           C
ATOM   1269  O   THR A 140       0.734   2.883  12.108  1.00  4.87           O
ATOM   1270  CB  THR A 140       3.015   1.766   9.939  1.00  5.70           C
ATOM   1271  OG1 THR A 140       3.196   0.990   8.758  1.00  7.35           O
ATOM   1272  CG2 THR A 140       3.169   0.848  11.133  1.00  6.06           C
HETATM 3258  O   HOH A2192       8.341  -2.987  13.103  1.00  3.44           O
HETATM 3269  O   HOH A2203       0.558   0.435   7.312  0.44 12.75           O
HETATM 3272  O   HOH A2206       4.011  -5.503   7.097  1.00  7.94           O
HETATM 3348  O   HOH A2282       6.791   2.057   9.696  0.44 13.20           O
HETATM 3349  O   HOH A2283       7.703   1.328   6.092  0.44 10.03           O
HETATM 3350  O   HOH A2284       8.499   0.709   5.936  0.56  6.26           O
HETATM 3375  O   HOH A2309       7.519  -1.331   4.340  1.00  8.93           O
HETATM 3381  O   HOH A2315       5.918   0.338   7.031  1.00 11.95           O
TER
"""

pdb_1exr_sym = """\
CRYST1   25.015   29.415   52.761  89.54  86.10  82.39 P 1
SCALE1      0.039976 -0.005341 -0.002731        0.00000
SCALE2      0.000000  0.034298  0.000034        0.00000
SCALE3      0.000000  0.000000  0.018997        0.00000
ATOM    518  N   GLY A  59      32.173  13.282  -5.408  1.00 14.69           N
ATOM    519  CA  GLY A  59      31.512  12.626  -6.522  1.00 16.07           C
ATOM    520  C   GLY A  59      32.329  12.708  -7.800  1.00 14.80           C
ATOM    521  O   GLY A  59      31.713  12.494  -8.863  1.00 17.34           O
ATOM    522  N   ASN A  60      33.616  12.979  -7.798  1.00 13.93           N
ATOM    523  CA  ASN A  60      34.395  13.043  -9.011  1.00 15.63           C
ATOM    524  C   ASN A  60      34.849  11.676  -9.428  1.00 14.22           C
ATOM    525  O   ASN A  60      35.497  11.598 -10.498  1.00 17.75           O
ATOM    526  CB  ASN A  60      35.533  14.043  -8.916  1.00 15.69           C
ATOM    527  CG  ASN A  60      36.616  13.666  -7.960  1.00  9.01           C
ATOM    528  OD1 ASN A  60      36.522  12.667  -7.198  1.00  9.98           O
ATOM    529  ND2 ASN A  60      37.638  14.489  -7.969  1.00 10.09           N
ATOM    530  N   GLY A  61      34.535  10.599  -8.680  1.00 11.70           N
ATOM    531  CA  GLY A  61      34.841   9.258  -9.131  1.00 11.82           C
ATOM    532  C   GLY A  61      36.078   8.633  -8.527  1.00  9.90           C
ATOM    533  O   GLY A  61      36.330   7.471  -8.747  1.00 11.11           O
ATOM    935  N   ARG A  94      34.136 -16.706  29.792  1.00 13.44           N
ATOM    936  CA  ARG A  94      35.261 -17.201  30.584  1.00 15.05           C
ATOM    937  C   ARG A  94      35.967 -18.368  29.930  1.00 15.88           C
ATOM    938  O   ARG A  94      37.204 -18.484  30.059  1.00 19.56           O
ATOM    939  CB  ARG A  94      34.779 -17.552  31.990  1.00 16.59           C
ATOM    940  CG  ARG A  94      34.216 -16.398  32.818  1.00 19.74           C
ATOM    941  CD  ARG A  94      33.862 -16.888  34.231  1.00 21.05           C
ATOM    942  NE  ARG A  94      35.044 -17.164  35.019  1.00 23.55           N
ATOM    943  CZ  ARG A  94      35.279 -18.159  35.849  1.00 24.35           C
ATOM    944  NH1 ARG A  94      34.386 -19.122  36.089  1.00 39.46           N
ATOM    945  NH2 ARG A  94      36.442 -18.178  36.481  1.00 31.24           N
ATOM   1005  N   ALA A 103      26.078 -19.251  33.669  1.00 11.12           N
ATOM   1006  CA  ALA A 103      27.044 -19.412  34.744  1.00 11.64           C
ATOM   1007  C   ALA A 103      28.186 -18.414  34.623  1.00  9.77           C
ATOM   1008  O   ALA A 103      28.620 -17.809  35.617  1.00 10.63           O
ATOM   1009  CB  ALA A 103      27.572 -20.847  34.756  1.00 14.26           C
ATOM   1010  N   GLU A 104      28.647 -18.214  33.408  1.00  9.68           N
ATOM   1011  CA  GLU A 104      29.733 -17.224  33.187  1.00  9.38           C
ATOM   1012  C   GLU A 104      29.282 -15.833  33.521  1.00  9.15           C
ATOM   1013  O   GLU A 104      29.995 -15.072  34.194  1.00 10.22           O
ATOM   1014  CB  GLU A 104      30.184 -17.298  31.725  1.00  9.64           C
ATOM   1015  CG  GLU A 104      30.993 -18.518  31.406  1.00 11.21           C
ATOM   1016  CD  GLU A 104      31.281 -18.734  29.941  1.00 12.08           C
ATOM   1017  OE1 GLU A 104      30.537 -18.227  29.084  1.00 11.43           O
ATOM   1018  OE2 GLU A 104      32.279 -19.456  29.640  1.00 14.06           O
ATOM   1027  N   ARG A 106      26.895 -14.883  35.508  1.00  9.22           N
ATOM   1028  CA  ARG A 106      26.719 -14.722  36.953  1.00  9.65           C
ATOM   1029  C   ARG A 106      28.056 -14.566  37.655  1.00  9.45           C
ATOM   1030  O   ARG A 106      28.210 -13.758  38.572  1.00 10.64           O
ATOM   1031  CB  ARG A 106      25.870 -15.881  37.488  1.00  9.92           C
ATOM   1032  CG  ARG A 106      25.440 -15.634  38.915  1.00 11.16           C
ATOM   1033  CD  ARG A 106      24.342 -16.590  39.346  1.00 11.19           C
ATOM   1034  NE  ARG A 106      24.571 -17.943  39.133  1.00 12.34           N
ATOM   1035  CZ  ARG A 106      25.209 -18.729  39.967  1.00 12.70           C
ATOM   1036  NH1 ARG A 106      25.783 -18.243  41.075  1.00 13.78           N
ATOM   1037  NH2 ARG A 106      25.363 -20.040  39.717  1.00 18.51           N
ATOM   1038  N  AHIS A 107      29.019 -15.368  37.207  0.34  8.71           N
ATOM   1039  CA AHIS A 107      30.341 -15.358  37.814  0.34 10.82           C
ATOM   1040  C  AHIS A 107      31.077 -14.031  37.693  0.34 12.08           C
ATOM   1041  O  AHIS A 107      31.706 -13.553  38.647  0.34 10.34           O
ATOM   1042  CB AHIS A 107      31.184 -16.457  37.142  0.34 11.82           C
ATOM   1043  CG AHIS A 107      32.441 -16.649  37.944  0.34 14.84           C
ATOM   1044  ND1AHIS A 107      33.504 -15.796  37.976  0.34 15.41           N
ATOM   1045  CD2AHIS A 107      32.745 -17.648  38.805  0.34 17.00           C
ATOM   1046  CE1AHIS A 107      34.432 -16.266  38.805  0.34 16.26           C
ATOM   1047  NE2AHIS A 107      33.988 -17.399  39.330  0.34 16.81           N
ATOM   1048  N  BHIS A 107      29.054 -15.396  37.312  0.66 10.87           N
ATOM   1049  CA BHIS A 107      30.352 -15.252  37.989  0.66 10.92           C
ATOM   1050  C  BHIS A 107      30.956 -13.870  37.775  0.66 11.22           C
ATOM   1051  O  BHIS A 107      31.380 -13.167  38.696  0.66 11.35           O
ATOM   1052  CB BHIS A 107      31.290 -16.355  37.492  0.66 11.21           C
ATOM   1053  CG BHIS A 107      32.739 -16.176  37.885  0.66 13.34           C
ATOM   1054  ND1BHIS A 107      33.159 -16.647  39.130  0.66 16.57           N
ATOM   1055  CD2BHIS A 107      33.826 -15.605  37.308  0.66 17.50           C
ATOM   1056  CE1BHIS A 107      34.453 -16.356  39.258  0.66 16.43           C
ATOM   1057  NE2BHIS A 107      34.900 -15.723  38.175  0.66 16.67           N
ATOM   1058  N   VAL A 108      30.977 -13.409  36.516  1.00 11.16           N
ATOM   1059  CA  VAL A 108      31.596 -12.123  36.256  1.00 11.53           C
ATOM   1060  C   VAL A 108      30.769 -11.023  36.930  1.00 11.56           C
ATOM   1061  O   VAL A 108      31.342 -10.062  37.513  1.00 12.76           O
ATOM   1062  CB  VAL A 108      31.797 -11.899  34.755  1.00 13.24           C
ATOM   1063  CG1 VAL A 108      32.393 -10.525  34.485  1.00 16.48           C
ATOM   1064  CG2 VAL A 108      32.707 -13.012  34.228  1.00 15.11           C
ATOM   1073  N   THR A 110      29.187 -11.175  39.754  1.00 11.08           N
ATOM   1074  CA  THR A 110      29.596 -11.157  41.171  1.00 11.80           C
ATOM   1075  C   THR A 110      30.881 -10.390  41.314  1.00 11.50           C
ATOM   1076  O   THR A 110      31.060  -9.583  42.258  1.00 12.59           O
ATOM   1077  CB  THR A 110      29.715 -12.619  41.665  1.00 14.17           C
ATOM   1078  OG1 THR A 110      28.415 -13.216  41.608  1.00 18.55           O
ATOM   1079  CG2 THR A 110      30.144 -12.649  43.113  1.00 21.50           C
ATOM   1080  N  AASN A 111      31.845 -10.634  40.457  0.72 11.72           N
ATOM   1081  CA AASN A 111      33.111  -9.834  40.563  0.72 10.91           C
ATOM   1082  C  AASN A 111      32.887  -8.368  40.312  0.72 10.52           C
ATOM   1083  O  AASN A 111      33.646  -7.556  40.907  0.72 12.14           O
ATOM   1084  CB AASN A 111      34.138 -10.489  39.672  0.72 12.73           C
ATOM   1085  CG AASN A 111      34.671 -11.785  40.267  0.72 15.99           C
ATOM   1086  OD1AASN A 111      35.190 -12.663  39.560  0.72 19.42           O
ATOM   1087  ND2AASN A 111      34.605 -12.004  41.573  0.72 15.05           N
ATOM   1088  N  BASN A 111      31.778 -10.642  40.375  0.28 11.41           N
ATOM   1089  CA BASN A 111      33.114 -10.070  40.212  0.28 11.93           C
ATOM   1090  C  BASN A 111      33.055  -8.580  39.879  0.28 11.88           C
ATOM   1091  O  BASN A 111      34.093  -7.890  39.997  0.28 15.67           O
ATOM   1092  CB BASN A 111      33.896 -10.846  39.145  0.28 15.73           C
ATOM   1093  CG BASN A 111      35.370 -10.896  39.477  0.28 19.65           C
ATOM   1094  OD1BASN A 111      35.791 -10.030  40.232  0.28 36.84           O
ATOM   1095  ND2BASN A 111      36.139 -11.837  38.963  0.28 29.48           N
ATOM   1166  N  AASP A 119      17.580 -16.656  38.457  0.34  8.94           N
ATOM   1167  CA AASP A 119      16.379 -16.069  37.854  0.34  8.53           C
ATOM   1168  C  AASP A 119      16.650 -14.679  37.279  0.34  8.66           C
ATOM   1169  O  AASP A 119      16.090 -14.294  36.242  0.34  6.57           O
ATOM   1170  CB AASP A 119      15.227 -15.989  38.848  0.34  8.93           C
ATOM   1171  CG AASP A 119      14.437 -17.290  38.953  0.34  9.03           C
ATOM   1172  OD1AASP A 119      14.661 -18.205  38.150  0.34 10.97           O
ATOM   1173  OD2AASP A 119      13.596 -17.271  39.874  0.34 10.12           O
ATOM   1174  N  BASP A 119      17.556 -16.687  38.419  0.66  8.99           N
ATOM   1175  CA BASP A 119      16.428 -16.100  37.662  0.66 10.11           C
ATOM   1176  C  BASP A 119      16.837 -14.731  37.126  0.66  8.54           C
ATOM   1177  O  BASP A 119      16.543 -14.410  35.961  0.66  8.09           O
ATOM   1178  CB BASP A 119      15.180 -16.037  38.509  0.66 11.81           C
ATOM   1179  CG BASP A 119      13.933 -15.929  37.672  0.66 18.19           C
ATOM   1180  OD1BASP A 119      13.105 -15.097  37.998  0.66 32.94           O
ATOM   1181  OD2BASP A 119      13.752 -16.725  36.715  0.66 27.56           O
TER
HETATM 1504  O   HOH A2031      28.414 -19.030  38.108  1.00 21.58           O
HETATM 1527  O   HOH A2054      31.784 -19.349  35.227  1.00 24.96           O
HETATM 1538  O   HOH A2065      36.943 -15.021  40.078  1.00 30.53           O
HETATM 1542  O   HOH A2069      31.525 -17.270  41.211  1.00 31.52           O
HETATM 1548  O   HOH A2075      11.848 -19.490  39.935  1.00 31.49           O
HETATM 1558  O   HOH A2085      14.146 -15.462  34.876  1.00 28.44           O
HETATM 1565  O   HOH A2092      29.078 -17.394  40.211  1.00 24.67           O
HETATM 1575  O   HOH A2102      32.697 -14.477  41.218  1.00 30.96           O
HETATM 1577  O   HOH A2104      11.478 -14.774  35.326  1.00 31.02           O
HETATM 1582  O   HOH A2109      10.155 -12.003  36.709  1.00 30.46           O
HETATM 1588  O   HOH A2115      31.012 -20.152  38.014  1.00 34.73           O
HETATM 1604  O   HOH A2131      14.288 -12.648  38.691  1.00 26.44           O
HETATM 1633  O   HOH A2160      31.707   8.248 -10.026  1.00 38.85           O
HETATM 1634  O   HOH A2161      12.461 -17.449  40.907  0.60 18.47           O
HETATM 1646  O   HOH A2173      11.730 -12.749  34.283  1.00 40.97           O
"""

if __name__ == "__main__":
  test_done_1exr()
  test_check_failure_restores()
  test_strict_restores_first()
  test_unpin_contract()
  test_restore_check_catches_omission()
  test_deadline()
  test_any_error_keeps_first()
  test_first_run_error_propagates()
  test_first_run_record_error()
  test_fallback_deadline()
  test_report_errors()
  test_memory_error_propagates()
  test_abort_in_rerun()
  test_diagnostic_errors()
  test_decide_error_second_his()
  test_set_flip_states()
  test_optimize()
  test_real_deadline()
  test_two_candidates()
  test_failure_after_first_his()
  test_blank_ring_rule_in_driver()
  test_check_written_marks()
  test_hint_and_report()
  test_blank_ring_split_h()
  test_metal_kept()
  test_multi_model()
  test_qmready()
  test_qmready_stderr()
  test_auto()
  test_flip_report_words()
  test_search_atom_info()
  test_blank_ring_last_pass()
  print("OK")
