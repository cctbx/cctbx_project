from __future__ import absolute_import, division, print_function
import io, json
from math import cos, sin, radians
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry
import iotbx.pdb
from mmtbx.validation import cbetadev
from mmtbx.validation import cbd_utils
from mmtbx.validation.regression.tst_cbetadev_05 import get_hierarchy, pdb_str

# per-model strength of the phi/psi correction on the a2.pdb fragment of
# tst_cbetadev_05 (residues 11-16 have corrections)

def run_cbetadev(hierarchy, out=None, quiet=True, **kwds):
  if out is None: out = null_out()
  return cbetadev.cbetadev(hierarchy, out=out, quiet=quiet, **kwds)

def only_strength(c):
  (s,) = c.phi_psi_correction.values()
  return s

def exercise_strength_record():
  c = run_cbetadev(get_hierarchy(), apply_phi_psi_correction=True)
  s = only_strength(c)
  assert s.n == 6 and s.fitted is not None and 0 <= s.applied <= 1
  assert run_cbetadev(get_hierarchy()).phi_psi_correction is None

def exercise_user_strength():
  h = get_hierarchy()
  base = run_cbetadev(h)
  for strength in [0., 0.5, 1.]:
    c = run_cbetadev(h, apply_phi_psi_correction=True,
                     phi_psi_correction_strength=strength)
    s = only_strength(c)
    assert s.user_set and s.applied == strength and not s.check_failed
    if strength == 0.:
      assert [r.deviation for r in c.results] == \
             [r.deviation for r in base.results]

def exercise_strength_scales():
  # the table vector c is v - (corrected at strength 1); at 0.5 the
  # corrected vector must be v - 0.5*c
  h = get_hierarchy()
  kw = dict(apply_phi_psi_correction=True)
  base = run_cbetadev(h)
  full = run_cbetadev(h, phi_psi_correction_strength=1., **kw)
  half = run_cbetadev(h, phi_psi_correction_strength=0.5, **kw)
  def xy(r):
    t = radians(r.dihedral_NABB)
    return r.deviation*cos(t), r.deviation*sin(t)
  n = 0
  for b, f, r in zip(base.results, full.results, half.results):
    assert b.id_str() == f.id_str() == r.id_str()
    (bx, by), (fx, fy) = xy(b), xy(f)
    c = (bx - fx, by - fy)
    if abs(c[0]) + abs(c[1]) < 1e-9:
      assert r.deviation == b.deviation
      continue
    n += 1
    assert approx_equal((r.deviation, r.dihedral_NABB),
      cbd_utils.apply_correction(b.deviation, b.dihedral_NABB, c, 0.5),
      eps=1e-6)
  assert n == 6, n

def exercise_outliers_only_same_strength():
  h = get_hierarchy()
  a = run_cbetadev(h, apply_phi_psi_correction=True)
  b = run_cbetadev(h, apply_phi_psi_correction=True, outliers_only=True)
  assert only_strength(a).fitted == only_strength(b).fitted
  assert only_strength(a).n == only_strength(b).n

def exercise_no_correctable_residues():
  # MLY has no table without use_parent_for_modified_residues
  h = get_hierarchy()
  for ag in h.atom_groups():
    ag.resname = 'MLY'
  out = io.StringIO()
  c = run_cbetadev(h, out=out, quiet=False, apply_phi_psi_correction=True)
  s = only_strength(c)
  assert s.n == 0 and s.fitted is None and s.applied == 0.
  assert 'no residues with a correction' in out.getvalue(), out.getvalue()
  js = json.loads(c.as_JSON())['summary_results']
  (v,) = js.values()
  assert v['phi_psi_correction_fitted_strength'] is None
  assert v['phi_psi_correction_residues'] == 0

def exercise_altlocs():
  # residue 13 with its CB split into A and B, B moved 0.1 A along x
  lines = []
  for line in pdb_str.splitlines():
    if line[12:16] == ' CB ' and line[22:26] == '  13':
      a = line[:16] + 'A' + line[17:54] + '  0.50' + line[60:]
      x = float(line[30:38]) + 0.1
      b = line[:16] + 'B' + line[17:30] + '%8.3f' % x + line[38:54] + \
          '  0.50' + line[60:]
      lines += [a, b]
    else:
      lines.append(line)
  h = iotbx.pdb.input(source_info=None, lines=lines).construct_hierarchy()
  c = run_cbetadev(h, apply_phi_psi_correction=True,
                   phi_psi_correction_strength=0.5)
  s = only_strength(c)
  # generate_protein_threes yields a three per conformer, keyed by the CB's
  # altloc, so conformer B gets its own phi/psi
  assert s.n == 7, s.n
  alts = [r for r in c.results if r.resseq.strip() == '13']
  assert sorted(r.altloc for r in alts) == ['A', 'B'], [r.altloc for r in alts]

def exercise_summary_and_json():
  h = get_hierarchy(n_models=2)
  out = io.StringIO()
  c = run_cbetadev(h, out=out, quiet=False, apply_phi_psi_correction=True,
                   phi_psi_correction_strength=0.5)
  text = out.getvalue()
  assert text.count('(user set)') == 2, text
  assert 'Correction strength model "1" : fitted' in text, text
  assert 'Correction strength model "2" : fitted' in text, text
  js = json.loads(c.as_JSON())['summary_results']
  assert len(js) == 2
  for v in js.values():
    assert v['phi_psi_correction_applied_strength'] == 0.5
    assert v['phi_psi_correction_user_set'] is True
    assert v['phi_psi_correction_check_failed'] is False
    assert v['phi_psi_correction_residues'] == 6
  # no correction, no fields
  js = json.loads(run_cbetadev(h).as_JSON())['summary_results']
  for v in js.values():
    assert not [k for k in v if k.startswith('phi_psi_correction')]

def exercise_format_summary():
  from libtbx import group_args
  def rec(**kw):
    d = dict(n=412, fitted=0.31, pulled=0.31, applied=0.31, user_set=False,
             check_failed=False, mean_uncorrected=0.061, mean_corrected=0.058)
    d.update(kw)
    return group_args(**d)
  f = cbetadev.format_correction_summary
  assert f('', rec()) == [
    '  Correction strength : fitted 0.31, applied 0.31 (412 residues)']
  assert f(' model "2"', rec(applied=1., user_set=True)) == [
    '  Correction strength model "2" : fitted 0.31, applied 1.00 (user set) '
    '(412 residues)']
  assert f('', rec(n=0, fitted=None, pulled=None, applied=0.)) == [
    '  Correction strength : no residues with a correction']
  lines = f('', rec(fitted=-0.42, pulled=-0.3, applied=0.))
  assert lines[1] == ('  WARNING: fitted strength -0.42 is below 0; the C-beta '
    'offsets in this model run against the high-resolution phi/psi pattern')
  lines = f('', rec(applied=0., check_failed=True, mean_corrected=0.065))
  assert lines[-1] == ('  Correction not applied: it would raise the mean '
    'deviation from 0.061 to 0.065 Angstrom')

def exercise_program_option():
  from iotbx.cli_parser import run_program
  from mmtbx.programs import cbetadev as cbetadev_program
  file_name = 'tst_cbetadev_06.pdb'
  with open(file_name, 'w') as f:
    f.write(pdb_str)
  try:
    run_program(program_class=cbetadev_program.Program,
                args=[file_name, 'phi_psi_correction_strength=1'],
                logger=null_out())
  except Sorry as e:
    assert str(e).find(
      'phi_psi_correction_strength requires apply_phi_psi_correction') > -1
  else:
    raise AssertionError('Sorry not raised')
  run_program(program_class=cbetadev_program.Program,
              args=[file_name, 'apply_phi_psi_correction=True',
                    'phi_psi_correction_strength=1'],
              logger=null_out())

def run():
  exercise_strength_record()
  exercise_user_strength()
  exercise_strength_scales()
  exercise_outliers_only_same_strength()
  exercise_no_correctable_residues()
  exercise_altlocs()
  exercise_summary_and_json()
  exercise_format_summary()
  exercise_program_option()
  print('OK')

if __name__ == '__main__':
  run()
