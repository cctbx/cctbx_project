from __future__ import absolute_import, division, print_function
import io, json
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry
import iotbx.pdb
from mmtbx.validation import cbetadev

# phi_psi_correction_region = positive: only residues in the positive phi
# region are corrected, and the strength comes from those residues alone

# 3noq chain A 10-14, backbone and CB: GLU 12 has phi 68.7, psi 25.7; PRO 11
# and VAL 13 have corrections at negative phi
pdb_str = """\
ATOM     78  N   PHE A  10       1.018  25.773  25.906  1.00  7.05           N
ATOM     79  CA  PHE A  10       1.266  26.317  27.227  1.00  6.77           C
ATOM     80  C   PHE A  10       2.741  26.229  27.573  1.00  5.97           C
ATOM     81  O   PHE A  10       3.483  25.403  27.044  1.00  6.02           O
ATOM     82  CB  PHE A  10       0.449  25.585  28.302  1.00  6.81           C
ATOM     89  N   PRO A  11       3.187  27.102  28.485  1.00  6.59           N
ATOM     90  CA  PRO A  11       4.592  27.045  28.938  1.00  6.67           C
ATOM     91  C   PRO A  11       4.963  25.654  29.424  1.00  5.60           C
ATOM     92  O   PRO A  11       4.219  25.024  30.171  1.00  6.09           O
ATOM     93  CB  PRO A  11       4.637  28.085  30.056  1.00  8.19           C
ATOM     96  N   GLU A  12       6.144  25.180  29.020  1.00  5.53           N
ATOM     97  CA  GLU A  12       6.663  23.872  29.432  1.00  5.74           C
ATOM     98  C   GLU A  12       5.862  22.715  28.829  1.00  5.27           C
ATOM     99  O   GLU A  12       5.905  21.597  29.338  1.00  5.87           O
ATOM    100  CB  GLU A  12       6.838  23.763  30.948  1.00  7.38           C
ATOM    105  N   VAL A  13       5.218  22.955  27.680  1.00  5.26           N
ATOM    106  CA  VAL A  13       4.665  21.861  26.891  1.00  5.07           C
ATOM    107  C   VAL A  13       5.795  20.934  26.440  1.00  4.61           C
ATOM    108  O   VAL A  13       6.918  21.352  26.172  1.00  5.20           O
ATOM    109  CB  VAL A  13       3.897  22.443  25.681  1.00  5.43           C
ATOM    112  N   GLN A  14       5.451  19.653  26.286  1.00  4.78           N
ATOM    113  CA  GLN A  14       6.303  18.683  25.617  1.00  5.03           C
ATOM    114  C   GLN A  14       6.040  18.876  24.131  1.00  4.88           C
ATOM    115  O   GLN A  14       4.990  18.476  23.636  1.00  5.66           O
ATOM    116  CB  GLN A  14       5.945  17.271  26.043  1.00  5.89           C
END
"""

def get_hierarchy():
  return iotbx.pdb.input(source_info=None,
    lines=pdb_str.splitlines()).construct_hierarchy()

def run_cbetadev(out=None, quiet=True, **kwds):
  if out is None: out = null_out()
  return cbetadev.cbetadev(get_hierarchy(), out=out, quiet=quiet, **kwds)

def by_resseq(c):
  return dict((r.resseq.strip(), r) for r in c.results)

def only_strength(c):
  (s,) = c.phi_psi_correction.values()
  return s

def exercise_only_positive_residue_changes():
  base = by_resseq(run_cbetadev())
  every = by_resseq(run_cbetadev(apply_phi_psi_correction=True,
                                 phi_psi_correction_strength=1.))
  pos = by_resseq(run_cbetadev(apply_phi_psi_correction=True,
                               phi_psi_correction_strength=1.,
                               phi_psi_correction_region='positive'))
  # the fragment must have corrected residues outside the region
  assert every['13'].deviation != base['13'].deviation
  assert pos['12'].deviation != base['12'].deviation
  assert approx_equal(pos['12'].deviation, every['12'].deviation)
  for resseq in base:
    if resseq == '12': continue
    assert pos[resseq].deviation == base[resseq].deviation, resseq

def exercise_strength_from_region_residues():
  c = run_cbetadev(apply_phi_psi_correction=True,
                   phi_psi_correction_region='positive')
  s = only_strength(c)
  assert s.n == 1
  # one residue against a prior of 10 residues at 0.5
  assert approx_equal(s.pulled, (s.fitted + 5)/11.)

def exercise_default_is_all():
  a = run_cbetadev(apply_phi_psi_correction=True)
  b = run_cbetadev(apply_phi_psi_correction=True,
                   phi_psi_correction_region='all')
  assert only_strength(a).n == only_strength(b).n > 1
  assert [r.deviation for r in a.results] == [r.deviation for r in b.results]

def exercise_bad_region():
  try:
    run_cbetadev(apply_phi_psi_correction=True,
                 phi_psi_correction_region='negative')
  except Sorry as e:
    assert 'phi_psi_correction_region' in str(e)
  else:
    raise AssertionError('Sorry not raised')

def exercise_json_and_summary():
  out = io.StringIO()
  c = run_cbetadev(out=out, quiet=False, apply_phi_psi_correction=True,
                   phi_psi_correction_region='positive')
  js = json.loads(c.as_JSON())['summary_results']['']
  assert js['phi_psi_correction_region'] == 'positive'
  assert 'positive phi only' in out.getvalue()
  js = json.loads(run_cbetadev(apply_phi_psi_correction=True).as_JSON())
  assert js['summary_results']['']['phi_psi_correction_region'] == 'all'
  out = io.StringIO()
  run_cbetadev(out=out, quiet=False, apply_phi_psi_correction=True)
  assert 'positive phi only' not in out.getvalue()

def exercise_program_option():
  from iotbx.cli_parser import run_program
  from mmtbx.programs import cbetadev as cbetadev_program
  import os, tempfile
  d = tempfile.mkdtemp()
  fn = os.path.join(d, 'frag.pdb')
  with open(fn, 'w') as f: f.write(pdb_str)
  try:
    run_program(program_class=cbetadev_program.Program,
                args=[fn, 'cbetadev.phi_psi_correction_region=positive'],
                logger=null_out())
  except Sorry as e:
    assert 'requires apply_phi_psi_correction' in str(e)
  else:
    raise AssertionError('Sorry not raised')
  run_program(program_class=cbetadev_program.Program,
              args=[fn, 'cbetadev.apply_phi_psi_correction=True',
                    'cbetadev.phi_psi_correction_region=positive'],
              logger=null_out())

def run():
  exercise_only_positive_residue_changes()
  exercise_strength_from_region_residues()
  exercise_default_is_all()
  exercise_bad_region()
  exercise_json_and_summary()
  exercise_program_option()
  print('OK')

if __name__ == '__main__':
  run()
