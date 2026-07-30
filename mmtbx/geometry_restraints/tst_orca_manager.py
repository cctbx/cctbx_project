"""Fast unit tests for mmtbx.geometry_restraints.orca_manager.

No ORCA binary is needed: every exercise checks the input string the manager
generates (optimisation keyword, %pal parallel line, %geom Tol overrides, the
NumFreq Hessian input) or the assembled command line. The one exercise that
drives get_hessian() stubs out run_cmd so ORCA is never launched.
"""
from __future__ import absolute_import, division, print_function
import os
import shutil
import tempfile

import iotbx.pdb
from scitbx.array_family import flex
from libtbx.utils import Sorry

from mmtbx.geometry_restraints import orca_manager

# A tiny Zn(II) site: two water donors (O) and the metal (ZN).
METAL_PDB = """\
HETATM    1  O   HOH A   1      59.041  72.758  32.575  1.00 20.00           O
HETATM    2  O   HOH A   2      57.856  75.840  29.342  1.00 20.00           O
HETATM    3 ZN    ZN A   3      59.483  73.605  29.897  1.00 20.00          ZN
"""


def make_manager(preamble, ligand_flags=(1, 0, 0), pdb_str=METAL_PDB,
                 method='HF', basis_set='def2-SVP', solvent_model='CPCM',
                 charge=2, multiplicity=1, nproc=1):
  """Build an orca_manager from a PDB string. The hierarchy is stashed on the
  manager so atom.parent() (used by id_str) stays valid."""
  hier = iotbx.pdb.input(source_info=None,
                         lines=pdb_str.splitlines()).construct_hierarchy()
  atoms = hier.atoms()
  assert atoms.size() == len(ligand_flags)
  m = orca_manager.orca_manager(atoms, method, basis_set, solvent_model,
                                charge, multiplicity, nproc, preamble=preamble)
  m.program_goal = 'opt'
  m.set_ligand_atoms(list(ligand_flags))
  m._tst_hier = hier  # keep the hierarchy (and atom.parent() chain) alive
  return m


def bang_line(text):
  """The single ORCA keyword line, e.g. '! HF def2-SVP CPCM LooseOpt '."""
  for line in text.splitlines():
    if line.startswith('! '):
      return line
  return ''


def exercise_parallel_and_optimisation():
  """%pal only when nproc>1; optimisation keyword is configurable."""
  # defaults: no %pal, LooseOpt
  out = make_manager('def').get_input_lines()
  assert '%pal' not in out, out
  assert 'LooseOpt' in bang_line(out), bang_line(out)

  # nproc>1 emits %pal in the header
  out4 = make_manager('p4', nproc=4).get_input_lines()
  assert '%pal nprocs 4 end' in out4, out4

  # optimisation='Opt' gives '! ... Opt' and no LooseOpt anywhere
  m = make_manager('opt')
  m.optimisation = 'Opt'
  out_opt = m.get_input_lines()
  assert 'LooseOpt' not in out_opt, out_opt
  assert ' Opt' in bang_line(out_opt), bang_line(out_opt)

  # gradients still map to EnGrad regardless of optimisation
  grad = make_manager('grad').get_input_lines(gradients=True)
  assert 'EnGrad' in bang_line(grad), bang_line(grad)
  assert 'LooseOpt' not in grad
  print('  parallel/optimisation OK')


def exercise_geom_convergence():
  """Tol* overrides land inside %geom, ahead of Constraints."""
  preset = orca_manager.GEOM_CONVERGENCE['looser']

  # named preset, with constraints present (ligand_flags [1,0,0])
  m = make_manager('geom')
  m.geom_convergence = 'looser'
  out = m.get_input_lines()
  gi, ti, ci = out.index('%geom'), out.index('TolMaxD'), out.index('Constraints')
  assert gi < ti < ci, out                     # Tol lines between %geom and Constraints
  for key in preset:
    assert '  %s ' % key in out, (key, out)
  assert 'TolMaxD 0.03' in out, out            # 3e-2 formats to 0.03

  # a plain dict is injected verbatim
  m.geom_convergence = {'TolMaxD': 0.05, 'TolE': 1e-3}
  out2 = m.get_input_lines()
  assert 'TolMaxD 0.05' in out2, out2
  assert 'TolE 0.001' in out2, out2
  g2, t2, c2 = out2.index('%geom'), out2.index('TolMaxD'), out2.index('Constraints')
  assert g2 < t2 < c2, out2

  # default (None) injects no Tol lines
  m.geom_convergence = None
  assert 'Tol' not in m.get_input_lines()

  # Tol overrides but no constraints: %geom block with no Constraints sub-block
  mnc = make_manager('nocon', ligand_flags=(1, 1, 1))
  mnc.geom_convergence = 'looser'
  out_nc = mnc.get_input_lines()
  assert '%geom' in out_nc, out_nc
  assert 'Constraints' not in out_nc, out_nc
  assert 'TolMaxD 0.03' in out_nc, out_nc
  print('  geom convergence OK')


def exercise_get_cmd():
  """self.exe wins; otherwise $PHENIX_ORCA; otherwise Sorry."""
  m = make_manager('cmd')
  saved = os.environ.pop('PHENIX_ORCA', None)
  try:
    # neither exe nor env -> Sorry
    m.exe = None
    try:
      m.get_cmd()
      raise AssertionError('expected Sorry when neither exe nor env is set')
    except Sorry:
      pass
    # explicit exe is used
    m.exe = '/opt/orca/orca'
    assert m.get_cmd() == '/opt/orca/orca orca_cmd.in', m.get_cmd()
    # env fallback when exe is None
    m.exe = None
    os.environ['PHENIX_ORCA'] = '/env/orca'
    assert m.get_cmd() == '/env/orca orca_cmd.in', m.get_cmd()
  finally:
    os.environ.pop('PHENIX_ORCA', None)
    if saved is not None:
      os.environ['PHENIX_ORCA'] = saved
  print('  get_cmd OK')


def drive_hessian(m, **kwds):
  """Run run_hessian with run_cmd stubbed (ORCA never launched): the stub writes
  a placeholder .hess so the method succeeds. Return (written .in, path)."""
  def fake_run_cmd(log=None):
    with open(m.get_hessian_filename(), 'w') as f:  # ORCA would write this
      f.write('$hessian\n')
  m.run_cmd = fake_run_cmd
  path = m.run_hessian(**kwds)
  with open(m.get_input_filename()) as f:
    written = f.read()
  return written, path


def exercise_hessian():
  """run_hessian writes a single NumFreq .in and returns the .hess path.

  Everything runs in a temp dir; run_cmd is stubbed so ORCA is never launched.
  """
  cwd = os.getcwd()
  tmp = tempfile.mkdtemp(prefix='tst_orca_')
  os.chdir(tmp)
  try:
    # partial Hessian: indices sorted + de-duplicated; one monolithic .in
    m = make_manager('hess')
    assert m.get_hessian_filename() == 'orca_hess.hess'
    inp, path = drive_hessian(m, hessian_atoms=[1, 0, 1])
    assert 'NumFreq' in inp, inp
    assert 'Partial_Hess {0 1}' in inp, inp
    assert '%maxcore 2048' in inp, inp
    assert '%pal' not in inp, inp
    assert '* xyz 2 1' in inp, inp                # charge 2, multiplicity 1
    assert path == os.path.abspath(m.get_hessian_filename()), path
    assert os.path.exists(path)

    # full Hessian: no %freq block
    full, _ = drive_hessian(make_manager('hfull'))
    assert 'NumFreq' in full and 'Partial_Hess' not in full, full

    # partial=False ignores hessian_atoms
    nofreq, _ = drive_hessian(make_manager('hno'), hessian_atoms=[0, 1], partial=False)
    assert 'Partial_Hess' not in nofreq, nofreq

    # nproc>1 -> %pal in the Hessian header too
    inp4, _ = drive_hessian(make_manager('h4', nproc=4), hessian_atoms=[0])
    assert '%pal nprocs 4 end' in inp4, inp4

    # sites_cart overrides the geometry
    sc = flex.vec3_double([(1, 2, 3)] * 3)
    inpsc, _ = drive_hessian(make_manager('hsc'), sites_cart=sc)
    assert '1.00000 2.00000 3.00000' in inpsc, inpsc

    # a run that produces no .hess raises Sorry
    mf = make_manager('hfail')
    mf.run_cmd = lambda log=None: None
    try:
      mf.run_hessian()
      raise AssertionError('expected Sorry when no .hess is produced')
    except Sorry:
      pass
  finally:
    os.chdir(cwd)
    shutil.rmtree(tmp, ignore_errors=True)
  print('  hessian OK')


def run():
  exercise_parallel_and_optimisation()
  exercise_geom_convergence()
  exercise_get_cmd()
  exercise_hessian()


if __name__ == '__main__':
  run()
  print('OK')
