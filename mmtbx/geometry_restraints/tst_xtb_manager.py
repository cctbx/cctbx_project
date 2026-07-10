"""Tests for mmtbx.geometry_restraints.xtb_manager.

Two tiers:

  * exercise_input_generation -- no xtb binary needed. Checks the xcontrol
    detailed input the manager writes for the opt / gradient / constrained-
    torsion modes, plus the option-parsing helpers and (when xtb is found) the
    assembled command line.

  * exercise_run -- needs a working xtb (located via xtb_manager.get_exe(), i.e.
    the build's bin dir or $PHENIX_XTB). Runs a tiny constrained optimisation, a
    single-point gradient and a torsion-constrained optimisation, and checks the
    energy / coordinate / gradient parsing. Skipped cleanly when xtb is absent.
"""
from __future__ import absolute_import, division, print_function
import os
import shutil
import tempfile

import iotbx.pdb
from cctbx import geometry_restraints
from libtbx.test_utils import approx_equal

from mmtbx.geometry_restraints import xtb_manager

# A water (donor O frozen, two H free) for opt / gradient / parsing.
WATER_PDB = """\
HETATM    1  O   HOH A   1      -0.000   0.000   0.000  1.00  0.00           O
HETATM    2  H1  HOH A   1       0.957   0.000   0.000  1.00  0.00           H
HETATM    3  H2  HOH A   1      -0.240   0.927   0.000  1.00  0.00           H
"""

# A propanol-like fragment with a real C-C-C-O torsion, for $constrain.
PROPANOL_PDB = """\
HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  C2  LIG A   1       1.520   0.000   0.000  1.00  0.00           C
HETATM    3  C3  LIG A   1       2.030   1.430   0.000  1.00  0.00           C
HETATM    4  O1  LIG A   1       3.450   1.450   0.300  1.00  0.00           O
HETATM    5  H1  LIG A   1      -0.380  -0.520   0.880  1.00  0.00           H
HETATM    6  H2  LIG A   1      -0.380  -0.520  -0.880  1.00  0.00           H
HETATM    7  H3  LIG A   1      -0.380   1.020   0.000  1.00  0.00           H
HETATM    8  H4  LIG A   1       1.880  -0.540   0.880  1.00  0.00           H
HETATM    9  H5  LIG A   1       1.880  -0.540  -0.880  1.00  0.00           H
HETATM   10  H6  LIG A   1       1.660   1.980   0.870  1.00  0.00           H
HETATM   11  H7  LIG A   1       1.660   1.980  -0.870  1.00  0.00           H
HETATM   12  H8  LIG A   1       3.800   2.350   0.300  1.00  0.00           H
"""


def make_manager(pdb_str, ligand_flags, preamble, method='GFN2-xTB',
                 solvent_model='ALPB=ether', charge=0, multiplicity=1, nproc=1):
  """Build an xtb_manager from a PDB string. The hierarchy is stashed on the
  manager so atom.parent() stays valid (get_torsion needs it)."""
  hier = iotbx.pdb.input(source_info=None,
                         lines=pdb_str.splitlines()).construct_hierarchy()
  atoms = hier.atoms()
  assert atoms.size() == len(ligand_flags)
  m = xtb_manager.xtb_manager(atoms, method, None, solvent_model, charge,
                              multiplicity, nproc, preamble=preamble)
  m.program_goal = 'opt'
  m.set_ligand_atoms(ligand_flags)
  m._tst_hier = hier  # keep the hierarchy (and atom.parent() chain) alive
  return m


def wrapped_dihedral(sites):
  d = geometry_restraints.dihedral(
    sites=sites, angle_ideal=0, weight=1, periodicity=1)
  return d.angle_model


def exercise_input_generation():
  """No xtb binary required: check the generated xcontrol input and helpers."""
  # --- option-parsing helpers ---
  assert xtb_manager._compress_index_ranges([1, 2, 3, 5, 6, 9]) == '1-3,5-6,9'
  assert xtb_manager._compress_index_ranges([]) == ''
  assert make_manager(WATER_PDB, [0, 1, 1], 'g2')._get_gfn() == 2
  assert make_manager(WATER_PDB, [0, 1, 1], 'g1', method='gfn1')._get_gfn() == 1
  assert make_manager(WATER_PDB, [0, 1, 1], 'gff', method='GFNFF')._get_gfn() == 'ff'
  assert make_manager(WATER_PDB, [0, 1, 1], 'gd', method=None)._get_gfn() == 2
  assert make_manager(WATER_PDB, [0, 1, 1], 's1')._get_solvent() == 'ether'
  assert make_manager(WATER_PDB, [0, 1, 1], 's2',
                      solvent_model='ether')._get_solvent() == 'ether'
  assert make_manager(WATER_PDB, [0, 1, 1], 's3',
                      solvent_model=None)._get_solvent() == ''

  m = make_manager(WATER_PDB, [0, 1, 1], 'inpgen')

  # --- optimisation input ---
  opt_lines = m.get_input_lines(optimise_ligand=True, optimise_h=True)
  for token in ('$chrg 0', '$spin 0', '$gfn', 'method=2', '$opt',
                'engine=inertial', 'maxcycle=1000', '$fix'):
    assert token in opt_lines, token
  # The donor O (atom 1, not a ligand atom and not H) is the only frozen atom.
  assert 'atoms: 1' in opt_lines, opt_lines

  # --- gradient input: a single point on the whole system, so no $opt/$fix ---
  grad_lines = m.get_input_lines(gradients=True)
  assert '$opt' not in grad_lines
  assert '$fix' not in grad_lines
  assert 'runtype=grad' in grad_lines
  assert '$gfn' in grad_lines

  # --- multiplicity -> unpaired-electron count ($spin = mult - 1) ---
  m5 = make_manager(WATER_PDB, [0, 1, 1], 'mult', multiplicity=5)
  assert '$spin 4' in m5.get_input_lines()

  # --- constrained torsions ---
  mt = make_manager(PROPANOL_PDB, [1] * 12, 'tors')
  tors_lines = mt.get_input_lines(constrain_torsions=True)
  for token in ('$constrain', 'force constant=1.0', 'dihedral:', ',auto'):
    assert token in tors_lines, token
  assert '$constrain' not in mt.get_input_lines(constrain_torsions=False)

  # --- command line (only when xtb is actually present) ---
  if xtb_manager.get_exe():
    cmd = m.get_cmd()
    for token in ('--gfn 2', '--alpb ether', '--opt', '-P 1',
                  '--namespace xtb_inpgen'):
      assert token in cmd, token
    assert '--grad' in m.get_cmd(gradients_only=True)
    m.robust = True
    rcmd = m.get_cmd()
    for token in ('--opt loose', '--etemp 1500', '--iterations 500'):
      assert token in rcmd, token
    m.robust = False
  print('  input generation OK')


def exercise_run():
  """Needs a working xtb. Runs tiny calculations and checks the parsers."""
  # --- constrained optimisation (donor O frozen, both H free) ---
  m = make_manager(WATER_PDB, [0, 1, 1], 'tst_xtb_opt')
  coords, buf = m.get_opt(file_read=False, redirect_output=True)
  assert buf.size() == 3 and coords.size() == 2          # buffer all; coords=H only
  assert approx_equal(buf[0], (0, 0, 0), eps=1e-3)       # frozen O stayed put
  energy, units = m.read_energy()
  assert units == 'Hartree'
  assert -6.0 < energy < -4.0, energy                    # GFN2 water ~ -5.08 Eh
  assert m.read_charge() == 0

  # --- reuse path: a fresh manager with the same preamble reads the output
  #     back (exercises check_file_read_safe + the xtb termination check). ---
  m_reuse = make_manager(WATER_PDB, [0, 1, 1], 'tst_xtb_opt')
  coords2, buf2 = m_reuse.get_opt(file_read=True, redirect_output=True)
  assert buf2.size() == 3
  assert approx_equal(list(buf2[1]), list(buf[1]), eps=1e-4)   # same H position

  # --- single-point gradient ---
  mg = make_manager(WATER_PDB, [0, 1, 1], 'tst_xtb_grad')
  mg.program_goal = 'gradients'
  ge, grad = mg.get_gradients(file_read=False)
  assert grad.size() == 3
  assert -6.0 < ge < -4.0, ge
  gmax = max(abs(c) for vec in grad for c in vec)
  assert gmax > 1e-2, gmax                               # genuinely non-zero

  # --- torsion-constrained optimisation ---
  mt = make_manager(PROPANOL_PDB, [1] * 12, 'tst_xtb_tors')
  mt.exclude_torsions_from_optimisation = True           # -> constrain_torsions
  start = [tuple(a.xyz) for a in mt.atoms]
  c_t, buf_t = mt.get_opt(file_read=False, redirect_output=True)
  assert buf_t.size() == 12
  log = mt.get_lines(mt.get_log_filename())
  assert 'constraining angle 1 2 3 4' in log, 'dihedral 1-2-3-4 not constrained'
  # The C1-C2-C3-O1 dihedral (atoms 0,1,2,3) should be held near its start value.
  before = wrapped_dihedral([start[0], start[1], start[2], start[3]])
  after = wrapped_dihedral([buf_t[0], buf_t[1], buf_t[2], buf_t[3]])
  drift = ((after - before + 180) % 360) - 180
  assert abs(drift) < 10.0, (before, after, drift)
  print('  run OK')


def run():
  exercise_input_generation()
  if xtb_manager.get_exe():
    cwd = os.getcwd()
    tmp = tempfile.mkdtemp(prefix='tst_xtb_')
    os.chdir(tmp)
    try:
      exercise_run()
    finally:
      os.chdir(cwd)
      shutil.rmtree(tmp, ignore_errors=True)
  else:
    print('  skipping run tests, xtb not found (set PHENIX_XTB)')


if __name__ == '__main__':
  run()
  print('OK')
