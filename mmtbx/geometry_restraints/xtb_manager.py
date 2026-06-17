"""QM manager for GFN-xTB constrained geometry optimisations.

Plugs into the ``base_qm_manager`` machinery (``get_opt`` / ``get_energy`` /
``get_gradients`` via ``write_input`` -> ``run_cmd`` -> ``read_*_output``), so
the quantum-restraints pipeline can drive xTB.

The tuning here targets transition-metal sites:

  * GFN2-xTB by default; the open shell is set by the unpaired-electron count
    (``$spin`` = multiplicity - 1), which fixes the orbital occupations in an
    otherwise spin-restricted run (no spin polarisation -- that needs xTB's
    opt-in ``--spinpol``, which requires the tblite backend and which this
    manager does not request).
  * ALPB implicit solvent by *named* solvent (default ``ether``, eps ~ 4.3,
    the closest common match to the eps = 4 protein-interior model).
  * The lbfgs optimiser (``$opt engine=lbfgs``) so a ``$fix`` block holds the
    frozen scaffold atoms *exactly* in Cartesian space on both of xTB's
    backends -- otherwise drift in the fixed atoms injects noise into the
    optimised region's geometry. The default ``rf`` optimiser drifts on both
    backends, and ``inertial`` (FIRE) -- exact on the native backend -- silently
    drops ``$fix`` on the tblite backend (which ``--spinpol`` forces); lbfgs is
    the only engine exact on both.
  * A generous geometry-optimisation cycle cap (xTB's automatic default is
    only ~2x the atom count, which a slow-converging case can exceed); xTB
    cycles are cheap.
  * A robustness preset (``robust=True``): a gentler optimiser level plus
    Fermi smearing and extra SCC iterations for open-shell clusters that
    otherwise crash/stall the SCF. NOTE its looser gradient leaves more
    residual play in bond lengths, so use it only when the default (tight)
    run fails.

xTB takes its runtype (opt vs gradient), GFN level and solvent on the command
line, while the detailed-input (xcontrol) file holds charge / spin / optimiser
/ frozen atoms. The detailed-input file is the one ``base_qm_manager`` tracks
as ``get_input_filename`` (so ``check_file_read_safe`` works on it); the actual
``.xyz`` xTB reads is written alongside, and the coordinates are also echoed as
``#`` comments into the detailed input so a geometry change is still caught by
the input diff. A per-job ``--namespace`` prefixes every xTB output file with
the preamble, so many (charge, spin, site) jobs can share a directory without
collisions.
"""
from __future__ import absolute_import, division, print_function
import os
import time
from itertools import groupby

from libtbx.utils import Sorry
from scitbx.array_family import flex
from libtbx import easy_run
from libtbx import Auto
import libtbx.load_env

from mmtbx.geometry_restraints import base_qm_manager

# Eh -> kcal/mol and bohr -> Angstrom, used to put gradients into kcal/mol/A.
harkcal = 627.50946900
bohrang = 0.52918

DEFAULT_XTB_GFN = 2
DEFAULT_XTB_SOLVENT = 'ether'
DEFAULT_XTB_MAXCYCLE = 1000
DEFAULT_XTB_THREADS = 1
# Harmonic force constant for $constrain dihedral restraints (Eh/rad^2). xTB has
# no exact dihedral constraint (only $fix is exact, and only on Cartesians), so
# torsions are held by a stiff restraint; 1.0 holds them to within ~0.2 deg.
DEFAULT_XTB_CONSTRAIN_FC = 1.0

# Robustness preset (applied with robust=True): a gentler optimiser level whose
# looser gradient lets a flat PES converge, plus Fermi smearing and extra SCC
# iterations for SCF stability on open-shell clusters.
XTB_ROBUST_OPTLEVEL = 'loose'
XTB_ROBUST_ETEMP = 1500
XTB_ROBUST_ITERATIONS = 500


def get_exe(verbose=False):
  """Locate the xtb executable: ``PHENIX_XTB`` env var, then the build's
  ``bin`` dir, then a Windows ``Library/bin/xtb.exe``. Returns the path or
  False."""
  env_exe = os.environ.get('PHENIX_XTB', False)
  bin_dir = libtbx.env.under_base('bin')
  exe_path = os.path.join(bin_dir, 'xtb')
  win_exe_path = os.path.join(libtbx.env.under_base('Library'), 'bin', 'xtb.exe')
  if verbose:
    print('"'*79)
    print(f'  Looking for xtb in env var PHENIX_XTB : {env_exe}')
    print(f'  Looking for xtb as {exe_path} : {os.path.exists(exe_path)}')
    print(f'  Looking for win xtb as {win_exe_path} : {os.path.exists(win_exe_path)}')
    print('"'*79)
  if env_exe:
    return env_exe
  for path in (exe_path, win_exe_path):
    if os.path.exists(path):
      return path
  return False


def _compress_index_ranges(indices):
  """[1,2,3,5,6] -> '1-3,5-6' (compact atom list for an xTB ``$fix`` block)."""
  # Consecutive integers share the same (value - position) key, so grouping on
  # it collects each contiguous run.
  parts = []
  for _, group in groupby(enumerate(sorted(set(indices))),
                          key=lambda pair: pair[1] - pair[0]):
    run = [value for _, value in group]
    parts.append(f'{run[0]}-{run[-1]}' if run[0] != run[-1] else f'{run[0]}')
  return ','.join(parts)


class xtb_manager(base_qm_manager.base_qm_manager):

  # Caller-overridable knobs.
  maxcycle = DEFAULT_XTB_MAXCYCLE
  threads = DEFAULT_XTB_THREADS
  robust = False

  error_lines = [
    'abnormal termination of xtb',
    '[ERROR]',
    '#ERROR!',
  ]

  # ---- filenames -----------------------------------------------------------

  def get_namespace(self):
    # Prefixes every xTB output file, so concurrent jobs in one dir don't clash.
    return f'xtb_{self.preamble}'

  def get_input_filename(self):
    # The xcontrol detailed-input file; this is what base_qm_manager tracks.
    return f'xtb_{self.preamble}.inp'

  def get_coordinate_input_filename(self):
    # The .xyz xTB reads (the geometry to optimise).
    return f'xtb_{self.preamble}.xyz'

  def get_coordinate_filename(self):
    # The optimised geometry xTB writes (namespaced xtbopt.xyz).
    return f'xtb_{self.preamble}.xtbopt.xyz'

  def get_gradient_filename(self):
    # Turbomole-format gradient file from a --grad run (namespaced).
    return f'xtb_{self.preamble}.gradient'

  def get_log_filename(self):
    return f'xtb_{self.preamble}.log'

  # ---- option helpers ------------------------------------------------------

  def _get_gfn(self):
    """GFN parametrisation from ``self.method``. Accepts forms like 'GFN2-xTB',
    'gfn2', '2' (-> int 2) or a GFN-FF request ('gfnff'/'gfn-ff' -> 'ff');
    defaults to ``DEFAULT_XTB_GFN``."""
    method = self.method
    if method in [None, Auto]:
      return DEFAULT_XTB_GFN
    method_str = str(method).lower()
    if 'ff' in method_str:
      return 'ff'
    digit = next((ch for ch in method_str if ch.isdigit()), None)
    return int(digit) if digit is not None else DEFAULT_XTB_GFN

  def _get_solvent(self):
    """ALPB solvent name from ``self.solvent_model``. Accepts a bare name
    ('ether'), or forms like 'ALPB=ether' / 'alpb ether'; empty disables
    solvation."""
    solvent_str = self.solvent_model
    if not solvent_str or solvent_str in [Auto]:
      return ''
    solvent_str = str(solvent_str).strip()
    if '=' in solvent_str:
      solvent_str = solvent_str.split('=')[-1].strip()
    elif ' ' in solvent_str:
      solvent_str = solvent_str.split()[-1].strip()
    return solvent_str

  def _get_threads(self):
    return self.nproc or self.threads or DEFAULT_XTB_THREADS

  # ---- input generation ----------------------------------------------------

  def _input_header(self, gradients_only=False):
    """xcontrol blocks that don't depend on the frozen-atom partition:
    charge/spin/GFN, and (for an optimisation) the lbfgs-optimiser ``$opt``
    block. A leading comment records the command-line-only settings (GFN,
    solvent, threads, robust) so a change in any of them is still detected by
    ``check_file_read_safe``."""
    gfn = self._get_gfn()
    solvent = self._get_solvent()
    charge = self.get_charge()  # base get_charge() already maps Auto/None -> 0
    n_unpaired = self.get_multiplicity() - 1  # xTB $spin = N(alpha) - N(beta)

    lines = [
      f'# xtb_manager detailed input ({self.preamble})',
      f'# gfn={gfn} solvent={solvent or "none"} threads={self._get_threads()} '
      f'robust={self.robust} runtype={"grad" if gradients_only else "opt"}',
      f'$chrg {int(charge)}',
      f'$spin {n_unpaired}',
    ]
    if gfn != 'ff':
      lines += ['$gfn', f'   method={gfn}', '$end']
    if not gradients_only:
      # lbfgs optimiser: the only engine that holds $fix *exactly* on both the
      # native and tblite backends. 'inertial' is exact on the native backend
      # but drifts the fixed atoms on tblite (which --spinpol forces); 'rf'
      # (default) drifts on both.
      lines += ['$opt', '   engine=lbfgs', f'   maxcycle={self.maxcycle}',
                '$end']
    return '\n'.join(lines) + '\n'

  def get_coordinate_lines(self, optimise_ligand=True, optimise_h=True,
                           constrain_torsions=False, include_fix=True):
    """The ``$fix`` block (every atom NOT selected for optimisation, frozen
    exactly) and, when ``constrain_torsions``, a ``$constrain`` block pinning
    the dihedral of each optimised heavy atom at its current value; followed by
    the geometry echoed as ``#`` comments. xTB ignores the comment lines, but
    including them means ``check_file_read_safe`` also detects a geometry change
    -- the coordinates xTB actually reads live in the separate .xyz the writer
    emits.

    ``include_fix=False`` suppresses the ``$fix`` block: a gradient evaluation
    wants forces on the whole system, and fixed atoms have their gradient
    components projected out.

    Torsions are restrained, not constrained exactly: xTB's only exact
    constraint is the Cartesian ``$fix``, so ``$constrain dihedral: ...,auto``
    holds each angle by a stiff harmonic restraint (force constant
    ``DEFAULT_XTB_CONSTRAIN_FC``).
    """
    frozen = []
    dihedrals = []
    if include_fix or constrain_torsions:
      for i, atom in enumerate(self.atoms):
        optimised = self._is_atom_for_opt(i, atom,
                                          optimise_ligand=optimise_ligand,
                                          optimise_h=optimise_h)
        if include_fix and not optimised:
          frozen.append(i + 1)  # $fix is 1-based
        if constrain_torsions and optimised and atom.element.strip() not in ('H', 'D'):
          try:
            torsion_atoms = self.get_torsion(i)  # 0-based [i,b,c,d], or None (HOH)
          except (IndexError, AttributeError):
            # get_torsion's greedy bond-walk can dead-end (e.g. at a terminal H)
            # before assembling a 4-atom chain; no well-defined torsion to pin.
            torsion_atoms = None
          if torsion_atoms:
            dihedral = tuple(j + 1 for j in torsion_atoms)  # $constrain is 1-based
            if dihedral not in dihedrals:
              dihedrals.append(dihedral)

    lines = []
    if frozen:
      lines += ['$fix', f'   atoms: {_compress_index_ranges(frozen)}', '$end']
    if dihedrals:
      lines += ['$constrain', f'   force constant={DEFAULT_XTB_CONSTRAIN_FC}']
      for d in dihedrals:
        lines.append(f'   dihedral: {",".join(str(j) for j in d)},auto')
      lines.append('$end')
    lines.append('# coordinates')
    for atom in self.atoms:
      x, y, z = atom.xyz
      lines.append(f'# {atom.element.strip():<2s} {x:14.6f} {y:14.6f} {z:14.6f}')
    return '\n'.join(lines) + '\n'

  def get_input_lines(self, optimise_ligand=True, optimise_h=True,
                      constrain_torsions=False, gradients=False):
    return self._input_header(gradients_only=gradients) + \
      self.get_coordinate_lines(optimise_ligand=optimise_ligand,
                                optimise_h=optimise_h,
                                constrain_torsions=constrain_torsions,
                                include_fix=not gradients,
                                )

  def get_xyz_lines(self):
    """The standard .xyz file xTB reads (full precision)."""
    lines = [str(len(self.atoms)), f'Generated by xtb_manager ({self.preamble})']
    for atom in self.atoms:
      elem = atom.element.strip()
      if not elem:
        # A blank element column would be silently misread; fail loudly.
        raise Sorry(f'atom {atom.quote()} has blank element column; '
                    'cannot write xTB xyz block')
      x, y, z = atom.xyz
      lines.append(f'{elem:<2s} {x:14.6f} {y:14.6f} {z:14.6f}')
    return '\n'.join(lines) + '\n'

  def write_input(self, content):
    # Write the detailed-input file base_qm_manager tracks, then materialise
    # the .xyz xTB actually reads. (_gradients_only, recovered from whether an
    # $opt block was emitted, tells run_cmd which runtype flag to use.)
    self._gradients_only = '$opt' not in content
    with open(self.get_input_filename(), 'w') as f:
      f.write(content)
    with open(self.get_coordinate_input_filename(), 'w') as f:
      f.write(self.get_xyz_lines())

  # ---- running -------------------------------------------------------------

  def get_cmd(self, gradients_only=False):
    exe = get_exe()
    if not exe:
      raise Sorry('xtb executable not found. Set the PHENIX_XTB environment '
                  'variable or install xtb in the build.')
    gfn = self._get_gfn()
    gfn_flag = '--gfnff' if gfn == 'ff' else f'--gfn {gfn}'
    solvent = self._get_solvent()
    solvent_flag = f' --alpb {solvent}' if solvent else ''
    threads = self._get_threads()
    if gradients_only:
      runtype = '--grad'
    elif self.robust:
      runtype = (f'--opt {XTB_ROBUST_OPTLEVEL} --etemp {XTB_ROBUST_ETEMP} '
                 f'--iterations {XTB_ROBUST_ITERATIONS}')
    else:
      runtype = '--opt'
    return (f'{exe} {self.get_coordinate_input_filename()} '
            f'--input {self.get_input_filename()} '
            f'--namespace {self.get_namespace()} '
            f'{gfn_flag}{solvent_flag} -P {threads} {runtype}')

  def _omp_env(self):
    # OpenMP settings applied to the run environment: pin the thread count (xTB
    # scales poorly on small clusters, so default 1 thread maximises throughput
    # when packing many jobs) and set xTB's recommended stack settings.
    threads = self._get_threads()
    return {
      'OMP_NUM_THREADS': f'{threads},1',
      'MKL_NUM_THREADS': str(threads),
      'OMP_STACKSIZE': '4G',
      'OMP_MAX_ACTIVE_LEVELS': '1',
    }

  def run_cmd(self, redirect_output=True, log=None):
    # xTB's success/error strings aren't in base_qm_manager.run_qm_cmd, so do a
    # self-contained run + check here. xTB output is always sent to the log file
    # because read_energy/read_charge parse it. The OpenMP environment is set on
    # os.environ (easy_run.go execs the command, so an inline VAR=val prefix
    # would be misread as the program name).
    t0 = time.time()
    cmd = self.get_cmd(gradients_only=getattr(self, '_gradients_only', False))
    log_filename = self.get_log_filename()
    full_cmd = f'{cmd} >& {log_filename}'
    print(f'  Starting : {full_cmd}', file=log)
    saved = {k: os.environ.get(k) for k in self._omp_env()}
    os.environ.update(self._omp_env())
    try:
      easy_run.go(full_cmd)
    finally:
      for k, v in saved.items():
        if v is None:
          os.environ.pop(k, None)
        else:
          os.environ[k] = v
    if not os.path.exists(log_filename):
      raise Sorry(f'xtb did not appear to write log file : {log_filename}')
    self._check_log(log_filename, log=log)
    self.times.append(time.time() - t0)

  def _check_log(self, log_filename, log=None):
    lines = self.get_lines(filename=log_filename)
    for marker in self.error_lines:
      if marker in lines:
        raise Sorry(f'xtb error in {log_filename}:\n  {marker}')
    # 'abnormal termination of xtb' contains 'normal termination of xtb' as a
    # substring, so it's caught by error_lines above before this success test.
    if 'normal termination of xtb' not in lines:
      raise Sorry(f'xtb does not seem to have terminated normally. '
                  f'Check {log_filename}')

  # ---- reading outputs -----------------------------------------------------

  def read_xyz_output(self):
    filename = self.get_coordinate_filename()
    if not os.path.exists(filename):
      raise Sorry(f'QM output filename not found: {filename}')
    with open(filename, 'r') as f:
      lines = f.read().splitlines()[2:]  # skip the atom-count and comment lines
    coords = flex.vec3_double()
    for line in lines:
      fields = line.split()
      if len(fields) >= 4:
        coords.append(tuple(float(v) for v in fields[1:4]))
    return coords

  def read_energy(self, filename=None):
    # |  TOTAL ENERGY    -5.079025732220 Eh  |  -- last occurrence is final.
    lines = self.get_lines(filename=filename)
    for line in lines.splitlines():
      if 'TOTAL ENERGY' in line:
        self.energy = float(line.replace('|', '').split()[2])
        self.units = 'Hartree'
    return self.energy, self.units

  def read_charge(self, filename=None):
    #  :  net charge    0  :
    lines = self.get_lines(filename=filename)
    for line in lines.splitlines():
      if 'net charge' in line:
        try:
          self.charge = int(round(float(line.replace(':', '').split()[2])))
        except (IndexError, ValueError):
          pass
    return self.charge

  def read_engrad_output(self, filename=None):
    """Parse the Turbomole-format ``gradient`` file from a --grad run. Returns
    (energy in Hartree, gradients as flex.vec3_double in kcal/mol/A)."""
    if filename is None:
      filename = self.get_gradient_filename()
    if not os.path.exists(filename):
      raise Sorry(f'xtb gradient file not found: {filename}')
    natoms = len(self.atoms)
    energy = None
    grad_rows = []
    with open(filename, 'r') as f:
      for line in f:
        stripped = line.strip()
        if not stripped or stripped.startswith('$'):
          continue
        if 'SCF energy' in stripped:
          # cycle = 1  SCF energy = -5.0788...  |dE/dxyz| = 0.017662
          tokens = stripped.replace('=', ' ').split()
          energy = float(tokens[tokens.index('energy') + 1])
          continue
        fields = stripped.split()
        # Coordinate lines carry a trailing element (4 tokens); gradient lines
        # are exactly 3 floats (Fortran D/E exponent notation).
        if len(fields) == 3:
          try:
            grad_rows.append([float(x.replace('D', 'E')) for x in fields])
          except ValueError:
            continue
    gradients = flex.vec3_double()
    for vals in grad_rows[-natoms:]:
      gradients.append([v * harkcal * bohrang for v in vals])
    self.energy = energy
    self.gradients = gradients
    return self.energy, self.gradients

  # ---- cleanup -------------------------------------------------------------

  def cleanup(self, level=None, verbose=False):
    if not self.preamble or level is None:
      return
    namespace = self.get_namespace()
    # Scratch/restart files with no downstream value.
    scratch_suffixes = ('.xtbrestart', '.xtbtopo.mol', '.wbo', '.charges',
                        '.engrad')
    keepers = (self.get_coordinate_filename(),  # optimised geometry
               self.get_log_filename(),
               self.get_input_filename(),
               self.get_coordinate_input_filename(),
               self.get_gradient_filename())
    for filename in os.listdir('.'):
      if not filename.startswith(namespace):
        continue
      if level == 'all':
        remove = filename not in keepers
      elif level == 'most':
        remove = filename.endswith(scratch_suffixes)
      else:
        remove = False
      if remove:
        if verbose:
          print(f'  removing {filename}')
        os.remove(filename)


if __name__ == '__main__':
  exe = get_exe(verbose=True)
  print(f'xtb executable {exe}')
