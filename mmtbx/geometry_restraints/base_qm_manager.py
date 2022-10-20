from __future__ import absolute_import, division, print_function
import os
from io import StringIO
import tempfile
import time

from libtbx.utils import Sorry
from scitbx.array_family import flex
from libtbx import adopt_init_args
from libtbx import easy_run

def loop_over_file(filename):
  f=open(filename, 'r')
  lines = f.read()
  del f
  for line in lines.splitlines():
    yield line

def process_qm_log_file(log_filename=None,
                        generator=None,
                        error_lines=None,
                        log=None,
                        ):
  if log_filename is not None: generator=loop_over_file(log_filename)
  error_line = None
  status = None
  for i, line in enumerate(generator):
    if line.find('GEOMETRY OPTIMIZATION CYCLE')>-1:
      cycle = int(line.split()[4])
      if cycle==1: print('  QM minimisation started', file=log)
    # if line.find('Max(Improp)')>-1:
    #   conv = process_orca_convergence(last_ten)
      # if cycle%10==0:
        # print('  Opt. cycle %s %s %s %0.1fmin' % (cycle, conv, i, (time.time()-t0)/60))
    # if line.find('OPTIMIZATION RUN DONE')>-1:
    #   status = True
    if line.find('ORCA TERMINATED NORMALLY')>-1:
      status = True
    if line.find('* JOB ENDED NORMALLY *')>-1:
      status = True
    if error_lines:
      for el in error_lines:
        if line.find(el)>-1:
          error_line = line
          break
    if error_line: break
  if error_line:
    raise Sorry(error_line)
  if not status:
    raise Sorry('QM does not seem to have converged. Check %s' % log_filename)
  return status

def run_qm_cmd(cmd,
               log_filename,
               error_lines=None,
               redirect_output=True,
               log=None,
               verbose=False,
               ):
  def loop_over_list(l):
    for line in l:
      yield l
  import time
  t0=time.time()
  if verbose: print('run_qm_cmd',cmd,log_filename)
  if redirect_output:
    cmd += ' >& %s' % log_filename

  print('  Starting : %s' % cmd, file=log)
  rc = easy_run.go(cmd)
  if not os.path.exists(log_filename):
    raise Sorry('QM program did not appear to write log file : %s' % log_filename)

  generator = loop_over_file(log_filename)
  # if redirect_output:
  #   generator = loop_over_file(log_filename)
  # else:
  #   generator = loop_over_list(rc.stdout_lines)

  status = process_qm_log_file(generator=generator,
                               error_lines=error_lines,
                               log=log)
  if rc.stderr_lines:
    print('stderr')
    for line in rc.stderr_lines:
      print(line)
    print('stdout')
    for line in rc.stdout_lines:
      print(line)
    assert 0
  return rc

class base_manager():
  def __init__(self,
               atoms,
               method,
               basis_set,
               solvent_model,
               charge,
               multiplicity,
               nproc=1,
               preamble=None,
               log=None,
               ):
    adopt_init_args(self, locals())
    self.times = flex.double()
    self.energies = {}
    if self.preamble is None:
      self.preamble = os.path.basename(tempfile.NamedTemporaryFile().name)
    # validate
    if self.basis_set is None: self.basis_set=''
    if self.solvent_model is None: self.solvent_model=''
    #
    self.ligand_atoms_array = None
    self.error_lines = []

  def __repr__(self):
    program = getattr(self, 'program', '')
    if program: program = ' - %s' % program
    outl = 'QI Manager%s\n' % program
    outl += ' charge: %s multiplicity: %s\n method: %s basis: "%s" solvent: "%s"\n' % (
      self.charge,
      self.multiplicity,
      self.method,
      self.basis_set,
      self.solvent_model,
      )
    residues = []
    for atom in self.atoms:
      outl += '  %s\n' % atom.quote()
      if atom.parent().id_str() not in residues:
        residues.append(atom.parent().id_str())
    outl += '  residues\n'
    outl += '\n  '.join(residues)
    return outl

  def get_charge(self): return self.charge

  def set_charge(self, charge): self.charge = charge

  def add_atoms(self, atoms, replace=False):
    if replace:
      self.atoms=atoms
    else:
      quotes = []
      for atom in self.atoms:
        quotes.append(atom.quote())
      for atom in atoms:
        assert atom.quote() not in quotes, 'found duplicate %s' % atom.quote()
      self.atoms.append(atoms)

  def set_ligand_atoms(self, selection_array):
    """Set the atoms that are ligand or the optimisation selection

    Args:
        selection_array (array): Selection array that specifies the atoms that
        will be optimised.
    """
    assert len(selection_array)==len(self.atoms)
    self.ligand_atoms_array = selection_array

  # def set_frozen_atoms(self, selection_array):
  #   """Summary

  #   Args:
  #       selection_array (TYPE): Description
  #   """
  #   assert len(selection_array)==len(self.atoms)
  #   self.freeze_a_ray = selection_array

  def get_opt(self,
              cleanup=False,
              file_read=False,
              coordinate_filename_ext='.xyz',
              log_filename_ext='.log',
              redirect_output=True,
              log=None):
    import random
    rc = []
    for atom in self.atoms:
      rc.append([])
      for i in range(3):
        rc[-1].append(atom.xyz[i]+(random.random()-0.5)/10)
    rc_buffer = rc
    tmp = []
    if self.ligand_atoms_array is not None:
      for sel, atom in zip(self.ligand_atoms_array, rc):
        if sel:
          tmp.append(atom)
      rc=tmp
    return flex.vec3_double(rc), flex.vec3_double(rc_buffer)

  def get_timings(self, energy=False):
    return '-'

class base_qm_manager(base_manager):

  def write_input(self, outl):
    f=open(self.get_input_filename(), 'w')
    f.write(outl)
    del f

  def check_file_read_safe(self, optimise_ligand=True, optimise_h=True):
    outl = self.get_input_lines(optimise_ligand=optimise_ligand,
                                optimise_h=optimise_h)
    filename = self.get_input_filename()
    lines=''
    if os.path.exists(filename):
      f=open(filename, 'r')
      lines = f.read()
      del f
      if not(outl==lines):
        print('differences '*5)
        print('proposed')
        print(outl)
        print('='*80)
        print(filename)
        print(lines)
        print('filename',filename)
        print('='*80)
        for i, (line1, line2) in enumerate(zip(outl.splitlines(),lines.splitlines())):
          if line1!=line2: print(' ! %s "%s" <> "%s"' % (i+1, line1, line2))
        print('='*80)
        raise Sorry('something has changed making the QM input files different')
    return outl==lines

  def opt_setup(self, optimise_ligand=True, optimise_h=True):
    outl = self.get_input_lines(optimise_ligand=optimise_ligand,
                                optimise_h=optimise_h)
    self.write_input(outl)

  def get_opt(self,
              optimise_h=True,
              cleanup=False,
              file_read=True,
              check_file_read_safe=True,
              coordinate_filename_ext='.xyz',
              log_filename_ext='.log',
              redirect_output=True,
              log=None):
    if self.program_goal in ['opt', 'strain']:
      optimise_ligand=True
    # elif self.program_goal in ['energy']:
    #   optimise_ligand=False

    coordinates = None
    rc=True
    if check_file_read_safe:
      rc = self.check_file_read_safe(optimise_ligand=optimise_ligand,
                                     optimise_h=optimise_h)
    if file_read and rc:
      filename = self.get_coordinate_filename()
      if os.path.exists(filename):
        lf = self.get_log_filename()
        if os.path.exists(lf):
          process_qm_log_file(lf, log=log)
        print('  Reading coordinates from %s\n' % filename, file=log)
        coordinates = self.read_xyz_output()
    if coordinates is None:
      self.opt_setup(optimise_ligand=optimise_ligand, optimise_h=optimise_h)
      self.run_cmd(redirect_output=redirect_output, log=log)
      coordinates = self.read_xyz_output()
    coordinates_buffer = coordinates
    if cleanup: self.cleanup(level=cleanup)
    if self.ligand_atoms_array is not None:
      tmp = []
      for sel, atom in zip(self.ligand_atoms_array, coordinates):
        if sel:
          tmp.append(atom)
      coordinates=tmp
    return flex.vec3_double(coordinates), flex.vec3_double(coordinates_buffer)

  def get_energy(self,
                 optimise_h=True,
                 cleanup=False,
                 file_read=True,
                 redirect_output=False,
                 log=StringIO(),
                 **kwds):
    energy=None
    old_preamble = self.preamble
    self.preamble += '_energy'
    optimise_ligand=False
    if file_read and self.check_file_read_safe(optimise_ligand=optimise_ligand,
                                               optimise_h=optimise_h):
      filename = self.get_log_filename()
      if os.path.exists(filename):
        if os.path.exists(filename):
          process_qm_log_file(filename, log=log)
        print('  Reading energy from %s\n' % filename, file=log)
        energy = self.read_energy()
    if energy is None:
      outl = self.get_input_lines(optimise_ligand=optimise_ligand,
                                  optimise_h=optimise_h)
      self.write_input(outl)
      self.run_cmd(redirect_output=redirect_output)
      energy = self.read_energy()
    if cleanup: self.cleanup(level=cleanup)
    print('  Current energy = %0.5f %s' % (self.energy, self.units), file=log)
    self.preamble = old_preamble
    return energy

  def get_strain(self,
                 cleanup=False,
                 file_read=True,
                 redirect_output=False,
                 log=StringIO(),
                 **kwds):
    old_preamble = self.preamble
    start_energy, junk = self.get_energy(optimise_h=True,
                                         redirect_output=redirect_output,
                                         cleanup=cleanup,
                                         log=log)
    self.preamble = old_preamble+'_strain'
    final_energy, junk = self.get_opt(redirect_output=redirect_output,
                                      cleanup=cleanup,
                                      log=log)
    final_energy, junk = self.read_energy()
    self.strain = start_energy-final_energy
    print('  Strain energy = %0.5f %s' % (self.strain, self.units), file=log)
    self.preamble = old_preamble
    return self.strain, self.units

  def get_timings(self, energy=None):
    if not self.times: return '-'
    f='  Timings : %0.2fs (%ss)' % (
      self.times[-1],
      self.times.format_mean(format='%.2f'))
    if energy:
      f+=' Energy : %0.6f' % energy
    return f

