from __future__ import absolute_import, division, print_function
import os
from io import StringIO
import tempfile
import time

from libtbx.utils import Sorry
from scitbx.array_family import flex
from libtbx import adopt_init_args
from libtbx import easy_run

def dist2(xyz1, xyz2):
  d2=0
  for i in range(3):
    d2 += (xyz2[i]-xyz1[i])**2
  return d2

def get_internal_coordinate_value(atom1, atom2, atom3=None, atom4=None):
  import math
  from cctbx import geometry_restraints
  if not atom3:
    d2 = dist2(atom1.xyz, atom2.xyz)
    return math.sqrt(d2)
  elif not atom4:
    sites = [atom1.xyz, atom2.xyz, atom3.xyz]
    ang = geometry_restraints.angle(
      sites=sites,
      angle_ideal=109,
      weight=1,
      )
    return ang.angle_model
  sites = [atom1.xyz, atom2.xyz, atom3.xyz, atom4.xyz]
  dih = geometry_restraints.dihedral(
      sites=sites,
      angle_ideal=0,
      weight=1.,
      periodicity=1)
  return dih.angle_model

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
                        verbose=False,
                        ):
  if log_filename is not None: generator=loop_over_file(log_filename)
  if error_lines is None:
    error_lines = {
      '* GEOMETRY IN ERROR.' : 'Check protonated ligand and/or protein pocket',
      }
  error_line = None
  status = None
  errors=None
  for i, line in enumerate(generator):
    if line.find('GEOMETRY OPTIMIZATION CYCLE')>-1:
      cycle = int(line.split()[4])
      if verbose and cycle==1: print('  QM minimisation started', file=log)
    # if line.find('Max(Improp)')>-1:
    #   conv = process_orca_convergence(last_ten)
      # if cycle%10==0:
        # print('  Opt. cycle %s %s %s %0.1fmin' % (cycle, conv, i, (time.time()-t0)/60))
    if line.find('FINAL HEAT OF FORMATION =')>-1:
      status = True
    if line.find('ORCA TERMINATED NORMALLY')>-1:
      status = True
    if line.find('* JOB ENDED NORMALLY *')>-1:
      status = True
    if error_lines:
      for el, ad in error_lines.items():
        if line.find(el)>-1:
          error_line = line
          advice = ad
          break
    if error_line: break
  if error_line:
    raise Sorry(f'{error_line}\nAdvice: {advice}')
  if not status:
    raise Sorry('QM does not seem to have converged. Check %s' % log_filename)
  return status, errors

def run_command(command):
    """
    execute <command> in a subprocess and check error code
    taken from ASE
    """
    from subprocess import Popen, PIPE
    if command == '':
        raise RuntimeError('no command for run_command :(')
    proc = Popen([command], shell=True, stderr=PIPE)
    proc.wait()
    exitcode = proc.returncode
    if exitcode != 0:
        error='%s exited with error code %i in %s' % (
                       command,exitcode,self.calc_dir)
        stdout,stderr = proc.communicate()
        print('shell output: ',stdout,stderr)
        raise RuntimeError(error)
    return 0

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
    for i, atom in enumerate(self.atoms):
      ann=''
      if self.ligand_atoms_array: ann=self.ligand_atoms_array[i]
      outl += '  %s %s\n' % (atom.quote(), ann)
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

  def get_energy(self, *args, **kwds): return 0, 'dirac'

  def read_energy(self, *args, **kwds): return 0, 'dirac'

  def read_charge(self, *args, **kwds): return 99

  def get_strain(self, *args, **kwds): return 0, 'dirac'

  def get_bound(self, *args, **kwds): return 0, 'dirac'

  def get_opt(self, *args, **kwds):
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

  def get_gradients(self, *args, **kwds):
    return []

  def get_timings(self, energy=False):
    return '-'

  def get_lines(self, filename=None):
    if filename is None:
      filename = self.get_log_filename()
    assert os.path.exists(filename), 'filename not found %s' % filename
    f=open(filename, 'r')
    lines=f.read()
    del f
    return lines

class base_qm_manager(base_manager):

  def write_input(self, outl):
    f=open(self.get_input_filename(), 'w')
    f.write(outl)
    del f

  def check_file_read_safe(self, optimise_ligand=True, optimise_h=True, constrain_torsions=False, gradients=False):
    outl = self.get_input_lines(optimise_ligand=optimise_ligand,
                                optimise_h=optimise_h,
                                constrain_torsions=constrain_torsions,
                                gradients=gradients,
                                )
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
        f, ext = os.path.splitext(filename)
        f=open('old%s' % ext, 'w')
        f.write(lines)
        del f
        f=open('new%s' % ext, 'w')
        f.write(outl)
        del f
        raise Sorry('something has changed making the QM input files different')
    return outl==lines

  def opt_setup(self, optimise_ligand=True, optimise_h=True, constrain_torsions=False):
    outl = self.get_input_lines(optimise_ligand=optimise_ligand,
                                optimise_h=optimise_h,
                                constrain_torsions=constrain_torsions,
                                )
    self.write_input(outl)

  def get_opt(self,
              optimise_h=True,
              constrain_torsions=False,
              cleanup=False,
              file_read=True,
              check_file_read_safe=True,
              coordinate_filename_ext='.xyz',
              log_filename_ext='.log',
              redirect_output=True,
              log=None,
              verbose=False):
    if self.program_goal in ['opt', 'strain']:
      optimise_ligand=True
    # elif self.program_goal in ['energy']:
    #   optimise_ligand=False
    constrain_torsions = self.exclude_torsions_from_optimisation

    coordinates = None
    rc=True
    if file_read and check_file_read_safe:
      if verbose: print('check_file_read_safe',check_file_read_safe)
      rc = self.check_file_read_safe(optimise_ligand=optimise_ligand,
                                     optimise_h=optimise_h,
                                     constrain_torsions=constrain_torsions,
                                     )
      if verbose: print('rc',rc)
    if file_read and rc:
      filename = self.get_coordinate_filename()
      if os.path.exists(filename):
        lf = self.get_log_filename()
        if os.path.exists(lf):
          process_qm_log_file(lf, log=log)
        if verbose: print('  Reading coordinates from %s\n' % filename)
        coordinates = self.read_xyz_output()
    if coordinates is None:
      if verbose: print('no coordinates')
      self.opt_setup(optimise_ligand=optimise_ligand,
                     optimise_h=optimise_h,
                     constrain_torsions=constrain_torsions,
                     )
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
                 constrain_torsions=False,
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
                                               optimise_h=optimise_h,
                                               constrain_torsions=constrain_torsions,
                                               ):
      filename = self.get_log_filename()
      if os.path.exists(filename):
        if os.path.exists(filename):
          process_qm_log_file(filename, log=log)
        # print('  Reading energy from %s\n' % filename, file=log)
        energy, units = self.read_energy()
    if energy is None:
      outl = self.get_input_lines(optimise_ligand=optimise_ligand,
                                  optimise_h=optimise_h,
                                  constrain_torsions=constrain_torsions,
                                  )
      self.write_input(outl)
      self.run_cmd(redirect_output=redirect_output)
      energy, units = self.read_energy()
    if cleanup: self.cleanup(level=cleanup)
    # print('  Current energy = %0.5f %s' % (self.energy, self.units), file=log)
    self.preamble = old_preamble
    return energy, units

  def get_strain(self,
                 cleanup=False,
                 file_read=True,
                 redirect_output=False,
                 log=StringIO(),
                 **kwds):
    #
    # Get the strain energy of a ligand. Only works? when this is applied to
    # the ligand.
    #
    old_preamble = self.preamble
    #
    # get energy with just the H opt
    #
    start_energy, junk = self.get_energy(optimise_h=True,
                                         redirect_output=redirect_output,
                                         cleanup=cleanup,
                                         file_read=file_read,
                                         log=log)
    self.preamble = old_preamble+'_strain'
    #
    # optimise all atoms positions and get energy
    #
    final_energy, units = self.get_opt(redirect_output=redirect_output,
                                       cleanup=cleanup,
                                       file_read=file_read,
                                       log=log)
    final_energy, units = self.read_energy()
    self.strain = start_energy-final_energy
    self.units = units
    # print('  Strain energy = %0.5f %s' % (self.strain, self.units), file=log)
    self.preamble = old_preamble
    return self.strain, self.units

  def get_something_energy( self,
                            preamble,
                            cleanup=False,
                            file_read=True,
                            redirect_output=False,
                            log=StringIO(),
                            verbose=False,
                            **kwds
                            ):
    if verbose:
      print('get %s' % preamble)
      print(self)
    old_preamble = self.preamble
    self.preamble += '_%s' % preamble
    energy, units = self.get_energy(optimise_h=True,
                                    redirect_output=redirect_output,
                                    cleanup=cleanup,
                                    file_read=file_read,
                                    log=log)
    self.preamble = old_preamble
    return energy, units

  def get_bound(self, **kwds):
    return self.get_something_energy('bound', **kwds)

  def get_pocket(self, **kwds):
    return self.get_something_energy('pocket', **kwds)

  def get_gradients(self,
                    cleanup=False,
                    file_read=True,
                    redirect_output=False,
                    log=StringIO(),
                    verbose=False,
                    **kwds):
    gradients=None
    old_preamble = self.preamble
    self.preamble += '_gradients'
    optimise_ligand=False
    filename = self.get_log_filename()
    if file_read and self.check_file_read_safe(optimise_ligand=optimise_ligand,
                                               optimise_h=False,
                                               # constrain_torsions=constrain_torsions,
                                               gradients=True,
                                               ):
      filename = self.get_log_filename()
      if os.path.exists(filename):
        if os.path.exists(filename):
          process_qm_log_file(filename, log=log)
        # print('  Reading energy and gradients from %s\n' % filename, file=log)
        energy, gradients = self.read_engrad_output()
    if gradients is None:
      outl = self.get_input_lines(optimise_ligand=optimise_ligand,
                                  optimise_h=False,
                                  # constrain_torsions=constrain_torsions,
                                  gradients=True,
                                  )
      self.write_input(outl)
      self.run_cmd(redirect_output=redirect_output)
      energy, gradients = self.read_engrad_output()
    if cleanup: self.cleanup(level=cleanup)
    # print('  Current energy = %0.5f %s' % (self.energy, self.units), file=log)
    self.preamble = old_preamble
    return energy, gradients

  def get_timings(self, energy=None):
    return '-'
    if not self.times:
      filename = self.get_log_filename()
      if os.path.exists(filename):
        f=open(filename, 'r')
        lines=f.read()
        del f
        for line in lines.splitlines():
          #TOTAL JOB TIME:          1622.77 SECONDS
          if line.find('TOTAL JOB TIME:')>-1:
            print(line)
            # print(energy)
            # if energy is None:
            #   rc=self.read_energy()
            #   print(rc)
            return '  Timings : %0.2fs' % float(line.split()[3])
      # assert 0
      return '-'
    f='  Timings : %0.2fs (%ss)' % (
      self.times[-1],
      self.times.format_mean(format='%.2f'))
    if energy:
      f+=' Energy : %0.6f' % energy
    return f

  def _is_atom_for_opt(self, i, atom, optimise_ligand=True, optimise_h=True):
    ligand_atom = self.ligand_atoms_array[i]
    if optimise_ligand:
      if ligand_atom:
        opt=1
      else:
        opt=0
    else:
      opt=0
    if optimise_h and atom.element_is_hydrogen():
      opt=1
    return opt

  def guess_bonds(self):
    bonds = []
    for i, atom1 in enumerate(self.atoms):
      bonds.append([])
      for j, atom2 in enumerate(self.atoms):
        if i==j: continue
        d2 = dist2(atom1.xyz, atom2.xyz)
        d2_limit=2.5
        if atom1.element_is_hydrogen() and atom2.element_is_hydrogen():
          continue
        elif atom1.element_is_hydrogen() or atom2.element_is_hydrogen():
          d2_limit=1.3
        elif atom1.element in ['P'] or atom2.element in ['P']:
          d2_limit=4
        if d2<d2_limit:
          bonds[i].append(j)
    # for i, atom1 in enumerate(self.atoms):
    #   print(i, atom1.quote())
    #   for j in bonds[i]:
    #     print('  - %s' % self.atoms[j].quote())
    return bonds

  def get_torsion(self, i):
    if not hasattr(self, 'bonds'):
      self.bonds = self.guess_bonds()
    if self.atoms[i].parent().resname in ['HOH']: return None
    rc = [i]
    next = i
    j=0
    while len(rc)<4:
      next_i=self.bonds[i][j]
      if next_i not in rc:
        rc.append(next_i)
        i=next_i
        j=0
      else:
        j+=1
    return rc

def main():
  # test different run methods
  from mmtbx.geometry_restraints import mopac_manager
  dat = ''' AM1 XYZ GEO-OK
 HOH.dat

 O -23.081267 1 18.356610 1 -21.628018 1
 H -22.165708 1 18.441448 1 -21.593653 1
 H -23.452909 1 18.384220 1 -20.699162 1

'''
  f=open('HOH.dat', 'w')
  f.write(dat)
  del f
  cmd = '%s HOH.dat' % mopac_manager.get_exe()
  print(cmd)
  for i, func in enumerate([os.system, run_command, run_qm_cmd]):
    print(func)
    if i==2:
      func(cmd, 'HOH.out', redirect_output=False)
    else:
      func(cmd)
    f=open('HOH.out', 'r')
    lines=f.read()
    del f
    print(lines[-150:])

if __name__ == '__main__':
  main()
