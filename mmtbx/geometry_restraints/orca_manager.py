from __future__ import absolute_import, division, print_function
import os
import time

from libtbx.utils import Sorry
from scitbx.array_family import flex
from mmtbx.geometry_restraints import base_qm_manager

harkcal = 627.50946900
bohrang = 0.52918

class orca_manager(base_qm_manager.base_qm_manager):

  error_lines = [
                  'ORCA finished by error termination in GSTEP',
                  '-> impossible',
                  'SCF NOT CONVERGED AFTER',
                  'SERIOUS PROBLEM IN SOSCF',
                ]
  def get_input_filename(self): return 'orca_%s.in' % self.preamble

  def get_coordinate_filename(self): return 'orca_%s.xyz' % self.preamble

  def get_log_filename(self): return 'orca_%s.log' % self.preamble

  def read_engrad_output(self):
    '''#
# Number of atoms
#
 5
#
# The current total energy in Eh
#
    -49.737578240166
#
# The current gradient in Eh/bohr
#
       0.009609074575
       0.007643624367
      -0.019142934602
       0.010258288141
      -0.020537435105
      -0.000346851479
       0.000773577750
       0.021293697927
       0.011393000407
      -0.018928466970
      -0.006660132835
       0.008456622796
      -0.001712473496
      -0.001739754355
      -0.000359837122
#
# The atomic numbers and current coordinates in Bohr
#
   8    59.0407136   72.7582356   32.5750991
   8    57.8558553   75.8403789   29.3417777
   8    58.8800869   71.4618835   28.1663680
   8    62.2022254   74.3474953   29.5553167
  16    59.4829095   73.6048329   29.8973572'''
    f=open('orca_%s.engrad' % self.preamble, 'r')
    lines = f.read()
    del f
    lines = lines.split('#')

    energy = None
    gradients = flex.vec3_double()
    for line in lines[6].splitlines():
      if len(line.strip()):
        energy = float(line)
        break
    tmp=[]
    for line in lines[9].splitlines():
      if len(line.strip()):
        tmp.append(float(line)*harkcal*bohrang)
        if len(tmp)==3:
          gradients.append(tmp)
          tmp=[]

    self.energy = energy
    self.gradients = gradients
    return self.energy, self.gradients

  def read_charge(self):
    filename = self.get_log_filename()
    f=open(filename, 'r')
    lines=f.readlines()
    del f
    #Sum of atomic charges:   -1.0000000
    for line in lines:
      if line.find('Sum of atomic charges:')>-1:
        if len(line.split())==5:
          self.charge = float(line.split()[-1])
        else:
          self.charge = 99
    return self.charge

  def read_energy(self):
    filename = self.get_log_filename()
    f=open(filename, 'r')
    lines=f.readlines()
    del f
    for line in lines:
      if line.find('FINAL SINGLE POINT ENERGY')>-1:
        if len(line.split())==5:
          self.energy = float(line.split()[-1])
        else:
          self.energy = 1e-9
        self.units = 'Hartree'
    return self.energy, self.units

  def read_xyz_output(self):
    filename = self.get_coordinate_filename()
    if not os.path.exists(filename):
      raise Sorry('QM output filename not found: %s' % filename)
    f=open(filename, 'r')
    lines = f.read()
    del f
    rc = flex.vec3_double()
    for i, line in enumerate(lines.splitlines()):
      if i>=2:
        tmp = line.split()
        rc.append((float(tmp[1]), float(tmp[2]), float(tmp[3])))
    return rc

  def get_cmd(self):
    cmd = '%s orca_%s.in' % (
      os.environ['PHENIX_ORCA'],
      self.preamble,
      )
    return cmd

  def run_cmd(self, redirect_output=True, log=None):
    t0=time.time()
    cmd = self.get_cmd()
    assert redirect_output
    base_qm_manager.run_qm_cmd(cmd,
                               'orca_%s.log' % self.preamble,
                               error_lines=self.error_lines,
                               redirect_output=redirect_output,
                               log=log,
                               )
    self.times.append(time.time()-t0)

  def _input_header(self, gradients_only=False):
    standard_options = '''%scf maxiter=500

SOSCFStart 0.00033 # Default value of orbital gradient is 0.0033. Here reduced by a factor of 10.

end
'''
    addiotnal_options=''
    ptr=1
    if gradients_only:
      ptr=2
    outl = '%s\n! %s %s %s %s %s\n\n' % (standard_options,
                                       self.method,
                                       self.basis_set,
                                       self.solvent_model,
                                       ['Opt', 'LooseOpt', 'EnGrad'][ptr],
                                       addiotnal_options,
                                       )
    return outl

  def get_coordinate_lines(self, optimise_ligand=True, optimise_h=True, constrain_torsions=False):
    outl = '* xyz %s %s\n' % (self.charge, self.multiplicity)
    for i, atom in enumerate(self.atoms):
      outl += ' %s %0.5f %0.5f %0.5f # %s %s\n' % (
        atom.element,
        atom.xyz[0],
        atom.xyz[1],
        atom.xyz[2],
        atom.id_str(),
        i,
        )
    outl += '*\n'
    return outl

  def get_input_lines(self, optimise_ligand=True, optimise_h=True, constrain_torsions=False, gradients=False):
    outl = self._input_header(gradients_only=gradients)
    outl += self.get_coordinate_lines(optimise_ligand=optimise_ligand,
                                      optimise_h=optimise_h,
                                      constrain_torsions=constrain_torsions,
                                      )
    freeze_outl = '''%geom
      Constraints
'''
    added = 0
    for i, (sel, atom) in enumerate(zip(self.ligand_atoms_array, self.atoms)):
      if optimise_h and atom.element in ['H', 'D']: continue
      opt = self._is_atom_for_opt(i,
                                  atom,
                                  optimise_ligand=optimise_ligand,
                                  optimise_h=optimise_h)
      tmp=''
      if constrain_torsions:
        if opt and atom.element not in ['H', 'D']:
          torsions = self.get_torsion(i)
          tmp = '{D %d %d %d %d C } # Constraining dihedral' % tuple(torsions)
          tmp += ' "%s %s %s" :' % (atom.parent().parent().parent().id,
                                    atom.parent().parent().resseq.strip(),
                                    atom.parent().resname,
                                    )
          for j in torsions: tmp += ' %s' % self.atoms[j].name.strip()
          tmp += '\n'
      if tmp:
        freeze_outl += tmp
      if sel and optimise_ligand: continue
      freeze_outl += '{C %d C} # Constraining xyz %s\n' % (i, atom.id_str())
      added+=1
    freeze_outl += 'end\nend\n'
    if added: outl+=freeze_outl
    assert outl.find('Opt')>-1 or outl.find('EnGrad')>-1
    return outl

  # def get_engrad(self):
  #   outl = '! %s %s %s EnGrad\n\n' % (self.method,
  #                                     self.basis_set,
  #                                     self.solvent_model)
  #   outl += self.get_coordinate_lines()
  #   if outl in self.energies:
  #     self.times.append(0)
  #     return self.energies[outl]
  #   self.write_input(outl)
  #   self.run_cmd()
  #   energy, gradients = self.read_engrad_output()
  #   self.print_timings(energy)
  #   self.energies[outl] = (energy, gradients)
  #   return energy, gradients

  def cleanup(self, level=None, verbose=False):
    if not self.preamble: return
    if level is None: return
    #
    tf = 'orca_%s.trj' % self.preamble
    if os.path.exists(tf):
      uf = 'orca_%s_trj.xyz' % self.preamble
      # print('rename',tf,uf)
      os.rename(tf, uf)
    most_keepers = ['.xyz', '.log', '.in', '.engrad', '.trj']
    for filename in os.listdir('.'):
      if filename.startswith('orca_%s' % self.preamble):
        if level=='most':
          name, ext = os.path.splitext(filename)
          if ext in most_keepers: continue
        if verbose: print('  removing',filename)
        os.remove(filename)

  def view(self, cmd, ext='.xyz'):
    # /Applications/Avogadro.app/Contents/MacOS/Avogadro
    print(cmd)
    tf = 'orca_%s' % self.preamble
    print(tf)
    filenames =[]
    for filename in os.listdir('.'):
      if filename.startswith(tf) and filename.endswith(ext):
        filenames.append(filename)
    filenames.sort()
    print(filenames)
    cmd += ' %s' % filenames[-1]
    easy_run.go(cmd)
