from __future__ import absolute_import, division, print_function
import os
from io import StringIO
import time

from scitbx.array_family import flex
from mmtbx.geometry_restraints import base_qm_manager

class mopac_manager(base_qm_manager.base_qm_manager):
  def get_coordinate_filename(self):
    return 'mopac_%s.arc' % self.preamble

  def get_log_filename(self):
    return 'mopac_%s.out' % self.preamble

  def _input_header(self):
    if self.nproc==0:
      nproc_str=''
    else:
      nproc_str='THREADS=%s' % self.nproc
    outl = '%s %s %s %s \n%s\n\n' % (
     self.method,
     self.basis_set,
     self.solvent_model,
     'CHARGE=%s %s' % (self.charge, nproc_str),
     self.preamble,
     )
    return outl

  def opt_setup(self, optimise_ligand=True, optimise_h=True):
    outl = self._input_header()
    assert self.multiplicity==1
    outl += self.get_coordinate_lines(optimise_ligand=optimise_ligand,
                                      optimise_h=optimise_h)
    self.write_input(outl)

  def get_coordinate_lines(self, optimise_ligand=True, optimise_h=True):
    outl = ''
    for i, atom in enumerate(self.atoms):
      ligand_atom = self.ligand_atoms_array[i]
      if optimise_ligand:
        if ligand_atom:
          opt=1
        else:
          opt=0
      else:
        opt=0
      if optimise_h and atom.element in ['H', 'D']:
        opt=1
      outl += ' %s %0.5f %d %0.5f %d %0.5f %d\n' % (
        atom.element,
        atom.xyz[0],
        opt,
        atom.xyz[1],
        opt,
        atom.xyz[2],
        opt,
        # atom.id_str(),
        # i,
        )
    outl += '\n'
    return outl

  def get_input_lines(self, optimise_ligand=True, optimise_h=True):
    outl = self._input_header()
    outl += self.get_coordinate_lines(optimise_ligand=optimise_ligand,
                                      optimise_h=optimise_h)
    return outl

  def get_input_filename(self):
    return 'mopac_%s.mop' % self.preamble

  def write_input(self, outl):
    f=open(self.get_input_filename(), 'w')
    f.write(outl)
    del f

  def read_xyz_output(self):
    filename = self.get_coordinate_filename()
    if not os.path.exists(filename):
      raise Sorry('QM output filename not found: %s' % filename)
    f=open(filename, 'r')
    lines = f.read()
    del f
    rc = flex.vec3_double()
    read_xyz = False
    for i, line in enumerate(lines.splitlines()):
      if read_xyz:
        tmp = line.split()
        if len(tmp)==7:
          rc.append((float(tmp[1]), float(tmp[3]), float(tmp[5])))
      if line.find('FINAL GEOMETRY OBTAINED')>-1:
        read_xyz=True
    return rc

  def get_cmd(self):
    cmd = '%s mopac_%s' % (
      os.environ['PHENIX_MOPAC'],
      self.preamble,
      )
    return cmd

  def run_cmd(self, redirect_output=True, log=None):
    t0=time.time()
    cmd = self.get_cmd()
    base_qm_manager.run_qm_cmd(cmd,
                               'mopac_%s.out' % self.preamble,
                               error_lines=self.error_lines,
                               redirect_output=redirect_output,
                               log=log,
                               )
    self.times.append(time.time()-t0)

  def get_energy(self,
                 optimise_h=True,
                 cleanup=False,
                 file_read=True,
                 coordinate_filename_ext=None, # not used
                 log_filename_ext=None, # not used
                 redirect_output=False, # not used
                 log=StringIO()):
    energy=None
    old_preamble = self.preamble
    self.preamble += '_energy'
    optimise_ligand=False
    if file_read and self.check_file_read_safe(optimise_ligand=optimise_ligand,
                                               optimise_h=optimise_h):
      filename = self.get_log_filename()
      if os.path.exists(filename):
        if os.path.exists(filename):
          base_qm_manager.process_qm_log_file(filename, log=log)
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

  def read_energy(self):
    filename = self.get_log_filename()
    f=open(filename, 'r')
    lines=f.readlines()
    del f
    '''FINAL HEAT OF FORMATION =       -132.17152 KCAL/MOL =    -553.00562 KJ/MOL





          TOTAL ENERGY            =      -1356.82771 EV'''
    for line in lines:
      if line.find('FINAL HEAT OF FORMATION')>-1:
        self.energy = float(line.split()[5])
        self.units = line.split()[6]
    return self.energy, None

  def get_strain(self,
                 cleanup=False,
                 file_read=True,
                 coordinate_filename_ext=None, # not used
                 log_filename_ext=None, # not used
                 redirect_output=False, # not used
                 log=StringIO()):
    old_preamble = self.preamble
    start_energy, junk = self.get_energy(optimise_h=True, cleanup=cleanup, log=log)
    self.preamble = old_preamble+'_strain'
    final_energy, junk = self.get_opt(redirect_output=False, cleanup=cleanup, log=log)
    final_energy, junk = self.read_energy()
    self.strain = start_energy-final_energy
    print('  Strain energy = %0.5f %s' % (self.strain, self.units), file=log)
    self.preamble = old_preamble
    return self.strain, self.units

  def cleanup(self, level=None, verbose=False):
    if level=='all':
      assert 0

