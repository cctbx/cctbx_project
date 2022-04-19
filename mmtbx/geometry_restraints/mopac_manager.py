from __future__ import absolute_import, division, print_function
import os
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

  def get_coordinate_lines(self, optimise_ligand=True, optimise_h=True, verbose=False):
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
      if verbose and ligand_atom: print(atom.quote())
    outl += '\n'
    return outl

  def get_input_lines(self, optimise_ligand=True, optimise_h=True):
    outl = self._input_header()
    outl += self.get_coordinate_lines(optimise_ligand=optimise_ligand,
                                      optimise_h=optimise_h)
    return outl

  def get_input_filename(self):
    return 'mopac_%s.mop' % self.preamble

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

  def cleanup(self, level=None, verbose=False):
    if level=='all':
      assert 0

