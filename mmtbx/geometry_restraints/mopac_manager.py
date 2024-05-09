from __future__ import absolute_import, division, print_function
import os
import time

from scitbx.array_family import flex
from mmtbx.geometry_restraints import base_qm_manager

import libtbx.load_env

from mmtbx.geometry_restraints.base_qm_manager import get_internal_coordinate_value

def get_exe():
  bin_dir = libtbx.env.under_base('bin')
  exe_path = os.path.join(bin_dir, 'mopac')
  bin_dir = libtbx.env.under_base('Library')
  win_exe_path = os.path.join(bin_dir,"bin", 'mopac.exe')
  if os.environ.get('PHENIX_MOPAC', False):
    return os.environ['PHENIX_MOPAC']
  elif os.path.exists(exe_path):
    return exe_path
  elif os.path.exists(win_exe_path):
    return win_exe_path
  return False

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
    outl = '%s %s %s %s DISP \n%s\n\n' % (
     self.method,
     self.basis_set,
     self.solvent_model,
     'CHARGE=%s %s' % (self.charge, nproc_str),
     self.preamble,
     )
    return outl

  def get_coordinate_lines(self,
                           optimise_ligand=True,
                           optimise_h=True,
                           constrain_torsions=False,
                           verbose=False):
    outl = ''
    for i, atom in enumerate(self.atoms):
      opt = self._is_atom_for_opt(i,
                                  atom,
                                  optimise_ligand=optimise_ligand,
                                  optimise_h=optimise_h)
      tmp = ''
      if constrain_torsions:
        def _get_B_A_D(torsions, round_d=False):
          b=get_internal_coordinate_value(self.atoms[torsions[0]],
                                          self.atoms[torsions[1]],
                                          )
          a=get_internal_coordinate_value(self.atoms[torsions[0]],
                                          self.atoms[torsions[1]],
                                          self.atoms[torsions[2]],
                                          )
          d=get_internal_coordinate_value(self.atoms[torsions[0]],
                                          self.atoms[torsions[1]],
                                          self.atoms[torsions[2]],
                                          self.atoms[torsions[3]],
                                          )
          if round_d:
            if abs(d)<45:
              d=0.
            elif abs(180-abs(d))<45:
              d=180.
            else:
              assert 0
          return b,a,d
        bad=None
        if opt and not atom.element_is_hydrogen():
          torsions = self.get_torsion(i)
          # C         1.5 0 120 1 180 1 2 5 3
          bad=list(_get_B_A_D(torsions))
          bad.insert(0, atom.element)
        if opt and atom.element_is_hydrogen():
          ligand_atom = self.ligand_atoms_array[i]
          if ligand_atom and atom.name in [' HD1', ' HE2']:
            torsions = self.get_torsion(i)
            bad = list(_get_B_A_D(torsions, round_d=True))
            bad.insert(0, atom.element)
        if bad:
          tmp = ' %s %12.5f 1 %12.5f 1 %12.5f 0' % tuple(bad)
          for j in torsions[1:]: tmp += ' %s' % (j+1)
          tmp += '\n'
      if tmp:
        outl += tmp
      else:
        outl += ' %s %0.5f %d %0.5f %d %0.5f %d\n' % (
          atom.element,
          atom.xyz[0],
          opt,
          atom.xyz[1],
          opt,
          atom.xyz[2],
          opt,
          )
      if verbose and ligand_atom: print(atom.quote())
    outl += '\n'
    return outl

  def get_input_lines(self, optimise_ligand=True, optimise_h=True, constrain_torsions=False):
    outl = self._input_header()
    outl += self.get_coordinate_lines(optimise_ligand=optimise_ligand,
                                      optimise_h=optimise_h,
                                      constrain_torsions=constrain_torsions,
                                      )
    return outl

  def get_input_filename(self):
    return 'mopac_%s.mop' % self.preamble

  def read_xyz_output(self):
    filename = self.get_coordinate_filename()
    filename = filename.replace('.arc', '.out')
    if not os.path.exists(filename):
      raise Sorry('QM output filename not found: %s' % filename)
    f=open(filename, 'r')
    lines = f.read()
    del f
    rc = flex.vec3_double()
    read_xyz = False
    i_done = 1e9
    for i, line in enumerate(lines.splitlines()):
      if i>i_done and not line.strip(): read_xyz=False
      if read_xyz:
        if line.find('NO.       ATOM           X           Y           Z')>-1: continue
        tmp = line.split()
        if len(tmp) in [5]:
          rc.append((float(tmp[2]), float(tmp[3]), float(tmp[4])))
      if line.find('CARTESIAN COORDINATES')>-1:
        read_xyz=True
        i_done=i+2
    return rc

  def get_cmd(self):
    cmd = '%s mopac_%s' % (
      get_exe(),
      self.preamble,
      )
    return cmd

  def run_cmd(self, redirect_output=False, log=None):
    t0=time.time()
    cmd = self.get_cmd()
    base_qm_manager.run_qm_cmd(cmd,
                               'mopac_%s.out' % self.preamble,
                               error_lines=self.error_lines,
                               redirect_output=redirect_output,
                               log=log,
                               )
    self.times.append(time.time()-t0)

  def read_energy(self, filename=None):
    lines = self.get_lines(filename=filename)
    # FINAL HEAT OF FORMATION =       -132.17152 KCAL/MOL =    -553.00562 KJ/MOL
    for line in lines.splitlines():
      if line.find('FINAL HEAT OF FORMATION = ')>-1:
        self.energy = float(line.split()[5])
        self.units = line.split()[6].lower()
    return self.energy, self.units

  def read_charge(self, filename=None):
    lines = self.get_lines(filename=filename)
    #  CHARGE ON SYSTEM =
    for line in lines.splitlines():
      if line.find('CHARGE ON SYSTEM = ')>-1:
        self.charge = float(line.split()[5])
        break
    return self.charge

  def read_gradients(self, filename=None):
    lines = self.get_lines(filename=filename)
    print(lines)
    assert 0

  def cleanup(self, level=None, verbose=False):
    if level=='all':
      assert 0

if __name__ == '__main__':
  exe = get_exe()
  print('mopac executable',exe)
