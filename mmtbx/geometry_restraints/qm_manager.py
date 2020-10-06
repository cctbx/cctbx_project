from __future__ import absolute_import, division, print_function
import os
from StringIO import StringIO
import time
import tempfile

# from libtbx.utils import Sorry
from scitbx.array_family import flex
from libtbx import adopt_init_args
from libtbx import easy_run
# from libtbx.str_utils import make_header

from cctbx.geometry_restraints.manager import manager as standard_manager

orca_master_phil_str = '''
  use_orca = True
    .type = bool
  method = PM3
    .type = str
  basis_set = None
    .type = str
  solvent_model = None
    .type = str
'''

harkcal = 627.50946900
bohrang = 0.52918

class orca_manager:
  def __init__(self,
               atoms,
               method,
               basis_set,
               solvent_model,
               charge,
               multiplicity,
               preamble=None,
               log=None,
               ):
    adopt_init_args(self, locals())
    self.times = flex.double()
    self.energies = {}
    if self.preamble is None:
      self.preamble = os.path.basename(tempfile.NamedTemporaryFile().name)

  def __repr__(self):
    outl = 'QI Manager\n'
    outl += ' charge: %s multiplicity: %s\n method: %s basis: "%s" solvent: "%s"\n' % (
      self.charge,
      self.multiplicity,
      self.method,
      self.basis_set,
      self.solvent_model,
      )
    for atom in self.atoms:
      outl += '  %s\n' % atom.quote()
    return outl

  def set_sites_cart(self, sites_cart):
    assert len(self.atoms)==len(sites_cart)
    for atom, site_cart in zip(self.atoms, sites_cart):
      atom.xyz = site_cart

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
    f=file('orca_%s.engrad' % self.preamble, 'rb')
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

  def read_xyz_output(self):
    f=file('orca_%s.xyz' % self.preamble, 'rb')
    lines = f.read()
    del f
    rc = flex.vec3_double()
    for i, line in enumerate(lines.splitlines()):
      if i>=2:
        tmp = line.split()
        rc.append((float(tmp[1]), float(tmp[2]), float(tmp[3])))
    return rc

  def write_input(self, outl):
    f=file('orca_%s.in' % self.preamble, 'wb')
    f.write(outl)
    del f

  def run_cmd(self):
    t0=time.time()
    cmd = '%s orca_%s.in' % (
      os.environ['PHENIX_ORCA'],
      self.preamble,
      )
    # print(cmd)
    rc = easy_run.go(cmd)
    if rc.stderr_lines:
      print('stderr')
      for line in rc.stderr_lines:
        print(line)
      print('stdout')
      for line in rc.stdout_lines:
        print(line)
      assert 0
    self.times.append(time.time()-t0)

  def get_coordinate_lines(self):
    outl = '* xyz %s %s\n' % (self.charge, self.multiplicity)
    for atom in self.atoms:
      outl += ' %s %0.5f %0.5f %0.5f\n' % (
        atom.element,
        atom.xyz[0],
        atom.xyz[1],
        atom.xyz[2],
        )
    outl += '*\n'
    return outl

  def print_timings(self, energy=None):
    print('  Timings : %0.2fs (%ss) Energy : %0.6f' % (
      self.times[-1],
      self.times.format_mean(format='%.2f'),
      energy))

  def get_engrad(self):
    outl = '! %s %s %s EnGrad\n\n' % (self.method,
                                      self.basis_set,
                                      self.solvent_model)
    outl += self.get_coordinate_lines()
    if outl in self.energies:
      self.times.append(0)
      return self.energies[outl]
    self.write_input(outl)
    self.run_cmd()
    energy, gradients = self.read_engrad_output()
    self.print_timings(energy)
    self.energies[outl] = (energy, gradients)
    return energy, gradients

  def get_opt(self):
    outl = '! %s %s %s Opt\n\n' % (self.method,
                                   self.basis_set,
                                   self.solvent_model)
    outl += self.get_coordinate_lines()
    self.write_input(outl)
    self.run_cmd()
    coordinates = self.read_xyz_output()
    return coordinates

class manager(standard_manager):
  def __init__(self,
               # pdb_hierarchy,
               params,
               energy_components=None,
               # gradients=None,
               # gradients_factory=None, #flex.vec3_double,
               log=StringIO()):
    # self.gradients_factory = gradients_factory
    adopt_init_args(self, locals(), exclude=["log"])
    self.validate()

  def validate(self):
    print(dir(self.params.qi))
    qi = self.params.qi
    assert qi.use_quantum_interface
    assert qi.selection
    if qi.orca.use_orca:
      print('Orca')

  def get_engrad(self, sites_cart):
    self.execution_manager.set_sites_cart(sites_cart)
    return self.execution_manager.get_engrad()

  def get_opt(self, sites_cart):
    self.execution_manager.set_sites_cart(sites_cart)
    return self.execution_manager.get_opt()

  def set_qm_info(self,
                  method,
                  basis_set,
                  solvent_model,
                  charge,
                  multiplicity,
                  ):
    adopt_init_args(self, locals())
    if self.basis_set is None:
      self.basis_set = ''
    if self.solvent_model is None:
      self.solvent_model = ''
    self.execution_manager = orca_manager( self.qm_atoms,
                                           self.method,
                                           self.basis_set,
                                           self.solvent_model,
                                           self.charge,
                                           self.multiplicity
                                           )

  def set_qm_atoms(self, qm_atoms):
    self.qm_atoms = qm_atoms
    self.qm_iseqs = []
    for atom in self.qm_atoms:
      self.qm_iseqs.append(atom.i_seq)

  def energies_sites(self,
                     sites_cart,
                     flags=None,
                     custom_nonbonded_function=None,
                     compute_gradients=False,
                     gradients=None,
                     disable_asu_cache=False,
                     normalization=False,
                     external_energy_function=None,
                     extension_objects=[],
                     site_labels=None,
                     log=None):
    result = standard_manager.energies_sites(
      self,
      sites_cart,
      flags=flags,
      custom_nonbonded_function=custom_nonbonded_function,
      compute_gradients=compute_gradients,
      gradients=gradients,
      disable_asu_cache=disable_asu_cache,
      normalization=normalization,
      external_energy_function=external_energy_function,
      extension_objects=extension_objects,
      site_labels=site_labels,
      )
    if compute_gradients:
      qm_sites_cart = []
      for i_seq in self.qm_iseqs:
        qm_sites_cart.append(sites_cart[i_seq])
      # coordinates = self.get_opt(qm_sites_cart)
      # print(list(coordinates))
      # assert 0
      energy, gradients = self.get_engrad(qm_sites_cart)
      for i_seq, gradient in zip(self.qm_iseqs, gradients):
        result.gradients[i_seq]=gradient
    return result

def digester(model,
             standard_geometry_restraints_manager,
             params,
             log=StringIO(),
             ):
  """Digest a standard GRM and create a QI GRM

  Args:
      model (TYPE): Model class
      standard_geometry_restraints_manager (TYPE): Standard GRM?
      params (TYPE): User parameteres
      log (TYPE, optional): Output log

  Returns:
      TYPE: Quantum Interface Restraints Manager

  """
  #
  # Digest
  #
  sgrm = standard_geometry_restraints_manager
  qi_grm = manager(params, log=log)
  for attr, value in vars(sgrm).items(): setattr(qi_grm, attr, value)
  qi_grm.standard_geometry_restraints_manager = sgrm
  #
  # Create selection lists
  #
  qm_atoms = []
  ligand_selection = model.selection(params.qi.selection)
  i_buffer_selection = []
  if params.qi.buffer>0.:
    buffer_selection_string = 'within(%s, %s)' % (params.qi.buffer,
                                                  params.qi.selection)
    buffer_selection = model.selection(buffer_selection_string)
    i_buffer_selection = buffer_selection.iselection()
  pdb_hierarchy = model.get_hierarchy()
  atoms = pdb_hierarchy.atoms()
  done = []
  for i, b in enumerate(ligand_selection.iselection()):
    qm_atoms.append(atoms[b])
    done.append(b)
  for i, b in enumerate(i_buffer_selection):
    if b in done: continue
    qm_atoms.append(atoms[b])
    done.append(b)
    ag = qm_atoms[-1].parent()
    assert len(ag.parent().atom_groups())==1
    for atom in ag.atoms():
      qm_atoms.append(atom)
      done.append(atom.i_seq)
  #
  # Add to QI GRM
  #
  qi_grm.set_qm_atoms(qm_atoms)
  qi_grm.set_qm_info(params.qi.orca.method,
                     params.qi.orca.basis_set,
                     params.qi.orca.solvent_model,
                     params.qi.charge,
                     params.qi.multiplicity,
                     )
  return qi_grm

def main():
  from iotbx import pdb
  pdb_lines = '''
HETATM   97  S   SO4 A  13      31.477  38.950  15.821  0.50 25.00           S
HETATM   98  O1  SO4 A  13      31.243  38.502  17.238  0.50 25.00           O
HETATM   99  O2  SO4 A  13      30.616  40.133  15.527  0.50 25.00           O
HETATM  100  O3  SO4 A  13      31.158  37.816  14.905  0.50 25.00           O
HETATM  101  O4  SO4 A  13      32.916  39.343  15.640  0.50 25.00           O
'''
  pdb_inp = pdb.input(lines=pdb_lines, source_info='lines')
  qi_grm = orca_manager(pdb_inp.atoms(),
                        'PM3',
                        '',
                        '',
                        -2,
                        1,
                        preamble='test',
                        )
  print(qi_grm)
  energy, gradients = qi_grm.get_engrad()
  print(energy, list(gradients))
  coordinates = qi_grm.get_opt()
  print(list(coordinates))

if __name__ == '__main__':
  main()
