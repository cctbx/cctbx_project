from __future__ import absolute_import, division, print_function
from StringIO import StringIO
# import sys, time

# from libtbx.utils import Sorry
# from scitbx.array_family import flex
from libtbx import adopt_init_args
# from libtbx.str_utils import make_header

from cctbx.geometry_restraints.manager import manager as standard_manager

orca_master_phil_str = '''
  use_orca = True
    .type = bool
  method = PM6
    .type = str
  basis_set = None
    .type = str
  solvent_model = CPCM
    .type = str
'''

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

  def set_qm_atoms(self, qm_atoms):
    self.qm_atoms = qm_atoms

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
      compute_gradients=False, #compute_gradients,
      gradients=gradients,
      disable_asu_cache=disable_asu_cache,
      normalization=normalization,
      external_energy_function=external_energy_function,
      extension_objects=extension_objects,
      site_labels=site_labels,
      )
    print(result)

    assert 0

def digester(model,
             standard_geometry_restraints_manager,
             params,
             log=StringIO(),
             ):
  sgrm = standard_geometry_restraints_manager
  qi_grm = manager(params, log=log)
  for attr, value in vars(sgrm).items(): setattr(qi_grm, attr, value)
  qi_grm.standard_geometry_restraints_manager = sgrm
  #
  qm_atoms = []
  ligand_selection = model.selection(params.qi.selection)
  buffer_selection_string = 'within(%s, %s)' % (params.qi.buffer,
                                                params.qi.selection)
  buffer_selection = model.selection(buffer_selection_string)
  pdb_hierarchy = model.get_hierarchy()
  atoms = pdb_hierarchy.atoms()
  done = []
  for i, b in enumerate(ligand_selection.iselection()):
    qm_atoms.append(atoms[b])
    done.append(b)
  for i, b in enumerate(buffer_selection.iselection()):
    if b in done: continue
    qm_atoms.append(atoms[b])
    done.append(b)
    ag = qm_atoms[-1].parent()
    assert len(ag.parent().atom_groups())==1
    for atom in ag.atoms():
      qm_atoms.append(atom)
      done.append(atom.i_seq)
  qi_grm.set_qm_atoms(qm_atoms)
  return qi_grm
