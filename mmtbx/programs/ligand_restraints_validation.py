from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate

# from cctbx.maptbx.box import shift_and_box_model

class Program(ProgramTemplate):

  description = '''
mmtbx.ligand_restraints_validation:

Usage examples:
  mmtbx.ligand_restraints_validation ligand.cif
  '''

  datatypes = ['phil', 'restraint']

  master_phil_str = """
  ligand_restraints_validation {
    input {
      restraints = None
        .type = path
      }
  }
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    pass

  # ---------------------------------------------------------------------------
  def get_energies_sites(self, model, use_hydrogens=False):
    assert model.restraints_manager is not None
    if(use_hydrogens):
      rm = model.restraints_manager
      sc = model.get_sites_cart()
    else:
      hd_sel = model.get_xray_structure().hd_selection()
      rm = model.restraints_manager.select(~hd_sel)
      sc = model.get_sites_cart().select(~hd_sel)
    energies_sites = \
      rm.geometry.energies_sites(
        sites_cart        = sc,
        compute_gradients = False)
    return energies_sites #.bond_deviations()[2]

  def processed_ligand_restaints(self, filename):
    from mmtbx.monomer_library.geostd_utils import get_as_hierarchy
    from mmtbx.monomer_library.geostd_utils import as_cif_object
    hierarchy=get_as_hierarchy(filename)
    hierarchy.sort_atoms_in_place()
    # hierarchy.show()
    starting_atoms=hierarchy.atoms()
    starting_xyz=starting_atoms.extract_xyz()
    if 0: hierarchy.write_pdb_file('test.pdb')
    #
    from mmtbx.model.model import manager
    m = manager(pdb_hierarchy=hierarchy,
                restraint_objects=[(filename, as_cif_object(filename))],
                )
    # m = shift_and_box_model(model = m, box_cushion = 5.)
    # for atom in m.get_hierarchy().atoms(): print(atom.format_atom_record())
    #
    # geom min
    #
    from mmtbx.refinement import geometry_minimization
    m.process(make_restraints=True)
    geometry_minimization.run2(
      restraints_manager = m.get_restraints_manager(),
      pdb_hierarchy = m.get_hierarchy(),
      correct_special_position_tolerance = 1.,
      bond        = True,
      nonbonded   = True,
      angle       = True,
      dihedral    = True,
      chirality   = True,
      planarity   = True,
      parallelity = True,
      )
    # m.set_sites_cart_from_hierarchy()
    # for atom in m.get_hierarchy().atoms(): print(atom.format_atom_record())
    h=m.get_hierarchy()
    if 0: h.write_pdb_file('test.pdb')
    ending_xyz=h.atoms().extract_xyz()
    self.results={}
    rc = starting_xyz.rms_difference(ending_xyz)
    print('\n  RMSD of starting and ending coordinates : %5.2f' % rc)
    self.results['RMSD']=rc
    es = self.get_energies_sites(m, use_hydrogens=True)
    # es.show()
    self.results['energies_sites']=es
    rc = m.restraints_as_geo()
    print(rc[:1000])
    # print(rc)

  def run(self, log=None):
    for filename in self.data_manager.get_restraint_names():
      self.processed_ligand_restaints(filename)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results
