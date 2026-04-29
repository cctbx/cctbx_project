"""Validate ligand restraints"""
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
    action {
      use_hydrogens = True
        .type = bool
    }
    output {
      write_pdb = None
        .type = path
      write_geo = None
        .type = path
      write_input_pdb = None
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
    params = self.params.ligand_restraints_validation
    if params.output.write_input_pdb:
      hierarchy.write_pdb_file(params.output.write_input_pdb)
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
    from six.moves import StringIO
    sio=StringIO()
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
      log         = sio,
      )
    # m.set_sites_cart_from_hierarchy()
    # for atom in m.get_hierarchy().atoms(): print(atom.format_atom_record())
    h=m.get_hierarchy()
    if 0: h.write_pdb_file('test.pdb')
    if params.output.write_pdb:
      h.write_pdb_file(params.output.write_pdb)
    ending_xyz=h.atoms().extract_xyz()

    result={}
    rc = starting_xyz.rms_difference(ending_xyz)
    print('\n  RMSD of starting and ending coordinates : %5.2f' % rc, file=self.logger)
    result['RMSD']=rc
    result['atoms']=len(ending_xyz)
    es = self.get_energies_sites(m, use_hydrogens=params.action.use_hydrogens)
    # es.show()
    result['energies_sites']=es
    print('\nGeometry Min. result')
    es.show(f=self.logger)
    if params.output.write_geo:
      gs = m.restraints_as_geo()
      f=open(params.output.write_geo, 'w')
      f.write(gs)
      del f
    return result

  def run(self, log=None):
    self.results={}
    for filename in self.data_manager.get_restraint_names():
      rc = self.processed_ligand_restaints(filename)
      self.results[filename]=rc

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results
