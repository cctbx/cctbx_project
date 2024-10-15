from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate

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
  def processed_ligand_restaints(self, filename):
    from mmtbx.monomer_library.geostd_utils import get_as_hierarchy
    hierarchy=get_as_hierarchy(filename)
    hierarchy.show()
    hierarchy.write_pdb_file('test.pdb')
    print(hierarchy)
    from mmtbx.model import model
    m = model(pdb_hierarchy=hierarchy)
    print(m)

  def run(self, log=None):
    for filename in self.data_manager.get_restraint_names():
      self.processed_ligand_restaints(filename)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results
