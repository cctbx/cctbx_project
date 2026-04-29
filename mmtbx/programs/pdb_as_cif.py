"""Convert PDB formatted model int mmCIF"""
from __future__ import absolute_import, division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import os

master_phil_str = '''
model_file_name = None
  .type = path
  .short_caption = Model file
  .multiple = False
  .help = Model file name
  .style = file_type:pdb bold input_file
'''

class Program(ProgramTemplate):
  datatypes = ['model', 'phil', 'restraint']
  master_phil_str = master_phil_str

  def validate(self):
    model = self.data_manager.get_model()
    hierarchy = model.get_hierarchy()
    pdb_atoms = hierarchy.atoms()
    pdb_atoms.set_chemical_element_simple_if_necessary()
    elements = pdb_atoms.extract_element().strip()
    if (not elements.all_ne("")):
      n_missing = elements.count("")
      raise RuntimeError("Missing element symbol for %d atoms." % n_missing)

  def run(self):
    file_name = self.data_manager.get_model_names()[0]
    print("Converting %s to mmCIF format." %file_name)
    basename = os.path.splitext(os.path.basename(file_name))[0]
    model = self.data_manager.get_model()
    txt = model.model_as_mmcif()
    print(txt)
    with open(basename+".cif", 'w') as f:
      f.write(txt)
    print("  wrote %s.cif" % basename)

