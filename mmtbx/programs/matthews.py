"""Calculate the expected Matthews coefficient given the crystal symmetry and
crystallized molecule(s)"""
from __future__ import absolute_import, division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
from iotbx import crystal_symmetry_from_any
import iotbx.bioinformatics
from cctbx import crystal
from libtbx.utils import Sorry

master_phil_str = """
space_group = None
  .type = space_group
unit_cell = None
  .type = unit_cell
n_residues = None
  .type = int(value_min=1)
  .optional = True
n_bases = None
  .type = int(value_min=1)
  .optional = True
"""

# =============================================================================

class Program(ProgramTemplate):

  description = '''
Calculate the expected Matthews coefficient given the crystal symmetry and
crystallized molecule(s).

phenix.matthews [data.hkl] [space_group] [unit_cel] [sequence] [n_residues] ...
'''

  datatypes = ['model', 'phil', 'miller_array', 'sequence']

  master_phil_str = master_phil_str

  # ---------------------------------------------------------------------------

  def validate(self):
    '''
    Make sure there is correct amount of inputs
    '''
    print('Validating inputs...\n', file=self.logger)
    #
    # number of files
    if len(self.data_manager.get_miller_array_names()) > 1:
      raise Sorry('Supply at most one reflection file.')
    if len(self.data_manager.get_model_names()) > 1:
      raise Sorry('Supply at most one model file.')
    if len(self.data_manager.get_sequence_names()) > 1:
      raise Sorry('Supply at most one model file.')
    # Space group & unit cell - either from data file or from command line
    if (self.params.space_group is None) and (self.params.unit_cell is None):
      if not self.data_manager.has_miller_arrays():
        raise Sorry('You must supply both a space group and a unit cell, ' +
          'either via the command line or a data file).')
    elif (self.params.space_group is None) or (self.params.unit_cell is None):
        raise Sorry('Supply the space group and a unit cell either via the '+
          'command line or a data file.')
    else:
      if self.data_manager.has_miller_arrays():
        raise Sorry('No need to use a data file if space group and unit cell '+
          'are supplied.')
    # composition: either from sequence, model or command line
    if self.data_manager.has_models() and self.data_manager.has_sequences():
      raise Sorry('Supply either a model file or a sequence.')
    if not self.data_manager.has_sequences() and not self.data_manager.has_models():
      if (self.params.n_residues is None) and (self.params.n_bases is None):
        raise Sorry('You must specify the composition of the crystallized '+
          'entity - either a sequence, a partial model, or the number of '+
          'protein residues or nucleic acid bases.')
    else :
      if (self.params.n_residues is not None) or (self.params.n_bases is not None):
        raise Sorry('You may only specify a sequence file OR a partial model '+
          'OR the numbers of residues and/or bases.')

  # ---------------------------------------------------------------------------

  def run(self):
    '''
    Calculate Matthews coefficient and print the results
    '''
    # space group
    print('Space group and unit cell:', file=self.logger)
    if (self.params.space_group is None) and (self.params.unit_cell is None):
      data_fn = self.data_manager.get_miller_array_names()[0]
      print('From data file ', data_fn, file=self.logger)
      cs = crystal_symmetry_from_any.extract_from(file_name=data_fn)
      if cs is not None:
        sg_from_file = cs.space_group_info()
        if sg_from_file is not None:
          self.params.space_group = cs.space_group_info()
        uc_from_file = cs.unit_cell()
        if uc_from_file is not None:
          self.params.unit_cell = uc_from_file
    else:
      print('From command line inputs', file=self.logger)
    print('Space group:', self.params.space_group, file=self.logger)
    print('Unit cell:', self.params.unit_cell, file=self.logger)

    print('\nComposition:', file=self.logger)

    # composition from sequence
    if self.data_manager.has_sequences():
      print('Using sequence file ',
        self.data_manager.get_default_sequence_name(), file=self.logger)
      assert (self.params.n_residues == self.params.n_bases == None)
      seq_object = self.data_manager.get_sequence()[0]
      n_residues, n_bases = iotbx.bioinformatics.composition_from_sequence(
        sequence=seq_object.sequence)
      self.params.n_residues = n_residues
      self.params.n_bases = n_bases
    #else :
    #  raise Sorry("No composition information could be obtained from the "+
    #    "sequence file.")

    # composition from model file
    if self.data_manager.has_models():
      print('Using model file ',
        self.data_manager.get_default_model_name(), file=self.logger)
      assert (self.params.n_residues == self.params.n_bases == None)
      model = self.data_manager.get_model()
      self.params.n_residues = 0
      self.params.n_bases = 0
      hierarchy = model.get_hierarchy()
      comp = hierarchy.composition()
      self.params.n_residues = comp.n_protein
      self.params.n_bases = comp.n_nucleotide
    if (self.params.n_residues):
        print('Number of residues: %d' % self.params.n_residues, file=self.logger)
    if (self.params.n_bases):
        print('Number of bases: %d' % self.params.n_bases, file=self.logger)

    symm = crystal.symmetry(
      space_group_info=self.params.space_group,
      unit_cell=self.params.unit_cell)
    from mmtbx.scaling import matthews
    result = matthews.matthews_rupp(
      crystal_symmetry=symm,
      n_residues=self.params.n_residues,
      n_bases=self.params.n_bases)
    result.show(out=self.logger)
    return result
