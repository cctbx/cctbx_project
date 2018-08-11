from __future__ import division, print_function
'''
'''

import iotbx.pdb
import mmtbx.model

from iotbx.file_reader import any_file
from iotbx.data_manager import DataManagerBase
from libtbx.utils import Sorry

# =============================================================================
# optional PHIL scope for specifying model usage
model_phil_str = '''
model
  .multiple = True
{
  file_name = None
    .type = path
  type = *x_ray neutron
    .type = choice(multi=False)
}
'''

# =============================================================================
class ModelDataManager(DataManagerBase):

  datatype = 'model'

  # ---------------------------------------------------------------------------
  # Models
  def add_model(self, filename, data):
    return self._add(ModelDataManager.datatype, filename, data)

  def set_default_model(self, filename):
    return self._set_default(ModelDataManager.datatype, filename)

  def get_model(self, filename=None):
    model = self._get(ModelDataManager.datatype, filename)
    if (self.supports('restraint')):
      restraint_objects = list()
      for filename in self.get_restraint_names():
        restraint_objects.append((filename, self.get_restraint(filename)))
      model.set_restraint_objects(restraint_objects)
    return model

  def get_model_names(self):
    return self._get_names(ModelDataManager.datatype)

  def get_default_model_name(self):
    return self._get_default_name(ModelDataManager.datatype)

  def remove_model(self, filename):
    return self._remove(ModelDataManager.datatype, filename)

  def has_models(self, expected_n=1, exact_count=False, raise_sorry=False):
    return self._has_data(ModelDataManager.datatype, expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_model_file(self, filename, pdb_interpretation_extract=None):
    # unique because any_file does not return a model object
    if (filename not in self.get_model_names()):
      a = any_file(filename)
      if (a.file_type != 'pdb'):
        raise Sorry('%s is not a recognized model file' % filename)
      else:
        model_in = iotbx.pdb.input(a.file_name)
        model = mmtbx.model.manager(
          model_input=model_in,
          pdb_interpretation_params=pdb_interpretation_extract,
          log=self.logger)
        self.add_model(filename, model)

  def process_model_str(self, label, model_str, pdb_interpretation_extract=None):
    model = mmtbx.model.manager(
      model_input=iotbx.pdb.input(source_info=None, lines=model_str),
      pdb_interpretation_params=pdb_interpretation_extract,
      log=self.logger)
    self.add_model(label, model)

  def write_model_file(self, filename, model_str, overwrite=False):
    self._write_text(ModelDataManager.datatype, filename,
                     model_str, overwrite=overwrite)

  def update_pdb_interpretation_for_model(
    self, filename, pdb_interpretation_extract):
    '''
    Pass PHIL extract to model class
    model class handles finding and matching of PHIL scopes
    '''
    self.get_model(filename=filename).set_pdb_interpretation_params(
      pdb_interpretation_extract)

# =============================================================================
# end
