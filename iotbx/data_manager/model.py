from __future__ import absolute_import, division, print_function
'''
'''

import iotbx.pdb
import mmtbx.model

from iotbx.file_reader import any_file
from iotbx.data_manager import DataManagerBase
from libtbx import Auto
from libtbx.utils import Sorry

# =============================================================================
class ModelDataManager(DataManagerBase):

  datatype = 'model'

  # ---------------------------------------------------------------------------
  # Models
  def add_model_phil_str(self):
    '''
    Add custom PHIL and storage for type
    '''

    # set up storage
    # self._model_types = dict()  # [filename] = type
    self._model_types = dict()
    self._default_model_type = 'x_ray'
    self._possible_model_types = ['x_ray', 'neutron', 'electron']

    # custom PHIL section
    custom_phil_str = '''
model
  .multiple = True
{
  file = None
    .type = path
    .short_caption = Model file
    .style = file_type:pdb input_file
  type = *%s
    .type = choice(multi=False)
}
''' % ' '.join(self._possible_model_types)

    # custom PHIL scope
    self._custom_model_phil = iotbx.phil.parse(custom_phil_str)

    return custom_phil_str

  def export_model_phil_extract(self):
    '''
    Export custom PHIL extract
    '''
    extract = list()
    filenames = self.get_model_names()
    for filename in filenames:
      item_extract = self._custom_model_phil.extract().model[0]
      item_extract.file = filename
      item_extract.type = self._model_types.get(
        filename, self._default_model_type)
      extract.append(item_extract)
    return extract

  def load_model_phil_extract(self, phil_extract):
    '''
    Load custom PHIL extract
    '''
    extract = phil_extract.data_manager.model
    for item_extract in extract:
      if ((not hasattr(item_extract, 'file')) or
          (not hasattr(item_extract, 'type'))):
        raise Sorry('This PHIL is not properly defined for the "model" datatype.\n There should be a parameter for the filename ("file") and type ("type").\n')

      # process file
      self.process_model_file(item_extract.file)
      self._model_types[item_extract.file] = item_extract.type

  def add_model(self, filename, data):
    return self._add(ModelDataManager.datatype, filename, data)

  def set_default_model_type(self, model_type):
    if (model_type not in self._possible_model_types):
      raise Sorry('Unrecognized model type, "%s," possible choices are %s.' %
                  (model_type, ', '.join(self._possible_model_types)))
    self._default_model_type = model_type

  def get_default_model_type(self):
    return self._default_model_type

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

  def set_model_type(self, filename=None, model_type=None):
    if (filename is None):
      filename = self.get_default_model_name()
    if (model_type is None):
      model_type = self._default_model_type
    elif (model_type not in self._possible_model_types):
      raise Sorry('Unrecognized model type, "%s," possible choices are %s.' %
                  (model_type, ', '.join(self._possible_model_types)))
    self._model_types[filename] = model_type

  def get_model_type(self, filename=None):
    if (filename is None):
      filename = self.get_default_model_name()
    return self._model_types.get(filename, self._default_model_type)

  def get_model_names(self, model_type=None):
    all_names = self._get_names(ModelDataManager.datatype)
    names = list()
    if (model_type is None):
      names = all_names
    else:
      for filename in all_names:
        if (model_type == self.get_model_type(filename)):
          names.append(filename)
    return names

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

  def get_default_output_model_filename(self):
    '''
    Function for returning the filename with extension. By default ".cif" will
    be used.
    '''
    filename = self.get_default_output_filename()
    if not (filename.endswith('.cif') or filename.endswith('.pdb')):
      filename += '.cif'
    return filename

  def write_model_file(self, model_str, filename=Auto, overwrite=Auto):
    if filename is Auto:
      filename = self.get_default_output_model_filename()
    self._write_text(ModelDataManager.datatype, model_str,
                     filename=filename, overwrite=overwrite)

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
