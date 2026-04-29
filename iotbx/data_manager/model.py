from __future__ import absolute_import, division, print_function
'''
'''

import iotbx.pdb
import mmtbx.model

from iotbx.file_reader import any_file
from iotbx.data_manager import DataManagerBase
from libtbx import Auto
from libtbx.utils import Sorry
import os

# =============================================================================
class ModelDataManager(DataManagerBase):

  datatype = 'model'

  # ---------------------------------------------------------------------------
  # Models
  def add_model_phil_str(self):
    '''
    Add custom PHIL and storage for type
    The type is a choice(multi=True) PHIL parameter, so it is always a list
    '''

    # set up storage
    # self._model_types = dict()  # [filename] = [type]
    self._model_types = dict()
    self._default_model_type = ['x_ray']
    self._possible_model_types = ['x_ray', 'neutron', 'electron', 'reference']

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
    .type = choice(multi=True)
    .short_caption = Model type(s)
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

  def load_model_phil_extract(self, phil_extract, process_files=True):
    '''
    Load custom PHIL extract
    '''
    extract = phil_extract.data_manager.model
    for item_extract in extract:
      if ((not hasattr(item_extract, 'file')) or
          (not hasattr(item_extract, 'type'))):
        raise Sorry('This PHIL is not properly defined for the "model" datatype.\n There should be a parameter for the filename ("file") and type ("type").\n')

      # process file
      if process_files and item_extract.file is not None:
        self.process_model_file(item_extract.file)
        self._model_types[item_extract.file] = item_extract.type

  def add_model(self, filename, data):
    return self._add(ModelDataManager.datatype, filename, data)

  def _is_valid_model_type(self, model_type):
    """
    Convenience function for checking if the model type is valid
    This will also check that model_type is a list to conform with the
    PHIL parameter

    Parameters
    ----------
    model_type: list
      The model_type(s) to check.

    Returns
    -------
    bool:
    """
    if not isinstance(model_type, list):
      raise Sorry('The model_type argument must be a list.')
    if len(model_type) == 0:
      return False
    valid = True
    for mt in model_type:
      valid = valid and (mt in self._possible_model_types)
    return valid

  def set_target_output_format(self, target_output_format):
   if target_output_format == 'mmcif': target_output_format='cif'
   if not target_output_format in ['cif','pdb']:
     raise  Sorry("Target output format (%s) not recognized, options are" %(
       target_output_format) + "pdb or cif")
   self._target_output_format = target_output_format

  def set_default_model_type(self, model_type):
    if not self._is_valid_model_type(model_type):
      raise Sorry('Unrecognized model type, "%s," possible choices are %s.' %
                  (model_type, ', '.join(self._possible_model_types)))
    self._default_model_type = model_type

  def get_default_model_type(self):
    return self._default_model_type

  def set_default_model(self, filename):
    return self._set_default(ModelDataManager.datatype, filename)

  def get_model(self, filename=None, model_type=None):
    """
    Retrieve a stored mmtbx.model.manager object

    If model_type is None and there is only one model type, then the
    model is returned. If there is more than one model type, then a
    Sorry is raised.

    If a model_type is specified when a model has more than one type, a
    copy of the model is returned if model_type is not the default type.

    Parameters
    ----------
    filename : str
        Optionally specify which model using its filepath
    model_type: str
        Optionally specify the type of the model
        The options are the same as for the scattering dictionary
        ["n_gaussian", "wk1995", "it1992", "electron", "neutron"] and
        "x_ray" which will default to "n_gaussian".

    Returns
    -------
    model
        The mmtbx.model.manager object

    """
    model = self._get(ModelDataManager.datatype, filename)
    if self.supports('restraint'):
      restraint_objects = list()
      for restraint_filename in self.get_restraint_names():
        restraint_objects.append((restraint_filename, self.get_restraint(restraint_filename)))
      model.set_restraint_objects(restraint_objects)
    if hasattr(model,'info'):  # save filename if possible
      if filename is None:
        filename = self.get_default_model_name()
      if filename:
        model.info().full_file_name = os.path.abspath(filename)
        model.info().file_name = os.path.split(filename)[-1]
    if model_type is None:
      if len(self.get_model_type(filename=filename)) > 1:
        raise Sorry('''
There is more than one model type, {}. You must specify one.
'''.format(self.get_model_type(filename=filename)))
    else:
      check_type = self.map_scattering_table_type(model_type)
      if check_type not in self.get_model_type(filename=filename):
        raise Sorry('''
The model type, {}, is not one of the types set for the model, {}.
The choices are {}.
'''.format(model_type, filename, self.get_model_type(filename=filename)))
      if len(self.get_model_type(filename=filename)) > 1 \
        and model_type not in self.get_default_model_type():
          model = model.deep_copy()

      # set scattering dictionary based on model type
      # if model_type == 'x_ray':
      #   model_type = 'n_gaussian'
      # model.setup_scattering_dictionaries(scattering_table=model_type)

    return model

  def set_model_type(self, filename=None, model_type=None):
    if (filename is None):
      filename = self.get_default_model_name()
    if (model_type is None):
      model_type = self._default_model_type
    elif not self._is_valid_model_type(model_type):
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
        if (model_type in self.get_model_type(filename)):
          names.append(filename)
    return names

  def get_default_model_name(self):
    return self._get_default_name(ModelDataManager.datatype)

  def remove_model(self, filename):
    return self._remove(ModelDataManager.datatype, filename)

  def has_models(self, expected_n=1, exact_count=False, raise_sorry=False):
    return self._has_data(ModelDataManager.datatype, expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_model_file(self, filename):
    """
    Parse a model file and store the mmtbx.model.manager object

    Parameters
    ----------
    filename : str
        The filepath as a string

    Returns
    -------
    filename : str
        The model filename added to the DataManager

    """
    # unique because any_file does not return a model object
    if (filename not in self.get_model_names()):
      a = any_file(filename)
      if (a.file_type != 'pdb'):
        raise Sorry('%s is not a recognized model file' % filename)
      else:
        model_in = a.file_content.input
        expand_with_mtrix = True  # default
        skip_ss_annotations = False
        if 'model_skip_expand_with_mtrix' in self.custom_options:
          expand_with_mtrix = False
        if 'model_skip_ss_annotations' in self.custom_options:
          skip_ss_annotations = True
        model = mmtbx.model.manager(
          model_input=model_in,
          expand_with_mtrix=expand_with_mtrix,
          skip_ss_annotations=skip_ss_annotations,
          log=self.logger)
        self.add_model(filename, model)
    return filename

  def process_model_str(self, label, model_str):
    model = mmtbx.model.manager(
      model_input=iotbx.pdb.input(source_info=None, lines=model_str),
      log=self.logger)
    self.add_model(label, model)

  def get_default_output_model_filename(self, extension=Auto):
    '''
    Function for returning the filename with extension. By default ".cif" will
    be used.
    '''
    filename = self.get_default_output_filename()
    if extension is Auto:
      extension = '.cif'
    if not (filename.endswith('.cif') or filename.endswith('.pdb')):
      filename += extension
    return filename

  def write_model_file(self, model_str, filename=Auto,
                       format=Auto, overwrite=Auto,
                       output_cs=True):
    '''
    Function for writing a model to file

    Parameters
    ----------
      model_str: str or mmtbx.model.manager object
        The string to be written or a model object. If a model object is
        provided, the format (PDB or mmCIF) of the original file is kept
        unless specified with format or target_output_format below.
        If a string is provided, the format must be specified as pdb or cif
      filename: str or Auto
        The output filename. If set to Auto, a default filename is
        generated based on params.output.prefix, params.output.suffix,
        and params.output.serial
      format: pdb or cif (mmcif treated as cif) or Auto.  If set to
         Auto, defaults to format of original file.
         If self._target_output_format is not None,
         always write model objects to this format if possible.
      overwrite: bool or Auto
        Overwrite filename if it exists. If set to Auto, the overwrite
        state of the DataManager is used.
      output_cs: bool
        Defines if crystal symmetry needs to be outputted. Passed directly
        into model_as_mmcif() and model_as_pdb()

    Returns
    -------
      filename: str
        The actual output filename. This may differ from the
        get_default_output_model_filename function since that sets the
        extension to cif by default. This function may alter the extension
        based on the desired format.
    '''


    if format == 'mmcif': format = 'cif'  # mmcif and cif are synonyms here

    if isinstance(model_str, mmtbx.model.manager):

      # Get overall preference for output format
      if (format is Auto) and hasattr(self,'_target_output_format') and (
           self._target_output_format is not None):
        format = self._target_output_format

      # Write as mmCIF if:
      #  1. format was supplied as 'cif' or
      #  2. format was Auto and target_output_format was set to 'cif', or
      #  3. format was Auto, no target_output_format set and this model was
      #       cif when read in, or
      #  4. model does not fit in PDB format
      if (format == 'cif') or (format is Auto and
            model_str.input_model_format_cif()) or (
          not model_str.get_hierarchy().fits_in_pdb_format()):
        extension = '.cif'
        format = 'cif'
        model_str = model_str.model_as_mmcif(output_cs=output_cs)
      else:
        extension = '.pdb'
        format = 'pdb'
        model_str = model_str.model_as_pdb(output_cs=output_cs)

    else:  # writing a string. Output format must be specified on input

      if format == 'cif':
        extension = '.cif'
      elif format == 'pdb':
        extension = '.pdb'
      else:  # no extension
        extension = Auto

    # Get the filename and add extension if necessary
    if filename is Auto:
      filename = self.get_default_output_model_filename(extension=extension)
    elif extension and (extension is not Auto) and (
        not extension_matches_ending(filename, extension)):
      other_extension = ".pdb" if extension == ".cif" else ".cif"
      fn,ext = os.path.splitext(filename)
      if ext == other_extension: # swap extension
        filename = fn + extension
      elif extension_matches_ending(filename, other_extension):
        filename = fn + "%s_%s" %(extension, ext.split("_")[1])
      else:
        filename += extension
    if model_str:
      return self._write_text(ModelDataManager.datatype, model_str,
                            filename=filename, overwrite=overwrite)
    else:
      return ''

def extension_matches_ending(filename, extension):
  if not extension:
    return True
  fn, ext = os.path.splitext(filename)
  if ext == extension:
    return True
  if "_" in ext and ext.split("_")[0] == extension:
    return True
 # =============================================================================
# end
