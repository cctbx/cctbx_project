'''
DataManager datatype for reading and writing YAML files

The PyYAML import is deliberately lazy (see _yaml_module): this module is
loaded by every default DataManager and ships inside cctbx-base / Olex2, where
PyYAML may be absent. Do NOT hoist "import yaml" to module scope.
'''

from iotbx.data_manager import DataManagerBase
from libtbx import Auto
from libtbx.utils import Sorry

# =============================================================================
class YamlDataManager(DataManagerBase):
  '''
  DataManager datatype for YAML files
  '''

  datatype = 'yaml'

  # ---------------------------------------------------------------------------
  def _yaml_module(self):
    '''
    Lazily import PyYAML.

    Returns
    -------
    module
        The imported yaml module

    Raises
    ------
    Sorry
        If PyYAML is not installed (with installation guidance), instead
        of a raw ImportError
    '''
    try:
      import yaml
    except ImportError:
      raise Sorry('PyYAML is required to read or write YAML files. '
                  'Install it with "libtbx.conda install pyyaml" or '
                  '"libtbx.pip install PyYAML".')
    return yaml

  # ---------------------------------------------------------------------------
  # YAML
  def add_yaml(self, filename, data):
    '''
    Add a parsed YAML object to the DataManager under a filename key.

    Parameters
    ----------
    filename : str
        The filename key for the data
    data : object
        The parsed Python object (e.g. dict or list)
    '''
    return self._add(YamlDataManager.datatype, filename, data)

  def set_default_yaml(self, filename):
    '''
    Set the default YAML object of the DataManager.

    Parameters
    ----------
    filename : str
        The filename key of the object to use as the default
    '''
    return self._set_default(YamlDataManager.datatype, filename)

  def get_yaml(self, filename=None):
    '''
    Return the parsed Python object for a YAML file.

    Parameters
    ----------
    filename : str, optional
        The filename key of the stored object. If None, the default
        YAML object is returned.

    Returns
    -------
    object
        The parsed Python object (dict, list, scalar, or None for an
        empty file)
    '''
    return self._get(YamlDataManager.datatype, filename)

  def get_yaml_names(self):
    '''
    Return the list of managed YAML filenames.

    Returns
    -------
    names : list of str
        The filename keys of the stored YAML objects
    '''
    return self._get_names(YamlDataManager.datatype)

  def get_default_yaml_name(self):
    '''
    Return the filename of the default YAML object.

    Returns
    -------
    filename : str or None
        The filename key of the default YAML object, or None if unset
    '''
    return self._get_default_name(YamlDataManager.datatype)

  def remove_yaml(self, filename):
    '''
    Remove a YAML object from the DataManager.

    Parameters
    ----------
    filename : str
        The filename key of the object to remove
    '''
    return self._remove(YamlDataManager.datatype, filename)

  def has_yamls(self, expected_n=1, exact_count=False, raise_sorry=False):
    '''
    Check if the DataManager has the expected number of YAML objects.

    Parameters
    ----------
    expected_n : int, optional
        The expected number of YAML objects
    exact_count : bool, optional
        If True, the number of objects must match expected_n exactly,
        otherwise at least expected_n objects are required.
    raise_sorry : bool, optional
        If True, raise Sorry instead of returning False

    Returns
    -------
    bool
        True if the expected number of YAML objects is available

    Raises
    ------
    Sorry
        If raise_sorry is True and the count check fails
    '''
    return self._has_data(YamlDataManager.datatype, expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_yaml_file(self, filename):
    '''
    Parse a YAML file (yaml.safe_load) and store the resulting Python object.

    An empty file parses to None and is stored as None.

    Parameters
    ----------
    filename : str
        The filepath of the YAML file

    Returns
    -------
    filename : str
        The filename added to the DataManager (the contract expected by
        _get and _set_default)

    Raises
    ------
    Sorry
        If the file does not exist, PyYAML is not installed, or the file
        is not valid YAML
    '''
    # _process_file parses once via read_file (file_type known) and stores the
    # object (an empty YAML file -> None); read_file surfaces the missing-PyYAML
    # / missing-file / malformed Sorry and transparently decompresses compressed
    # YAML (.yaml.gz/.xz/.bz2/...) - yaml has no compression-aware reader of its own
    return self._process_file(YamlDataManager.datatype, filename)

  def get_default_output_yaml_filename(self):
    '''
    Return the default output filename with a .yaml extension.

    Returns
    -------
    filename : str
        The default output filename, with ".yaml" appended if missing
    '''
    filename = self.get_default_output_filename()
    if not filename.endswith('.yaml'):
      filename += '.yaml'
    return filename

  def write_yaml_file(self, yaml_object, filename=Auto, overwrite=Auto):
    '''
    Write YAML to file.

    Parameters
    ----------
    yaml_object : object or str
        A Python object (serialized with yaml.safe_dump) or a
        pre-formatted string (written verbatim)
    filename : str, optional
        The output filename. If Auto, the default output filename with a
        ".yaml" extension is used.
    overwrite : bool, optional
        If True, an existing file is overwritten. If Auto, the overwrite
        setting of the DataManager is used.

    Returns
    -------
    filename : str
        The filename of the written file

    Raises
    ------
    Sorry
        If PyYAML is not installed and serialization is required, if the
        file already exists and overwrite is False, or if writing fails
    '''
    if filename is Auto:
      filename = self.get_default_output_yaml_filename()
    if isinstance(yaml_object, str):
      text = yaml_object
    else:
      yaml = self._yaml_module()
      text = yaml.safe_dump(yaml_object, default_flow_style=False)
    return self._write_text(YamlDataManager.datatype, text,
                            filename=filename, overwrite=overwrite)

# =============================================================================
# end
