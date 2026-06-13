from __future__ import absolute_import, division, print_function
'''
DataManager datatype for reading and writing JSON files
'''

import json

from iotbx.data_manager import DataManagerBase
from libtbx import Auto

# =============================================================================
class JsonDataManager(DataManagerBase):
  '''
  DataManager datatype for JSON files
  '''

  datatype = 'json'

  # ---------------------------------------------------------------------------
  # JSON
  def add_json(self, filename, data):
    '''
    Add a parsed JSON object to the DataManager under a filename key.

    Parameters
    ----------
    filename : str
        The filename key for the data
    data : object
        The parsed Python object (e.g. dict or list)
    '''
    return self._add(JsonDataManager.datatype, filename, data)

  def set_default_json(self, filename):
    '''
    Set the default JSON object of the DataManager.

    Parameters
    ----------
    filename : str
        The filename key of the object to use as the default
    '''
    return self._set_default(JsonDataManager.datatype, filename)

  def get_json(self, filename=None):
    '''
    Return the parsed Python object for a JSON file.

    Parameters
    ----------
    filename : str, optional
        The filename key of the stored object. If None, the default
        JSON object is returned.

    Returns
    -------
    object
        The parsed Python object (dict, list, or scalar)
    '''
    return self._get(JsonDataManager.datatype, filename)

  def get_json_names(self):
    '''
    Return the list of managed JSON filenames.

    Returns
    -------
    names : list of str
        The filename keys of the stored JSON objects
    '''
    return self._get_names(JsonDataManager.datatype)

  def get_default_json_name(self):
    '''
    Return the filename of the default JSON object.

    Returns
    -------
    filename : str or None
        The filename key of the default JSON object, or None if unset
    '''
    return self._get_default_name(JsonDataManager.datatype)

  def remove_json(self, filename):
    '''
    Remove a JSON object from the DataManager.

    Parameters
    ----------
    filename : str
        The filename key of the object to remove
    '''
    return self._remove(JsonDataManager.datatype, filename)

  def has_jsons(self, expected_n=1, exact_count=False, raise_sorry=False):
    '''
    Check if the DataManager has the expected number of JSON objects.

    Parameters
    ----------
    expected_n : int, optional
        The expected number of JSON objects
    exact_count : bool, optional
        If True, the number of objects must match expected_n exactly,
        otherwise at least expected_n objects are required.
    raise_sorry : bool, optional
        If True, raise Sorry instead of returning False

    Returns
    -------
    bool
        True if the expected number of JSON objects is available

    Raises
    ------
    Sorry
        If raise_sorry is True and the count check fails
    '''
    return self._has_data(JsonDataManager.datatype, expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_json_file(self, filename):
    '''
    Parse a JSON file and store the resulting Python object.

    Parameters
    ----------
    filename : str
        The filepath of the JSON file

    Returns
    -------
    filename : str
        The filename added to the DataManager (the contract expected by
        _get and _set_default)

    Raises
    ------
    Sorry
        If the file does not exist or is not valid JSON
    '''
    # _process_file parses once via read_file (file_type known) and stores the
    # object; read_file also transparently decompresses compressed JSON
    # (.json.gz/.xz/.bz2/...) - json has no compression-aware reader of its own
    return self._process_file(JsonDataManager.datatype, filename)

  def get_default_output_json_filename(self):
    '''
    Return the default output filename with a .json extension.

    Returns
    -------
    filename : str
        The default output filename, with ".json" appended if missing
    '''
    filename = self.get_default_output_filename()
    if not filename.endswith('.json'):
      filename += '.json'
    return filename

  def write_json_file(self, json_object, filename=Auto, overwrite=Auto):
    '''
    Write JSON to file.

    Parameters
    ----------
    json_object : object or str
        A Python object (serialized with json.dumps) or a pre-formatted
        string (written verbatim)
    filename : str, optional
        The output filename. If Auto, the default output filename with a
        ".json" extension is used.
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
        If the file already exists and overwrite is False, or if writing
        fails
    '''
    if filename is Auto:
      filename = self.get_default_output_json_filename()
    if isinstance(json_object, str):
      text = json_object
    else:
      text = json.dumps(json_object, indent=2)
    return self._write_text(JsonDataManager.datatype, text,
                            filename=filename, overwrite=overwrite)

# =============================================================================
# end
