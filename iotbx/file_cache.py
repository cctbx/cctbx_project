
# This module is used by the PHENIX GUI to manage file objects and associated
# phil parameters.

from libtbx.utils import hashlib_md5
from libtbx import adopt_init_args
import os

class manager (object) :
  file_type = None
  file_type_label = None
  def __init__ (self,
                allowed_param_names=None,
                allowed_multiple_params=None,
                debug=False,
                auto_reload_files=True,
                use_md5_sum=False) :
    adopt_init_args(self, locals())
    assert ((allowed_param_names is None) or
            (isinstance(allowed_param_names, list)))
    assert ((allowed_multiple_params is None) or
            (isinstance(allowed_param_names, list)))
    self.clear_cache()
    self._param_callbacks = {}

  def clear_cache (self) :
    self._file_mtimes = {}
    self._file_md5sums = {}
    self._cached_input_files = {}
    self._param_files = {}
    self.clear_format_specific_cache()

  def clear_format_specific_cache (self) :
    pass

  def set_param_callback (self, file_param_name, callback_handler) :
    assert hasattr(callback_handler, "__call__")
    self._param_callbacks[file_param_name] = callback_handler

  def add_file_callback (self, file_name) :
    pass

  def remove_file_callback (self, file_name) :
    pass

  def input_files (self) :
    for (file_name, file_object) in self._cached_input_files.iteritems() :
      yield (file_name, file_object)

  def save_file (self, input_file=None, file_name=None) :
    if (input_file is None) :
      from iotbx import file_reader
      input_file = file_reader.any_file(file_name)
    input_file.assert_file_type(self.file_type)
    file_name = input_file.file_name
    self._cached_input_files[file_name] = input_file
    self.add_file_callback(file_name)
    return self.save_other_file_data(input_file)

  def add_file (self, *args, **kwds) :
    return self.save_file(*args, **kwds)

  def save_other_file_data (self, input_file) :
    return None

  def file_is_modified (self, file_name) :
    if (not self.auto_reload_files) :
      return False
    elif (not file_name in self._cached_input_files) :
      return True
    elif self.use_md5_sum :
      file_records = open(file_name).read()
      m = hashlib_md5(file_records)
      old_md5sum = self._file_md5sums.get(file_name, None)
      if (old_md5sum is None) or (old_md5sum != m) :
        self._file_md5sums[file_name] = m
        return True
    else :
      mtime = os.path.getmtime(file_name)
      old_mtime = self._file_mtimes.get(file_name, None)
      if (old_mtime is None) or (old_mtime < mtime) :
        self._file_mtimes[file_name] = mtime
        return True
    return False

  def remove_file (self, file_name) :
    if (file_name in self._cached_input_files) :
      self._cached_input_files.pop(file_name)
      if self.use_md5_sum :
        self._file_md5sums.pop(file_name)
      else :
        self._file_mtimes.pop(file_name)
      if (self.allowed_param_names is not None) :
        for param_name, param_file in self._param_files.iteritems() :
          if self.allow_multiple(param_name) :
            if param_file.contains(file_name) :
              parma_file.remove(file_name)
              if (len(param_file) == 0 ) :
                self._param_files.pop(param_name)
              break
          elif (param_file == file_name) :
            self._param_files.pop(param_name)
            break
      self.remove_file_callback(file_name)

  def get_file (self, file_name=None, file_param_name=None) :
    from iotbx import file_reader
    if (file_name is None) and (file_param_name is not None) :
      file_name = self._param_files.get(file_param_name)
      if (isinstance(file_name, list)) :
        return file_name
    if (file_name is None) :
      return None
    assert os.path.isfile(file_name)
    if (file_name in self._cached_input_files) :
      if self.file_is_modified(file_name) :
        input_file = file_reader.any_file(file_name)
        self.save_file(input_file)
      return self._cached_input_files[file_name]
    return None

  def allow_multiple (self, param_name) :
    if (self.allowed_multiple_params is not None) :
      return (param_name in self.allowed_multiple_params)
    return False

  def get_current_file_names (self) :
    file_names = self._cached_input_files.keys()
    file_names.sort()
    return file_names

  def set_param_file (self, file_name, file_param_name, input_file=None,
      run_callback=True) :
    if self.allowed_param_names is not None :
      if not file_param_name in self.allowed_param_names :
        raise KeyError("Unrecognized input file parameter %s."%file_param_name)
    if (file_name is None) or (file_name == "") or (file_name == "None") :
      self._param_files.pop(file_param_name)
    else :
      if (input_file is not None) :
        self.save_file(input_file)
      elif (self.get_file(file_name) is None) :
        from iotbx import file_reader
        input_file = file_reader.any_file(file_name)
        self.save_file(input_file)
      if self.allow_multiple(file_param_name) :
        if (not file_param_name in self._param_files) :
          self._param_files[file_param_name] = []
        self._param_files[file_param_name].append(file_name)
      else :
          self._param_files[file_param_name] = file_name
    if run_callback :
      callback = self._param_callbacks.get(file_param_name, None)
      if (callback is not None) :
        callback(file_name)

  def unset_param_file (self, file_name, file_param_name, run_callback=True) :
    if (self.allow_multiple(file_param_name) and
        (file_param_name in self._param_files)) :
      param_files = self._param_files.get(file_param_name)
      if (file_name in param_files) :
        param_files.remove(file_name)
      if (len(param_files) == 0) :
        self._param_files.pop(file_param_name)
    else :
      if (self._param_files[file_param_name] == file_name) :
        self._param_files.pop(file_param_name)
    if run_callback :
      callback = self._param_callbacks.get(file_param_name, None)
      if (callback is not None) :
        callback(None)

  def get_file_params (self, file_name) :
    params = []
    for file_param_name in self._param_files :
      param_files = self._param_files[file_param_name]
      if isinstance(param_files, list) and (file_name in param_files) :
        params.append(file_param_name)
      elif (param_files == file_name) :
        params.append(file_param_name)
    return params

  def get_param_files (self, file_param_name) :
    file_name = self._param_files.get(file_param_name)
    if isinstance(file_name, list) :
      return file_name
    else :
      return [ file_name ]

  def get_file_type_label (self, file_name=None, input_file=None) :
    return self.file_type_label
