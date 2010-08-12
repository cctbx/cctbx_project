
from libtbx.utils import hashlib_md5
import os

class manager (object) :
  file_type = None
  def __init__ (self,
                allowed_param_names=None,
                allowed_multiple_params=None,
                debug=False,
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

  def save_file (self, input_file) :
    input_file.assert_file_type(self.file_type)
    file_name = input_file.file_name
    self._cached_input_files[file_name] = input_file
    return self.save_other_file_data(input_file)

  def save_other_file_data (self, input_file) :
    return None

  def file_is_modified (self, file_name) :
    if self.use_md5_sum :
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
        if self.allow_multiple(file_param_name) :
          for param_name, file_list in self._param_files.iteritems() :
            if file_list.contains(file_name) :
              file_list.remove(file_name)
              break
        else :
          for param_name, param_file in self._param_files.iteritems() :
            if (param_file == file_name) :
              self._param_files.pop(param_name)
              break

  def get_file (self, file_name) :
    from iotbx import file_reader
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

  def set_param_file (self, file_name, file_param_name, input_file=None,
      run_callback=True) :
    if self.allowed_param_names is not None :
      if not file_param_name in self.allowed_param_names :
        raise KeyError("Unrecognized input file parameter %s."%file_param_name)
    if (file_name is None) or (file_name == "") or (file_name == "None") :
      self._param_files.pop(file_param_name)
    else :
      input_file = self.get_file(file_name)
      if (input_file is None) :
        self.save_file(file_name)
      if (self.allow_multiple(file_param_name) and
          (file_param_name in self._param_files)) :
          self._param_files[file_param_name].append(file_name)
      else :
          self._param_files[file_param_name] = file_name
    if run_callback :
      callback = self._param_callbacks.get(file_param_name, None)
      if (callback is not None) :
        callback(file_name)
