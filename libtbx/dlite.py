"light-weight, simple source_path, target_path dependency management"
from __future__ import absolute_import, division, print_function

from libtbx import easy_pickle
from libtbx.utils import Sorry
import hashlib
import time
import sys, os

def try_loading_db(file_name):
  '''
  This function tries to load an existing pickle file. If the pickle file
  cannot be loaded, an empty target_db object is returned. This helps
  simplify the rebuilding of chem_data databases when switching Python
  versions.

  Parameters
  ----------
  file_name : str
      The filename for the database file

  Returns
  -------
  target_db : target_db
      If the file exists and the pickle file can be loaded, the target_db
      object. If the pickle protocol is not supported, the file is removed
      and an empty target_db object is returned.
  '''
  db = None
  try:
    db = target_db(file_name)
  except (Sorry, ValueError) as s:
    if 'unsupported pickle protocol' in str(s):
      os.remove(file_name)
      db = target_db(file_name)
  assert db is not None
  return db

class node_info(object):

  def __init__(self, path):
    self.path = path
    self.mtime = None
    self.md5 = None

  def full_path(self, path_prefix=None):
    if (path_prefix is None): return self.path
    return os.path.join(path_prefix, self.path)

  def current_mtime(self, path_prefix=None):
    full_path = self.full_path(path_prefix=path_prefix)
    if (not os.path.exists(full_path)): return None
    return os.path.getmtime(full_path)

  def current_md5(self, path_prefix=None):
    full_path = self.full_path(path_prefix=path_prefix)
    if (not os.path.exists(full_path)): return None
    m = hashlib.md5()
    with open(full_path, "rb") as f:
      m.update(f.read())
    return m.hexdigest()

  def has_changed(self, path_prefix=None, mtime_resolution=2):
    old_mtime = self.mtime
    if (old_mtime is None): return True
    self.mtime = self.current_mtime(path_prefix=path_prefix)
    if (self.mtime == old_mtime
        and time.time() > old_mtime + mtime_resolution): return False
    if (self.md5 is None): return True
    old_md5 = self.md5
    self.md5 = self.current_md5(path_prefix=path_prefix)
    return self.md5 != old_md5

class pair_info(object):

  def __init__(self, source_path, target_path, needs_update=True):
    self.source = node_info(path=source_path)
    self.target = node_info(path=target_path)
    self.needs_update = needs_update

  def eval_needs_update(self, source_path=None, path_prefix=None):
    if (source_path != self.source.path):
      self.source = node_info(path=source_path)
      self.needs_update = True
    elif (not self.needs_update):
      if (   self.source.has_changed(path_prefix=path_prefix)
          or self.target.has_changed(path_prefix=path_prefix)):
        self.needs_update = True
    return self.needs_update

  def start_building_target(self, path_prefix=None):
    if (self.source.mtime is None):
      self.source.mtime = self.source.current_mtime(path_prefix=path_prefix)
    if (self.source.md5 is None):
      self.source.md5 = self.source.current_md5(path_prefix=path_prefix)

  def done_building_target(self, path_prefix=None):
    self.target.mtime = self.target.current_mtime(path_prefix=path_prefix)
    self.target.md5 = self.target.current_md5(path_prefix=path_prefix)
    self.needs_update = False

class target_db(object):

  def __init__(self, file_name, file_name_during_write=None):
    self.file_name = file_name
    if (file_name_during_write is None and self.file_name is not None):
      self.file_name_during_write = self.file_name + ".new"
    else:
      self.file_name_during_write = file_name_during_write
    if (self.file_name is None
        or not os.path.exists(self.file_name)):
      self.pair_infos = {}
    else:
      self.pair_infos = easy_pickle.load(file_name=self.file_name)

  def write(self):
    assert self.file_name is not None
    easy_pickle.dump(file_name=self.file_name_during_write, obj=self.pair_infos)
    if (os.path.exists(self.file_name)):
      os.remove(self.file_name)
    os.rename(self.file_name_during_write, self.file_name)

  def pair_info(self, source_path, target_path, path_prefix=None):
    result = self.pair_infos.get(target_path)
    if (result is None):
      result = pair_info(source_path=source_path, target_path=target_path)
      self.pair_infos[target_path] = result
    else:
      result.eval_needs_update(
        source_path=source_path, path_prefix=path_prefix)
    return result

  def show(self, out=None):
    if (out is None): out = sys.stdout
    for pair_info in self.pair_infos.values():
      for attr in ["source", "target"]:
        node = getattr(pair_info, attr)
        print(attr+":", node.path, "mtime:", node.mtime,\
                                           "md5:", node.md5, file=out)
      print("-"*79, file=out)
