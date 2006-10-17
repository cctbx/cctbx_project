"light-weight, simple source_path, target_path dependency management"

from libtbx import easy_pickle
import md5
import time
import sys, os

class node_info(object):

  def __init__(self, path):
    self.path = path
    self.mtime = None
    self.md5 = None

  def path_exists(self):
    return os.path.exists(self.path)

  def current_mtime(self):
    if (not self.path_exists()): return None
    return os.path.getmtime(self.path)

  def current_md5(self):
    if (not self.path_exists()): return None
    m = md5.new()
    m.update(open(self.path, "rb").read())
    return m.hexdigest()

  def has_changed(self, mtime_resolution=2):
    old_mtime = self.mtime
    if (old_mtime is None): return True
    self.mtime = self.current_mtime()
    if (self.mtime == old_mtime
        and time.time() > old_mtime + mtime_resolution): return False
    if (self.md5 is None): return True
    old_md5 = self.md5
    self.md5 = self.current_md5()
    return self.md5 != old_md5

class pair_info(object):

  def __init__(self, source_path, target_path, needs_update=True):
    self.source = node_info(path=source_path)
    self.target = node_info(path=target_path)
    self.needs_update = needs_update

  def eval_needs_update(self, source_path=None):
    if (source_path != self.source.path):
      self.source = node_info(path=source_path)
      self.needs_update = True
    elif (not self.needs_update):
      if (   self.source.has_changed()
          or self.target.has_changed()):
        self.needs_update = True
    return self.needs_update

  def start_building_target(self):
    if (self.source.mtime is None):
      self.source.mtime = self.source.current_mtime()
    if (self.source.md5 is None):
      self.source.md5 = self.source.current_md5()

  def done_building_target(self):
    self.target.mtime = self.target.current_mtime()
    self.target.md5 = self.target.current_md5()
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

  def pair_info(self, source_path, target_path):
    result = self.pair_infos.get(target_path)
    if (result is None):
      result = pair_info(source_path=source_path, target_path=target_path)
      self.pair_infos[target_path] = result
    else:
      result.eval_needs_update(source_path=source_path)
    return result

  def show(self, out=None):
    if (out is None): out = sys.stdout
    for pair_info in self.pair_infos.values():
      for attr in ["source", "target"]:
        node = getattr(pair_info, attr)
        print >> out, attr+":", node.path, "mtime:", node.mtime,\
                                           "md5:", node.md5
      print >> out, "-"*79
