import os, pickle

class env:

  def __init__(self):
    libtbx_build = os.environ["LIBTBX_BUILD"]
    file_name = os.path.join(libtbx_build, "libtbx_env")
    try:
      f = open(file_name)
    except:
      raise RuntimeError("Cannot open libtbx environment file: " + file_name)
    try:
      e = pickle.load(f)
    except:
      raise RuntimeError("Corrupt libtbx environment file: " + file_name)
    self.__dict__.update(e)
    assert self.LIBTBX_BUILD == libtbx_build

  def dist_path(self, package_name):
    return self.dist_paths[package_name.upper() + "_DIST"]

  def current_working_directory_is_libtbx_build(self):
    return os.path.normpath(os.getcwd()) == self.LIBTBX_BUILD

class _build_options:

  def set(self, **kw):
    assert env().current_working_directory_is_libtbx_build()
    self.__dict__.update(kw)

build_options = _build_options()
