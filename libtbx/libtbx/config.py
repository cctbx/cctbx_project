import sys, os, pickle

class UserError(Exception): pass

class env:

  def __init__(self):
    libtbx_build = os.environ["LIBTBX_BUILD"]
    file_name = os.path.join(libtbx_build, "libtbx_env")
    try:
      f = open(file_name)
    except:
      raise UserError("Cannot open libtbx environment file: " + file_name)
    try:
      e = pickle.load(f)
    except:
      raise UserError("Corrupt libtbx environment file: " + file_name)
    self.__dict__.update(e)
    assert self.LIBTBX_BUILD == libtbx_build
    if (not self.python_version_major_minor == sys.version_info[:2]):
      raise UserError("Python version incompatible with this build."
       + " Version used to configure: %d.%d." % self.python_version_major_minor
       + " Version used now: %d.%d." % sys.version_info[:2])

  def dist_name(self, package_name):
    return package_name.upper() + "_DIST"

  def has_dist(self, package_name):
    return self.dist_name(package_name) in self.dist_paths

  def dist_path(self, package_name):
    return self.dist_paths[self.dist_name(package_name)]

  def current_working_directory_is_libtbx_build(self):
    return os.path.normpath(os.getcwd()) == self.LIBTBX_BUILD

class _build_options:

  def set(self, **kw):
    assert env().current_working_directory_is_libtbx_build()
    self.__dict__.update(kw)

build_options = _build_options()
