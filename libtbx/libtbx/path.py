import os

def norm_join(*args):
  return os.path.normpath(apply(os.path.join, args))

def create_target_dir(target_file):
  target_dir = os.path.split(target_file)[0]
  if (not os.path.isdir(target_dir)):
    os.makedirs(target_dir)

_shortpath_bat = None

def abs_path_short(abs_path):
  if (os.name != "nt"): return abs_path
  global _shortpath_bat
  if (_shortpath_bat is None):
    _shortpath_bat = "shortpath.bat"
    libtbx_build = os.environ.get("LIBTBX_BUILD", None)
    if (libtbx_build is not None):
      _shortpath_bat = os.path.join(libtbx_build, _shortpath_bat)
  assert os.path.exists(_shortpath_bat)
  return os.popen('call "%s" "%s"' %
    (_shortpath_bat, abs_path), "r").readline().rstrip()

def abs_path_clean(abs_path):
  if (os.name != "nt"): return abs_path
  short = abs_path_short(abs_path).split(os.sep)
  orig = abs_path.split(os.sep)
  clean = []
  for o,s in zip(orig, short):
    if (o.find(" ") < 0):
      clean.append(o)
    else:
      clean.append(s)
  return os.sep.join(clean)
