from __future__ import absolute_import, division, print_function
import os
import sys
import boost_adaptbx.boost.python as bp
try:
  ext = bp.import_ext('fast_linalg_ext')
  env = ext.env
  try_to_initialise = True
except ImportError:
  env = lambda: None
  env.initialised = False
  try_to_initialise = False

def find_dlls(lib_dirs):
  files = []
  for lib_dir in lib_dirs:
    if not os.path.exists(lib_dir):
      continue
    files += [os.path.join(lib_dir, x)\
      for x in os.listdir(lib_dir)\
        if 'openblas' in x and (x.endswith('.dll') or x.endswith('.dylib') or x.endswith('.so') or '.so.' in x)]
  return files

def find_old_layout_libs():
  import numpy, scipy
  dir_name = ".dylibs" if sys.platform == 'darwin' else ".libs"
  dirs = [os.path.join(os.path.dirname(scipy.__file__), dir_name),
          os.path.join(os.path.dirname(numpy.__file__), dir_name)]
  return find_dlls(dirs)

def find_new_layout_libs():
  import numpy, scipy
  from pathlib import Path
  dir_name = ".dylibs" if sys.platform == 'darwin' else ".libs"
  dirs = [os.path.join(Path(os.path.dirname(scipy.__file__)).parent.absolute(), "scipy%s" %dir_name),
          os.path.join(Path(os.path.dirname(numpy.__file__)).parent.absolute(), "numpy%s" %dir_name)]
  return find_dlls(dirs)

def find_extra_dll_libs():
  """ OpenBLAS in scipy/extra-dll, an older Windows layout

  Covered by neither find_old_layout_libs (scipy/.libs) nor
  find_new_layout_libs (a sibling scipy.libs).
  """
  import scipy
  return find_dlls([os.path.join(os.path.dirname(scipy.__file__),
                                 "extra-dll")])

def candidate_libraries():
  """ Shared libraries to try, in order of cost

  A generator, so that the numpy and scipy imports done by the finders are only
  paid for if the plain name on the system search path does not work.
  """
  yield "openblas.dll" if sys.platform[:3] == "win" else "libopenblas.so"
  for finder in (find_new_layout_libs, find_old_layout_libs,
                 find_extra_dll_libs):
    try:
      for path in finder():
        yield path
    except Exception:
      continue

# initialise() has two distinct failure modes: it may raise, or it may return
# leaving initialised False when the library was found but its symbols did not
# resolve. Both mean this candidate is unusable, so both move on to the next.
if not env.initialised and try_to_initialise:
  tried = 0
  for lib_file in candidate_libraries():
    tried += 1
    try:
      env.initialise(lib_file)
    except Exception:
      continue
    if env.initialised:
      print("Successfully initialised OpenBlas from %s:" %lib_file)
      print(env.build_config)
      break
  if not env.initialised:
    print("Could not initialise OpenBlas: %d candidate(s) tried" %tried)
