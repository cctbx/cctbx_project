from __future__ import absolute_import, division, print_function
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
  dirs = [os.path.join(os.path.dirname(scipy.__file__), ".libs"),
          os.path.join(os.path.dirname(numpy.__file__), ".libs")]
  return find_dlls(dirs)

def find_new_layout_libs():
  import numpy, scipy
  from pathlib import Path
  dirs = [os.path.join(Path(os.path.dirname(scipy.__file__)).parent.absolute(), "scipy.libs"),
          os.path.join(Path(os.path.dirname(numpy.__file__)).parent.absolute(), "numpy.libs")]
  return find_dlls(dirs)

if not env.initialised and try_to_initialise:
  if sys.platform[:3] == "win":
    lib_path = "openblas.dll"
  else:
    lib_path = "libopenblas.so"
  try:
    env.initialise(lib_path)
    if env.initialised:
      print("Successfully initialised OpenBlas at %s:" %lib_path)
      print(env.build_config)
    else:
      print("Located OpenBlas but could not initialise")
  except Exception:
    import os
    try:
      files = find_new_layout_libs()
      if not files:
        files = find_old_layout_libs()
      if not files:
        print("Could not locate usable OpenBlas")
      for lib_file in files:
        try:
          env.initialise(lib_file.encode("utf-8"))
        except:
          continue
        if env.initialised:
          print("Successfully initialised OpenBlas from %s:" %lib_file)
          print(env.build_config)
          break
    except Exception as e:
      print("Could not initialise OpenBlas: %s" %e)
