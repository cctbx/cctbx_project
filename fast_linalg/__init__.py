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
    if sys.platform[:3] == "win":
      import os
      try:
        import scipy
        scipy_path = os.path.dirname(scipy.__file__)
        lib_dirs = ["extra-dll", ".libs"]
        files = None
        for lib_dir in lib_dirs:
          scipy_libs_path = os.path.join(scipy_path, lib_dir)
          if not os.path.exists(scipy_libs_path):
            continue
          files = [os.path.join(x, scipy_libs_path)\
            for x in os.listdir(scipy_libs_path) if 'openblas' in x and 'dll' in x]
          break
        if files:
          env.initialise(files[0].encode("utf-8"))
          if env.initialised:
            print("Successfully initialised SciPy OpenBlas:")
            print(env.build_config)
      except Exception as e:
        print("Could not initialise OpenBlas: %s" %e)
    else:
      print("Could not initialise OpenBlas")
