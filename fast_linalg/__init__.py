import sys
import boost.python
ext = boost.python.import_ext('fast_linalg_ext')
env = ext.env
if not env.initialised:
  if sys.platform[:3] == "win":
    lib_path = "openblas.dll"
  else:
    lib_path = "libopenblas.so"
  print("Trying OpenBlas at: %s" %lib_path)
  try:
    env.initialise(lib_path)
    if env.initialised:
      print("Successfully initialised OpenBlas:")
      print(env.build_config)
    else:
      print("Located OpenBlas but could not initialise")
  except:
    print("Could not initialise OpenBlas")