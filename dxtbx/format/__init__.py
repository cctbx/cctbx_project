from __future__ import division

# Pick up Dectris's LZ4 and bitshuffle compression plugins
def setup_hdf5_plugin_path():
  import os, libtbx.load_env
  plugin_path = os.path.abspath(os.path.join(abs(libtbx.env.python_exe), "..", "..", "lib"))
  try:
    os.environ['HDF5_PLUGIN_PATH'] = os.environ['HDF5_PLUGIN_PATH'] + ":" + plugin_path
  except Exception:
    os.environ['HDF5_PLUGIN_PATH'] = plugin_path
setup_hdf5_plugin_path()
