from __future__ import absolute_import, division

# Pick up Dectris's LZ4 and bitshuffle compression plugins
def setup_hdf5_plugin_path():
    import os, libtbx.load_env

    plugin_path = libtbx.env.under_base("lib/plugins")
    try:
        os.environ["HDF5_PLUGIN_PATH"] = (
            plugin_path + ":" + os.environ["HDF5_PLUGIN_PATH"]
        )
    except Exception:
        os.environ["HDF5_PLUGIN_PATH"] = plugin_path


setup_hdf5_plugin_path()
