from __future__ import absolute_import, division

import os

import libtbx.load_env


def setup_hdf5_plugin_path():
    """Sets up the plugin path for HDF5.

    Pick up Dectris's LZ4 and bitshuffle compression plugins
    """
    plugin_path = libtbx.env.under_base(os.path.join("lib", "plugins"))
    try:
        os.environ["HDF5_PLUGIN_PATH"] = (
            plugin_path + os.pathsep + os.environ["HDF5_PLUGIN_PATH"]
        )
    except Exception:
        os.environ["HDF5_PLUGIN_PATH"] = plugin_path


setup_hdf5_plugin_path()
