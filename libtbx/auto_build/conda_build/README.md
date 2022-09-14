# Scripts for conda-build

1) `install_build.py` copies or links complete module directories to a location
2) `update_libtbx_env.py` updates the bookkeeping for modules installed into `$PREFIX` and updates the dispatchers.
3) `write_env_files.py` creates the <program>_env.{c}sh or bat files for setting
up a user's installation

# Not used
1) `create_setup.py` uses `setuptools` to partially find and copy the Python files in `cctbx_project`

The installation process is still being developed, but generally, the
`install_build.py` script copies (or links) the necessary files into the
`$PREFIX` directory. The `update_libtbx_env.py` script is called by
`libtbx.python` to do the initial copy of the environment. A second call
to the script by `$PREFIX/bin/python` will update the dispatchers.
