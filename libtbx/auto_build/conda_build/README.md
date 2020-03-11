# Scripts for conda-build

1) `create_setup.py` uses `setuptools` to partially find and copy the Python files in `cctbx_project`
2) `install_modules.py` copies or links complete module directories to a location
3) `update_libtbx_env.py` updates the bookkeeping for modules installed into `$PREFIX` and updates the dispatchers.
