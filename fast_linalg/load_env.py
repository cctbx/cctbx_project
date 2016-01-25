import fast_linalg

# The following import cannot live in fast_linalg/__init__.py because
# the code in the latter is executed when fast_linalg.precompiled is imported
# during the processing of libtbx_refresh.py, at which time fast_linalg_ext does
# not exist yet.
import boost.python
ext = boost.python.import_ext('fast_linalg_ext')
fast_linalg.env = fast_linalg_ext.env
