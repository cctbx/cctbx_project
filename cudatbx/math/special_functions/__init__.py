import boost.python
import libtbx.load_env
if (libtbx.env.build_options.enable_cuda):
  cudatbx_special_functions_ext = boost.python.import_ext("cudatbx_special_functions_ext")
  from cudatbx_special_functions_ext import *
