import boost.python
import libtbx.load_env
if (libtbx.env.build_options.enable_cuda):
  cudatbx_scattering_ext = boost.python.import_ext("cudatbx_scattering_ext")
  from cudatbx_scattering_ext import *
