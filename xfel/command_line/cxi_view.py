# LIBTBX_SET_DISPATCHER_NAME cxi.view
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

from xfel.cxi.display_spots import view_raw_image
import sys,os

if (__name__ == "__main__"):
  files = [arg for arg in sys.argv[1:] if os.path.isfile(arg)]
  arguments = [arg for arg in sys.argv[1:] if not os.path.isfile(arg)]
  for file in files:
    view_raw_image(file, *arguments, **({'display':True}))
