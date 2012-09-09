from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME libtbx.ipython

import IPython
if (__name__ == "__main__"):
  try:
    IPython.Shell.start().mainloop()
  except AttributeError:
    from IPython import embed
    embed()
