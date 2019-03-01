from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME libtbx.ipython
import re
import sys

from IPython import start_ipython

if __name__ == '__main__':
      sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
      sys.exit(start_ipython())

