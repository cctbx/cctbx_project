from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME sphinx.build

import sys
try:
  # try importing scipy.linalg before any cctbx modules, otherwise we
  # sometimes get a segmentation fault/core dump if it is imported after.
  # scipy.linalg is a dependency of e.g. xfel/clustering
  import scipy.linalg # import dependency
except ImportError:
  pass

if __name__ == '__main__':
    try:
        from sphinx import main
        sys.exit(main(sys.argv))
    except ImportError:
        from sphinx.cmd.build import main
        sys.exit(main())
