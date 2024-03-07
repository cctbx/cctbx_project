# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.preference

from __future__ import division
import sys
from xfel.util.preference import params_from_phil, phil_scope, run, message


if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(phil_scope, sys.argv[1:])
  run(params)
