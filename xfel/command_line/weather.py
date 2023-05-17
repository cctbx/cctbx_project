from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.weather

from xfel.util.weather import params_from_phil, run, message
import sys

if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(sys.argv[1:])
  run(params)
