# LIBTBX_SET_DISPATCHER_NAME ptw
from __future__ import absolute_import, division, print_function

import sys

from pkg_resources import DistributionNotFound, load_entry_point

# modify sys.argv so the command line help shows the right executable name
sys.argv[0] = 'ptw'

if __name__ == '__main__':
  try:
    ptw = load_entry_point('pytest-watch>=4.1.0', 'console_scripts', 'ptw')

  except DistributionNotFound:
    # Install package if necessary
    import pip
    pip.main(['install', 'pytest-watch>=4.1.0'])
    ptw = load_entry_point('pytest-watch>=4.1.0', 'console_scripts', 'ptw')

  sys.exit(ptw())
