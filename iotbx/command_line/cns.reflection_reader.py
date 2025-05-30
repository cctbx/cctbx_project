"""Read a CNS reflection file"""
from __future__ import absolute_import, division, print_function
from iotbx.cns import reflection_reader
import sys
if (__name__ == "__main__"):
  reflection_reader.run(sys.argv[1:])
