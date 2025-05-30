"""Read a CNS SDB file"""
from __future__ import absolute_import, division, print_function
from iotbx.cns import sdb_reader
import sys
if (__name__ == "__main__"):
  sdb_reader.run(sys.argv[1:])
