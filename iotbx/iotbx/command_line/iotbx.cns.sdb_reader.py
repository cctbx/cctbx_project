#! /usr/bin/env python
from iotbx.cns import sdb_reader
import sys
sdb_reader.run(sys.argv[1:])
