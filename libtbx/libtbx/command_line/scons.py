#! /usr/bin/env python

import sys, os

dist_root = os.environ["LIBTBX_DIST_ROOT"]
engine = os.path.normpath(os.path.join(dist_root, "scons/engine"))
if (not os.path.isdir(engine)):
  engine = os.path.normpath(os.path.join(dist_root, "scons/src/engine"))
sys.path.insert(0, engine)
try: import SCons
except: del sys.path[0]

from SCons import Script
Script.main()
