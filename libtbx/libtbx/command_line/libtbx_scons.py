#! /usr/bin/env python
import sys, os
try: libtbx_scons = os.environ["LIBTBX_SCONS"]
except: pass
else:
  if (libtbx_scons != "default"):
    sys.path.insert(0, libtbx_scons)
from SCons import Script
Script.main()
