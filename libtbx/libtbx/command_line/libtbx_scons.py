#! /usr/bin/env python
import sys, os
libtbx_scons = os.environ["LIBTBX_SCONS"]
sys.path.insert(0, libtbx_scons)
from SCons import Script
Script.main()
