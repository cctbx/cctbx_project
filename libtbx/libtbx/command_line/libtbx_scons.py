import sys, os

libtbx_dist = os.environ["LIBTBX_DIST"]
engine = os.path.normpath(os.path.join(libtbx_dist, "../scons/engine"))
if (not os.path.isdir(engine)):
  engine = os.path.normpath(os.path.join(libtbx_dist, "../scons/src/engine"))
sys.path.insert(0, engine)
try: import SCons
except: del sys.path[0]

from SCons import Script
Script.main()
