import libtbx.load_env
import sys, os

engine = os.path.join(libtbx.env.scons_dist_path, "engine")
if (not os.path.isdir(engine)):
  engine = os.path.join(libtbx.env.scons_dist_path, "src", "engine")
if (os.path.isdir(engine)):
  sys.path.insert(0, engine)
  try: import SCons
  except: del sys.path[0]

def show_times():
  t = os.times()
  usr_plus_sys = t[0] + t[1]
  try: ticks = sys.gettickeraccumulation()
  except: ticks = None
  s = "usr+sys time: %.2f" % usr_plus_sys
  if (ticks is not None):
    s += ", ticks: %d" % ticks
    if (ticks != 0):
      s += ", micro-seconds/tick: %.3f" % (usr_plus_sys*1.e6/ticks)
  print s

import atexit
atexit.register(show_times)

from SCons import Script
Script.main()
