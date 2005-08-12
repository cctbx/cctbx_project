import libtbx.load_env
import time
import sys, os

time_start = time.time()

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
  s = "usr+sys time: %.2f seconds" % usr_plus_sys
  if (ticks is not None):
    s += ", ticks: %d" % ticks
    if (ticks != 0):
      s += ", micro-seconds/tick: %.3f" % (usr_plus_sys*1.e6/ticks)
  print s
  wall_clock_time = time.time() - time_start
  print "wall clock time:",
  if (wall_clock_time < 120):
    print "%.2f seconds" % wall_clock_time
  else:
    m = int(wall_clock_time / 60 + 1.e-6)
    s = wall_clock_time - m * 60
    print "%d minutes %.2f seconds (%.2f seconds total)" % (
      m, s, wall_clock_time)

def run():
  import atexit
  atexit.register(show_times)

  from SCons import Script
  Script.main()

if (__name__ == "__main__"):
  run()
