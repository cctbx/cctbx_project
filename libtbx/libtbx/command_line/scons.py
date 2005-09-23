from libtbx.utils import Sorry
from libtbx.str_utils import show_string
import libtbx.load_env
import time
import sys, os

time_start = time.time()

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

def find_scons_engine_path():
  join = os.path.join
  isdir = os.path.isdir
  if (libtbx.env.scons_dist_path is not None):
    result = join(libtbx.env.scons_dist_path, "engine")
    if (isdir(result)): return result
    result = join(libtbx.env.scons_dist_path, "src", "engine")
    if (isdir(result)): return result
  for path in libtbx.env.repository_paths:
    result = join(path, "scons", "engine")
    if (isdir(result)): return result
    result = join(path, "scons", "src", "engine")
    if (isdir(result)): return result
  return None

def run():
  engine_path = find_scons_engine_path()
  if (engine_path is not None):
    sys.path.insert(0, engine_path)
    try: import SCons
    except ImportError: del sys.path[0]
  try: import SCons.Script
  except ImportError:
    msg = ["SCons is not available.",
      "  A possible solution is to unpack a SCons distribution in",
      "  one of these directories:"]
    for path in libtbx.env.repository_paths:
      msg.append("    " + show_string(path))
    msg.extend([
      "  SCons distributions are available at this location:",
      "    http://www.scons.org/",
      "  It may be necessary to rename the unpacked distribution, e.g.:",
      "    mv scons-0.96.1 scons"])
    raise Sorry("\n".join(msg))
  import atexit
  atexit.register(show_times)
  SCons.Script.main()

if (__name__ == "__main__"):
  run()
