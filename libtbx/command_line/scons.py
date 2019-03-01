from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry, show_times_at_exit
from libtbx.str_utils import show_string
import libtbx.load_env
import sys, os

def find_scons_engine_path():
  join = os.path.join
  isdir = os.path.isdir
  if (libtbx.env.scons_dist_path is not None):
    result = libtbx.env.scons_dist_path / "engine"
    if result.isdir(): return abs(result)
    result = libtbx.env.scons_dist_path / "src" / "engine"
    if result.isdir(): return abs(result)
  for path in libtbx.env.repository_paths:
    result = path / "scons" / "engine"
    if result.isdir(): return abs(result)
    result = path / "scons" / "src" /"engine"
    if result.isdir(): return abs(result)
    # More recent scons-local packages don't use 'engine' directory
    result = path / "scons"
    if result.isdir(): return abs(result)
  return None

def dummy_fetch_win32_parallel_msg():
  pass

def run():
  debug_import = "--debug=import" in sys.argv[1:]
  def show_traceback():
    if (debug_import):
      import traceback
      print(file=sys.stderr)
      traceback.print_exc()
      print(file=sys.stderr)
  engine_path = find_scons_engine_path()
  if (engine_path is not None):
    sys.path.insert(0, engine_path)
    try: import SCons
    except ImportError:
      show_traceback()
      del sys.path[0]
  try: import SCons.Script
  except ImportError:
    show_traceback()
    msg = ["SCons is not available.",
      "  A possible solution is to unpack a SCons distribution in",
      "  one of these directories:"]
    for path in libtbx.env.repository_paths:
      msg.append("    " + show_string(abs(path)))
    msg.extend([
      "  SCons distributions are available at this location:",
      "    http://www.scons.org/",
      "  It may be necessary to rename the unpacked distribution, e.g.:",
      "    mv scons-0.96.1 scons"])
    raise Sorry("\n".join(msg))
  import SCons.Script.Main
  if (hasattr(SCons.Script.Main, "fetch_win32_parallel_msg")):
    SCons.Script.Main.fetch_win32_parallel_msg = dummy_fetch_win32_parallel_msg
  show_times_at_exit()
  SCons.Script.main()

if (__name__ == "__main__"):
  run()
