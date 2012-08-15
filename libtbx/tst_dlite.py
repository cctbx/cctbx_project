from __future__ import division
from libtbx import dlite
from libtbx.utils import null_out
import time
import sys, os

def update_target(out):
  target_db = dlite.target_db(file_name="tmp_dlite")
  pair_info = target_db.pair_info(
    source_path="tmp_source",
    target_path="tmp_target")
  update_done = False
  if (pair_info.needs_update):
    pair_info.start_building_target()
    open(pair_info.target.path, "w").write(
      open(pair_info.source.path).read().upper())
    pair_info.done_building_target()
    update_done = True
  target_db.write()
  target_db.show(out=out)
  return update_done

def exercise(args, mtime_resolution=2):
  if ("--verbose" in args):
    out = None
  else:
    out = null_out()
  if (os.path.exists("tmp_dlite")):
    os.remove("tmp_dlite")
  print >> open("tmp_source", "w"), "a"
  assert update_target(out=out)
  assert not update_target(out=out)
  print >> open("tmp_source", "w"), "b"
  assert update_target(out=out)
  assert not update_target(out=out)
  time.sleep(mtime_resolution)
  assert not update_target(out=out)
  print >> open("tmp_source", "w"), "b"
  assert not update_target(out=out)
  print >> open("tmp_target", "w"), "B"
  assert not update_target(out=out)
  print >> open("tmp_target", "w"), "C"
  assert update_target(out=out)
  assert not update_target(out=out)
  os.remove("tmp_target")
  assert update_target(out=out)
  assert not update_target(out=out)
  print "OK"

if (__name__ == "__main__"):
  exercise(args=sys.argv[1:])
