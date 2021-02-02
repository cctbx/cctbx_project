from __future__ import absolute_import, division, print_function
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
    with open(pair_info.target.path, "w") as fw, open(pair_info.source.path) as fr:
      fw.write(fr.read().upper())
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
  with open("tmp_source", "w") as f:
    print("a", file=f)
  assert update_target(out=out)
  assert not update_target(out=out)
  with open("tmp_source", "w") as f:
    print("b", file=f)
  assert update_target(out=out)
  assert not update_target(out=out)
  time.sleep(mtime_resolution)
  assert not update_target(out=out)
  with open("tmp_source", "w") as f:
    print("b", file=f)
  assert not update_target(out=out)
  with open("tmp_target", "w") as f:
    print("B", file=f)
  assert not update_target(out=out)
  with open("tmp_target", "w") as f:
    print("C", file=f)
  assert update_target(out=out)
  assert not update_target(out=out)
  os.remove("tmp_target")
  assert update_target(out=out)
  assert not update_target(out=out)
  print("OK")

if (__name__ == "__main__"):
  exercise(args=sys.argv[1:])
