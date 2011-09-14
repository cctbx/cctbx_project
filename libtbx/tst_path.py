def exercise_move_old_create_new_directory():
  from libtbx.path import move_old_create_new_directory as mocnd
  import os
  mocnd("tmp_mocnd")
  assert len(os.listdir("tmp_mocnd")) == 0
  for i in xrange(3):
    mocnd("tmp_mocnd/a")
  assert sorted(os.listdir("tmp_mocnd")) == ["a", "a_001", "a_002"]
  open("tmp_mocnd/a_23", "w")
  mocnd("tmp_mocnd/a")
  assert sorted(os.listdir("tmp_mocnd")) == [
    "a", "a_001", "a_002", "a_024", "a_23"]
  open("tmp_mocnd/a_log", "w")
  mocnd("tmp_mocnd/a")
  assert sorted(os.listdir("tmp_mocnd")) == [
    "a", "a_001", "a_002", "a_024", "a_025", "a_23", "a_log"]
  mocnd("tmp_mocnd/b")
  assert sorted(os.listdir("tmp_mocnd")) == [
    "a", "a_001", "a_002", "a_024", "a_025", "a_23", "a_log", "b"]
  mocnd("tmp_mocnd/b", serial_sep="", serial_fmt="%d")
  assert sorted(os.listdir("tmp_mocnd")) == [
    "a", "a_001", "a_002", "a_024", "a_025", "a_23", "a_log", "b", "b1"]

def run(args):
  assert len(args) == 0
  exercise_move_old_create_new_directory()
  from libtbx.path import random_new_directory_name
  assert len(random_new_directory_name()) == len("tmp_dir_00000000")
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
