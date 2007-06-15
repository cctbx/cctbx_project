from libtbx import easy_run
import libtbx.load_env
import sys, os

def run(args):
  assert len(args) == 0
  if (os.name == "nt"):
    exe_suffix = ".exe"
  else:
    exe_suffix = ""
  command = libtbx.env.under_build("scitbx/error/tst_error"+exe_suffix)
  assert os.path.isfile(command)
  command_out = libtbx.easy_run.fully_buffered(
    command=command).raise_if_errors().stdout_lines
  assert command_out[0].startswith("scitbx InternalError: ")
  assert command_out[1:] == ["  x = 1.1", "  n = 1", "  j = 2"]
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
