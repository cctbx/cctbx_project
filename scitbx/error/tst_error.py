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
  assert command.find('"') < 0
  command_out = libtbx.easy_run.fully_buffered(
    command='"'+command+'"').raise_if_errors().stdout_lines
  assert command_out[0].startswith("scitbx Internal Error: ")
  assert command_out[1:4] == ["  x = 1.1", "  n = 1", "  j = 2"]
  assert command_out[4].startswith(
    "Control flow passes through branch that should be unreachable: ")
  assert len(command_out) == 5
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
