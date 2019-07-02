from __future__ import absolute_import, division, print_function
import sys, os, time, string

def run_string(cmd_and_args):
  t0 = time.time()
  os.system(cmd_and_args)
  t1 = time.time()
  return t1 - t0

def run_argv(argv):
  return run_string(string.join(argv))

if (__name__ == "__main__"):
  print("u+s:", run_argv(sys.argv[1:]))
